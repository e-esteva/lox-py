#!/usr/bin/env python3
"""
LoxCode Alignment Tool - A Python utility for processing paired-end sequencing LoxCode data.

Refactored into loxcode_align.core for importable use.
"""

import glob
import gzip
import os
import sys
from collections import defaultdict
from typing import List, Dict, Tuple, Optional

try:
    import Levenshtein
    HAS_LEVENSHTEIN = True
except ImportError:
    HAS_LEVENSHTEIN = False
    print("Warning: python-Levenshtein not found. Using basic edit distance implementation.")

try:
    from Bio import SeqIO
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

# --- your functions unchanged ---
def edit_distance(s1: str, s2: str) -> int:
    """Calculate edit distance using Levenshtein if available, otherwise use basic implementation."""
    if HAS_LEVENSHTEIN:
        return Levenshtein.distance(s1, s2)
    else:
        # Basic edit distance implementation
        if len(s1) < len(s2):
            return edit_distance(s2, s1)
        
        if len(s2) == 0:
            return len(s1)
        
        previous_row = list(range(len(s2) + 1))
        for i, c1 in enumerate(s1):
            current_row = [i + 1]
            for j, c2 in enumerate(s2):
                insertions = previous_row[j + 1] + 1
                deletions = current_row[j] + 1
                substitutions = previous_row[j] + (c1 != c2)
                current_row.append(min(insertions, deletions, substitutions))
            previous_row = current_row
        
        return previous_row[-1]


def read_fastq_manual(file_handle) -> Optional[Tuple[str, str, str, str]]:
    """Manually read a FASTQ record (4 lines)."""
    try:
        header = file_handle.readline().strip()
        if not header:
            return None
        sequence = file_handle.readline().strip()
        plus = file_handle.readline().strip()
        quality = file_handle.readline().strip()
        
        # Validate we have all required lines
        if not sequence or not plus or not quality:
            return None
            
        return header, sequence, plus, quality
    except Exception as e:
        return None


def open_file(filename: str):
    """Open file handle, automatically detecting gzip compression."""
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')


def read_fastq_records(filename: str):
    """Generator to read FASTQ records, supporting both regular and gzipped files."""
    if HAS_BIOPYTHON:
        try:
            with open_file(filename) as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    yield str(record.seq)
        except Exception as e:
            print(f"BioPython parsing failed for {filename}: {e}, falling back to manual parsing")
            # Fall back to manual parsing
            try:
                with open_file(filename) as f:
                    while True:
                        record = read_fastq_manual(f)
                        if record is None:
                            break
                        yield record[1]  # sequence
            except Exception as e2:
                print(f"Manual parsing also failed for {filename}: {e2}")
                return
    else:
        try:
            with open_file(filename) as f:
                while True:
                    record = read_fastq_manual(f)
                    if record is None:
                        break
                    yield record[1]  # sequence
        except Exception as e:
            print(f"Error reading {filename}: {e}")
            return


def consensus(a: List[int], b: List[int]) -> List[int]:
    """
    Create consensus sequence from two barcode vectors.
    """
    max_overlap = 0
    max_offset = 0
    asize, bsize = len(a), len(b)
    
    # Find best overlap
    for offset in range(-bsize, asize):
        match = 0
        not_match = 0
        
        for k in range(len(b)):
            if k + offset < 0 or k + offset >= asize:
                continue
            
            if a[k + offset] == b[k] and b[k] != 0:
                match += 1
            if a[k + offset] != b[k] and a[k + offset] != 0 and b[k] != 0:
                not_match += 1
        
        if match > max_overlap and not_match == 0:
            max_overlap = match
            max_offset = offset
    
    # Handle special case for non-overlapping sequences
    if max_overlap == 0:
        if (len(a) + len(b) == 13 and 
            a.count(0) == 0 and b.count(0) == 0):
            # Check for intersection
            a_sorted = sorted(a)
            b_sorted = sorted(b)
            intersection = set(a_sorted) & set(b_sorted)
            if len(intersection) > 0:
                return []
            
            # Return concatenated sequence
            return a + b
        else:
            return []
    
    # Build consensus
    consensus_seq = []
    
    # Add prefix from sequence a
    for i in range(max_offset):
        consensus_seq.append(a[i])
    
    # Add overlapping region
    for k in range(len(b)):
        if k + max_offset < len(a):
            if a[k + max_offset] == b[k]:
                consensus_seq.append(b[k])
            elif a[k + max_offset] != 0 and b[k] == 0:
                consensus_seq.append(a[k + max_offset])
            elif a[k + max_offset] == 0 and b[k] != 0:
                consensus_seq.append(b[k])
            else:
                consensus_seq.append(b[k])  # Default to b in case of conflict
        else:
            consensus_seq.append(b[k])
    
    # Handle negative offset
    if max_offset < 0:
        consensus_seq = consensus_seq[-max_offset:]
    
    return consensus_seq


def align_read(read: str, elements: Dict[str, int]) -> List[int]:
    """
    Align a read and extract LoxCode elements.
    """
    element_vector = []
    
    # LoxP sequences
    LOX_F = "ATAACTTCGTATAATGTATGCTATACGAAGTTAT"
    LOX_R = "ATAACTTCGTATAGCATACATTATACGAAGTTAT"
    
    # Parse read structure
    lox_sites = []
    element_seqs = []
    
    l = 0
    s = 8
    i = 0
    
    while True:
        if 34 + l + s > len(read):
            break
            
        lox_sites.append(read[l:l+34])
        element_seqs.append(read[34+l:34+l+s])
        
        if i % 2 == 0:
            l += 42
            s = 14
        else:
            l += 48
            s = 8
        i += 1
    
    # Process each LoxP site and element
    for k in range(len(lox_sites)):
        # Check LoxP site quality
        dist_f = edit_distance(LOX_F, lox_sites[k])
        dist_r = edit_distance(LOX_R, lox_sites[k])
        
        if min(dist_f, dist_r) < 10:
            # Try exact match first
            element_seq = element_seqs[k]
            if element_seq in elements and len(element_seq) == len(element_seq):
                element_vector.append(elements[element_seq])
                continue
            
            # Try approximate match
            found = False
            exclude_elements = {1, 3, 5, 2, 8}
            
            for ele_seq, ele_id in elements.items():
                dist = edit_distance(element_seq, ele_seq)
                if dist == 1 and ele_id not in exclude_elements:
                    element_vector.append(ele_id)
                    found = True
                    break
            
            if not found:
                element_vector.append(0)
        else:
            element_vector.append(0)
    
    return element_vector


def process_file_pair(r1_file: str, r2_file: str, output_file: str, elements: Dict[str, int]) -> float:
    """
    Process a pair of FASTQ files and generate barcode counts using streaming approach.
    """
    lox_counts = defaultdict(lambda: {'first_seen': 0, 'count': 0})
    total_reads = 0
    barcode_reads = 0
    
    try:
        # Check if files exist
        if not os.path.exists(r1_file):
            print(f"R1 file not found: {r1_file}")
            return 0.0
        if not os.path.exists(r2_file):
            print(f"R2 file not found: {r2_file}")
            return 0.0
            
        print(f"Processing {os.path.basename(r1_file)} and {os.path.basename(r2_file)}...")
        
        # Stream through files - equivalent to C++ approach
        r1_generator = read_fastq_records(r1_file)
        r2_generator = read_fastq_records(r2_file)
        
        # Process read pairs one at a time (streaming)
        while True:
            try:
                r1_seq = next(r1_generator)
                r2_seq = next(r2_generator)
            except StopIteration:
                # One or both files have ended - stop processing
                break
            
            total_reads += 1
            
            # Progress reporting for large files
            if total_reads % 100000 == 0:
                print(f"Processed {total_reads:,} read pairs...")
            
            try:
                # Find primer sequences
                r1_start = r1_seq.find("TACCGAGCTCGAATTTGCAC")
                if r1_start == -1:
                    continue
                r1_trim = r1_seq[r1_start + 20:]
                
                r2_start = r2_seq.find("GCGCCTGGATGAATTCGTGT")
                if r2_start == -1:
                    continue
                r2_trim = r2_seq[r2_start + 20:]
                
                # Align reads
                r1_elements = align_read(r1_trim, elements)
                r2_elements = align_read(r2_trim, elements)
                
                # Process R2 (reverse and negate)
                r2_elements.reverse()
                r2_elements = [-x for x in r2_elements]
                
                # Create consensus
                consensus_elements = consensus(r1_elements, r2_elements)
                
                # Validate consensus
                if (len(consensus_elements) > 0 and 
                    consensus_elements.count(0) == 0 and 
                    len(set(consensus_elements)) == len(consensus_elements)):
                    
                    barcode_reads += 1
                    consensus_tuple = tuple(consensus_elements)
                    
                    if consensus_tuple not in lox_counts:
                        lox_counts[consensus_tuple] = {'first_seen': barcode_reads, 'count': 1}
                    else:
                        lox_counts[consensus_tuple]['count'] += 1
                        
            except Exception as e:
                print(f"Error processing read pair {total_reads}: {e}")
                continue
    
    except Exception as e:
        print(f"Error processing files {r1_file}, {r2_file}: {e}")
        import traceback
        traceback.print_exc()
        return 0.0
    
    # Write output
    try:
        with open(output_file, 'w') as f:
            for barcode, data in lox_counts.items():
                f.write(f"{data['count']},")
                f.write(" ".join(map(str, barcode)))
                f.write("\n")
    except Exception as e:
        print(f"Error writing output file {output_file}: {e}")
        return 0.0
    
    percentage = 100 * barcode_reads / total_reads if total_reads > 0 else 0
    print(f"{os.path.basename(output_file)}: {percentage:.2f}% ({barcode_reads:,}/{total_reads:,} reads)")
    
    return percentage


# --- end unchanged functions ---

# LoxCode elements dictionary (moved global)
ELEMENTS = {
    "ACTCCGCA": 1, "TCCAGAATTTGTAT": 2, "ACATCCAC": 3, "AAAGGAATTTCTCC": 4,
    "ATTTCCTC": 5, "GCCCGAATTTTTTC": 6, "GCTACTGG": 7, "ATGAGAATTTATGG": 8,
    "AACTAGAA": 9, "TGCAGAATTTCCTC": 10, "CGACACTT": 11, "AACGGAATTTTCAA": 12,
    "CGTGTTTG": 13, "ACACGAATTCATCC": 14, "GTGCAAATTCGAGA": -15, "CAAACACG": -13,
    "TTGAAAATTCCGTT": -12, "AAGTGTCG": -11, "GAGGAAATTCTGCA": -10, "TTCTAGTT": -9,
    "CCATAAATTCTCAT": -8, "CCAGTAGC": -7, "GAAAAAATTCGGGC": -6, "GAGGAAAT": -5,
    "GGAGAAATTCCTTT": -4, "GTGGATGT": -3, "ATACAAATTCTGGA": -2, "TGCGGAGT": -1
}

def process_all_pairs(input_folder, r1_pattern, r2_pattern, output_folder):
    """Find and process all R1/R2 file pairs in a folder."""
    r1_patterns = [
        os.path.join(input_folder, f"*{r1_pattern}"),
        os.path.join(input_folder, f"*{r1_pattern}.gz")
    ]
    r2_patterns = [
        os.path.join(input_folder, f"*{r2_pattern}"),
        os.path.join(input_folder, f"*{r2_pattern}.gz")
    ]

    r1_files, r2_files = [], []
    for pattern in r1_patterns:
        r1_files.extend(glob.glob(pattern))
    for pattern in r2_patterns:
        r2_files.extend(glob.glob(pattern))

    r1_files = sorted(r1_files)
    r2_files = sorted(r2_files)

    if len(r1_files) == 0 and len(r2_files) == 0:
        raise FileNotFoundError(f"No FASTQ files found in {input_folder}")

    if len(r1_files) != len(r2_files):
        raise ValueError("Number of R1 and R2 files is unequal")

    for r1_file, r2_file in zip(r1_files, r2_files):
        base_name = os.path.basename(r1_file)
        if base_name.endswith('.gz'):
            base_name = base_name[:-3]  # strip .gz
        base_name = base_name.replace("R1_001.fastq", "")
        output_file = os.path.join(output_folder, f"{base_name}loxcode_counts.csv")

        process_file_pair(r1_file, r2_file, output_file, ELEMENTS)

