#!/usr/bin/env python3

import sys
import os
import subprocess
from typing import List, Tuple
import time
import glob

def read_file_list(path: str) -> List[str]:
    with open(path, 'r') as f:
        return [line.strip() for line in f if line.strip()]

def read_expected_results(path: str) -> List[int]:
    with open(path, 'r') as f:
        return [int(line.strip()) for line in f if line.strip()]

def cleanup_files():
    """Remove all output and pathway files."""
    for f in glob.glob("*Out") + glob.glob("*Pathway"):
        try:
            os.remove(f)
        except:
            pass

def run_single_test(executable: str, mol_file: str) -> Tuple[int, float]:
    try:
        cleanup_files()
        start_time = time.time()
        result = subprocess.run([executable, mol_file], 
                              capture_output=True, 
                              text=True,
                              timeout=3600)
        
        base_name = os.path.splitext(mol_file)[0]
        output_file = f"{base_name}Out"
        
        if not os.path.exists(output_file):
            raise RuntimeError(f"Output file not found: {output_file}")
            
        with open(output_file, 'r') as f:
            lines = f.readlines()
            ai_line = next((line for line in lines if "has assembly index:" in line), None)
            if not ai_line:
                raise ValueError("Assembly index not found in output")
            ai = int(ai_line.split(":")[-1].strip())
            
            time_line = next((line for line in lines if "time to completion:" in line), None)
            runtime = float(time_line.split(":")[-1].strip()) if time_line else time.time() - start_time

        cleanup_files()
        return ai, runtime

    except Exception as e:
        cleanup_files()
        return -1, -1

def main():
    if len(sys.argv) != 4:
        print("Usage: test_assembly_index.py <executable> <mol_file_list> <expected_results>")
        sys.exit(1)

    executable = sys.argv[1]
    mol_files = read_file_list(sys.argv[2])
    expected_results = read_expected_results(sys.argv[3])

    if len(mol_files) != len(expected_results):
        print("Error: Number of test files doesn't match number of expected results")
        sys.exit(1)

    print(f"\nRunning {len(mol_files)} tests...\n")
    print(f"{'Molecule':<40} {'Expected':<10} {'Got':<10} {'Time(s)':<10} {'Status'}")
    print("-" * 80)

    total_passed = total_failed = total_errors = 0

    for mol_file, expected in zip(mol_files, expected_results):
        ai, runtime = run_single_test(executable, mol_file)
        
        if ai == -1:
            status, total_errors = "ERROR", total_errors + 1
        elif ai == expected:
            status, total_passed = "PASS", total_passed + 1
        else:
            status, total_failed = "FAIL", total_failed + 1

        print(f"{os.path.basename(mol_file):<40} {expected:<10} {ai:<10} {runtime:< 10.2f} {status}")

    print(f"\nSummary:")
    print(f"Passed: {total_passed}")
    print(f"Failed: {total_failed}")
    print(f"Errors: {total_errors}")
    print(f"Total:  {len(mol_files)}")

    cleanup_files()
    sys.exit(1 if total_failed > 0 or total_errors > 0 else 0)

if __name__ == "__main__":
    main()