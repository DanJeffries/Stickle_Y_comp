import re
import sys

def parse_file(filename):
    results = {}
    record_count = 0
    record_data = ""

    with open(filename) as file:
        for line in file:
            line = line.strip()
            
            if re.match(r"^\d", line):
                record_count += 1
                line_parts = line.split()
                div=line_parts[1]
                chromosome = line_parts[4]
                start_position = line_parts[5]
                end_position = line_parts[6]
                element = line_parts[9] if line_parts[8] == "C" else line_parts[8]
                
            elif re.match(r"^Kimura", line):
                line_parts = line.split()
                k2p_value = line_parts[4]
                record_data = f"{chromosome}\t{start_position}\t{end_position}\t{element}\t{k2p_value}\t{div}"
                if record_count > 0:
                    results[record_count] = record_data
            else:
                continue

    return results

def display_results(results):
    print("No.\tChr\tStart\tEnd\tFamily\tK2P\tDIV")
    for key, value in results.items():
        print(f"{key}\t{value}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <input_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    parsed_results = parse_file(input_file)
    display_results(parsed_results)

