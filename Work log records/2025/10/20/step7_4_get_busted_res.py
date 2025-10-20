#!/usr/bin/env python3
import json
import os
import glob
import argparse
from typing import Dict, Any

def extract_busted_results(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract key results from a BUSTED JSON.
    """
    results = {}
    if "test results" not in data:
        return results  # Not BUSTED
    
    results["p_value"] = data["test results"].get("p-value", None)
    results["lrt"] = data["test results"].get("LRT", None)
    
    # Background test
    bg_test = data.get("test results background", {})
    results["bg_p_value"] = bg_test.get("p-value", None)
    results["bg_lrt"] = bg_test.get("LRT", None)
    
    # Shared distributions test
    shared_test = data.get("test results shared distributions", {})
    results["shared_p_value"] = shared_test.get("p-value", None)
    results["shared_lrt"] = shared_test.get("LRT", None)
    
    # Omega distributions from Unconstrained model
    unconstrained = data.get("fits", {}).get("Unconstrained model", {})
    rate_dist = unconstrained.get("Rate Distributions", {})
    
    # Test branches omega
    test_omega = rate_dist.get("Test", {})
    results["test_omega_2"] = test_omega.get("2", {}).get("omega", None)
    results["test_prop_2"] = test_omega.get("2", {}).get("proportion", None)
    
    # Background omega
    bg_omega = rate_dist.get("Background", {})
    results["bg_omega_2"] = bg_omega.get("2", {}).get("omega", None)
    results["bg_prop_2"] = bg_omega.get("2", {}).get("proportion", None)
    
    return results

def main():
    parser = argparse.ArgumentParser(description="Extract BUSTED results from JSON files in a folder.")
    parser.add_argument("--input-folder", required=True, help="Path to the input folder containing JSON files.")
    parser.add_argument("--output", required=True, help="Path to the output TXT file.")
    
    args = parser.parse_args()
    
    summary = []
    
    # Find all JSON files
    json_files = glob.glob(os.path.join(args.input_folder, "*.json"))
    
    for file_path in json_files:
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
            
            # Get gene name from filename prefix (without .json)
            filename = os.path.basename(file_path)
            gene_name = os.path.splitext(filename)[0]
            
            # Extract BUSTED results
            busted_results = extract_busted_results(data)
            
            if busted_results:
                row = [gene_name]
                row.extend([str(busted_results.get(key, "NA")) for key in ["p_value", "lrt", "bg_p_value", "bg_lrt", "shared_p_value", "shared_lrt", "test_omega_2", "test_prop_2", "bg_omega_2", "bg_prop_2"]])
                summary.append("\t".join(row))
            else:
                print(f"Warning: Not a BUSTED JSON for {filename}")
                
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
    
    # Write header
    header = ["Gene", "p_value", "LRT", "bg_p_value", "bg_LRT", "shared_p_value", "shared_LRT", "test_omega_2", "test_prop_2", "bg_omega_2", "bg_prop_2"]
    
    with open(args.output, 'w') as f:
        f.write("\t".join(header) + "\n")
        for row in summary:
            f.write(row + "\n")
    
    print(f"Extraction complete. Results saved to {args.output}")
    print(f"Processed {len(summary)} BUSTED analyses.")

if __name__ == "__main__":
    main()
