import json
import sys

def modify_json_file(file_path):
    """
    Loads a JSON file, adds a 'sample_key' to nested dictionaries,
    and saves it back to the same file.
    """
    try:
        # Open the file for reading
        with open(file_path, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        return
    except json.JSONDecodeError:
        print(f"Error: The file '{file_path}' does not contain valid JSON.")
        return

    # Iterate over the outer dictionary (e.g., "blocco1", "blocco2")
    for blocco_key, samples in data.items():
        if isinstance(samples, dict):
            # Iterate over the inner dictionary (e.g., "c26STAT3", "sham")
            for sample_name, details in samples.items():
                if isinstance(details, dict):
                    # Add the 'sample_key' field
                    details['sample_key'] = f"{blocco_key}_{sample_name}"

    try:
        # Open the same file for writing, which will overwrite it
        with open(file_path, 'w') as f:
            json.dump(data, f, indent=4)
        print(f"Successfully modified and saved '{file_path}'")
    except IOError as e:
        print(f"Error: Could not write to file '{file_path}'. Reason: {e}")

if __name__ == "__main__":
    # Ensure a file path is provided as a command-line argument
    if len(sys.argv) != 2:
        print("Usage: python modify_json.py <path_to_your_json_file>")
    else:
        file_path = sys.argv[1]
        modify_json_file(file_path)
