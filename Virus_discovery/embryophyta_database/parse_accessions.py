import json
import csv

def extract_accessions(json_path, tsv_path):
    """
    Parses a JSON catalog file to extract accession numbers and writes them
    to a TSV file.

    Args:
        json_path (str): The path to the input JSON file.
        tsv_path (str): The path for the output TSV file.
    """
    accession_numbers = []

    try:
        # Step 1: Read and parse the JSON file
        with open(json_path, 'r') as f:
            data = json.load(f)

        # Step 2: Extract accession numbers from the 'assemblies' list
        # We iterate through each item and get the value of the 'accession' key if it exists.
        if 'assemblies' in data and isinstance(data['assemblies'], list):
            for assembly in data['assemblies']:
                if 'accession' in assembly:
                    accession_numbers.append(assembly['accession'])
        
        # Step 3: Write the extracted numbers to a TSV file
        with open(tsv_path, 'w', newline='') as tsv_file:
            writer = csv.writer(tsv_file, delimiter='\t')
            # Write a header row for clarity
            writer.writerow(['accession_number'])
            # Write each accession number as a new row
            for acc in accession_numbers:
                writer.writerow([acc])
        
        print(f"Successfully extracted {len(accession_numbers)} accession numbers.")
        print(f"Output saved to: {tsv_path}")

    except FileNotFoundError:
        print(f"Error: The file '{json_path}' was not found.")
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from '{json_path}'. Please check the file format.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == '__main__':
    # Define the input and output file names
    json_file = 'dataset_catalog.json'
    tsv_file = 'accession_numbers.tsv'
    
    # Run the extraction function
    extract_accessions(json_file, tsv_file)