#!/usr/bin/env python3
"""
Example usage of the eqdsk2vmec_converter package.
"""

from source import convert_eqdsk_to_vmec

def main():
    # Example usage
    filename = 'g02895743.geq'  # Replace with your GEQDSK file
    
    print(f"Converting {filename} to VMEC input format...")
    
    try:
        convert_eqdsk_to_vmec(filename)
        print("Conversion completed successfully!")
    except FileNotFoundError:
        print(f"Error: File {filename} not found.")
        print("Please provide a valid GEQDSK file.")
    except Exception as e:
        print(f"Error during conversion: {e}")

if __name__ == "__main__":
    main()

