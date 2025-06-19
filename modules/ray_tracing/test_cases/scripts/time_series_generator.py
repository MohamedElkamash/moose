import csv
import os

# Input CSV file (no header)
input_file = 'output_files/developed_1.csv'  # replace with your filename
output_dir = 'developed_1'

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

with open(input_file, 'r') as infile:
    reader = csv.reader(infile)
    
    for i, row in enumerate(reader):
        if i > 102:
            break


        if len(row) != 3:
            print(f"Skipping invalid row {i}: {row}")
            continue
        
        x, y, z = row
        filename = f"dev_fast_{int(i)}.csv"  # Ensure clean float formatting
        filepath = os.path.join(output_dir, filename)

        with open(filepath, 'w', newline='') as outfile:
            writer = csv.writer(outfile, delimiter=',')
            writer.writerow(['x', 'y', 'z'])
            writer.writerow([x, y, z])

print(f"Generated files in: {output_dir}")