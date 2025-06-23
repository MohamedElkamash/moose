import csv
import os

# Input CSV file (no header)
input_file = '../particle_position.csv'  # replace with your filename
output_dir = 'couette'

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)
j = 0

with open(input_file, 'r') as infile:
    reader = csv.reader(infile)
    
    for i, row in enumerate(reader):
        if i > 2019:
            break


        # if len(row) != 3:
        #     print(f"Skipping invalid row {i}: {row}")
        #     continue
        
        if i%10 == 0:
            t, x, y, z, vx, vy, vz = row
            filename = f"c_{int(i/10)}.csv"  # Ensure clean float formatting
            filepath = os.path.join(output_dir, filename)

            with open(filepath, 'w', newline='') as outfile:
                writer = csv.writer(outfile, delimiter=',')
                writer.writerow(['x', 'y', 'z'])
                writer.writerow([x, y, z])
print(f"Generated files in: {output_dir}")