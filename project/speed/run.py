import numpy as np
import itertools
import os
import subprocess
import re
from concurrent.futures import ThreadPoolExecutor

# Configuration
root_dir = '../..'
input_path = './speed.py'
input_temp_dir = './temp'
exec_name = "speed.out"
build_dir = f'{root_dir}/build'
exec_path = f'{root_dir}/build/bin/{exec_name}'
run_dir = '.'
# 異なる初期条件の並列スレッド数
n_thread = 1
# １回の計算の並列コア数
n_proc = 1

# Regular expression pattern to identify lines that need to be evaluated
pattern = re.compile(r'(\w+)\s*:\s*(\d+\.?\d*e?-?\d*)\s*->\s*(\d+\.?\d*e?-?\d*),\s*(n|na)=(\d+)')

# Read the input file
with open(input_path, 'r') as f:
    lines = f.readlines()

# Parse information for each variable
variables = {}
other_lines = []
all_lines = []

for line in lines:
    match = pattern.match(line.strip())
    if match:
        key, range_start, range_end, seq_type, n = match.groups()
        range_start, range_end = float(range_start), float(range_end)
        n = int(n)
        if seq_type == 'n':
            values = np.geomspace(range_start, range_end, n)
        elif seq_type == 'na':
            values = np.linspace(range_start, range_end, n)
        variables[key] = values
        all_lines.append((key, None))  # Placeholder for the variable
    else:
        other_lines.append(line.strip())
        all_lines.append((None, line.strip()))  # Preserve non-variable lines

# Generate all combinations
keys = list(variables.keys())
combinations = list(itertools.product(*variables.values()))

# Create output directory if it does not exist
os.makedirs(input_temp_dir, exist_ok=True)

# Create a new file for each combination
output_filenames = []
for i, combination in enumerate(combinations):
    output_filename = os.path.join(input_temp_dir, f'{i:06d}.txt')
    output_filenames.append(output_filename)
    with open(output_filename, 'w') as f:
        comb_dict = dict(zip(keys, combination))
        for key, line in all_lines:
            if key is None:
                f.write(f'{line}\n')
            else:
                f.write(f'{key} = {comb_dict[key]:.6e}\n')

# Execute the Make command
result = subprocess.run(['make', exec_name], cwd=build_dir)

# Check if the compilation was successful
if result.returncode != 0:
    print("Compilation failed, skipping execution.")
else:
    # Function to execute a command
    def run_command(filename):
        command = ['mpirun', '-np', str(n_proc), exec_path, input_path]
        subprocess.run(command, cwd=run_dir)

    # Execute in parallel
    with ThreadPoolExecutor(max_workers=n_thread) as executor:
        executor.map(run_command, output_filenames)
