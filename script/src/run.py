import numpy as np
import itertools
import os
import subprocess
import re
from concurrent.futures import ThreadPoolExecutor

# Configuration
conf_in_path = '../../input/ic/0.py'
conf_out_dir = '../../input/ic_temp'
build_dir = '../../build'
exec_dir = '../../src'
# 異なる初期条件の並列スレッド数
n_thread = 1
# １回の計算の並列コア数
n_proc = 4

# Regular expression pattern to identify lines that need to be evaluated
pattern = re.compile(r'(\w+)\s*:\s*(\d+\.?\d*e?-?\d*)\s*->\s*(\d+\.?\d*e?-?\d*),\s*(n|na)=(\d+)')

# Read the input file
with open(conf_in_path, 'r') as f:
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
os.makedirs(conf_out_dir, exist_ok=True)

# Create a new file for each combination
output_filenames = []
for i, combination in enumerate(combinations):
    output_filename = os.path.join(conf_out_dir, f'{i:06d}.txt')
    output_filenames.append(output_filename)
    with open(output_filename, 'w') as f:
        comb_dict = dict(zip(keys, combination))
        for key, line in all_lines:
            if key is None:
                f.write(f'{line}\n')
            else:
                f.write(f'{key} = {comb_dict[key]:.6e}\n')

# Execute the Make command
subprocess.run(['make'], cwd=build_dir)

# Function to execute a command
def run_command(filename):
    command = ['mpirun', '-np', str(n_proc), '../build/bin/cff.out', f'../input/ic_temp/{filename}']
    subprocess.run(command, cwd=exec_dir)

# Execute in parallel
with ThreadPoolExecutor(max_workers=n_thread) as executor:
    executor.map(run_command, output_filenames)
