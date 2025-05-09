import numpy as np
import itertools
import os
import subprocess
import re
from multiprocessing import Process
from datetime import datetime
import sys

# Configuration
ws_root_path = '/home/shusuke/workspace/elasticity'  # This must be an absolute path
input_path = 'debug.py'
exec_name = 'debug.out'
build_path = f'{ws_root_path}/build'
exec_path = f'{ws_root_path}/build/bin/{exec_name}'

run_dir = 'runs'  # Base directory for individual runs

# Number of threads and processes
n_thread = 10  # Number of parallel runs
n_proc = 1     # Number of MPI processes per run

args = sys.argv
if 1 <= len(args) and args[0].isdigit(): n_proc = args[0]

# Regular expression pattern to identify lines that need to be evaluated
pattern = re.compile(r'(\w+)\s*:\s*(\d+\.?\d*e?-?\d*)\s*->\s*(\d+\.?\d*e?-?\d*),\s*(n|na)=(\d+)')

# Execute the Make command
print('make')
result = subprocess.run(['make', exec_name], cwd=build_path)

# Check if the compilation was successful
if result.returncode != 0:
    print("Compilation failed.")
    exit(-1)

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

# Create the run directory if it does not exist
os.makedirs(run_dir, exist_ok=True)

# Get the current date and time for directory names
timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

# Create a new file and unique directory for each combination
output_filenames = []
run_directories = []
for i, combination in enumerate(combinations):
    run_subdir = os.path.join(run_dir, f'{timestamp}-{i:04d}')
    os.makedirs(run_subdir, exist_ok=True)
    os.makedirs(os.path.join(run_subdir, 'dump'), exist_ok=True)
    os.makedirs(os.path.join(run_subdir, 'debug', 'distr'), exist_ok=True)
    run_directories.append(run_subdir)

    output_filename = os.path.join(run_subdir, 'input.txt')
    output_filenames.append(output_filename)
    
    with open(output_filename, 'w') as f:
        comb_dict = dict(zip(keys, combination))
        for key, line in all_lines:
            if key is None:
                f.write(f'{line}\n')
            else:
                f.write(f'{key} = {comb_dict[key]:.6e}\n')

from concurrent.futures import ThreadPoolExecutor

def run_command_in_background(run_dir):
    # Temporary directory with '~' prefix
    temp_dir = f"{os.path.dirname(run_dir)}/~{os.path.basename(run_dir)}"
    os.rename(run_dir, temp_dir)

    stdout_file = os.path.join(temp_dir, 'stdout.log')
    stderr_file = os.path.join(temp_dir, 'stderr.log')
    command = ['mpirun', '-np', str(n_proc), exec_path, 'input.txt']

    try:
        # Start the process in the background using Popen with session detachment
        with open(stdout_file, 'w') as stdout, open(stderr_file, 'w') as stderr:
            process = subprocess.Popen(
                command,
                cwd=temp_dir,
                stdout=stdout,
                stderr=stderr,
                # stdin='/dev/null',
                start_new_session=True  # Ensures detachment from the terminal
            )

        # Monitor the process in a separate thread
        def monitor_process():
            try:
                process.wait()  # Wait for the process to complete
                if process.returncode != 0:
                    raise subprocess.CalledProcessError(process.returncode, command)
            except Exception:
                # Rename directory to indicate failure
                error_dir = f"{os.path.dirname(temp_dir)}/!{os.path.basename(temp_dir)[1:]}"
                os.rename(temp_dir, error_dir)
            else:
                # Restore original directory name after success
                os.rename(temp_dir, run_dir)

        # Submit the monitoring task to a thread
        with ThreadPoolExecutor(max_workers=1) as executor:
            executor.submit(monitor_process)

    except Exception:
        pass  # Fail silently to avoid terminal output

# Execute in parallel processes without blocking terminal
with ThreadPoolExecutor(max_workers=n_thread) as executor:
    for run_dir in run_directories:
        executor.submit(run_command_in_background, run_dir)

# All processes submitted. Terminal remains available.
print("All processes started in background.")
