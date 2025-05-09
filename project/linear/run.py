import numpy as np
import itertools
import os
import subprocess
import re
from multiprocessing import cpu_count
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
import sys
import time
import psutil

# Configuration
name = 'linear'
ws_root = '/home/shusuke/workspace/cfe'  # This must be an absolute path
input_path = f'{name}.py'
exec_name = f'{name}.out'
build_path = f'{ws_root}/build'
exec_path = f'{ws_root}/build/bin/{exec_name}'

run_dir = 'runs'  # Base directory for individual runs

queue_prefix = '*'
error_prefix = '!'
busy_prefix = '~'

# Number of threads and processes
n_cores = psutil.cpu_count(logical=False)  # Use physical cores
n_procs = 4  # Number of MPI processes per run

args = sys.argv
if len(args) > 1 and args[1].isdigit():
    n_procs = int(args[1])

# Ensure we do not exceed hardware limits
# max_parallel_jobs = max(1, n_cores // n_procs)  # Ensure at least 1 job is allowed
max_parallel_jobs = 1

print(f"CPU Cores Available: {n_cores}")
print(f"Max Parallel Jobs: {max_parallel_jobs}")
print(f"Running up to {max_parallel_jobs} parallel runs with {n_procs} MPI processes each.")

# Regular expression pattern to identify lines that need to be evaluated
pattern = re.compile(r'(\w+)\s*:\s*(\d+\.?\d*e?-?\d*)\s*->\s*(\d+\.?\d*e?-?\d*),\s*(n|na)=(\d+)')

# Execute the Make command
print('make')
result = subprocess.run(['make', exec_name], cwd=build_path, stdout=sys.stdout, stderr=sys.stderr)

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
        values = np.geomspace(range_start, range_end, n) if seq_type == 'n' else np.linspace(range_start, range_end, n)
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
    base_dir = os.path.join(run_dir, f'{timestamp}_{i:04d}')
    waiting_dir = os.path.join(run_dir, f'{queue_prefix}{timestamp}_{i:04d}')
    
    os.makedirs(waiting_dir, exist_ok=True)
    os.makedirs(os.path.join(waiting_dir, 'dump'), exist_ok=True)
    os.makedirs(os.path.join(waiting_dir, 'debug', 'distr'), exist_ok=True)
    run_directories.append(waiting_dir)  # Store as waiting dir initially

    output_filename = os.path.join(waiting_dir, 'input.txt')
    output_filenames.append(output_filename)
    
    with open(output_filename, 'w') as f:
        comb_dict = dict(zip(keys, combination))
        f.writelines([
            f'{key} = {comb_dict[key]:.6e}\n' if key else f'{line}\n' for key, line in all_lines
        ])

def run_command_in_background(waiting_dir):
    """Runs an MPI job in a separate process and ensures proper handling of failures."""
    base_name = waiting_dir.replace(queue_prefix, '', 1)  # Remove waiting prefix
    running_dir = os.path.join(run_dir, f'{busy_prefix}{os.path.basename(base_name)}')
    
    os.rename(waiting_dir, running_dir)  # Mark as running

    stdout_file = os.path.join(running_dir, 'stdout.log')
    stderr_file = os.path.join(running_dir, 'stderr.log')
    command = ['mpirun', '-np', str(n_procs), exec_path, 'input.txt']

    try:
        with open(stdout_file, 'w') as stdout, open(stderr_file, 'w') as stderr:
            process = subprocess.Popen(
                command,
                cwd=running_dir,
                stdout=stdout,
                stderr=stderr,
                start_new_session=True
            )

        process.wait()
        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, command)

        os.rename(running_dir, base_name)  # Restore directory name after success
    except Exception as e:
        error_dir = os.path.join(run_dir, f'{error_prefix}{os.path.basename(base_name)}')
        os.rename(running_dir, error_dir)
        print(f"Job failed in {running_dir}: {str(e)}")
        with open(stderr_file) as f:
            print(f"Error log:\n{f.read()}")

# Execute in parallel processes while respecting CPU limits
print("Starting parallel execution...")

running_jobs = []
with ThreadPoolExecutor(max_workers=max_parallel_jobs) as executor:
    for waiting_dir in run_directories:
        # Wait if we are running too many jobs
        while len(running_jobs) >= max_parallel_jobs:
            for future in as_completed(running_jobs):
                running_jobs.remove(future)
                try:
                    future.result()  # Retrieve any exceptions
                except Exception as e:
                    print(f"Error in execution: {str(e)}")
            time.sleep(1)  # Small delay to prevent excessive CPU polling

        # Submit a new job
        future = executor.submit(run_command_in_background, waiting_dir)
        running_jobs.append(future)

    # Wait for remaining jobs to finish
    for future in as_completed(running_jobs):
        try:
            future.result()
        except Exception as e:
            print(f"Error in execution: {str(e)}")

print("All processes completed.")