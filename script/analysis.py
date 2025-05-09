import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def read_data(file_path):
    """Read CSV data from the given file path using pandas and extract time."""
    try:
        with open(file_path, "r") as file:
            file.readline()  # skip 1st line
            file_metadata = file.readline().strip().split('\t')
            time = float(file_metadata[1])  # Extract time from the second element

        # Read the rest of the data, skipping the first three rows
        data = pd.read_csv(file_path, skiprows=2, sep='\t')
        return data, time

    except Exception as e:
        print(f"Error in reading {file_path}: {e}")
        return pd.DataFrame(), 0.0
    
def filter_files(data_dir, max_files=None):
    # List all CSV files in the data folder, sorted by the integer suffix in their filenames
    data_files = sorted(
        [os.path.join(data_dir, f) for f in os.listdir(data_dir) if f.endswith(".tsv")],
        key=lambda x: int(os.path.splitext(os.path.basename(x))[0].split("_")[-1]),
    )

    if not data_files:
        print("No files found in the specified directory.")
        exit(-1)

    # If limiting how many files to load, select a subset
    if max_files and len(data_files) > max_files:
        step = max(1, len(data_files) // max_files)
        selected_files = data_files[::step][:max_files]
    else:
        selected_files = data_files
        
    return selected_files

def collect_tsv_data(base_dir):
    combined_data = pd.DataFrame()
    
    # 任意のディレクトリからTSVファイルを探す
    for dir_name in os.listdir(base_dir):  # base_dir内の全ディレクトリを探索
        if not dir_name.startswith('~~'):
            dir_path = os.path.join(base_dir, dir_name)
            if os.path.isdir(dir_path):
                tsv_file = os.path.join(dir_path, 'result.tsv')
                if os.path.exists(tsv_file):
                    # TSVファイルを読み込む（最初の1行＋ヘッダのみ）
                    df = pd.read_csv(tsv_file, sep='\t', nrows=1)
                    combined_data = pd.concat([combined_data, df], ignore_index=True)

    return combined_data
