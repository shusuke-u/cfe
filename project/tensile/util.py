import os
import pandas as pd

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