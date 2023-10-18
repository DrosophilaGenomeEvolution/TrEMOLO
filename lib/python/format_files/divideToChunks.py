import pandas as pd
import numpy as np
import os
from multiprocessing import Pool, cpu_count, current_process

def save_chunk(chunk_data, folder):
    pid = current_process().pid
    idx = chunk_data['idx']
    filename = os.path.join(folder, f"tmp_chunk_{idx}.csv")
    chunk_data['data'].to_csv(filename, index=False)
    print(f"Chunk {idx} saved to {filename}")

def main(input_file, num_splits, output_folder):
    # Lecture du fichier csv
    df = pd.read_csv(input_file)

    # Créer le dossier de sortie s'il n'existe pas
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Divisez le DataFrame en sous-ensembles
    chunks = [(idx, sub_df) for idx, sub_df in enumerate(np.array_split(df, num_splits), 1)]

    # Utilisez le multiprocessing pour sauvegarder chaque morceau
    with Pool(processes=min(cpu_count(), num_splits)) as pool:
        pool.starmap(save_chunk, [({"idx": idx, "data": data}, output_folder) for idx, data in chunks])

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Divise un fichier CSV en plusieurs morceaux et les sauvegarde dans un dossier spécifié.")
    parser.add_argument('input_file', type=str, help='Chemin vers le fichier CSV d\'entrée.')
    parser.add_argument('num_splits', type=int, help='Nombre de morceaux à créer.')
    parser.add_argument('output_folder', type=str, help='Dossier où sauvegarder les morceaux de CSV.')
    
    args = parser.parse_args()
    
    main(args.input_file, args.num_splits, args.output_folder)
