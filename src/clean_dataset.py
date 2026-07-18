import sys
import numpy as np
from csv_utils import build_output_path, read_matrix_csv, write_matrix_csv


def select_top_variance_genes(data, gene_count):
    """Palier 1, Principe : A ce palier, on veut enlever les gènes (colonnes de D) qui ne
    fluctuent pas beaucoup entre observations. Pour ce faire : calculer la variance de chaque
    gène et garder uniquement les N gènes avec la plus forte variance.
    """
    variances = np.var(data, axis=0)
    highest_variance_indices = np.argsort(variances)[::-1][:gene_count]

    return highest_variance_indices


def clean_dataset(datadir, dataname, gene_count):
    """Palier 1, Paramètres d'entrée : datadir dataname K, où datadir est le chemin vers le
    dossier contenant les données, dataname le nom du fichier csv contenant les données brutes
    et K est le nombre de gènes a conserver.
    """
    input_path = build_output_path(datadir, dataname)
    header, row_ids, data = read_matrix_csv(input_path)

    selected_indices = select_top_variance_genes(data, gene_count)
    cleaned_header = [header[index] for index in selected_indices]
    cleaned_data = data[:, selected_indices]

    output_filename = f"{dataname}_cleaned_{gene_count}.csv"
    output_path = build_output_path(datadir, output_filename)
    output_decimals = 5

    write_matrix_csv(output_path, cleaned_header, row_ids, cleaned_data, output_decimals)


def main():
    """Palier 1, Paramètres d'entrée : datadir dataname K"""
    datadir = sys.argv[1]
    dataname = sys.argv[2]
    gene_count = int(sys.argv[3])

    clean_dataset(datadir, dataname, gene_count)


if __name__ == "__main__":
    main()
