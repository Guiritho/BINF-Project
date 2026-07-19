import sys
import numpy as np
from csv_utils import build_output_path, read_group_assignments, read_matrix_csv


def select_group_gene_indices(group_assignments, group_id):
    """Palier 4, Principe : Pour avoir une représentation génétique de chaque groupe, on va
    calculer son eigengene. Pour un groupe G : 1. Extraire l'ensemble K des positions des gènes
    présents dans G.
    """
    indices = [index for index, assigned_id in enumerate(group_assignments) if assigned_id == group_id]

    return indices


def normalize_group_matrix(data_subset):
    """Palier 4, Principe, point 3 : Centrer et normaliser D(G) : D(G) = (D(G) - E[D(G)]) /
    Var[D(G)], où E et Var sont respectivement la moyenne et variance calculées sur tous les
    éléments de la matrice.
    """
    # E[D(G)]
    mean_value = np.mean(data_subset)

    # Var[D(G)]
    variance_value = np.var(data_subset)

    # D(G) = (D(G) - E[D(G)]) / Var[D(G)]
    normalized = (data_subset - mean_value) / variance_value

    return normalized


def compute_eigengene(data_subset):
    """Palier 4, Principe, points 4-5 : Calculer la décomposition en valeurs singulières
    (SVD) de D(G) : D(G) = U x X x V^T, où U, X et V sont trois matrices. Vous pouvez utiliser la
    fonction svd de numpy.linalg. L'eigengene du groupe G correspond à la première colonne de la
    matrice V (Attention : on parle ici de V et pas V^T).
    """
    normalized = normalize_group_matrix(data_subset)

    # 4
    _, _, transposed_right_singular = np.linalg.svd(normalized)

    # 5
    eigengene_component_index = 0
    raw_eigengene = transposed_right_singular[eigengene_component_index, :]
    eigengene_round_decimals = 3
    eigengene = np.round(raw_eigengene, eigengene_round_decimals)

    return eigengene


def write_eigen_csv(filepath, eigengenes_by_group):
    """Palier 4, Format de sortie : Un fichier csv nommé [grpname]_eigen.csv où chaque ligne
    correspond à un groupe, trié par ordre d'identifiants, et les colonnes sont les valeurs du
    vecteur d'eigengene (arrondies à 3 décimales) séparées par des virgules. Pour chaque ligne,
    la première valeur est l'id du groupe. Notez que le nombre de valeurs va varier entre chaque
    groupe.
    """
    eigengene_round_decimals = 3
    value_format = f"%.{eigengene_round_decimals}f"

    with open(filepath, "w", newline="") as text_file:
        for group_id, eigengene in eigengenes_by_group:
            formatted_values = [value_format % value for value in eigengene]
            group_id_text = str(group_id)
            output_row = [group_id_text] + formatted_values
            row_text = ",".join(output_row)
            line = f"{row_text}\n"
            text_file.write(line)


def run_extract_eigen(datadir, dataname, groupname):
    """Palier 4, Paramètres d'entrée : datadir dataname groupname, où datadir est le chemin
    vers le dossier contenant les données, dataname le nom du fichier csv contenant les données
    d'expression génétique (du palier 1) et groupname le nom du fichier csv contenant les
    informations de groupes pour chaque gène (au format du palier précédent). Principe, point 2 :
    Récupérer la sous-matrice d'expression de ces gènes D(G) = D_{*,K}.
    """
    data_path = build_output_path(datadir, dataname)
    _, _, data = read_matrix_csv(data_path)

    group_path = build_output_path(datadir, groupname)
    group_assignments = read_group_assignments(group_path)

    distinct_groups = set(group_assignments)
    sorted_groups = sorted(distinct_groups)

    eigengenes_by_group = []
    for group_id in sorted_groups:
        indices = select_group_gene_indices(group_assignments, group_id)
        data_subset = data[:, indices]
        eigengene = compute_eigengene(data_subset)
        eigengenes_by_group.append((group_id, eigengene))

    output_filename = f"{groupname}_eigen.csv"
    output_path = build_output_path(datadir, output_filename)

    write_eigen_csv(output_path, eigengenes_by_group)


def main():
    """Palier 4, Paramètres d'entrée : datadir dataname groupname"""
    datadir = sys.argv[1]
    dataname = sys.argv[2]
    groupname = sys.argv[3]

    run_extract_eigen(datadir, dataname, groupname)


if __name__ == "__main__":
    main()
