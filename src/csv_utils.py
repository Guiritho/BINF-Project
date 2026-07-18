import csv
from pathlib import Path
import numpy as np


def read_matrix_csv(filepath):
    """Palier 1, Format d'entrée : Le ficher csv d'entrées est au format suivant : la première
    ligne contient les identifiants de chaque colonne (la première colonne de cette ligne est
    vide), ensuite pour chaque ligne, la première colonne contient son identifiant et les
    colonnes suivantes une valeur d'expression pour chaque gène. Les données d'expression
    génétiques sont structurées sous la forme d'une matrice D de taille M x N où M est le nombre
    d'observations (échantillons, cellules ou patients) et N le nombre de gènes mesurés.
    """
    with open(filepath, newline="") as csv_file:
        reader = csv.reader(csv_file)
        header_row = next(reader)
        gene_names = header_row[1:]

        row_ids = []
        data_rows = []
        for row in reader:
            row_ids.append(row[0])
            data_rows.append(row[1:])

    data = np.array(data_rows, dtype=float)

    return gene_names, row_ids, data


def write_matrix_csv(filepath, header, row_ids, data, decimals):
    """Palier 1, Format de sortie : Vous devez créer un fichier au format csv nommé
    [dataname]_cleaned_[K].csv (où [dataname] et [K] sont à remplacer par leurs valeur de
    paramètres d'entrées associés) situé dans le dossier donné par le paramètre datadir. Ce
    fichier a le même format que le fichier d'entrée mais avec les valeurs arrondies à 5
    décimales.
    """
    rounded_data = np.round(data, decimals)
    value_format = f"%.{decimals}f"

    with open(filepath, "w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        header_list = list(header)
        header_row = [""] + header_list
        writer.writerow(header_row)
        for row_id, row_values in zip(row_ids, rounded_data):
            formatted_values = [value_format % value for value in row_values]
            output_row = [row_id] + formatted_values
            writer.writerow(output_row)


def write_matrix_only_csv(filepath, data, decimals):
    """Palier 2, Format de sortie : Pour l'option TOM uniquement : un fichier au format csv
    nommé [dataname]_TOM_[beta].csv (où [dataname] et [beta] sont à remplacer par la valeur des
    paramètres d'entrées associés) situé dans le dossier donné par le paramètre datadir. Ce
    fichier contient la matrice de recouvrement (TOM, définie ci-dessous), les valeurs d'une
    même ligne sont séparées par des virgules et arrondies à 5 décimales. Note : à cette étape,
    on conserve uniquement les données de la matrice, pas les noms des lignes et colonnes.
    """
    rounded_data = np.round(data, decimals)
    value_format = f"%.{decimals}f"

    with open(filepath, "w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        for row_values in rounded_data:
            formatted_values = [value_format % value for value in row_values]
            writer.writerow(formatted_values)


def read_matrix_only_csv(filepath):
    """Palier 3, Paramètres d'entrée : dataname le nom du fichier csv contenant les données
    (au format du palier précédent).
    """
    data = np.loadtxt(filepath, delimiter=",")

    return data


def read_group_assignments(filepath):
    """Palier 3, Format de sortie : Ce fichier contient sur chaque ligne k l'id du groupe
    auquel appartient le kème gène de la matrice TOM. L'ordre doit être le même que celui des
    lignes et colonnes de la matrice d'entrée.
    """
    with open(filepath) as text_file:
        group_ids = [int(line.strip()) for line in text_file if line.strip()]

    return group_ids


def write_group_assignments(filepath, group_ids):
    """Palier 3, Format de sortie : Ce fichier contient sur chaque ligne k l'id du groupe
    auquel appartient le kème gène de la matrice TOM. L'ordre doit être le même que celui des
    lignes et colonnes de la matrice d'entrée.
    """
    with open(filepath, "w", newline="") as text_file:
        for group_id in group_ids:
            line = f"{group_id}\n"
            text_file.write(line)


def build_output_path(datadir, filename):
    """datadir est le chemin vers le dossier contenant les données. Les fichiers produits par
    les 4 programmes sont situés dans le dossier donné par le paramètre datadir.
    """
    output_path = Path(datadir) / filename

    return output_path
