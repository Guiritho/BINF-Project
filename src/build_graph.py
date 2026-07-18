import sys
import numpy as np
from scipy.stats import linregress
from csv_utils import build_output_path, read_matrix_csv, write_matrix_only_csv


def compute_correlation_adjacency(data):
    """Palier 2, Étape 1 : Matrice de corrélation Pour construire le réseau d'expression
    génétique, on va calculer la corrélation entre l'expression de chaque gène dans les
    observations. Pour une matrice de données D et deux gènes d'indices i et j, la corrélation
    est donnée par :

        corr_{i,j} = CovD_{*,i},D_{*,j} / Var[D_{*,i}] x Var[D_{*,j}]

    avec Cov la covariance, Var la variance et le symbole * en indice indique que l'on
    considère toute la ligne / colonne. Ici vous pouvez utiliser la fonction corrcoef de numpy.
    On transforme ensuite cette corrélation dans [-1, 1] en une valeur dans [0, 1] que l'on
    conserve dans une matrice A de taille N x N où N est le nombre de gènes (colonnes) de D :

        a_{i,j} = corr_{i,j}+1 / 2.
    """
    correlation = np.corrcoef(data, rowvar=False)
    adjacency = (correlation + 1) / 2

    return correlation, adjacency


def compute_connectivity(alpha):
    """Palier 2, Étape 2, point 1 : Détermination du seuil mou β Plutôt que de biner la matrice de
    corrélation, on va utiliser un seuil "mou" β et définir une nouvelle matrice α :

        α_{i,j} = (a_{i,j})^β.

    Pour déterminer β, on fait l'hypothèse que notre réseau doit être "scale-free", une
    propriété observée dans la plupart des réseaux biologiques. Dans un tel réseau, la
    distribution des poids de connectivité k suit une loi de puissance :

        p(k) ~ k^-β

    où p représente la densité de probabilité et β est le paramètre de la distribution. On va
    alors chercher la valeur positive entière minimum de β tel que la densité de probabilité
    observée p^{~} soit la plus proche d'une loi de puissance. Pour ce faire, voici les étapes
    à suivre : 1. Calculer pour chaque gène sa connectivité c :

        c_{i} = Σ_{j=0}^{N-1} α_{i,j}
    """
    connectivity = alpha.sum(axis=1)

    return connectivity


def compute_power_law_r_squared(connectivity, gene_count):
    """Palier 2, Étape 2, points 2 à 6 : Générer la distribution empirique p^{~} via un
    histogramme, vous pouvez utiliser numpy.histogram. Récupérer les coordonnées b des paquets
    (bins) et les valeurs o associées. Utiliser ceil(sqrt(N)) paquets, où N est le nombre de
    gènes et ceil la fonction arrondissant à l'entier supérieur. Si besoin, convertissez les
    valeurs des paquets b_{i} pour qu'elles représentent les centres et non les bords de chaque
    paquet. Enlever toutes les paires (b_{i}, o_{i}), 0 <= i < N, pour lesquelles o_{i} = 0. Une
    loi de puissance implique une relation linéaire en axes log-log. On va donc approximer
    notre distribution par la droite :

        log10(o) = u x log10(b) + v

    en utilisant, par exemple, linregress de scipy.stats. Les coefficients (u, v) en eux même
    ne nous intéressent pas, ce qu'on veut c'est récupérer le coefficient de détermination
    associé à cette valeur de β : R²_{β}.
    """
    # 2
    gene_count_sqrt = np.sqrt(gene_count)
    gene_count_sqrt_ceiled = np.ceil(gene_count_sqrt)
    bin_count = int(gene_count_sqrt_ceiled)
    counts, edges = np.histogram(connectivity, bins=bin_count)

    # 3
    centers = (edges[:-1] + edges[1:]) / 2

    # 4
    non_empty_bins = counts > 0
    counts = counts[non_empty_bins]
    centers = centers[non_empty_bins]

    # 5
    log_centers = np.log10(centers)
    log_counts = np.log10(counts)
    regression = linregress(log_centers, log_counts)

    # 6
    r_squared = regression.rvalue ** 2

    return r_squared


def find_optimal_beta(adjacency, gene_count):
    """Palier 2, Étape 2, point 7 : Répéter ces étapes pour toutes les valeurs entières
    0 < β <= 30 et retourner la plus petite valeur de β tel que R²_{β} > 0.8.
    """
    min_beta = 1
    max_beta = 30
    beta_upper_bound = max_beta + 1
    for beta in range(min_beta, beta_upper_bound):
        # On récupère la connectivité (point 1)
        alpha = adjacency ** beta
        connectivity = compute_connectivity(alpha)

        # On récupère un nouveau beta (points 2 à 6)
        r_squared = compute_power_law_r_squared(connectivity, gene_count)
        r_squared_threshold = 0.8

        if r_squared > r_squared_threshold:
            return beta

    return max_beta


def compute_link_overlap(alpha_tilde):
    """Palier 2, Étape 3 : Détermination de la matrice de connectivité Finalement, la matrice
    de connectivité finale du graphe est donnée par la matrice de topologie recouvrante
    (Topological Overlap Matrix, TOM) T définie par :

        t_{i,j} = |α_{i,j} + Σ_{u!=i,j} α^{~}_{i,u} x α^{~}_{u,j}|
                  / (min(c_{i}, c_{j}) + 1 - |α_{i,j}|),

    avec |.| la valeur absolue, α^{~}_{i,j} = α_{i,j} x sign(corr_{i,j}) et sign la fonction
    donnant le signe.
    """
    # Somme de tous les u (u=i et u=j inclus)
    # Σ_{u=0}^{N-1} α^{~}_{i,u} x α^{~}_{u,j}
    full_sum = alpha_tilde @ alpha_tilde

    # Les deux termes à retirer
    # u=i (α^{~}_{i,i} x α^{~}_{i,j})
    # u=j (α^{~}_{i,j} x α^{~}_{j,j})
    diagonal = np.diag(alpha_tilde)
    self_terms = alpha_tilde * (diagonal[:, None] + diagonal[None, :])

    # Sur la diagonale (i=j), u=i et u=j sont le même terme
    # du coup ça veut dire que self_terms le compte deux fois
    # donc on soustrait α^{~}_{i,i}^2 une fois pour le retirer une fois
    squared_diagonal = diagonal ** 2
    diagonal_correction = np.diag(squared_diagonal)
    self_terms -= diagonal_correction

    # On retire tous les termes u=i et u=j
    # Σ_{u!=i,j} α^{~}_{i,u} x α^{~}_{u,j}
    link_overlap = full_sum - self_terms

    return link_overlap


def compute_tom(adjacency, correlation, beta):
    """Palier 2, Étape 3 : Détermination de la matrice de connectivité Finalement, la matrice
    de connectivité finale du graphe est donnée par la matrice de topologie recouvrante
    (Topological Overlap Matrix, TOM) T définie par :

        t_{i,j} = |α_{i,j} + Σ_{u!=i,j} α^{~}_{i,u} x α^{~}_{u,j}|
                  / (min(c_{i}, c_{j}) + 1 - |α_{i,j}|),

    avec |.| la valeur absolue, α^{~}_{i,j} = α_{i,j} x sign(corr_{i,j}) et sign la fonction
    donnant le signe.
    """
    # α_{i,j}
    alpha = adjacency ** beta

    # α^{~}_{i,j} = α_{i,j} x sign(corr_{i,j})
    alpha_tilde = alpha * np.sign(correlation)

    # c_{i} et c_{j}
    connectivity = compute_connectivity(alpha)

    # Σ_{u!=i,j} α^{~}_{i,u} x α^{~}_{u,j}
    link_overlap = compute_link_overlap(alpha_tilde)

    # min(c_{i}, c_{j}
    min_connectivity = np.minimum.outer(connectivity, connectivity)

    # α_{i,j} + Σ_{u!=i,j} α^{~}_{i,u} x α^{~}_{u,j}
    alpha_plus_overlap = alpha + link_overlap

    # |α_{i,j} + Σ_{u!=i,j} α^{~}_{i,u} x α^{~}_{u,j}|
    numerator = np.abs(alpha_plus_overlap)

    # |α_{i,j}|
    alpha_abs = np.abs(alpha)

    # min(c_{i}, c_{j}) + 1 - |α_{i,j}|
    denominator = min_connectivity + 1 - alpha_abs

    # t_{i,j} = |α_{i,j} + Σ_{u!=i,j} α^{~}_{i,u} x α^{~}_{u,j}|
    #           / (min(c_{i}, c_{j}) + 1 - |α_{i,j}|)
    tom = numerator / denominator

    return tom


def run_pow(datadir, dataname):
    """Palier 2 : POW : le programme ne prend pas d'autre paramètre et affiche sur la sortie
    standard valeur opitmale de β calculée (cf ci-dessous).
    """
    input_path = build_output_path(datadir, dataname)
    _, _, data = read_matrix_csv(input_path)
    gene_count = data.shape[1]

    _, adjacency = compute_correlation_adjacency(data)
    beta = find_optimal_beta(adjacency, gene_count)

    print(beta)


def run_tom(datadir, dataname, beta):
    """Palier 2 : TOM : le programme utilise aussi le paramètre beta et calcule la matrice
    TOM basée sur celui-ci (cf. ci-dessous).
    """
    input_path = build_output_path(datadir, dataname)
    _, _, data = read_matrix_csv(input_path)

    correlation, adjacency = compute_correlation_adjacency(data)
    tom = compute_tom(adjacency, correlation, beta)

    output_filename = f"{dataname}_TOM_{beta}.csv"
    output_path = build_output_path(datadir, output_filename)
    output_decimals = 5

    write_matrix_only_csv(output_path, tom, output_decimals)


def main():
    """Palier 2, Paramètres d'entrée : datadir dataname option [beta], où datadir est le chemin
    vers le dossier contenant les données, dataname le nom du fichier csv contenant les données
    (au format du palier précédent). Les options possibles sont : POW : le programme ne prend
    pas d'autre paramètre et affiche sur la sortie standard valeur opitmale de β calculée (cf
    ci-dessous). TOM : le programme utilise aussi le paramètre beta et calcule la matrice TOM
    basée sur celui-ci (cf. ci-dessous).
    """
    datadir = sys.argv[1]
    dataname = sys.argv[2]
    option = sys.argv[3]
    is_tom_option = option == "TOM"
    has_beta_argument = len(sys.argv) > 4

    if is_tom_option and has_beta_argument:
        beta = int(sys.argv[4])
        run_tom(datadir, dataname, beta)
    else:
        run_pow(datadir, dataname)


if __name__ == "__main__":
    main()
