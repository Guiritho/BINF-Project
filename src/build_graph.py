import sys


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
    pass


def compute_connectivity(alpha):
    """Palier 2, Étape 2 : Détermination du seuil mou β Plutôt que de biner la matrice de
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
    pass


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
    pass


def find_optimal_beta(adjacency, gene_count):
    """Palier 2, Étape 2, point 7 : Répéter ces étapes pour toutes les valeurs entières
    0 < β <= 30 et retourner la plus petite valeur de β tel que R²_{β} > 0.8.
    """
    pass


def compute_link_overlap(alpha_tilde):
    """Palier 2, Étape 3 : Détermination de la matrice de connectivité Finalement, la matrice
    de connectivité finale du graphe est donnée par la matrice de topologie recouvrante
    (Topological Overlap Matrix, TOM) T définie par :

        t_{i,j} = |α_{i,j} + Σ_{u!=i,j} α^{~}_{i,u} x α^{~}_{u,j}|
                  / (min(c_{i}, c_{j}) + 1 - |α_{i,j}|),

    avec |.| la valeur absolue, α^{~}_{i,j} = α_{i,j} x sign(corr_{i,j}) et sign la fonction
    donnant le signe.
    """
    pass


def compute_tom(adjacency, correlation, beta):
    """Palier 2, Étape 3 : Détermination de la matrice de connectivité Finalement, la matrice
    de connectivité finale du graphe est donnée par la matrice de topologie recouvrante
    (Topological Overlap Matrix, TOM) T définie par :

        t_{i,j} = |α_{i,j} + Σ_{u!=i,j} α^{~}_{i,u} x α^{~}_{u,j}|
                  / (min(c_{i}, c_{j}) + 1 - |α_{i,j}|),

    avec |.| la valeur absolue, α^{~}_{i,j} = α_{i,j} x sign(corr_{i,j}) et sign la fonction
    donnant le signe.
    """
    pass


def run_pow(datadir, dataname):
    """Palier 2 : POW : le programme ne prend pas d'autre paramètre et affiche sur la sortie
    standard valeur opitmale de β calculée (cf ci-dessous).
    """
    pass


def run_tom(datadir, dataname, beta):
    """Palier 2 : TOM : le programme utilise aussi le paramètre beta et calcule la matrice
    TOM basée sur celui-ci (cf. ci-dessous).
    """
    pass


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
    has_beta_argument = len(sys.argv) > 4

    if option == "TOM" and has_beta_argument:
        beta = int(sys.argv[4])
        run_tom(datadir, dataname, beta)
    else:
        run_pow(datadir, dataname)


if __name__ == "__main__":
    main()
