import sys


def select_group_gene_indices(group_assignments, group_id):
    """Palier 4, Principe : Pour avoir une représentation génétique de chaque groupe, on va
    calculer son eigengene. Pour un groupe G : 1. Extraire l'ensemble K des positions des gènes
    présents dans G.
    """
    pass


def normalize_group_matrix(data_subset):
    """Palier 4, Principe, point 3 : Centrer et normaliser D(G) : D(G) = (D(G) - E[D(G)]) /
    Var[D(G)], où E et Var sont respectivement la moyenne et variance calculées sur tous les
    éléments de la matrice.
    """
    pass


def compute_eigengene(data_subset):
    """Palier 4, Principe, points 4-5 : Calculer la décomposition en valeurs singulières
    (SVD) de D(G) : D(G) = U x X x V^T, où U, X et V sont trois matrices. Vous pouvez utiliser la
    fonction svd de numpy.linalg. L'eigengene du groupe G correspond à la première colonne de la
    matrice V (Attention : on parle ici de V et pas V^T).
    """
    pass


def write_eigen_csv(filepath, eigengenes_by_group):
    """Palier 4, Format de sortie : Un fichier csv nommé [grpname]_eigen.csv où chaque ligne
    correspond à un groupe, trié par ordre d'identifiants, et les colonnes sont les valeurs du
    vecteur d'eigengene (arrondies à 3 décimales) séparées par des virgules. Pour chaque ligne,
    la première valeur est l'id du groupe. Notez que le nombre de valeurs va varier entre chaque
    groupe.
    """
    pass


def run_extract_eigen(datadir, dataname, groupname):
    """Palier 4, Paramètres d'entrée : datadir dataname groupname, où datadir est le chemin
    vers le dossier contenant les données, dataname le nom du fichier csv contenant les données
    d'expression génétique (du palier 1) et groupname le nom du fichier csv contenant les
    informations de groupes pour chaque gène (au format du palier précédent). Principe, point 2 :
    Récupérer la sous-matrice d'expression de ces gènes D(G) = D_{*,K}.
    """
    pass


def main():
    """Palier 4, Paramètres d'entrée : datadir dataname groupname"""
    datadir = sys.argv[1]
    dataname = sys.argv[2]
    groupname = sys.argv[3]

    run_extract_eigen(datadir, dataname, groupname)


if __name__ == "__main__":
    main()
