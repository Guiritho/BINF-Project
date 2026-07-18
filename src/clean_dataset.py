import sys


def select_top_variance_genes(data, gene_count):
    """Palier 1, Principe : A ce palier, on veut enlever les gènes (colonnes de D) qui ne
    fluctuent pas beaucoup entre observations. Pour ce faire : calculer la variance de chaque
    gène et garder uniquement les N gènes avec la plus forte variance.
    """
    pass


def clean_dataset(datadir, dataname, gene_count):
    """Palier 1, Paramètres d'entrée : datadir dataname K, où datadir est le chemin vers le
    dossier contenant les données, dataname le nom du fichier csv contenant les données brutes
    et K est le nombre de gènes a conserver.
    """
    pass


def main():
    """Palier 1, Paramètres d'entrée : datadir dataname K"""
    datadir = sys.argv[1]
    dataname = sys.argv[2]
    gene_count = int(sys.argv[3])

    clean_dataset(datadir, dataname, gene_count)


if __name__ == "__main__":
    main()
