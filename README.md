# Gene Expression Module Analysis

Four-step pipeline: select the most variable genes, build a TOM matrix, identify gene modules, and compute their eigengenes.

## Dependencies

- Python 3.13.13
- NumPy 2.2.6
- SciPy 1.17.1

```bash
python3 -m pip install numpy==2.2.6 scipy==1.17.1
```

Run the commands from the project root. `<datadir>` is the data directory and `<dataname>` is the CSV filename, including its extension.

## Input format

The input CSV contains gene names on the first row, an observation identifier in the first column, and numeric expression values in the remaining columns:

```text
,GENE_A,GENE_B
OBS_1,1.2,3.4
```

## Usage

### Step 1: select genes

Keeps the `K` genes with the highest variance.

```bash
python3 src/clean_dataset.py <datadir> <dataname> <K>
# Output: <datadir>/<dataname>_cleaned_<K>.csv
```

### Step 2: build the TOM matrix

Automatically determines `beta` (`1 <= beta <= 30`, `R^2 > 0.8`):

```bash
python3 src/build_graph.py <datadir> <dataname> POW
# Output: prints beta
```

Then builds the TOM matrix:

```bash
python3 src/build_graph.py <datadir> <dataname> TOM <beta>
# Output: <datadir>/<dataname>_TOM_<beta>.csv
```

### Step 3: identify modules

```bash
python3 src/find_modules.py <datadir> <dataname> static <tau> <l>
# Output: <datadir>/<dataname>_static_<l>_<tau>.csv

python3 src/find_modules.py <datadir> <dataname> adaptive <tau>
# Output: <datadir>/<dataname>_adaptive_<tau>.csv

python3 src/find_modules.py <datadir> <dataname> dynamic <tau>
# Output: <datadir>/<dataname>_dynamic_<tau>.csv
```

`tau` is an integer and `l` is a real number. Each output file contains one group identifier per gene, in the same order as the TOM matrix.

### Step 4: extract eigengenes

```bash
python3 src/extract_eigen.py <datadir> <dataname> <groupname>
# Output: <datadir>/<groupname>_eigen.csv
```

`<dataname>` is the expression file from step 1 and `<groupname>` is a file produced in step 3. Eigengene values are rounded to three decimal places.

## Tests

```bash
# Quick test
python3 tests/run_tests.py --quick --no-color

# Full test suite
python3 tests/run_tests.py --no-color

# Example: run step 1 only
python3 tests/run_tests.py --palier 1 --datasets small --gene-counts 50 --no-color
```

Test datasets and reference files are stored in `data/tests/`. Source files are in `src/` and the test runner is `tests/run_tests.py`.

## Auteurs / Authors

- Thomas Morin
- Clovis Thoréton

---

# Analyse de modules d'expression génétique

Pipeline en quatre paliers : sélection des gènes variables, construction d'une matrice TOM, recherche de modules, puis calcul de leurs eigengenes.

## Dépendances

- Python 3.13.13
- NumPy 2.2.6
- SciPy 1.17.1

```bash
python3 -m pip install numpy==2.2.6 scipy==1.17.1
```

Les commandes sont à exécuter depuis la racine du projet. `<datadir>` est le dossier des données et `<dataname>` le nom du fichier CSV, extension comprise.

## Format d'entrée

Le fichier CSV contient les noms des gènes sur la première ligne, un identifiant d'observation dans la première colonne, puis les valeurs numériques d'expression :

```text
,GENE_A,GENE_B
OBS_1,1.2,3.4
```

## Utilisation

### Palier 1 : sélection des gènes

Conserve les `K` gènes de plus forte variance.

```bash
python3 src/clean_dataset.py <datadir> <dataname> <K>
# Sortie : <datadir>/<dataname>_cleaned_<K>.csv
```

### Palier 2 : matrice TOM

Calcule automatiquement `beta` (`1 <= beta <= 30`, `R^2 > 0.8`) :

```bash
python3 src/build_graph.py <datadir> <dataname> POW
# Sortie : affiche beta
```

Construit ensuite la matrice TOM :

```bash
python3 src/build_graph.py <datadir> <dataname> TOM <beta>
# Sortie : <datadir>/<dataname>_TOM_<beta>.csv
```

### Palier 3 : modules

```bash
python3 src/find_modules.py <datadir> <dataname> static <tau> <l>
# Sortie : <datadir>/<dataname>_static_<l>_<tau>.csv

python3 src/find_modules.py <datadir> <dataname> adaptive <tau>
# Sortie : <datadir>/<dataname>_adaptive_<tau>.csv

python3 src/find_modules.py <datadir> <dataname> dynamic <tau>
# Sortie : <datadir>/<dataname>_dynamic_<tau>.csv
```

`tau` est un entier ; `l` est un nombre réel. Les fichiers produits contiennent un identifiant de groupe par gène, dans l'ordre des lignes et colonnes de la matrice TOM.

### Palier 4 : eigengenes

```bash
python3 src/extract_eigen.py <datadir> <dataname> <groupname>
# Sortie : <datadir>/<groupname>_eigen.csv
```

`<dataname>` est le fichier d'expression du palier 1 et `<groupname>` un fichier produit au palier 3. Les eigengenes sont arrondis à trois décimales.

## Tests

```bash
# Test rapide
python3 tests/run_tests.py --quick --no-color

# Tous les tests
python3 tests/run_tests.py --no-color

# Exemple : tester uniquement le palier 1
python3 tests/run_tests.py --palier 1 --datasets small --gene-counts 50 --no-color
```

Les jeux de données et fichiers de référence sont dans `data/tests/`. Le code source se trouve dans `src/` et le lanceur de tests dans `tests/run_tests.py`.

## Auteurs

- Thomas Morin
- Clovis Thoréton
