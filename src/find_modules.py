import sys
import numpy as np
from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import squareform
from csv_utils import build_output_path, read_matrix_only_csv, write_group_assignments


def build_dendrogram(tom_matrix):
    """Palier 3, Principe : Maintenant que l'on a un graphe, sous forme de matrice de
    connectivité, on va chercher les ensembles de gènes très corrélés entre eux. Cela va nous
    permettre de construire des "modules" de gènes se comportant de manière similaire [3].
    Étape 1 : construction de l'arbre initial On commence par calculer un arbre de distances
    entre gènes appelé dendrogramme. Vous pouvez utiliser la fonction linkage de scipy.hierarchy
    en utilisant les options method="average" et optimal_ordering=True sur 1 - T où T est la
    matrice de topologie recouvrante obtenue à la fin du palier précédent (Attention : le format
    d'entrée de linkage est une matrice sous forme compressée). Pour obtenir un arbre T à partir
    de la matrice de liaisons on peut utiliser la fonction to_tree de scipy.cluster.hierarchy.
    """
    distance_matrix = 1 - tom_matrix
    np.fill_diagonal(distance_matrix, 0)
    condensed_distance = squareform(distance_matrix, checks=False)
    linkage_matrix = linkage(condensed_distance, method="average", optimal_ordering=True)
    tree = to_tree(linkage_matrix)

    return tree


def extract_heights_and_leaf_map(tree):
    """Palier 3, Étape 1 : Une fois qu'on a cet arbre, on veut obtenir la liste des distances
    associées à chaque noeud interne dans l'ordre de T (c'est à dire un parcours profondeur main
    gauche infixe). On appelle cette liste la liste des hauteurs h.
    """
    heights = []
    leaf_parent_index = {}
    stack = []
    node = tree

    while stack or node is not None:
        while node is not None:
            stack.append(node)
            node = node.left
        node = stack.pop()

        if not node.is_leaf():
            parent_index = len(heights)
            heights.append(node.dist)
            if node.left.is_leaf():
                leaf_parent_index[node.left.id] = parent_index
            if node.right.is_leaf():
                leaf_parent_index[node.right.id] = parent_index

        node = node.right

    return heights, leaf_parent_index


def measure_backward_run(h_bar, position):
    """Palier 3, Étape 2 : sa taille t_{k} qui correspond au nombres de positions successives
    < k de même signe que h^{-}_{k}, et son point de cassure (breakpoint) b_{k} = k - t_{k}.
    Renvoie donc la taille t_{k} et la position b_{k} associée à la transition descendante en k.
    """
    sign_at_position = np.sign(h_bar[position])
    size = 0
    current = position
    while current - 1 >= 0 and np.sign(h_bar[current - 1]) == sign_at_position:
        size += 1
        current -= 1

    return size, current


def measure_forward_run(h_bar, position):
    """Palier 3, Étape 2 : symétrique de measure_backward_run pour une transition montante
    (h^{-} passant de négatif à positif). Compte les positions successives > k de même signe
    que h^{-}_{k}, afin de mesurer la taille d'un groupe qui débute à cette transition même
    lorsqu'aucune transition descendante ne le referme.
    """
    sign_at_position = np.sign(h_bar[position])
    size = 0
    current = position
    height_count = len(h_bar)
    while current + 1 < height_count and np.sign(h_bar[current + 1]) == sign_at_position:
        size += 1
        current += 1

    return size


def find_transition_breakpoints(heights, l_value, tau):
    """Palier 3, Étape 2 : A partir du tableau des hauteurs h et d'un niveau de hauteur l, on
    défini le tableau :

        h^{-} = h - l

    dans ce tableau, on cherche les points transitions correspondant à l'ensemble des positions
    k tel que h^{-}_{k} > 0 et h^{-}_{k+1} < 0. A chaque transition, on associe deux éléments
    (cf. Figure 1) : sa taille t_{k} qui correspond au nombres de positions successives < k de
    même signe que h^{-}_{k}, et son point de cassure (breakpoint) b_{k} = k - t_{k}, la
    première position où la valeur de h^{-} devient du même signe que h^{-}_{k}. Si on atteint
    le début de la liste sans trouver une telle position, alors on choisira b_{k} = 0.
    """
    # h^{-} = h - l
    h_bar = np.array(heights) - l_value
    height_count = len(h_bar)
    breakpoints = {0, height_count}

    position_upper_bound = height_count - 1
    for position in range(position_upper_bound):
        next_position = position + 1
        # h^{-}_{k} > 0 et h^{-}_{k+1} < 0
        if h_bar[position] > 0 and h_bar[next_position] < 0:
            # t_{k}, b_{k} = k - t_{k}
            run_size, breakpoint_position = measure_backward_run(h_bar, position)
            # t_{k} > τ
            if run_size > tau:
                breakpoints.add(breakpoint_position)
        # h^{-}_{k} < 0 et h^{-}_{k+1} > 0
        elif h_bar[position] < 0 and h_bar[next_position] > 0:
            run_size = measure_forward_run(h_bar, next_position)
            if run_size > tau:
                breakpoints.add(next_position)

    sorted_breakpoints = sorted(breakpoints)

    return sorted_breakpoints


def assign_groups(breakpoints, height_count):
    """Palier 3, Étape 2 : Pour trouver les groupes, on va chercher tous les points de
    cassures dans h^{-}, puis garder uniquement les points significatifs, dont la taille de la
    transition associée est > τ. Tous les éléments entre deux points de cassures significatifs
    font parti du même groupe. Attention : Pour éviter de perdre des éléments, on va toujours
    considérer que la position 0 est le premier point de cassure significatif et |h^{-}| - 1 le
    dernier.
    """
    group_of_position = [0] * height_count
    breakpoint_pairs = zip(breakpoints, breakpoints[1:])
    for group_id, (start, end) in enumerate(breakpoint_pairs):
        for position in range(start, end):
            group_of_position[position] = group_id

    return group_of_position


def find_contiguous_runs(group_ids):
    """Palier 3, Étape 2 : regroupe les positions consécutives partageant le même id de
    groupe, nécessaire pour repérer les groupes de taille 1 dans merge_lone_dendrogram_leaves.
    """
    runs = []
    start = 0
    for position in range(1, len(group_ids) + 1):
        if position == len(group_ids) or group_ids[position] != group_ids[start]:
            runs.append((start, position))
            start = position

    return runs


def merge_lone_dendrogram_leaves(group_ids):
    """Palier 3, Étape 2 : post-traitement de la découpe statique (cf. article [3], TreecutCore,
    section 2.1.1 du supplément), donc appliqué à chaque découpe de base et hérité par les
    découpes adaptative et dynamique. Fusionne chaque groupe de taille 1 avec son voisin dans
    l'ordre du dendrogramme (celui de droite, ou celui de gauche en dernière position), répété
    jusqu'à ce qu'il ne reste plus de groupe de taille 1.
    """
    group_ids = list(group_ids)
    singleton_found = True
    while singleton_found:
        singleton_found = False
        runs = find_contiguous_runs(group_ids)
        for run_index, (start, end) in enumerate(runs):
            if end - start != 1:
                continue
            neighbor_index = run_index + 1 if run_index + 1 < len(runs) else run_index - 1
            neighbor_start, _ = runs[neighbor_index]
            group_ids[start] = group_ids[neighbor_start]
            singleton_found = True
            break

    return group_ids


def renumber_consecutively(group_ids):
    """Palier 3, Étape 2 : réattribue des identifiants consécutifs (0, 1, 2, ...) dans l'ordre
    d'apparition, pour donner un étiquetage canonique après fusion des groupes de taille 1.
    """
    remap = {}
    renumbered = []
    for group_id in group_ids:
        if group_id not in remap:
            remap[group_id] = len(remap)
        renumbered.append(remap[group_id])

    return renumbered


def static_cut(heights, l_value, tau):
    """Palier 3, Étape 2 : Découpe statique de l'arbre L'arbre généré par linkage donne une
    hiérarchie de distances entre noeuds mais ne fourni pas de groupes. Pour les obtenir, on va
    découper l'arbre à un certain niveau de hauteur. Le post-traitement de l'article [3]
    (TreecutCore, section 2.1.1 du supplément) fait partie de cette découpe de base : les groupes
    de taille 1 sont rattachés à un groupe voisin.
    """
    breakpoints = find_transition_breakpoints(heights, l_value, tau)
    height_count = len(heights)
    group_ids = assign_groups(breakpoints, height_count)
    merged_group_ids = merge_lone_dendrogram_leaves(group_ids)
    renumbered_group_ids = renumber_consecutively(merged_group_ids)

    return renumbered_group_ids


def adaptive_cut(heights, tau):
    """Palier 3, Étape 3 : Découpe adaptative La méthode de découpe statique n'indique pas
    comment trouver le seuil de découpe l. Pour éviter les problèmes liés au choix d'un seuil
    arbitraire, on va calculer des valeurs de seuils adaptées pour chaque groupe. Soit un groupe
    H composées de n distances h1, ... , hn alors on a :

        l_{m} = (1/n) Σ_{i=1}^n h_{i}
        l_{u} = (1/2)[l_{m} + max(h1, ... , hn)]
        l_{d} = (1/2)[l_{m} + min(h1, ... , hn)]

    A partir de ces quantités, on définie la procédure de découpe adaptative pour un groupe de
    hauteurs H = h1, ... , hn : 1. Calculer l_{m}, l_{d}, l_{u} 2. G = Découpe statique de H
    avec le seuil l_{m} 3. G = G' 4. Si |G| == 1 alors G' = Découpe statique de H avec le seuil
    l_{d} 5. Si G == G' alors G' = Découpe statique de H avec le seuil l_{u} 6. retourner G' où
    |G| représente le nombre de groupes différents dans G.
    """
    # 1
    average_weight = 0.5
    # l_{m}
    mean_height = np.mean(heights)
    # l_{u}
    upper_threshold = average_weight * (mean_height + np.max(heights))
    # l_{d}
    lower_threshold = average_weight * (mean_height + np.min(heights))

    # 2
    groups = static_cut(heights, mean_height, tau)

    # 3
    refined_groups = groups

    # 4
    distinct_group_count = len(set(groups))
    if distinct_group_count == 1:
        refined_groups = static_cut(heights, lower_threshold, tau)

    # 5
    if groups == refined_groups:
        refined_groups = static_cut(heights, upper_threshold, tau)

    # 6
    return refined_groups


def group_positions(group_ids):
    """Palier 3, Étape 4, point 2 : Appeler récursivement la procédure de découpe adaptative
    (étape 3) sur H(k) et tous les sous groupes ainsi générés jusqu'à convergence.
    """
    positions_by_group = {}
    for position, group_id in enumerate(group_ids):
        positions_by_group.setdefault(group_id, []).append(position)

    return positions_by_group


def merge_recursive_split(heights, tau, first_split_groups):
    """Palier 3, Étape 4, point 2 : Pour chaque groupe H(k) ∈ Ω : Appeler récursivement la
    procédure de découpe adaptative (étape 3) sur H(k) et tous les sous groupes ainsi générés
    jusqu'à convergence, quand plus aucun nouveau groupe n'est formé.
    """
    final_groups = [0] * len(heights)
    next_id = 0

    positions_by_group = group_positions(first_split_groups)
    for positions in positions_by_group.values():
        sub_heights = [heights[position] for position in positions]
        refined_groups = refine_adaptive_recursive(sub_heights, tau)
        sub_groups = renumber_consecutively(refined_groups)
        for position, sub_group_id in zip(positions, sub_groups):
            final_groups[position] = next_id + sub_group_id
        distinct_sub_groups = set(sub_groups)
        next_id += len(distinct_sub_groups)

    return final_groups


def refine_adaptive_recursive(heights, tau):
    """Palier 3, Étape 4, point 2 : Pour chaque groupe H(k) ∈ Ω : Appeler récursivement la
    procédure de découpe adaptative (étape 3) sur H(k) et tous les sous groupes ainsi générés
    jusqu'à convergence, quand plus aucun nouveau groupe n'est formé.
    """
    if len(heights) <= 1:
        single_element_groups = [0] * len(heights)
        return single_element_groups

    groups = adaptive_cut(heights, tau)
    distinct_group_count = len(set(groups))
    if distinct_group_count == 1:
        return groups

    refined_groups = merge_recursive_split(heights, tau, groups)

    return refined_groups


def dynamic_cut(heights, tau):
    """Palier 3, Étape 4 : Découpe dynamique de l'arbre Finalement, à partir des deux
    procédures précédentes on peut définir une procédure de découpe dynamique d'un tableau de
    hauteurs H = h1, ... hn : 1. Déterminer un premier ensemble de groupes Ω = H(1), ...,
    H(K) avec une découpe statique utilisant une haute valeur de seuil l = 0.99 x max(H).
    """
    if len(heights) == 0:
        return []

    initial_height_ratio = 0.99
    initial_threshold = initial_height_ratio * max(heights)
    initial_groups = static_cut(heights, initial_threshold, tau)
    final_groups = merge_recursive_split(heights, tau, initial_groups)

    return final_groups


def cut_heights(option, heights, tau, l_value):
    """Palier 3, Paramètres d'entrée : option peut être soit static, adaptive ou dynamic, avec :
    Pour l'option static, on donne le paramètre tau (type int) correspondant au seuil de
    significativité des points de cassures et le paramètre l (type float) correspondant au
    seuil de découpe souhaité et . Pour l'option dynamic, on donne le paramètre tau (type int)
    correspondant au seuil de significativité des points de cassures. Pour l'option adaptive on
    donne le paramètre tau (type int) correspondant au seuil de significativité des points de
    cassures.
    """
    if option == "static":
        group_ids = static_cut(heights, l_value, tau)
    elif option == "adaptive":
        group_ids = adaptive_cut(heights, tau)
    else:
        group_ids = dynamic_cut(heights, tau)

    return group_ids


def propagate_to_leaves(group_of_node, leaf_parent_index, gene_count):
    """Palier 3 : Finalement, quelque soit la méthode, une fois le tableau des groupes obtenu
    il faut propager les groupes aux feuilles de l'arbre. Le groupe d'une feuille correspondant
    au groupe de son noeud parent. Attention : l'ordre dans le tableau final doit être le même
    que celui de la matrice d'entrée.
    """
    output = [0] * gene_count
    for leaf_id, parent_index in leaf_parent_index.items():
        output[leaf_id] = group_of_node[parent_index]

    return output


def build_output_filename(dataname, option, tau, l_value):
    """Palier 3, Format de sortie : Un fichier texte nommé [dataname]_[option]_[l]_[tau].csv
    dans le cas statique ou dataname_[option]_[tau].csv dans les autres cas (les éléments entre
    crochets sont à remplacés par les valeurs des paramètres associés). Ce fichier contient sur
    chaque ligne k l'id du groupe auquel appartient le kème gène de la matrice TOM. L'ordre doit
    être le même que celui des lignes et colonnes de la matrice d'entrée.
    """
    if option == "static":
        output_filename = f"{dataname}_{option}_{l_value}_{tau}.csv"
    else:
        output_filename = f"{dataname}_{option}_{tau}.csv"

    return output_filename


def run_find_modules(datadir, dataname, option, cut_params):
    """Palier 3, Paramètres d'entrée : datadir dataname option tau [l], où datadir est le
    chemin vers le dossier contenant les données, dataname le nom du fichier csv contenant les
    données (au format du palier précédent), option peut être soit static, adaptive ou dynamic,
    avec : Pour l'option static, on donne le paramètre tau (type int) correspondant au seuil de
    significativité des points de cassures et le paramètre l (type float) correspondant au
    seuil de découpe souhaité et . Pour l'option dynamic, on donne le paramètre tau (type int)
    correspondant au seuil de significativité des points de cassures. Pour l'option adaptive on
    donne le paramètre tau (type int) correspondant au seuil de significativité des points de
    cassures.
    """
    tau, l_value = cut_params
    input_path = build_output_path(datadir, dataname)
    tom_matrix = read_matrix_only_csv(input_path)
    gene_count = tom_matrix.shape[0]

    tree = build_dendrogram(tom_matrix)
    heights, leaf_parent_index = extract_heights_and_leaf_map(tree)

    group_of_node = cut_heights(option, heights, tau, l_value)
    groups = propagate_to_leaves(group_of_node, leaf_parent_index, gene_count)

    output_filename = build_output_filename(dataname, option, tau, l_value)
    output_path = build_output_path(datadir, output_filename)
    write_group_assignments(output_path, groups)


def main():
    """Palier 3, Paramètres d'entrée : datadir dataname option tau [l]"""
    datadir = sys.argv[1]
    dataname = sys.argv[2]
    option = sys.argv[3]
    tau = int(sys.argv[4])
    l_value = float(sys.argv[5]) if option == "static" else None

    run_find_modules(datadir, dataname, option, (tau, l_value))


if __name__ == "__main__":
    main()
