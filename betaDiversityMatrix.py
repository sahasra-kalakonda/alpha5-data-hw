from helperFunctions import jaccard_distance, bray_curtis_distance, jaccard_binary_distance


def beta_diversity_matrix(all_maps: dict[str, dict[str, int]], distance_metric: str) -> tuple[list[str], list[list[float]]]:
    """
    Multiple Beta Diversity Problem.

    Takes a map of frequency maps along with a distance metric ("Bray-Curtis"
    or "Jaccard") and returns a sorted list of the sample names, along with a
    distance matrix whose (i, j) entry is the distance between sample i and
    sample j using the chosen metric.

    Parameters:
    - all_maps: dict[str, dict[str, int]]
    - distance_metric: str ("Bray-Curtis" or "Jaccard")

    Returns:
    - tuple[list[str], list[list[float]]]: sorted sample names and distance matrix
    """
    sample_names = sorted(all_maps.keys())  # fixed, reproducible ordering
    n = len(sample_names)

    # pick the distance function once, and reject anything we don't recognize
    # (rather than silently falling back to a default metric)
    if distance_metric == "Bray-Curtis":
        distance_func = bray_curtis_distance
    elif distance_metric == "Jaccard":
        distance_func = jaccard_distance              # abundance-weighted (Ruzicka)
    elif distance_metric == "Jaccard-binary":
        distance_func = jaccard_binary_distance       # classic presence/absence
    else:
        raise ValueError(f"Error: unknown distance_metric '{distance_metric}' given to beta_diversity_matrix()")

    # start with an n x n matrix of zeros
    matrix: list[list[float]] = [[0.0 for _ in range(n)] for _ in range(n)]

    # both Bray-Curtis and Jaccard are symmetric (distance(A, B) == distance(B, A)),
    # so we only need to compute the upper triangle (i <= j) and mirror each
    # result into the lower triangle, roughly halving the number of distance calls
    for i in range(n):
        for j in range(i, n):
            map1 = all_maps[sample_names[i]]
            map2 = all_maps[sample_names[j]]
            distance = distance_func(map1, map2)
            matrix[i][j] = distance
            matrix[j][i] = distance

    return sample_names, matrix
