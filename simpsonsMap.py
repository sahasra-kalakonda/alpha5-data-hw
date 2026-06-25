from helperFunctions import simpsons_index


def simpsons_map(all_maps: dict[str, dict[str, int]]) -> dict[str, float]:
    """
    Multiple Evenness Problem.

    Takes a map of frequency maps and returns a map whose values are the
    GINI-SIMPSON index 1 - D of each sample (D = sum of (n/N)^2). This is an
    evenness/diversity measure: higher means the community is more even (two random
    reads are more likely to be different species). It is NOT the raw Simpson's index
    D itself -- see helperFunctions.simpsons_index for the distinction.

    Parameters:
    - all_maps: dict[str, dict[str, int]]

    Returns:
    - dict[str, float]: map of sample names to their Gini-Simpson index (1 - D)
    """
    result: dict[str, float] = {}  # sample name -> Gini-Simpson index (1 - D)

    for sample_name, sample_map in all_maps.items():
        result[sample_name] = 1.0 - simpsons_index(sample_map)

    return result
