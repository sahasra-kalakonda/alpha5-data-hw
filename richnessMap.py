from helperFunctions import richness


def richness_map(all_maps: dict[str, dict[str, int]]) -> dict[str, int]:
    """
    Multiple Richness Problem.

    Takes a map of frequency maps and returns a map whose values are the
    richness of each sample.

    Parameters:
    - all_maps: dict[str, dict[str, int]]

    Returns:
    - dict[str, int]: map of sample names to their richness values
    """
    result: dict[str, int] = {}  # sample name -> richness

    for sample_name, sample_map in all_maps.items():
        result[sample_name] = richness(sample_map)

    return result
