from helperFunctions import simpsons_index


def simpsons_map(all_maps: dict[str, dict[str, int]]) -> dict[str, float]:
    """
    Multiple Evenness Problem.

    Takes a map of frequency maps and returns a map whose values are
    Simpson's index of each sample.

    Parameters:
    - all_maps: dict[str, dict[str, int]]

    Returns:
    - dict[str, float]: map of sample names to their Simpson's index values
    """
    result: dict[str, float] = {}  # sample name -> Simpson's index

    for sample_name, sample_map in all_maps.items():
        result[sample_name] = simpsons_index(sample_map)

    return result
