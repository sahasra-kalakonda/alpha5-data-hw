def richness(sample: dict[str, int]) -> int:
    """
    Returns the richness of the frequency map — the number of keys
    corresponding to nonzero values.

    Parameters:
    - sample: dict[str, int]

    Returns:
    - int: the richness of the sample
    """
    c = 0  # running count of nonzero species

    for key in sample:
        val = sample[key]
        if val < 0:
            raise ValueError("Error: negative value in frequency map given to richness()")
        if val > 0:
            c += 1

    return c


def jaccard_distance(map1: dict[str, int], map2: dict[str, int]) -> float:
    """
    Returns the Jaccard distance between two frequency maps.

    Parameters:
    - map1: dict[str, int]
    - map2: dict[str, int]

    Returns:
    - float: Jaccard distance between map1 and map2
    """
    c = sum_of_minima(map1, map2)
    t = sum_of_maxima(map1, map2)

    return 1 - (float(c)/float(t))


def sum_of_maxima(map1: dict[str, int], map2: dict[str, int]) -> int:
    """
    Returns the sum of per-key maxima across the union of two frequency maps
    (i.e., the sum used in the Jaccard distance denominator).

    Parameters:
    - map1: dict[str, int]
    - map2: dict[str, int]

    Returns:
    - int: sum of per-key maxima
    """
    c = 0  # running sum

    for key in map2:
        if key in map1:
            c += max_2(map1[key], map2[key])
        else:
            c += map2[key]

    for key in map1:
        if key not in map2:
            c += map1[key]

    # panic if c is equal to zero since we don't want 0/0
    if c == 0:
        raise ValueError("Error: no species common to two maps given to sum_of_maxima")

    return c


def max_2(n1: int, n2: int) -> int:
    """Returns the larger of two integers."""
    if n1 < n2:
        return n2
    else:
        return n1


def frequency_map(patterns: list[str]) -> dict[str, int]:
    """
    Produces the frequency map of a collection of patterns.

    Parameters:
    - patterns: list[str]

    Returns:
    - dict[str, int]: frequency map of strings to their counts in patterns
    """
    freq_map: dict[str, int] = {}
    for pattern in patterns:
        freq_map[pattern] = freq_map.get(pattern, 0) + 1
    return freq_map


def simpsons_index(sample: dict[str, int]) -> float:
    """
    Returns Simpson's index: the sum of (n/N)^2 where n is the count for a
    given species and N is the total number of individuals.

    Parameters:
    - sample: dict[str, int]

    Returns:
    - float: Simpson's index for the sample
    """
    total = sample_total(sample)
    simpson = 0.0

    if total == 0:
        raise ValueError("Error: Empty frequency map given to simpsons_index()!")

    for key in sample:
        val = sample[key]
        x = float(val) / float(total)
        simpson += x * x

    return simpson


def sample_total(freq_map: dict[str, int]) -> int:
    """Returns the sum of all counts in a frequency map."""
    total = 0  # running sum
    for key in freq_map:
        total += freq_map[key]
    return total


def bray_curtis_distance(map1: dict[str, int], map2: dict[str, int]) -> float:
    """
    Returns the Bray-Curtis distance between two frequency maps.

    Parameters:
    - map1: dict[str, int]
    - map2: dict[str, int]

    Returns:
    - float: Bray-Curtis distance between map1 and map2
    """
    c = sum_of_minima(map1, map2)
    s1 = sample_total(map1)
    s2 = sample_total(map2)

    if s1 == 0 or s2 == 0:
        raise ValueError("Error: sample given to bray_curtis_distance() has no positive values.")

    av = average(float(s1), float(s2))
    return 1 - (float(c)/av)


def average(x: float, y: float) -> float:
    """Returns the average of two floats."""
    return (x + y) / 2.0


def sum_of_minima(map1: dict[str, int], map2: dict[str, int]) -> int:
    """
    Returns the sum of per-key minima across the intersection of two frequency maps.

    Parameters:
    - map1: dict[str, int]
    - map2: dict[str, int]

    Returns:
    - int: sum of per-key minima
    """
    c = 0  # running sum

    for key in map1:
        if key in map2:
            c += min_2(map1[key], map2[key])

    return c


def min_2(n1: int, n2: int) -> int:
    """Returns the smaller of two integers."""
    if n1 < n2:
        return n1
    else:
        return n2


def main() -> None:
    pass


if __name__ == "__main__":
    main()
