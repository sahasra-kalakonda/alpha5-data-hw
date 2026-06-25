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
    Returns the (abundance-weighted) Jaccard distance between two frequency maps.

    NOTE ON THE METRIC: because we use the COUNTS in the min/max sums, this is the
    abundance-weighted Jaccard, also known as the Ruzicka distance:

        d = 1 - sum_i min(a_i, b_i) / sum_i max(a_i, b_i)

    It only reduces to the classic presence/absence Jaccard (|A intersect B| / |A union B|)
    if every count is 1. Keep that in mind when comparing it to Bray-Curtis: both of our
    distances are abundance-based, so they are expected to agree closely.

    Parameters:
    - map1: dict[str, int]
    - map2: dict[str, int]

    Returns:
    - float: abundance-weighted Jaccard (Ruzicka) distance between map1 and map2
    """
    c = sum_of_minima(map1, map2)
    t = sum_of_maxima(map1, map2)

    return 1 - (float(c)/float(t))


def jaccard_binary_distance(map1: dict[str, int], map2: dict[str, int]) -> float:
    """
    Returns the CLASSIC presence/absence Jaccard distance between two frequency maps:

        d = 1 - |A intersect B| / |A union B|

    where A and B are the SETS of OTUs present (count > 0) in each sample. Unlike
    jaccard_distance (which is abundance-weighted Ruzicka), this ignores how many reads
    each OTU has -- it only asks WHICH species are shared. We use it as a genuinely
    independent check on Bray-Curtis: BC and Ruzicka share the same min-overlap numerator,
    but this binary metric does not, so agreement here is real evidence, not a tautology.

    Parameters:
    - map1: dict[str, int]
    - map2: dict[str, int]

    Returns:
    - float: presence/absence Jaccard distance between map1 and map2
    """
    a = {k for k, v in map1.items() if v > 0}
    b = {k for k, v in map2.items() if v > 0}
    union = a | b
    if not union:
        raise ValueError("Error: both samples empty in jaccard_binary_distance (cannot divide by 0)")
    return 1.0 - len(a & b) / len(union)


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

    # panic if c is equal to zero since we don't want 0/0. Note: because this sum is over
    # the UNION of both maps, c == 0 happens ONLY when BOTH maps are empty -- two disjoint
    # but non-empty samples still give c > 0.
    if c == 0:
        raise ValueError("Error: both frequency maps are empty in sum_of_maxima (cannot divide by 0)")

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
    Returns Simpson's index D = sum of (n/N)^2, where n is the count for a given
    species and N is the total number of individuals.

    INTERPRETATION: D is the probability that two reads drawn at random are the SAME
    species (a dominance/concentration measure) -- HIGHER D means LESS diverse. The
    complementary Gini-Simpson index 1 - D (computed in simpsonsMap.simpsons_map) is the
    probability they are DIFFERENT species -- HIGHER means MORE even/diverse. We report
    1 - D in our tables and plots because "higher = more diverse" reads more naturally.

    Parameters:
    - sample: dict[str, int]

    Returns:
    - float: Simpson's index D (dominance) for the sample
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
