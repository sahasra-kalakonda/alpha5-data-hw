import os


def read_samples_from_directory(directory: str) -> dict[str, dict[str, int]]:
    """
    Reads all .txt files in the given directory as frequency maps.

    Parameters:
    - directory: str

    Returns:
    - dict[str, dict[str, int]]: map of sample names to frequency maps
    """
    all_maps: dict[str, dict[str, int]] = {}
    for filename in os.listdir(directory):
        if filename.endswith(".txt"):
            # Remove the file extension (".txt") to obtain the name of the sample
            sample_name = filename.replace(".txt", "")
            all_maps[sample_name] = read_freq_map_from_file(os.path.join(directory, filename))
    return all_maps


def read_freq_map_from_file(filename: str) -> dict[str, int]:
    """
    Reads a file with one string per line and returns its frequency map.

    Parameters:
    - filename: str

    Returns:
    - dict[str, int]: frequency map of strings in the file
    """
    freq_map: dict[str, int] = {}
    with open(filename, 'r') as file:
        for line in file:
            val = line.strip()
            freq_map[val] = freq_map.get(val, 0) + 1
    return freq_map


def main() -> None:
    pass


if __name__ == "__main__":
    main()
