import functools
import itertools
import random
import re
import subprocess
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from tqdm.auto import tqdm

from structphy.executables import find_executable

RMSD_re = re.compile(r"RMSD=\W+([+-]?([0-9]*[.])?[0-9]+),")
TMscores_re = re.compile(r"TM-score=\W+([+-]?([0-9]*[.])?[0-9]+) \(")
identical_percent_re = re.compile(r"Seq_ID=n_identical/n_aligned=\W+([+-]?([0-9]*[.])?[0-9]+)\W")


# Wrapper for the tmalign binary
@functools.cache
def TMalign(pdb1: Path, pdb2: Path, tmalign_path: Path | str):

    command = [str(tmalign_path), str(pdb1), str(pdb2)]
    process = subprocess.run(command, capture_output=True, text=True, check=True)
    output = process.stdout

    rmsd_match = RMSD_re.search(output)
    identical_match = identical_percent_re.search(output)
    TMscores = TMscores_re.findall(output)
    if rmsd_match is None or identical_match is None or len(TMscores) < 2:
        message = f"Could not parse TM-align output for {pdb1} and {pdb2}."
        raise ValueError(message)
    RMSD = float(rmsd_match.group(1))
    TMscore_a = float(TMscores[0][0])
    TMscore_b = float(TMscores[1][0])
    identical_of_aligned = float(identical_match.group(1))

    return {
        "RMSD": RMSD,
        "TMscore_a": TMscore_a,
        "TMscore_b": TMscore_b,
        "identical_of_aligned": identical_of_aligned,
        "pdb_a": pdb1,
        "pdb_b": pdb2,
    }


# Wrapper required for multiprocessing
def tmalign_wrapper(args):
    return TMalign(args[0], args[1], args[2])


def generate_matrix_from_bootstraps(structure_files: list[Path], n_threads: int) -> pd.DataFrame:
    executable = find_executable("TMalign")

    # Get all combinations of structures, and add the TMalign binary location to the argument
    all_structure_combinations = [
        (first, second, executable) for first, second in itertools.combinations(structure_files, r=2)
    ]

    # Run all tmaligns multithreaded
    with Pool(n_threads) as pool:
        tm_results = list(
            tqdm(
                pool.imap_unordered(tmalign_wrapper, all_structure_combinations),
                total=len(all_structure_combinations),
                leave=False,
                desc="Current bootstrap",
                ascii=True,
            )
        )

    # Take the maximum score between two proteins as the re
    tm_scores = [
        {
            "max_score": max(tm_result["TMscore_a"], tm_result["TMscore_b"]),
            "pdb_a": tm_result["pdb_a"].name.split(".")[0],
            "pdb_b": tm_result["pdb_b"].name.split(".")[0],
        }
        for tm_result in tm_results
    ]

    # Load scores into a dataframe with columns ['pdb_a', 'pdb_b', 'max_score']
    tm_df = pd.DataFrame.from_records(tm_scores)

    # Make a copy with pdb_a and pdb_b swapped to fill out both triangles of the matrix
    tm_copy = tm_df.copy()
    a_col = tm_copy["pdb_a"].copy()
    tm_copy["pdb_a"] = tm_copy["pdb_b"]
    tm_copy["pdb_b"] = a_col
    tm_copy = tm_copy[tm_copy["pdb_a"] != tm_copy["pdb_b"]]
    tm_mix = pd.concat([tm_df, tm_copy])

    # Pivot and sort the score records to produce the square similarity matrix
    similarity_matrix = tm_mix.pivot_table(columns="pdb_b", index="pdb_a", values="max_score")
    similarity_matrix = similarity_matrix.reindex(sorted(similarity_matrix.columns), axis=1)
    similarity_matrix = similarity_matrix.sort_index()
    similarity_matrix.index.name = None
    similarity_matrix.columns.name = None

    # Make sure all entries are filled
    # nan_indices = np.where(pd.isnull(piv))
    # assert (nan_indices[0] == nan_indices[1]).all(), 'NaNs found outside the diagonal'

    # 1-similarity = distance
    # Convert similarity given by tmalign to a distance
    # Fill diagonal with 0s for self similarity
    distance_matrix = 1 - similarity_matrix
    distance_matrix = distance_matrix.fillna(0)
    return distance_matrix


def generate_bootstrap_matrices_from_structures(
    structure_files: list[Path], n_threads: int, n_bootstraps: int, *, seed: int | None = None
) -> list[pd.DataFrame]:

    # Set up ids_dict to sample a set of bootstrap structures for each matrix
    # Of form {id: [Path(id#0), Path(id#1), ...], ...}
    variants: dict[str, list[Path]] = {}
    for path in sorted(structure_files):
        variants.setdefault(path.stem.split("#")[0], []).append(path)

    rng = random.Random(seed)
    bootstrap_matrices = []
    for _ in tqdm(range(n_bootstraps), desc="Total bootstraps", ascii=True, position=0):
        chosen = [rng.choice(paths) for paths in variants.values()]
        bootstrap_matrices.append(generate_matrix_from_bootstraps(chosen, n_threads=n_threads))

    return bootstrap_matrices


def make_fake_outgroups(distance_matrices: list[pd.DataFrame], fake_outgroup: str) -> list[pd.DataFrame]:

    faked_distance_matrices = []
    for distance_matrix in distance_matrices:
        distance_matrix[fake_outgroup] = 1_000_000.0
        distance_matrix.loc[fake_outgroup] = 1_000_000.0
        distance_matrix.loc[fake_outgroup, fake_outgroup] = 0.0
        faked_distance_matrices.append(distance_matrix)

    return faked_distance_matrices
