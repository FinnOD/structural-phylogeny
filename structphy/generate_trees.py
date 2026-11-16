import subprocess
from multiprocessing import Pool
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd
from tqdm.auto import tqdm

from structphy.executables import find_executable


def fastme(phylip_matrix: str) -> str:
    executable = find_executable("fastme")
    with TemporaryDirectory() as directory:
        workdir = Path(directory)
        (workdir / "matrix").write_text(phylip_matrix)
        subprocess.run(
            [executable, "-i", "matrix", "-I", "info", "-O", "matrix.out", "-o", "tree", "-s", "-m", "NJ"],
            cwd=workdir,
            check=True,
            capture_output=True,
            text=True,
        )
        return (workdir / "tree").read_text().strip()


def convert_matrix_to_phylip(distance_df: pd.DataFrame) -> str:
    rows = [str(len(distance_df))]
    for name, row in distance_df.iterrows():
        label = str(name).split("#")[0].replace(" ", "_")
        values = " ".join(f"{value:.8f}" for value in row)
        rows.append(f"{label} {values}")
    return "\n".join(rows) + "\n"


def matrix_to_fastme_newick(distance_df: pd.DataFrame) -> str:
    return fastme(convert_matrix_to_phylip(distance_df))


def matrices_to_fastme_newick(distance_dfs: list[pd.DataFrame], n_threads: int) -> list[str]:
    with Pool(n_threads) as pool:
        return list(
            tqdm(
                pool.imap(matrix_to_fastme_newick, distance_dfs),
                total=len(distance_dfs),
                desc="Building trees",
                ascii=True,
            )
        )
