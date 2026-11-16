from io import StringIO
from pathlib import Path


def fasta_to_dict(fasta_file: Path) -> tuple[dict[str, str], dict[str, str]]:
    return parse_fasta(fasta_file.read_text(encoding="utf-8-sig"))


def parse_fasta(text: str) -> tuple[dict[str, str], dict[str, str]]:
    """Return sequences keyed by full header, with and without alignment markup."""
    records: dict[str, list[str]] = {}
    parts: list[str] | None = None
    for line_number, raw_line in enumerate(StringIO(text, newline=None), start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            header = line[1:].strip()
            if not header:
                message = f"Empty FASTA header at line {line_number}."
                raise ValueError(message)
            if header in records:
                message = f"Duplicate FASTA header {header!r} at line {line_number}."
                raise ValueError(message)
            parts = []
            records[header] = parts
        else:
            if parts is None:
                message = f"Expected a FASTA header before sequence data at line {line_number}."
                raise ValueError(message)
            if ">" in line:
                message = f"Unexpected '>' in FASTA sequence data at line {line_number}."
                raise ValueError(message)
            parts.append("".join(line.split()))

    if not records:
        message = "The input FASTA contains no sequences."
        raise ValueError(message)
    sequences = {}
    for header, parts in records.items():
        if not parts:
            message = f"FASTA record {header!r} contains no sequence."
            raise ValueError(message)
        sequences[header] = "".join(parts)

    ungapped = {header: seq.replace("-", "").upper() for header, seq in sequences.items()}
    return sequences, ungapped


def fasta_dict_to_bootstrap_string(fasta_dict: dict[str, str], n_bootstraps: int) -> str:

    bootstrap_fasta_out = ""
    for id, seq in fasta_dict.items():
        for i in range(n_bootstraps):
            bootstrap_fasta_out += f">{id}#{i}\n"
            bootstrap_fasta_out += seq + "\n\n"

    return bootstrap_fasta_out
