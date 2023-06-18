from pathlib import Path
from typing import Tuple, Dict
import re

def fasta_to_dict(fasta_file: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    
    with open(fasta_file, 'r') as f:
            fasta_read = f.read()
    id_seq = ([id_seq.strip() for id_seq in fasta_read.split('>') if id_seq.strip()])
    id_seq_items = [[split for split in id_seq_pair.split('\n') if split] for id_seq_pair in id_seq]
    fasta_dict_full = {id_seq_pair[0]:id_seq_pair[1] for id_seq_pair in id_seq_items}

    # Keep gaps ('-') and inserts (lowercase) in fasta_dict
    # Remove the gaps and uppercase inserts in fasta_dict_no_gaps
    fasta_dict_no_gaps = {id: seq.replace('-', '').upper() for id, seq in fasta_dict_full.items()}
    
    return fasta_dict_full, fasta_dict_no_gaps

def fasta_dict_to_bootstrap_string(fasta_dict: Dict[str, str], n_bootstraps: int) -> str:

    bootstrap_fasta_out = ''
    for id, seq in fasta_dict.items():
        for i in range(n_bootstraps):
            bootstrap_fasta_out += f'>{id}#{i}\n'
            bootstrap_fasta_out += seq + '\n\n'
    
    return bootstrap_fasta_out