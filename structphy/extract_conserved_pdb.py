from pathlib import Path

aa_dict = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def remove_inserts_from_structure(pdb_filename: Path, manual_alignment: str, output_pdb_filename: Path):

  with open(pdb_filename) as f:
    full_pdb_lines = f.readlines()

  # Check first to assure the alignment is correct between sequence and PDB
  lines = [x for x in full_pdb_lines if 'ATOM' in x]
  pdb_seq = {int(x[22:30]):aa_dict[x[17:20]] for x in lines}
  fasta_seq = manual_alignment.replace('-', '').strip()
  pdb_seq_string = "".join(pdb_seq[key] for key in sorted(pdb_seq))
  assert (pdb_seq_string == fasta_seq.upper())


  out_lines = []
  header_finished = False
  for line in full_pdb_lines:
    # Keep the header
    if ('ATOM' not in line):
      if not header_finished:
        out_lines.append(line)
      continue
    header_finished = True

    # Extract PDB atom features using columns of PDB file
    pdb_atm_num = int(line[6:11])
    pdb_res_num = int(line[22:30])
    pdb_res_long = line[17:20]
    pdb_res = aa_dict[pdb_res_long]

    fasta_res = fasta_seq[pdb_res_num-1]
    assert fasta_res.upper() == pdb_res
    if fasta_res.isupper():
      out_lines.append(line)

  ter_line = f'TER    {pdb_atm_num}      {pdb_res_long} A {pdb_res_num}\n'
  end_line = 'END\n'
  out_lines.append(ter_line)
  out_lines.append(end_line)
  
  with open(output_pdb_filename, 'w') as f:
    f.writelines(out_lines)
