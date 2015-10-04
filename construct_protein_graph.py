#!/usr/bin/env python

"""Construct Protein Graph.


Usage:
  construct_protein_graph.py <pdb_file>
  construct_protein_graph.py (-h | --help)
  construct_protein_graph.py --version

"""

import os
import re
from docopt import docopt
from Bio.PDB import PDBParser


def construct_protein_graph(structure, atom_type, delta):
    graph = []
    for residue_a in structure.get_residues():
        for residue_b in structure.get_residues():
            if residue_a.id[0] != ' ' or residue_b.id[0] != ' ': # HETATM
                continue
            if residue_a.id[1] == residue_b.id[1]:
                continue
            if residue_a.child_dict[atom_type] - residue_b.child_dict[atom_type] <= delta:
                graph.append((residue_a.id[1], residue_b.id[1]))
    return graph


if __name__ == "__main__":
    arguments = docopt(__doc__)
    pdb_file  = arguments['<pdb_file>']
    delta = 5
    atom_type = 'CA'

    parser = PDBParser()
    pdbid = re.match(r'(.+)\.pdb', os.path.basename(pdb_file)).group(1)
    structure = parser.get_structure(pdbid, pdb_file)
    protein_graph = construct_protein_graph(structure, atom_type, delta)
    for edge in protein_graph:
        print edge[0], edge[1]
