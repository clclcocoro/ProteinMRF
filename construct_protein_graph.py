#!/usr/bin/env python

"""Construct Protein Graph.


Usage:
  construct_protein_graph.py <binding_residue_file> <pdb_file_dir>
  construct_protein_graph.py (-h | --help)
  construct_protein_graph.py --version

"""

import os
import re
from docopt import docopt
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from Bio import pairwise2


def construct_protein_graph(structure, pdbidx_to_bindresidx, atom_type, delta):
    graph = []
    for a_i, residue_a in enumerate(structure.get_residues()):
        for b_i, residue_b in enumerate(structure.get_residues()):
            if residue_a.id[0] != ' ' or residue_b.id[0] != ' ': # HETATM
                continue
            if residue_a.id[1] == residue_b.id[1]:
                continue
            if residue_a.child_dict[atom_type] - residue_b.child_dict[atom_type] <= delta:
                if a_i in pdbidx_to_bindresidx and b_i in pdbidx_to_bindresidx:
                    graph.append((pdbidx_to_bindresidx[a_i], pdbidx_to_bindresidx[b_i]))
    return graph


def get_sequence_from_pdb_structure(structure):
    sequence = ''
    for residue in structure.get_residues():
        if residue.id[0] == ' ':
            sequence += seq1(residue.resname)
    return sequence


def get_sequence_from_bindres_file(bindres_file):
    with open(bindres_file) as fp:
        for i, line in enumerate(fp):
            if i == 1:
                return line.rstrip()


def get_index_map(sequence_bindres, sequence_pdb):
    alignment = pairwise2.align.globalmd(sequence_bindres, sequence_pdb, 0, -1, -100, -100, -0.5, -0.1)
    pdbidx_to_bindresidx = {}
    index_of_pdb = 0
    index_of_bindres = 0
    for amino_bindres, amino_pdb in zip(alignment[0][0], alignment[0][1]):
        if amino_bindres != '-' and amino_pdb != '-':
            pdbidx_to_bindresidx[index_of_pdb] = index_of_bindres
        if amino_bindres != '-':
            index_of_bindres += 1
        if amino_pdb != '-':
            index_of_pdb += 1
    return pdbidx_to_bindresidx


if __name__ == "__main__":
    arguments = docopt(__doc__)
    bindres_file  = arguments['<binding_residue_file>']
    pdb_file_dir  = arguments['<pdb_file_dir>']
    delta = 5
    atom_type = 'CA'

    parser = PDBParser()
    pdbid_chain = re.match(r'(.+)\.bindres', os.path.basename(bindres_file)).group(1)
    if len(pdbid_chain) == 5:
        pdbid = pdbid_chain[:4]
        chain = pdbid_chain[4]
    elif len(pdbid_chain) == 6:
        pdbid = pdbid_chain[:4]
        chain = pdbid_chain[5]
    else:
        raise ValueError("bindres_file is invalid {}".format(bindres_file))

    sequence_bindres = get_sequence_from_bindres_file(bindres_file)
    structure = parser.get_structure(pdbid, "{}/{}.pdb".format(pdb_file_dir, pdbid))
    structure_of_chain = structure[0][chain]
    sequence_pdb = get_sequence_from_pdb_structure(structure_of_chain)
    pdbidx_to_bindresidx = get_index_map(sequence_bindres, sequence_pdb)
    protein_graph = construct_protein_graph(structure_of_chain, pdbidx_to_bindresidx, atom_type, delta)
    for edge in protein_graph:
        print edge[0], edge[1]
