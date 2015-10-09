#!/usr/bin/env python

import sys
sys.path.append("..")

import unittest
import construct_protein_graph
from Bio.PDB import PDBParser

class TestConstructProteinGraph(unittest.TestCase):

    def test_get_sequence_from_bindres_file(self):
        bindres_file = "./test1.bindres"
        sequence = construct_protein_graph.get_sequence_from_bindres_file(bindres_file)
        self.assertEqual("VNIKTNPFK", sequence)

    def test_get_sequence_from_pdb_structure(self):
        pdb_file = "./test.pdb"
        p = PDBParser()
        structure = p.get_structure('test', pdb_file)
        structure_of_chain = structure[0]['A']
        sequence = construct_protein_graph.get_sequence_from_pdb_structure(structure_of_chain)
        self.assertEqual("VNIKTNPFK", sequence)

    def test_get_index_map(self):
        bindres_file = "./test1.bindres"
        sequence_bindres = construct_protein_graph.get_sequence_from_bindres_file(bindres_file)
        pdb_file = "./test.pdb"
        p = PDBParser()
        structure = p.get_structure('test', pdb_file)
        structure_of_chain = structure[0]['A']
        sequence_pdb = construct_protein_graph.get_sequence_from_pdb_structure(structure_of_chain)
        pdbidx_to_bindresidx = construct_protein_graph.get_index_map(sequence_bindres, sequence_pdb)
        correct_bindres_idx = [i for i in xrange(len(sequence_pdb))]
        for i in xrange(len(sequence_pdb)):
            self.assertEqual(correct_bindres_idx[i], pdbidx_to_bindresidx[i])
            
        bindres_file = "./test2.bindres"
        sequence_bindres = construct_protein_graph.get_sequence_from_bindres_file(bindres_file)
        pdbidx_to_bindresidx = construct_protein_graph.get_index_map(sequence_bindres, sequence_pdb)
        correct_bindres_idx = [i for i in xrange(len(sequence_pdb))]
        for i in xrange(len(sequence_pdb)):
            self.assertEqual(correct_bindres_idx[i], pdbidx_to_bindresidx[i])
 
        bindres_file = "./test3.bindres"
        sequence_bindres = construct_protein_graph.get_sequence_from_bindres_file(bindres_file)
        pdbidx_to_bindresidx = construct_protein_graph.get_index_map(sequence_bindres, sequence_pdb)
        correct_bindres_idx = [i+1 for i in xrange(len(sequence_pdb))]
        for i in xrange(len(sequence_pdb)):
            self.assertEqual(correct_bindres_idx[i], pdbidx_to_bindresidx[i])
 
        bindres_file = "./test4.bindres"
        sequence_bindres = construct_protein_graph.get_sequence_from_bindres_file(bindres_file)
        pdbidx_to_bindresidx = construct_protein_graph.get_index_map(sequence_bindres, sequence_pdb)
        correct_bindres_idx = [i+2 for i in xrange(len(sequence_pdb))]
        for i in xrange(len(sequence_pdb)):
            self.assertEqual(correct_bindres_idx[i], pdbidx_to_bindresidx[i])

        bindres_file = "./test5.bindres"
        sequence_bindres = construct_protein_graph.get_sequence_from_bindres_file(bindres_file)
        pdbidx_to_bindresidx = construct_protein_graph.get_index_map(sequence_bindres, sequence_pdb)
        correct_bindres_idx = [1, 2, 3, 7, 8, 9, 11, 12, 13]
        for i in xrange(len(sequence_pdb)):
            self.assertEqual(correct_bindres_idx[i], pdbidx_to_bindresidx[i])

        bindres_file = "./test6.bindres"
        sequence_bindres = construct_protein_graph.get_sequence_from_bindres_file(bindres_file)
        pdbidx_to_bindresidx = construct_protein_graph.get_index_map(sequence_bindres, sequence_pdb)
        correct_bindres_idx = [i for i in xrange(len(sequence_bindres))]
        for i in xrange(len(sequence_pdb)):
            if i == 0:
                self.assertTrue(not i in pdbidx_to_bindresidx)
            else:
                self.assertEqual(correct_bindres_idx[i-1], pdbidx_to_bindresidx[i])

    def test_construct_protein_graph(self):
        bindres_file = "./test5.bindres"
        sequence_bindres = construct_protein_graph.get_sequence_from_bindres_file(bindres_file)
        pdb_file = "./test.pdb"
        p = PDBParser()
        structure = p.get_structure('test', pdb_file)
        structure_of_chain = structure[0]['A']
        sequence_pdb = construct_protein_graph.get_sequence_from_pdb_structure(structure_of_chain)
        pdbidx_to_bindresidx = construct_protein_graph.get_index_map(sequence_bindres, sequence_pdb)
        atom_type = 'CA'
        delta = 5
        graph = construct_protein_graph.construct_protein_graph(structure_of_chain, pdbidx_to_bindresidx, atom_type, delta)
        correct_edge = [(1, 2), (1, 3), (1, 7), (2, 1), (3, 1), (7, 1), (8, 9), (8, 11), (9, 8), (9, 11), (11, 8), (11, 9)]
        for i in xrange(len(graph)):
            self.assertEqual(correct_edge[i], graph[i])
 

if __name__ == "__main__":
    unittest.main()
