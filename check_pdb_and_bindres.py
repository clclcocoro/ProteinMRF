#!/usr/bin/env python

import os
import re
from Bio.PDB import PDBParser

bindres_dir = '/Users/clclcocoro/work/atp_bindingPred/bindres/atp227'
pdb_dir = '/Users/clclcocoro/work/atp_bindingPred/pdb/atp227'
p = PDBParser()
for bindres_file in os.listdir(bindres_dir):
    pdbid_chain = re.match(r'(.+)\.bindres', bindres_file).group(1)
    pdbid = pdbid_chain[:4]
    chain = pdbid_chain[4]
    seqlen_bindres = 0
    with open('{}/{}'.format(bindres_dir, bindres_file)) as fp_b:
        for line in fp_b:
            if line != '' and not re.match('>', line):
                seqlen_bindres = len(line.rstrip())
    p = PDBParser()
    s = p.get_structure(pdbid, '{}/pdb{}.ent'.format(pdb_dir, pdbid.lower()))
    seqlen_pdb = 0
    for res in s[0][chain].get_residues():
        if res.id[0] == ' ': # ATOM Record, not HETATM
            seqlen_pdb += 1
    if seqlen_bindres != seqlen_pdb:
        print "Doesn't match sequence length pdbid_chain {} seqlen_bindres {} seqlen_pdb {}".format(pdbid_chain, seqlen_bindres, seqlen_pdb)
    else:
        print "Match         pdbid_chain {}".format(pdbid_chain)
