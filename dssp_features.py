#!/usr/bin/env python

"""
Extract DSSP features directly from a DSSP output file with this script.
The DSSP output file can be attained by running the dssp executable file with
the target pdb file as input.
You will need both the pdb file and the dssp output file
"""

import os
import math
import Bio
from Bio.PDB import *

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "MENSAdb"

def pdb_parser(input_pdb):

    """
    Parses pdb from .pdb file
    """
    parser = PDBParser()
    pdb_name = input_pdb[0:-4]
    structure = parser.get_structure(pdb_name, input_pdb)
    return structure, pdb_name

def DSSP_gap_dictionary():

    """
    DSSP files are not easily separable by default character.
    As such, it can be accurately split by space values, which are stated below
    """
    dssp_gaps = {"0":[0,5], "1":[5,10], "2":[10,12], "3": [12,14], "4":[14,22], "5":[22,33],
                 "6":[34,38],"7":[38,50], "8":[50,61], "9":[61,72], "10":[72,83], "11":[83,91],
                  "12":[91,97], "13":[97,103], "14":[103,109], "15":[109,115], "16":[115,122],
                   "17":[122,129], "18":[129, 136], "19":[136, 150]}
    return dssp_gaps

def round_number(input_number, round_to = 2):

    """
    Simple function to round a number to "round_to"
    decimal houses. Change "round_to" to modify the default number of decimal houses
    """
    if len(str(input_number)) == 0:
        input_number = 0.0
    return round(float(input_number), round_to)

def run_for_chain(input_file, input_target_chain, feature_number):

    feature_gaps = DSSP_gap_dictionary()
    feature_residues = {}
    residues_count = 0
    useful = False
    to_break = [7,8,9,10]
    for row in input_file:
        if useful == True:
            if row[feature_gaps["2"][0]:feature_gaps["2"][-1]].replace(" ","") == input_target_chain:
                residues_count = residues_count + 1
                if feature_number in to_break:
                    feature_to_store = row[feature_gaps[str(feature_number)][0]:feature_gaps[str(feature_number)][1]].replace(" ","").split(",")
                    average_float = (float(feature_to_store[0]) + float(feature_to_store[1])) / float(2)
                    feature_value = round_number(average_float)
                    feature_residues[residues_count] = feature_value
                else:
                    feature_value = round_number(row[feature_gaps[str(feature_number)][0]:feature_gaps[str(feature_number)][1]].replace(" ",""))
                    feature_residues[residues_count] = feature_value
        if row[feature_gaps["0"][0]:feature_gaps["0"][-1]].replace(" ","") == "#":
            useful = True
    return feature_residues
        

def DSSP_features(input_pdb, feature_number, dssp_termination = "_dssp.txt", pdb_chain = None):

    """
    Retrieves the features 0-13 described bellow, from bioython structure object
    - 0            DSSP index 
    - 1            Amino acid number
    - 2            Amino acid code
    - 3            Chain 
    - 4            Secondary Structure
    - 5            BP
    - 6            ASA
    - 7            NH-->O_1_relidx 
    - 8            O-->NH_1_relidx 
    - 9            NH-->O_1_energy 
    - 10           O-->NH_1_energy 
    - 11           TCO 
    - 12           KAPPA 
    - 13           Alpha 
    - 14           Phi 
    - 15           Psi 
    - 16           X-CA
    - 17           Y-CA
    - 18           Z-CA

    Change the "dssp_termination" argument if it does not correspond to your own
    """
    structure = pdb_parser(input_pdb)[0]
    dssp_name = "output/" + input_pdb[0:-4] + dssp_termination
    opened_file = open(dssp_name, "r").readlines()
    chain_SS_sequences = []
    useful = False
    
    output_features = run_for_chain(opened_file, pdb_chain , feature_number)
    return output_features


