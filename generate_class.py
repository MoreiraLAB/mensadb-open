#!/usr/bin/env python

"""
Use this script to assess the interfacial classification of protein residues
"""

import os
import sys

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "MENSAdb"

def process_vmd_output(input_file_name):

    """
    Process the output of vmd
    """
    import re
    opened_file = open(input_file_name, "r").readlines()
    no_par_file = opened_file[0].replace("\n","")
    split_file = no_par_file.split("} {")
    res_dict = {}
    for res in split_file:
        processed_res = re.sub('[}{]', '', res)
        res_dict[processed_res.split(" ")[0]] = processed_res.split(" ")[1]           
    return res_dict

def run_generate_classifier(input_pdb):

    """
    Run the remaining scripts to attain the necessary files
    """
    hydrogens_command = "python mensa_class/hydrogens.py " + input_pdb
    os.system(hydrogens_command)
    surface_A_command = "python mensa_class/surface_A.py " + input_pdb[0:-4] + "_complex.pdb"
    os.system(surface_A_command)
    surface_B_command = "python mensa_class/surface_B.py " + input_pdb[0:-4] + "_complex.pdb"
    os.system(surface_B_command)
    interface_command_A = "python mensa_class/interface_A.py " + input_pdb[0:-4] + "_complex.pdb"
    os.system(interface_command_A)
    interface_command_B = "python mensa_class/interface_B.py " + input_pdb[0:-4] + "_complex.pdb"
    os.system(interface_command_B)

def process_data(input_pdb, chains):

    """
    Process the VMD output files
    """
    interface_A = process_vmd_output("output/interface_output_A")
    interface_B = process_vmd_output("output/interface_output_B")
    surface_A_name = "output/surface_output_" + input_pdb[0:-4] + "_complex_" + chains[0]
    surface_A = process_vmd_output(surface_A_name)
    surface_B_name = "output/surface_output_" + input_pdb[0:-4] + "_complex_" + chains[1]
    surface_B = process_vmd_output(surface_B_name)
    return interface_A, surface_A, interface_B, surface_B

def process_pdb(input_pdb):

    """
    Process the pdb in a position based manner
    """
    processed_pdb = []
    for row in open(input_pdb):
        try:
            identifier = row[0:4]
            res_name = row[17:20].replace(" ","")
            res_number = row[22:26].replace(" ","")
            chain = row[21].replace(" ","")
            temp = float(row[60:66].replace(" ",""))
            processed_pdb.append([res_number, res_name, chain])
        except:
            continue
    return processed_pdb
            
def assess_prot(input_processed_pdb, input_features_list, chains_names):

    """
    Check which of the residues are:
    - interfacial (3)
    - non - interfacial (2)
    - non - surface (0)
    """
    output_dict_A = {}
    output_dict_B = {}
    for row in input_processed_pdb:
        res_number = row[0]
        chain = row[2]
        if chain == chains_names[0]:
            if (res_number in input_features_list[0].keys()) and (res_number in input_features_list[1].keys()):
                output_dict_A[res_number] = "3"
            if (res_number not in input_features_list[0].keys()) and (res_number in input_features_list[1].keys()):
                output_dict_A[res_number] = "2"
            if (res_number not in input_features_list[0].keys()) and (res_number not in input_features_list[1].keys()):
                output_dict_A[res_number] = "0"
            if (res_number in input_features_list[0].keys()) and (res_number not in input_features_list[1].keys()):
                output_dict_A[res_number] = "0"
        if chain == chains_names[1]:
            if (res_number in input_features_list[2].keys()) and (res_number in input_features_list[3].keys()):
                output_dict_B[res_number] = "3"
            if (res_number not in input_features_list[2].keys()) and (res_number in input_features_list[3].keys()):
                output_dict_B[res_number] = "2"
            if (res_number not in input_features_list[2].keys()) and (res_number not in input_features_list[3].keys()):
                output_dict_B[res_number] = "0"
            if (res_number in input_features_list[2].keys()) and (res_number not in input_features_list[3].keys()):
                output_dict_B[res_number] = "0"
    return output_dict_A, output_dict_B


def generate_output(input_file):

    """
    Run all the functions for the pdb
    """
    chains_names = input_file.split("_")[1][0:2]
    run_generate_classifier(input_file)
    processed_info = process_data(input_file, chains_names)
    processed_pdb = process_pdb(input_file)
    class_results = assess_prot(processed_pdb, processed_info, chains_names)
    return class_results

"""
Input your file as input file, the program will then use VMD to:
- add hydrogens
- calculate interface
- calculate surface
- compare to the input pdb and retrieve full dictionaries
"""
