#!/usr/bin/env python

"""
This script takes as input a .pdb file and returns all information listed
Usage: generate_outputs(input_pdb).joint_call(autodock, autodock_2)
The output is a list of outputs:
[class_values, dssp_chain_A, dssp_chain_B, dssp_chain_comp_A, dssp_chain_comp_B, binana_features, jsd_values_A, jsd_values_B]
"""

import os
import sys
import Bio
from Bio.PDB import *
from Bio.PDB.Polypeptide import three_to_one

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "MENSAdb"

def move_to_folder(folder = "output"):

    """
    Move remaining files to folder
    """
    for files in os.listdir(os.getcwd()):
        if files.endswith(".tcl") or files.endswith(".pdb") or files.endswith(".fasta") or files.endswith(".tpl"):
            new_file = folder + "/" + files
            os.rename(files, new_file)

def pdb_to_fasta(pdb_input):

    """
    Write a fasta file for the processed pdb
    """
    p = PDBParser(PERMISSIVE=1)
    structure = p.get_structure(pdb_input, pdb_input)
    file_name = pdb_input[0:-4] + ".fasta"
    fasta_file = open(file_name, 'w')
    for model in structure:
        for chain in model:
            seq = list()
            chainID = chain.get_id()

            for residue in chain:
                if is_aa(residue.get_resname(), standard=True):
                    seq.append(three_to_one(residue.get_resname()))
                else:
                    seq.append("X")
            chain_line = ">Chain_" + chainID + "\n" + str("".join(seq)) + "\n" + "\n"
            fasta_file.write(chain_line)

    fasta_file.close()

def move_to_folder(folder = "output"):

    """
    Move remaining files to folder
    """

    for files in os.listdir(os.getcwd()):
        if files.endswith(".tcl") or files.endswith(".pdb") or files.endsiwth(".fasta"):
            new_file = folder + "/" + files
            os.rename(files, new_file)

class generate_outputs:

    def __init__(self, input_pdb, autodock, autodock_2, psiblast_path, nr_path):

        """
        Define commonly used variables
        """
        self.home = os.getcwd()
        self.autodock = autodock #"MGL=/opt/mgltools_x86_64Linux2_1.5.6/bin\n"
        self.autodock_2 = autodock_2 #"ADT=/opt/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24\n"
        self.psiblast_path = psiblast_path #"/services/psiblast"
        self.nr_path = nr_path #"nr"
        self.pdb_path = input_pdb
        self.pdb_name = input_pdb[0:-4]
        self.complex_name = input_pdb[0:-4] + "_complex"
        self.complex_name_A = self.complex_name + "_A.pdb"
        self.complex_name_B = self.complex_name + "_B.pdb"

    def call_dssp(self):

        """
        Call dssp for both chains
        """
        import dssp_features
        output_name_A = "output/" + self.complex_name + "_A_dssp.txt"
        dssp_command_A = "dssp -i " + self.complex_name_A + " > " + output_name_A
        os.system(dssp_command_A)
        dssp_chain_A = dssp_features.DSSP_features(self.complex_name_A, 6, pdb_chain = "A")

        output_name_B = "output/" + self.complex_name + "_B_dssp.txt"
        dssp_command_B = "dssp -i " + self.complex_name_B + " > " + output_name_B
        os.system(dssp_command_B)
        dssp_chain_B = dssp_features.DSSP_features(self.complex_name_B, 6, pdb_chain = "B")

        complex_pdb = self.complex_name + ".pdb"
        output_name = "output/" + self.complex_name + "_dssp.txt"
        dssp_command = "dssp -i " + complex_pdb + " > " + output_name
        os.system(dssp_command)
        dssp_chain_comp_A = dssp_features.DSSP_features(complex_pdb, 6, pdb_chain = "A")
        dssp_chain_comp_B = dssp_features.DSSP_features(complex_pdb, 6, pdb_chain = "B")
        return dssp_chain_A, dssp_chain_B, dssp_chain_comp_A, dssp_chain_comp_B

    def write_binana_sh(self, autodock_path, autodock_path_2, command_chain_A, command_chain_B):

        """
        Generate pdbqts from the pdb files
        """
        opened_sh = open("pdbqt_generator.sh","w")
        opened_sh.write(autodock_path)
        opened_sh.write(autodock_path_2)
        opened_sh.write(command_chain_A)
        opened_sh.write(command_chain_B)
        opened_sh.close()


    def call_binana(self, autodock_path, autodock_path_2):

        """
        Prepare pdbqt generationg and run binana
        """
        import process_binana
        autodock_chain_A = "$MGL/pythonsh $ADT/prepare_ligand4.py -l " + self.complex_name_A + " -A hydrogens -o output/" + self.complex_name + "_A.pdbqt\n"
        autodock_chain_B = "$MGL/pythonsh $ADT/prepare_receptor4.py -r " + self.complex_name_B + " -A hydrogens -o output/" + self.complex_name + "_B.pdbqt\n"
        self.write_binana_sh(autodock_path, autodock_path_2, autodock_chain_A, autodock_chain_B)
        os.system("sh pdbqt_generator.sh")
        binana_output = "binana_" + self.pdb_path
        binana_command = "python services/binana_1_2_0.py -receptor output/" + self.pdb_name + "_complex_B.pdbqt -ligand output/" \
                            + self.pdb_name + "_complex_A.pdbqt -output_file " + binana_output
        os.system(binana_command)
        binana_output = process_binana.retrieve_features(binana_output)
        return binana_output

    def call_pssm(self):

        """
        Call psiblast, modify psiblast_path and nr_path according to your own paths
        """
        pdb_to_fasta(self.complex_name_A)
        psiblast_string_A = self.psiblast_path + " -query " + self.pdb_name + "_complex_A.fasta" \
                         + " -evalue 0.001 -num_iterations 2 -db " + self.nr_path + " -outfmt 5 -out " \
                         + "output/pssm_" + self.pdb_name + "_A.txt -out_ascii_pssm output/pssm_" + self.pdb_name + "_A.pssm -num_threads 6"
        os.system(psiblast_string_A)

        pdb_to_fasta(self.complex_name_B)
        psiblast_string_B = self.psiblast_path + " -query " + self.pdb_name + "_complex_B.fasta" \
                         + " -evalue 0.001 -num_iterations 2 -db " + self.nr_path + " -outfmt 5 -out " \
                         + "output/pssm_" + self.pdb_name + "_B.txt -out_ascii_pssm output/pssm_" + self.pdb_name + "_B.pssm -num_threads 6"
        os.system(psiblast_string_B)

    def call_class(self):

        """
        Retrieve class values
        """
        import generate_class
        return generate_class.generate_output(self.pdb_path)

    def joint_call(self, autodock, autodock_2):

        """
        Call class, dssp, binana, pssm (jsd) functions and retrieve its values
        Output_order = class, dssp chain A, dssp chain B, binana features, jsd values chain A, jsd values chain B
        """
        import features_pssm
        class_values = self.call_class()
        dssp_chain_A, dssp_chain_B, dssp_chain_comp_A, dssp_chain_comp_B = self.call_dssp()
        binana_features = self.call_binana(self.autodock, self.autodock_2)
        self.call_pssm()
        pssm_output_A = "output/pssm_" + self.pdb_name + "_A.pssm"
        jsd_values_A = features_pssm.joint_call(pssm_output_A)
        pssm_output_B = "output/pssm_" + self.pdb_name + "_B.pssm"
        jsd_values_B = features_pssm.joint_call(pssm_output_B)
        move_to_folder()
        return class_values, dssp_chain_A, dssp_chain_B, dssp_chain_comp_A, dssp_chain_comp_B, binana_features, jsd_values_A, jsd_values_B

'''
home = os.getcwd()
autodock = "MGL=/opt/mgltools_x86_64Linux2_1.5.6/bin\n"
autodock_2 = "ADT=/opt/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24\n"
psiblast_path = "services/psiblast"
nr_path = "nr"
'''
