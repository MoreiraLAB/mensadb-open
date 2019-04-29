#!/usr/bin/env python

"""
Extract jsd value from pssm output file
"""

import csv
import os
import math
import pandas as pd
import numpy as np
import string
import dit
from dit.divergences import jensen_shannon_divergence
from numpy import apply_along_axis

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "MENSA"

def open_pssm_file(pssm_name):

    """
    Process pssm table
    """
    opened_file = open(pssm_name, "r").readlines()
    start = False
    stop = False
    for row in opened_file:
        new_row = row.split()
        if stop == False:
            try:
                if new_row[1] == "Lambda":
                    stop = True
            except:
                continue
        if stop == True:
            break
        if start == True and new_row[0] not in string.ascii_uppercase:
            yield new_row
        if start == False:
            if new_row[0] == "Last":
                start = True


def retrieve_col(pssm_file, target_col = -2):

    """
    Retrieve the jsd value from the pssm table
    """
    processed_file = open_pssm_file(pssm_file)

    output_dict = {}
    for row in processed_file:
        output_dict[row[0]] = row[target_row]
    return output_dict

def prepare_dataset(input_list):

    """
    Prepare dataset excluding "X" residues
    """
    res_list = []
    probs_list = []
    for row in input_list:
        if row[1] != "X":
            res_list.append(row[1])
            probs_list.append([float(i) for i in row[22:42]])
    return res_list, probs_list

def scale_table(input_table):

    """
    Scale probability table and ensure that the row sum is 100%
    """
    output_table = []
    for row in input_table:
        percent_row = [float(i/100) for i in row]
        total = sum(percent_row)
        if total != 1:
            percent_row = apply_along_axis(lambda x: x / total, 0, percent_row).tolist()
        output_table.append(percent_row)
    return output_table

def calculate_jsd(input_values, input_probabilities):

    """
    Calculated Jensen-Shannon Divergence upon the table
    """
    processed_probabilities = scale_table(input_probabilities)
    x = dit.ScalarDistribution(amino_list, input_values, sample_space = amino_list, sort = True)
    jsd_values = []
    for row in processed_probabilities:
        y = dit.ScalarDistribution(amino_list, row, sample_space = amino_list, sort = True)
        jsd_values.append(jensen_shannon_divergence([x,y]))
    return jsd_values


def generate_dict(input_list):

    """
    Generate dictionary from list. list values are dictionary keys
    """
    empty_dict = {}
    for row in input_list:
        empty_dict[row] = 0
    return empty_dict

def calc_prob(input_dict, input_len):

    """
    Calculate amino acid probability occurence
    """
    output_list = []
    for res in amino_list:
        output_list.append(float(input_dict[res])/float(input_len))
    return output_list

def probability_generator(input_ordered_res):

    """
    Export probability
    """
    amino_dict = generate_dict(amino_list)
    count = 0
    for res in input_ordered_res:
        count += 1
        amino_dict[res] += 1
    return calc_prob(amino_dict, count)

def joint_call(input_pssm_name):

    """
    Call the above functions
    """
    processed_pssm = list(open_pssm_file(input_pssm_name))
    residues, probabilities = prepare_dataset(processed_pssm)
    overall_probs = probability_generator(residues)
    jsd_output = calculate_jsd(overall_probs, probabilities)
    return jsd_output


amino_list = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

