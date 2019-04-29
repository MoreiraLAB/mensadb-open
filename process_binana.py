#!/usr/bin/env python

"""
Extract information from BINANA output files. You will require the pdb output of BINANA
"""

import os
import csv

home = os.getcwd()
__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "MENSAdb"

def open_file(file_name):

	"""
	Open and process the input pdb file
	"""
	opened_file = open(file_name, "r").readlines()
	processed_file = []
	for lines in opened_file:
		new_line = lines.replace("\n","")
		new_line = lines.replace("/","")
		new_line = new_line.split()
		processed_file.append(new_line)
	return processed_file

def find_sections(input_processed_file, header_sections):


	"""
	Locate the feature sections on the input processed_file
	"""
	count = -1
	entered_section = False
	empty_sections = []
	non_empty_sections = {}
	for line in input_processed_file:
		try:
			if (entered_section == True) and (len(line) == 1) and (line[0] == "REMARK"):
				entered_section = False
			if len(line) != 0:
				if line[0] == "Raw":
					entered_section = True
					if count < 7:
						count += 1
			if len(line) != 0 and entered_section == True:
				if header_sections[count] not in non_empty_sections:
					non_empty_sections[header_sections[count]] = []
					non_empty_sections[header_sections[count]].append(line)
				else:
					non_empty_sections[header_sections[count]].append(line)
		except:
			if header_sections[count] not in empty_sections:
				empty_sections.append(header_sections[count])
		if (len(line) != 0) and (len(line) > 1):
			if (line[1] == "The") and (count == 7):
				break
	return non_empty_sections


def read_feat_dict(input_dict, feature_name, starting_point, chain_name_position, original_chains = "AB", SB = False):

	import re
	"""
	Process the file to yield the information for both chains
	"""
	res_dict_A = {}
	res_dict_B = {}
	reverse_start = chain_name_position
	reverse_start_string = starting_point
	if SB == True:
		reverse_start = chain_name_position - 1
		reverse_start_string = starting_point - 1
	for i in range(0,len(input_dict[feature_name])):
		if input_dict[feature_name][i][0] == 'REMARK':
			res_list_1 = input_dict[feature_name][i][1].strip("'")
			if res_list_1[chain_name_position] == original_chains[1]:
				count_1 = 0
				for breakpoint in res_list_1:
					if breakpoint == ")":
						break
					count_1 += 1
				res_number_1 = int(res_list_1[starting_point:count_1])
				if res_number_1 not in res_dict_B:
					res_dict_B[res_number_1] = 1
				else:
					res_dict_B[res_number_1] = res_dict_B[res_number_1] + 1
			if res_list_1[chain_name_position] == original_chains[0]:
				count_1 = 0
				for breakpoint in res_list_1:
					if breakpoint == ")":
						break
					count_1 += 1
				res_number_1 = int(res_list_1[starting_point:count_1])
				if res_number_1 not in res_dict_A:
					res_dict_A[res_number_1] = 1
				else:
					res_dict_A[res_number_1] = res_dict_A[res_number_1] + 1
			res_list_2 = input_dict[feature_name][i][-1].strip("'")
			if res_list_2[reverse_start] == original_chains[1]:
				count_2 = 0
				for breakpoint in res_list_2:
					if breakpoint == ")":
						break
					count_2 += 1
				res_number_2 = int(re.sub('[)()]','',str(res_list_2[reverse_start_string:count_2])))
				if res_number_2 not in res_dict_B:
					res_dict_B[res_number_2] = 1
				else:
					res_dict_B[res_number_2] = res_dict_B[res_number_2] + 1
			if res_list_2[reverse_start] == original_chains[0]:
				count_2 = 0
				for breakpoint in res_list_2:
					if breakpoint == ")":
						break
					count_2 += 1
				res_number_2 = int(re.sub('[)()]','',res_list_2[reverse_start_string:count_2]))
				if res_number_2 not in res_dict_A:
					res_dict_A[res_number_2] = 1
				else:
					res_dict_A[res_number_2] = res_dict_A[res_number_2] + 1

	return res_dict_A, res_dict_B

def retrieve_features(input_file_name):

	"""
	Evaluate file for all possible sections.
	The different features list designates features which have different processing
	Call "write_csv_from_dict" to write a csv file for each feature
	Uncomment the "write_csv_from_dict" calls to write a csv file for each feature
	"""
	possible_sections = ["below_2.5A", "below_4A", "HB", "hydrophobic_contacts", "pipi_stack", "T_stack", "cat_pi","SB"]
	different_features = ["pipi_stack", "T_stack", "cat_pi","SB"]
	opened_file = open_file(input_file_name)
	feature_sections = find_sections(opened_file, possible_sections)
	chains = input_file_name.split("_")[2].split(".")[0]
	complex_name = input_file_name.split("_")[1]
	features_output = {}
	for feature in possible_sections:
		if feature != "SB":
			dict_A, dict_B = read_feat_dict(feature_sections, feature, 6, 0, original_chains = chains)
			features_output[feature] = [dict_A, dict_B]
		if feature in different_features:
			dict_A, dict_B = read_feat_dict(feature_sections, feature, 7, 1, original_chains = chains, SB = True)
			features_output[feature] = [dict_A, dict_B]
	return features_output
