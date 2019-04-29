#conda install -c conda-forge vmd-python


import vmd
from vmd import molecule
from vmd import evaltcl
import os
import sys

def process_vmd_output(input_file_name):

    import re
    opened_file = open(input_file_name, "r").readlines()
    no_par_file = opened_file[0].replace("\n","")
    split_file = no_par_file.split("} {")
    res_dict = {}
    for res in split_file:
        processed_res = re.sub('[}{]', '', res)
        res_dict[processed_res.split(" ")[0]] = processed_res.split(" ")[1]
    return res_dict

def detect_chains(input_pdb):

    opened_pdb = open(input_pdb, "r").readlines()
    chains = []
    for value in opened_pdb:
        try:
            if value[21] not in chains:
                chains.append(value[21])
        except:
            continue
    return chains

def add_hydrogens(input_pdb):

    output_string = 'mol new ' + input_pdb + '\n' \
     + 'set pdb_name ' + '"' + input_pdb[0:-4] + '"' + '\n' \
     + 'set protein [atomselect top all]' + '\n' \
     + 'set chains [lsort -unique [$protein get chain]]' + '\n' \
     + 'foreach chain $chains {' + '\n' \
     + 'set sel [atomselect top "protein and chain $chain"]' + '\n' \
     + 'set current_chain_name "${pdb_name}_${chain}.pdb"' + '\n' \
     + "$sel writepdb $current_chain_name" + '\n' \
     + '}' + '\n' \
     + 'package require psfgen' + '\n' \
     + 'topology mensa_class/top_na.inp' + '\n' \
     + 'alias residue HIS HSD' + '\n' \
     + 'alias residue HOH TIP3' + '\n' \
     + 'alias residue ZN ZN2' + '\n' \
     + 'alias atom ILE CD1 CD' + '\n' \
     + 'alias atom HOH O OH2' + '\n' \
     + 'pdbalias residue DG GUA' + '\n' \
     + 'pdbalias residue DC CYT' + '\n' \
     + 'pdbalias residue DA ADE' + '\n' \
     + 'pdbalias residue DT THY' + '\n' \
     + 'foreach bp { GUA CYT ADE THY URA } {' + '\n' \
     + 'pdbalias atom $bp "O5\*" O5' + "'" + '\n' \
     + 'pdbalias atom $bp "C5\*" C5' + "'" + '\n' \
     + 'pdbalias atom $bp "O4\*" O4' + "'" + '\n' \
     + 'pdbalias atom $bp "C4\*" C4' + "'" + '\n' \
     + 'pdbalias atom $bp "C3\*" C3' + "'" + '\n' \
     + 'pdbalias atom $bp "O3\*" O3' + "'" + '\n' \
     + 'pdbalias atom $bp "C2\*" C2' + "'" + '\n' \
     + 'pdbalias atom $bp "O2\*" O2' + "'" + '\n' \
     + 'pdbalias atom $bp "C1\*" C1' + "'" + '\n' \
     + '}' + '\n' \
     + 'foreach monomer $chains {' + '\n' \
     + 'set current_mol_name "${pdb_name}_${monomer}.pdb"'  + '\n' \
     + 'pdbalias residue HIS HSD' + '\n' \
     + 'pdbalias residue HOH TIP3' + '\n' \
     + 'pdbalias residue ZN ZN2' + '\n' \
     + 'pdbalias atom ILE CD1 CD' + '\n' \
     + 'pdbalias atom HOH O OH2' + '\n' \
     + 'pdbalias residue DG GUA' + '\n' \
     + 'pdbalias residue DC CYT' + '\n' \
     + 'pdbalias residue DA ADE' + '\n' \
     + 'pdbalias residue DT THY' + '\n' \
     + 'segment $monomer {' + '\n' \
     + 'pdb $current_mol_name'  + '\n' \
     + '}' + '\n' \
     + 'coordpdb $current_mol_name $monomer' + '\n' \
     + '}' + '\n' \
     + 'guesscoord' + '\n' \
     + 'set output_name "${pdb_name}_complex.pdb"' + '\n' \
     + 'writepdb $output_name' + '\n' \
     + 'mol new "${pdb_name}_complex.pdb"' + '\n' \
     + 'set pdb_name ' + '"${pdb_name}_complex"' + '\n' \
     + 'set protein [atomselect top all]' + '\n' \
     + 'set chains [lsort -unique [$protein get chain]]' + '\n' \
     + 'foreach chain $chains {' + '\n' \
     + 'set sel [atomselect top "protein and chain $chain"]' + '\n' \
     + 'set current_chain_name "${pdb_name}_${chain}.pdb"' + '\n' \
     + "$sel writepdb $current_chain_name" + '\n' \
     + '}' + '\n' \
     + 'quit' + '\n' \
     + 'exit'
    opened_file = open("add_H.tcl.tpl", "w")
    opened_file.write(output_string)
    return output_string

def run_class_pdb(input_pdb_name):
    add_hydrogens(input_pdb_name)
    evaltcl("play add_H.tcl.tpl")

home = os.getcwd()
input_file = home + "/" + sys.argv[1]
run_class_pdb(input_file)
