#conda install -c conda-forge vmd-python


import vmd
from vmd import molecule
from vmd import evaltcl
import os
import sys

home = os.getcwd()
output_folder = home + "/output"

def detect_chains(input_pdb):

    opened_pdb = open(input_pdb, "r").readlines()
    chains = []
    for value in opened_pdb[1:]:
        try:
            if value[21] not in chains:
                chains.append(value[21])
        except:
            continue
    return chains
        
class surface:

    def __init__(self, pdb_path, first = True):

        self.input_pdb_path = pdb_path
        self.chains = detect_chains(self.input_pdb_path)
        if first == True:
            self.current_chain = self.chains[0]
        if first == False:
            self.current_chain = self.chains[1]

    def surface_template(self, chain):

        output_name = output_folder + "/surface_output_" + str(self.input_pdb_path).split("/")[-1][0:-4] + "_" + chain
        output_string = 'mol new ' + self.input_pdb_path + '\n' \
         + 'set allsel [atomselect top "all and chain ' + chain + '"]' + '\n' \
         + 'set chain ' + chain + '\n' \
         + 'set tot_sasa [dict create ARG 241 TRP 259 TYR 229' \
         + ' LYS 211 PHE 218 MET 204 GLN 189 HIS 194 HSD 194 HSE 194 HSP 194 GLU' \
         + ' 183 LEU 180 ILE 182 ASN 158 ASP 151 CYS 140 VAL 160 THR 146 PRO 143 ' \
         + 'SER 122 ALA 113 GLY 85]' + '\n' \
         + 'set residlist [lsort -unique [$allsel get resid]]' + '\n' \
         + 'set surf_list [list]' + '\n' \
         + 'foreach r $residlist {' + '\n' \
         + 'set sel [atomselect top "resid $r and chain $chain"]' + '\n' \
         + 'set temp_rsasa [measure sasa 1.4 $allsel -restrict $sel]' + '\n' \
         + 'set temp_name [lsort -unique [$sel get resname]]' + '\n' \
         + 'set temp_id [lsort -unique [$sel get resid]]' + '\n' \
         + 'set temp_tot [dict get $tot_sasa $temp_name]' + '\n' \
         + 'set rsasa [expr $temp_rsasa/$temp_tot]' + '\n' \
         + 'if {$rsasa > 0.2} {lappend surf_list "$temp_id $temp_name"}' + '\n' \
         + '}' + '\n' \
         + 'set fileId [open ' + output_name + ' w]' + '\n' \
         + 'puts $fileId $surf_list' + '\n' \
         + 'close $fileId' + '\n' \
         + 'quit' + '\n' \
         + 'exit' + '\n'
        return output_string


    def generate_surface_temp(self):

        tcl_name = output_folder + "/surface_output_" + str(self.input_pdb_path).split("/")[-1][0:-4] + "_" + str(self.current_chain) + ".tcl"
        opened_file = open(tcl_name, "w")
        writeable_string = surface(self.input_pdb_path).surface_template(chain = str(self.current_chain))
        opened_file.write(writeable_string)

    def run_read_output_surf(self):
        
        tcl_name = output_folder + "/surface_output_" + str(self.input_pdb_path).split("/")[-1][0:-4] + "_" + str(self.current_chain) + ".tcl"
        run_command = "play "+ tcl_name
        evaltcl(run_command)

def run_class_pdb(input_pdb_name):

    surface(input_pdb_name).generate_surface_temp()
    surface_chain = surface(input_pdb_name).run_read_output_surf()

home = os.getcwd()
print(sys.argv)
input_file = home + "/" + sys.argv[1]
run_class_pdb(input_file)



