mol new /Users/pedrofilipe/Documents/GitHub/mensadb/example_files/1a0t_PQ.pdb
set pdb_name "/Users/pedrofilipe/Documents/GitHub/mensadb/example_files/1a0t_PQ"
set protein [atomselect top all]
set chains [lsort -unique [$protein get chain]]
foreach chain $chains {
set sel [atomselect top "protein and chain $chain"]
set current_chain_name "${pdb_name}_${chain}.pdb"
$sel writepdb $current_chain_name
}
package require psfgen
topology top_na.inp
alias residue HIS HSD
alias residue HOH TIP3
alias residue ZN ZN2
alias atom ILE CD1 CD
alias atom HOH O OH2
pdbalias residue DG GUA
pdbalias residue DC CYT
pdbalias residue DA ADE
pdbalias residue DT THY
foreach bp { GUA CYT ADE THY URA } {
pdbalias atom $bp "O5\*" O5'
pdbalias atom $bp "C5\*" C5'
pdbalias atom $bp "O4\*" O4'
pdbalias atom $bp "C4\*" C4'
pdbalias atom $bp "C3\*" C3'
pdbalias atom $bp "O3\*" O3'
pdbalias atom $bp "C2\*" C2'
pdbalias atom $bp "O2\*" O2'
pdbalias atom $bp "C1\*" C1'
}
foreach monomer $chains {
set current_mol_name "${pdb_name}_${monomer}.pdb"
pdbalias residue HIS HSD
pdbalias residue HOH TIP3
pdbalias residue ZN ZN2
pdbalias atom ILE CD1 CD
pdbalias atom HOH O OH2
pdbalias residue DG GUA
pdbalias residue DC CYT
pdbalias residue DA ADE
pdbalias residue DT THY
segment $monomer {
pdb $current_mol_name
}
coordpdb $current_mol_name $monomer
}
guesscoord
set output_name "${pdb_name}_complex.pdb"
writepdb $output_name
quit
exit