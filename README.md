# MENSAdb

#####################################################################################################################
The features extracted for the database of membrane protein dimer analysis can be replicated through this repository.
#####################################################################################################################

- INSTALLATION REQUIREMENTS: First you will need to run *python setup.py* to install all the dependencies necessary for feature extraction.

- Before feature extraction you should perform a pre-processing of the pdb files. For that you need to:
		- Trim non-transmembrane residues
		- Remove heteroatoms
		- Mutate exotic amino acids
		- Model incomplete structures
		- Dimer extraction from the structure Files
		- Add hydrogens

To see additional details in how to perform data pre-processing please see our review - "Structural Characterization of Membrane Protein Dimers" published in Methods in Molecular Biology - Protein Supersecondary Structures (https://www.springer.com/us/book/9781493991600).

#############################################################
A - Obtain all the features using a single pdb file as input.
#############################################################

- **call_all.py** deploys all the above as well as the need libraries to attain the output files. Run: generate_outputs(your_pdb).joint_call(autodock, autodock_2). The output features are: [class_values, dssp_chain_A, dssp_chain_B, dssp_chain_comp_A, dssp_chain_comp_B, binana_features, jsd_values_A, jsd_values_B]

######################################################################
B - Obtain each feature individually using a single pdb file as input.
######################################################################

- **dssp_features.py** extract the features from a dssp output file. Also requires the corresponding pdb file. To attain the dssp output file use the DSSP executable. For that, first check if you have dssp binary available in your $PATH. This can be done by running: *which dssp* or *which mkdssp*. If you don't have you can install it using CONDA *conda install -c salilab dssp*. To attain the dssp output file try running: *dssp -i [pdb_name.pdb] > [output_name.txt]* or *mkdssp -i [pdb_name.pdb] > [output_name.txt]*. To attain DSSP features, you can run: *python dssp_features.py*.

	  - DSSP index
    - Amino acid number
    - Amino acid code
    - Chain
    - Secondary Structure
    - BP
    - ASA
    - NH-->O_1_relidx
    - O-->NH_1_relidx
    - NH-->O_1_energy
    - O-->NH_1_energy
    - TCO
    - KAPPA
    - Alpha
    - Phi
    - Psi
    - X-CA
    - Y-CA
    - Z-CA

- **features_pssm.py** extract the pssm "jsd" features from psi-blast output file. To retrieve the pssm files needed you will require the psiblast local installation, the non-redundant (nr) database and your input file, with this, run: *psiblast -query [fasta_file.fasta] -evalue 0.001 -num_iterations 2 -db [nr] -outfmt 5 -out pssm_output_name.txt -out_ascii_pssm [output_name.pssm] -num_threads 6"*. Running this step can be very time-consuming, depending on the computer and the protein. To attain PSSM "jsd" features output, you can run: *python features_pssm.py*.

- **process_binana.py** extract the features from the BINding ANAlyser (BINANA - to download go to http://rocce-vm0.ucsd.edu/data/sw/hosted/binana/#download) output file. To attain the BINANA output, you can run: *python binana_1_2_0.py -receptor /path/to/receptor.pdbqt -ligand /path/to/ligand.pdbqt -output_file /path/to/output.pdb*, as stated in the website of this software. To use this command, you will need their binana_1_2_0.py script, as well as the ".pdbqt" input files. To attain the selected features from the BINANA output, you can run: *python process_binana.py*. A single csv will be written for each of the possible features. These features are relate to a dimer, specifically.

	- Below 2.5 Angstrom residues
	- Below 4 Angstrom residues
	- Hydrogen Bonds
	- Hydrophobic contacts
	- Pi-pi bond stack
	- T - stack
	- Cation - pi interaction
	- Salt-bridges

- **generate_class.py** uses vmd to extract the interfacial and surface classification for each residue. Makes use of 5 other scripts that are located on the "mensa_class" folder. To use these scripts, it is required the installation of python based vmd. This can be done with: *conda install -c conda-forge vmd-python*. The whole code can be run with *generate_outputs(input_pdb).joint_call(autodock, autodock_2)*. Check the path list and replace with your locations. The possible classes are:

	- non-interface and non-surface: 0
	- non-interface and surface: 2
	- interface and surface: 3

Please cite:
