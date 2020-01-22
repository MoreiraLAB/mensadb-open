from mensadb_fetcher import generate_outputs
import pandas as pd
import os
from statistics import mean
import sys
from exampledicr import getdicio
import numpy as np

def convDic(A, last_aa):
    '''
    Converts aminoacid numeration from the numeration contained within
    the feature outputs to a relative numeration includin the whole
    protein instead of a numeration by chain.
    '''
    B = {}
    n = 1
    for entry in A.keys():
        if entry == np.number:
            B[n + int(last_aa)] = A[int(entry)]
        else:
            B[n + int(last_aa)] = A[entry]
        n += 1
    return B

def turnRelative(A, chain, posA, posB):
    '''
    Converts aminoacid numeration from the numeration contained within
    the pbd file to a relative numeration according to the interface
    defined in the literature.
    '''
    B = {}
    if chain == 'A':
        for q in A.keys():
            B[posA[int(q)]] = A[q]
    else:
        for q in A.keys():
            B[posB[int(q)]] = A[q]
    return B

class Run:
    '''
    Runs all dependencies of mensadb and produces a csv file combining
    all the information. This csv file also includes a scraped version
    of the pdb file of the protein in analysis.
    '''

    def __init__(self, pdbid, chains, autodock_1, autodock_2, database = 'nr'):
    '''
    Defines the starting attributes of the run class.
    When calling Run, a pdb code (corresponding to the file) and the
    corresponding chains should be given.
    '''
        self.pdbid      = pdbid
        self.chains     = chains
        self.file       = self.pdbid + '_' + self.chains + '.pdb'
        self.autodock_1 = autodock_1
        self.autodock_2 = autodock_2
        self.psiblast   = 'psiblast'
        self.nr         = database

    def runer(self):
        '''
        Runs all the dependencies and loads them as class attributes.
        Obtained characteristics:
        - aminoacid classification
        - monomer and complex solvent accessible surface areas
        - interactions at each aminoacid
        - conservation
        '''
        protein = generate_outputs(self.file, self.autodock_1, self.autodock_2, self.psiblast, self.nr)
        self.all_features = protein.joint_call()

    def getSequence(self):
        '''
        Reads the corresponding pdb file, obtaining a relationship
        between all aminoacids and its characteristics.
        This method also scraped raw b-factor data.
        Obtained characteristics:
        - aminoacid numeration
        - chains
        - b-factor
        '''
        self.runer()
        file = open(self.file, 'r').readlines()

        self.seq      = {}
        self.mon      = {}
        self.orchains = {}
        self.position = {}
        self.bfactor  = {}

        n = 1

        BF      = []
        b = file[0][22:28].replace(' ', '')

        for line in file:

            typ     = line[0:6].replace(' ', '')
            atom    = line[11:16].replace(' ', '')
            aa      = line[17:20].replace(' ', '')
            chain   = line[20:22].replace(' ', '')
            pos     = line[22:28].replace(' ', '')
            bf      = line[60:66].replace(' ', '')

            if pos == b and typ == 'ATOM':
                BF.append(float(bf))
            elif pos != b and typ == 'ATOM':
                self.bfactor[n-1] = mean(BF)
                BF = []
                BF.append(float(bf))
                b = pos
            else:
                continue

            if typ == 'ATOM' and atom == 'CA':
                self.seq[n] = aa
                self.mon[n] = chain
                self.position[n] = pos
                if chain == 'A' or chain == self.chains[0]:
                    self.orchains[n] = self.chains[0]
                elif chain == 'B' or chain == self.chains[1]:
                    self.orchains[n] = self.chains[1]
                n += 1

            else:
                continue

        self.bfactor[n-1] = mean(BF)

        self.last_aa = 0

        for entry in self.mon.keys():
            if self.mon[entry] == 'A' or self.mon[entry] == self.chains[0]:
                self.last_aa = entry
            else:
                break

    def getFeatures(self):
        '''
        Correlates all features and their aminoacid position.
        '''
        self.getSequence()

        self.revposA = {}
        self.revposB = {}

        for w in self.position.keys():
            if w <= self.last_aa:
                self.revposA[int(self.position[w])] = w
            else:
                self.revposB[int(self.position[w])] = w

        self.dssp_asa_A         = self.all_features[1]
        dssp_asa_pB             = self.all_features[2]

        self.dssp_asa_B         = {}

        for entryB in dssp_asa_pB.keys():
            self.dssp_asa_B[entryB + self.last_aa] = dssp_asa_pB[entryB]

        self.dssp_asa_mon       = self.dssp_asa_A
        self.dssp_asa_mon.update(self.dssp_asa_B)
        self.dssp_asa_complex_A = self.all_features[3]
        dssp_asa_complex_pB = self.all_features[4]
        self.dssp_asa_complex_B = {}
        for entrypB in dssp_asa_complex_pB.keys():
            self.dssp_asa_complex_B[entrypB + self.last_aa] = dssp_asa_complex_pB[entrypB]
        self.dssp_asa_complex   = self.dssp_asa_complex_A
        self.dssp_asa_complex.update(self.dssp_asa_complex_B)

        self.class_a        = turnRelative(self.all_features[0][0], 'A', self.revposA, self.revposB)
        self.classy         = {}

        vmd_class_a         = convDic(self.all_features[0][1], self.last_aa)
        self.class_a.update(vmd_class_a)

        self.binana_b25     = turnRelative(self.all_features[5]['below_2.5A'][0], 'A', self.revposA, self.revposB)
        binana_b251         = turnRelative(self.all_features[5]['below_2.5A'][1], 'B', self.revposA, self.revposB)

        self.binana_b4      = turnRelative(self.all_features[5]['below_4A'][0], 'A', self.revposA, self.revposB)
        binana_b41          = turnRelative(self.all_features[5]['below_4A'][1], 'B', self.revposA, self.revposB)

        self.binana_hydroge = turnRelative(self.all_features[5]['HB'][0], 'A', self.revposA, self.revposB)
        binana_hydroge1     = turnRelative(self.all_features[5]['HB'][1], 'B', self.revposA, self.revposB)

        self.binana_hydroph = turnRelative(self.all_features[5]['hydrophobic_contacts'][0], 'A', self.revposA, self.revposB)
        binana_hydroph1     = turnRelative(self.all_features[5]['hydrophobic_contacts'][1], 'B', self.revposA, self.revposB)

        self.binana_pis     = turnRelative(self.all_features[5]['pipi_stack'][0], 'A', self.revposA, self.revposB)
        binana_pis1         = turnRelative(self.all_features[5]['pipi_stack'][1], 'B', self.revposA, self.revposB)

        self.binana_ts      = turnRelative(self.all_features[5]['T_stack'][0], 'A', self.revposA, self.revposB)
        binana_ts1          = turnRelative(self.all_features[5]['T_stack'][1], 'B', self.revposA, self.revposB)

        self.binana_salt    = turnRelative(self.all_features[5]['SB'][0], 'A', self.revposA, self.revposB)
        binana_salt1        = turnRelative(self.all_features[5]['SB'][1], 'B', self.revposA, self.revposB)

        self.binana_cation_pi= turnRelative(self.all_features[5]['cat_pi'][0], 'A', self.revposA, self.revposB)
        binana_cation_pi1    = turnRelative(self.all_features[5]['SB'][1], 'B', self.revposA, self.revposB)

        self.binana_b25.update(binana_b251)
        self.binana_b4.update(binana_b41)
        self.binana_hydroge.update(binana_hydroge1)
        self.binana_hydroph.update(binana_hydroph1)
        self.binana_pis.update(binana_pis1)
        self.binana_ts.update(binana_ts1)
        self.binana_salt.update(binana_salt1)
        self.binana_cation_pi.update(binana_cation_pi1)

        print(self.binana_b25, self.binana_b4, self.binana_hydroge, sep = '\n')

        self.jsd = self.all_features[6] + self.all_features[7]

    def buildDataFrame(self):
        '''
        Builds a dataframe using the data acquired so far.
        '''
        self.getFeatures()

        available_data = [
            self.dssp_asa_mon,
            self.dssp_asa_complex,
            self.class_a,
            self.binana_b25,
            self.binana_b4,
            self.binana_hydroge,
            self.binana_hydroph,
            self.binana_pis,
            self.binana_ts,
            self.binana_salt,
            self.binana_cation_pi
            ]

        for feat in available_data:
            true_calc = feat.keys()

            for posit in self.position.keys():
                if posit not in true_calc:
                    feat[posit] = float(0)
                else:
                    continue

        header = ['pdb_id',
            'original_chains',
            'chain_position',
            'residue_name',
            'residue_number',
            'classifier',
            'dssp_rel_asa_mon',
            'dssp_rel_asa',
            'dssp_delta_asa',
            'dssp_rel2_asa',
            'b_factor_mean',
            'binana_dist_2_5',
            'binana_dist_4_0',
            'binana_h_bonds',
            'binana_hydrophobic',
            'binana_s_bridges',
            'binana_pi_pi',
            'binana_t_stacking',
            'binana_cation_pi',
            'pssm_jsd']

        self.dataframe = pd.DataFrame(columns = header,
            index = list(self.position.keys()))

        for num in self.position.keys():
            print(num)
            fb_env_mean = []

            for wind_ele in range(num - 5, num + 6):
                try:
                    fb_env_mean.append(self.bfactor[wind_ele])
                except:
                    continue

            if num <= self.last_aa:
                if self.class_a[num] == '3':
                    if self.dssp_asa_mon[num] != 0.0:
                        row = {'pdb_id': self.pdbid,
                            'original_chains': self.chains,
                            'chain_position': 'A',
                            'residue_name': self.seq[num],
                            'residue_number': self.position[num],
                            'classifier': self.class_a[num],
                            'dssp_rel_asa_mon': self.dssp_asa_mon[num],
                            'dssp_rel_asa': self.dssp_asa_complex[num],
                            'dssp_delta_asa': float(self.dssp_asa_mon[num])-float(self.dssp_asa_complex[num]),
                            'dssp_rel2_asa': (float(self.dssp_asa_mon[num])-float(self.dssp_asa_complex[num]))/self.dssp_asa_mon[num],
                            'b_factor_mean': self.bfactor[num],
                            'binana_dist_2_5': self.binana_b25[num],
                            'binana_dist_4_0': self.binana_b4[num],
                            'binana_h_bonds': self.binana_hydroge[num],
                            'binana_hydrophobic': self.binana_hydroph[num],
                            'binana_s_bridges':self.binana_salt[num],
                            'binana_pi_pi': self.binana_pis[num],
                            'binana_t_stacking': self.binana_ts[num],
                            'binana_cation_pi': self.binana_cation_pi[num],
                            'pssm_jsd': self.jsd[num-1]}

                    else:
                        row = {'pdb_id': self.pdbid,
                            'original_chains': self.chains,
                            'chain_position': 'A',
                            'residue_name': self.seq[num],
                            'residue_number': self.position[num],
                            'classifier': self.class_a[num],
                            'dssp_rel_asa_mon': self.dssp_asa_mon[num],
                            'dssp_rel_asa': self.dssp_asa_complex[num],
                            'dssp_delta_asa': float(self.dssp_asa_mon[num])-float(self.dssp_asa_complex[num]),
                            'dssp_rel2_asa': 0.0,
                            'b_factor_mean': self.bfactor[num],
                            'binana_dist_2_5': self.binana_b25[num],
                            'binana_dist_4_0': self.binana_b4[num],
                            'binana_h_bonds': self.binana_hydroge[num],
                            'binana_hydrophobic': self.binana_hydroph[num],
                            'binana_s_bridges':self.binana_salt[num],
                            'binana_pi_pi': self.binana_pis[num],
                            'binana_t_stacking': self.binana_ts[num],
                            'binana_cation_pi': self.binana_cation_pi[num],
                            'pssm_jsd': self.jsd[num-1]}

                else:
                    if self.dssp_asa_mon[num] != 0.0:
                        row = {'pdb_id': self.pdbid,
                            'original_chains': self.chains,
                            'chain_position': 'A',
                            'residue_name': self.seq[num],
                            'residue_number': self.position[num],
                            'classifier': self.class_a[num],
                            'dssp_rel_asa_mon': self.dssp_asa_mon[num],
                            'dssp_rel_asa': self.dssp_asa_complex[num],
                            'dssp_delta_asa': float(self.dssp_asa_mon[num])-float(self.dssp_asa_complex[num]),
                            'dssp_rel2_asa': (float(self.dssp_asa_mon[num])-float(self.dssp_asa_complex[num]))/self.dssp_asa_mon[num],
                            'b_factor_mean': self.bfactor[num],
                            'binana_dist_2_5': float(0),
                            'binana_dist_4_0': float(0),
                            'binana_h_bonds': float(0),
                            'binana_hydrophobic': float(0),
                            'binana_s_bridges':float(0),
                            'binana_pi_pi': float(0),
                            'binana_t_stacking': float(0),
                            'binana_cation_pi': float(0),
                            'pssm_jsd': self.jsd[num-1]}

                    else:
                        row = {'pdb_id': self.pdbid,
                            'original_chains': self.chains,
                            'chain_position': 'A',
                            'residue_name': self.seq[num],
                            'residue_number': self.position[num],
                            'classifier': self.class_a[num],
                            'dssp_rel_asa_mon': self.dssp_asa_mon[num],
                            'dssp_rel_asa': self.dssp_asa_complex[num],
                            'dssp_delta_asa': float(self.dssp_asa_mon[num])-float(self.dssp_asa_complex[num]),
                            'dssp_rel2_asa': 0.0,
                            'b_factor_mean': self.bfactor[num],
                            'binana_dist_2_5': float(0),
                            'binana_dist_4_0': float(0),
                            'binana_h_bonds': float(0),
                            'binana_hydrophobic': float(0),
                            'binana_s_bridges': float(0),
                            'binana_pi_pi': float(0),
                            'binana_t_stacking': float(0),
                            'binana_cation_pi': float(0),
                            'pssm_jsd': self.jsd[num-1]}

            else:
                if self.class_a[num] == '3':
                    if self.dssp_asa_mon[num] != 0.0:
                        row = {'pdb_id': self.pdbid,
                            'original_chains': self.chains,
                            'chain_position': 'B',
                            'residue_name': self.seq[num],
                            'residue_number': self.position[num],
                            'classifier': self.class_a[num],
                            'dssp_rel_asa_mon': self.dssp_asa_mon[num],
                            'dssp_rel_asa': self.dssp_asa_complex[num],
                            'dssp_delta_asa': float(self.dssp_asa_mon[num])-float(self.dssp_asa_complex[num]),
                            'dssp_rel2_asa': (float(self.dssp_asa_mon[num])-float(self.dssp_asa_complex[num]))/self.dssp_asa_mon[num],
                            'b_factor_mean': self.bfactor[num],
                            'binana_dist_2_5': self.binana_b25[num],
                            'binana_dist_4_0': self.binana_b4[num],
                            'binana_h_bonds': self.binana_hydroge[num],
                            'binana_hydrophobic': self.binana_hydroph[num],
                            'binana_s_bridges':self.binana_salt[num],
                            'binana_pi_pi': self.binana_pis[num],
                            'binana_t_stacking': self.binana_ts[num],
                            'binana_cation_pi': self.binana_cation_pi[num],
                            'pssm_jsd': self.jsd[num-1]}

                    else:
                        row = {'pdb_id': self.pdbid,
                            'original_chains': self.chains,
                            'chain_position': 'B',
                            'residue_name': self.seq[num],
                            'residue_number': self.position[num],
                            'classifier': self.class_a[num],
                            'dssp_rel_asa_mon': self.dssp_asa_mon[num],
                            'dssp_rel_asa': self.dssp_asa_complex[num],
                            'dssp_delta_asa': float(self.dssp_asa_mon[num])-float(self.dssp_asa_complex[num]),
                            'dssp_rel2_asa': 0.0,
                            'b_factor_mean': self.bfactor[num],
                            'binana_dist_2_5': self.binana_b25[num],
                            'binana_dist_4_0': self.binana_b4[num],
                            'binana_h_bonds': self.binana_hydroge[num],
                            'binana_hydrophobic': self.binana_hydroph[num],
                            'binana_s_bridges':self.binana_salt[num],
                            'binana_pi_pi': self.binana_pis[num],
                            'binana_t_stacking': self.binana_ts[num],
                            'binana_cation_pi': self.binana_cation_pi[num],
                            'pssm_jsd': self.jsd[num-1]}

                else:
                    if self.dssp_asa_mon[num] != 0.0:
                        row = {'pdb_id': self.pdbid,
                            'original_chains': self.chains,
                            'chain_position': 'B',
                            'residue_name': self.seq[num],
                            'residue_number': self.position[num],
                            'classifier': self.class_a[num],
                            'dssp_rel_asa_mon': self.dssp_asa_mon[num],
                            'dssp_rel_asa': self.dssp_asa_complex[num],
                            'dssp_delta_asa': float(self.dssp_asa_mon[num])-float(self.dssp_asa_complex[num]),
                            'dssp_rel2_asa': (float(self.dssp_asa_mon[num])-float(self.dssp_asa_complex[num]))/self.dssp_asa_mon[num],
                            'b_factor_mean': self.bfactor[num],
                            'binana_dist_2_5': float(0),
                            'binana_dist_4_0': float(0),
                            'binana_h_bonds': float(0),
                            'binana_hydrophobic': float(0),
                            'binana_s_bridges':float(0),
                            'binana_pi_pi': float(0),
                            'binana_t_stacking': float(0),
                            'binana_cation_pi': float(0),
                            'pssm_jsd': self.jsd[num-1]}

                    else:
                        row = {'pdb_id': self.pdbid,
                            'original_chains': self.chains,
                            'chain_position': 'B',
                            'residue_name': self.seq[num],
                            'residue_number': self.position[num],
                            'classifier': self.class_a[num],
                            'dssp_rel_asa_mon': self.dssp_asa_mon[num],
                            'dssp_rel_asa': self.dssp_asa_complex[num],
                            'dssp_delta_asa': float(self.dssp_asa_mon[num])-float(self.dssp_asa_complex[num]),
                            'dssp_rel2_asa': 0.0,
                            'b_factor_mean': self.bfactor[num],
                            'binana_dist_2_5': float(0),
                            'binana_dist_4_0': float(0),
                            'binana_h_bonds': float(0),
                            'binana_hydrophobic': float(0),
                            'binana_s_bridges': float(0),
                            'binana_pi_pi': float(0),
                            'binana_t_stacking': float(0),
                            'binana_cation_pi': float(0),
                            'pssm_jsd': self.jsd[num-1]}



            self.dataframe.loc[num] = pd.Series(row)

    def correction(self, do = 'make'):
    '''
    Applies a correction factor to all features, excluding accessible
    surface area. This factor was precalculated before, and loaded as 
    norms. This step can be omited from the pipeline by changin the 'do'
    variable form "make" to "stay".
    This method also creates an output csv file containing all the
    information of interest.
    '''
        self.buildDataFrame()
        if do == 'make':
            dic = {}

            home = os.getcwd()
            m_docs = home + '/data/normalizators.csv'
            norms = pd.read_csv(m_docs).set_index('Residues')

            n_core = norms['Non-surface']
            n_surface = norms['Non-interfacial Surface']
            n_interface = norms['Interfacial Surface']

            allfeatures = self.dataframe

            corr_features = ['b_factor_mean',
                'binana_dist_2_5',
                'binana_dist_4_0',
                'binana_h_bonds',
                'binana_hydrophobic',
                'binana_s_bridges',
                'binana_pi_pi',
                'binana_t_stacking',
                'binana_cation_pi',
                'pssm_jsd']

            for index, row in allfeatures.iterrows():
                if row['classifier'] == '0':
                    n_val = n_core.loc[n_core.index == row['residue_name']]
                    w = row.copy().drop(corr_features)
                    z = row[corr_features].multiply(n_val.values[0])
                    y = pd.concat([w,z])
                    dic[index] = y
                elif row['classifier'] == '2':
                    n_val = n_surface.loc[n_surface.index == row['residue_name']]
                    w = row.copy().drop(corr_features)
                    z = row[corr_features].multiply(n_val.values[0])
                    y = pd.concat([w,z])
                    dic[index] = y
                elif row['classifier'] == '3':
                    n_val = n_interface.loc[n_interface.index == row['residue_name']]
                    w = row.copy().drop(corr_features)
                    z = row[corr_features].multiply(n_val.values[0])
                    y = pd.concat([w,z])
                    dic[index] = y
                else:
                    continue

            f_df = pd.DataFrame.from_dict(dic).transpose()
            out_name = self.pdbid + "_mensadb.csv"
            f_df.to_csv(out_name)
            self.normalized = f_df

        elif do == 'make':
            self.normalized = self.dataframe

prot = Run(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
prot.correction()
