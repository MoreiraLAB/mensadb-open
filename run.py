from process_binana import retrieve_features as binana_features
from dssp_features import DSSP_features
from generate_class import run as classifier
import pandas as pd
from statistics import mean

def convDic(A, last_aa):
	B = {}
	for entry in A.keys():
		B[entry + last_aa] = A[entry]
	return B

class Run:

	def __init__(self, pdbid, chains):
		self.pdbid   = pdbid
		self.chains  = chains
		self.file    = self.pdbid + '_' + self.chains + '.pdb'

	def generateClass(self):
		self.classes = classifier(self.file)
		self.complex = pdbid + '_' + chains + '.pdb'
		self.binana  = 'binana_' + self.pdbid + '_' + chains + '.pdb'
		self.dssp_A  = pdbid + '_' + self.chains + '_A.pdb'
		self.dssp_B  = pdbid + '_' + self.chains + '_B.pdb'

	def getSequence(self):
		self.generateClass()
		file = open(self.complex, 'r').readlines()
		self.seq      = {}
		self.mon      = {}
		self.orchains = {}
		self.position = {}
		self.bfactor  = {}

		n = 1

		atoms   = [file[0][22:26].replace(' ', '')]
		BF      = []

		for line in file:

			bf      = line[61:67].replace(' ', '')
			chain   = line[21:23].replace(' ', '')
			typ     = line[0:6].replace(' ', '')
			atom    = line[11:16].replace(' ', '')
			aa      = line[17:20].replace(' ', '')
			pos     = line[22:26].replace(' ', '')

			if pos != atoms[-1] and typ == 'ATOM':
				atoms = []
				atoms.append(pos)
				self.bfactor[int(n)] = mean(BF)
				self.seq[int(n)] = aa
				self.mon[int(n)] = chain
				self.position[int(n)] = pos
				if chain == 'A':
					self.orchains[int(n)] = self.chains[0]
				elif chain == 'B':
					self.orchains[int(n)] = self.chains[1]
				n += 1
				print(aa, BF)
				BF = []
				BF.append(float(bf))
			elif typ == 'ATOM':
				atoms.append(pos)
				BF.append(float(bf))
			else:
				continue

		self.last_aa = 0

		for entry in self.mon.keys():
			if self.mon[entry] == 'A':
				self.last_aa = entry
			else:
				break

	def getDSSPFeatures(self):
		n_features = {
			0: 'DSSP index',
			1: 'Amino acid number',
			2: 'Amino acid code',
			3: 'Chain',
			4: 'Secondary Structure',
			5: 'BP',
			6: 'ASA',
			7: 'NH-->O_1_relidx',
			8: 'O-->NH_1_relidx',
			9: 'NH-->O_1_energy',
			10: 'O-->NH_1_energy',
			11: 'TCO',
			12: 'KAPPA',
			13: 'Alpha',
			14: 'Phi',
			15: 'Psi',
			16: 'X-CA',
			17: 'Y-CA',
			18: 'Z-CA'
			}

		self.dssp_asa_A         = DSSP_features(self.dssp_A, 6, pdb_chain = 'A')
		dssp_asa_pB             = DSSP_features(self.dssp_B, 6, pdb_chain = 'B')
		self.dssp_asa_B         = {}
		for entryB in dssp_asa_pB.keys():
			self.dssp_asa_B[entryB + self.last_aa] = dssp_asa_pB[entryB]
		self.dssp_asa_mon       = self.dssp_asa_A
		self.dssp_asa_mon.update(self.dssp_asa_B)
		self.dssp_asa_complex_A = DSSP_features(self.complex, 6, pdb_chain = 'A')
		dssp_asa_complex_pB = DSSP_features(self.complex, 6, pdb_chain = 'B')
		self.dssp_asa_complex_B = {}
		for entrypB in dssp_asa_complex_pB.keys():
			self.dssp_asa_complex_B[entrypB + self.last_aa] = dssp_asa_complex_pB[entrypB]
		self.dssp_asa_complex   = self.dssp_asa_complex_A
		self.dssp_asa_complex.update(self.dssp_asa_complex_B)

	def getBinanaFeatures(self):
		self.binana_feat    = binana_features(self.binana)

		self.binana_b25     = self.binana_feat['below_2.5A'][0]
		binana_b251         = convDic(self.binana_feat['below_2.5A'][1], self.last_aa)

		self.binana_b4      = self.binana_feat['below_4A'][0]
		binana_b41          = convDic(self.binana_feat['below_4A'][1], self.last_aa)

		self.binana_hydroge = self.binana_feat['HB'][0]
		binana_hydroge1     = convDic(self.binana_feat['HB'][1], self.last_aa)

		self.binana_hydroph = self.binana_feat['hydrophobic_contacts'][0]
		binana_hydroph1     = convDic(self.binana_feat['hydrophobic_contacts'][1], self.last_aa)

		self.binana_pis     = self.binana_feat['pipi_stack'][0]
		binana_pis1         = convDic(self.binana_feat['pipi_stack'][1], self.last_aa)

		self.binana_ts      = self.binana_feat['T_stack'][0]
		binana_ts1          = convDic(self.binana_feat['T_stack'][1], self.last_aa)

		self.binana_salt    = self.binana_feat['SB'][0]
		binana_salt1        = convDic(self.binana_feat['SB'][1], self.last_aa)

		self.binana_b25.update(binana_b251)
		self.binana_b4.update(binana_b41)
		self.binana_hydroge.update(binana_hydroge1)
		self.binana_hydroph.update(binana_hydroph1)
		self.binana_pis.update(binana_pis1)
		self.binana_ts.update(binana_ts1)
		self.binana_salt.update(binana_salt1)

	def buildDataFrame(self):
		self.getSequence()
		self.getDSSPFeatures()
		self.getBinanaFeatures()

		available_data = [
			self.dssp_asa_mon,
			self.dssp_asa_complex,
			self.binana_b25,
			self.binana_b4,
			self.binana_hydroge,
			self.binana_hydroph,
			self.binana_pis,
			self.binana_ts,
			self.binana_salt
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
			'dssp_asa_mon',
			'dssp_asa_cpx',
			'dssp_asa_delta',
			'dssp_asa_asa',
			'b_factor_mean',
			'b_factor_env',
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
			fb_env_mean = []

			for wind_ele in range(num - 5, num + 6):
				try:
					fb_env_mean.append(self.bfactor[wind_ele])
				except:
					continue

			if self.dssp_asa_mon[num] != 0.0:
				row = {'pdb_id': self.pdbid,
					'original_chains': self.chains,
					'chain_position': self.mon[num],
					'residue_name': self.seq[num],
					'residue_number': self.position[num],
					'classifier': 2,
					'dssp_asa_mon': self.dssp_asa_mon[num],
					'dssp_asa_cpx': self.dssp_asa_complex[num],
					'dssp_asa_delta': float(self.dssp_asa_mon[num])-float(self.dssp_asa_complex[num]),
					'dssp_asa_asa': (float(self.dssp_asa_mon[num])-float(self.dssp_asa_complex[num]))/self.dssp_asa_mon[num],
					'b_factor_mean': self.bfactor[num],
					'b_factor_env': mean(fb_env_mean),
					'binana_dist_2_5': self.binana_b25[num],
					'binana_dist_4_0': self.binana_b4[num],
					'binana_h_bonds': self.binana_hydroge[num],
					'binana_hydrophobic': self.binana_hydroph[num],
					'binana_s_bridges':self.binana_salt[num],
					'binana_pi_pi': self.binana_pis[num],
					'binana_t_stacking': self.binana_ts[num],
					'binana_cation_pi': 0.0,
					'pssm_jsd': 0.0}

			else:
				row = {'pdb_id': self.pdbid,
					'original_chains': self.chains,
					'chain_position': self.mon[num],
					'residue_name': self.seq[num],
					'residue_number': self.position[num],
					'classifier': 3,
					'dssp_asa_mon': self.dssp_asa_mon[num],
					'dssp_asa_cpx': self.dssp_asa_complex[num],
					'dssp_asa_delta': float(self.dssp_asa_mon[num])-float(self.dssp_asa_complex[num]),
					'dssp_asa_asa': 0.0,
					'b_factor_mean': self.bfactor[num],
					'b_factor_env': mean(fb_env_mean),
					'binana_dist_2_5': self.binana_b25[num],
					'binana_dist_4_0': self.binana_b4[num],
					'binana_h_bonds': self.binana_hydroge[num],
					'binana_hydrophobic': self.binana_hydroph[num],
					'binana_s_bridges':self.binana_salt[num],
					'binana_pi_pi': self.binana_pis[num],
					'binana_t_stacking': self.binana_ts[num],
					'binana_cation_pi': 0.0,
					'pssm_jsd': 0.0}

			self.dataframe.loc[num] = pd.Series(row)

	def correction(self, what = 'make'):
		self.buildDataFrame()
		self.dataframe.to_csv('dataframe.csv')

		if what == 'make':
			dic = {}

			norms = pd.read_csv('normalizators.csv').set_index('Residues')
			sander = pd.read_csv('sander_cons.csv', sep = ';').set_index('residue')

			n_core = norms['Non-surface']
			n_surface = norms['Non-interfacial Surface']
			n_interface = norms['Interfacial Surface']

			allfeatures = self.dataframe

			sander_features = ['dssp_asa_mon',
				'dssp_asa_cpx',
				'dssp_asa_delta',
				'dssp_asa_asa']

			corr_features = ['b_factor_mean',
				'b_factor_env',
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
				if row['classifier'] == 0:
					n_val = n_core.loc[n_core.index == row['residue_name']]
					s_val = sander.loc[sander.index == row['residue_name']]
					w = row.copy().drop(corr_features + sander_features)
					z = row[corr_features].multiply(n_val.values[0])
					x = row[sander_features].multiply(s_val.values[0][0])
					y = pd.concat([w,x,z])
					dic[index] = y
				elif row['classifier'] == 2:
					n_val = n_surface.loc[n_surface.index == row['residue_name']]
					s_val = sander.loc[sander.index == row['residue_name']]
					w = row.copy().drop(corr_features + sander_features)
					z = row[corr_features].multiply(n_val.values[0])
					x = row[sander_features].multiply(s_val.values[0][0])
					y = pd.concat([w,x,z])
					dic[index] = y
				elif row['classifier'] == 3:
					n_val = n_interface.loc[n_interface.index == row['residue_name']]
					s_val = sander.loc[sander.index == row['residue_name']]
					w = row.copy().drop(corr_features + sander_features)
					z = row[corr_features].multiply(n_val.values[0])
					x = row[sander_features].multiply(s_val.values[0][0])
					y = pd.concat([w,x,z])
					dic[index] = y
				else:
					continue

			f_df = pd.DataFrame.from_dict(dic).transpose()
			f_df.to_csv('test.csv')
			self.normalized = f_df

		elif what == 'stay':
			self.normalized = self.dataframe

prot = Run('1a0t', 'PQ')
prot.correction()
