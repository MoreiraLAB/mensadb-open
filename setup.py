import os
import urllib.request

biopython = 'easy_install -f http://biopython.org/DIST/ biopython'
numpy = 'python -m pip install --user numpy scipy matplotlib ipython jupyter pandas sympy nose'
dssp = 'conda install -c salilab dssp'
pandas = 'python3 -m pip install --upgrade pandas'
dit = 'pip install dit'
vmd = 'conda install -c conda-forge vmd-python'

url = 'https://downloads.sourceforge.net/project/binana/binana_1_2_0.zip?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbinana%2Ffiles%2Flatest%2Fdownload&ts=1553616248'
urllib.request.urlretrieve(url, 'binana.py')

installs = [numpy, biopython, pandas, dit, vmd]

for i in installs:
	os.system(i)
