import sys
import os
module_path = os.path.abspath(os.path.join('../src'))
if module_path not in sys.path:
    sys.path.append(module_path)
from demographicmodel import DemographicModel

'''
Read a fastNeutrino model file and output the flags to run simulate_joint_sfs.py
'''

script, fn = sys.argv
model = DemographicModel(fn)
model.rescale()
model.print_msprime_flags()
