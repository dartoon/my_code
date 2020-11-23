import sys
from gsf import gsf

# import argparse
# parser = argparse.ArgumentParser(description='Run gsf.')
# parser.add_argument('parfile', metavar='parfile', type=str, help='Configuration file.')
# parser.add_argument('fplt', metavar='fplt', type=int, help='Flag for run (int: 0,1,2,3).')
# parser.add_argument('--id', default=None, help='Manual input for object ID.')
# args = parser.parse_args()

gsf.run_gsf_all('sample.input', 1, idman=None)

#Tested Using Z/Z_sum = 0.1, 0.2, 0.4, T/Gyr = 0.1, 0.2, 0.4
#Find: inferred M/L are very similar (<0.1 dex)
