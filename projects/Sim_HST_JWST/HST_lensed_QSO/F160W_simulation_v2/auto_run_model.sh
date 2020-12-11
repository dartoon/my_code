#!/bin/bash
#for idx in {0..7}
#for idx in {8..15}
#for idx in {16..23}
#for idx in {24..31}
#for idx in {32..39}
for idx in {40..37}
	do
		echo "python 1_model_qso_inp.py ${idx}"
		python 1_model_qso_inp.py ${idx}
	done