#!/bin/bash
for idx in 0 
#{0..15}
#for idx in {16..31}
#for idx in {31..47}
	do
		echo "python 1_model_qso_inp.py ${idx}"
		python 1_model_qso_inp.py ${idx}
	done