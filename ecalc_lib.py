#!/usr/bin/python3
import sys
import os
import numpy as np
import csv
import pandas as pd
import tempfile
import subprocess
import signal


def get_ecalc(file):
	os.environ['ECALCDATA'] = '/serv/www/html_lilian/libs/ecalc-master/data'
	os.environ['BIOPLIB_DEPRECATED_QUIET']='silence'
	os.environ['DATADIR']='/serv/www/html_lilian/libs/data'
	signal.signal(signal.SIGSEGV, signal.SIG_IGN)
	energy_list=[]
	pdh_out = open("./out.pdh", "a")
	pdh_out.close()
	command='pdbgetchain L,H {}| /serv/www/html_lilian/libs/mutmodel -R | /serv/www/html_lilian/libs/bioptools-master/src/pdbchain | /serv/www/html_lilian/libs/bioptools-master/src/pdbhadd -c | /serv/www/html_lilian/libs/bioptools-master/src/pdbcter -c | /serv/www/html_lilian/libs/bioptools-master/src/pdbrenum > {}'.format(file, "./out.pdh")
	pdh_out = subprocess.check_output(command, shell=True)

	command='/serv/www/html_lilian/libs/ecalc-master/src/ecalc -p {}'.format("./out.pdh")
	ecalc_out=subprocess.check_output(command, shell=True)
	ecalc_out=str(ecalc_out, 'utf-8')	
					
	for line in ecalc_out.split('\n'):
		if 'Base structure energy' in line:
			energy=line.split(' = ')[1]
			energy=energy.replace('/n','')
			energy=float(energy)
			energy_list+=[energy]
	
	energy_df = pd.DataFrame(data=energy_list, columns=['energy'], index=None)

	os.remove("./out.pdh")
	return energy_df

if __name__ == '__main__':
	energy_df=get_ecalc(sys.argv[1])
	print(energy_df)