#!/usr/bin/python3

#########################################################################################
#Purpose of Programme:
#1: take each modelled loop structure and run loop flapper programme
#2: calculate different features of each flapped structure: i.e. ecalc (only VdW, electrostats,..), accessibility, packing quality, etc. 
#   build machine learning model 
#3: output the optimal angle of the structure
#4: assess quality of loop flapper programm by calculating the angle of the original PDB file and comparing with prediction of optimal angle

#export PYTHONPATH=/home/lilian/lib/python3.6/site-packages/
import sys
import os
import signal
import numpy as np
import csv
import pandas as pd
import itertools
import math
from biopandas.pdb import PandasPdb
from scipy.spatial.transform import Rotation
import subprocess
import accessibility_lib
import protrusion_lib
from functools import reduce
import shutil
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LinearRegression
import seaborn as sns
from sklearn.model_selection import KFold
from sklearn.metrics import matthews_corrcoef
import save_RMS_lib
import ast
from zipfile import ZipFile
import glob
import subprocess
#PART 1: Reading the input PDB model

def input_process(pdb_file):
	ppdb = PandasPdb().read_pdb(os.path.join(pdb_file))
	ppdb=ppdb.df['ATOM']

	#get index of start-atom and stop-atom of loop
	index_95 = ppdb.loc[(ppdb['residue_number'] == 95)& (ppdb['chain_id']=="H")].index.tolist()[0]
	index_102 = ppdb.loc[(ppdb['residue_number'] == 102)&(ppdb['chain_id']=="H")].index.tolist()[-1]


	#get coordinates of full structure
	full_structure=ppdb[["x_coord", "y_coord", "z_coord"]]
	full_structure=full_structure.to_numpy()

	#get base 95 point
	ppdb95=ppdb.loc[ppdb['residue_number'] == 95]
	ppdb95=ppdb95.loc[ppdb['chain_id'] == "H"]
	ppdb95=ppdb95.loc[ppdb['atom_name'] == "CA"]
	ppdb95=ppdb95[["x_coord","y_coord", "z_coord"]]
	ppdb95=ppdb95.to_numpy()

	#get base 102
	ppdb102=ppdb.loc[ppdb['residue_number'] == 102]
	ppdb102=ppdb102.loc[ppdb['chain_id'] == "H"]
	ppdb102=ppdb102.loc[ppdb['atom_name'] == "CA"]
	ppdb102=ppdb102[["x_coord","y_coord", "z_coord"]]
	ppdb102=ppdb102.to_numpy()
	return (ppdb, index_95,index_102, full_structure, ppdb95, ppdb102)



#shift the whole structure so that the base point H95 is at the origin; then rotate so second basepoint H102 as on x-axis.
def move_to_origin(full_structure, ppdb95, ppdb102):
	#move entire structure so that resH95 is at origin
	for i in range(0, np.shape(full_structure)[0]):
		#for x
		full_structure[i, 0]=full_structure[i, 0]-ppdb95[0,0]
		#for y
		full_structure[i, 1]=full_structure[i, 1]-ppdb95[0,1]
		#for z
		full_structure[i, 2]=full_structure[i, 2]-ppdb95[0,2]
	#define new H102 base-point position after transformation
	ppdb102[0,0]=ppdb102[0,0]-ppdb95[0,0]
	ppdb102[0,1]=ppdb102[0,1]-ppdb95[0,1]
	ppdb102[0,2]=ppdb102[0,2]-ppdb95[0,2]

	#align vectors from origin to res H102 and vector from origin along x axis
	rotation_origin, rmsd=Rotation.align_vectors(ppdb102,[[1,0,0]])
	rotated_full_structure=rotation_origin.apply(full_structure)
	rotated_full_df = pd.DataFrame(data=rotated_full_structure, columns=["x_coord", "y_coord", "z_coord"])

	return(ppdb95, ppdb102, rotated_full_df)


#PART2 flap the loop and make PDB files with altered angles

#rotate the loop structure by one degree around the x-axis
#outputs the rotated_full_df that has a changed loop structure, flapped by 1 degree
def flap_loop(rotated_full_df, index_95,index_102, angle_diff):
	r = Rotation.from_euler('x', angle_diff, degrees=True)
	for i in range(index_95, index_102+1):
		a_loop_vector=rotated_full_df.iloc[[i]]
		new_loop_vector=r.apply(a_loop_vector)[0]
		rotated_full_df.loc[i] = new_loop_vector

	return (rotated_full_df)

#Creates a PDB file in the save_path directory; name of file will be flapped_structure and the difference in angle to the original appended. 
#takes the structure at it's original position with the rotated loop as input
def pdb_integrate(original_rotation_structure, pdb_file, angle_diff, save_path):
	new_structure=pd.DataFrame(data=original_rotation_structure, index=None, columns=["x_coord","y_coord","z_coord"])
	ppdb = PandasPdb().read_pdb(os.path.join(pdb_file))
	ppdb.df['ATOM'][['x_coord','y_coord','z_coord']] = new_structure
	path=os.path.join(save_path,"flapped_structure_"+str(angle_diff)+".pdb")
	ppdb.to_pdb(path=path, 
			records=None, 
			gz=False, 
			append_newline=True)
	return (None)

#full loop flapping process, flap the loop from 
def full_flap (rotated_full_df, index_95, index_102, rotation_structure, ppdb95, pdb_file, save_path, max_angle_diff):
	for direction in ["clockwise", "counterclockwise"]:
		if direction=="clockwise":
			angles=range(0, int(max_angle_diff)+1,1)
		else:
			angles=range(0,(-1)*int(max_angle_diff)-1,-1)
		for angle_diff in angles:
			flapped_full_df=flap_loop(rotated_full_df, index_95,index_102, angle_diff)
			#return to initial position (i.e reverse of function move_to_origin)
			#first reverse rotation to place H102 on x-axis
			reverse_rotation=rotation_structure.inv()
			original_rotation_structure=reverse_rotation.apply(flapped_full_df)
			#then reverse transformation to place H95 at origin
			for i in range(0, np.shape(original_rotation_structure)[0]):
				#for x
				original_rotation_structure[i, 0]=original_rotation_structure[i, 0]+ppdb95[0,0]
				#for y
				original_rotation_structure[i, 1]=original_rotation_structure[i, 1]+ppdb95[0,1]
				#for z
				original_rotation_structure[i, 2]=original_rotation_structure[i, 2]+ppdb95[0,2]
			pdb_integrate(original_rotation_structure, pdb_file, angle_diff, save_path)
	return None

def run(pdb_file, save_path, max_angle_diff):
	#process the input pdb model file, get the index of the base residues
	(ppdb, index_95,index_102, full_structure, ppdb95, ppdb102)=input_process(pdb_file)

	#move the full structure so that H95 is on the origin and H102 is on the x-axis 
	(ppdb95, ppdb102, rotated_full_df)=move_to_origin(full_structure, ppdb95, ppdb102)
	
	#align vectors from origin to res H102 and vector from origin along x axis   #why again???
	rotation_structure, rmsd=Rotation.align_vectors(ppdb102,[[1,0,0]])
	rotated_full_structure=rotation_structure.apply(full_structure)
	rotated_full_df = pd.DataFrame(data=rotated_full_structure, columns=["x_coord", "y_coord", "z_coord"])  

	#flap the loop, output pdb files for all angles from -max_angle_diff to +max_angle_diff; 
	#file of 0 degree angle difference should correspond to original file, use as sanity check
	full_flap (rotated_full_df, index_95, index_102, rotation_structure, ppdb95, pdb_file, save_path, max_angle_diff)
	return None

#PART3 FEATURE CALCULATION: Calculate features of flapped PDB structures

def get_ecalc(save_path):
	here="/home/lilian/loop_flapper"#os.path.dirname(__file__)
	qualiloop_path="/home/lilian/sync_project/WWW/qualiloop"
	run_gromacs_path=os.path.join(here, "run_gromacs.pm")
	os.environ['ECALCDATA'] = os.path.join(here,"ecalc_data")
	os.environ['BIOPLIB_DEPRECATED_QUIET']='silence'
	#os.environ['DATADIR']=os.path.join(qualiloop_path,"data")
	mutmodel_path=os.path.join(here,"mutmodel")
	#print(mutmodel_path)
	pdbchain_path=os.path.join(qualiloop_path,"bin/bioptools/bioptools-1.9/src/pdbchain")
	pdbhadd_path=os.path.join(qualiloop_path,"bin/bioptools/bioptools-1.9/src/pdbhadd")
	pdbcter_path=os.path.join(qualiloop_path,"bin/bioptools/bioptools-1.9/src/pdbcter")
	pdbrenum_path=os.path.join(qualiloop_path,"bin/bioptools/bioptools-1.9/src/pdbrenum")
	#ecalc_path=os.path.join(qualiloop_path,"bin/ecalc-master/src/ecalc")
	#input(ecalc_path)
	#input(os.path.abspath(qualiloop_path))

	#signal.signal(signal.SIGSEGV, signal.SIG_IGN)
	energy_list=[]
	for file in os.listdir(save_path):
		if file.endswith(".pdb"):
			filename=os.path.splitext(os.path.basename(file))[0]
			file=os.path.join(save_path,file)
			#command='{}|{} -R | {} | {} -c | {} > {}'.format(file,mutmodel_path,pdbchain_path,pdbhadd_path,pdbrenum_path, os.path.join(save_path, filename+'_pdh'+".pdb"))
			#p = subprocess.run(command, shell=True)
			#command='{} -p {}'.format("ecalc", os.path.join(save_path, filename+"_pdh"+'.pdb'))
			#ecalc_out=subprocess.check_output(command, shell=True)
			#ecalc_out=str(ecalc_out, 'utf-8')	
			command='perl {} {} /home/lilian/loop_flapper/tmp'.format(run_gromacs_path, file)#os.path.join(save_path, filename+"_pdh"+'.pdb'))
			energy_out=subprocess.check_output(command, shell=True)
			#energy_out = open(energy_out_path, "r")
			print(energy_out)
			energy_out=str(energy_out, 'utf-8')	
			for line in energy_out.split('\n'):
				if 'Potential Energy' in line:
					energy=line.split(' = ')[1]
					energy=energy.replace('/n','')
					try:
						energy=float(energy)
					except:
						energy=np.nan
					energy_list+=[[str(file), energy]]
			#os.remove(os.path.join(save_path, filename+"_pdh"+'.pdb'))
		print(energy)
	energy_df = pd.DataFrame(data=energy_list, columns=['file', 'energy'], index=None)
	return energy_df

def accessibility (save_path):
	accessibility_df=pd.DataFrame()
	for file in os.listdir(save_path):
		if file.endswith(".pdb"):
			filename=os.path.splitext(os.path.basename(file))[0]
			file=os.path.join(save_path,file)
			if accessibility_df.empty != True:
				new_line=accessibility_lib.pdbsolv_run(file)
				new_line["file"]=str(file)
				accessibility_df=pd.concat([new_line, accessibility_df])
			else:
				accessibility_df=accessibility_lib.pdbsolv_run(file)
				accessibility_df["file"]=str(file)
	return (accessibility_df)

def get_protrusion(pdb_file):
	print(pdb_file)
	protrusion_df= protrusion_lib.calc_protrusion(pdb_file)
	print(protrusion_df)
	protrusion_res= protrusion_df["tip_pos"].loc[0]
	try:
		protrusion_res=int(protrusion_res)
	except:
		pass

	ppdb = PandasPdb().read_pdb(os.path.join(pdb_file))
	ppdb=ppdb.df['ATOM']
	protrusion_res = ppdb.loc[(ppdb['residue_number'] == protrusion_res)& (ppdb['chain_id']=="H")& (ppdb['atom_name'] == "CA")]
	A=protrusion_res[["x_coord","y_coord","z_coord"]]
	try:
		A=A.to_numpy()[0]
	except:
		A=np.nan
	return A,protrusion_df

def full_protrusion(save_path):
	protrusion_df=pd.DataFrame()
	for file in os.listdir(save_path):
		file_path=os.path.join(save_path, file)
		A,protrusion_row=get_protrusion(file_path)
		protrusion_row["file"]=file_path
		protrusion_df=pd.concat([protrusion_row, protrusion_df])
	return protrusion_df






def get_RMSD_actual(save_path, filename, actual_pdbs):
	actual=os.path.join(actual_pdbs,filename+".pdb")
	RMSD_num_df=pd.DataFrame()
	#command='pdbchain -c L,H {}| pdbgetchain L,H | {} -c | {} -c | {} > {}'.format(actual, "pdbhadd", "pdbcter", "pdbrenum", os.path.join(save_path, filename+'actual'+".pdb"))
	command='pdbchain -c L,H {}| pdbgetchain L,H > {}'.format(actual, os.path.join(save_path, os.path.join(filename+'actual'+".pdb")))
	print(command)
	try:
		p = subprocess.check_output(command, shell=True)
	except:
		return pd.DataFrame()
	p_out=str(p, 'utf-8') 
	print(p_out)
	for file in os.listdir(save_path):
		RMSD_row=pd.DataFrame()
		RMSD_row["file"]=[str(file)]
		file=os.path.join(save_path, file)
		actual=os.path.join(save_path, os.path.join(filename+'actual'+".pdb"))
		command="/home/lilian/sync_project/WWW/qualiloop/bin/ProFit_V3.3/src/profit -f {} {} {}".format("/home/lilian/loop_flapper/RMS_calcglobal_AA", actual, file)  #can use -h to include hetatoms
		profit_out=subprocess.check_output(command, shell=True)
		profit_out=str(profit_out, 'utf-8') 
		print(profit_out)
		counter=0
		
		for line in profit_out.split('\n'):
			counter+=1
			if "RMS:" in line:
				RMS=line.split("RMS: ")[1]
				print(RMS)
				RMS=float(RMS)
				RMSD_row["globalallatoms"]=[RMS]
			else:
				pass
		RMSD_num_df=pd.concat([RMSD_row, RMSD_num_df])
	print(RMSD_num_df)
	#input()
	#for file in os.listdir(save_path):
	#	os.remove(os.path.join(save_path, file))
	#get best angle according to lowest RMSD value
	try:
		all_RMSD_df_cp=RMSD_num_df[["globalallatoms", "file"]]
	except:
		all_RMSD_df_cp=pd.DataFrame()
		return pd.DataFrame()

	best_RMSD = min(all_RMSD_df_cp["globalallatoms"].to_list())  # minimum value
	best_angle=all_RMSD_df_cp.loc[(all_RMSD_df_cp['globalallatoms'] == best_RMSD, "file")]
	#input(best_angle)
	for i in range(best_angle.size):
		best_angle_ele=best_angle[i]
		if "actual" in str(best_angle_ele):
			print("here")
			new_list=all_RMSD_df_cp["globalallatoms"].to_list()
			new_list.remove(best_RMSD)
			#input(new_list)
			best_RMSD = min(new_list)  # minimum value
			best_angle=all_RMSD_df_cp.loc[(all_RMSD_df_cp['globalallatoms'] == best_RMSD, "file")].tolist()
			best_angle=[x.replace("flapped_structure_","").replace(".pdb","") for x in best_angle]
		#input(best_angle)
	try:
		if best_angle.empty ==True:
			best_angle=np.nan
	except:
		if len(best_angle)==0:
			best_angle=np.nan

	best_angle_df=pd.DataFrame()
	try:
		best_angle_df["best_angle"]=[best_angle[0]]
		best_angle_df["best_angle_others"]=[best_angle]
	except:
		best_angle_df["best_angle"]=[best_angle]
		best_angle_df["best_angle_others"]=[[np.nan,np.nan]]

	
	best_angle_df["file"]=filename
	
	#input(best_angle_df)

	return best_angle_df



def get_original_angle(protrusion_point):
	#A=maximal protruion point, point B is defined as point A projected onto the Y-Z Plane 
	#As the rotation is only about the X-axis, the X-coordinate can be disregarded
	B=[0,0,0]
	angle=math.atan2(protrusion_point[2], protrusion_point[1]) - math.atan2(B[2], B[1])
	return angle

def get_angle_df(save_path):
	angle_df=pd.DataFrame()
	for file in os.listdir(save_path):
		if file.endswith(".pdb"):
			filename=os.path.splitext(os.path.basename(file))[0]
			file=os.path.join(save_path,file)
			if angle_df.empty != True:
				A=get_protrusion(file)
				angle=get_original_angle(A)
				new_line=pd.DataFrame(data={"angle":[angle]})
				new_line["file"]=str(file)
				angle_df=pd.concat([new_line, angle_df])
			else:
				A=get_protrusion(file)
				angle=get_original_angle(A)
				angle_df=pd.DataFrame(data={"angle":[angle]})
				angle_df["file"]=str(file)
				
	return (angle_df)



def merge_features(save_path, filename, actual_pdbs, max_angle_diff):
	accessibility_df=accessibility(save_path)
	print("ACCESSIBILITY;", accessibility_df)
	energy_df=get_ecalc(save_path)
	print("ENERGY;",energy_df)
	#angle_df=get_angle_df(save_path)
	#print("ANGLE;", angle_df)
	protrusion_df=full_protrusion(save_path)
	print("PROTRUSION;",protrusion_df)

	if accessibility_df.empty!=False or energy_df.empty!=False or protrusion_df.empty!=False:
		print("WARNING")
		return df_final, df_file_angles
	#df_final = reduce(lambda left,right: pd.merge(left,right,on='file'), feature_df_list)
	df_features=pd.merge(accessibility_df, energy_df, on="file")
	df_features=pd.merge(df_features,protrusion_df, on="file")
	print(df_features)
	row_df=df_to_row(df_features, max_angle_diff, save_path)

	best_angle_df=get_RMSD_actual(save_path, filename, actual_pdbs)
	if best_angle_df.empty==True:
		row_df=pd.DataFrame()
	final_df=pd.concat([best_angle_df, row_df], axis=1)
	for i in final_df.columns:
		print(i)
	try:
		final_df.drop(labels="file",axis=1)
	except:
		pass
	return final_df

def df_to_row(df_final, max_angle_diff, save_path):
	column_names=[]
	full_row=pd.DataFrame()
	df_final["file"].to_list()[0]
	for i in df_final.columns:
		if i=="file":
			pass
		else:
			list1=[x.replace(save_path,"") for x in df_final["file"].to_list()]
			list1=[x.replace(".pdb","") for x in list1]
			for a in list1:
			#range(((-1)*int(max_angle_diff)), (int(max_angle_diff)+1)):
				column_names+=[i+a]
			row_df=df_final[[i]].T
			row_df.columns = column_names
			full_row.reset_index(drop=True, inplace=True)
			row_df.reset_index(drop=True, inplace=True)
			full_row=pd.concat([row_df,full_row], axis=1)
			column_names=[]
	#input(full_row)
	return full_row

#do ML prediction using this feature table, try to predict the angle associated with the features
def maḱe_train_test_df(model_pdbs, save_path, max_angle_diff, actual_pdbs):
	#try:
	#	shutil.rmtree(save_path)
	#except:
	#	os.mkdir(save_path)
	full_feature_df=pd.DataFrame()
	mode_name=os.path.basename(model_pdbs)
	for file in os.listdir(model_pdbs):
		if file.endswith(".pdb") or file.endswith(".pdb.model"):
			filename=os.path.splitext(os.path.basename(file))[0]
			if filename.endswith(".pdb"):
				filename=os.path.splitext(os.path.basename(filename))[0]
			file=os.path.join(model_pdbs,file)
			if os.stat(file).st_size == 0:
				continue
		else:
			continue

		print("FILENAME:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",filename)
		#create loop flapping structures for file
		if not os.path.exists(save_path):
			os.mkdir(save_path)
		run(file, save_path, max_angle_diff)
		#make individual df of this file 
		final_df=merge_features(save_path, filename, actual_pdbs, max_angle_diff)
		full_feature_df=pd.concat([full_feature_df, final_df])
		#if len(full_feature_df[["best_angle"]].to_numpy().tolist())>=3:
		#	break

		full_feature_df.to_csv("full_feature_df"+mode_name+".csv")
	print(full_feature_df)
	full_feature_df = full_feature_df.dropna(subset=["best_angle"])
	full_feature_df.to_csv("full_feature_df.csv")
	zipObj = ZipFile(os.path.join(save_path, "flapper_structures.zip"), "w")
	files = glob.glob(os.path.join(save_path))
	for f in files:
		zipObj.write(f)
	zipObj.close()
	

	return full_feature_df


def train_model(full_feature_df):
	full_feature_df.dropna(subset = ["best_angle"], inplace=True)
	msk = np.random.rand(len(full_feature_df)) < 0.8
	train_df = full_feature_df[msk]
	test_df = full_feature_df[~msk]

	del_cols=[]
	for i in train_df.columns:
		if "tip_res" in i:
			del_cols.append(i)
	del_cols.append("file")
	train_df.drop(labels=del_cols,axis=1, inplace=True)
	test_df.drop(labels=del_cols,axis=1, inplace=True)
	train_df.drop(labels=["best_angle_others"],axis=1, inplace=True)
	train_y=train_df[["best_angle"]].to_numpy()
	train_df.drop(labels=["best_angle"],axis=1, inplace=True)
	train_X=train_df.to_numpy()
	print(train_X)
	print(train_y)
	
	#model=KNeighborsClassifier()
	#model= LogisticRegression()
	#model = RandomForestClassifier(n_estimators=200)
	model=LinearRegression()
	model.fit(train_X, train_y)
	print("score",model.score(train_X, train_y))
	#for i in train_df.columns:
		#sns_plot=sns.pairplot(full_feature_df, x_vars=i, y_vars="best_angle", size=7, aspect=0.7, kind='reg')
		#sns_plot.savefig('output'+i+'.png')
	num_folds = 10
	acc_per_fold = []
	MCC_per_fold=[]
	fold_no = 1
	kfold = KFold(n_splits=num_folds, shuffle=True)
	fold_dic={}
	for train, test in kfold.split(train_X, train_y):
		print('------------------------------------------------------------------------')
		print(f'Training for fold {fold_no} ...')
		#model=KNeighborsClassifier()
		model=LinearRegression()
		#model=LogisticRegression()
		#model = RandomForestClassifier(n_estimators=200)
		model.fit(train_X[train], train_y[train])
		rf_prediction = model.predict(train_X[test])
		#MCC_of_fold=matthews_corrcoef(train_y[test],rf_prediction, sample_weight=None )
		scores = model.score(train_X[test], train_y[test])
		acc_per_fold.append(scores * 100)
		#MCC_per_fold.append(MCC_of_fold)
		#fold_dic.update({MCC_of_fold:fold_no})
		fold_no = fold_no + 1
	# == Provide average scores ==
	joblib.dump(model, os.path.join(os.path.dirname(__file__),"random_forest_model"+".pkl")) 
	print('Score per fold')
	for i in range(0, len(acc_per_fold)):
	  print(f'> Fold {i+1} -  Accuracy: {acc_per_fold[i]}%')
	print('Average scores for all folds:')
	print(f'> Accuracy: {np.mean(acc_per_fold)} (+- {np.std(acc_per_fold)})')
	return (rf_prediction)


def param_parser():
	conf_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,add_help=False)
	conf_parser.add_argument("-config",dest="config", help="Specify config file", metavar="FILE")
	args, remaining_argv = conf_parser.parse_known_args()
	print("args",args)

	settings = {}
	if args.config:
		config = configparser.SafeConfigParser()
		config.read([args.config])
		print(config.items("model_1"))
		settings.update(dict(config.items("model_1")))
	#input(settings)

	parser = argparse.ArgumentParser(epilog="For any comments or inquiries please contact zcbtlm0@ucl.ac.uk",#fromfile_prefix_chars='@', 
	  formatter_class=argparse.RawDescriptionHelpFormatter,parents=[conf_parser],#, argparse.ArgumentDefaultsHelpFormatter
	  description=('''\
		 Please do not mess up this text!
		 --------------------------------
			   I have indented it
			   exactly the way
			   I want it
		 ''')) #textwrap.dedent
	parser.set_defaults(**settings)
	verbosity = parser.add_mutually_exclusive_group()
	verbosity.add_argument("-v", "--verbose", action="store_true")
	verbosity.add_argument("-q", "--quiet", action="store_true")

	parser.add_argument('-save_path',default=None, type=Path,help='save the flapped structures and other generated files in a folder with this name.')
	parser.add_argument('-max_angle',type=int, help='', default=None)
	parser.add_argument('-pdb_dir',type=Path, required=True, help='path to file with pdb files to be flapped')


	args2 = vars(parser.parse_args())

	settings.update({k: v for k, v in args2.items() if v is not None})
	settings.update({k: v for k, v in args2.items() if v is None and k not in settings.keys()}) 
	for key,val in settings.items():
		try:
			val=eval(val)
			settings.update({key:val})
		except:
			pass
	train_df, test_df=maḱe_train_test_df(settings["pdb_dir"], settings["save_path"], settings["max_angle"], actual_pdbs)
	train_model(train_df, test_df)


if __name__ == '__main__':
	model_pdbs=sys.argv[1]
	actual_pdbs=sys.argv[2]
	save_path=sys.argv[3]
	max_angle_diff=sys.argv[4]

	try:
		readmode=sys.argv[5]
		readmode=read_csv(readmode)
		print(readmode)
	except:
		pass
	
	full_feature_df=maḱe_train_test_df(model_pdbs, save_path, max_angle_diff, actual_pdbs)
	rf_prediction=train_model(full_feature_df)
	print(rf_prediction)






"""FILENAME:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2WUC_1
ACCESSIBILITY;     Access   Relacc    Scacc  Screlacc  Access_avg  Relacc_avg  Scacc_avg  Screlacc_avg                                           file
0  304.985  152.715  272.533   164.202   38.123125   19.089375  34.066625     20.525250   ./flapped_full_test/flapped_structure_-4.pdb
0  200.181  100.423  179.056   107.426   25.022625   12.552875  22.382000     13.428250    ./flapped_full_test/flapped_structure_1.pdb
0  311.697  151.379  276.544   158.413   38.962125   18.922375  34.568000     19.801625   ./flapped_full_test/flapped_structure_-1.pdb
0  205.414  102.948  184.405   110.045   25.676750   12.868500  23.050625     13.755625    ./flapped_full_test/flapped_structure_2.pdb
0  198.990   99.561  177.582   106.084   24.873750   12.445125  22.197750     13.260500  ./flapped_full_test/flapped_structure_-10.pdb
0  307.757  152.378  274.596   162.238   38.469625   19.047250  34.324500     20.279750   ./flapped_full_test/flapped_structure_-3.pdb
0  260.552  132.455  234.320   145.275   32.569000   16.556875  29.290000     18.159375   ./flapped_full_test/flapped_structure_-7.pdb
0  224.186  112.187  202.630   119.764   28.023250   14.023375  25.328750     14.970500    ./flapped_full_test/flapped_structure_4.pdb
0  312.420  151.203  276.961   157.914   39.052500   18.900375  34.620125     19.739250   ./flapped_full_test/flapped_structure_10.pdb
0  291.503  149.291  259.776   161.041   36.437875   18.661375  32.472000     20.130125   ./flapped_full_test/flapped_structure_-5.pdb
0  272.092  141.563  242.472   154.500   34.011500   17.695375  30.309000     19.312500   ./flapped_full_test/flapped_structure_-6.pdb
0  261.671  133.535  235.147   146.584   32.708875   16.691875  29.393375     18.323000    ./flapped_full_test/flapped_structure_7.pdb
0  280.390  145.607  249.428   158.002   35.048750   18.200875  31.178500     19.750250    ./flapped_full_test/flapped_structure_8.pdb
0  252.074  125.013  229.903   138.007   31.509250   15.626625  28.737875     17.250875    ./flapped_full_test/flapped_structure_6.pdb
0  213.274  106.501  192.799   113.938   26.659250   13.312625  24.099875     14.242250    ./flapped_full_test/flapped_structure_3.pdb
0  312.420  151.203  276.961   157.914   39.052500   18.900375  34.620125     19.739250    ./flapped_full_test/flapped_structure_0.pdb
0  236.716  118.067  213.425   126.871   29.589500   14.758375  26.678125     15.858875    ./flapped_full_test/flapped_structure_5.pdb
0  224.186  112.187  202.630   119.764   28.023250   14.023375  25.328750     14.970500   ./flapped_full_test/flapped_structure_-9.pdb
0  304.985  152.715  272.533   164.202   38.123125   19.089375  34.066625     20.525250    ./flapped_full_test/flapped_structure_9.pdb
0  311.211  152.594  276.498   160.619   38.901375   19.074250  34.562250     20.077375   ./flapped_full_test/flapped_structure_-2.pdb
0  249.548  123.984  226.208   134.901   31.193500   15.498000  28.276000     16.862625   ./flapped_full_test/flapped_structure_-8.pdb
Error: (mutmodel) 
790 hydrogens were added.
Traceback (most recent call last):
  File "flapper_full.py", line 355, in <module>
	train_df, test_df=maḱe_train_test_df(actual_pdbs, save_path, max_angle_diff)
  File "flapper_full.py", line 281, in maḱe_train_test_df
	df_final, df_file_angles=merge_features(save_path)
  File "flapper_full.py", line 237, in merge_features
	energy_df=get_ecalc(save_path)
  File "flapper_full.py", line 165, in get_ecalc
	ecalc_out=subprocess.check_output(command, shell=True)
  File "/usr/lib64/python3.6/subprocess.py", line 356, in check_output
	**kwargs).stdout
  File "/usr/lib64/python3.6/subprocess.py", line 438, in run
	output=stdout, stderr=stderr)
subprocess.CalledProcessError: Command './libs/ecalc-master/src/ecalc -p ./flapped_full_test/flapped_structure_-8.pdh' died with <Signals.SIGSEGV: 11>.
[lilian@serv1 html_lilian]$"""



'''5TIH_1
Too many chains to add disulphide topology
Failure in scanning for disulphides
1026 hydrogens were added.
Too many chains to add disulphide topology
Failure in scanning for disulphides'''

'''
6B0A_1
728 hydrogens were added.
Too many chains to add disulphide topology
Failure in scanning for disulphides
731 hydrogens were added.'''
