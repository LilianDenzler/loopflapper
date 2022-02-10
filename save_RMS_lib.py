#!/usr/bin/env python
#python3 save_RMS ~/git/CDR_H3_quality/RMS_calc ~/sync_project/actual_PDBs_NR ~/sync_project/abYmod_structures ~/sync_project/Feature

import sys
import os
import os.path
import csv
import pandas as pd
import subprocess
import numpy as np
import glob


def get_file_paths(RMSD_mode):
    qualiloop_path=os.path.abspath(os.path.join((os.path.dirname(__file__)), os.pardir))
    try:
        profit_path=os.path.join(qualiloop_path,"bin/ProFit_V3.3/profit")
    except:
        profit_path=glob.glob('**/profit', recursive=True)
        if len(profit_path)==0:
            raise ValueError('no profit file was found')
        elif len(profit_path)>=1:
            profit_path=profit_path[0]

    try:
        RMS_calc_path=os.path.join(qualiloop_path,"bin/RMS_calc"+RMSD_mode)
    except:
        RMS_calc_path=glob.glob('**/RMS_calc'+RMSD_mode, recursive=True)
        if len(RMS_calc_path)==0:
            raise ValueError('no RMS_calc file was found')
        elif len(RMS_calc_path)>=1:
            RMS_calc_path=RMS_calc_path[0]
    try:
        RMS_feature_path=os.path.join(qualiloop_path,"bin/data/RMS_feature_csv.csv")
    except:
        RMS_feature_path=glob.glob('**/data/RMS_feature_csv.csv', recursive=True)
        if len(RMS_feature_path)==0:
            RMS_feature_path=None
            print('no RMS_feature_path file was found')
        elif len(RMS_feature_path)>=1:
            RMS_feature_path=RMS_feature_path[0]

    return profit_path, RMS_calc_path, RMS_feature_path

def pass_Profit_commands(file, actual_file, mode):
    data=[]
    filesize = os.path.getsize(file)
    if filesize == 0:
        print(file, "empty")
        return pd.DataFrame([[np.nan]*len(mode)], columns=mode)
    for RMSD_mode in mode:
        profit_path="profit"
        RMS_calc=os.path.join(os.path.dirname(__file__),"RMS_calc"+RMSD_mode)
        command="{} -f {} -h {} {}".format(profit_path,RMS_calc, actual_file, file)
        profit_out=subprocess.check_output(command, shell=True)
        
        profit_out=str(profit_out, 'utf-8') 
        print(profit_out)
        list_write=[]
        for line in profit_out.split('\n'):
            columns=line.split("    ")
            if len(columns)< 2 and 'RMS:' in line:
                line=line.replace(" ","")
                line=line.replace("\n","")
                list_write.append(float(line.split(":")[1]))
        if len(list_write)==2:
            list_write=[list_write[1]]
        elif len(list_write)==3 and RMSD_mode=="global_AA":
            list_write=[list_write[1]]
        elif len(list_write)!=1:
            print("EERRROORR")
            list_write=[np.nan]
            #because if the RMSD_mode is global_CA the FIT command prints the RMSD for the whole
            #of the fitted structure. The second RMSD value printed is the RMSd over the loop only.
        data.extend(list_write)
    RMSD_num_df = pd.DataFrame([data], columns=mode)
    return (RMSD_num_df)


def run_num(file, actual_file, mode):
   
    RMSD_num_df=pass_Profit_commands(file, actual_file, mode)

    return RMSD_num_df
    

def RMSD_binary(RMSD_mode, file, actual_file, threshold, RMSD_num_df=None):
    columns=[RMSD_mode+"_bin"+str(threshold)]
    data=[]
    input(RMSD_mode)
    print("RMSD_num_df",RMSD_num_df)
    input()
    if RMSD_num_df is not None:
        RMSD_num_col=RMSD_num_df[RMSD_mode].tolist()
        input(RMSD_num_col)
    else: 
        RMSD_num_df=run_num(file, actual_file, [RMSD_mode])
        RMSD_num_col=RMSD_num_df.tolist()

    for a in RMSD_num_col:
        if a is None:
            print(a)
            a=np.nan
        a=float(a)
        if a>=threshold:
            a=1
        elif a==np.nan:
            a=np.nan
        else:
            a=0
        write=[a]
        data.append(write)
    df_bin = pd.DataFrame(data, columns=columns)
    return (df_bin)
    

def get_all_bins(RMSD_mode, file,actual_file, threshold_list, RMSD_num_df=None):
    df_all_bins=pd.DataFrame()
    if RMSD_num_df is None:
        filename=os.path.splitext(os.path.basename(file))[0]
        filename=filename.replace(".pdb", "")
        ID_df=pd.DataFrame()
        ID_df["ID"]=[filename]
    else:
        ID_df=RMSD_num_df["ID"].to_frame()
    for i in threshold_list:
        df_bin=RMSD_binary(RMSD_mode, file, actual_file, i, RMSD_num_df=RMSD_num_df)
        if df_bin.empty==True:
            continue
        else:
            #print(df_all_bins, df_bin,ID_df)
            df_all_bins = pd.concat([df_all_bins, df_bin], axis=1)
    #print(df_all_bins.shape, ID_df.shape)
    df_all_bins = pd.concat([df_all_bins,ID_df], axis=1)
        
    return df_all_bins


def RMSD_nom(RMSD_mode, threshlow,threshup, file, actual_file,counter=None, RMSD_num_df=None):
    if RMSD_num_df is None:
        filename=os.path.splitext(os.path.basename(file))[0]
        filename=filename.replace(".pdb", "")
        ID_df=pd.DataFrame()
        ID_df["ID"]=[filename]
    else:
        ID_df=RMSD_num_df["ID"].to_frame()
    if not counter:
        counter=""
    columns=[RMSD_mode+str("_nom")+str(counter)]
    data=[]
    values=np.arange(0,len(threshup))

    try:
        RMSD_col=RMSD_num_df[RMSD_mode].tolist()
    except:
        RMSD_num_df=run_num(file, actual_file, [RMSD_mode])
        RMSD_num_df=pd.concat([RMSD_num_df,ID_df], axis=1)
        try:
            RMSD_col= RMSD_num_df.tolist()
        except:
            input("ERR1")
            return pd.DataFrame()

    
    col=[]
    for a in RMSD_col:
        a=float(a)
        start=a
        for (threshold_low,threshold_up,val) in list(zip(threshlow,threshup,values)):
            if a is None or a==np.nan:
                a=np.nan
                break
            if a>threshold_low and a<=threshold_up:
                #print("YES",threshold_low, a, threshold_up)
                a=val
                break
            else:
                pass
                #print("NO",threshold_low, a, threshold_up)
        try:
            write=int(a)
        except:
            if a==np.nan:
                write=np.nan
        #print(write)
        if a==start:
            write=None

        col.append(write)
    
    df_nom = pd.DataFrame(col, columns=columns)
    df_nom=pd.concat([df_nom, ID_df],axis=1)
    #print("!!!!!!!!!")
    #print(df_nom)
    return (df_nom)
   