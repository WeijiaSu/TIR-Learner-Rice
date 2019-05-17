import pandas as pd
import os


desired_width=500

pd.set_option('display.width', desired_width)

pd.set_option('display.max_columns',20)

os.chdir("/Users/weijiasu/Dropbox/Research/BioinformaticsProject/Rice")

def split(s):
    return s.split("_")[1]+"_"+s.split("_")[2]+"_"+s.split("_")[3]


def FlankingBlast(file):
    f=pd.read_table(file,header=None,sep="\t",comment="#")
    f=f.loc[f[11]>=50]
    f["coor"]=f[0].apply(lambda x: split(x))
    coor=list(f["coor"])
    return list(set(coor))

def TArepeats(s):
    t=s.upper().count("T")
    a=s.upper().count("A")
    ta=t+a
    if (ta>len(s)*0.7):
        return True
    else:
        return False

def tir(s):
    return s.split("_")[0].split(":")[1]

def tsd(s):
    return s.split("_")[3].split(":")[1]

def RemovePutativeDTA(file,output):
    f=pd.read_table(file,header=None,sep="\t")
    print(f)
    f["tir"]=f[8].apply(lambda x:tir(x))
    print(f.shape)
    n=f.shape[0]
    f=f.loc[(f["tir"]!="CGGACAGTCCG") & (f["tir"]!="CGGACTGTCCG") & (f["tir"]!="GGACAGTCCGG") & (f["tir"]!="GGACTGTCCGG")]
    print(f.shape)
    print(f.shape[0]/n)
    f.to_csv(output,header=None,index=None,sep="\t")

############################################################################## Remove Entries with Homology in Flanking Sequences ##################################
print("############################################################################## Remove Putative DTA ##################################")

for genome in ["."]:
    for dataset in ["result_module_1", "result_module_1","result_module_1"]:
        file = RemovePutativeDTA("%s/%s/%s.gff3"%(genome,dataset,dataset),"%s/%s/%s.gff3" % (genome, dataset, dataset) )
        print("Removing Putative DTA in %s"%(genome))
print("########################################################################## Finished : Remove Putative DTA ##############################")

######################################################################################################################################################################

######################################################################################################################################################################


def renameSource(file):
    f = pd.read_table(file, header=None, sep="\t")
    f[1]="IRF_HC"
    f.to_csv(file,header=None,index=None,sep="\t")
    return f

def removeDupinSingle(file):
    f=pd.read_table(file,header=None,sep="\t")
    f=f.sort_values([0,3,4],ascending=[True,True,True])
    f=f.drop_duplicates([0,3,4],keep="last")
    return f

def RemoveTA(f):
    f["TIR"]=f[8].apply(lambda x:tir(x))
    f["TA"]=f["TIR"].apply(lambda x:TArepeats(x))
    sub=f.loc[f["TA"]==False]
    return sub[[0,1,2,3,4,5,6,7,8]]

def combineHBIRFhb(f1,f2):
    f=f1.append(f2,ignore_index=True)
    f=f.sort_values([0,3,4],ascending=[True,True,True])
    f = f.drop_duplicates([0, 3, 4], keep="first")
    return f

def combineAll(f1,f2,out):
    f = f1.append(f2, ignore_index=True)
    f.to_csv(out,header=None,index=None,sep="\t")
    return f

############################################################################## Process and Combining three gff files ################################################
print("############################################################################## Processing and Combining three gff files ##################################")

for genome in ["."]:
    os.chdir("%s/combine"%(genome))
    cp="cp ../result_module_1/result_module_1.gff3 ."
    os.system(cp)
    cp="cp ../result_module_2/result_module_2.gff3 ."
    os.system(cp)
    cp="cp ../result_module_3/result_module_3.gff3 ."
    os.system(cp)
    renameSource("result_module_2.gff3")
    f_hb=removeDupinSingle("result_module_1.gff3")
    f_irfref=removeDupinSingle("result_module_2.gff3")
    f_sb=removeDupinSingle("result_module_3.gff3")
    f_irfref=RemoveTA(f_irfref)
    f_sb = RemoveTA(f_sb)
    comHBIRF=combineHBIRFhb(f_hb,f_irfref)
    comAll=combineAll(comHBIRF, f_sb, "%combined.gff3")
print("#######################################################################Finished: Processing and Combining three gff files ##################################")
######################################################################################################################################################################

def ProcessGff(file,output):
    f=pd.read_table(file,header=None,sep="\t")
    f["pri"]=0
    mask=f[2]=="DTM"
    f.loc[mask, "pri"] = 1
    mask=f[2]=="DTC"
    f.loc[mask, "pri"] = 2
    mask = f[2] == "DTA"
    f.loc[mask, "pri"] = 3
    mask = f[2] == "DTT"
    f.loc[mask, "pri"] = 4
    mask = f[2] == "DTH"
    f.loc[mask, "pri"] = 5
    f["copy3"] = f[3]
    f["copy4"] = f[4]
    f.copy3 = f.copy3.shift(1).fillna(value=0).astype("int64")
    f.copy4 = f.copy4.shift(1).fillna(value=0).astype("int64")
    mask = ((f[3] == f["copy3"]) & (f[4] == f["copy4"]))
    f.loc[mask, 1] = "Both"
    f = f.sort_values([0, 3, 4, "pri",1], ascending=[True,True, True, True, True])
    f = f.drop_duplicates([0, 3, 4,"pri"], keep="first")
    f = f.sort_values([0, 3, 4,1, "pri"], ascending=[True,True,True, True, True])
    f=f.drop_duplicates([0,3,4],keep="first")
    f["length"]=f[4]-f[3]+1
    f=f[[0,1,2,3,4,5,6,7,8,"pri","length"]]
    f.to_csv(output,header=None,index=None,sep="\t")

############################################################################## Preparing for Removing Overlaps ################################################
print("############################################################################## Preparing for Removing Overlaps ##################################")
for genome in ["."]:
    os.chdir("%s/combine"%(genome))
    ProcessGff("combined.gff3", "combine_all_process.gff3")
    f = pd.read_table("combine_all_process.gff3", header=None, sep="\t")
    for i in range(1,13):
        chr=f.loc[f[0]==i]
        chr.to_csv("combine_all_process_%s.gff3"%(i),header=None,index=None,sep="\t")
print("############################################################################## Finished: Removing Overlaps ##################################")
######################################################################################################################################################################

def splitInfor(s,i):
    l=s.split("_")
    return float(l[i])


def ProcessSelect(file):
    f = pd.read_table(file, header=None, sep="\t")
    f = f.sort_values([0,3, 4, 9,10], ascending=[True, True, True, True,True])

    f["copy3"]=f[3]
    f["copy4"]=f[4]
    f["copy10"]=f[10]
    f["TIRp"] = f[8].apply(lambda x: splitInfor(x, 2))
    f["TSDp"] = f[8].apply(lambda x: splitInfor(x, 5))
    f["copy9"] = f[9]
    f["copyTIRp"] = f["TIRp"]
    f["copyTSDp"] = f["TSDp"]
    f.copyTIRp = f.TIRp.shift(1).fillna(value=0)
    f.copyTSDp = f.TSDp.shift(1).fillna(value=0)
    f.copy3=f.copy3.shift(1).fillna(value=0).astype("int64")
    f.copy4 = f.copy4.shift(1).fillna(value=0).astype("int64")
    f.copy10 = f.copy10.shift(1).fillna(value=0).astype("int64")
    f.copy9 = f.copy9.shift(1).fillna(value=0).astype("int64")
    f.copyTIRp = f.TIRp.shift(1).fillna(value=0)
    f.copyTSDp = f.TSDp.shift(1).fillna(value=0)
    f["3_copy4"]=f[3]-f["copy4"]
    f["4_copy4"]=f[4]-f["copy4"]
    f["len-len"]=f[10]-f["copy10"]
    f["pri-pri"]=f[9]-f["copy9"]
    f["tir-tir"]=f["TIRp"]-f["copyTIRp"]
    f["tsd-tsd"]=f["TSDp"]-f["copyTSDp"]
    f=f.sort_values([0,3,4,1],ascending=[True,True,True,True])
    return f



def getRemoveList(f):
    removeList=[]
    overlap=(((f[3]-f["copy3"]<=30) & (f[3]-f["copy3"]>=0)) | ((f["copy4"]-f[4]<=30) & (f["copy4"]-f[4]>=0) ))
    r1=f.loc[(f["3_copy4"]<=0)&(f["4_copy4"]<=0) & (overlap==True)]
    if r1.shape[0]==0:
        pass
    removeList=removeList+list(r1.index.values)
    r2=f.loc[(f["3_copy4"]<=0)&(f["4_copy4"]>=0)]
    if r2.shape[0]==0:
        pass
    for index, row in r2.iterrows():
        if(row["pri-pri"]>0):
            removeList.append(index)
        elif(row["pri-pri"]<0):
            removeList.append(index-1)
        else:
            if(row["tir-tir"]<0):
                removeList.append(index)
            elif(row["tir-tir"]>0):
                removeList.append(index-1)
            else:
                if (row["tsd-tsd"] < 0):
                    removeList.append(index)
                elif (row["tsd-tsd"] > 0):
                    removeList.append(index-1)
                else:
                    if(row["len-len"] <= 0):
                        removeList.append(index)
                    else:
                        removeList.append(index - 1)
    return removeList


def CheckOverlap(file):
    re_open = pd.read_table(file, header=None, sep="\t")
    re_open["copy3"] = re_open[3]
    re_open["copy4"] = re_open[4]
    re_open["copy10"] = re_open[10]
    re_open.copy3 = re_open.copy3.shift(1).fillna(value=0).astype("int64")
    re_open.copy4 = re_open.copy4.shift(1).fillna(value=0).astype("int64")
    re_open.copy10 = re_open.copy10.shift(1).fillna(value=0).astype("int64")
    re_open["3_copy4"] = re_open[3] - re_open["copy4"]
    re_open["4_copy4"] = re_open[4] - re_open["copy4"]
    re_open["len-len"] = re_open[10] - re_open["copy10"]
    overlap = (((re_open[3]-re_open["copy3"]<=30) & (re_open[3]-re_open["copy3"]>=0)) | ((re_open["copy4"]-re_open[4]<=30) & (re_open["copy4"]-re_open[4]>=0) ))
    r1 = re_open.loc[(re_open["3_copy4"] <= 0) & (re_open["4_copy4"] <= 0) & (overlap==True)]
    r2 = re_open.loc[(re_open["3_copy4"] <= 0) & (re_open["4_copy4"] >= 0)]
    if (r1.shape[0]!=0 or r2.shape[0]!=0):
        return False
    else:
        return True

def deleteOverlap(file,output):
    f=ProcessSelect(file)
    l=getRemoveList(f)
    if len(l)!=0:
        newf = f.drop(f.index[l])
        newf = newf[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]
        newf.to_csv(output, header=None, index=None, sep="\t")
        deleteOverlap(output, output)
    else:
        newf = f[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]
        newf.to_csv(output, header=None, index=None, sep="\t")
        return newf

############################################################################## Removing Overlaps ################################################
print("############################################################################## Removing Overlaps ##################################")

for genome in ["."]:
    os.chdir("%s/combine"%(genome))
    for i in range(1, 13):
        newf = deleteOverlap("combine_all_process_%s.gff3"%(i), "combine_all_process_%s_Oremoved.gff3.txt" % (i))
    cat="cat *_combine_all_process_*_Oremoved.gff3.txt > FinalAnn.gff3"
    os.system(cat)
    final=pd.read_table("FinalAnn.gff3",header=0,sep=",")
    print("Final Size:")
    print(final.shape)
print("############################################################################## Finished: Removing Overlaps ##################################")
######################################################################################################################################################################


def RemoveOerlap(file,outname):
    removeList=[]
    f=pd.read_table(file,header=None,sep="\t")
    f = f.sort_values([0, 3, 4, 9, 10], ascending=[True, True, True, True, True])

    f["copy3"] = f[3]
    f["copy4"] = f[4]
    f["copy10"] = f[10]
    f.copy3 = f.copy3.shift(1).fillna(value=0).astype("int64")
    f.copy4 = f.copy4.shift(1).fillna(value=0).astype("int64")
    f.copy10 = f.copy10.shift(1).fillna(value=0).astype("int64")
    f["3_copy4"] = f[3] - f["copy4"]
    f["4_copy4"] = f[4] - f["copy4"]
    f["len-len"] = f[10] - f["copy10"]
    f["copy_0"]=f[0]
    f["0_copy0"]=f[0]-f["copy_0"]

    f = f.sort_values([0, 3, 4, 1], ascending=[True, True, True, True])
    sub1=f.loc[(f["3_copy4"]<=0)&(f["len-len"]<=0)&(f["0_copy0"]==0)]
    removeList = removeList + list(sub1.index.values)
    sub2=f.loc[(f["3_copy4"]<=0)&(f["len-len"]>0)&(f["0_copy0"]==0)]
    removeList = removeList + [i-1 for i in list(sub2.index.values)]
    f=f.drop(f.index[removeList])
    f=f.drop(["copy3","copy4","copy10","3_copy4","4_copy4","len-len"],axis=1)
    f.to_csv(outname,header=None,index=None,sep="\t")
    return removeList

############################################################################## Deleting internal copies ################################################
print("############################################################################## Deleting internal copies ##################################")
for genome in ["."]:
    os.chdir("%s/combine"%(genome))
    RemoveOerlap("FinalAnn.gff3", "FinalAnn_noInternal.gff3")
    f=pd.read_table("FinalAnn_noInternal.gff3"),header=None,sep="\t")

print("############################################################################## Finished: Deleting internal copies ##################################")
######################################################################################################################################################################
