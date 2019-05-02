import pandas as pd



def mo(s):
    Keyword={}
    Keyword["DTA"]=["hat","dta"]
    Keyword["DTC"]=["cact" , "enspm" , "dtc"]
    Keyword["DTH"]=["harbinger" , "pif" , "dth","tourist"]
    Keyword["DTM"]=["mudr" ,"mule", "mu" , "mutator" , "dtm"]
    Keyword["DTT"]=["tc1" , "mariner" , "dtt","stow"]
    for fam in ["DTA","DTC","DTH","DTM","DTT"]:
        if any(c in s.lower() for c in Keyword[fam]):
            return s+"_"+fam

def modify(file):
    f=pd.read_table(file,header=None,sep="\t")
    f[1]=f[1].apply(lambda x : mo(x))
    f=f.to_csv(file,header=None,index=None,sep="\t")


for chr in range(1,13):
    modify("%s_IRF_withhomology_80_unique"%(chr))
