#!/usr/bin/env python
import csv
import glob
import os
import re

details={
    "def_res_class":{"Susceptible":"S","Resistant":"R"}
}

def build_table(table_corresp,dir_tables_res_data,data):
    dict_corresp={}
    inp = csv.DictReader(open(table_corresp,"r"),delimiter=",",quotechar='"')
    ids_without_biosample={}
    for i in inp:
        biosample=i["BioSample Accession"]
        if biosample=="":
            ids_without_biosample[i["Genome ID"]]=1
        else:
            dict_corresp[i["Genome ID"]]=biosample
    print("* [INFO] There are {} isolates with a biosample.".format(len(dict_corresp.keys())))
    print("* [WARNING] There are {} isolates without a biosample. They will be discarded.".format(len(ids_without_biosample.keys())))
    files_to_parse=glob.glob(dir_tables_res_data+"/PATRIC_amr*")
    data_to_write=[]
    count_total_entries=0
    count_discarded=0
    for current_file in files_to_parse:
        inp = csv.DictReader(open(current_file,"r"),delimiter=",")
        for i in inp:
            count_total_entries+=1
            #entry: <xref> <tags> <antb> <exp_type> <conc> <conc_units> <mic_summary> <mic_conc_tested> <res_class> <method> <medium> <device> <doi> <who_compliant>
            entry=[""]*13
            if not i["Genome ID"] in dict_corresp:
                count_discarded+=1
                print(i)
                continue
            if i["Resistant Phenotype"]=="" or i["Resistant Phenotype"]=="Intermediate":
                continue
            #I define the entries that are MICs
            if i["Laboratory Typing Method"]=="MIC":
                entry[0]=dict_corresp[i["Genome ID"]]
                entry[1]="patric"
                entry[2]=i["Antibiotic"].upper()
                if entry[2]=="RIFAMPIN":
                    entry[2]="RIFAMPICIN"
                elif entry[2]=="PARA-AMINOSALICYLIC ACID":
                    entry[2]="PARA-AMINOSALICYLIC_ACID"
                entry[3]="MIC"
                entry[5]=i["Measurement Unit"]
                entry[6]=i["Measurement"]
            if i["Laboratory Typing Method"]=="MGIT":
                entry[0]=dict_corresp[i["Genome ID"]]
                entry[1]="patric"
                entry[2]=i["Antibiotic"].upper()
                if entry[2]=="RIFAMPIN":
                    entry[2]="RIFAMPICIN"
                elif entry[2]=="PARA-AMINOSALICYLIC ACID":
                    entry[2]="PARA-AMINOSALICYLIC_ACID"
                entry[3]="DST"
                conc_tested=i["Laboratory Typing Method Version"]
                if not conc_tested=="":
                    explode=conc_tested.split(" ")
                    entry[4]=explode[0]
                    entry[5]=explode[1]
                entry[8]=data["def_res_class"][i["Resistant Phenotype"]]
            if i["Laboratory Typing Method"]=="agar proportion method":
                entry[0]=dict_corresp[i["Genome ID"]]
                entry[1]="patric"
                entry[2]=i["Antibiotic"].upper()
                if entry[2]=="RIFAMPIN":
                    entry[2]="RIFAMPICIN"
                elif entry[2]=="PARA-AMINOSALICYLIC ACID":
                    entry[2]="PARA-AMINOSALICYLIC_ACID"
                entry[3]="DST"
                conc_tested=i["Laboratory Typing Method Version"]
                ugmL=re.compile("ug/ml$")
                if ugmL.search(conc_tested):
                    entry[5]="μg/mL"
                    entry[4]=conc_tested.replace("ug/ml","")
                else:
                    entry[4]=conc_tested
                entry[8]=data["def_res_class"][i["Resistant Phenotype"]]
                entry[9]="agar proportion method"
            if i["Laboratory Typing Method"]=="LJM":
                entry[0]=dict_corresp[i["Genome ID"]]
                entry[1]="patric"
                entry[2]=i["Antibiotic"].upper()
                if entry[2]=="RIFAMPIN":
                    entry[2]="RIFAMPICIN"
                elif entry[2]=="PARA-AMINOSALICYLIC ACID":
                    entry[2]="PARA-AMINOSALICYLIC_ACID"
                entry[3]="DST"
                conc_tested=i["Laboratory Typing Method Version"]
                ugmL=re.compile("ug/ml$")
                if ugmL.search(conc_tested):
                    entry[5]="μg/mL"
                    entry[4]=conc_tested.replace("ug/ml","")
                else:
                    entry[4]=conc_tested
                entry[8]=data["def_res_class"][i["Resistant Phenotype"]]
                entry[10]="Löwenstein-Jensen"
            if i["Laboratory Typing Method"]=="":
                entry[0]=dict_corresp[i["Genome ID"]]
                entry[1]="patric"
                entry[2]=i["Antibiotic"].upper()
                if entry[2]=="RIFAMPIN":
                    entry[2]="RIFAMPICIN"
                elif entry[2]=="PARA-AMINOSALICYLIC ACID":
                    entry[2]="PARA-AMINOSALICYLIC_ACID"
                entry[3]="UNKNOWN"
                entry[8]=data["def_res_class"][i["Resistant Phenotype"]]
            if i["Laboratory Typing Method"]=="Computational Prediction":
                continue
            data_to_write.append(entry)
    print("* [INFO] total_count_entries={}; discarded={}".format(count_total_entries,count_discarded))
    return(data_to_write)

def parse_rows_take_decisions(list_rows_table, file_out):
    """Takes the list of lists generated by build table and generates a tsv file with the following fields: <biosample> <antibiotic> <res_class> <tags>"""
    with open(file_out,"w") as outf:
        for row in list_rows_table:
            if row[3] in ["UNKNOWN","DST"]:
                outf.write("{}\t{}\t{}\t{}\n".format(row[0],row[2],row[8],row[1]))

rows=build_table("/home/lf61/mfarhat/rollingDB/tables/jupyter/sources/resistance_data/patric/PATRIC_genome.csv","/home/lf61/mfarhat/rollingDB/tables/jupyter/sources/resistance_data/patric/",details)
parse_rows_take_decisions(rows, "./patric/patric.res")

