#!/usr/bin/env python
import csv
import re
from collections import *
import sys

sys.path.insert(0,'./')
from utils.critical_concentrations_who import *

def recursive_defaultdict():
    return defaultdict(recursive_defaultdict)

details={
    "antb":["INH","RIF","RFB","EMB","STR","ETA","CIP","CYS","CAP","KAN","OFX","PAS","PZA","AMI","MOXI","PRO","CLO","LEVO","CLAR","GATI","AMOXCLAV","LIN"],
    "abbrev_antb":{
        "INH":"ISONIAZID",
        "RIF":"RIFAMPICIN",
        "RFB":"RIFABUTIN",
        "EMB":"ETHAMBUTOL",
        "STR":"STREPTOMYCIN",
        "ETA":"ETHIONAMIDE",
        "CIP":"CIPROFLOXACIN",
        "CYS":"CYCLOSERINE",
        "CAP":"CAPREOMYCIN",
        "KAN":"KANAMYCIN",
        "OFX":"OFLOXACIN",
        "PAS":"PARA-AMINOSALICYLIC_ACID",
        "PZA":"PYRAZINAMIDE",
        "AMI":"AMIKACIN",
        "MOXI":"MOXIFLOXACIN",
        "PRO":"PROTHIONAMIDE",
        "CLO":"CLOFAZIMINE",
        "LEVO":"LEVOFLOXACIN",
        "CLAR":"CLARITHROMYCIN",
        "GATI":"GATIFLOXACIN",
        "AMOXCLAV":"AMOXICILLIN_CLAVULANATE",
        "LIN":"LINEZOLID",
    }
}

def resolve_aliases(table_mic_data: str, table_strain_identification: str):
    """Takes a """
    dict_aliases=recursive_defaultdict()
    dict_corresp=recursive_defaultdict()
    # I get the IDs of the internal strains present on the table_strain_identification
    list_ids_table_strain_identification=[]
    list_ids_table_mic_data=[]
    inp = csv.DictReader(open(table_strain_identification,"r"),delimiter="\t")
    for i in inp:
        if i["internal_xref"]!="":
            list_ids_table_strain_identification.append(i["internal_xref"])
    inp = csv.DictReader(open(table_mic_data,"r"),delimiter="\t")
    count_entries=0
    for i in inp:
        count_entries+=1
        prim_id=i["ID"].replace("-","").replace("_","").replace(" ","")
        list_ids_table_mic_data.append(prim_id)
        alt_id=i["Alt ID"].replace("-","").replace("_","").replace(" ","")
        dict_aliases[prim_id]=prim_id
        dict_aliases[alt_id]=prim_id
    for id_table_identification in list_ids_table_strain_identification:
        if id_table_identification in dict_aliases:
            dict_corresp[dict_aliases[id_table_identification]]=id_table_identification
    intersection=set(list_ids_table_strain_identification).intersection(dict_corresp.keys())
    not_matched=set(list_ids_table_mic_data)-set(dict_corresp.keys())
    print("[INFO] I tried to resolve the aliases. Here are the stats:")
    print("  * I found a match for {} isolates ({:.1f}%)".format(len(intersection),len(intersection)/count_entries*100))
    print("  * I did NOT find a match for {} isolates ({:.1f}%)".format(len(not_matched),len(not_matched)/count_entries*100))
    missing_from_genomic_data=set(list_ids_table_mic_data) -set(not_matched) - set(intersection)
    print("  * Isolates we did not pushed yet through megapipe: {} ({:.1f}%)".format(len(missing_from_genomic_data),len(missing_from_genomic_data)/count_entries*100))
    return(dict_corresp, intersection, not_matched,missing_from_genomic_data)

def build_table(table_mic_data: str,data: dict,ids_ok: list, dict_corresp: defaultdict) -> list:
    data_to_write=[]
    inp = csv.DictReader(open(table_mic_data,"r"),delimiter="\t")
    entries_skipped=0
    num_entries=0
    for i in inp:
        for j in details["antb"]:
        #entry: <xref> <tags> <antb> <exp_type> <conc> <conc_units> <mic_summary> <mic_conc_tested> <res_class> <method> <medium> <device> <doi> <who_compliant>
            num_entries+=1
            entry=[""]*13
            # I define the antibiotic
            entry[2]=details["abbrev_antb"][j]
            # I define the value of the MIC
            value_mic=i[j]
            if value_mic=="NA":
                continue
            elif value_mic in ["r","s"]:
                entry[3]="UNKNOWN"
                entry[8]=value_mic.upper()
            else:
                entry[6]=value_mic
                # These are MIC data
                entry[3]="MIC"
                # conc (units)
                entry[5]="μg/ml"
            # I define the ID of the isolate
            if dict_corresp[i["ID"].replace("-","").replace("_","").replace(" ","")]:
                id_isolate=dict_corresp[i["ID"].replace("-","").replace("_","").replace(" ","")]
            else:
                entries_skipped+=1
                continue
            entry[0]=id_isolate
            # I define the tag
            if(i["Source Lab"]=="MSLI"):
                tag="pools"
                if entry[2] in ["ISONIAZID","RIFAMPICIN","ETHAMBUTOL","STREPTOMYCIN","KANAMYCIN","CAPREOMYCIN","ETHIONAMIDE","CYCLOSERINE","PARA-AMINOSALICYLIC_ACID","AMIKACIN","LEVOFLOXACIN","OFLOXACIN","CIPROFLOXACIN"]:
                    entry[9]="indirect proportion method"
                    entry[10]="Middlebrook 7H10"
                elif entry[2]=="PYRAZINAMIDE":
                    entry[9]="BACTEC MGIT"
            elif (i["Source Lab"]=="RIVM"):
                tag="pools"
                entry[9]="BACTEC MGIT"
            # SES = Socios en Salud
            elif(i["Source Lab"]=="SES"):
                tag="cetr"
                if entry[2]=="PYRAZINAMIDE":
                    entry[9]="BACTEC MGIT"
                else:
                    entry[10]="Middlebrook 7H10"
            else:
                if entry[2] in ["ISONIAZID","RIFAMPICIN","ETHAMBUTOL","STREPTOMYCIN","PARA-AMINOSALICYLIC_ACID"]:
                    entry[10]="Löwenstein-Jensen"
                elif entry[2] in ["OFLOXACIN","KANAMYCIN","CAPREOMYCIN"]:
                    entry[10]="Middlebrook 7H11"
                entry[9]="indirect proportion method"
                tag="tdr"
            entry[1]=tag
            data_to_write.append(entry)
    print("[INFO] {}/{} entries skipped ({:.1f}%)".format(entries_skipped,num_entries,entries_skipped/num_entries*100))
    return(data_to_write)

def parse_rows_take_decisions(list_rows_table, file_out, csv_critical_conc):
    """Takes the list of lists generated by build table and generates a tsv file with the following fields: <biosample> <antibiotic> <res_class> <tags>"""
    # I load the critical concentrations data
    thresholds=get_who_thresholds(csv_critical_conc)
    # I define the media conversions
    media_conv={"Middlebrook 7H10":"m7h10","Middlebrook 7H11":"m7h11","Löwenstein-Jensen":"lj"}
    method_conv={"BACTEC MGIT":"mgit960"}
    # Here I will store the list of the discarded entries
    list_discarded=[]
    count_all_entries=0
    count_discarded=0
    reasons_discarded={"no_media_method":0,"conc_tested_not_allows_decision":0,"new_regex_needed":0}
    interm=re.compile("([0-9].*[0-9]*)-([0-9].*[0-9]*)")
    great_eq_than=re.compile(">=([0-9]+.*[0-9]*)")
    great_than=re.compile(">([0-9]+.*[0-9]*)")
    less_than=re.compile("<=*([0-9]+.*[0-9]*)")
    numb=re.compile("(^[0-9]$)")
    with open(file_out,"w") as outf:
        for row in list_rows_table:
            count_all_entries+=1
            if row[3] in ["UNKNOWN","DST"]:
                outf.write("{}\t{}\t{}\t{}\n".format(row[0],row[2],row[8],row[1]))
            elif row[3] in ["MIC"]:
                res_interm=re.findall(interm,row[6])
                res_great_eq_than=re.findall(great_eq_than,row[6])
                res_great_than=re.findall(great_than,row[6])
                res_less_than=re.findall(less_than,row[6])
                res_numb=re.findall(numb,row[6])
                if res_interm:
                    #print("interm: {}".format(row[6]))
                    inf_lim=float((res_interm[0][0]))
                    sup_lim=float((res_interm[0][1]))
                    if row[10] in media_conv:
                        medium=media_conv[row[10]]
                    elif row[9] in method_conv:
                        medium=method_conv[row[9]]
                    else:
                        count_discarded+=1
                        list_discarded.append(row)
                        reasons_discarded["no_media_method"]+=1
                        continue
                    who_threshold=thresholds[row[2]][medium]
                    if inf_lim > who_threshold:
                        outf.write("{}\t{}\t{}\t{}\n".format(row[0],row[2],"R",row[1]))
                    elif sup_lim <= who_threshold:
                        outf.write("{}\t{}\t{}\t{}\n".format(row[0],row[2],"S",row[1]))
                    else:
                        count_discarded+=1
                        list_discarded.append(row)
                        reasons_discarded["conc_tested_not_allows_decision"]+=1
                elif res_great_eq_than:
                    #print("great_eq_than: {}".format(row[6]))
                    lim=float((res_great_eq_than[0]))
                    if row[10] in media_conv:
                        medium=media_conv[row[10]]
                    elif row[9] in method_conv:
                        medium=method_conv[row[9]]
                    else:
                        count_discarded+=1
                        list_discarded.append(row)
                        reasons_discarded["no_media_method"]+=1
                        continue
                    who_threshold=thresholds[row[2]][medium]
                    if lim > who_threshold:
                        outf.write("{}\t{}\t{}\t{}\n".format(row[0],row[2],"R",row[1]))
                    else:
                        count_discarded+=1
                        list_discarded.append(row)
                        reasons_discarded["conc_tested_not_allows_decision"]+=1
                elif res_great_than:
                    #print("great_than: {}".format(row[6]))
                    lim=float((res_great_than[0]))
                    if row[10] in media_conv:
                        medium=media_conv[row[10]]
                    elif row[9] in method_conv:
                        medium=method_conv[row[9]]
                    else:
                        count_discarded+=1
                        list_discarded.append(row)
                        reasons_discarded["no_media_method"]+=1
                        continue
                    who_threshold=thresholds[row[2]][medium]
                    if lim >= who_threshold:
                        outf.write("{}\t{}\t{}\t{}\n".format(row[0],row[2],"R",row[1]))
                    else:
                        count_discarded+=1
                        list_discarded.append(row)
                        reasons_discarded["conc_tested_not_allows_decision"]+=1
                elif res_less_than:
                    #print("less_than: {}".format(row[6]))
                    lim=float((res_less_than[0]))
                    if row[10] in media_conv:
                        medium=media_conv[row[10]]
                    elif row[9] in method_conv:
                        medium=method_conv[row[9]]
                    else:
                        count_discarded+=1
                        list_discarded.append(row)
                        reasons_discarded["no_media_method"]+=1
                        continue
                    who_threshold=thresholds[row[2]][medium]
                    if lim <= who_threshold:
                        outf.write("{}\t{}\t{}\t{}\n".format(row[0],row[2],"S",row[1]))
                    else:
                        count_discarded+=1
                        list_discarded.append(row)
                        reasons_discarded["conc_tested_not_allows_decision"]+=1
                elif res_numb:
                    #print("numb: {}".format(row[6]))
                    lim=float((res_numb[0]))
                    if row[10] in media_conv:
                        medium=media_conv[row[10]]
                    elif row[9] in method_conv:
                        medium=method_conv[row[9]]
                    else:
                        count_discarded+=1
                        list_discarded.append(row)
                        reasons_discarded["no_media_method"]+=1
                        continue
                    who_threshold=thresholds[row[2]][medium]
                    if lim > who_threshold:
                        outf.write("{}\t{}\t{}\t{}\n".format(row[0],row[2],"R",row[1]))
                    elif lim <= who_threshold:
                        outf.write("{}\t{}\t{}\t{}\n".format(row[0],row[2],"S",row[1]))
                    else:
                        count_discarded+=1
                        list_discarded.append(row)
                        reasons_discarded["conc_tested_not_allows_decision"]+=1
                else:
                    reasons_discarded["new_regex_needed"]+=1
                    list_discarded.append(row)
                    print("unknown regex: {} // {} {} {} {}".format(row[6],row[0],row[2],row[8],row[1]))
    print("[INFO] Total number of entries: {}".format(count_all_entries))
    print("[INFO] {} entries were discarded ({:.1f}%)".format(count_discarded,count_discarded/count_all_entries*100))
    for reason in reasons_discarded:
        print("* {}: {}".format(reason,reasons_discarded[reason]))
    return(list_discarded)

(d_corresp,list_ids_ok, list_ids_not_matched,list_ids_missing_from_genomic_data)=resolve_aliases("../jupyter/sources/resistance_data/Resistance_data_internal_strains-20180514.tsv", "/home/lf61/mfarhat/rollingDB/tables/table_strain_identification_data6.tsv")
dat=build_table("../jupyter/sources/resistance_data/Resistance_data_internal_strains-20180514.tsv",details,list_ids_ok,d_corresp)
discarded=parse_rows_take_decisions(dat,"./internal/internal.res","/n/data1/hms/dbmi/farhat/rollingDB/tables/jupyter/sources/critical_concentrations/criticalConcentrations.csv")
for entry in discarded:
    print(entry)


