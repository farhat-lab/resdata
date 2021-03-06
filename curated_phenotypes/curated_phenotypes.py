#!/usr/bin/env python
import csv
import re
from collections import *
import sys
import requests

sys.path.insert(0,'./')
from utils.critical_concentrations_who import *

def recursive_defaultdict():
    return defaultdict(recursive_defaultdict)

details={
    "antb_conv":{
        "Oflaxacin":"OFLOXACIN",
        "Rifampin":"RIFAMPICIN",
        "Para-Aminosalicylate Sodium":"PARA-AMINOSALICYLIC_ACID",
        "Para-aminosalicylic acid":"PARA-AMINOSALICYLIC_ACID"
    },
    # I will use this dictionary to create the tags for Sonny's data
    "tag_conv":{
        "Casali et al., 2014":"casali2014",
        "Biek et al., 2012":"biek2012",
        "Blouin et al., 2012":"blouin2012",
        "Bryant et al., 2013":"bryant2013",
        "Chatterjee et al., 2017":"chatterjee2017",
        "Clark et al., 2013":"clark2013",
        "Cohen et al., 2015":"cohen2015",
        "Comas et al., 2013":"comas2013",
        "Gardy et al., 2011":"gardy2011",
        "Guerra-Assuncao et al., 2014":"guerra2014",
        "Lieberman et al., 2016":"lieberman2016",
        "Manson et al., 2017":"manson2017",
        "Merker et al., 2015":"merker2015",
        "Perdigao et al., 2014":"perdigao2014",
        "Phelan et al., 2016":"phelan2016",
        "TBARC.Belarus":"tbarc.belarus",
        "TBARC.CDRC":"tbarc.cdrc",
        "TBARC.India":"tbarc.india",
        "TBARC.Iran":"tbarc.iran",
        "TBARC.Moldova":"tbarc.moldova",
        "TBARC.MRC":"tbarc.mrc",
        "TBARC.Romania":"tbarc.romania",
        "TBARC.Sweden":"tbarc.sweden",
        "Walker et al., 2013":"walker2013",
        "Walker et al., 2015":"walker2015",
        "Winglee et al., 2015":"winglee2015",
        "Zhang et al., 2013":"zhang2013"
    }
}

def getNCBIIdType(my_id):
    #bioproject
    re_bp=re.compile("^PRJ")
    if(len(re_bp.findall(my_id))!=0):
        return "bioproject"
    #biosample
    re_bs=re.compile("^SAM")
    if(len(re_bs.findall(my_id))!=0):
        return "biosample"
    #run
    re_r=re.compile("^ERR")
    re_r2=re.compile("^SRR")
    if((len(re_r.findall(my_id))!=0) or (len(re_r2.findall(my_id))!=0)):
        return "sra"
    return "unknown"

def trace_request(**kw):
    request = "http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db={db}&rettype=runinfo&term={term}".format(**kw)
    response = requests.get(url=request,stream=True)
    response.encoding = 'utf-8'
    return response

def build_table(table: str,data: dict) -> list:
    data_to_write=[]
    list_discarded=[]
    reasons_discarded={"no_antb":0,"no_res_class":0}
    inp = csv.DictReader(open(table,"r"),delimiter="\t")
    entries_discarded=0
    entries_parsed=0
    num_entries=0
    for i in inp:
    #entry: <xref> <tags> <antb> <exp_type> <conc> <conc_units> <mic_summary> <mic_conc_tested> <res_class> <method> <medium> <device> <doi> <who_compliant>
        num_entries+=1
        entry=[""]*13
        # I define the antibiotic
        if i["MSDRUG"] in ["","ND"]:
            entries_discarded+=1
            reasons_discarded["no_antb"]+=1
            list_discarded.append(i)
            continue
        elif i["MSDRUG"] in data["antb_conv"]:
            entry[2]=data["antb_conv"][i["MSDRUG"]]
        else:
            entry[2]=i["MSDRUG"].upper()
        # I define the tag
        tag=data["tag_conv"][i["Study"]]
        entry[1]=tag
        # main part to define the table
        if tag in ["blouin2012","bryant2013","chatterjee2017","clark2013","cohen2015","gardy2011","lieberman2016","merker2015","tbarc.belarus","tbarc.moldova","tbarc.romania","tbarc.sweden","walker2013","walker2015","zhang2013"]:
            res_class=i["Sensitive (S), Resistant (R), Intermediate (I), or No Data (ND)"]
            if res_class in ["ND","","-"]:
                continue
            else:
                entry[3]="UNKNOWN"
                entry[8]=res_class
        elif tag in ["biek2012","casali2014","comas2013","guerra2014","manson2017","perdigao2014","phelan2016","tbarc.india","tbarc.cdrc","tbarc.iran","tbarc.mrc","winglee2015"]:
            continue
        else:
            print("[INFO] Unknown tag: {}".format(tag))
        # I need to get the biosample
        if tag=="bryant2013":
            entry[0]=i["Accession Number"]
        else:
            r = trace_request(db=getNCBIIdType(i["Accession Number"]), term=i["Accession Number"])
            reader=csv.DictReader(r.iter_lines(decode_unicode=True))
            biosamples=[]
            for ent in reader:
                try:
                    #print(ent["BioSample"])
                    biosamples.append(ent["BioSample"])
                    if(len(set(biosamples))>1):
                        print("[ERROR] More than one biosample for {} [tag: {}]. Skipping this NCBI ID".format(i["Accession Number"],tag))
                        #print(r.text)
                        continue
                except:
                    print("[ERROR] I  cannot find a biosample for {}".format(i["Accession Number"]))
                    #print(r.text)
                    continue
            entry[0] = next(iter(set(biosamples)))
        data_to_write.append(entry)
    #print("[INFO] {}/{} entries skipped ({:.1f}%)".format(entries_skipped,num_entries,entries_skipped/num_entries*100))
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


dat=build_table("./sources/resistance_data/CuratedPhenotypesUpdated_2018-05-14.csv",details)
discarded=parse_rows_take_decisions(dat,"./curated_phenotypes/curated_phenotypes.res","./sources/critical_concentrations/criticalConcentrations.csv")
for entry in discarded:
    print(entry)

