#!/usr/bin/env python
import csv
from collections import *

def recursive_defaultdict():
    return defaultdict(recursive_defaultdict)

def get_collisions(files: list) -> dict:
    """Takes a list of res files as input and returns a dictionary"""
    #the format of a res file is: <biosample> <antb> <res_class> <tags>
    d_res_class=recursive_defaultdict()
    d_tags=recursive_defaultdict()
    collisions={}
    for current_file in files:
        inp = csv.reader(open(current_file,"r"), delimiter='\t')
        for i,entry in enumerate(inp):
            if len(entry) != 4:
                print("[WARNING] Invalid format. I will skip line {}.".format(i))
                continue
            biosample=entry[0]
            antb=entry[1]
            res_class=entry[2]
            tags=entry[3]
            #I get the current value (if it exists!)
            if d_res_class[biosample][antb]:
                res_class_from_dict=d_res_class[biosample][antb]
                if(res_class_from_dict != res_class):
                    d_tags[biosample][antb]["tags"]+=","+tags
                    collisions[biosample+"@"+antb]=d_tags[biosample][antb]["tags"]
            else:
                d_res_class[biosample][antb]=res_class
                d_tags[biosample][antb]["tags"]=tags
    print("[INFO] I found {} collisions".format(len(collisions.keys())))
    return(collisions)

def generate_summary_from_res(res_files: list,out_file) -> None:
    """Takes a list of res files as input and builds a summary table for the resistance data (out_file). Does not check for collisions"""
    #the format of a res file is: <biosample> <antb> <res_class> <tags>
    d_res_class=recursive_defaultdict()
    d_antb={}
    for current_file in res_files:
        inp = csv.reader(open(current_file,"r"), delimiter='\t')
        for entry in inp:
            biosample=entry[0]
            antb=entry[1]
            d_antb[antb]=1
            res_class=entry[2]
            tags=entry[3]
            d_res_class[biosample][antb]=res_class
    with open(out_file, "w") as outf:
        outf.write("Isolate\t{}\n".format("\t".join(sorted(d_antb.keys()))))
        for biosample in sorted(d_res_class.keys()):
            entry=[biosample]
            for antb in sorted(d_antb.keys()):
                try:
                    res_class=d_res_class[biosample][antb]
                    if res_class in ["R","S"]:
                        entry.append(res_class)
                    else:
                        entry.append("")
                except:
                    entry.append("")
                    pass
            outf.write("\t".join(entry)+"\n")
    print("[INFO] I generated the report! ({})".format(out_file))
    return()

def generate_summary_from_res_remove_collisions(res_files: list,out_file: str) -> None:
    """Takes a list of res files as input and builds a summary table for the resistance data (out_file). It checks for collisions and removes them"""
    #the format of a res file is: <biosample> <antb> <res_class> <tags>
    d_res_class=recursive_defaultdict()
    d_antb={}
    for current_file in res_files:
        inp = csv.reader(open(current_file,"r"), delimiter='\t')
        for entry in inp:
            biosample=entry[0]
            antb=entry[1]
            d_antb[antb]=1
            res_class=entry[2]
            tags=entry[3]
            d_res_class[biosample][antb]=res_class
    collisions=get_collisions(res_files)
    for collision in collisions:
        explode_collision=collision.split("@")
        # I set the entry to ""
        d_res_class[explode_collision[0]][explode_collision[1]]=""
    with open(out_file, "w") as outf:
        outf.write("Isolate\t{}\n".format("\t".join(sorted(d_antb.keys()))))
        for biosample in sorted(d_res_class.keys()):
            entry=[biosample]
            for antb in sorted(d_antb.keys()):
                try:
                    res_class=d_res_class[biosample][antb]
                    if res_class in ["R","S"]:
                        entry.append(res_class)
                    else:
                        entry.append("")
                except:
                    entry.append("")
                    pass
            outf.write("\t".join(entry)+"\n")
    print("[INFO] I generated the report! ({})".format(out_file))
    return()
