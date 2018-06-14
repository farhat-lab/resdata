#!/usr/bin/env python


def write_table(rows,file_out,fields):
    with open(file_out,"w") as outf:
        outf.write("\t".join(fields)+"\n")
        for row in rows: 
            outf.write("\t".join(row)+"\n")
