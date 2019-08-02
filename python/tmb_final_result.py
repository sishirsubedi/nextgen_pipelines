import sys
import os
import pandas as pd
import numpy as np

def get_seq_stats (f1,f2):
    cols = ["Total-Reads","Q20","Duplicate","Coverage","Coverage-2X","Coverage-10X","Coverage-100X"]
    in_file = open(f1,"r")
    out_file = open(f2,"a")
    for line in in_file:
        if line[0] != "#":
            vals = line.split(',')
            if vals[0] in cols:
                out_file.write("Tumor_"+vals[0]+","+vals[1]+'\n')
                out_file.write("Normal_"+vals[0]+","+vals[2]+'\n')
    in_file.close()
    out_file.close()


def get_breadth_coverage (f1,f2):
    in_file = open(f1,"r")
    out_file = open(f2,"a")
    for line in in_file:
        vals = line.split(' ')
        out_file.write("Tumor_Breadth-Coverage"+","+vals[0]+'\n')
    in_file.close()
    out_file.close()


def get_variant_results (f1,f2):
    cols = ["varscan-strelka","varscan-mutect","mutect-strelka","varscan-strelka-mutect"]
    in_file = open(f1,"r")
    out_file = open(f2,"a")
    for line in in_file:
        if line[0] != "#":
            vals = line.split(',')
            if vals[0] in cols:
                out_file.write(vals[0]+","+vals[1]+'\n')
    in_file.close()
    out_file.close()

def get_tmb_results (f1,f2):
    in_file = open(f1,"r")
    out_file = open(f2,"w")
    for line in in_file:
        vals = line.split(',')
        out_file.write("SampleID"+","+vals[0]+'\n')
        out_file.write("Sample"+","+vals[1]+'\n')
        out_file.write("TMB-Total-Variants"+","+vals[2]+'\n')
        out_file.write("TMB-Score"+","+vals[3]+'\n')
        out_file.write("TMB-Group"+","+vals[4]+'\n')
    in_file.close()
    out_file.close()


file_seq_stats = sys.argv[1]  # alignment result
file_breadth_coverage = sys.argv[2] # tumor breadth coverage
file_variant_results = sys.argv[3] # pair variant results
file_tmb_result = sys.argv[4] # tmb results
file_final_result = sys.argv[5]  # out file

print(file_seq_stats,
file_breadth_coverage,
file_variant_results,
file_tmb_result,
file_final_result
)

get_tmb_results(file_tmb_result,file_final_result)
get_variant_results(file_variant_results,file_final_result)
get_seq_stats(file_seq_stats,file_final_result)
get_breadth_coverage(file_breadth_coverage,file_final_result)
os.system("sed -i \'/^$/d\' %s" %file_final_result)
