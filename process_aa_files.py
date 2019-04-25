import sys
import os
import string
from os import listdir
import time
from DataExtractionSubmodules import RunInterProScan
process_tag = sys.argv[1]
#this script requires a folder containing amino acid files that contains the "process_tag"
#amino acid files are read from the folder and processed by InterProScan, after which they are deleted
#the results are printed to a file labeled "interproscan_" followed by a time stamp
def run_interproscan(infile, outputfilename):
	InterProScanExecutablePath = "/home/u7/pgnelson/my_extra/my_interproscan/interproscan-5.33-72.0/interproscan.sh"
	InterProScanOutputPath = "/home/u7/pgnelson/my_extra/out.txt"
	interpro_out = RunInterProScan(infile, InterProScanExecutablePath, InterProScanOutputPath)
	outputfile = open(outputfilename , "a")
	for key in interpro_out:
		num_entries = len(interpro_out[key]['PfamUID'])
		for entry in range(0, num_entries):
			print(key, interpro_out[key]['PfamUID'][entry], interpro_out[key]['PfamStart'][entry], interpro_out[key]['PfamStop'][entry], interpro_out[key]['Eval'][entry], sep = ",")
			print(key, interpro_out[key]['PfamUID'][entry], interpro_out[key]['PfamStart'][entry], interpro_out[key]['PfamStop'][entry], interpro_out[key]['Eval'][entry], sep = ",", file = outputfile)
	outputfile.close()

os.chdir("/home/u7/pgnelson/my_extra")
directories = os.listdir('.')
for directory in directories:
	if directory.find('aa_files') > -1 and directory.find(process_tag) > -1 and os.path.isdir(directory):
		os.chdir(directory)
		break
files  = os.listdir()
outputfilename = "interproscan_"+str(time.time())+".txt"
num_files_processed = 0
for filename in files:
	if filename.find(".fasta")>0 and filename.find("done")==-1:
		run_interproscan(filename, outputfilename)
		os.remove(filename)
		num_files_processed +=1
