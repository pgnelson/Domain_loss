import sys, json, os, csv, time, shutil, ftputil, string, re, datetime
import numpy as np
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
#This script was written by Sara Willis



def RunInterProScan(FastaFilename, InterProScanExecutablePath, InterProScanOutputPath):
	# InterProScan is run only using Pfam as the reference database using the input specifications.
	# Output is redirected to a tab-delimited file.
	inter_out = open(InterProScanOutputPath, 'w')
	inter_out.close()
	os.system('%s --applications Pfam --input %s --outfile %s -f TSV' %(InterProScanExecutablePath,FastaFilename,InterProScanOutputPath))

	# The results dictionary is defined
	InterProScanResultsDictionary = {}
	# The InterProScan output file is opened for parsing
	try:
		with open(InterProScanOutputPath, 'r') as f:
			reader = csv.reader(f, delimiter = '\t')
			# Each row is parsed for the relevant data
			for row in reader:
				# Accession of the sequence that was searched against the Pfam database
				ProteinAccession = row[0]
				# Information about the protein. Will not be used
				SequenceMD5Digest = row[1]
				# Length of the searched sequence
				SequenceLength = row[2]
				# More descriptions that will not be used
				Analysis = row[3]
				# Pfam hit ID
				SignatureAccession = row[4]
				# Description of Pfam
				SignatureDescription = row[5]
				# Location of Pfam domain start in protein coordinates
				StartLocation = row[6]
				# Location of Pfam domain stop in protein coordinates
				StopLocation = row[7]
				# Evalue associated with Pfam hit
				Score = row[8]
				# More descriptions that won't be used
				Status = row[9]
				# Date -- will not be used
				Date = row[10]

				'''
				Each entry in the InterProScanResultsDictionary will have nested dictionary with the following format:

				{Protein Accession: {PfamUID: [Pfam_1, ..., Pfam_n], PfamStart: [Pfam_1_start, ..., Pfam_n_start], PfamStop: [Pfam_1_stop, ..., Pfam_n_stop], Eval: [Pfam_1_eval, ..., Pfam_n_eval]}}
				'''
				# If the searched protein is not already in the dictionary, it is added along with all
				# the relevant information
				if ProteinAccession not in InterProScanResultsDictionary:
					# Each entry has its relevant information added as a length-one list so that further Pfam
					# data can be appended to the lists
					InterProScanResultsDictionary[ProteinAccession] = {'PfamUID':[SignatureAccession], 'PfamStart':[StartLocation], 'PfamStop':[StopLocation], 'Eval': [Score]}
				else:
					# If the searched protein is already in the dictionary, then the relevant information is
					# appended to the preexisting lists
					InterProScanResultsDictionary[ProteinAccession]['PfamUID'].append(SignatureAccession)
					InterProScanResultsDictionary[ProteinAccession]['PfamStart'].append(StartLocation)
					InterProScanResultsDictionary[ProteinAccession]['PfamStop'].append(StopLocation)
					InterProScanResultsDictionary[ProteinAccession]['Eval'].append(Score)
	except:
		InterProScanResultsDictionary = {}
	# Once the file is completely parsed, it is removed
	if os.path.exists(InterProScanOutputPath) == True:
		os.remove(InterProScanOutputPath)
	# The user is then returned the dictionary
	return InterProScanResultsDictionary
