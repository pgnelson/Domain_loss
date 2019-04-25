from ftplib import FTP
import gzip
import sys
import os
import string
#this script downloads the newest annotated genome of a given species from genbank and returns the file name
def download_species(species_name):
	num_chromosomes = 0
	assembled_chromosomes = 0
	outfilename = ""
	with FTP("ftp.ncbi.nih.gov") as ftp:
		ftp.login()
		ftp.cwd("/genomes/refseq")
		curr_dir = ftp.pwd()
		kingdoms = ["plant", "fungi", "invertebrate", "protozoa", "vertebrate_mammalian", "vertebrate_other"]
		for kingdom in kingdoms:
			try:
				ftp.cwd(kingdom)
				ftp.cwd(species_name + "/latest_assembly_versions")
				break
			except:
				if kingdom == "vertebrate_other":
					print("species not found")
					return("")
				ftp.cwd("/genomes/refseq")
		files = []
		ftp.dir(files.append)
		print(kingdom, len(files))
		latest_version_num = 0
		latest_version = ""
		for file in files:
			filename = file.split()[8]
			version_num=filename.split("_")
			version_num  = float(version_num[1].lstrip("0"))
			if version_num > latest_version_num:
				print(version_num)
				latest_version_num = version_num
				latest_version = filename
		ftp.cwd(latest_version)
		filenames = ftp.nlst()
		genbank_filename = ""
		for filename in filenames:
			print(filename)
			if filename.find("genomic.gbff.gz")>0 and filename.find("rna")==-1:
				print(filename)
				genbank_filename = filename
				file = open(filename, 'wb')
				ftp.retrbinary('RETR '+ filename, file.write)
				file.close()
		ftp.quit()
		f = gzip.open(genbank_filename, 'rt')
		file_content = f.read()
		outfilename = species_name+".gb"
		out = open(outfilename, "w")
		out.write(file_content)
		f.close()
		os.remove(genbank_filename)
	return(outfilename)
