from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from DataExtractionSubmodules import RunInterProScan
import os
import time
from get_genomes_from_genbank import download_species
from Bio import SeqIO
import re
import sys
import random
import psutil
process_id = str(sys.argv[1])
#if a file exists named "Log_/process_id", the script will read the file and start with the species, contig, and amino acid number indicated by the log file
#otherwise it will start with the first species in the "species_list" below and create a log file.

#write_aas takes a list of amino acid sequences and writes them to a file such that each file contains 20000 sequences.
#When a folder contains 25 files it creats a new folder. It also updates the log file with the position of the last amino acid written.
def write_aas(aas, species, num_dir_generated, num_aa_files_generated, log_dir, num_seq_per_analysis, aa_processed_from_contig, index, process_id, aa_out_dir):
	num_aa_files_printed_writing = 0
	while(len(aas) >= num_seq_per_analysis):
		#working_dir = '/Users/GavinLeighton/downloads/generating_working/'
		working_dir = '/home/u7/pgnelson/my_extra/generating_working/'
	#	if psutil.disk_usage('/home/u7/pgnelson/my_extra')[3] > 90:
	#		time.sleep(60)
		if os.path.isdir(working_dir) == False:
			os.mkdir(working_dir)
		aa_filename = working_dir + "aa_file_" + species + "_"+ str(num_dir_generated) + "_" + str(random.randint(0,10**6)) + ".fasta"
		try:
			aafile = open(aa_filename, "a")
		except:
			time.sleep(1)
			aafile = open(aa_filename, "a")
		for i in range(num_seq_per_analysis):
			print(aas[i][0], aas[i][1], sep = "\n", file = aafile)
		aafile.close()
		aas = list(aas[num_seq_per_analysis:])
		num_aa_files_printed_writing += num_seq_per_analysis
		aa_processed_from_contig += num_seq_per_analysis
		log = open(log_dir + "/log_" + process_id + ".txt", "w")
		print(species, index, aa_processed_from_contig, sep = ",", file = log)
		log.close()
		num_aa_files_generated +=1
		if num_aa_files_generated >= 25:
			curr_output_dir = aa_out_dir + "/aa_files_" + species +"_"+ str(num_dir_generated) + "/"
			if os.path.isdir(curr_output_dir) == False:
				os.rename(working_dir, curr_output_dir)
			num_dir_generated +=1
			num_aa_files_generated = 0
	return(num_dir_generated, num_aa_files_generated, num_aa_files_printed_writing, aa_processed_from_contig)

#process_contig takes a contig from a genbank genome file, removes annoted coding sequences, uses aa_w_Ns_cut_out to translate the remainder, and passes the translated peptides to write_aas to be written to file
def process_contig(contig, species, contig_aa_list, aa_files_todo_before_printing, writing_aas, num_dir_generated, num_aa_files_generated, num_seq_per_analysis, num_aacids, index, process_id, aa_out_dir):
	num_aa_files_printed_contig = 0
	num_aacids_contig = 0
	cd_locations = []
	at = 0
	cg = 0
	at_count_coding	= 0
	cg_count_coding = 0
	seq_id = species + "_" + contig.id
	num_cds = 0
	aa_processed_from_contig = 0
	writing_aas =1
	for feature in contig.features:
		if feature.type == "CDS":
			cdlocation = [feature.location.nofuzzy_start,feature.location.nofuzzy_end]
			cd_locations.append(cdlocation)
			num_cds +=1
	cd_locations = sorted(cd_locations,key=lambda x:x[0])
	start = 0
	if len(cd_locations) > 0: #removes cds, translates, and writes aas
		for location in cd_locations:
			end = location[0]
			at_count_coding += contig.seq[location[0]:location[1]].count("A")+contig.seq[location[0]:location[1]].count("T")+contig.seq[location[0]:location[1]].count("a")+contig.seq[location[0]:location[1]].count("t")
			cg_count_coding += contig.seq[location[0]:location[1]].count("C")+contig.seq[location[0]:location[1]].count("G")+contig.seq[location[0]:location[1]].count("c")+contig.seq[location[0]:location[1]].count("g")
			if len(contig.seq[start:end])>0:
				fragment = contig.seq[start:end]
				aa_w_Ns_cut_out_return = aa_w_Ns_cut_out(fragment, seq_id)
				#returns(new_aas, at, cg, num_new_aacids)
				contig_aa_list = contig_aa_list + aa_w_Ns_cut_out_return[0]
				at += aa_w_Ns_cut_out_return[1]
				cg += aa_w_Ns_cut_out_return[2]
				num_aacids_contig += aa_w_Ns_cut_out_return[3]
			start = location[1]
			if writing_aas == 0:
				if len(contig_aa_list) > aa_files_todo_before_printing:
					contig_aa_list = list(contig_aa_list[aa_files_todo_before_printing:])
					aa_files_todo_before_printing = 0
					writing_aas == 1
				else:
					aa_files_todo_before_printing -=len(contig_aa_list)
					contig_aa_list = []
			if len(contig_aa_list) > 20000 and writing_aas == 1:
				print_pos = int(len(contig_aa_list)/20000)*20000
				aa_to_write = contig_aa_list[:print_pos]
				contig_aa_list = contig_aa_list[print_pos:]
				write_output = write_aas(aa_to_write, species, num_dir_generated, num_aa_files_generated, log_dir, num_seq_per_analysis, aa_processed_from_contig, index, process_id, aa_out_dir)
				#returns(num_dir_generated, num_aa_files_generated, num_aa_files_printed_writing, aa_processed_from_contig)
				num_dir_generated = write_output[0]
				num_aa_files_generated = write_output[1]
				num_aa_files_printed_contig += write_output[2]
				aa_processed_from_contig = write_output[3]
	elif len(contig.seq)>0 :#in case there are no cds
		aa_w_Ns_cut_out_return =  aa_w_Ns_cut_out(contig.seq, seq_id)
		contig_aa_list = contig_aa_list + aa_w_Ns_cut_out_return[0]
		at += aa_w_Ns_cut_out_return[1]
		cg += aa_w_Ns_cut_out_return[2]
		num_aacids_contig += aa_w_Ns_cut_out_return[3]
		if len(contig_aa_list) > 20000 and writing_aas == 1:
			print_pos = int(len(contig_aa_list)/20000)*20000
			aa_to_write = contig_aa_list[:print_pos]
			contig_aa_list = contig_aa_list[print_pos:]
			write_output = write_aas(aa_to_write, species, num_dir_generated, num_aa_files_generated, log_dir, num_seq_per_analysis, aa_processed_from_contig, index, process_id, aa_out_dir)
			num_dir_generated = write_output[0]
			num_aa_files_generated = write_output[1]
			num_aa_files_printed_contig += write_output[2]
			aa_processed_from_contig = write_output[3]

	return(contig_aa_list, at, cg, at_count_coding, cg_count_coding, num_dir_generated, num_aa_files_generated, num_aa_files_printed_contig , aa_processed_from_contig, num_aacids_contig)

#aa_w_Ns_cut_out takes a nucleic acid sequence, cuts out any runs of "N">9, and translates the remainder using get_AA_runs and returns the resulting peptide list
def aa_w_Ns_cut_out(nsequence, seq_id):
	at = 0
	cg = 0
	frag_start = 0
	num_new_aacids = 0
	new_aas = []
	while str(nsequence[frag_start:]).find("NNNNNNNNN")>0:
		n_start_loc = str(nsequence[frag_start:]).find("NNNNNNNNN") + frag_start
		frag_aas =  get_AA_runs(nsequence[frag_start:n_start_loc-1], seq_id, frag_start)
		num_new_aacids += len(frag_aas)
		new_aas =  new_aas + frag_aas
		at = nsequence[frag_start:n_start_loc-1].count("A") + nsequence[frag_start:n_start_loc-1].count("T")+nsequence[frag_start:n_start_loc-1].count("a") + nsequence[frag_start:n_start_loc-1].count("t")
		cg = nsequence[frag_start:n_start_loc-1].count("C") + nsequence[frag_start:n_start_loc-1].count("G")+nsequence[frag_start:n_start_loc-1].count("c") + nsequence[frag_start:n_start_loc-1].count("g")
		try:
			frag_start = re.search(r'[^N]', str(nsequence[n_start_loc:])).start() + n_start_loc
		except:
			frag_start = len(nsequence) -1
	if len(nsequence[frag_start:])>0 :
		at = nsequence[frag_start:].count("A") + nsequence[frag_start:].count("T")+nsequence[frag_start:].count("a") + nsequence[frag_start:].count("t")
		cg = nsequence[frag_start:].count("C") + nsequence[frag_start:].count("G")+nsequence[frag_start:].count("c") + nsequence[frag_start:].count("g")
		frag_aas = get_AA_runs(nsequence[frag_start:], seq_id, frag_start)
		num_new_aacids += len(frag_aas)
		new_aas =  new_aas + frag_aas
	return(new_aas, at, cg, num_new_aacids)
#get_AA_runs translates an nucleic acid sequence and makes a list of the resulting peptides (each peptide is a run of amino acids with no stop codons)
def get_AA_runs(sequence, id, pos):
	processed_list = []
	for frame in range(0,3):
		aa_sequence = sequence[frame:].translate()
		aa_list = aa_sequence.split("*")
		position = pos + frame
		for peptide in aa_list:
			label = ">" +id+ "_pos_" + str(position) + "_frame_" + str(frame) + "_forward"
			if len(peptide) >= 10:
				processed_list.append([label, str(peptide)])
			position += len(peptide)*3+3
	rev = sequence.reverse_complement()
	for frame in range(0,3):
		position = pos + len(rev)
		aa_sequence = rev[frame:].translate()
		aa_list = aa_sequence.split("*")
		for peptide in aa_list:
			position -= len(peptide)*3+3
			label = ">"+id+ "_pos_" + str(position) + "_frame_" + str(frame) + "_rev"
			if len(peptide) >= 20:
				processed_list.append([label, str(peptide)])
	return(processed_list)
#process_species will download a genome from genbank, and generate files containing amino acid seqeunces from the genome
def process_species(species_name, curr_contig, curr_aas, process_id, aa_out_dir):
	num_aa_files_printed = 0
	num_aacids_printed = 0
	start_time = time.time()
	count_curr_species_aa_files = 0
	directories = os.listdir('.')
	num_dir_generated =0
	num_aa_files_generated = 0
	log_dir = "/home/u7/pgnelson"#os.getcwd()
	genbank_data  = download_species(species_name)
	num_aacids = 0
	at_count = 0
	cg_count = 0
	at_count_coding = 0
	cg_count_coding = 0
	aa_list = []
	num_seq_per_analysis = 20000

	if len(genbank_data)>0:
		for index, record in enumerate(SeqIO.parse(genbank_data, "genbank")):
			print(index, num_aacids, curr_contig)
			if index >= curr_contig:
				writing_aas = 1
				if index == curr_contig:
					writing_aas = 0
				else:
					curr_aas = 0
			#	if psutil.disk_usage('/home/u7/pgnelson/my_extra')[3] > 80:
			#		time.sleep(5*60)
				contig_results = process_contig(record, species, aa_list, curr_aas, writing_aas, num_dir_generated, num_aa_files_generated, num_seq_per_analysis, num_aacids, index, process_id, aa_out_dir)
				#returns(contig_aa_list, at, cg, at_count_coding, cg_count_coding, num_dir_generated, num_aa_files_generated, num_aa_files_printed , aa_processed_from_contig, num_aacids)

				aa_list = contig_results[0]
				at_count += contig_results[1]
				cg_count += contig_results[2]
				at_count_coding += contig_results[3]
				cg_count_coding += contig_results[4]
				num_dir_generated = contig_results[5]
				num_aa_files_generated = contig_results[6]
				num_aacids_printed += contig_results[7]
				aa_processed_from_contig = contig_results[8]
				num_aacids += contig_results[9]

		if len(aa_list) > 0:
			write_output = write_aas(aa_list, species, num_dir_generated, num_aa_files_generated, log_dir, num_seq_per_analysis, aa_processed_from_contig, index, process_id, aa_out_dir)
			num_dir_generated = write_output[0]
			num_aa_files_generated = write_output[1]
			num_aa_files_printed += write_output[2]
			aa_processed_from_contig = write_output[3]
		num_aacids_file = open("num_aacids.txt", "a")
		curr_output_dir = aa_out_dir + "/aa_files_" + species +"_"+ str(num_dir_generated) + "/"
		#working_dir = '/home/u7/pgnelson/my_extra/generating_working/'
		working_dir = '/home/u7/pgnelson/my_extra/generating_working/'

		if os.path.isdir(working_dir):
			if os.path.isdir(curr_output_dir):
				num_dir_generated +=1
				curr_output_dir = aa_out_dir + "/aa_files_" + species +"_"+ str(num_dir_generated) + "/"
			os.rename(working_dir, curr_output_dir)
		num_dir_generated +=1
		num_aa_files_generated = 0
		print (species_name, num_aacids_printed, num_aacids, at_count, cg_count, at_count_coding, cg_count_coding, sep = "\t", file = num_aacids_file)
		print (species_name, num_aacids_printed, num_aacids, at_count, cg_count, at_count_coding, cg_count_coding, sep = "\t")
		num_aacids_file.close()#

		os.remove(genbank_data)
		#os.remove(log_dir + "/log_" + process_id + ".txt")
	else:
		num_aacids_file = open("num_aacids.txt", "a")
		print (species_name, "not found", file = num_aacids_file)
		num_aacids_file.close()#
	return()


species_list=["Saccharomyces_cerevisiae","Arabidopsis_thaliana","Ailuropoda_melanoleuca","Anas_platyrhynchos","Anolis_carolinensis","Aotus_nancymaae","Bos_taurus","Caenorhabditis_elegans","Callithrix_jacchus","Canis_lupus","Capra_hircus","Carlito_syrichta","Cavia_aperea","Cavia_porcellus","Cebus_capucinus","Cercocebus_atys","Chinchilla_lanigera","Chlorocebus_sabaeus","Ciona_intestinalis","Colobus_angolensis","Cricetulus_griseus","Cricetulus_griseus","Danio_rerio","Dasypus_novemcinctus","Dipodomys_ordii","Drosophila_melanogaster","Eptatretus_burgeri","Equus_caballus","Felis_catus","Ficedula_albicollis","Fukomys_damarensis","Gadus_morhua","Gallus_gallus","Gasterosteus_aculeatus","Gorilla_gorilla","Heterocephalus_glaber","Homo_sapiens","Ictidomys_tridecemlineatus","Jaculus_jaculus","Latimeria_chalumnae","Lepisosteus_oculatus","Loxodonta_africana","Macaca_fascicularis","Macaca_mulatta","Macaca_nemestrina","Mandrillus_leucophaeus","Meleagris_gallopavo","Mesocricetus_auratus","Microcebus_murinus","Microtus_ochrogaster","Monodelphis_domestica","Mus_caroli","Mus_musculus","Mus_pahari","Mus_spretus","Mustela_putorius","Myotis_lucifugus","Nannospalax_galili","Nomascus_leucogenys","Octodon_degus","Ornithorhynchus_anatinus","Oryctolagus_cuniculus","Oryzias_latipes","Otolemur_garnettii","Ovis_aries","Pan_paniscus","Pan_troglodytes","Panthera_pardus","Panthera_tigris","Papio_anubis","Pelodiscus_sinensis","Peromyscus_maniculatus","Petromyzon_marinus","Propithecus_coquereli","Rattus_norvegicus","Rhinopithecus_bieti","Rhinopithecus_roxellana","Saimiri_boliviensis","Sarcophilus_harrisii","Sus_scrofa","Taeniopygia_guttata","Takifugu_rubripes","Tetraodon_nigroviridis","Xenopus_tropicalis","Xiphophorus_maculatus","Aegilops_tauschii","Amborella_trichopoda","Arabidopsis_lyrata","Beta_vulgaris","Brachypodium_distachyon","Brassica_napus","Brassica_oleracea","Brassica_rapa","Cajanus_cajan","Chlamydomonas_reinhardtii","Chondrus_crispus","Cucumis_sativus","Cyanidioschyzon_merolae","Daucus_carota","Galdieria_sulphuraria","Glycine_max","Gossypium_raimondii","Helianthus_annuus","Hordeum_vulgare","Lupinus_angustifolius","Manihot_esculenta","Medicago_truncatula","Musa_acuminata","Oryza_barthii","Oryza_brachyantha","Oryza_glaberrima","Oryza_meridionalis","Oryza_punctata","Oryza_rufipogon","Ostreococcus_lucimarinus","Phaseolus_vulgaris","Physcomitrella_patens","Populus_trichocarpa","Prunus_persica","Selaginella_moellendorffii","Solanum_tuberosum","Sorghum_bicolor","Theobroma_cacao","Trifolium_pratense","Triticum_aestivum","Triticum_urartu","Vitis_vinifera","Zea_mays","Acyrthosiphon_pisum","Aedes_aegypti","Amphimedon_queenslandica","Anopheles_darlingi","Anopheles_gambiae","Apis_mellifera","Atta_cephalotes","Belgica_antarctica","Bombyx_mori","Brugia_malayi","Caenorhabditis_brenneri","Caenorhabditis_briggsae","Caenorhabditis_japonica","Caenorhabditis_remanei","Capitella_teleta","Crassostrea_gigas","Culex_quinquefasciatus","Danaus_plexippus","Daphnia_pulex","Dendroctonus_ponderosae","Drosophila_ananassae","Drosophila_erecta","Drosophila_grimshawi","Drosophila_mojavensis","Drosophila_persimilis","Drosophila_pseudoobscura","Drosophila_sechellia","Drosophila_simulans","Drosophila_virilis","Drosophila_willistoni","Drosophila_yakuba","Heliconius_melpomene","Helobdella_robusta","Lepeophtheirus_salmonis","Lucilia_cuprina","Mayetiola_destructor","Melitaea_cinxia","Mnemiopsis_leidyi","Nematostella_vectensis","Octopus_bimaculoides","Pediculus_humanus","Pristionchus_pacificus","Rhodnius_prolixus","Solenopsis_invicta","Stegodyphus_mimosarum","Strigamia_maritima","Tribolium_castaneum","Trichoplax_adhaerens","Zootermopsis_nevadensis","Bigelowiella_natans","Dictyostelium_discoideum","Emiliania_huxleyi","Entamoeba_histolytica","Giardia_intestinalis","Hyaloperonospora_arabidopsidis","Leishmania_major","Phaeodactylum_tricornutum","Phytophthora_infestans","Phytophthora_parasitica","Phytophthora_ramorum","Phytophthora_sojae","Plasmodium_berghei","Plasmodium_chabaudi","Plasmodium_falciparum","Plasmodium_knowlesi","Plasmodium_vivax","Thalassiosira_pseudonana","Toxoplasma_gondii","Trypanosoma_brucei","Aspergillus_clavatus","Aspergillus_oryzae","Aspergillus_terreus","Fusarium_fujikuroi","Schizosaccharomyces_pombe","Ustilago_maydis","Volvox_carteri","Acropora_digitifera","Aedes_albopictus","Apis_cerana","Apis_dorsata","Apis_florea","Aplysia_californica","Bactrocera_dorsalis","Bactrocera_oleae","Bemisia_tabaci","Bicyclus_anynana","Branchiostoma_belcheri","Branchiostoma_floridae","Camponotus_floridanus","Ceratina_calcarata","Ceratitis_capitata","Ceratosolen_solmsi","Cimex_lectularius","Crassostrea_virginica","Cyphomyrmex_costatus","Drosophila_arizonae","Drosophila_biarmipes","Drosophila_bipectinata","Drosophila_busckii","Drosophila_elegans","Drosophila_eugracilis","Drosophila_ficusphila","Drosophila_kikkawai","Drosophila_miranda","Drosophila_navojoa","Drosophila_obscura","Drosophila_rhopaloa","Drosophila_serrata","Drosophila_suzukii","Drosophila_takahashii","Dufourea_novaeangliae","Echinococcus_granulosus","Eufriesea_mexicana","Eurytemora_affinis","Metaseiulus_occidentalis","Habropoda_laboriosa","Hydra_vulgaris","Limulus_polyphemus","Linepithema_humile","Mizuhopecten_yessoensis","Monomorium_pharaonis","Musca_domestica","Myzus_persicae","Nicrophorus_vespilloides","Nilaparvata_lugens","Onthophagus_taurus","Montastraea_faveolata","Orussus_abietinus","Papilio_machaon","Papilio_polytes","Papilio_xuthus","Parasteatoda_tepidariorum","Pieris_rapae","Pogonomyrmex_barbatus","Priapulus_caudatus","Pseudomyrmex_gracilis","Saccoglossus_kowalevskii","Stomoxys_calcitrans","Trachymyrmex_cornetzi","Trachymyrmex_septentrionalis","Trachymyrmex_zeteki","Varroa_destructor","Vollenhovia_emeryi","Wasmannia_auropunctata","Bactrocera_cucurbitae","Ananas_comosus","Camelina_sativa","Capsella_rubella","Capsicum_annuum","Carica_papaya","Cicer_arietinum","Citrus_clementina","Citrus_sinensis","Cucumis_melo","Elaeis_guineensis","Erythranthe_guttata","Eucalyptus_grandis","Eutrema_salsugineum","Fragaria_vesca","Gossypium_arboreum","Gossypium_hirsutum","Hevea_brasiliensis","Ipomoea_nil","Jatropha_curcas","Juglans_regia","Lactuca_sativa","Malus_domestica","Micromonas_pusilla","Momordica_charantia","Morus_notabilis","Nelumbo_nucifera","Nicotiana_sylvestris","Nicotiana_tabacum","Nicotiana_tomentosiformis","Oryza_sativa","Ostreococcus_tauri","Phoenix_dactylifera","Populus_euphratica","Prunus_avium","Prunus_mume","Pyrus_x_bretschneideri","Raphanus_sativus","Ricinus_communis","Sesamum_indicum","Solanum_pennellii","Spinacia_oleracea","Tarenaya_hassleriana","Vigna_radiata","Acinonyx_jubatus","Balaenoptera_acutorostrata","Bison_bison","Bos_indicus","Bos_mutus","Bubalus_bubalis","Camelus_bactrianus","Camelus_dromedarius","Camelus_ferus","Castor_canadensis","Ceratotherium_simum","Chrysochloris_asiatica","Condylura_cristata","Echinops_telfairi","Elephantulus_edwardii","Eptesicus_fuscus","Equus_asinus","Equus_przewalskii","Erinaceus_europaeus","Hipposideros_armiger","Leptonychotes_weddellii","Lipotes_vexillifer","Manis_javanica","Marmota_marmota","Meriones_unguiculatus","Miniopterus_natalensis","Myotis_brandtii","Myotis_davidii","Monachus_monachus","Ochotona_princeps","Odobenus_rosmarus","Odocoileus_virginianus","Orcinus_orca","Orycteropus_afer","Pantholops_hodgsonii","Phascolarctos_cinereus","Physeter_catodon","Pongo_abelii","Pteropus_alecto","Pteropus_vampyrus","Rhinolophus_sinicus","Rousettus_aegyptiacus","Sorex_araneus","Trichechus_manatus","Tursiops_truncatus","Ursus_maritimus","Vicugna_pacos","Acanthisitta_chloris","Acanthochromis_polyacanthus","Alligator_mississippiensis","Alligator_sinensis","Anser_cygnoides","Caprimulgus_carolinensis","Apaloderma_vittatum","Aptenodytes_forsteri","Aquila_chrysaetos","Balearica_regulorum","Boleophthalmus_pectinirostris","Buceros_rhinoceros","Calidris_pugnax","Callorhinchus_milii","Calypte_anna","Cariama_cristata","Chaetura_pelagica","Charadrius_vociferus","Chelonia_mydas","Chlamydotis_macqueenii","Chrysemys_picta","Clupea_harengus","Colius_striatus","Columba_livia","Corvus_brachyrhynchos","Coturnix_japonica","Crocodylus_porosus","Cuculus_canorus","Cyprinodon_variegatus","Cyprinus_carpio","Egretta_garzetta","Esox_lucius","Eurypyga_helias","Falco_cherrug","Falco_peregrinus","Fulmarus_glacialis","Fundulus_heteroclitus","Gavia_stellata","Gavialis_gangeticus","Gekko_japonicus","Geospiza_fortis","Haliaeetus_albicilla","Haliaeetus_leucocephalus","Haplochromis_burtoni","Hippocampus_comes","Ictalurus_punctatus","Labrus_bergylta","Larimichthys_crocea","Lates_calcarifer","Lepidothrix_coronata","Leptosomus_discolor","Lonchura_striata","Manacus_vitellinus","Maylandia_zebra","Melopsittacus_undulatus","Merops_nubicus","Mesitornis_unicolor","Nanorana_parkeri","Neolamprologus_brichardi","Nestor_notabilis","Nipponia_nippon","Nothobranchius_furzeri","Notothenia_coriiceps","Numida_meleagris","Oncorhynchus_mykiss","Opisthocomus_hoazin","Parus_major","Pelecanus_crispus","Phaethon_lepturus","Phalacrocorax_carbo","Picoides_pubescens","Poecilia_mexicana","Pogona_vitticeps","Protobothrops_mucrosquamatus","Pseudopodoces_humilis","Pterocles_gutturalis","Pundamilia_nyererei","Pygocentrus_nattereri","Pygoscelis_adeliae","Salmo_salar","Scleropages_formosus","Serinus_canaria","Seriola_dumerili","Sinocyclocheilus_anshuiensis","Sinocyclocheilus_grahami","Sinocyclocheilus_rhinocerous","Stegastes_partitus","Struthio_camelus","Sturnus_vulgaris","Tauraco_erythrolophus","Thamnophis_sirtalis","Tyto_alba","Xenopus_laevis","Zonotrichia_albicollis","Ursus_arctos","Dendrobium_catenatum"]
# skipped Limulus_polyphemus
log_dir = os.getcwd()#"/home/u7/pgnelson"#os.getcwd()
#aa_out_dir = "/Users/GavinLeighton/downloads"
aa_out_dir = "/home/u7/pgnelson/my_extra"
curr_species = ""
curr_contig = 0
curr_aas = 0
try:
	log = open(log_dir + "/log_" + process_id + ".txt", "r")
	curr_log = str(log.read())
	log_items = curr_log.split(",")
	curr_species = log_items[0]
	curr_contig = int(log_items[1])
	curr_aas = int(log_items[2])
	print(curr_species, curr_contig, curr_aas)
	log.close()
except:
	curr_species = species_list[0]
	curr_contig = 0
	curr_aas = 0

#curr_species = ""
#curr_contig =0
#curr_aas = 0
for species in species_list:
	if species == curr_species or curr_species == "":
		curr_species = ""
		print(species)
		log = open(log_dir + "/log_" + process_id + ".txt", "w")
		print(species, 0, 0, sep = ",", file = log)
		log.close()
		process_species(species, curr_contig, curr_aas, process_id, aa_out_dir)
