#!/usr/bin/env python

"""
NetworKIN(tm), (C) 2005,2006,2007.
Drs Rune Linding & Lars Juhl Jensen

Released under stay the *#(@(#@)(%)(@!!!*$(# away license, until we published all our papers!
I.e. its NOT released, if you did not obtain this software from any of the above you will be
legally prosecuted.

Usage: ./networkin.py Organism FastaFile SitesFile

If no sites file is given NetworKIN will predict on all T/S/Y residues 
in the given sequences.
"""

### Coding strategy & TODO
# columns: 1-5, 10: substrate, position, id, networkin score, tree,...,substrate name
###

# Changelog
# 21.10.06: tested blast -a 8 option, no differences on proteome
# 21.10.06: Filtercode broken
#			-fixed
# 28.01.07: Testing autophosphorylation (self == 1 in update script)
# 30.07.07: Working on v1.5 milestone
# 07.05.08: Initiated 2.5 w. scaling factor
#

import sys, os, subprocess, re, tempfile, random, operator
import threading
from optparse import OptionParser
import multiprocessing
import csv
from string import *

#debugging
import time

# Weighting parameter, 0=only motif, 1=only STRING
# estimated while benchmarking an then hardcoded here
# feel free to play with it, but it is your own responsibility
ALPHAS = {"9606": 0.85, "4932": 0.65} # 4932 Saccharomyces cerevisiae

global options
# Temporary directory to store the files
#tempfile.tempdir= '/tmp'

# Location of NetworKIN input files (string network, alias file etc.)
#datadir = sys.argv[0].rsplit("/", 1)[0]+'/data'
#global DATADIR
#DATADIR = ""

# Number of threads used for NetPhorest and BLAST
# setting it to a high number can also help in case NetworKIN uses too much memory
# as NetPhorest result files will be read one after the other to save memory
# -> the more files, the less memory usage
#global NUMBER_OF_PROCESSES
#NUMBER_OF_PROCESSES = "1";

# should temporary files of netphorest be zipped?
# saves some diskspace and makes it sometimes faster (depending on the CPU/harddisk speed)
#global SAVE_DISKSPACE
#SAVE_DISKSPACE = "";

# Should the analysis be limited to specific trees?
# just comment it out if you want to predict on all
#limitTrees = ['SH2']
#limitTrees = ['KIN']

################################################################################
#                                                                              #
#                             the code starts here                             #
#                                                                              #
################################################################################

# Run system binary
'''
def myPopen(cmd):
	try:
		pipe = subprocess.Popen(cmd, shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout = pipe.stdout.readlines()
	except:
		sys.stderr.write('ERROR executing: '+repr(cmd)+'\n')
		sys.exit()
	else:
		print(stdout)
		return stdout 
'''	
### chat-gpt's myPopen:
def myPopen(cmd):
    try:
        pipe = subprocess.Popen(cmd, shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = pipe.communicate()
        # Decode the byte strings to utf-8
        stdout = stdout.decode('utf-8')
        stderr = stderr.decode('utf-8')
        if pipe.returncode != 0:
            # If the command failed, print the stderr and raise an exception
            error_message = 'ERROR executing: ' + repr(cmd) + '\n' + stderr
            sys.stderr.write(error_message + '\n')
            raise subprocess.CalledProcessError(pipe.returncode, cmd, output=error_message)
        else:
            # If the command succeeded, return the stdout
            return stdout
    except subprocess.CalledProcessError as e:
        # Handle the subprocess.CalledProcessError exception
        print("Error: Command '{}' returned non-zero exit status {}".format(e.cmd, e.returncode))
        return e.output  # Return the error message
    except Exception as e:
        # Handle other exceptions
        error_message = 'ERROR executing: ' + repr(cmd) + '\n' + str(e) + '\n'
        sys.stderr.write(error_message)
        return error_message

# Read sequences from fasta file
def readFasta(fastafile):
	id_seq = {}
	aminoacids = re.compile('[^ACDEFGHIKLMNPQRSTVWYXB]')
	data = fastafile.readlines()
	fastafile.close()
	for line in data:
		if line[0] != ';':
			if line[0] == '>':
				line = line.strip()
				id = line[1:].split(' ', 1)[0]
#				id = line[1:-1].split(' ', 1)[0]
				seq = ''
			else:
				seq += aminoacids.sub('', line)
				if len(seq) > 0:
					id_seq[id] = seq
		else:
			pass
	#print(id_seq)
	return id_seq

# Read phosphorylation sites from tsv file
# id -> position -> residue
def readPhosphoSites(sitefile):
	id_pos_res = {}
	if sitefile:
		data = sitefile.readlines()
		sitefile.close()
		for line in data:
			tokens = line.split('\t')
			id = tokens[0]
			pos = int(tokens[1])
			try:
				res = tokens[2][:-1]
			except:
				res = ""
			if id in id_pos_res:
				id_pos_res[id][pos] = res
			else:
				id_pos_res[id] = {pos: res}
	return id_pos_res

#Alias hashes
def readAliasFiles(organism, datadir):
	alias_hash = {}
	desc_hash = {}
	name_hash={}
	'''
	# Read alias db (alias_best.v9.0):
	try:
		alias_db = myPopen('gzip -cd %s/%s.alias_best.tsv.gz'%(datadir, organism))
		for line in alias_db:
			(taxID, seqID, alias) = line.strip().split('\t')[:3]
			alias_hash[seqID] = alias
	except:
		sys.stderr.write("No aliases available for organism: '%s'\n"%organism)
	'''
	# Read alias db (protein.aliases.v12.0):
	try:
		alias_db = myPopen('gzip -cd %s/%s.protein.aliases.v12.0.txt.gz' % (datadir, organism))

		alias_db = alias_db.split('\n')
		#print(alias_db)
		for line in alias_db:

			if line.startswith('#'):
				continue
			(taxID, seqID, alias) = line.split('\t')
			#print((taxID, seqID, alias))
			alias_hash[seqID] = alias
	except:
		sys.stderr.write("No aliases available for organism: '%s'\n"%organism)

	# Read desc db
	try:
		desc_db = myPopen('gzip -cd %s/%s.text_best.v9.0.tsv.gz'%(datadir, organism))
		#desc_db = myPopen('gzip -cd %s/%s.protein.info.v12.0.txt.gz'%(datadir, organism))
		desc_db = desc_db.split('\n')
		for line in desc_db:
			#print(line)
			(taxID, seqID, desc) = line.split('\t')[:3]
			desc_hash[seqID] = desc
	except:
		sys.stderr.write("No descriptions available for organism: '%s'\n"%organism)

	try:
		#desc_db = myPopen('gzip -cd %s/%s.text_best.v9.0.tsv.gz'%(datadir, organism))
		name_db = myPopen('gzip -cd %s/%s.protein.info.v12.0.txt.gz'%(datadir, organism))
		name_db = name_db.split('\n')
		for line in name_db:
			#print(line)
			(string_protein_id,preferred_name,protein_size,annotation) = line.split('\t')
			name_hash[string_protein_id] = preferred_name
	except:
		sys.stderr.write("No names available for organism: '%s'\n"%organism)

	#print(alias_hash)
	return alias_hash, desc_hash , name_hash
'''
# Run Netphorest
def runNetPhorest(id_seq, id_pos_res, save_diskspace, number_of_processes):
	id_pos_tree_pred = {}

	#check how many sequences we actually have
	number_of_sequences = 0
	if id_pos_res == {}:
		number_of_sequences = len(id_seq)
	else:
		number_of_sequences = len(id_pos_res)

	if number_of_sequences < number_of_processes:
		number_of_processes = number_of_sequences
	#print(number_of_sequences)
	# use multiple instances of netphorest
	file_in = list(range(number_of_processes))
	file_out = list(range(number_of_processes))

	# create filehandles
	for i in range(number_of_processes):
		file_in[i] = tempfile.NamedTemporaryFile(mode='w+t')
		file_out[i] = tempfile.NamedTemporaryFile(mode='w+t')

	# distribute data into different files
	line_counter = 0;
	if id_pos_res == {}:
		for id in id_seq:
			line_counter = line_counter + 1
			number = line_counter % number_of_processes
			file_in[number].write(">%s\n%s\n"%(id, id_seq[id]))
		#file = open(args[1], 'r')
		#for line in file:
		#	if( re.search('^>',line) ):
		#		line_counter = line_counter + 1
		#	number = line_counter % number_of_processes
		#	file_in[number].write(line)
	else:
		for id in id_seq:
			if id in id_pos_res:
				line_counter = line_counter + 1
				number = line_counter % number_of_processes
				file_in[number].write(">%s\n%s\n"%(id, id_seq[id]))
	for i in range(number_of_processes):
		file_in[i].flush()

	# run NetPhorest for each file in parallel
	if options.verbose:
		sys.stderr.write("\n")
	for i in range(number_of_processes):
		if options.verbose:
			sys.stderr.write(netphorest_bin+' < '+file_in[i].name+ ' > '+file_out[i].name+'\n');
		if save_diskspace:
			arg = (netphorest_bin+' < '+file_in[i].name+ '| gzip -9 > '+file_out[i].name,)
		else:
			arg = (netphorest_bin+' < '+file_in[i].name+ ' > '+file_out[i].name,)
		threading.Thread(target=myPopen, args=arg).start()
	if options.verbose:
		sys.stderr.write("Running on %s sequences\n"%line_counter)

	# wait for threads to finish
	while (threading.active_count() > 1):
		sys.stderr.write('.')
		time.sleep(5)

	return file_out
'''

# new runNetPhorest:

def runNetPhorest(id_seq, id_pos_res, save_diskspace, number_of_processes):
	#print(netphorest_bin)
	id_pos_tree_pred = {}
	number_of_sequences = len(id_seq) if not id_pos_res else len(id_pos_res)
	number_of_processes = min(number_of_sequences, number_of_processes)

	file_in = [tempfile.NamedTemporaryFile(mode='w+t') for _ in range(number_of_processes)]
	file_out = [tempfile.NamedTemporaryFile(mode='w+t') for _ in range(number_of_processes)]

	# Distribute data into different files
	line_counter = 0
	for id in id_seq:
		if not id_pos_res or id in id_pos_res:
			line_counter += 1
			number = line_counter % number_of_processes
			file_in[number].write(">%s\n%s\n" % (id, id_seq[id]))

	for f in file_in:
		f.flush()

	# Define function to run NetPhorest
	def run_netphorest_single(input_file, output_file):
		netphorest_command = f"{netphorest_bin} < {input_file.name}"
		if save_diskspace:
			netphorest_command += " | gzip -9"
		with open(output_file.name, 'w') as out:
			subprocess.run(netphorest_command, shell=True, stdout=out)

	# Run NetPhorest for each file in parallel
	processes = []
	for i in range(number_of_processes):
		p = multiprocessing.Process(target=run_netphorest_single, args=(file_in[i], file_out[i]))
		p.start()
		processes.append(p)

	if options.verbose:
		sys.stderr.write(f"Running on {line_counter} sequences\n")

	# Wait for processes to finish
	for p in processes:
		p.join()

	return file_out


# Parse the NetPhorest output
# expected format:
# Name  Position        Residue Peptide Method    Orgnism    Tree    Classifier      Posterior       Prior
# O00151  2       T       ----MtTQQID     nn    human    KIN     CDK2_CDK3_group    0.040050    0.028284
def parseNetphorestFile(filename, id_pos_res, save_diskspace):

	if save_diskspace:
		command = "gzip -cd %s"%(filename)
	else:
		command = "cat %s"%(filename)

	if options.verbose:
		sys.stderr.write('Parsing NetPhorest result file: "' + filename + '"\n');
	try:
		netphorest_results = myPopen(command)
	except:
		sys.stderr.write("Going to sleep for 1 hour to give you time to debug:\nCrashed with '%s'\n"%command)
		#time.sleep(3600)

	id_pos_tree_pred = {}
	netphorest_results=netphorest_results.split('\n')
	netphorest_results.pop(-1)
	for line in netphorest_results:
		if(re.match('^#',line)):
			continue
		tokens = line.split('\t')
		id = tokens[0]
		pos = int(tokens[1])
		# NetPhorest2 introduces organism column
		(res, peptide, method, organism, tree, pred) = tokens[2:8]
		try:
			if tree not in limitTrees:
				continue
		except:
			pass
		try:
			score = float(tokens[8])
		except:
			continue
		if (id in id_pos_res and pos in id_pos_res[id]) or id_pos_res == {}:
			if id in id_pos_tree_pred:
				if pos in id_pos_tree_pred[id]:
					if tree in id_pos_tree_pred[id][pos]:
						id_pos_tree_pred[id][pos][tree][pred] = (res,peptide,score)
					else:
						id_pos_tree_pred[id][pos][tree] = { pred: (res,peptide,score) }
				else:
					id_pos_tree_pred[id][pos]= { tree: { pred: (res,peptide,score) } }					
			else:
				id_pos_tree_pred[id]= { pos: { tree: { pred: (res,peptide,score) } } }
		else:
			pass
	else:
		pass
	return id_pos_tree_pred

# Map incoming peptides to STRING sequences
def mapPeptides2STRING(blastDir, organism, fastafilename, id_pos_res, id_seq, number_of_processes, datadir):
	sys.stderr.write("Mapping using blast\n")
	incoming2string = {}
	string2incoming = {}

	# Speedup, only blast sequences with site specified
	with open('tmp/blast_tmpfile.txt','w') as blast_tmpfile:
		#blast_tmpfile = tempfile.NamedTemporaryFile()
		if id_pos_res == {}:
			for id in id_seq:
				blast_tmpfile.write('>'+id+'\n'+id_seq[id]+'\n')
		else:
			for id in id_pos_res:
				try:
					blast_tmpfile.write('>'+id+'\n'+id_seq[id]+'\n')
				except:
					sys.stderr.write("No sequence available for '%s'\n"%id)
		blast_tmpfile.flush()

		blastDB = "%s/%s.protein.sequences.v12.0.fa"%(datadir,organism)

		# Check if blast database is actually initialized, if not: do it
		#if not os.path.isfile(blastDB+'.pin'):
		#	command = "%s/formatdb -i %s"%(blastDir.rsplit("/", 1)[0], blastDB)
		#	sys.stderr.write("Looks like blast database is not initialized, trying to run:\n%s\n"%command)
		#	myPopen(command)

		#command = "%s -a %s -p blastp -e 1e-10 -m 8 -d %s -i %s | sort -k12nr"%(blastDir, number_of_processes, blastDB, blast_tmpfile.name)
		#command = 'blastp -evalue 1e-10 -outfmt 6 -db NetworKIN_CODE_v3.0/data/9606.protein.sequences.fa -query tmp/tmpvkhu68f5 | sort -k12nr'
	command = f'blastp -evalue 1e-10 -outfmt 6 -db {blastDB} -query tmp/blast_tmpfile.txt | sort -k12nr'

	blast_out = myPopen(command)
	#print(command)
	#print(blast_out)
	blast_out = blast_out.split('\n')
	blast_out.pop(-1)
	#print(blast_out[-2])
	for line in blast_out:
		#print(line)
		tokens = line.split('\t')
		incoming = tokens[0]

		if incoming not in incoming2string:
			# get rid of organism prefix
			string = tokens[1].replace("%s."%organism, "")
			identity = float(tokens[2])
			evalue = float(tokens[-2])

			if string in string2incoming:
				sys.stderr.write('Best hit for '+incoming+' is not reciprocal\n')
			if identity < 90:
				sys.stderr.write('Best hit for '+incoming+' has only '+str(round(identity, 2))+' %identity\n')
			if evalue > 1e-40:
				sys.stderr.write('Best hit for '+incoming+' has high E-value '+str(round(evalue, 2))+' \n')
			if incoming in incoming2string:
				incoming2string[incoming][string] = True
			else:
				incoming2string[incoming] = { string: True }
			if string in string2incoming:
				string2incoming[string][incoming] = True
			else:
				string2incoming[string] = { incoming: True }
		else:
			pass
	#print(incoming2string, string2incoming)
	return incoming2string, string2incoming

# Random mapping of incoming identifiers to string
def mapRandom(id_seq):
	sys.stderr.write("Mapping random\n")
	import random

	incoming2string = {}
	string2incoming = {}

	stringIDs = []
	file = open('%s/%s.protein.sequences.fa'%(datadir,organism), 'r')
	data = file.readlines()
	file.close()

	for line in data:
		if line[0] == '>':
			name = line[1:-1]
			stringIDs.append(name)
			string2incoming[name] = {}
	max = len(stringIDs) - 1

	for incoming in id_seq:
		if incoming not in incoming2string:
			int = random.randint(0, max)
			string = stringIDs[int]
			incoming2string[incoming] = { string: True }
			if string in string2incoming:
				string2incoming[string][incoming] = True
			else:
				string2incoming[string] = { incoming: True }
	return incoming2string, string2incoming

# In case we run NetworKIN on the same sequence set we use in STRING, we can skip the mapping by blast
def mapOne2one(id_seq):
	sys.stderr.write("Mapping one2one\n")
	incoming2string = {}
	string2incoming = {}
	for incoming in id_seq:
		incoming2string[incoming] = {incoming: True}
		string2incoming[incoming] = {incoming: True}
	return incoming2string, string2incoming

# In case we have a better mapping then we can expect from blasting, we can use an external file to do so
# file format:
# incoming ID -> STRING ID
def mapFromFile(filename):
	sys.stderr.write("Mapping using external mapping file\n")
	incoming2string = {}
	string2incoming = {}

	command = "cat %s"%(filename)
	try:
		mappingFile = myPopen(command)
	except:
		sys.stderr.write("Going to sleep, crashed with '%s'\n"%command)
		#time.sleep(3600)

	for line in mappingFile:
		if(re.match('^#',line)):
			continue
	
		line = line.strip()
		tokens = line.split('\t')

		incoming = tokens[0]
		string = tokens[1]

		if incoming in incoming2string:
			incoming2string[incoming][string] = True
		else:
			incoming2string[incoming] = { string: True }
		if string in string2incoming:
			string2incoming[string][incoming] = True
		else:
			string2incoming[string] = { incoming: True }

	return incoming2string, string2incoming

# Load the precalculated STRING network file
def loadSTRINGdata(string2incoming, datadir, number_of_processes):
	# bestpath ( <= v9.0):
	#command = 'gzip -cd %s/%s.bestpath.tsv.gz'%(datadir, organism)
	# protein.links.txt:
	command = 'gzip -cd %s/%s.protein.links.v12.0.txt.gz'%(datadir, organism)
	try:
		data = myPopen(command)
	except:
		sys.stderr.write("Error loading STRING data using '%s', sleeping fo 1h.\n"%command)
		#time.sleep(3600)

	tree_pred_string_data = {}
	data=data.split('\n')
	data.pop(-1)

	for line in data:
		#print(line)
		#
		# for bestpath :
		
		try:
			(tree, pred, name, string1, string2, stringscore, null, path) = line.split('\t')
		# path to itself,  we will miss the path information
		except:
			(tree, pred, name, string1, string2, stringscore, null) = line.split('\t')
			path = ""
		
		'''
		# for protein.links.txt:
		line = line.split()
		if line[0].startswith('protein'):
			continue
		#(string1, string2, stringscore) = line.split('\t')
		string1=line[0]
		string2=line[1]
		stringscore=line[2]
		#print((string1, string2, stringscore))
		string1=string1.replace(f'{organism}.','')#[len(str(organism)+'.'):]
		string2=string2.replace(f'{organism}.','')#[len(str(organism)+'.'):]
		stringscore=float(stringscore)/1000
		tree='notdef'
		pred='notdef'
		name='notdef'
		path='notdef'
		'''
		'''
		try:
			if tree not in limitTrees:
				continue
		except:
			pass
		'''

		if string2 in string2incoming.keys():
			#print(string2)
			if tree in tree_pred_string_data:
				if pred in tree_pred_string_data[tree]:
					if string2 in tree_pred_string_data[tree][pred]:
						tree_pred_string_data[tree][pred][string2][string1] = {"_name": name}
					else:
						tree_pred_string_data[tree][pred][string2] = {string1: {"_name": name}}
				else:
					tree_pred_string_data[tree][pred] = {string2: {string1: {"_name": name}}}
			else:
				tree_pred_string_data[tree] = {pred: {string2: {string1: {"_name": name}}}}
			tree_pred_string_data[tree][pred][string2][string1]["_score"] = float(stringscore)
			tree_pred_string_data[tree][pred][string2][string1]["_path"] = path
		else:
			pass
	#print(tree_pred_string_data)
	return tree_pred_string_data

def printResult(id_pos_tree_pred, tree_pred_string_data, incoming2string, string_alias, string_desc,string_name, organism, mode,res_dir, fasta_file):
	ALPHA = ALPHAS[organism]
	csv_filename = os.path.join(res_dir, f'{fasta_file}.result.csv')
	with open(csv_filename, 'w', newline='') as csvfile:
		writer = csv.writer(csvfile, delimiter=',')
		'''
		writer.writerow(
			['Name', 'Position',  'id', 'NetworKIN score','Tree','NetPhorest Group',
			 'NetPhorest score', 'STRING score', 'Input STRING ID', 'Substrate STRING ID', 'Input description',
			 'Substrate description', 'Input Name', 'Substrate Name', 'Sequence window', 'Nodes'])
		'''
		writer.writerow(
			['substrate',	'position',	'id',	'networkin_score',	'tree',
			 'netphorest_group',	'netphorest_score',	'string_identifier',	'string_score',
			 'substrate_name',	'sequence',	'string_path'])

	# For each ID in NetPhorest
	for id in id_pos_tree_pred:
		# We have a mapping to STRING
		if id in incoming2string:
			# For each predicted position
			for pos in id_pos_tree_pred[id]:

				# For each of the trees (KIN, SH@ etc.)
				for tree in id_pos_tree_pred[id][pos]:
					#print(tree)
					score_results = {}
					# For each single classifier
					for pred in id_pos_tree_pred[id][pos][tree]:
						#print(pred)
						# For each mapped sequence
						for string1 in incoming2string[id]:
							if string1 in string_alias:
								bestName1 = string_alias[string1]
							else:
								bestName1 = ''
							if string1 in string_desc:
								desc1 = string_desc[string1]
							else:
								desc1 = ''

							if tree in tree_pred_string_data and pred in tree_pred_string_data[tree] and string1 in tree_pred_string_data[tree][pred]:

							#if string1 in tree_pred_string_data['notdef']['notdef']:
								(res, peptide, netphorestScore) = id_pos_tree_pred[id][pos][tree][pred]
								for string2 in tree_pred_string_data[tree][pred][string1]:
								#for string2 in tree_pred_string_data['notdef']['notdef'][string1]:
									if string2 in string_alias:
										bestName2 = string_alias[string2]
									else:
										bestName2 = ''
									if string2 in string_desc:
										desc2 = string_desc[string2]
									else:
										desc2 = ''
									stringScore = tree_pred_string_data[tree][pred][string1][string2]["_score"]
									path = tree_pred_string_data[tree][pred][string1][string2]["_path"]
									name = tree_pred_string_data[tree][pred][string1][string2]["_name"]
									#stringScore = tree_pred_string_data['notdef']['notdef'][string1][string2]["_score"]
									#path = tree_pred_string_data['notdef']['notdef'][string1][string2]["_path"]
									#name = tree_pred_string_data['notdef']['notdef'][string1][string2]["_name"]
									networkinScore = pow(stringScore, ALPHA)*pow(netphorestScore, 1-ALPHA)

									# NetworKIN result
									'''
									result = id+'\t'+res+str(pos)+'\t'+tree+'\t'+pred+'\t'+name+'\t'+ \
									str(round(networkinScore,4))+'\t'+str(round(netphorestScore,4))+'\t'+str(round(stringScore,4))+'\t'+ \
									string1+'\t'+string2+'\t'+bestName1+'\t'+bestName2+'\t'+desc1+'\t'+desc2+'\t'+peptide+'\t'+path+'\n'
									'''
									try:
										result = string1+'\t'+str(pos)+'\t'+string_name[organism+'.'+string2]+'\t'+str(round(networkinScore,4))+'\t'+tree+'\t'+ \
										pred+'\t'+str(round(netphorestScore,4))+'\t'+string2+'\t'+str(round(stringScore,4))+'\t'+string_name[organism+'.'+string1]+'\t'+ \
										peptide+'\t'+path+'\n'
									except:
										result = string1 + '\t' + str(pos) + '\t' + name + '\t' + str(round(networkinScore, 4)) + '\t' + tree + '\t' + \
										pred + '\t' + str(round(netphorestScore, 4)) + '\t' + string2 + '\t' + str(round(stringScore, 4)) + '\t' + bestName2 + '\t' + \
										peptide + '\t' + path + '\n'
									

									'''
									result = id+'\t'+str(pos)+str(round(networkinScore,4))+'\t'+tree+'\t'+pred+'\t'+name+'\t'+ \
									'\t'+str(round(netphorestScore,4))+'\t'+str(round(stringScore,4))+'\t'+ \
									string1+'\t'+string2+'\t'+bestName1+'\t'+bestName2+'\t'+desc1+'\t'+desc2+'\t'+peptide+'\t'+path+'\n'
									'''

									with open(csv_filename, 'a', newline='') as csvfile:
										writer = csv.writer(csvfile, delimiter=',')
										writer.writerow(result.split('\t'))

									if networkinScore not in score_results:
										score_results[networkinScore] = []

									score_results[networkinScore].append(result)
					'''
					if mode == 'network':
						highestScore = sorted(score_results.keys(), reverse=True)[0]

						if len(score_results[highestScore]) > 1:
							index = random.randint(0,len(score_results[highestScore])-1)
							sys.stdout.write((score_results[highestScore][index]))
						else:
							sys.stdout.write((score_results[highestScore][0]))
						pass
					else:
						for score in sorted(score_results.keys(), reverse=True):
							sys.stdout.write("".join(score_results[score]))
					'''
		else:
			pass
	return

#MAIN
def Main():

	#for key, value in os.environ.items():
	#	print(f"{key}: {value}")

	sys.stderr.write("Reading fasta input file\n")
	id_seq = readFasta(fastafile)
	if options.verbose:
		sys.stderr.write("%s sequences loaded\n"%len(id_seq.keys()))

	if sitesfile:
		sys.stderr.write("Reading phosphosite file\n")
		id_pos_res = readPhosphoSites(sitesfile)
	else:
		id_pos_res = {}

	sys.stderr.write("Loading aliases and descriptions\n")
	(string_alias, string_desc, string_name) = readAliasFiles(args[0], options.datadir);

	# Default way of mapping using BLAST
	incoming2string, string2incoming = mapPeptides2STRING(blastDir, organism, fastafile.name, id_pos_res, id_seq, options.threads, options.datadir)

	# Hack for random mapping to proteins
	#incoming2string, string2incoming = mapRandom(id_seq)

	# Use if a mapping file for the input can be provided
	#incoming2string, string2incoming = mapFromFile("/home/red1/hhorn/projects2/2012_03_22_Jesper/ensmusp_ensp.tsv")

	# Used if only ensembl of the same version used
	#incoming2string, string2incoming = mapOne2one(id_seq)

	# Load the STRING network data
	sys.stderr.write("Loading STRING network\n")
	tree_pred_string_data = loadSTRINGdata(string2incoming, options.datadir, options.threads)


	# Run NetPhorest
	sys.stderr.write("Running NetPhorest")
	netphorestTmpFiles = runNetPhorest(id_seq, id_pos_res, options.compress, options.threads)
	sys.stderr.write('\n')

	# Writing result to STDOUT
	sys.stderr.write("Writing results\n")
	sys.stdout.write("#Name\tPosition\tTree\tNetPhorest Group\tKinase/Phospho-binding domain\t\
									 NetworKIN score\tNetPhorest score\tSTRING score\t\
									 Input STRING ID\tSubstrate STRING ID\tInput description\tSubstrate description\tInput Name\tSubstrate Name\tSequence window\tNodes\n")

	res_dir='results'
	if not os.path.exists(res_dir):
		os.makedirs(res_dir)

	#print(string_name)
	for i in range(len(netphorestTmpFiles)):
		id_pos_tree_pred = parseNetphorestFile(netphorestTmpFiles[i].name, id_pos_res, options.compress)
		#tree_pred_string_data['notdef']={'notdef':'notdef'}
  		
		printResult(id_pos_tree_pred, tree_pred_string_data, incoming2string, string_alias, string_desc, string_name, args[0], options.mode,res_dir,args[1])
 
	return

if __name__ == '__main__':
	#BLAST
	try:
		blastDir = os.environ['BLAST_PATH']
		
	except:
		blastDir=""
	#NETPHOREST
	try:
		netphorest_bin = "" #os.environ['NETPHOREST_PATH']
	except:
		netphorest_bin=""

	usage = "usage: %prog [options] organism FASTA-file [sites-file]"
	parser = OptionParser(usage=usage, version="%prog 3.0")
	parser.add_option("-n", "--netphorest", dest="netphorest_bin", default=netphorest_bin,
										help="set the location of the NetPhorest binary, overwrites the 'NETPHOREST_PATH' environmental variable. [ENV: %default]")
	parser.add_option("-b", "--blast", dest="blast", default=blastDir,
										help="set the directory for the BLAST binaries (formatdb and blastall), overwrites the 'BLAST_PATH' environmental variable. [ENV: %default]")
	parser.add_option("-m", "--mode", dest="mode", default=False,
										help="if set to 'network', gives only one best scoring result for each site. In case of multiple candidate kinases with the same core, the selection hapens randomly. [default: %default]")
	parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
										help="print out everything [default: %default]")
	parser.add_option("-t", "--threads", dest="threads", default=1, type="int",
										help="number of available threads/CPUs. Also leads to less memory usage, as result files are read sequentially [default: %default]")
	parser.add_option("-c", "--compress", dest="compress", default=True,
										help="compress temporary result files, saves discspace [default: %default]")
	parser.add_option("-d", "--data", dest="datadir", default=sys.argv[0].rsplit("/", 1)[0]+'/data',
										help="location for the additional files like the pre-computed STRING network, STRING sequence database etc. [default: %default]")
	parser.add_option("--tmp", dest="tmpdir", default="/tmp",
										help="location for the temporary files [default: %default]")
	#python3 NetworKIN_CODE_v3.0/NetworKIN.py -n /netphorest -b /ncbi-blast-2.15.0+/bin -d NetworKIN_CODE_v3.0/data 9606 test.fas
	#python3 NetworKIN_CODE_v3.0/NetworKIN.py -n netphorest/netphorest -d NetworKIN_CODE_v3.0/data 9606 test1.fas
	#python3 NetworKIN_CODE_v3.0/NetworKIN.py -n netphorest/netphorest -d NetworKIN_CODE_v3.0/data 9606 cured_morpho_seqs_v2.fa
	global options
	(options, args) = parser.parse_args()

	tempfile.tempdir= options.tmpdir

	#ORGANISM
	try:
		organism = int(args[0])
	except:
		parser.error("Organism not defined!")

	#SEQUENCE FILE
	try:
		fastafile = open(args[1], 'r')
	except:
		sys.stderr.write("%s"%args)
		parser.error("FASTA-file not defined!")

	#SITES FILE
	try:
		sitesfile = open(args[2], 'r')
	except:
		sitesfile = False

	#BLAST
	if options.blast:
		blastDir = options.blast
	
	#NETPHOREST
	if options.netphorest_bin:
		netphorest_bin = options.netphorest_bin
		
	# Show runtime parameters
	if(options.verbose):
		sys.stderr.write('\nPredicting using parameters as follows:\nOrganism:\t%s\nFastaFile:\t%s\n'%(organism, args[1]))
		sys.stderr.write('Threads:\t%s\nCompress:\t%s\n'%(options.threads, options.compress))
		if sitesfile:
			sys.stderr.write('Sitesfile:\t%s\n'%sitesfile)
		else:
			sys.stderr.write('No sites-file given, predicting on all S,T,Y residues.\n')
		if options.mode:
			sys.stderr.write('Mode:\t\t%s\n'%options.mode)
		sys.stderr.write("Blast dir: %s" % blastDir)
		sys.stderr.write("NetPhorest binary: %s" % netphorest_bin)
		sys.stderr.write('\n')
	Main()
