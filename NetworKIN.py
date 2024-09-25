#!/usr/bin/env python3

"""
NetworKIN(tm), (C) 2005,2006,2007,2013.
Drs Rune Linding, Lars Juhl Jensen, Heiko Horn & Jinho Kim

Usage: ./networkin.py Organism FastaFile SitesFile

If no sites file is given NetworKIN will predict on all T/S/Y residues
in the given sequences.
"""

### Coding strategy & TODO
# FIX: filtering when no sites given
# 2) add code for prediction of downstream recognition site (14-3-3 first)
# 3) add code for predicting downstream binding module (14-3-3 first)
# 4) finnish up output format and CLI options
# 5) add code for doing 'pathway' analysis
# 6) predict phenotypes
###

# Changelog
# 21.10.06: tested blast -a 8 option, no differences on proteome
# 21.10.06: Filtercode broken
#			-fixed
# 28.01.07: Testing autophosphorylation (self == 1 in update script)
# 30.07.07: Working on v1.5 milestone
# 07.05.08: Initiated 2.5 w. scaling factor
# 11.07.13: New scoring scheme (Bayesian)
# 03.02.14: Several minor bug fixes. blast sorting bug fix (k12nr -> k12gr)
#


###   python3 NetworKIN.py -n netphorest/netphorest -b /mnt/c/Program Files/NCBI/blast-2.16.0+/bin/blastp.exe -d data 9606 test1.fas test1.tsv
###   mv /home/fandrich/projekte/networkin_bayesian/results/KinomeXplorer_all_predictions_v2.csv /home/basar/Signaling_Group/
###   mv /home/fandrich/projekte/networkin_bayesian/results/cured_morpho_seqs_v2.fa.result.csv /home/basar/Signaling_Group/


import sys, os, subprocess, re, tempfile, random, operator
import threading
from optparse import OptionParser
import multiprocessing
import glob
import csv
from string import *
from likelihood import ReadConversionTableBin
from likelihood import ConvertScore2L
import platform

# debugging
import time

# Weighting parameter, 0=only motif, 1=only STRING
# estimated while benchmarking an then hardcoded here
# feel free to play with it, but it is your own responsibility
ALPHAS = {"9606": 0.85, "4932": 0.65}
dSpeciesName = {"9606": "human", "4932": "yeast"}
dPenalty = {"9606": {"hub penalty": 100, "length penalty": 800}, "4932": {"hub penalty": 170, "length penalty": 1000}}

NETWORKIN_SITE_FILE = 1
PROTEOME_DISCOVERER_SITE_FILE = 2
MAX_QUANT_DIRECT_OUTPUT_FILE = 3
RUNES_SITE_FILE = 4

global options


# Temporary directory to store the files
# tempfile.tempdir= '/tmp'

# Location of NetworKIN input files (string network, alias file etc.)
# datadir = sys.argv[0].rsplit("/", 1)[0]+'/data'
# global DATADIR
# DATADIR = ""

# Number of threads used for NetPhorest and BLAST
# setting it to a high number can also help in case NetworKIN uses too much memory
# as NetPhorest result files will be read one after the other to save memory
# -> the more files, the less memory usage
# global NUMBER_OF_PROCESSES
# NUMBER_OF_PROCESSES = "1";

# should temporary files of netphorest be zipped?
# saves some diskspace and makes it sometimes faster (depending on the CPU/harddisk speed)
# global SAVE_DISKSPACE
# SAVE_DISKSPACE = "";

# Should the analysis be limited to specific trees?
# just comment it out if you want to predict on all
# limitTrees = ['SH2']
# limitTrees = ['KIN', "SH2", "PTP"]
# limitTrees = ['KIN']

################################################################################
#                                                                              #
#                             the code starts here                             #
#                                                                              #
################################################################################

class CSheet(list):
    pass


'''
# Run system binary
def myPopen(cmd):
	try:
		pipe = subprocess.Popen(cmd, shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout = pipe.stdout.readlines()
	except:
		sys.stderr.write('ERROR executing: '+`cmd`+'\n')
		sys.exit()
	else:
		return stdout 
'''


### chat-gpt's myPopen:
def myPopen(cmd):
    if platform.system() == 'Windows':
        sys.stderr.write('WINDOWS\n')
        # Wrap the command to run in a Unix-like shell using WSL
        command = f'wsl bash -c "{cmd}"'
    else:
        command = cmd
    try:
        pipe = subprocess.Popen(command, shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = pipe.communicate()
        # Decode the byte strings to utf-8
        stdout = stdout.decode('utf-8')
        stderr = stderr.decode('utf-8')
        if pipe.returncode != 0:
            # If the command failed, print the stderr and raise an exception
            error_message = 'ERROR executing: ' + repr(command) + '\n' + stderr
            sys.stderr.write(error_message + '\n')
            raise subprocess.CalledProcessError(pipe.returncode, command, output=error_message)
        else:
            # If the command succeeded, return the stdout
            return stdout
    except subprocess.CalledProcessError as e:
        # Handle the subprocess.CalledProcessError exception
        print("Error: Command '{}' returned non-zero exit status {}".format(e.cmd, e.returncode))
        return e.output  # Return the error message
    except Exception as e:
        # Handle other exceptions
        error_message = 'ERROR executing: ' + repr(command) + '\n' + str(e) + '\n'
        sys.stderr.write(error_message)
        return error_message


# Read sequences from fasta file
def readFasta(fastafile):
    id_seq = {}
    aminoacids = re.compile('[^ACDEFGHIKLMNPQRSTVWYXB]')
    data = fastafile.readlines()
    fastafile.close()
    seq = ''
    for line in data:
        if line[0] != ';':
            if line[0] == '>':
                line = line.strip()
                # id = line[1:].split(' ', 1)[0]
                id = line[1:].split('_', 1)[0]
                #				id = line[1:-1].split(' ', 1)[0]
                seq = ''
            else:
                seq += aminoacids.sub('', line)
                if len(seq) > 0:
                    id_seq[id] = seq
        else:
            continue
    return id_seq


def CheckInputType(sitesfile):
    f = open(sitesfile, 'rU')
    line = f.readline()
    f.close()

    tokens = line.split()
    if len(tokens) == 3:
        return NETWORKIN_SITE_FILE
    elif len(tokens) == 2:
        return PROTEOME_DISCOVERER_SITE_FILE
    elif tokens[0] == "Proteins" and tokens[4] == "Leading":
        return MAX_QUANT_DIRECT_OUTPUT_FILE
    elif tokens[1] == "phospho":
        return RUNES_SITE_FILE
    else:
        sys.stderr("Unknown format of site file")
        sys.exit()


# Read phosphorylation sites from tsv file
# id -> position -> residue
def readPhosphoSites(sitesfile):
    id_pos_res = {}
    f = open(sitesfile, 'rU')
    if f:
        data = f.readlines()
        f.close()
        for line in data:
            tokens = line.split('\t')
            id = tokens[0]
            try:
                pos = int(tokens[1])
            except:
                sys.stderr.write(line)
                raise
            try:
                res = tokens[2][:-1]
            except:
                res = ""
            if id in id_pos_res:
                id_pos_res[id][pos] = res
            else:
                id_pos_res[id] = {pos: res}
    else:
        sys.stderr.write("Could not open site file: %s" % sitesfile)
        sys.exit()

    return id_pos_res


def readPhosphoSitesProteomeDiscoverer(fastafile, sitesfile):
    # This function takes a tab-separated list of ProteinID \t Peptide, in ProteomeDiscoverer format (i.e. phosphosites listed as
    # lowercase), maps the peptide to the full-length protein sequence and conducts the NetworKIN search.
    #

    # store all fasta sequences
    fasta = open(fastafile).readlines()
    fastadict = {}
    for line in fasta:
        if line[0] == ">":  # use this as each new entry in FASTA file starts with '>'
            ensp = line.split()[0].strip(">")[
                   0:15]  # strip the '>' symbol from the line before copying value into dictionary
            # print ensp
            fastadict[ensp] = ""
        # retrieve the AA sequence belonging to individual ESPN identifiers and store as value for ESPN key
        else:
            seq = line.strip()
            fastadict[ensp] = seq

    # store all peptides in dictionary, with parent protein as key
    peptides = open(sitesfile).readlines()
    peptidedict = {}
    for line in peptides:
        line = line.strip()
        tokens = line.split("\t")
        protID = tokens[0]
        peptide = tokens[1]
        if protID in peptidedict.iterkeys():
            if peptide in peptidedict[protID]:
                pass
            else:
                peptidedict[protID][peptide] = ""
        else:
            peptidedict[protID] = {peptide: ""}

    # now map peptide onto full length sequence, then get absolute phosphosite locations
    id_pos_res = {}
    for protID in peptidedict.iterkeys():
        for peptide in peptidedict[protID]:

            # first make peptide upper case to match to FASTA sequence
            UPPERpeptide = peptide.upper()
            # print UPPERpeptide

            # get parent protein sequence
            sequence = fastadict[protID]
            peptideindex = sequence.index(UPPERpeptide)
            # print peptideindex

            x = 0
            for letter in peptide:
                if letter.islower():
                    if letter == "s" or letter == "t" or letter == "y":

                        phoslocation = peptideindex + x + 1
                        phosresidue = sequence[phoslocation - 1]

                        if protID in id_pos_res.iterkeys():
                            id_pos_res[protID][phoslocation] = phosresidue
                        else:
                            id_pos_res[protID] = {phoslocation: phosresidue}

                    else:
                        # non phosho modification
                        pass
                else:
                    pass

                x += 1
    return id_pos_res


# return a list of dictionaries
# Usage l[index][column name]
def ReadSheet(fname, offset=0):
    l = CSheet()

    f = open(fname, 'rU')

    for i in range(offset):
        f.readline()

    columns = f.readline().strip().split('\t')

    l.columns = columns
    for line in f.readlines():
        instance = {}
        l.append(instance)
        fields = line.strip().split('\t')

        if len(fields) > len(columns):
            sys.stderr.write("Error in data file")
            sys.stderr.write(columns)
            sys.stderr.write(fields)
            raise

        for i in range(len(columns)):
            try:
                instance[columns[i]] = fields[i]
            except IndexError:
                instance[columns[i]] = ''

    # sys.stderr.write("No. of entries in file input: %d" % len(l))
    f.close()

    return l


def readPhosphoSitesMaxQuant(fname, only_leading=False):
    id_pos_res = {}

    phosphosites = ReadSheet(fname)

    for site in phosphosites:
        Ids = site["Proteins"].split(';')
        positions = [int(x) for x in site["Positions within proteins"].split(';')]
        aa = site["Amino acid"]

        leading_protein_ids = site["Leading proteins"].split(';')

        for i in range(len(Ids)):
            Id = Ids[i]
            pos = positions[i]

            if only_leading and not Id in leading_protein_ids:
                continue

            if Id in id_pos_res:
                id_pos_res[Id][pos] = aa
            else:
                id_pos_res[Id] = {pos: aa}

    return id_pos_res

def readRunessitesfile(sitesfile):
	id_pos_res = {}
	f = open(sitesfile, 'rU')
	if f:
		data = f.readlines()
		f.close()
		for line in data:
			tokens = line.split(' ')
			id = tokens[3]
			res_pos = tokens[2]
			try:
				pos = int(res_pos[2:])
			except:
				sys.stderr.write(line)
				raise
			try:
				res = res_pos[0]
			except:
				res = ""
			if id in id_pos_res:
				id_pos_res[id][pos] = res
			else:
				id_pos_res[id] = {pos: res}
	else:
		sys.stderr.write("Could not open site file: %s" % sitesfile)
		sys.exit()

	return id_pos_res


'''
#Alias hashes
def readAliasFiles(organism, datadir):
	alias_hash = {}
	desc_hash = {}

	# Read alias db
	try:
		alias_db = myPopen('gzip -cd %s/%s.alias_best.tsv.gz'%(datadir, organism))
		for line in alias_db:
			(taxID, seqID, alias) = line.strip().split('\t')[:3]
			alias_hash[seqID] = alias
	except:
		sys.stderr.write("No aliases available for organism: '%s'\n"%organism)

	# Read desc db
	try:
		desc_db = myPopen('gzip -cd %s/%s.text_best.tsv.gz'%(datadir, organism))
		for line in desc_db:
			(taxID, seqID, desc) = line.split('\t')[:3]
			desc_hash[seqID] = desc
	except:
		sys.stderr.write("No descriptions available for organism: '%s'\n"%organism)

	return alias_hash, desc_hash
'''


###  my Alias hashes
def readAliasFiles(organism, datadir):
    alias_hash = {}
    desc_hash = {}
    name_hash = {}
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
        for line in alias_db:
            if line.startswith('#'):
                continue
            (seqID, alias, source) = line.split('\t')
            if source == ('Ensembl_HGNC'):
                alias_hash[seqID[5:]] = alias
    except:
        sys.stderr.write("No aliases available for organism: '%s'\n" % organism)

    # Read desc db
    try:
        desc_db = myPopen('gzip -cd %s/%s.text_best.v9.0.tsv.gz' % (datadir, organism))
        # desc_db = myPopen('gzip -cd %s/%s.protein.info.v12.0.txt.gz'%(datadir, organism))
        desc_db = desc_db.split('\n')
        for line in desc_db:
            # print(line)
            (taxID, seqID, desc) = line.split('\t')[:3]
            desc_hash[seqID] = desc
    except:
        sys.stderr.write("No descriptions available for organism: '%s'\n" % organism)

    try:
        # desc_db = myPopen('gzip -cd %s/%s.text_best.v9.0.tsv.gz'%(datadir, organism))
        name_db = myPopen('gzip -cd %s/%s.protein.info.v12.0.txt.gz' % (datadir, organism))
        name_db = name_db.split('\n')
        for line in name_db:
            # print(line)
            (string_protein_id, preferred_name, protein_size, annotation) = line.split('\t')
            name_hash[string_protein_id] = preferred_name
    except:
        sys.stderr.write("No names available for organism: '%s'\n" % organism)

    # print(alias_hash)
    return alias_hash, desc_hash, name_hash


# Run Netphorest
def runNetPhorest(id_seq, id_pos_res, save_diskspace, number_of_processes, number_of_active_processes=3, fast=False,
                  leave_intermediates=False):
    if number_of_active_processes < 2:
        raise ValueError("Number of maximum threads is less than 2")
    id_pos_tree_pred = {}

    # check how many sequences we actually have
    number_of_sequences = 0
    if id_pos_res == {}:
        number_of_sequences = len(id_seq)
    else:
        number_of_sequences = len(id_pos_res)

    if number_of_sequences < number_of_processes:
        number_of_processes = number_of_sequences

    # use multiple instances of netphorest
    file_in = list(range(number_of_processes))
    file_out = list(range(number_of_processes))

    class CDummy:
        def __init__(self, name):
            self.name = name

    # create filehandles
    for i in range(number_of_processes):
        file_in[i] = tempfile.NamedTemporaryFile()

        if fast or leave_intermediates:
            if save_diskspace:
                file_out[i] = CDummy("%s.%s.gz" % (fn_netphorest_output, str(i)))  # jhkim

            else:
                file_out[i] = CDummy("%s.%s.txt" % (fn_netphorest_output, str(i)))
        else:
            file_out[i] = tempfile.NamedTemporaryFile()

    # distribute data into different files
    line_counter = 0;
    if id_pos_res == {}:
        for id in id_seq:
            line_counter = line_counter + 1
            number = line_counter % number_of_processes
            file_in[number].write((">%s\n%s\n" % (id, id_seq[id])).encode('utf-8'))
    # file = open(args[1], 'r')
    # for line in file:
    #	if( re.search('^>',line) ):
    #		line_counter = line_counter + 1
    #	number = line_counter % number_of_processes
    #	file_in[number].write(line)
    else:
        for id in id_seq:
            if id in id_pos_res:
                line_counter = line_counter + 1
                number = line_counter % number_of_processes
                file_in[number].write((">%s\n%s\n" % (id, id_seq[id])).encode('utf-8'))
    for i in range(number_of_processes):
        file_in[i].flush()

    # run NetPhorest for each file in parallel
    if options.verbose:
        sys.stderr.write("\n")

    if options.verbose:
        sys.stderr.write("Running on %s sequences\n" % line_counter)

    for i in range(number_of_processes):
        if fast and os.path.isfile(file_out[i].name):
            if options.verbose:
                sys.stderr.write("%s is already exist.\n" % file_out[i].name)
            continue

        while threading.active_count() >= number_of_active_processes:
            sys.stderr.write('a')
        # time.sleep(5)

        if options.verbose:
            sys.stderr.write(netphorest_bin + ' < ' + file_in[i].name + ' > ' + file_out[i].name + '\n');
        if save_diskspace:
            arg = (netphorest_bin + ' < ' + file_in[i].name + '| gzip -9 > ' + file_out[i].name,)
        else:
            arg = (netphorest_bin + ' < ' + file_in[i].name + ' > ' + file_out[i].name,)
        # print(arg)
        #####
        process = subprocess.Popen([netphorest_bin], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        output, error = process.communicate(input=open(file_in[i].name, 'rb').read())

        # Decode the output from bytes to string (assuming it's in UTF-8 encoding)
        output_str = output.decode('utf-8')

        # Print the output
        '''
		print("Output:")
		print(output_str)
		print(len(output_str.split('\n')))
		'''
        #####

        threading.Thread(target=myPopen, args=arg).start()

    # wait for threads to finish
    while (threading.active_count() > 1):
        pass
    # sys.stderr.write('.')
    # time.sleep(5)

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
'''


# Parse the NetPhorest output
# expected format:
# Name  Position        Residue Peptide Method    Orgnism    Tree    Classifier      Posterior       Prior
# O00151  2       T       ----MtTQQID     nn    human    KIN     CDK2_CDK3_group    0.040050    0.028284
def parseNetphorestFile(filename, id_pos_res, save_diskspace):
    if save_diskspace:
        command = "gzip -cd %s" % (filename)
    else:
        command = "cat %s" % (filename)

    if options.verbose:
        sys.stderr.write('Parsing NetPhorest result file: "' + filename + '"\n');
    try:
        netphorest_results = myPopen(command)
    except:
        sys.stderr.write("Going to sleep for 1 hour to give you time to debug:\nCrashed with '%s'\n" % command)
        time.sleep(3600)

    id_pos_tree_pred = {}

    # print(netphorest_results)
    netphorest_results = netphorest_results.split('\n')
    for line in netphorest_results:
        if (re.match('^#', line)) or (line == ''):
            continue
        tokens = line.split('\t')
        id = tokens[0]
        pos = int(tokens[1])
        # NetPhorest2 introduces organism column
        try:
            (res, peptide, method, organism, tree, pred) = tokens[2:8]
        except ValueError:
            sys.stderr.write(tokens)
            raise
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
                        id_pos_tree_pred[id][pos][tree][pred] = (res, peptide, score)
                    else:
                        id_pos_tree_pred[id][pos][tree] = {pred: (res, peptide, score)}
                else:
                    id_pos_tree_pred[id][pos] = {tree: {pred: (res, peptide, score)}}
            else:
                id_pos_tree_pred[id] = {pos: {tree: {pred: (res, peptide, score)}}}
        else:
            pass
    else:
        pass
    return id_pos_tree_pred


def ReadLines(fname):
    f = open(fname)
    lines = f.readlines()
    f.close()
    return lines


def WriteString(fname, s):
    f = open(fname, 'w')
    f.write(s)
    f.close()


### my mapPeptidestoString:

def mapPeptides2STRING(blastDir, organism, fastafilename, id_pos_res, id_seq, number_of_processes, datadir):
    sys.stderr.write("Mapping using blast\n")
    incoming2string = {}
    string2incoming = {}

    with open('tmp/blast_tmpfile.txt', 'w') as blast_tmpfile:
        if id_pos_res == {}:
            for id in id_seq:
                blast_tmpfile.write('>' + id + '\n' + id_seq[id] + '\n')
        else:
            for id in id_pos_res:
                try:
                    blast_tmpfile.write('>' + id + '\n' + id_seq[id] + '\n')
                except:
                    sys.stderr.write("No sequence available for '%s'\n" % id)
        blast_tmpfile.flush()

        blastDB = "%s/%s.protein.sequences.v12.0.fa" % (datadir, organism)

    # Check if blast database is actually initialized, if not: do it
    if not os.path.isfile(blastDB+'.pin'):
        command = f"{blastDir.rsplit('/', 1)[0]}/makeblastdb -in {blastDB} -out {blastDB} -parse_seqids -dbtype prot"
        sys.stderr.write(f"Looks like blast database is not initialized, trying to run:\n{command}\n")
        myPopen(command)

    if platform.system() == 'Windows':
        sys.stderr.write('WINDOWS\n')
        # Wrap the command to run in a Unix-like shell using WSL
        command_windows = f'blastp -evalue 1e-10 -outfmt 6 -db {blastDB} -query tmp/blast_tmpfile.txt | sort -k12nr'
        # command = f'wsl bash -c "{command_windows}"'
        command = command_windows
    else:
        # Use the Unix-like command as is
        command = f'blastp -evalue 1e-10 -outfmt 6 -db {blastDB} -query tmp/blast_tmpfile.txt | sort -k12nr'

    blast_out = myPopen(command)
    # print(command)
    # print(blast_out)
    blast_out = blast_out.split('\n')
    blast_out.pop(-1)
    # print(blast_out[-2])
    for line in blast_out:
        # print(line)
        tokens = line.split('\t')
        incoming = tokens[0]

        if incoming not in incoming2string:
            # get rid of organism prefix
            string = tokens[1].replace("%s." % organism, "")
            identity = float(tokens[2])
            evalue = float(tokens[-2])

            if string in string2incoming:
                sys.stderr.write('Best hit for ' + incoming + ' is not reciprocal\n')
            if identity < 90:
                sys.stderr.write('Best hit for ' + incoming + ' has only ' + str(round(identity, 2)) + ' %identity\n')
            if evalue > 1e-40:
                sys.stderr.write('Best hit for ' + incoming + ' has high E-value ' + str(round(evalue, 2)) + ' \n')
            if incoming in incoming2string:
                incoming2string[incoming][string] = True
            else:
                incoming2string[incoming] = {string: True}
            if string in string2incoming:
                string2incoming[string][incoming] = True
            else:
                string2incoming[string] = {incoming: True}
        else:
            pass
    # print(incoming2string, string2incoming)
    return incoming2string, string2incoming


# Random mapping of incoming identifiers to string
def mapRandom(id_seq):
    sys.stderr.write("Mapping random\n")
    import random

    incoming2string = {}
    string2incoming = {}

    stringIDs = []
    file = open('%s/%s.protein.sequences.fa' % (datadir, organism), 'r')
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
            incoming2string[incoming] = {string: True}
            if string in string2incoming:
                string2incoming[string][incoming] = True
            else:
                string2incoming[string] = {incoming: True}
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

    command = "cat %s" % (filename)
    try:
        mappingFile = myPopen(command)
    except:
        sys.stderr.write("Going to sleep, crashed with '%s'\n" % command)
        time.sleep(3600)

    for line in mappingFile:
        if (re.match('^#', line)):
            continue

        line = line.strip()
        tokens = line.split('\t')

        incoming = tokens[0]
        string = tokens[1]

        if incoming in incoming2string:
            incoming2string[incoming][string] = True
        else:
            incoming2string[incoming] = {string: True}
        if string in string2incoming:
            string2incoming[string][incoming] = True
        else:
            string2incoming[string] = {incoming: True}

    return incoming2string, string2incoming


# Load the precalculated STRING network file
def loadSTRINGdata(string2incoming, datadir, number_of_processes):
    # command = 'gzip -cd %s/%s.string_data.tsv.gz'%(datadir, organism)

    # fn_bestpath = "%s/%s.string_000_%04d_%04d.tsv.gz" % (os.path.join(datadir, "string_data"), organism, dPenalty[organism]["hub penalty"], dPenalty[organism]["length penalty"])
    fn_bestpath = "%s/%s.links.v12.0.tsv.gz" % (os.path.join(datadir, "string_data"), organism)

    if not os.path.isfile(fn_bestpath):
        sys.stderr.write("Best path file does not exist: %s" % fn_bestpath)

    '''	    
	command = "gzip -cd %s" % fn_bestpath

	try:
		data = myPopen(command)
	except:
		sys.stderr.write("Error loading STRING data using '%s', sleeping fo 1h.\n"%command)
		time.sleep(3600)
	'''

    tree_pred_string_data = {}

    # for memory efficiency
    import gzip
    f = gzip.open(fn_bestpath)
    f = f.read().decode('utf-8')
    f = f.split('\n')
    for line in f:
        # print(line)
        if line == '':
            continue
        # line = line.strip()
        tokens = line.split('\t')
        # print(tokens)
        if len(tokens) == 9:
            (tree, group, name, string1, string2, stringscore, stringscore_indirect, path) = tokens
        elif len(tokens) == 6:
            (tree, group, name, string1, string2, stringscore) = tokens
            stringscore = float(stringscore) / 1000
            string1 = string1.replace(f'{organism}.', '')
            string2 = string2.replace(f'{organism}.', '')
        '''
		if len(tokens) == 8:
			(tree, group, name, string1, string2, stringscore, stringscore_indirect, path) = tokens
		elif len(tokens) == 7:
			(tree, group, name, string1, string2, stringscore, stringscore_indirect) = tokens
			path = ""	# path to itself,  we will miss the path information
		elif len(tokens) == 6:
			(name, string1, string2, stringscore, stringscore_indirect, path) = tokens
		elif len(tokens) == 5:
			(name, string1, string2, stringscore, stringscore_indirect) = tokens
			path = ""	# path to itself,  we will miss the path information
		'''

        if string2 in string2incoming:
            if string2 in tree_pred_string_data:
                tree_pred_string_data[string2][string1] = {"_name": name}
            else:
                tree_pred_string_data[string2] = {string1: {"_name": name}}

            if options.path == "direct":
                tree_pred_string_data[string2][string1]["_score"] = float(stringscore)
            elif options.path == "indirect":
                tree_pred_string_data[string2][string1]["_score"] = float(stringscore_indirect)  # Use indirect path
            else:
                raise "Path information should be either direct or indirect."
            if len(tokens) == 9:
                tree_pred_string_data[string2][string1]["_path"] = path
            elif len(tokens) == 6:
                tree_pred_string_data[string2][string1]["_path"] = 'notdef'
        else:
            pass

    # f.close()
    return tree_pred_string_data


def InsertValueIntoMultiLevelDict(d, keys, value):
    for i in range(len(keys) - 1):
        # if not d.has_key(keys[i]):
        if not (keys[i] in d.keys()):
            d[keys[i]] = {}
        d = d[keys[i]]

    if not (keys[-1] in d.keys()):
        d[keys[-1]] = []
    d[keys[-1]].append(value)


def ReadGroup2DomainMap(path_group2domain_map):
    map_group2domain = {}  # KIN   group   name

    f = open(path_group2domain_map, "rU")
    #	KinomeXplorer_all_predictions_v2.csv
    for line in f.readlines():
        tokens = line.split()
        InsertValueIntoMultiLevelDict(map_group2domain, tokens[:2], tokens[2])

    f.close()

    return map_group2domain


def SetValueIntoMultiLevelDict(d, keys, value):
    for i in range(len(keys) - 1):
        # if not d.has_key(keys[i]):
        if keys[i] not in d:
            d[keys[i]] = {}
        d = d[keys[i]]
    if (keys[-1] in d) and type(d[keys[-1]]) != type(value):
        sys.stderr.write("Causion: multi-dict already has value and try to assign a value of different type")
        pass
    if (keys[-1] in d):
        if d[keys[-1]] != value:
            # sys.stderr.write("This operation replaces a value (%s)" % " ".join(map(lambda(x):str(x), keys)))
            sys.stderr.write("This operation replaces a value ({})".format(' '.join(map(str, keys))))
    d[keys[-1]] = value


def printResult(id_pos_tree_pred, tree_pred_string_data, incoming2string, string_alias, string_desc, string_names,
                organism, mode, dir_likelihood_conversion_tbl, map_group2domain, res_dir, fasta_file):
    ALPHA = ALPHAS[organism]
    species = dSpeciesName[organism]

    fas = open(fasta_file, 'r')
    csv_filename = os.path.join(res_dir, f'{fas.name}.result.csv')
    with open(csv_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow([ 'Name'
                        , 'Position'
                        , 'Tree'
                        , 'NetPhorest Group'
                        , 'Kinase/Phosphatase/Phospho-binding domain'
                        , 'NetworKIN score'
                        , 'NetPhorest probability'
                        , 'STRING score'
                        , 'Target STRING ID'
                        , 'Kinase STRING ID'
                        , 'Target Name'
                        , 'Kinase Name'
                        , 'Target description'
                        , 'Kinase description'
                        , 'Peptide sequence window'
                        , 'Intermediate nodes'])

        dLRConvTbl = {}
        for fname in glob.glob(os.path.join(dir_likelihood_conversion_tbl, "conversion_tbl_*_smooth*")):
            netphorest_or_string, species_of_conversion_table, tree, player_name = \
            re.findall("conversion_tbl_([a-z]+)_smooth_([a-z]+)_([A-Z0-9]+)_([a-zA-Z0-9_/-]+)",
                       os.path.basename(os.path.splitext(fname)[0]))[0]
            # species, tree, player_name = os.path.basename(os.path.splitext(fname)[0]).rsplit('_', 3)[1:]

            if species_of_conversion_table != species:
                continue

            conversion_tbl = ReadConversionTableBin(fname)
            SetValueIntoMultiLevelDict(dLRConvTbl,
                                       [species_of_conversion_table, tree, player_name, netphorest_or_string],
                                       conversion_tbl)

            if options.verbose:
                sys.stderr.write("Conversion table %s %s %s %s\n" % (
                species_of_conversion_table, tree, player_name, netphorest_or_string))

        # For each ID in NetPhorest
        for id in id_pos_tree_pred:
            # We have a mapping to STRING
            if id in incoming2string:
                # For each predicted position
                for pos in id_pos_tree_pred[id]:
                    # For each of the trees (KIN, SH@ etc.)
                    for tree in id_pos_tree_pred[id][pos]:
                        score_results = {}
                        # For each single classifier
                        for pred in id_pos_tree_pred[id][pos][tree]:
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
                                # print(tree_pred_string_data)
                                if string1 in tree_pred_string_data:
                                    (res, peptide, netphorestScore) = id_pos_tree_pred[id][pos][tree][pred]
                                    for string2 in tree_pred_string_data[string1]:

                                        if string2 in string_alias:
                                            bestName2 = string_alias[string2]
                                        else:
                                            bestName2 = ''
                                        if string2 in string_desc:
                                            desc2 = string_desc[string2]
                                        else:
                                            desc2 = ''
                                        stringScore = tree_pred_string_data[string1][string2]["_score"]
                                        path = tree_pred_string_data[string1][string2]["_path"]
                                        name = tree_pred_string_data[string1][string2]["_name"]  # string2 = kinase

                                        # sys.stderr.write("%s %s %s\n" % (tree, pred, name))
                                        if (tree not in map_group2domain) or (pred not in map_group2domain[tree]) or (
                                                name not in map_group2domain[tree][pred]):
                                            if options.string_for_uncovered:
                                                if species == "human":
                                                    if tree in ["1433", "BRCT", "WW", "PTB", "WD40", "FHA"]:
                                                        conversion_tbl_string = dLRConvTbl[species]["SH2"]["general"][
                                                            "string"]
                                                    else:
                                                        conversion_tbl_string = dLRConvTbl[species][tree]["general"][
                                                            "string"]
                                                elif species == "yeast":
                                                    conversion_tbl_string = dLRConvTbl[species][tree]["general"][
                                                        "string"]
                                                else:
                                                    raise "This species is not supported"

                                                likelihood_netphorest = 1
                                                likelihood_string = ConvertScore2L(stringScore, conversion_tbl_string)
                                                unified_likelihood = likelihood_netphorest * likelihood_string
                                                networkinScore = unified_likelihood

                                                # NetworKIN result
                                                result = id + '\t' + res + str(
                                                    pos) + '\t' + tree + '\t' + pred + '\t' + name + '\t' + \
                                                         format(networkinScore, '.4f') + '\t' + "N/A" + '\t' + format(
                                                    stringScore, '.4f') + '\t' + \
                                                         string1 + '\t' + string2 + '\t' + bestName1 + '\t' + bestName2 + '\t' + desc1 + '\t' + desc2 + '\t' + peptide + '\t' + path + '\n'

                                                with open(csv_filename, 'a', newline='') as csvfile:
                                                    writer = csv.writer(csvfile, delimiter=',')
                                                    writer.writerow(result.split('\t'))
                                            else:
                                                continue
                                        else:
                                            # sys.stderr.write("%s %s\n" % (string1, string2))

                                            if species == "human":
                                                if tree in ["1433", "BRCT", "WW", "PTB", "WD40", "FHA"]:
                                                    conversion_tbl_netphorest = dLRConvTbl[species]["SH2"]["general"][
                                                        "netphorest"]
                                                    conversion_tbl_string = dLRConvTbl[species]["SH2"]["general"][
                                                        "string"]
                                                else:
                                                    if (name in dLRConvTbl[species][tree]):
                                                        conversion_tbl_netphorest = dLRConvTbl[species][tree][name][
                                                            "netphorest"]
                                                        conversion_tbl_string = dLRConvTbl[species][tree][name][
                                                            "string"]
                                                    else:
                                                        conversion_tbl_netphorest = \
                                                        dLRConvTbl[species][tree]["general"]["netphorest"]
                                                        conversion_tbl_string = dLRConvTbl[species][tree]["general"][
                                                            "string"]
                                            elif species == "yeast":
                                                if (name in dLRConvTbl[species][tree]):
                                                    conversion_tbl_netphorest = dLRConvTbl[species][tree][name][
                                                        "netphorest"]
                                                    conversion_tbl_string = dLRConvTbl[species][tree][name]["string"]
                                                else:
                                                    conversion_tbl_netphorest = dLRConvTbl[species][tree]["general"][
                                                        "netphorest"]
                                                    conversion_tbl_string = dLRConvTbl[species][tree]["general"][
                                                        "string"]
                                            else:
                                                raise "This species is not supported"

                                            likelihood_netphorest = ConvertScore2L(netphorestScore,
                                                                                   conversion_tbl_netphorest)
                                            likelihood_string = ConvertScore2L(stringScore, conversion_tbl_string)
                                            unified_likelihood = likelihood_netphorest * likelihood_string
                                            networkinScore = unified_likelihood
                                            # networkinScore = pow(stringScore, ALPHA)*pow(netphorestScore, 1-ALPHA)

                                            # NetworKIN result

                                            try:
                                                result = string_names[organism + '.' + id] + '\t' + res + str(
                                                    pos) + '\t' + tree + '\t' + pred + '\t' + name + '\t' + \
                                                         format(networkinScore, '.4f') + '\t' + format(netphorestScore,
                                                                                                       '.4f') + '\t' + format(
                                                    stringScore, '.4f') + '\t' + \
                                                         string1 + '\t' + string2 + '\t' + bestName1 + '\t' + bestName2 + '\t' + desc1 + '\t' + desc2 + '\t' + peptide + '\t' + path + '\n'
                                            except:
                                                result = id + '\t' + res + str(
                                                    pos) + '\t' + tree + '\t' + pred + '\t' + name + '\t' + \
                                                         format(networkinScore, '.4f') + '\t' + format(netphorestScore,
                                                                                                       '.4f') + '\t' + format(
                                                    stringScore, '.4f') + '\t' + \
                                                         string1 + '\t' + string2 + '\t' + bestName1 + '\t' + bestName2 + '\t' + desc1 + '\t' + desc2 + '\t' + peptide + '\t' + path + '\n'

                                            with open(csv_filename, 'a', newline='') as csvfile:
                                                writer = csv.writer(csvfile, delimiter=',')
                                                writer.writerow(result.split('\t'))

                                        if networkinScore not in score_results:
                                            score_results[networkinScore] = []

                                        score_results[networkinScore].append(result)
                        if mode == 'network':
                            highestScore = sorted(score_results.keys(), reverse=True)[0]

                            if len(score_results[highestScore]) > 1:
                                index = random.randint(0, len(score_results[highestScore]) - 1)
                            # sys.stdout.write((score_results[highestScore][index]))
                            else:
                                # sys.stdout.write((score_results[highestScore][0]))
                                pass
                        else:
                            for score in sorted(score_results.keys(), reverse=True):
                                # sys.stdout.write("".join(score_results[score]))
                                pass
            else:
                pass
    return


# MAIN
def Main():
    res_dir = 'results'
    if not os.path.exists(res_dir):
        os.makedirs(res_dir)

    tmpdir = 'tmp'
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)

    sys.stderr.write("Reading fasta input file\n")
    id_seq = readFasta(fastafile)
    if options.verbose:
        sys.stderr.write("%s sequences loaded\n" % len(id_seq.keys()))

    if sitesfile:
        sys.stderr.write("Reading phosphosite file\n")
        input_type = CheckInputType(sitesfile)
        if input_type == NETWORKIN_SITE_FILE:
            id_pos_res = readPhosphoSites(sitesfile)
        elif input_type == PROTEOME_DISCOVERER_SITE_FILE:
            id_pos_res = readPhosphoSitesProteomeDiscoverer(fn_fasta, sitesfile)
        elif input_type == MAX_QUANT_DIRECT_OUTPUT_FILE:
            id_pos_res = readPhosphoSitesMaxQuant(sitesfile)
        elif input_type == RUNES_SITE_FILE:
            id_pos_res = readRunessitesfile(sitesfile)
    else:
        id_pos_res = {}

    sys.stderr.write("Loading aliases and descriptions\n")
    (string_alias, string_desc, string_names) = readAliasFiles(args[0], options.datadir);

    if organism == "9606":
        path_group2domain_map = os.path.join(options.datadir, "group_human_protein_name_map.tsv")
    elif organism == "4932":
        path_group2domain_map = os.path.join(options.datadir, "group_yeast_KIN.tsv")

    map_group2domain = ReadGroup2DomainMap(path_group2domain_map)

    # Default way of mapping using BLAST
    # incoming2string, string2incoming = mapPeptides2STRING(blastDir, organism, fastafile.name, id_pos_res, id_seq, options.threads, options.datadir, options.fast, options.leave)
    incoming2string, string2incoming = mapPeptides2STRING(blastDir, organism, fastafile.name, id_pos_res, id_seq,
                                                          options.threads, options.datadir)

    # Hack for random mapping to proteins
    # incoming2string, string2incoming = mapRandom(id_seq)

    # Use if a mapping file for the input can be provided
    # incoming2string, string2incoming = mapFromFile("/home/red1/hhorn/projects2/2012_03_22_Jesper/ensmusp_ensp.tsv")

    # Used if only ensembl of the same version used
    # incoming2string, string2incoming = mapOne2one(id_seq)

    # Load the STRING network data
    sys.stderr.write("Loading STRING network\n")
    tree_pred_string_data = loadSTRINGdata(string2incoming, options.datadir, options.threads)
    # Run NetPhorest
    sys.stderr.write("Running NetPhorest\n")
    # netphorestTmpFiles = runNetPhorest(id_seq, id_pos_res, options.compress, options.threads, options.active_threads, options.fast, options.leave)
    netphorestTmpFiles = runNetPhorest(id_seq, id_pos_res, options.compress, options.threads)
    # sys.stderr.write('\n')

    # Writing result to STDOUT
    sys.stderr.write("Writing results\n")
    sys.stdout.write(
        "#Name\tPosition\tTree\tNetPhorest Group\tKinase/Phosphatase/Phospho-binding domain\tNetworKIN score\tNetPhorest probability\tSTRING score\tTarget STRING ID\tKinase/Phosphatase/Phospho-binding domain STRING ID\tTarget description\tKinase/Phosphatase/Phospho-binding domain description\tTarget Name\tKinase/Phosphatase/Phospho-binding domain Name\tPeptide sequence window\tIntermediate nodes\n")
    for i in range(len(netphorestTmpFiles)):
        id_pos_tree_pred = parseNetphorestFile(netphorestTmpFiles[i].name, id_pos_res, options.compress)
        if options.path == "direct":
            dir_likelihood_conversion_tbl = os.path.join(options.datadir, "likelihood_conversion_table_direct")
        elif options.path == "indirect":
            dir_likelihood_conversion_tbl = os.path.join(options.datadir, "likelihood_conversion_table_indirect")
        else:
            raise "Path information should be either direct or indirect."
        # print(tree_pred_string_data)
        printResult(id_pos_tree_pred, tree_pred_string_data, incoming2string, string_alias, string_desc, string_names,
                    args[0], options.mode, dir_likelihood_conversion_tbl, map_group2domain, res_dir, fn_fasta)

    return


if __name__ == '__main__':
    # BLAST
    try:
        blastDir = os.environ['BLAST_PATH']
    except:
        blastDir = ""
    # NETPHOREST
    try:
        netphorest_bin = os.environ['NETPHOREST_PATH']
    except:
        netphorest_bin = ""

    usage = "usage: %prog [options] organism FASTA-file [sites-file]"
    parser = OptionParser(usage=usage, version="%prog 3.0")
    parser.add_option("-n", "--netphorest", dest="netphorest_bin", default=netphorest_bin,
                      help="set the location of the NetPhorest binary, overwrites the 'NETPHOREST_PATH' environmental variable. [ENV: %default]")
    parser.add_option("-b", "--blast", dest="blast", default=blastDir,
                      help="set the directory for the BLAST binaries (formatdb and blastall), overwrites the 'BLAST_PATH' environmental variable. [ENV: %default]")
    parser.add_option("-m", "--mode", dest="mode", default=False,
                      help="if set to 'network', gives only one best scoring result for each site. In case of multiple candidate kinases with the same core, the selection hapens randomly. [default: %default]")
    parser.add_option("-p", "--path", dest="path", default="direct",
                      help="NetworKIN uses both direct and indirect paths. Otherwise, it uses only indirect paths. [default: %default]")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
                      help="print out everything [default: %default]")
    parser.add_option("-f", "--fast", dest="fast", default=False, action="store_true",
                      help="Speed up by using the intermediate files of previous run [default: %default]")
    parser.add_option("-l", "--leave", dest="leave", default=False, action="store_true",
                      help="leave intermediate files [default: %default]")
    parser.add_option("-u", "--uncovered", dest="string_for_uncovered", default=False, action="store_true",
                      help="Use STRING likelihood for uncovered Kinases [default: %default]")
    parser.add_option("-t", "--threads", dest="threads", default=1, type="int",
                      help="number of available threads/CPUs. Also leads to less memory usage, as result files are read sequentially [default: %default]")
    parser.add_option("--nt", dest="active_threads", default=2, type="int",
                      help="number of active threads at a time")

    parser.add_option("-c", "--compress", dest="compress", default=True,
                      help="compress temporary result files, saves discspace [default: %default]")
    parser.add_option("-d", "--data", dest="datadir",
                      default=os.path.join(os.path.split(os.path.realpath(sys.argv[0]))[0], 'data'),
                      help="location for the additional files like the pre-computed STRING network, STRING sequence database etc. [default: %default]")
    # if os.environ.has_key("TMPDIR"):
    #	parser.add_option("--tmp", dest="tmpdir", default=os.environ["TMPDIR"],
    #									help="location for the temporary files [default: %default]")
    # else:
    #	print >> sys.stderr, "TMPDIR environmental variable is not defined. Please define this variable or specify tmp directory by using --tmp command line option"
    #	sys.exit()

    global options
    (options, args) = parser.parse_args()

    if options.active_threads < 2:
        parser.error("Number of active thread (--nt) is less than 2")
    # tempfile.tempdir= options.tmpdir

    # ORGANISM
    try:
        organism = args[0]
    except:
        parser.error("Organism not defined!")

    # SEQUENCE FILE
    try:
        fn_fasta = args[1]
        fn_blast_output = "%s.%s.blast.out" % (fn_fasta, organism)
        fn_netphorest_output = "%s.%s.netphorest.out" % (fn_fasta, organism)
        fastafile = open(fn_fasta, 'rU')
    except:
        sys.stderr.write("%s" % args)
        parser.error("FASTA-file not defined!")

    # SITES FILE
    try:
        sitesfile = args[2]
    except:
        sitesfile = False

    # BLAST
    if options.blast:
        blastDir = options.blast

    # NETPHOREST
    if options.netphorest_bin:
        netphorest_bin = options.netphorest_bin

    # Show runtime parameters
    if (options.verbose):
        sys.stderr.write(
            '\nPredicting using parameters as follows:\nOrganism:\t%s\nFastaFile:\t%s\n' % (organism, fn_fasta))
        sys.stderr.write('Threads:\t%s\nCompress:\t%s\n' % (options.threads, options.compress))
        if options.string_for_uncovered:
            sys.stderr.write("Use STRING likelihood when the kinases is not covered NetPhorest.\n")
        if sitesfile:
            sys.stderr.write('Sitesfile:\t%s\n' % sitesfile)
        else:
            sys.stderr.write('No sites-file given, predicting on all S,T,Y residues.\n')
        if options.mode:
            sys.stderr.write('Mode:\t\t%s\n' % options.mode)
        sys.stderr.write("Blast dir: %s" % blastDir)
        sys.stderr.write("NetPhorest binary: %s" % netphorest_bin)
        sys.stderr.write('\n')
    Main()
