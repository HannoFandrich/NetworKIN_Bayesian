import sys, os, subprocess, re, tempfile, random, operator
import platform
from itertools import chain


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

def readAliasFiles(organism, datadir):
    name_hash = {}
    try:
        # desc_db = myPopen('gzip -cd %s/%s.text_best.v9.0.tsv.gz'%(datadir, organism))
        name_db = myPopen('gzip -cd %s/%s.protein.aliases.v12.0.txt.gz' % (datadir, organism))
        name_db = name_db.split('\n')
        name_hash={}
        for line in name_db:
            # print(line)
            (seqID, alias, source) = line.split('\t')
            key=seqID[5:]
            if key in name_hash.keys():
                if alias in name_hash[key]:
                    continue
                else:
                    name_hash[key].append(alias)
            else:
                name_hash[key]=[alias]
    except:
        sys.stderr.write("No names available for organism: '%s'\n" % organism)

    return  name_hash

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
    f = open(path_group2domain_map, "r")
    names=list(name_hash.values())
    names=list(chain(list(chain(*names[1:]))))
    c=0
    cl=0
    not_mapped=[]
    with open("hanno_group_human_protein_name_map.tsv", "w") as output_file:
        for line in f.readlines():
            cl+=1
            tokens = line.split()
            name = tokens[2]
            name_og = tokens[2]
            if name in names:
                c+=1
                new_line = f"{tokens[0]}\t{tokens[1]}\t{name_og}\t{name}\n"
                output_file.write(new_line)
            else:
                if len(name)==3:
                    name = name.upper()
                if name[:4]=='1433':
                    name=name.replace('1433','14-3-3')
                    '''
                if name[:5]=='PIK3R':
                    name=name[:6] # don't know what these subdivisions are
                if name[:4]=='PTPN' and '_' in name:
                    name=name[:-2] # again don't know what these subdivisions are
                    '''
                if '_' in name:
                    name=name[:-2]
                name = name.replace('CK1', 'CSNK1')
                name = name.replace('CaMKII', 'CaMK2')
                name = name.replace('CTDSPL1', 'CTDSPL')
                name = name.replace('CTPTP', 'CTDSP')  # sacco et al.
                name = name.replace('DMPK1', 'DMPK')
                name=name.replace('HH498','TNNI3K') # https://www.guidetopharmacology.org/GRAC/FamilyDisplayForward?familyId=521
                name=name.replace('HSER','GUCY2C') # http://kinase.com/web/current/kinbase/gene/5631
                name=name.replace('JAKb','JAK') # realy dont know what JAKb is supposed to be
                name=name.replace('KHS2', 'MAP4K3') # http://kinase.com/web/current/kinbase/gene/5656
                name=name.replace('LOC440917','14-3-3epsilon') # discontinued https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=DetailsSearch&Term=440817
                name=name.replace('MDSP','DUSP13') # https://www.uniprot.org/uniprotkb/Q6B8I1/entry
                name = name.replace('p70S6K', 'PS6K')  # https://www.uniprot.org/uniprotkb/Q94533/entry
                name=name.replace('PKG1cGKI','PRKG1') # Compare bestpathv9.0, aliasesv9.0 and aliasesv12.0
                name = name.replace('PKG2cGKII', 'PRKG2') # Compare bestpathv9.0, aliasesv9.0 and aliasesv12.0
                name=name.replace('RSK5','RSKL') # https://dlb.idrblab.net/ttd/data/target-ppi/details/t50224
                #name=name.replace('SACM1L','SAC1')
                name=name.replace('SgK424','TEX14') # http://www.kinhub.org/kinases.html
                name = name.replace('STLK3','STK39') # https://www.phosphosite.org/proteinAction.action?id=758&showAllSites=true
                name = name.replace('STLK5','LYK5') # http://www.phosphonet.ca/?search=Q7RTN6
                name = name.replace('STLK6','ALS2CR2') # http://www.phosphonet.ca/kinasepredictor.aspx?uni=Q8N7Q3&ps=T136
                name = name.replace('TAO3','TAOK3') # https://www.uniprot.org/uniprotkb/Q9H2K8/entry
                name=name.replace('TRCPbeta','BTRC') # https://en.wikipedia.org/wiki/BTRC_(gene)



                name = name.replace('Domain2_', '')
                name = name.replace('Domain2','')
                name = name.replace('Aurora', 'Aurora kinase ')
                name = name.replace('eta', ' eta')
                name = name.replace('alpha',' alpha')
                name = name.replace('b eta', ' beta')
                name = name.replace('gamma', ' gamma')
                name = name.replace('delta', ' delta')
                name = name.replace('epsilon', ' epsilon')
                name = name.replace('z eta', ' zeta')
                name = name.replace('th eta', ' theta')
                name = name.replace('iota', ' iota')
                name = name.replace('kappa', ' kappa')
                name = name.replace('sigma', ' sigma')
                if name in names:
                    c += 1
                    new_line = f"{tokens[0]}\t{tokens[1]}\t{name_og}\t{name}\n"
                    output_file.write(new_line)
                else:
                    #name = name.replace('CK1','CKI')

                    name = name.replace('_',' ')
                    name = name.replace(' alpha', '-alpha')
                    name = name.replace(' beta', '-beta')
                    name = name.replace(' gamma', '-gamma')
                    name = name.replace(' delta', '-delta')
                    name = name.replace(' epsilon', '-epsilon')
                    name = name.replace(' zeta', '-zeta')
                    name = name.replace(' eta', '-eta')
                    name = name.replace(' theta', '-theta')
                    name = name.replace(' iota', '-iota')
                    name = name.replace(' kappa', '-kappa')
                    name = name.replace(' sigma', '-sigma')
                    name=name.replace('PKC-theta','nPKC-theta') # https://www.uniprot.org/uniprotkb/Q04759/entry
                    name = name.replace('PKC-zeta', 'nPKC-zeta')  # https://www.uniprot.org/uniprotkb/O19111/entry
                    if name in names:
                        c += 1
                        new_line = f"{tokens[0]}\t{tokens[1]}\t{name_og}\t{name}\n"
                        output_file.write(new_line)
                    else:
                        name = name.replace('_', ' ')
                        name = name.replace('-alpha', ' C-alpha')
                        name = name.replace('C-alpha2', 'C-alpha 2')
                        name = name.replace('-beta', ' C-beta')
                        name = name.replace('-gamma', ' C-gamma')
                        name = name.replace('C-gamma2', 'C-gamma 2')
                        name = name.replace('C-gamma3', 'C-gamma 3')
                        name = name.replace('-delta', ' C-delta')
                        name = name.replace('-epsilon', ' C-epsilon')
                        name = name.replace('-zeta', ' C-zeta')
                        name = name.replace('-eta', ' C-eta')
                        name = name.replace('-theta', ' C-theta')
                        name = name.replace('-iota', ' C-iota')
                        name = name.replace('-kappa', ' C-kappa')
                        name = name.replace('AlphaK','Alpha kinase ')
                        name = name.replace('ACTR','ACVR')
                        name = name.replace('AUXI', 'AUXI_HUMAN')
                        if name in names:
                            c += 1
                            new_line = f"{tokens[0]}\t{tokens[1]}\t{name_og}\t{name}\n"
                            output_file.write(new_line)
                        else:
                            name =name.replace('K ',' kinase ')
                            if name in names:
                                c += 1
                            else:

                                name = name.replace('II', '2')
                                name = name.replace('III', '3')

                                name = name.replace(' C-alpha', 'A')
                                name = name.replace(' C-beta', 'B')
                                name = name.replace(' C-gamma', 'G')
                                name = name.replace(' C-delta', 'D')
                                name = name.replace(' C-epsilon', 'E')
                                name = name.replace(' C-zeta', 'Z')
                                name = name.replace(' C-eta', 'L')
                                name = name.replace(' C-theta', 'Q')
                                name = name.replace(' C-iota', 'I')
                                '''
                                name = name.replace(' C-kappa', 'J')
                                '''
                                name=name.replace('A 2','A2')
                                name=name.replace('G 2','G2')
                                name = name.replace('G 3', 'G3')

                                name = name.replace('CSNK1A2', 'CSNK1A1L')  # https://www.phosphosite.org/proteinAction?id=3753&showAllSites=true

                                name = name.upper()
                                if name in names:
                                    c += 1
                                    new_line = f"{tokens[0]}\t{tokens[1]}\t{name_og}\t{name}\n"
                                    output_file.write(new_line)
                                else:
                                    new_line = f"{tokens[0]}\t{tokens[1]}\t{name_og}\t{name_og}\n"
                                    name=name_og
                                    not_mapped.append(name)
                                    output_file.write(new_line)
            InsertValueIntoMultiLevelDict(map_group2domain, tokens[:2], name)
        output_file.close()
    f.close()
    print(not_mapped)
    return map_group2domain
name_hash = readAliasFiles('9606','.')
map_group2domain = ReadGroup2DomainMap("string_data/group_human_protein_name_map.tsv")