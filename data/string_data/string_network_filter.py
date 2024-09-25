import gzip
import csv
import shutil
import os
### Create .links file with in the same format as old bestpath.v9
### needs to be executed everytime new data from string is added

# Open the new protein links file
with gzip.open('9606.protein.links.v12.0.txt.gz', 'rb') as f:
    content = f.read().decode('utf-8')
lines = content.split('\n')
interaction_list=[]
for line in lines:
    # Split each line by tab ('\t') to get columns
    columns = line.split(' ')
    interaction_list.append(columns)

print(len(interaction_list))


# Get Protein names from .info file
name_hash={}
with open('9606.protein.info.v12.0.txt', 'r') as file:
    content = file.read()#.decode('utf-8')
    name_db = content.split('\n')
    name_db.pop(-1)
    for line in name_db:
        #print(line)
        (string_protein_id,preferred_name,protein_size,annotation) = line.split('\t')
        name_hash[string_protein_id] = preferred_name


# Get the corresponding Group for every Kinase from the group.tsv file
data_dict = {}
with open('group_human_protein_name_map.tsv', 'r') as file:
    for line in file:
        try:
            columns = line.strip().split('\t')
            KIN= columns[0]
            group= columns[1]
            kinase = columns[2]
            data_dict[kinase] = [KIN, group]
        except:
            continue


# write it all in a new .links file
file_name = '9606.links.v12.0.tsv'
with open(file_name, 'w', newline='') as file:
    writer = csv.writer(file, delimiter='\t')
    c=0
    for i in interaction_list:
        try:
            kinase=name_hash[i[0]]
            if kinase in data_dict:
                KIN=data_dict[kinase][0]
                group=data_dict[kinase][1]
                writer.writerow([KIN,group,kinase,i[0],i[1],i[2]])
                c+=1
        except:
            continue

# Gzip the file
with open(file_name, 'rb') as f_in:
    with gzip.open(file_name + '.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

# Delete the original unzipped .tsv file
os.remove(file_name)
