import pandas as pd

df = pd.read_csv('MS_Gaussian_updated_09032023.csv')
df = df.dropna(subset=['site'])
# Create pairs of adjacent values from the 'ENSP' and 'site' columns
pairs = list(zip(df['ENSP'], df['site']))
unique_pairs = list(set(pairs))
unique_pairs_df = pd.DataFrame(unique_pairs, columns=['ENSP', 'site'])

dft=pd.read_csv('data/string_data/9606.protein.info.v12.0.txt', delimiter='\t')
dft['#string_protein_id'] = dft['#string_protein_id'].str.slice(start=5)

ENSPS=[]
fastafile = open('cured_morpho_seqs_v2.fa','rU')
data = fastafile.readlines()
for line in data:
    if line[0] == '>':
        ENSPS.append(line[1:16])

with open('MS_Gaussian_updated_09032023.tsv', 'w') as f:
    for e in ENSPS:
        #print(e)
        if e in unique_pairs_df['ENSP'].values:
            #print(unique_pairs_df[unique_pairs_df['ENSP']==e])
            matching_rows = unique_pairs_df[unique_pairs_df['ENSP'] == e]
            for _, row in matching_rows.iterrows():
                first_column = row['ENSP']
                second_column = row['site']
                formatted_row = f"{first_column}\t{second_column[2:]}\t{second_column[0]}"
                # Write the formatted row to the file
                if second_column[0] in ['S','T','Y']:
                    f.write(formatted_row + '\n')
    