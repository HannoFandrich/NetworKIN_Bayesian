import pandas as pd

# Read the CSV file
csv="cured_morpho_seqs_v2.fa.result.csv"
csv = pd.read_csv(csv)
print(csv)
print(csv.loc[(csv['Target Name']=='MAP2K2') & (csv['Kinase Name']=='MAPK3')])
#print(csv.loc[(csv['Target Name']=='AKT1') & (csv['Kinase Name']=='MAPK3')])