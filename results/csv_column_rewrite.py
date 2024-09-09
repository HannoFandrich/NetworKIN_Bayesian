import pandas as pd

# Read the CSV file
old_csv="cured_morpho_seqs_v2.fa.result.csv"
<<<<<<< HEAD
new_csv='KinomeXplorer_all_predictions_v3.csv'
=======
new_csv='KinomeXplorer_all_predictions_v2.csv'
>>>>>>> origin/Hanno
df = pd.read_csv(old_csv)

new_df = df[['Target STRING ID',
             'Position',
             'Name',
             'NetworKIN score',
             'Tree',
             'NetPhorest Group',
             'NetPhorest probability',
             'Kinase/Phosphatase/Phospho-binding domain STRING ID',
             'STRING score',
             'Target Name',
             'Peptide sequence window',
             'Intermediate nodes']]

new_df.rename(columns={'Target STRING ID': '#substrate'}, inplace=True)
new_df.rename(columns={'Position': 'position'}, inplace=True)
new_df.rename(columns={'Name': 'id'}, inplace=True)
new_df.rename(columns={'NetworKIN score': 'networkin_score'}, inplace=True)
new_df.rename(columns={'Tree': 'tree'}, inplace=True)
new_df.rename(columns={'NetPhorest Group': 'netphorest_group'}, inplace=True)
new_df.rename(columns={'NetPhorest probability': 'netphorest_score'}, inplace=True)
new_df.rename(columns={'Kinase/Phosphatase/Phospho-binding domain STRING ID': 'string_identifier'}, inplace=True)
new_df.rename(columns={'STRING score': 'string_score'}, inplace=True)
new_df.rename(columns={'Target Name': 'substrate_name'}, inplace=True)
new_df.rename(columns={'Peptide sequence window': 'sequence'}, inplace=True)
new_df.rename(columns={'Intermediate nodes': 'string_path'}, inplace=True)
new_df = new_df[new_df['networkin_score'] >= 0.00001]
new_df['Iteration'] = 0
print(df.head)
print(new_df.head)

new_df.to_csv(new_csv, index=False)
