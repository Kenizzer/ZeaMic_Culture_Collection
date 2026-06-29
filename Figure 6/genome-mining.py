## Code by: Anna Garrell 2026

import os
import pandas as pd

os.chdir('./USDA-genomemining')

# Replace 'tab_files_folder' with the folder where your .tab files are located
tab_files_folder = 'USDA-genomemining/tabfile'
# Replace 'summary_folder' with the folder where you want to save the summary files for each individual tab file
summary_folder = 'USDA-genomemining/results'
# Replace 'txt_files_folder' with the folder where you want to save the text files
txt_files_folder= 'USDA-genomemining/annotatedgenomestxtfiles'
# Replace 'all_genes_folder' with the folder where you want to save the combined summary file
all_genes_folder = 'USDA-genomemining/results'


# Get a list of all .tab files in the folder
tab_files = [f for f in os.listdir(tab_files_folder) if f.endswith('.tab')]

# Make an empty list to save the summaries of the individual tab files
summary_dfs = [[]]*len(tab_files)
i = 0 # this index is for saving the individual tab files

for tab_file in tab_files:
    
    # Read the tab file
    tab_df = pd.read_csv(os.path.join(tab_files_folder, tab_file), sep='\t')

    # Create a DataFrame to store the summary
    summary_df = pd.DataFrame()
    
    summary_df['Function'] = tab_df.iloc[:, 1] # Assuming the KO column is the second column (index 1)
    
    summary_df['KO'] = tab_df.iloc[:, 2]  # Assuming the KO column is the third column (index 2)
        # Add in a column that specifies the tab file where KO is found
            # use df.insert to put the new column as the first column
    summary_df.insert(0, 'Original_Tab_File', tab_file.split('.')[0])
            
    txt_files = [f for f in os.listdir(txt_files_folder) if f.endswith('.txt')]

    for txt_file in txt_files:
        # Read the text file
        txt_df = pd.read_csv(os.path.join(txt_files_folder, txt_file),
                             sep='\t', header=None, skiprows = 3, skipfooter=3, 
                             engine = 'python')
        
        txt_df2 = pd.read_csv(os.path.join(txt_files_folder, txt_file), sep='\t', header=None, skiprows=5)

        # Get the header (file name) without the extension
        header = os.path.splitext(txt_file)[0]

        # Initialize a list to store the results
        results = []

        # Iterate through the rows of the tab DataFrame
        for index, row in tab_df.iterrows():
            # Extract the search term from the row
            search_term = row['KO']

            # Search for the search term in the 12th column of the text DataFrame
            found = txt_df.iloc[:, 11].str.contains(search_term, case=False).any()

            # Append 1 if found, 0 if not found
            results.append(1 if found else 0)

        # Add the results to the summary DataFrame with the corresponding header
        summary_df[header] = results

    # Save the summary dataframe of each tab file to a CSV file
    summary_name = os.path.splitext(tab_file)[0] + '_summary.csv'
    summary_path = os.path.join(summary_folder, summary_name)
    summary_df.to_csv(summary_path, index=False)
    
    # Save the summary of each tab file in a list to create the full summary file
    summary_dfs[i] = summary_df
    i += 1
    

# Combine all of the summary files into one file
all_genes = pd.concat(summary_dfs)

    # save the combined summary
all_genes.to_csv(all_genes_folder + '/all_genes_summary.csv', index=False)

# in the tab file; add in a column for the function of the KO in the same tab file as the original KO
    # allways the second column always named function