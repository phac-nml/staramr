import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
import sys

# This file updates the resfinder and pointfinder tsv files to the proper format
# To run the program, use the command python update_staramr_db.py <file> <resfinder,pointfinder>
# Requires pandas, sys and xlrd>=1.0.0 as dependencies, create a new conda environment for this script

def get_accession(value):
    split_value = value.split('_')
    accession_value = ""
    if(len(split_value) == 2):
        # Means a gene does not have a variant number, needs to be added in
        accession_value = split_value[1]
    elif(len(split_value) == 3):
        accession_value = split_value[2]
    else:
        accession_value = "_".join(split_value[2:4])
    return accession_value

def get_gene_name(value):
    split_value = value.split('_')
    gene_value = "_".join(split_value[0:2])
    
    return gene_value

def format_drug_name(value):
    # Convert string into lowercase first
    value = value.lower()
    
    # There's 2 formatting cases we have to deal with
    # some values in the drug column had a combination of commas and spaces
    # need to format so we seperate values by commas with no spaces
    if(" " in value and "," in value):
        value = value.replace(" ", "")
    elif(" " in value):
        value = value.replace(" ", ",")

    # Check if ciprofloxacinI/R is in there
    if("ciprofloxacini/r" in value):
        value = value.replace("i/r", " I/R")

    # Check if the word acid is in there
    index = value.find('acid')
    
    # Replace the comma with a space ex clavulanic,acid to clavulanic acid
    if(index > 0):
        index = index-1
        value = value[:index] + " " + value[index + 1:]
    
    split_value = value.split(',')
    # If there's more than one drug type
    if(len(split_value) > 1):
        drug_value = ",".join(split_value)
    else:
        drug_value = value
        
    return drug_value

def get_point_gene_name(value):
    if("ampCprom" in value):
        gene_value = 'ampC_promoter_size_53bp'
    else:
        value = value.replace("(", "=")
        split_value = value.split('=')
        gene_value = split_value[0]
    
    return gene_value

def get_codon_pos(value):
    value = value.replace("(", "=")
    value = value.replace(")", "")
    value = value.replace("-", "=-")
    split_value = value.split('=')
    codon_value = split_value[1]
    
    return codon_value

def format_resfinder(filename):
  columns = ['Class', 'Gene' ,'Accession', 'Drug']
  final_frame = pd.DataFrame(columns=columns)
  df = pd.ExcelFile(filename)
  sheet_names = df.sheet_names

  dfs = pd.read_excel(filename, sheet_name=None)

  for sheet in sheet_names:
    dfs = pd.read_excel(filename, sheet_name=sheet, header=None)
    class_name = sheet.lower()

    if(class_name == 'b-lactam'):
      class_name = 'beta-lactam'

    if(len(dfs.columns) >= 3):
      dfs = dfs[dfs.columns[0:2]]
      
    # ResFinder needs the following columns: [Class,Gene,Accession,Drug]
    dfs.columns = ['Gene','Drug']
    
    # Parse the Accession string into it's own column
    dfs['Accession'] = dfs.apply(lambda x: get_accession(x['Gene']),axis=1)
    dfs['Gene'] = dfs.apply(lambda x: get_gene_name(x['Gene']),axis=1)
    dfs['Class'] = class_name
    dfs = dfs.reindex(columns=columns)

    # If there are empty values in the drug column, replace it with None for now
    dfs['Drug'] = dfs['Drug'].where((pd.notnull(dfs['Drug'])), 'None')
    dfs['Drug'] = dfs.apply(lambda x: format_drug_name(x['Drug']),axis=1)

    final_frame = final_frame.append(dfs)
    
  # Remove the rows that contains none in them
  none_rows = final_frame.loc[final_frame['Drug'] == 'none']

  if not none_rows.empty:
    print("The following genes have been removed for no drug value:\n")
    print(none_rows[columns[0:3]].to_string(index=False))

    final_frame = final_frame[final_frame['Drug'].map(lambda x: str(x)!="none")]

  return final_frame

def format_pointfinder(filename):
  columns = ['Organism', 'Gene', 'Codon Pos.', 'Drug']
  final_frame = pd.DataFrame(columns=columns)
  df = pd.ExcelFile(filename)
  sheet_names = df.sheet_names

  dfs = pd.read_excel(filename, sheet_name=None)

  for sheet in sheet_names:
    dfs = pd.read_excel(filename, sheet_name=sheet, header=None)
    class_name = sheet.lower()
    
    # PointFinder needs the following: [Organism, Gene, Codon Pos. ,Drug]
    dfs.columns = ['Gene','Drug']
    
    dfs['Codon Pos.'] = dfs.apply(lambda x: get_codon_pos(x['Gene']),axis=1)
    dfs['Gene'] = dfs.apply(lambda x: get_point_gene_name(x['Gene']),axis=1)
    dfs['Organism'] = class_name
    dfs = dfs.reindex(columns=columns)

     # If there are empty values in the drug column, replace it with None for now
    dfs['Drug'] = dfs['Drug'].where((pd.notnull(dfs['Drug'])), 'None')
    dfs['Drug'] = dfs.apply(lambda x: format_drug_name(x['Drug']),axis=1)

    final_frame = final_frame.append(dfs)

  # Remove the rows that contains none in them
  none_rows = final_frame.loc[final_frame['Drug'] == 'none']

  if not none_rows.empty:
    print("The following genes have been removed for no drug value:\n")
    print(none_rows[columns[0:3]].to_string(index=False))

    final_frame = final_frame[final_frame['Drug'].map(lambda x: str(x)!="none")]

  return final_frame

# Dynamically read in a file
command = sys.argv

# Exit the program if the arguments is not 2
if(len(command) <= 2):
  print("Error:  <file name> <resfinder, pointfinder>")
  exit()
else:
  # If the filename has spaces in them, use [1:] to get a 
  # list of the rest of the arguments and turn it into a string
  filename = command[1]
  database = command[2].lower()

  if(database != 'resfinder' and database != 'pointfinder'):
    print("Error:  <file name> <resfinder, pointfinder>")
    exit()
  else:

    if database == 'resfinder':
      final_frame = format_resfinder(filename)
      file_name = 'ARG_drug_key_resfinder.tsv'
    else:
      final_frame = format_pointfinder(filename)
      file_name = 'ARG_drug_key_pointfinder.tsv'
    
    # Convert the file to a tsv file
    print('\n{} has been created.\n'.format(file_name))
    final_frame.to_csv(file_name, sep='\t', encoding='utf-8', index=False)