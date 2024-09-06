import os
import argparse
import subprocess
import fileinput
import sys
import pandas as pd
import shutil

# Argument parser setup
parser = argparse.ArgumentParser(description='Process populations file and create separate files for each population.')
parser.add_argument("-wd", "--workingdir", required=True, help='The working directory where the populations file is located and where the output files will be created.')
parser.add_argument("-pop", "--popfile", required=True, help='populations tsv file')
parser.add_argument("-vcf", "--variants", required=True, help='vcf file used for population statistics, must match IDs in populations file')
parser.add_argument("-tjd", "--tajimawind", required=True, help='integer representing bp windows for Tajimas D calculation')
parser.add_argument("-maf", "--maf", required=True, help='Minor allele frequency filter (i.e. .05)')
parser.add_argument("-thin", "--thin", required=True, help='Filter SNPs in vcf file by distance to account for linked (i.e. 200 for 200bp)')
args = parser.parse_args()

#Define function to get column averages
def append_column_averages(csv_file):

    # Determine the delimiter based on the file extension
    _, file_extension = os.path.splitext(csv_file)

    if file_extension.lower() in ('.tsv', '.txt', 'Tajima.D'):
        delimiter = '\t'
    else:
        delimiter = ','  # Default to comma for CSV and other cases
    
    # Load the data
    df = pd.read_csv(csv_file, sep=delimiter)
    
    # Convert all columns to numeric, errors='ignore' will leave non-numeric columns unchanged
    df = df.apply(pd.to_numeric, errors='ignore')
    
    # Calculate the mean of each numeric column
    averages = df.mean()
    
    # Set the first entry of averages to "Averages"
    if not averages.empty:
        first_column = df.columns[0] if df.columns.size > 0 else None
        if first_column:
            averages[first_column] = "Averages"
    
    # Append the averages as a new row at the bottom of the DataFrame
    df = df.append(averages, ignore_index=True)
    
    # Save the updated DataFrame back to the same CSV file
    df.to_csv(csv_file, index=False, sep=delimiter)


#Define function to replaces lines
def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)

#function to replace characters
def replace_in_file(file_path, old_string, new_string):
    # Read the contents of the file
    with open(file_path, 'r') as file:
        file_contents = file.read()

    # Replace the old string with the new string
    updated_contents = file_contents.replace(old_string, new_string)

    # Write the updated contents back to the file
    with open(file_path, 'w') as file:
        file.write(updated_contents)

# Change to the working directory
os.chdir(args.workingdir)

# Initialize a dictionary to store population data
population_dict = {}

# Parse population file
with open(os.path.join(args.workingdir, args.popfile), 'r') as populations:
    for line in populations:
        if line.strip():  # Skip empty lines
            id, pop = line.strip().split('\t')
            if pop not in population_dict:
                population_dict[pop] = []
            population_dict[pop].append(id)

#Write the output files for each population
for pop, ids in population_dict.items():
    with open(os.path.join(args.workingdir, f"{pop}_poplist.txt"), 'w') as pop_file:
        for id in ids:
            pop_file.write(f"{id}\n")

#Get PIXY Pi (In Progress)
#bgzip vcf
#subprocess.call(["bgzip %s > " % (args.variants, args.variants+'.gz')], shell=True)
#subprocess.call(["tabix %s" % (args.variants)], shell=True)

#Run Pixy
#subprocess.call(["pixy --vcf %s --populations %s --window_size %s" % (args.workingdir+args.variants, args.popfile, args.window)], shell=True)

#run vcftools Pi calculation
for popfiles in os.listdir(args.workingdir):
    if popfiles.endswith('_poplist.txt'):
        subprocess.call(["vcftools --vcf %s --keep %s --site-pi --out %s" % (args.variants, popfiles, popfiles[:-12])], shell=True)
        subprocess.call(["vcftools --vcf %s --keep %s --TajimaD %s --out %s" % (args.variants, popfiles, args.tajimawind, popfiles[:-12])], shell=True)
    else:
        continue

#Get Averages:
#Get Pi averages
for file in os.listdir(args.workingdir):
    if file.endswith("sites.pi"):
        #replace tab with comma
        replace_in_file(file, '\t', ',')
        #replace nan with 0
        replace_in_file(file, ' nan', '0')
        append_column_averages(file)
        #rename file to have csv extension
        os.rename(file, file[:-9]+'_Pi.csv')

#Make filtered vcf file
mafstr=args.maf.replace('.', '')
subprocess.call(["vcftools --vcf %s --maf %s --thin %s --recode --recode-INFO-all --out maf%s_thin%s_snps" % (args.variants, args.maf, args.thin, mafstr, args.thin)],shell =True)
filtvcf = "maf"+mafstr+"_thin"+args.thin+"_snps.recode.vcf"

# #Get per individual heterozygosity, hardy weinberg for FIS, and Tajima's D
for popfiles in os.listdir(args.workingdir):
    if popfiles.endswith('_poplist.txt'):
        subprocess.call(["vcftools --vcf %s --keep %s --het --out %s" % (filtvcf, popfiles, popfiles[:-12])], shell=True)
        subprocess.call(["vcftools --vcf %s --keep %s --hardy --out %s" % (filtvcf, popfiles, popfiles[:-12])], shell=True)
    else:
        continue

os.chdir(args.workingdir)

#process hardy output for Obs/Exp Hom/Het and calculating FIS, others

for file in os.listdir(args.workingdir):
    if file.endswith(".hwe"):
        with open(file, 'r') as infile:
            for line in infile:
                if 'P_HET_EXCESS' in line:
                        print(line)
                        name ='CHR,POS,OBS_HOM1,OBS_HET,OBS_HOM2,EXP_HOM1,EXP_HET,EXP_HOM2,hiSq_HWE,P_HWE,P_HET_DEFICIT,P_HET_EXCESS' + '\n'
                        print(name)
                        #replace header
                        replaceAll(file, line, name)
                        #replace / with comma
                        replace_in_file(file, '/', ',')
                        #replace tab with comma
                        replace_in_file(file, '\t', ',')
                        #replace -nan to 0
                        replace_in_file(file, '-nan', '0')
                        # get averages
                        append_column_averages(file)
                        #rename file to have csv extension
                        os.rename(file, file[:-4]+'_HWE.csv')

#Convert raw values of Hobserved Hetobserved Hexpected into proportions

def getprop(directory):
    # Iterate over each file in the directory
    for filename in os.listdir(directory):
        if filename.endswith('_HWE.csv'):
            file_path = os.path.join(directory, filename)
            
            # Read the CSV file into a DataFrame
            df = pd.read_csv(file_path)
            
            # Check if the required columns are present
            if 'OBS_HOM1' not in df.columns or 'OBS_HOM2' not in df.columns or 'OBS_HET' not in df.columns or 'EXP_HET' not in df.columns:
                print(f"Skipping {filename}: Required columns not found.")
                continue
            
            # Create the 'OBS_HOM_Total' column
            df['OBS_HOM_Total'] = df['OBS_HOM1'] + df['OBS_HOM2']
            
            # Calculate the sum of 'OBS_HOM_Total' and 'OBS_HET' for each row
            sum_hom_het = df['OBS_HOM_Total'] + df['OBS_HET']
            
            # Create the 'OBS_HOM_prop' column
            df['OBS_HOM_prop'] = df['OBS_HOM_Total'] / sum_hom_het
            
            # Create the 'OBS_HET_prop' column
            df['OBS_HET_prop'] = df['OBS_HET'] / sum_hom_het
            
            # Create the 'EXP_HET_prop' column
            df['EXP_HET_prop'] = df['EXP_HET'] / sum_hom_het
            
            # Save the modified DataFrame back to the same file or a new file

            output_file = os.path.join(directory, f"prop_{filename}")
            df.to_csv(output_file, index=False)
            print(f"Processed {filename}, saved as {output_file}")

getprop(args.workingdir)
#Calculate FIS from HWE output
#Define function to calculate FIS=(He-Ho)/He; if He = 0; output is 0

def calculate_fis(directory):
    # Iterate over each file in the directory
    for filename in os.listdir(directory):
        if filename.startswith('prop_'):
            file_path = os.path.join(directory, filename)
            
            # Read the CSV file into a DataFrame
            df = pd.read_csv(file_path)
            
            # Check if the required columns are present
            if 'EXP_HET' not in df.columns or 'OBS_HET' not in df.columns:
                print(f"Skipping {filename}: Required columns not found.")
                continue
            
            # Exclude the last row from the FIS calculation
            data_to_calculate = df.iloc[:-1]  # Skip the last row for FIS calculation
            
            # Calculate FIS and handle division by zero
            df['FIS'] = data_to_calculate.apply(
                lambda row: (row['EXP_HET_prop'] - row['OBS_HET_prop']) / row['EXP_HET_prop']
                if row['EXP_HET_prop'] != 0 else 0, axis=1)
            
            # Calculate the average FIS, treating NaN as 0
            average_fis = df['FIS'].fillna(0).mean()  # Fill NaN with 0 before calculating mean
            
            # Output the average FIS value in the excluded last row
            df.at[len(df)-1, 'FIS'] = average_fis  # Place the average in the last row's FIS column
            
            fileout=filename[5:]

            # Output the updated DataFrame back to the same file or a new file
            output_file = os.path.join(directory, f"FIS_{fileout}")
            df.to_csv(output_file, index=False)
            print(f"Processed {filename}, saved as {output_file}")

#Get FIS
calculate_fis(args.workingdir)

#remove intermediate files
for file in os.listdir(args.workingdir):
    if file.startswith('prop_'):
        os.remove(file)

#Get Tajima's D averages
for file in os.listdir(args.workingdir):
    if file.endswith("Tajima.D"):
        #replace tab with comma
        replace_in_file(file, '\t', ',')
        #replace nan with 0
        replace_in_file(file, ' nan', '0')
        append_column_averages(file)
        #rename file to have csv extension
        os.rename(file, file[:-9]+'_TjD.csv')

#compile population statistics into one file
def compile_popstats(directory, output_file):
    # Initialize a list to store the compiled data
    compiled_data = []
    
    # Identify all unique populations based on the filenames
    population_names = set(
        filename.split('_')[0] 
        for filename in os.listdir(directory) 
        if filename.endswith(('_Pi.csv', '_TjD.csv', '_HWE.csv')) and not filename.startswith(('FIS_','prop_'))
    )
    
    # Iterate over each population
    for population in population_names:
        # Initialize a dictionary to store values for this population
        population_data = {'population': population}
        
        # Process the _Pi.csv file
        pi_file = os.path.join(directory, f"{population}_Pi.csv")
        if os.path.exists(pi_file):
            df_pi = pd.read_csv(pi_file)
            population_data['Pi'] = df_pi['PI'].iloc[-1]  # Get the last value in the PI column
        else:
            population_data['Pi'] = None
        
        # Process the _TjD.csv file
        tjd_file = os.path.join(directory, f"{population}_TjD.csv")
        if os.path.exists(tjd_file):
            df_tjd = pd.read_csv(tjd_file)
            population_data['TjD'] = df_tjd['TajimaD'].iloc[-1]  # Get the last value in the TajimaD column
        else:
            population_data['TjD'] = None
        
        # Process the _HWE.csv file
        hwe_file = os.path.join(directory, f"FIS_{population}_HWE.csv")
        if os.path.exists(hwe_file):
            df_hwe = pd.read_csv(hwe_file)
            population_data['FIS'] = df_hwe['FIS'].iloc[-1]  # Get the last value in the FIS column
            population_data['EXP_HET'] = df_hwe['EXP_HET_prop'].iloc[-1]  # Get the last value in the EXP_HET column
            population_data['OBS_HET'] = df_hwe['OBS_HET_prop'].iloc[-1]  # Get the last value in the OBS_HET column
        else:
            population_data['FIS'] = None
            population_data['EXP_HET'] = None
            population_data['OBS_HET'] = None
        
        # Append the population data to the compiled list
        compiled_data.append(population_data)
    
    # Convert the compiled data to a DataFrame
    compiled_df = pd.DataFrame(compiled_data)
    
    # Save the compiled DataFrame to a CSV file
    compiled_df.to_csv(output_file, index=False)
    print(f"Compiled data saved to {output_file}")

#Compilte population statistics
outname='populations_statistics_maf'+mafstr+'_thin'+args.thin+'.csv'
compile_popstats(args.workingdir, outname)


#Concatenate individual 'het' dataframes into one csv:
import os
import pandas as pd

def concatenate_het_files(directory, output_file='ind_population_stats.csv'):
    # Initialize an empty list to store dataframes
    dataframes = []
    
    # Iterate over each file in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.het'):
            file_path = os.path.join(directory, filename)
            
            # Read the TSV file into a DataFrame
            df = pd.read_csv(file_path, sep='\t')
            
            # Extract the population name from the filename (excluding the '.het' suffix)
            population_name = filename.replace('.het', '')
            
            # Add the 'population' column with the extracted population name
            df['population'] = population_name
            
            # Append the dataframe to the list
            dataframes.append(df)
    
    # Concatenate all the dataframes into one, ensuring headers are not repeated
    combined_df = pd.concat(dataframes, ignore_index=True)
    
    # Save the combined dataframe as a CSV file
    combined_df.to_csv(os.path.join(directory, output_file), index=False)
    print(f"Combined data saved as {output_file}")

# Concatenate het files into single csv
concatenate_het_files(args.workingdir)


#move intermediate output files to popstat_files directoy
def move_popstats(directory):
    # Define the new folder name
    new_folder = os.path.join(directory, 'popstat_files')
    
    # Create the new folder if it doesn't exist
    if not os.path.exists(new_folder):
        os.makedirs(new_folder)
    
    # List of extensions to move
    extensions = (
        'poplist.txt', 'Pi.csv', '.het', '.hwe', 'sites.pi', '.log', '_HWE.csv', 'TjD.csv'
    )
    
    # Iterate over all files in the directory
    for filename in os.listdir(directory):
        # Skip the popstat_files folder itself
        if filename == 'popstat_files':
            continue
        
        # Get the full file path
        file_path = os.path.join(directory, filename)
        
        # Check if the file ends with any of the specified extensions
        if os.path.isfile(file_path) and filename.endswith(extensions):
            shutil.move(file_path, os.path.join(new_folder, filename))
            print(f"Moved {filename} to {new_folder}")


move_popstats(args.workingdir)



