### Sketch for computing DICE scores between genomes, from input FASTA files
# Richard Mayne, 2019
# v1

# HEADER /
import pandas as pd
import os
#import psutil

# For dev use (from python-resources)
#import resource

### globalVars
df2 = []
df3 = []
bigdict = {}
subdict = {}
qry_gsubdict = {}
sbj_gsubdict = {}

def scrape_data():
    
    ### Gather data from all .txt files in defined pathway
    for file in filelist:
        filename = os.fsdecode(file)
        if filename.endswith( ('.txt') ):
            df = pd.read_table(file, header=None)
            #assign column names
            default_outfmt6_cols = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'.strip().split(' ')
            df.columns = default_outfmt6_cols
            #filter columns in dataframe
            df_filtered = df[(df['pident'] >= 30) & (df['length'] >= 30) & (df['evalue'] <= 0.01)]
            #handle empty entries       
            if df_filtered.empty:
                Query = df.iloc[0,0]
                Subject = df.iloc[0,1]
                Total = 0
            else:
                Query = df_filtered.iloc[0,0]
                Subject = df_filtered.iloc[0,1] 
                Total = df_filtered['bitscore'].sum()
            df2.append((Query, Subject, Total))

def assign_to_dicts():
    for i in df2.index:
        qry_genome = df2.iloc[i,0]
        sbj_genome = df2.iloc[i,1]
        score = df2.iloc[i,2]

    #if the query genome exists
        if qry_genome in bigdict:
            subdict = bigdict[qry_genome]
            subdict[sbj_genome] = score
            bigdict[qry_genome] = subdict

        #else assign to dictionary
        else:
            subdict = {sbj_genome : score}
            bigdict[qry_genome] = subdict

def analyse():        
    # for each key in the list of genomes
    for i in genomes:
    # get  the qry_genome subdict from the bigdict
        qry_gsubdict = bigdict[i]
    # for each sbj_genome in the subdict
        for key in qry_gsubdict.keys():
    # get the sbj_genome subdict from the bigdict
            sbj_gsubdict = bigdict[key]
    # calculate the distance (need to convert the values to floats to allow operators)
            x = 1 - ((float(qry_gsubdict.get(key))) + (float(sbj_gsubdict.get(i))))/((float(qry_gsubdict.get(i))) + (float(sbj_gsubdict.get(key))))
            #print(i, key, x)


            df3.append((i, key, x))

# HEADER End \


### Folder Select
topDir = os.getcwd()

## Folder Selection Option 1
#lowDir = topDir + '/tblastx'

## Folder Selection Option 2
print('Please type the name of the input directory containing FASTA output')
usrResponse = str(input(''))
lowDir = (topDir + '//' + usrResponse)

os.chdir(lowDir)
filelist = os.listdir(lowDir)

### Scrape call
scrape_data()

#convert df2 list to dataframe 
df2 = pd.DataFrame(list(df2),columns=['Query','Subject','Total'])

assign_to_dicts()

# create a list of the keys as 'genomes'
genomes = list(bigdict.keys())

analyse()

df3 = pd.DataFrame(list(df3),columns=['Query','Subject','Distance'])

os.chdir(topDir)
try:
    os.mkdir('Results')
except:
    print('I cant overwrite your results folder, move it somewhere then try again')
    
os.chdir('.//Results')
df3.to_csv(r'genome_comparison.txt', index = None, header=True)
print('Im done')


# Resource monitor for dev purposes
#process = psutil.Process(os.getpid())
#print(process.memory_info().rss)
