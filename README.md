# Building S. aureus spectral library  


## Create spectral library from 111 open search results.
For each search result file, processed by xinteract (PeptideProphet/iProphet), create a spectral library . Note that we have used IRT peptide defined in file cirt_peps.txt. The iRT peptides are from the following paper. 

>Parker S J, Rost H, Rosenberger G, et al. Identification of a Set of Conserved Eukaryotic Internal Retention Time Standards for Data-independent Acquisition Mass Spectrometry*[S][J]. Molecular & Cellular Proteomics, 2015, 14(10): 2800-2813. 
>
>https://www.sciencedirect.com/science/article/pii/S1535947620326335?via%3Dihub



```
/path/to/tpp51/bin/spectrast -c_IRT./cirt_peps.txt -c_IRR -cP0.0 -cq0.01 -cIHCD -co -cNXYZ interact-XYZ.ipro.pep.xml
```

## Remove the PSMs with a delta mass grreater than 0.02
For each spectral library, create another library with only delta mass close to zero Da. 

```
spectrast -c_MDF0.02 -cN${name}_md0 ${name}.splib
```

## Merge 111 mdf0 spectral libraries
Combine the 111 libraries with option -cJU -cAC

```bash
./merge_111_mdf0_splib.bash

```

## Create library of PSMs with delta mass greater than 0.02
Twelve delta mass values are selected. Next, we will choose PSMs with those PTMs, but localized on a high abundent site. The modified sites are collected from pep.xml files.

### Convert pep.xml into tsv format
In this step, we use an inhouse tool to convert pep.xml file into a tsv file. in total we have 111 tsv files. 

### Build Peptide-deltaMass vs localization information.
Using 111 tsv files, genearte a table of localization information for each (peptide, deltamass) pair. Each row is a distinct (peptide, deltamass ) pair, each column is a data file. The content in cell of row i and column j  of the table, is the distribution of the delta mass on different sites. Assume the  the content of cell of row i, column j is  :

1,4; 2,9

Then it means, for the given peptide, delta mass (PTM) of row i,  there are 13 PSMs identified in data file j, where 4 PSMs got first site modified, 9 PSMs got second site modified. 

### Summarize total frequency 
For each row, calculate the total frequency for each site. 

### Generate AA frequency
Considering the amino acid of each site, the site of 12 PTMs (delta masses) can be summaried as 12 histograms. Each is "sum" of the amino frequency indicated by all the PSMs with that PTM. 

### Filter PSM with amino acid site frequency of a PTM
PSMs with a >15% site frequency are collected from the 111 files to generate the final library. The selected PTMs and the corresponding sites are used in the spectrast.usermods file. 

### Generate TSV file
```python
# First, let's get the list of PSMs. 
import os
import numpy as np
import matplotlib.pyplot as plt 

# get list of tsv table.
def getListOfPepMassDiffFileNames():
    table_list = '/data/wulong/data/jordyLibrary/opensearch/mzXMLs/peptide_massdiff_loc_111.list'
    filepath= os.path.split(table_list)[0]

    filenames = []
    with open(table_list,'r') as f:
        filenames = f.readlines()
        filenames = [os.path.join(filepath, os.path.split(each.strip())[1]) for each in filenames]
    return filenames

# final format is 
# 
#
# mzMLFile        scanNumber      peptideID       probability     scores  protein
# --------------------------------------------------------------------------------------------------------------
# /path/to/data/DrugCombination_S1D_S1D-1.mzXML 6       K.HYYEQPK.F/3   0.9799  0.00024664660000
#     Q2FXA4
# /path/to/data/DrugCombination_S1D_S1D-1.mzXML 8       R.LNFSHGSHEEHK.G/2      0.9241  0.0020882680
# 0000        Q2FXM9
# /path/to/data/DrugCombination_S1D_S1D-1.mzXML 9       K.HYYEQPK.F/3   0.9108  0.00172646400000
#     Q2FXA4

def tokenize_modified_peptides(peptide):
    pep_token=[]
    if peptide[0] != 'n':
        pep_token.append('n')

    k=0
    while k < len(peptide):
        if peptide[k]=='[':
            # find first left bracket
            pep_token[-1] += peptide[k]
            while peptide[k] !=']':
                k+=1
                pep_token[-1] += peptide[k]
            # token done
            k+=1
        else:
            pep_token.append(peptide[k])
            k+=1

    # print(f'tokenize of peptide {peptide}\t {pep_token}\n of row {row}')
    return pep_token 

class modInfo:
    def __init__(self, name, oldToken, massDiff, newToken):
        self.name=name
        self.oldToken= oldToken
        self.massDiff = massDiff 
        self.newToken = newToken 



def getModList():

    # C,Carbamidomethyl: 5007
#                   K,Acetyl: 6886
#                   M,Oxidation: 3260
#                   S,Phospho: 698
#                   W,Oxidation: 9090 ; W,doubly_oxidization: 9997
#                   Y,Oxidation: 6107
#                   n,Acetyl: 16 ; n,Addition_of_Lysine: 39365 ; n,Carbamidomethyl: 20791 ; n,Carbamyl: 6061 ; n,Met-Loss: 528 
    modList={}
    modList[16]=modInfo("Oxidation",['W','Y'], 16, ['W[202]', 'Y[179]'])
    modList[32]=modInfo("Dioxidation",['W'], 32, ['W[218]'])
    modList[80]=modInfo("Phospho",['S','T'], 80, ['S[167]','T[181]'])
    modList[42]=modInfo("Acetyl",['K'], 42, ['K[170]'])
    modList[-131]=modInfo("Met-Loss",['n'],-131,['n[-129]'])
    modList[57]=modInfo("Carbamidomethyl",['n'],57,['n[58]'])
    modList[43]=modInfo("Carbamyl",['n'],43,['n[44]'])
    modList[128]=modInfo("Lys",['n'],128,['n[129]'])
    return modList 

def replaceToken(modinfo,pep_token,loc):
    site_i = -1
    site_OK=False 
    for k in loc:
        for m in range(len(modinfo.oldToken)):
            ost, nst=modinfo.oldToken[m], modinfo.newToken[m]
            # print(f"old site {ost} to new site {nst}")
            offset = 0
            if ost == 'n':
                offset = -1

            if pep_token[k+offset] == ost:
                # matche to first site.
                # great!
                pep_token[k+offset] = nst
                site_OK=True
                site_i = k
                break
        if site_OK:
            break 
    return site_OK, pep_token


filelist = getListOfPepMassDiffFileNames()
modList = getModList()
print(f"file number {len(filelist)}")
GoodDeltaMass = [ 42.010565, 57.021464, 43.005814,  31.989829, 128.094963, -131.040485, 15.994915, 79.966331]
selectedPTM=np.array(GoodDeltaMass)
mass_tol = 0.02
output_table=[]

for i, eachfile in enumerate(filelist):
    # if i > 7:
    #     break
    # for each file, search for the dm of 128 and site of n term. print out it!
    print(f"..{i}..",end='')
    lines = []
    with open(eachfile, 'r') as f:
        lines = f.readlines()
    
    lines=[each.rstrip('\n').split('\t') for each in lines]
    data=lines[1:]
    header = lines[0]
    mod_pep_i = header.index('modified_peptide')
    massDiff_i = header.index('massDiff')
    loc_i = header.index('localization')

    # print(f"{header}\n{data[0]}")
    # print(f"{data[0][mod_pep_i]}\t{data[0][massDiff_i]}")
    for j in range(len(data)):
        row = data[j]
        peptide = row[mod_pep_i]
        dm = float(row[massDiff_i])
        res = selectedPTM[(selectedPTM<dm+mass_tol) & (selectedPTM > dm-mass_tol)]
        if res.shape[0] == 0 :
            # selectedPTM not found, go for next ; 
            continue;
        # Yes, find selected PTM 
        loc = row[loc_i]
        if loc == '':
            # print(f'no localzliation from MSFragger. sip this PSM for {dm}\n{row}')
            continue
        loc = loc.split('_')
        loc = [int(each) for each in loc]
        # Yes, find localization information
        # process each case
        dm_int = round(dm)
        if dm_int not in modList.keys():
            print(f"invalid dm {dm_int} {dm}")
            continue 
        else:
            pep_token = tokenize_modified_peptides(peptide)
            site_OK, pep_token = replaceToken(modList[dm_int],pep_token[:], loc)            
            
            if not site_OK:
                continue

            if pep_token[0] == 'n':
                pep_token[0] = ''
            
            peptide=''.join(pep_token)
            # print(f"newPeptide is {peptide}")
            output_str=f"/data/wulong/data/jordyLibrary/opensearch/mzXMLs/{'.'.join(row[0].split('.')[:-3])}.mzXML\t{row[1]}\t{peptide}/{row[4]}\t{row[9]}\t{1.0}\t{row[11]}"
            output_table.append(output_str)

print(f".Done")

        

# prepare to save the table
# print(f"final row \n{output_table[-1]}")
finaltable = 'selected_ptm_for_tsv_libimportor.tsv'
with open(finaltable,'w') as f:
    f.write('\n'.join(output_table))

print(f'{finaltable} saved! row num {len(output_table)}')
```


### Create the user mod file
To use the TSVImporter function, we need to provide the following usermods file. 
```usermods
n|+128.094963|Lys
n|+43.005814
n|+57.021464
n|-131.040485|Met-Loss
W|+15.994915
Y|+15.994915
W|+31.989829|Dioxidization
S|+79.966331
T|+79.966331
K|+42.010565
```

### Generate library from tsv file
Generate the library file and refresh the database. Making consensus. 

```bash
./spectrast  -Mspectrast.usermods -cD/path/to/db/<seq-db-name>.fasta -cIHCD -cP0.0 -cNselected_ptm selected_ptm_for_tsv_libimportor.tsv 
./spectrast -cJU -cAC -cD/path/to/db/<seq-db-name>.fasta -cNselected_ptm_consensus selected_ptm.splib
```

## Combine the spectral library, with PTM (deltaMass > 0.02), and without PTM (deltaMass<0.02>)

```bash
./spectrast -Mspectrast.usermods -cJU -cAC -cNs_aureus_111_mod_unmod s_aureus_111_mdf0.splib selected_ptm_consensus.splib
```
Quality control with SpectraST
```bash
splib_name=s_aureus_111_mod_unmod;  ./spectrast -Mspectrast.usermods -cAQ -cN${splib_name}_Q ${splib_name}.splib
```

## The S. aureus spectral library
There are two tar.gz files, the one named s_aureus_111_mod_unmod_Q.tar.gz contains 10 DECOY PSMs. For searching purpose, use the other one, s_aureus_111_mod_unmod_Q_nodecoy.tar.gz, where the DECOY PSMs are removed. 


