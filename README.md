# Building S. aureus spectral library  

## Create spectral library from 111 open search results.
For each search result file, processed by xinteract (PeptideProphet/iProphet), create a spectral library . Note that we have used IRT peptide defined in file cirt_peps.txt. The iRT peptides are from the following paper. 

>Parker S J, Rost H, Rosenberger G, et al. Identification of a Set of Conserved Eukaryotic Internal Retention Time Standards for Data-independent Acquisition Mass Spectrometry*[S][J]. Molecular & Cellular Proteomics, 2015, 14(10): 2800-2813. 
>
>https://www.sciencedirect.com/science/article/pii/S1535947620326335?via%3Dihub


```
spectrast -c_IRT./cirt_peps.txt -c_IRR -cP0.0 -cq0.01 -cIHCD -co -cNXYZ interact-XYZ.ipro.pep.xml
```

## Remove the PSMs with a delta mass greater than 0.02
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
In this step, we use an inhouse tool  (named xmlparser) to convert pep.xml file into a tsv file. in total we have 111 tsv files. 
```bash
xmlparser --inputpepxml interact-XYZ.ipro.pep.xml
```

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
```bash
python ./generate_tsv_for_spectrast.py
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

## Combine the spectral library, with PTM (deltaMass > 0.02), and without PTM (deltaMass<0.02)

```bash
./spectrast -Mspectrast.usermods -cJU -cAC -cNs_aureus_111_mod_unmod s_aureus_111_mdf0.splib selected_ptm_consensus.splib
```
Quality control with SpectraST
```bash
splib_name=s_aureus_111_mod_unmod;  ./spectrast -Mspectrast.usermods -cAQ -cN${splib_name}_Q ${splib_name}.splib
```

## The S. aureus spectral library
There are two tar.gz files, the one named s_aureus_111_mod_unmod_Q.tar.gz contains 10 DECOY PSMs. For searching purpose, use the other one, s_aureus_111_mod_unmod_Q_nodecoy.tar.gz, where the DECOY PSMs are removed. 

```bash
name=s_aureus_111_mod_unmod_Q
spectrast -cf'Protein!~DECOY' -cN${name}_nodecoy ${name}.splib
```


