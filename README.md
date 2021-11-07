# s_aureus_libs

Building spectral library from 111 data files. 


## Create spectral library from 111 open search results.
For each search result file, processed by xinteract (PeptideProphet/iProphet), create a spectral library .

```
/path/to/tpp51/bin/spectrast -cP0.0 -cq0.01 -cIHCD -co -cNXYZ interact-XYZ.ipro.pep.xml
```

## Remove the PSMs with a delta mass grreater than 0.02
For each spectral library, create another library with only delta mass close to zero Da. 

```
spectrast -c_MDF0.02 -cN${name}_md0 ${name}.splib
```

## Merge 111 mdf0 spectral libraries
Combine the 111 libraries with option -cJU -cAC

```
./merge_111_mdf0_splib.bash

```

## Create library of PSMs with delta mass greater than 0.02
Twelve delta mass values are selected. 
