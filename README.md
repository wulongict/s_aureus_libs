# s_aureus_libs

Building spectral library from 111 data files. 


## Create spectral library from 111 open search results.
For each search result file, processed by xinteract (PeptideProphet/iProphet), create a spectral library .

```
/path/to/tpp51/bin/spectrast -cP0.0 -cq0.01 -cIHCD -co -cNXYZ ../mzXMLs/interact-XYZ.ipro.pep.xml
```


