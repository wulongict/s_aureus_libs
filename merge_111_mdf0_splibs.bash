files=`ls splib | grep mdf0`
filenum=`echo $files | wc -w `
/tools/tpp51/bin/spectrast -cJU -cAC -cNs_aureus_111_mdf0 $files 