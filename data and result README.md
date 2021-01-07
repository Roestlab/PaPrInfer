# Location of the Data and Results

Due to the size of the data files, they are not 
stored on github. As well, the nature of sqlite 
files makes the result and the data exist in the 
same file

There are 4 data set that is intended for this
program, but only one of them was actually
in use as of 2021/01/01. 

1. The one that is in use is on the lab google drive. 
(Project > diaPASEF)
Name:
20180911_TIMS2_12-2_AnBr_SA_diaPASEF_Test10_42eV_1_A1_01_2927_test_lib_overridegroupid.osw

2. The other data set that is also on the google drive is
not used because it does not have
degenerate peptides

3. The other one 20180911_TIMS2_12-2_AnBr_SA_diaPASEF_Test10_42eV_1_A1_01_2927_test_lib_hela_maxlib.osw
does not have pyprophet scoring and has a empty 
feature_ms2 table (so running pyprophet is not possible)
and this program uses the q-value of the peptides 
which comes from pyprophet scoring

4. The last data set is at graham cluster /project/6011811/frankmak/tims/Mar2018_helaBenchmark/data/openswath/ProteinIsoform/merged.osw
