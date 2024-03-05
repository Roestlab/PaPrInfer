# Location of the Data and Results

Due to the size of the data files, they are not 
stored on GitHub. As well, the nature of sqlite 
files makes the result and the data exist in the 
same file

There are 4 data set that is intended for this
program, but only one of them was actually
in use as of 2021/01/01. 

1. The one that is in use is on the lab Google Drive. 
(Project > diaPASEF)
Name:


    20180911_TIMS2_12-2_AnBr_SA_diaPASEF_Test10_42eV_1_A1_01_2927_test_lib_overridegroupid.osw

2. The other data set that is also on the Google Drive is
not used because it does not have
degenerate peptides

3. The other one does not have pyprophet scoring and has an empty 
feature_ms2 table (so running pyprophet is not possible)
and IDPicker uses the q-value of the peptides 
which comes from pyprophet scoring.


    20180911_TIMS2_12-2_AnBr_SA_diaPASEF_Test10_42eV_1_A1_01_2927_test_lib_hela_maxlib.osw


4. The last data set is at graham cluster /project/6011811/frankmak/tims/Mar2018_helaBenchmark/data/openswath/ProteinIsoform/merged.osw

Update 2022/09/26

Since 1 has problems with epifany, I use another
set of data (5 and 6)

5. This is the dataset using swissprot database Top two are the same, just one in tsv and one in osw, most of the time, I only used
    the osw one.


    projects/def-hroest/data/diaPASEFManuscriptOutput/merged_-1.osw
   
    projects/def-hroest/data/diaPASEFManuscriptOutput/merged_-1.tsv
   
    projects/def-hroest/data/diaPASEFManuscriptOutput/pyprophet_export.tsv 


7. This is the dataset using uniprot/trembl database


    projects/def-hroest/data/josh_kai_2022/2021-12-23-proteinIsoformRslts.osw
   
    projects/def-hroest/data/josh_kai_2022/peptide.tsv

peptide.tsv is one of the output files from msfragger, which is needed
for adding non-razor protein.
