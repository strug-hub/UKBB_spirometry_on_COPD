#!/bin/bash

pruned1KGplusUKBB="/hpf/largeprojects/struglis/cfcentre/strug/Naim/uk_biobank/Ethnic_PCA/12-unrelated_1KG_plus_UKBiobank_pruned_set_mind0.1"
pcadir="data/intermediate_files/pca"
copddata="data/clean/ukbb_spiro_and_geno_qc_v2.csv"

[ ! -d "$pcadir" ] && mkdir "$pcadir"

# The 1000 Genomes and UKBB genotype arrays have been combined and pruned before, so create symbolic links to them
ln -s "${pruned1KGplusUKBB}.bed" "$pcadir"
ln -s "${pruned1KGplusUKBB}.bim" "$pcadir"
ln -s "${pruned1KGplusUKBB}.fam" "$pcadir"
ln -s "${pruned1KGplusUKBB}.log" "$pcadir"

# intermediate files:
keepfile="${pcadir}/ukbb_gold2-4_copd_cases.txt"
plinkout="${pcadir}/13-ukbb_copd_set"
#kingout="${pcadir}/14-copdset_kinshipmat"
flashpcaout="${pcadir}/15-ukbb_copd_flashpca2"

sed '1d' "$copddata" |awk -F "," '{if($10=="TRUE") print $1"\t"$1}' > "$keepfile"
plink --bfile "$pruned1KGplusUKBB" --keep "$keepfile" --make-bed --out "$plinkout"
#king -b "${plinkout}.bed" --kinship --prefix "$kingout"

# Run PCA for the set with COPD (i.e. GOLD2-4)
#Rscript ~/scripts/R_scripts/PC-AiR_v1_empty_kinFile2.R "$plinkout" "$kingout" "$pcairout"
# fails due to very large memory requirement; run flashPCA2 instead

Rscript ~/scripts/run_flashPCA2.R "$plinkout" "$flashpcaout"



# Next, we want to run PCA for an analysis of all UKBB with good spirometry/genotype data:
ukbbdataQCd="data/clean/ukbb_spiro_and_geno_qc_v2.csv"
keepfile="${pcadir}/ukbb_good_spirometry_and_genoQC_set.txt"
plinkout="${pcadir}/16-ukbb_spiro_set"
#kingout="${pcadir}/17-ukbbspiro_kinshipmat"
flashpcaout="${pcadir}/18-ukbb_ukbbspiro_flashpca2"

sed '1d' "$ukbbdataQCd" |awk -F "," '{print $1"\t"$1}' > "$keepfile"
plink --bfile "$pruned1KGplusUKBB" --keep "$keepfile" --make-bed --out "$plinkout"
#king -b "${plinkout}.bed" --kinship --prefix "$kingout" # very large matrix and not using it in PCA step below

# Run PCA for the set of UKBB individuals that underwent spirometry and genotyping QC
# Rscript ~/scripts/R_scripts/PC-AiR_v1_empty_kinFile2.R "$plinkout" "$kingout" "$pcairout"
# fails due to very large memory requirement; run flashPCA2 instead

Rscript ~/scripts/run_flashPCA2.R "$plinkout" "$flashpcaout"



#===================================
# PCA for GOLD3-4 subset
#===================================
# intermediate files:
keepfile="${pcadir}/ukbb_gold3-4_copd_cases.txt"
plinkout="${pcadir}/19-ukbb_gold34_set"
flashpcaout="${pcadir}/20-ukbb_gold34_flashpca2"

sed '1d' "$copddata" |awk -F "," '{if($11==3 || $11==4) print $1"\t"$1}' > "$keepfile"
plink --bfile "$pruned1KGplusUKBB" --keep "$keepfile" --make-bed --out "$plinkout" --memory 15000 

Rscript ~/scripts/run_flashPCA2.R "$plinkout" "$flashpcaout"


#===================================
# PCA for GOLD2 subset
#===================================
# intermediate files:
keepfile="${pcadir}/ukbb_gold2_cases.txt"
plinkout="${pcadir}/21-ukbb_gold2_set"
flashpcaout="${pcadir}/22-ukbb_gold2_flashpca2"

sed '1d' "$copddata" |awk -F "," '{if($11==2) print $1"\t"$1}' > "$keepfile"
plink --bfile "$pruned1KGplusUKBB" --keep "$keepfile" --make-bed --out "$plinkout" --memory 15000 

Rscript ~/scripts/run_flashPCA2.R "$plinkout" "$flashpcaout"
