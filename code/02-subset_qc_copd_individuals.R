library(data.table)
library(RNOmni)
library(httr)
library(jsonlite)
source(".apikey")

# Load data
dat <- fread(file.path("data","intermediate_files","ukb24727_spirometry.tab"), header=T, stringsAsFactors=F) # 502,543
genotyped_fam <- fread("/hpf/largeprojects/struglis/datasets/uk_biobank_40946/genotypes/ukb40946_cal_chr1_v2_s488295.fam", header=F) # 488,377

x <- merge(genotyped_fam[,'V1',drop=F], dat, by.x="V1", by.y="f.eid") # 488,295
# phenotype and genotype files' EID's overlap for 488,295 individuals

dat <- x # 488,295
setnames(dat,"V1","f.eid")

# Step 1. From the genotyped individuals, identify those that fail
# 1. Remove failed QC individuals (column 8236, UDI=22010-0.0; poor heterozygosity/missingness)
dat <- x[!(f.22010.0.0 %in% 1)] # 487,826; 469 removed

# 2. Remove individuals with no genetic sex information 
# or genetic sex not matching reported gender
dat <- dat[!is.na(f.22001.0.0)] # 487,826;; none removed
dat <- dat[-which(f.22001.0.0 != f.31.0.0)] # 487,449; 377 removed

# 3. Remove related individuals
relinds <- fread(file.path("data","intermediate_files", "set_of_related_ind_to_rm.txt"), header=F, stringsAsFactors=F) # 36,100
# check all fid==iid
# for(i in 1:nrow(relinds)) {
#   if(relinds[i,1] != relinds[i,2]) print(relinds[i,])
# }
#V1      V2
#1: VN061 HG02061 --> not in the dataset
relinds <- relinds[V2!="HG02061"]
relinds$V1 <- as.integer(relinds$V1)
relinds$V2 <- as.integer(relinds$V2) # 36,099
dat <- dat[!(f.eid %in% relinds$V2)] # 451,445; 36,004 removed

# 4. Remove samples with sex chromosome aneuploidy
dat <- dat[-which(f.22019.0.0 %in% 1)] # 451,014; 431 removed

# 5. Remove non-Caucasian (column 8193; 22006-0.0)
dat <- dat[f.22006.0.0==1] # 377,006; 74,008 removed


#============================
# Spirometry QC
#============================

## Acceptability of blows
acceptable_blows <- c("ACCEPT", "BELOW6SEC ACCEPT", "BELOW6SEC")
acceptable_blows_mat_idx <- apply(dat[ ,grep("20031.0", colnames(dat)),with=F ], 2, function(x) ifelse(is.na(x) | (x %in% acceptable_blows), 1, 0)) # 797,669 acceptable blows from 377,006 participants

## Assess start of blow quality
source("code/back_ev_calculation.R")
blow_curves_mat <- dat[, grep("3066.0", colnames(dat)), with=F]
acceptable_blow_start <- matrix(rep(NA, nrow(blow_curves_mat)*ncol(blow_curves_mat)), ncol=ncol(blow_curves_mat))
for(i in 1:nrow(blow_curves_mat)) {
  if((i %% 10000) == 0) print(paste0(i,"/",nrow(acceptable_blow_start)," ",format(round(i/nrow(acceptable_blow_start)*100,digits = 2),nsmall=1),"% done"))
  for(j in 1:ncol(blow_curves_mat)) {
    acceptable_blow_start[i,j] <- blow_start_quality(blow_curves_mat[i,j,with=F])
  }
}
fwrite(acceptable_blow_start, "data/intermediate_files/acceptable_blow_starts.csv"
       ,quote=F, row.names=F, col.names=T)
# 752,112 acceptable blow starts from 377,006 participants


x <- unname(acceptable_blows_mat_idx * acceptable_blow_start)
x <- ifelse(x==0,NA,x) # 550,563 acceptable blows and blow starts in 377,006 participants
participants_to_rm <- rowSums(is.na(x))==3 # 91,951 marked for removal


eid <- dat[,'f.eid']
fvc_mat <- dat[, grep("3062.0", colnames(dat)), with=F]
fev_mat <- dat[, grep("3063.0", colnames(dat)), with=F]
pef_mat <- dat[, grep("3064.0", colnames(dat)), with=F]

fvc_mat <- fvc_mat * x
fvc_mat <- fvc_mat[-which(participants_to_rm),]
fev_mat <- fev_mat * x
fev_mat <- fev_mat[-which(participants_to_rm),]
pef_mat <- pef_mat * x
pef_mat <- pef_mat[-which(participants_to_rm),]
eid <- eid[-which(participants_to_rm),] # 285,055


## Get best measures
fvc_best <- apply(fvc_mat,1,max,na.rm=T)
fev_best <- apply(fev_mat,1,max,na.rm=T)
pef_best <- apply(pef_mat,1,max,na.rm=T)

## Assess reproducibility of measures
### best measures have to be within 250mL from any other blow (including unacceptable blows)
nullify <- function(arow) {
  idx <- which(arow %in% max(arow,na.rm=T))[1]
  arow[idx] <- NA
  return(arow)
}
fvc_mat2 <- dat[-which(participants_to_rm), grep("3062.0", colnames(dat)), with=F]
fev_mat2 <- dat[-which(participants_to_rm), grep("3063.0", colnames(dat)), with=F]
fvc_best_colidx <- apply(fvc_mat2,1,function(x) which(x %in% max(x,na.rm=T))[1])
fev_best_colidx <- apply(fev_mat2,1,function(x) which(x %in% max(x,na.rm=T))[1])
fvc_mat2 <- t(apply(fvc_mat2, 1, function(x) nullify(x)))
fev_mat2 <- t(apply(fev_mat2, 1, function(x) nullify(x)))


fvc_reproducible <- abs(fvc_mat2 - fvc_best)
fvc_reproducible <- apply(fvc_reproducible, 1, function(x) any(x<0.25,na.rm=T)) # 268,863
fev_reproducible <- abs(fev_mat2 - fev_best)
fev_reproducible <- apply(fev_reproducible, 1, function(x) any(x<0.25,na.rm=T)) # 275,347

fev_and_fvc_reproducible <- fev_reproducible & fvc_reproducible # 264,873

fvc_best_reproducible <- fvc_best[fev_and_fvc_reproducible] # 264,873
fev_best_reproducible <- fev_best[fev_and_fvc_reproducible] # 264,873
fev_fvc_ratio <- fev_best_reproducible / fvc_best_reproducible # 264,873
eid <- eid[fev_and_fvc_reproducible,] # 264,873

dat <- dat[-which(participants_to_rm),] # 285,055 with acceptable blows
dat <- dat[fev_and_fvc_reproducible,] # 264,873 with reproducible blows

spiro_qc_df <- cbind(eid, fvc_best_reproducible, fev_best_reproducible, pef_best, fev_fvc_ratio) # 264,873
# 40,369 have FEV1/FVC ratio < 0.7

dat <- cbind(dat, spiro_qc_df[,-"f.eid",with=F])


# 5. Save the current dataset to calculate FEV1pp on GLI calculator
gli_query_data <- dat[, grep("eid|21022.0|50.0|22001.0|22006.0|fev_best_reproducible|fvc_best_reproducible|pef_best", colnames(dat)), with=F]

setnames(gli_query_data,
         c("f.eid", "f.50.0.0", "f.21022.0.0", "f.22001.0.0", "f.22006.0.0", "fvc_best_reproducible", "fev_best_reproducible", "pef_best"),
         c("eid", "height","age","sex","ethnic","fvc","fev1", "pef"))
setcolorder(gli_query_data,
            c("eid","age","height","sex","ethnic","fev1","fvc","pef"))
gli_query_data[,sex := ifelse(sex==0,'F',ifelse(sex==1,'M',NA))]
fwrite(gli_query_data
       ,file.path("data", "intermediate_files","gli_calc_data_v2.csv")
       ,quote=F, row.names=F, col.names=T, sep=",")
# 264,640 out of the 264,873 have complete data to calculate FEV1pp


# 7. Use GLI Calculator API to calculate FEV1pp
to_query <- na.omit(gli_query_data)
#boxplot(to_query$fev1)
#boxplot(to_query$fvc) - one clear outlier
#boxplot(to_query$fev1/to_query$fvc)
# remove one extreme outlier
gli_query_data <- gli_query_data[fvc<10]
to_query <- to_query[fvc<10]



rest_url <- "https://gli-api.ersnet.org/public/"
fev1pp <- NULL
results <- NULL
for(i in 1:nrow(to_query)) {
  print(paste0("Acquiring sample ",i," of ",nrow(to_query),
               " ",format(round((i/nrow(to_query))*100
                                ,digits = 2),nsmall=2), "% done"))
  response <- GET(paste0(rest_url,"type/spiro"
                         ,"/age/",to_query$age[i]
                         ,"/height/",to_query$height[i]
                         ,"/sex/",tolower(to_query$sex[i])
                         ,"/ethnic/",to_query$ethnic[i]
                         ,"/fev1/",to_query$fev1[i]
                         ,"/fvc/",to_query$fvc[i]
  )
  , add_headers("x-api-key" = apikey))
  if(response$status_code==200){
    spiro <- content(response, type="application/json")
    fev1pp <- as.numeric(spiro$fev1_pp)
    results <- rbind(results, c(eid=to_query$eid[i], fev1pp=fev1pp))
  } else {
    message <- content(response, type="application/json")$message
    print(paste0("API status code ",response$status_code,": ",message))
    break
  }
}
fwrite(results, "data/intermediate_files/fev1pp_results_v2.csv"
       ,quote=F, row.names=F, col.names=T)
results <- data.table(results)
if(nrow(results) != nrow(gli_query_data)) {
  stop("Only partial FEV1pp results were obtained")
}
results <- data.table(results)
results <- merge(to_query[,-c("fev1pp")], results, by='eid')

gli_query_data <- results # # 264,639

hasCOPD <- with(gli_query_data, (fev1/fvc)<0.7 & fev1pp<80) # 22,199 have COPD as per GOLD2-4 definitions
GOLDlevel <- with(gli_query_data, ifelse((fev1/fvc)<0.7, ifelse(fev1pp>=80,1,
                                                                ifelse(fev1pp>=50,2,
                                                                       ifelse(fev1pp>=30,3,
                                                                              ifelse(fev1pp<30,4,NA)))),NA))
#table(GOLDlevel)
gli_query_data <- cbind(gli_query_data,hasCOPD,GOLDlevel)
samplefile <- fread("/hpf/largeprojects/struglis/datasets/uk_biobank_40946/imputation/sample_files/ukb40946_imp_chr1_v3_s487324.sample")
samplefile <- samplefile[-1,] # 487,409
samplefile$index <- 1:nrow(samplefile)
samplefile$samplename <- paste0("(anonymous_sample_", samplefile$index, ")")
samplefile <- samplefile[,-c("missing","sex","ID_2")]
setnames(samplefile, "ID_1","eid")
i <- which(samplefile$eid %in% gli_query_data$eid) # 264,273 out of 264,873
gli_query_data_samples <- paste0("(anonymous_sample_", i, ")")

x <- merge(gli_query_data, samplefile, by="eid")
gli_query_data <- x # 264,273 (with 22,176 GOLD2-4-defined COPD)
fwrite(gli_query_data, file.path("data","clean","ukbb_spiro_and_geno_qc_v2.csv")
       , quote=F, row.names=F, col.names=T, sep=",")
