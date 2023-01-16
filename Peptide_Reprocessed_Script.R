#Script for processing proteome discover export data for the identification of bitter peptides in cheese.

options( scipen = 1 )
options( digits = 5 )
##Install and load required packages#######
install.packages("matrixTests")
install.packages("gridExtra")
install.packages("gtools")
install.packages("Matrix")
install.packages("ggplot")
install.packages("dplyr")
install.packages("stringr")
library(gtools)
library(matrixTests)
library(dplyr)
library(matrixStats)
library(ggplot2) 
library(ggrepel)
library(tidyr)
library(Matrix)
library(gridExtra)
library(stringr)
####Import raw data####
#Imports literature reference database
lit_ref <- read.csv("Supplemental_Information_CSV_Database.csv", header = TRUE, stringsAsFactors = FALSE, fileEncoding = "latin1", comment.char = '#')	
lit_ref$Groups <-as.factor(lit_ref$Groups)

#Imports peptide data list
data1 <- read.table("Kuhfeld_Bitter_Peptide_Reprocess-(1)_PeptideGroups.txt", header = TRUE, stringsAsFactors = FALSE, fileEncoding = "latin1",)


#oreder of samples and bitterness grouping in data
#adds sample ID to peptide abundance columns after seperating

#imports file name list and extracts sample name, ID, and grouping
#inpute files line 14 and 15 column file name was manually edited from "_E_8_7_a_" & "_E_8_7_b_" to "_E_8_7a_" & "_E_8_7b_" to be consistent with the formating for the threshold samples
sample_names_df <- read.table("Kuhfeld_Bitter_Peptide_Reprocess-(1)_InputFiles.txt", header = TRUE, stringsAsFactors = FALSE)
sample_names_ext_name <- str_extract(sample_names_df$File.Name, "([A-Z]_[0-9]_[0-9][a-z]?)")
sample_names_ext_rename <- str_replace(sample_names_ext_name, "(?<=[0-9])[_]", ".")
sample_names_df$sample_names <-sample_names_ext_rename
sample_names_ext_ID <-sample_names_df$File.ID
sample_names <- sample_names_ext_rename

#imports sample overview to bring in bitterness and grouping

samdf_import <- read.csv("Sample_Overview.csv", header = TRUE, stringsAsFactors = FALSE)

for (i in 1:nrow(samdf_import)) {
  if (samdf_import[i,'Mean.bitterness.value'] == "ND") {
    samdf_import[i,'Mean.bitterness.value'] = 0.5} # MBI LOQ = 0.5  
  }
colnames(samdf_import)[colnames(samdf_import) == 'Mean.bitterness.value'] <- 'MBI'
samdf_import$MBI <- as.numeric(samdf_import$MBI)
colnames(samdf_import)[colnames(samdf_import) == 'Age'] <- 'age'
colnames(samdf_import)[colnames(samdf_import) == 'Bitterness.grouping'] <- 'grouping'
colnames(samdf_import)[colnames(samdf_import) == 'Categorical.grouping'] <- 'bitterness'
colnames(samdf_import)[colnames(samdf_import) == 'Sample.ID'] <- 'sample_names'
sample_merg_df <- merge(samdf_import, sample_names_df, by.x = "sample_names", all.x = TRUE, sort = FALSE)
samdf <- subset(sample_merg_df, select = c(sample_names, age, MBI, grouping, bitterness, File.ID))
colnames(samdf)[colnames(samdf) == 'File.ID'] <- 'sample_names_ext_ID'

#MBI <- c(0.7,0,0,0,1.2,1.4,1.4,2.1,2.1,2.4,2.9,2.9,3.0,3.4)  #mean bitterness intensity of cheese from NCSU
#age <- c(0.2,0.2,0.2,0.2,3.3,2.9,6.0,7.3,6.2,5.6,7.2,5.7,8.7,8.7)  #age of cheese from creamery (years)
#samdf <- data.frame(MBI,age,sample_names,sample_names_ext_ID)

# samdf$grouping <- NA
# samdf$bitterness <-NA
# for (i in 1:nrow(samdf)){
#   temp_sample_name <-unlist(strsplit (x = samdf[i,'sample_names'], split = ""))
#   if ("T" %in% temp_sample_name)
#     {samdf[i,'grouping'] <- "Threshold" 
#     samdf[i,'bitterness'] <- "Non-bitter"}
#   if ("L" %in% temp_sample_name)
#     {samdf[i,'grouping'] <- "Low" 
#     samdf[i,'bitterness'] <- "Non-bitter"}
#   if ("M" %in% temp_sample_name)
#     {samdf[i,'grouping'] <- "Moderate" 
#     samdf[i,'bitterness'] <- "Bitter"}
#   if ("E" %in% temp_sample_name)
#     {samdf[i,'grouping'] <- "Extreme" 
#     samdf[i,'bitterness'] <- "Bitter"}
# }

data1_names <-colnames(data1)
data1_names <- str_replace(data1_names, "Abundances.Normalized", "ANI") #replaces abundance normalized imputed data with abrevation "ANI"
data1_names <- str_replace(data1_names, "Abundances.Origin", "ADI") #replaces abundance origin with abrevation "ADI" for detected or imputed
data1_names <- str_replace(data1_names, "Abundance", "A") #replaces abundance  with abrevation "A" for data with out normilization or imputation
data1_names <- str_replace(data1_names, ".Sample.....?.?.?.?.?.?", "") #removes sample and grouping infromation

data1_names_new <- list()
for (i in data1_names)  {
  temp_data1_name <-unlist(strsplit (x = i, split = "[/.]"))
  if(temp_data1_name[2] %in% sample_names_ext_ID)
    {
      temp_data2_name <- samdf$sample_names[samdf$sample_names_ext_ID==temp_data1_name[2]]
      temp_comb_name <- paste0(temp_data1_name[1],':',temp_data2_name)
      data1_names_new[[length(data1_names_new) + 1]] <- temp_comb_name    # Append new list element
  }
  else
  {
    data1_names_new[[length(data1_names_new) + 1]] <- i   # Append new list element
  }
}


####data cruching#### 
#relabels peptide abudnaces columns with sample names and merges to create new labled dataframe(datam)
datam <- data1
colnames(datam) <- data1_names_new

#counts the number of real MS detections in the samples
datam$count <- NA
for (i in 1:nrow(datam)){
  datam[i,'count'] <- (sum(c(datam[i,grep("^ADI:",data1_names_new)] == "Det")))
}

####Q value Function####
amino_acids <- rbind(
  'Q' = -100,
  'N' = -10,
  'G' = 0,
  'S' = 40,
  'T' = 440,
  'H' = 500,
  'D' = 540,
  'E' = 550,
  'R' = 730,
  'A' = 730,
  'C' = 0,
  'M' = 1300,
  'K' = 1500,
  'V' = 1690,
  'L' = 2420,
  'P' = 2620,
  'F' = 2650,
  'Y' = 2870,
  'I' = 2970,
  'W' = 3000)
amino_acids <- amino_acids[sort(unique(rownames(amino_acids))), ]

Q_value_results <- c() 
for (i in 1:nrow(datam)){
  peptide <-unlist(strsplit (x = datam[i,"Sequence"], split = ""))
  peptide_vec <- table(peptide)
  Q_value <- sum(amino_acids[names(peptide_vec)]*peptide_vec)/length(peptide)
  Q_value_results <-rbind(Q_value_results,c(datam[i,"Sequence"],Q_value))
}
datam$Q_value <- as.numeric(Q_value_results[,2])

####determines standard difference mean#######
nb_sample_list <-paste0("ANI:",samdf$sample_names[samdf$bitterness == "Non-bitter"])
b_sample_list <-paste0("ANI:",samdf$sample_names[samdf$bitterness == "Bitter"])

datam$mean_NB <- rowMeans(as.matrix(datam[nb_sample_list]), na.rm =TRUE)
datam$mean_B <- rowMeans(as.matrix(datam[b_sample_list]), na.rm =TRUE)
datam$diff_means <- datam$mean_B - datam$mean_NB 
datam$sd_NB <- rowSds(as.matrix(datam[nb_sample_list]), na.rm =TRUE)
datam$sd_B <- rowSds(as.matrix(datam[b_sample_list]), na.rm =TRUE)
n_NB <- sum(samdf$bitterness == "Non-bitter")
n_B <-  sum(samdf$bitterness == "Bitter")
datam$pooled_sd <- sqrt(((n_NB-1)*datam$sd_NB^2 + (n_B-1)*datam$sd_B^2)/ (n_NB+n_B-2))
datam$stand_diff_mean <- datam$diff_means / datam$pooled_sd

# #determines standard difference mean
# data1$mean_T_L <- rowMeans(as.matrix(data1[,c(2,3,4,5,6,7,8)]), na.rm =TRUE)
# data1$mean_M_E <- rowMeans(as.matrix(data1[,c(9,10,11,12,13,14,15)]), na.rm =TRUE)
# data1$diff_means <- data1$mean_M_E - data1$mean_T_L 
# data1$sd_T_L <- rowSds(as.matrix(data1[,c(2,3,4,5,6,7,8)]), na.rm =TRUE)
# data1$sd_M_E <- rowSds(as.matrix(data1[,c(9,10,11,12,13,14,15)]), na.rm =TRUE)
# n_T_L <- 7
#n_M_E <- 7
#data1$pooled_sd <- sqrt(((n_T_L-1)*data1$sd_T_L^2 + (n_M_E-1)*data1$sd_M_E^2)/ (n_T_L+n_M_E-2))
#data1$stand_diff_mean <- data1$diff_means / data1$pooled_sd

####determins p value and fold change of T_L and M_E#####
datam$pvalue <- row_t_welch((datam[nb_sample_list]), (datam[b_sample_list]))$pvalue
datam$pvalue <- as.numeric(datam$pvalue, na.rm = TRUE)
datam$foldchange <-(datam$mean_B/datam$mean_NB)
datam$foldchange <- as.numeric(datam$foldchange, na.rm = TRUE)
datam$logfoldchange <-cbind(log2(datam$foldchange))
datam$logfoldchange <- as.numeric(datam$logfoldchange, na.rm =TRUE)
datam$logpvalue <- cbind(-log10(datam$pvalue))
datam$logpvalue <- as.numeric(datam$logpvalue, na.rm =TRUE)



####converts casein to greek symbol#####
colnames(datam)[colnames(datam) == 'Master.Protein.Accessions'] <- 'Protein'
Beta<-"\U03B2"
BetaA1<- str_c(Beta,"A1")
BetaA2<- str_c(Beta, "A2")
Alphas1<-paste0("\U03B1","s1")
Alphas2<-paste0("\U03B1","s2")
Kappa<-paste("\U03BA")
#regular expressions and for loop to deal duplicate data due to beta casein P02666 without genetic varriant A1 & A2
  datam$Positions.in.Proteins <-gsub("; P02666A2[ ][//[].....?.?.?.?.?","",as.character(datam$Positions.in.Proteins))
for (i in 1:nrow(datam))  {
  temp_pnp <-unlist(strsplit (x=datam[i,'Positions.in.Proteins'], split = " "))
  if (datam[i,'Protein'] == "P02666A1; P02666A2")
    {
    datam[i,'Positions.in.Proteins'] <-paste("P02666",temp_pnp[2])
    datam[i,'Protein'] <- "P02666"
    }
}

#labels position in protein
datam$Positions.in.Proteins <-gsub("P02666",Beta,as.character(datam$Positions.in.Proteins))
datam$Positions.in.Proteins <-gsub("P02666A1",BetaA1,as.character(datam$Positions.in.Proteins))
datam$Positions.in.Proteins <-gsub("P02666A2",BetaA2,as.character(datam$Positions.in.Proteins))
datam$Positions.in.Proteins <-gsub("P02662",Alphas1,as.character(datam$Positions.in.Proteins))
datam$Positions.in.Proteins <-gsub("P02663",Alphas2,as.character(datam$Positions.in.Proteins))
datam$Positions.in.Proteins <-gsub("P02668",Kappa,as.character(datam$Positions.in.Proteins))

#labels protein
datam$Protein[datam$Protein == "P02666"] <- Beta
datam$Protein[datam$Protein == "P02666A1"] <- BetaA1
datam$Protein[datam$Protein == "P02666A2"] <- BetaA2
datam$Protein[datam$Protein == "P02662"] <- Alphas1
datam$Protein[datam$Protein == "P02663"] <- Alphas2
datam$Protein[datam$Protein == "P02668"] <- Kappa

####Filters out data based on count, molecular weight and Q-value####
#number of rows or orginal data set pre-filter
nrow(datam)
      #[1] 2100

#filters for count >=4 (only includes peptides with atleast four instrement responces) creates new data frame called datamfc (merged, filtered, count)
datamfc <- data.frame(subset(datam, datam$count >= 4))
nrow(datamfc)
      #[1] 1140

#filters for molecular weight <= 3000  creates new data frame called datafcmw (merged, filtered, count, molecular weight)
datamfcmw <- data.frame(subset(datamfc, datamfc$Theo.MHplus.in.Da<= 3000))
nrow(datamfcmw)
      #[1] 1105

#filters for Q-Value >=1200  creates new data frame called data_filtered (filtered, count(>4), molecular weight (<=3000),Q-Value (>=1200))
data_filtered <- data.frame(subset(datamfcmw, datamfcmw$Q_value >= 1200))
nrow(data_filtered)
      #[1] 872

#joins sequence and modification columns and checks for duplicates
data_filtered$combined_sequence <-paste0(data_filtered$Sequence,data_filtered$Modifications)
n_distinct(data_filtered$combined_sequence)
      #[1] 872 
nrow(data_filtered) == n_distinct(data_filtered$combined_sequence)
      #[1] TRUE, all values are unique and there are no duplicates

####Calculates linear correlations####
pdf("Bar,box, and scatter plot of peptides abudance to cheese bitterness.pdf",height=8,width=12)
cor.results_A <- c()
cor.results_B <- c()
par(mfrow=c(2,2))
sample_list_ANI <-paste0("ANI.",samdf$sample_names)
for (i in 1:nrow(data_filtered)) {
  data_cor <- data.frame(b=samdf$MBI,age=samdf$age,grp=factor(samdf$grouping,levels=c("Threshold", "Low", "Moderate", "Extreme")),a=c(unlist(data_filtered[i,sample_list_ANI])) )
  cor_B <- cor(data_cor[,c("a","b")])[2,1]
  cor.results_B <-rbind(cor.results_B,c(data_filtered[i,1],cor_B))
  #data_cor <- data.frame(b=bitter,age=age,grp_1E=factor(grp_1E,levels=c("T","L","M","E")),a=c(unlist(data_filtered[i,sample_list_ANI])) )
  cor_A <- cor(data_cor[,c("a","age")])[2,1]
  cor.results_A <-rbind(cor.results_A,c(data_filtered[i,1],cor_A))
#  boxplot(a~grp,data=data_cor,main=paste(data_filtered[i,'Protein'],data_filtered[i,'Positions.in.Proteins'],data_filtered[i,'Instrument_Run']),xlab="bitterness grouping",ylab="Relative abundance")
  dat <- matrix(data_filtered[i,sample_list_ANI],1,14)
  colnames(dat) <- sample_names_ext_rename
#  barplot( dat, main="Relative abundance bar plot", xlab="Sample ID" , ylab="Relative abundance", cex.names=.5)
#  plot(age,data_filtered[i,sample_list_ANI],xlab="Cheese Age (years)",ylab="Relative abundance",main=paste("Relative abundance vs. cheese age","corr=",round(cor_A,3)))
  lm_A <- lm(a ~ age,data=data_cor)
#  abline(lm_A,)
#  plot(samdf$MBI,data_filtered[i,samdf$sample_names],xlab="Mean bitterness score",ylab="Relative abundance",main=paste("Relative abundance vs. cheese bitterness","corr=",round(cor_B,3)))
  lm_B <- lm(a ~ b,data=data_cor)
#  abline(lm_B)
}
dev.off()
#combines peptide, corr_age, and Corr_bitterness into one table  
data_filtered$Corr_age<-as.numeric(cor.results_A[,2])
data_filtered$Corr_bitter<-as.numeric(cor.results_B[,2])

####Creates data frame for export and ML analysis####
#data_rf <- (data_filtered[sample_names_ext_rename])

#generates
#counter <- 2
#rf_row_name_list <- "OUT1"
#while(counter<=nrow(data_rf)){
#  rf_row_name_list <- append(rf_row_name_list, str_c("OUT", counter, sep =""))
#  counter <- counter + 1
#  }

#colnames(data_rf) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6","Sample7","Sample8","Sample9","Sample10","Sample11","Sample12","Sample13","Sample14")
#rownames(data_rf) <- c(rf_row_name_list)
#saveRDS(data_rf, file = "data_rf.Rdata")
#write.csv(data_rf, "data_rf.csv")
#write.csv(samdf, "sample_cats.csv")

####compares linear relationship of age and bitterness to peptide abundance####
cor_age=c(unlist(data_filtered[ ,'Corr_age']))
cor_age=as.numeric(cor_age)
cor_bitter=c(unlist(data_filtered[ ,'Corr_bitter']))
cor_bitter=as.numeric(cor_bitter)
t.test(cor_bitter, cor_age, paired = TRUE, alternative = "two.sided")
# 	Paired t-test
# 
# data:  cor_bitter and cor_age
# t = -7, df = 871, p-value = 4e-12
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.030 -0.017
# sample estimates:
# mean of the differences 
#                  -0.023 
mean(cor_age, na.rm = TRUE)
      #[1] 0.033
mean(cor_bitter, na.rm = TRUE)
      #[1] 0.0098
#plot(cor_age,cor_bitter)

#####generates bitter peptide canadites list####
#data_filtered <- data_filtered[order(data_filtered$stand_diff_mean, decreasing = TRUE),]
#removes bulky extraq columns
A_sample_list <-paste0("A.",samdf$sample_names)
data_filtered[A_sample_list] <-NULL
ADI_sample_list <-paste0("ADI.",samdf$sample_names)
data_filtered[ADI_sample_list] <-NULL
ANI_sample_list <-paste0("ANI.",samdf$sample_names)
data_filtered[ANI_sample_list] <-NULL
data_filtered["Sequence"] <-NULL
data_filtered["Protein"] <-NULL
data_filtered["Peptide.Groups.Peptide.Group.ID"] <-NULL
data_filtered["Modifications"] <-NULL
data_filtered["logpvalue"] <-NULL
data_filtered["logfoldchange"] <-NULL
data_filtered["Quan.Info"] <-NULL

#orders peptide based on rank average meam of stand diff mean and cor to bitter
data_filtered <- data_filtered[order(data_filtered$stand_diff_mean, decreasing = TRUE),]
data_filtered$SMD_order <- order(data_filtered$stand_diff_mean, decreasing = TRUE)
data_filtered <-  data_filtered[order(data_filtered$Corr_bitter, decreasing = TRUE),]
data_filtered$Cor_order <- order(data_filtered$Corr_bitter, decreasing = TRUE)
data_filtered$RAM-c()
#data_filtered$SAW<-c()
#dfmax <- nrow(data_filtered)
for(i in 1:nrow(data_filtered)){
#  data_filtered[i,'SAW'] <-  1-((data_filtered[i,'SMD_order']/dfmax)*.5) + ((data_filtered[i,'Cor_order']/dfmax)*0.5)
  data_filtered[i,'RAM'] <-  ((data_filtered[i,'SMD_order']) + (data_filtered[i,'Cor_order']))/2
    }
data_filtered <-  data_filtered[order(data_filtered$RAM, decreasing = FALSE),]
data_filtered$RAM_order <- order(data_filtered$RAM, decreasing = FALSE)


#matches peptide list against literature reference database
colnames(lit_ref)[colnames(lit_ref) == 'Peptide.sequence'] <- 'combined_sequence'
lit_ref_seq <- lit_ref$combined_sequence
length(lit_ref_seq)
    #[1] 235
lit_ref_short <- subset(lit_ref, select = c(combined_sequence,Groups,Bitterness.Threshold.value..umol.L., Mean.Bitterness.Intensity..1.15.))

data_filtered_final <-full_join(data_filtered, lit_ref_short, by = "combined_sequence")

data_filtered_final <- data_filtered_final[,c("combined_sequence", "Positions.in.Proteins", "RAM_order", "RAM", "SMD_order", "stand_diff_mean", "Cor_order", "Corr_bitter", "pvalue", "foldchange",
                                  "Corr_age","mean_NB", "mean_B", "diff_means",
                                  "sd_NB", "sd_B", "pooled_sd", "Groups","Bitterness.Threshold.value..umol.L.", "Mean.Bitterness.Intensity..1.15.","Sequence.Length", "count", "Q_value","Theo.MHplus.in.Da","Top.Apex.RT.in.min")]

order(data_filtered$Corr_bitter, decreasing = TRUE)
#standard mean difference (SMD) >= 1.0
smd_grt_one_seq <- data_filtered$combined_sequence[data_filtered$stand_diff_mean >=1.00]
length(smd_grt_one_seq)
      #[1] 97
data_top_smd <- subset(data_filtered, combined_sequence %in% smd_grt_one_seq)


data_lit_match_seq <- data_lit_match_df$combined_sequence
data_lit_match_lr <- subset(lit_ref, combined_sequence %in% data_lit_match_seq, select = c(combined_sequence,Groups,Bitterness.Threshold.value..umol.L., Mean.Bitterness.Intensity..1.15.))
data_lit_match <-full_join(data_lit_match_df, data_lit_match_lr)

nrow(data_lit_match)
    #[1] 56

#top 54 bitterneess correlation values (MBI to peptide abundance) >= 0.505
bit_cor_top_97<- head(data_filtered[order(data_filtered$Corr_bitter, decreasing = T),], n=97)
bit_cor_top_97_seq <- bit_cor_top_97$combined_sequence
length(bit_cor_top_97_seq)
    #[1] 97

#top 54 peptides from machine learning



#Joins list together
data_bit_candidates_lit_smd <- full_join(data_top_smd, data_lit_match)
nrow(data_bit_candidates_lit_smd)
    #[1] 144

data_bit_candidates_lit_smd_cor <- full_join(data_bit_candidates_lit_smd, bit_cor_top_97)
nrow(data_bit_candidates_lit_smd_cor)
    #[1] 178

data_bit_candidates_lit_smd_cor$top_SMD[data_bit_candidates_lit_smd_cor$combined_sequence %in% smd_grt_one_seq] <- "SMD"
data_bit_candidates_lit_smd_cor$top_cor[data_bit_candidates_lit_smd_cor$combined_sequence %in% bit_cor_top_97_seq] <- "Cor"
data_bit_candidates_lit_smd_cor$lit_match[data_bit_candidates_lit_smd_cor$combined_sequence %in% lit_ref_seq] <- "Lit"
data_bit_candidates_lit_smd_cor <- data_bit_candidates_lit_smd_cor[,c("combined_sequence", "Positions.in.Proteins", "top_SMD", "top_cor", "lit_match", "stand_diff_mean","Corr_bitter","pvalue","foldchange",
                                  "Corr_age","mean_NB", "mean_B", "diff_means",
                                  "sd_NB", "sd_B", "pooled_sd", "Sequence.Length", "count", "Q_value","Theo.MHplus.in.Da","Top.Apex.RT.in.min", "Bitterness.Threshold.value..umol.L.", "Mean.Bitterness.Intensity..1.15.")]
#selects peptides that match 3 criteria
bpc_3match <- subset(data_bit_candidates_lit_smd_cor, (top_SMD == "SMD" & top_cor == "Cor" & lit_match == "Lit"), 
                     select = c(combined_sequence, Positions.in.Proteins, count, Theo.MHplus.in.Da , Q_value, stand_diff_mean, 
                                Corr_bitter, Bitterness.Threshold.value..umol.L., Mean.Bitterness.Intensity..1.15.))
# pruchased on the list YPFPGPIHN[S], YPFPGPIHNS,PFPGPIHNS
# RINKKIEKF	11.8 out of 15 on MBI scale in wrong peptie gourp

options(scipen = -1, digits = 3); View(data_filtered)
options(scipen = -1, digits = 2); View(data_bit_candidates_lit_smd_cor)


#examined data of incorectly labels top 20 pepetides
top_20<- read.csv("Top_20_List.csv", header = TRUE,)
grp_top_20=c(unlist(top_20$Sequence))
data_wrong_20<-data.frame(subset(data_filtered,data_filtered$Sequence%in%grp_top_20))
