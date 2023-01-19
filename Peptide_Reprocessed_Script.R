#Script for processing proteome discover export data for the identification of bitter peptides in cheese.

options( scipen = 1 )
options( digits = 5 )



####Install and load required packages#######
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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



####Import raw data, creates important datafROMes(samdf = sample information, lit_ref = literature reference database, datarn = insturment data with new names) and renames columns####
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Imports literature reference database
lit_ref <- read.csv("Supplemental_Information_CSV_Database.csv", header = TRUE, stringsAsFactors = FALSE, fileEncoding = "latin1", comment.char = '#')	

#Imports peptide data list
data1 <- read.table("Kuhfeld_Bitter_Peptide_Reprocess-(1)_PeptideGroups.txt", header = TRUE, stringsAsFactors = FALSE, fileEncoding = "latin1",)

#imports file name list and extracts sample name, ID, and grouping
sample_names_df <- read.table("Kuhfeld_Bitter_Peptide_Reprocess-(1)_InputFiles.txt", header = TRUE, stringsAsFactors = FALSE) #input files line 14 and 15 column file name was manually edited from "_E_8_7_a_" & "_E_8_7_b_" to "_E_8_7a_" & "_E_8_7b_" to be consistent with the formating for the threshold samples
sample_names_ext_name <- str_extract(sample_names_df$File.Name, "([A-Z]_[0-9]_[0-9][a-z]?)") # this is a regular expression that extracts the sample name from the File.ID column. ie. "E_5.7"
sample_names_ext_rename <- str_replace(sample_names_ext_name, "(?<=[0-9])[_]", ".") # this is a regular expression that replaces the "_" between numbers in sample name with a "."
sample_names_df$sample_names <-sample_names_ext_rename #creates new column with extracted renamed names
sample_names_ext_ID <-sample_names_df$File.ID #creates list with file IDs

#imports sample overview to bring in age, Mean bitterness intensity = MBI, bitterness grouping = bitterness and catagorical grouping = grouping
samdf_import <- read.csv("Sample_Overview.csv", header = TRUE, stringsAsFactors = FALSE)

#for loop that replaces the non-detect values in MBI with LOQ (0.5)
for (i in 1:nrow(samdf_import)) {
  if (samdf_import[i,'Mean.bitterness.value'] == "ND") {
    samdf_import[i,'Mean.bitterness.value'] = 0.5} # the LOQ of 0.5 is subsituted for non-detected(ND) results to add in correlation analysis  
}

#renames columns and turns MBI to numeric
colnames(samdf_import)[colnames(samdf_import) == 'Mean.bitterness.value'] <- 'MBI'
samdf_import$MBI <- as.numeric(samdf_import$MBI)
colnames(samdf_import)[colnames(samdf_import) == 'Age'] <- 'age'
colnames(samdf_import)[colnames(samdf_import) == 'Bitterness.grouping'] <- 'grouping'
colnames(samdf_import)[colnames(samdf_import) == 'Categorical.grouping'] <- 'bitterness'
colnames(samdf_import)[colnames(samdf_import) == 'Sample.ID'] <- 'sample_names'

#creates a final version of samdf(sample data frame) which merges the sample overview and instrument input files and removes the extra columns, renames File.ID
sample_merg_df <- merge(samdf_import, sample_names_df, by.x = "sample_names", all.x = TRUE, sort = FALSE)
samdf <- subset(sample_merg_df, select = c(sample_names, age, MBI, grouping, bitterness, File.ID))
colnames(samdf)[colnames(samdf) == 'File.ID'] <- 'sample_names_ext_ID'

#shortens names of columns in data1 (instrument results) dataframe
data1_names <-colnames(data1)
data1_names <- str_replace(data1_names, "Abundances.Normalized", "ANI") #replaces abundance normalized imputed data with abrevation "ANI"
data1_names <- str_replace(data1_names, "Abundances.Origin", "ADI") #replaces abundance origin with abrevation "ADI" for detected or imputed
data1_names <- str_replace(data1_names, "Abundance", "A") #replaces abundance  with abrevation "A" for data with out normilization or imputation
data1_names <- str_replace(data1_names, ".Sample.....?.?.?.?.?.?", "") #removes sample and grouping infromation

#loop that replaces file id ie. "F1" with sample name ie. "T_0.2a"
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

#creates new dataframe datarn (insturment data re-named) and renames columns with shortend sample names
datarn <- data1
colnames(datarn) <- data1_names_new



####data cruching#### 
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#counts the number of real Mass spec detections in the samples where "Det" = detected value and "Imp" = imputed and saves as column in df
datarn$count <- NA
for (i in 1:nrow(datarn)){
  datarn[i,'count'] <- (sum(c(datarn[i,grep("^ADI.",data1_names_new)] == "Det")))
}

####Q value Function####
#uses ney 1979 formula for determing hydrophobicity or Q-value of a peptide based on free energy transfer of amino acids which are the lines below assigning amino acids there energy transfer value
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

#loop tha calculates Q-Value, 
  #creates a vector with count of each amino acid in the peptide, 
  #then multiplies this by the energy value of each peptide, 
  #then sums these values
  #then divided by the sequence length of the peptide 
  #finally saves values into dataframe as new column
Q_value_results <- c() 
for (i in 1:nrow(datarn)){
  peptide <-unlist(strsplit (x = datarn[i,"Sequence"], split = ""))
  peptide_vec <- table(peptide)
  Q_value <- sum(amino_acids[names(peptide_vec)]*peptide_vec)/length(peptide)
  Q_value_results <-rbind(Q_value_results,c(datarn[i,"Sequence"],Q_value))
}
datarn$Q_value <- as.numeric(Q_value_results[,2])

####determines standard difference mean#######
#calculates Standard mean difference (Cohen's D test) which is a measure of effect size. It is the mean divided by the standard deviation of a difference between two random values each from one of two groups
nb_sample_list <-paste0("ANI:",samdf$sample_names[samdf$bitterness == "Non-bitter"]) #first group non-bitter peptides [1] "ANI:T_0.2a" "ANI:T_0.2b" "ANI:T_0.2c" "ANI:T_0.2d" "ANI:L_3.3"  "ANI:L_2.9"  "ANI:L_6.0" 
b_sample_list <-paste0("ANI:",samdf$sample_names[samdf$bitterness == "Bitter"]) #second group bitter peptides "ANI:M_7.3"  "ANI:M_6.2"  "ANI:M_5.6"  "ANI:E_7.2"  "ANI:E_8.7a" "ANI:E_5.7"  "ANI:E_8.7b"

#calculates means of each group and finds difference
datarn$mean_NB <- rowMeans(as.matrix(datarn[nb_sample_list]), na.rm =TRUE)
datarn$mean_B <- rowMeans(as.matrix(datarn[b_sample_list]), na.rm =TRUE)
datarn$diff_means <- datarn$mean_B - datarn$mean_NB 

#calculates standard dev. of each group and number of samples in each group 
datarn$sd_NB <- rowSds(as.matrix(datarn[nb_sample_list]), na.rm =TRUE)
datarn$sd_B <- rowSds(as.matrix(datarn[b_sample_list]), na.rm =TRUE)
n_NB <- sum(samdf$bitterness == "Non-bitter")
n_B <-  sum(samdf$bitterness == "Bitter")

#determines the pooled stand dev and finally the standard mean difference
datarn$pooled_sd <- sqrt(((n_B-1)*datarn$sd_B^2 + (n_NB-1)*datarn$sd_NB^2)/ (n_NB+n_B-2))
datarn$stand_diff_mean <- datarn$diff_means / datarn$pooled_sd

####determins p value and fold change of bitter vs non-bitter group, calculates log values for later volcano ploting and finally converts values to numeric form
datarn$pvalue <- row_t_welch((datarn[nb_sample_list]), (datarn[b_sample_list]))$pvalue
datarn$pvalue <- as.numeric(datarn$pvalue, na.rm = TRUE)
datarn$foldchange <-(datarn$mean_B/datarn$mean_NB)
datarn$foldchange <- as.numeric(datarn$foldchange, na.rm = TRUE)
datarn$logfoldchange <-cbind(log2(datarn$foldchange))
datarn$logfoldchange <- as.numeric(datarn$logfoldchange, na.rm =TRUE)
datarn$logpvalue <- cbind(-log10(datarn$pvalue))
datarn$logpvalue <- as.numeric(datarn$logpvalue, na.rm =TRUE)



####converts casein to greek symbol and renames columns#####
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
colnames(datarn)[colnames(datarn) == 'Master.Protein.Accessions'] <- 'Protein' #renames column for easier reference
#references unicode for alpha, beta, and kappa greek symbols
Beta<-"\U03B2"
BetaA1<- str_c(Beta,"A1")
BetaA2<- str_c(Beta, "A2")
Alphas1<-paste0("\U03B1","s1")
Alphas2<-paste0("\U03B1","s2")
Kappa<-paste("\U03BA")

#regular expressions and for loop to deal duplicate data due to beta casein P02666 genetic varriant A1 & A2
  datarn$Positions.in.Proteins <-gsub("; P02666A2[ ][//[].....?.?.?.?.?","",as.character(datarn$Positions.in.Proteins)) # regular expresion to removes second instance of beta casein in position in protein column

#loop that labels any beta casein without overlap with the genetic varriant site 67 as just beta casein (P02666) in both protein and position in protein columns.
  for (i in 1:nrow(datarn))  {
  temp_pnp <-unlist(strsplit (x=datarn[i,'Positions.in.Proteins'], split = " "))
  if (datarn[i,'Protein'] == "P02666A1; P02666A2")
    {
    datarn[i,'Positions.in.Proteins'] <-paste("P02666",temp_pnp[2])
    datarn[i,'Protein'] <- "P02666"
    }
}

#labels position in protein
datarn$Positions.in.Proteins <-gsub("P02666",Beta,as.character(datarn$Positions.in.Proteins)) #all beta casein with out a 67 position overlap
datarn$Positions.in.Proteins <-gsub("P02666A1",BetaA1,as.character(datarn$Positions.in.Proteins)) #beta casein with P67 position overlap labled as A1
datarn$Positions.in.Proteins <-gsub("P02666A2",BetaA2,as.character(datarn$Positions.in.Proteins)) #beta casein with H67 position overlap labled as A2
datarn$Positions.in.Proteins <-gsub("P02662",Alphas1,as.character(datarn$Positions.in.Proteins))
datarn$Positions.in.Proteins <-gsub("P02663",Alphas2,as.character(datarn$Positions.in.Proteins))
datarn$Positions.in.Proteins <-gsub("P02668",Kappa,as.character(datarn$Positions.in.Proteins))

#labels protein
datarn$Protein[datarn$Protein == "P02666"] <- Beta #all beta casein with out a 67 position overlap
datarn$Protein[datarn$Protein == "P02666A1"] <- BetaA1 #beta casein with P67 position overlap labled as A1
datarn$Protein[datarn$Protein == "P02666A2"] <- BetaA2 #beta casein with H67 position overlap labled as A2
datarn$Protein[datarn$Protein == "P02662"] <- Alphas1
datarn$Protein[datarn$Protein == "P02663"] <- Alphas2
datarn$Protein[datarn$Protein == "P02668"] <- Kappa



####Filters out data based on count, molecular weight and Q-value####
  #to reduce possible bitter peptide candidates to only those peptides with enough instrument data and literature based molecular weight and Q-value selction criter
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#number of rows or orginal data set pre-filter
nrow(datarn)
      #[1] 2100

#filters for count >=4 (only includes peptides with atleast four instrement responces) creates new data frame called datamfc (re-named, filtered, count)
datarnfc <- data.frame(subset(datarn, datarn$count >= 4))
nrow(datarnfc)
      #[1] 1140

#filters for molecular weight <= 3000  creates new data frame called datafcmw (filtered, count, molecular weight)
datafcmw <- data.frame(subset(datarnfc, datarnfc$Theo.MHplus.in.Da<= 3000))
nrow(datafcmw)
      #[1] 1105

#filters for Q-Value >=1200  creates new data frame called data_filtered (filtered, count(>4), molecular weight (<=3000), Q-Value (>=1200))
data_filtered <- data.frame(subset(datafcmw, datafcmw$Q_value >= 1200))
nrow(data_filtered)
      #[1] 872

#joins sequence and modification columns and checks for duplicates
data_filtered$combined_sequence <-paste0(data_filtered$Sequence,data_filtered$Modifications) #this combines the sequence ie. "APKHKEMPFPKYP" with modification "1xOxidation [M7]" to create new column combined_sequence ie"APKHKEMPFPKYP1xOxidation [M7]"

#tests to see if there is duplicate data in data frame
n_distinct(data_filtered$combined_sequence)
      #[1] 872 
nrow(data_filtered) == n_distinct(data_filtered$combined_sequence)
      #[1] TRUE, all values are unique and there are no duplicates



####Calculates linear correlations####
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cor.results_A <- c()
cor.results_B <- c()
sample_list_ANI <-paste0("ANI.",samdf$sample_names) #list of sample names used in correlation
#loop that calculates linear regression of insturment inensity to age and bitterness, generates plots for each pepetide
pdf("Bar,box, and scatter plot of peptides abudance to cheese bitterness.pdf",height=8,width=12)
par(mfrow=c(2,2))
for (i in 1:nrow(data_filtered)) {
  data_cor <- data.frame(b=samdf$MBI,age=samdf$age,grp=factor(samdf$grouping,levels=c("Threshold", "Low", "Moderate", "Extreme")),a=c(unlist(data_filtered[i,sample_list_ANI])) ) #creates internal dataframe, groups based on bitterness
  cor_B <- cor(data_cor[,c("a","b")])[2,1] #corelates instrumental data to bitteress
  cor.results_B <-rbind(cor.results_B,c(data_filtered[i,1],cor_B)) #saves data from above
  cor_A <- cor(data_cor[,c("a","age")])[2,1] #corelates instrumental data to age
  cor.results_A <-rbind(cor.results_A,c(data_filtered[i,1],cor_A)) #saves data from above
  boxplot(a~grp,data=data_cor,main=paste(data_filtered[i,'combined_sequence'],"SMD:",data_filtered[i,'stand_diff_mean']),xlab="bitterness grouping",ylab="Relative abundance") #generate box plots based on bitterness grouping for each pepetide
  dat <- matrix(data_filtered[i,sample_list_ANI],1,14) #creates internal matrix for bar plot
  colnames(dat) <- sample_names_ext_rename 
  barplot( dat, main="Relative abundance bar plot", xlab="Sample ID" , ylab="Relative abundance", cex.names=.5) #generates bar plot showing relative abundance for each sample
  plot(samdf$age,data_filtered[i,sample_list_ANI],xlab="Cheese Age (years)",ylab="Relative abundance",main=paste("Relative abundance vs. cheese age","corr=",round(cor_A,3))) #creates scrater plot with correlation age of relative abundance to age
  lm_A <- lm(a ~ age,data=data_cor)
  abline(lm_A,)
  plot(samdf$MBI,data_filtered[i,sample_list_ANI],xlab="Mean bitterness score",ylab="Relative abundance",main=paste("Relative abundance vs. cheese bitterness","corr=",round(cor_B,3))) #creates scrater plot with correlation age of relative abundance to bitterness
  lm_B <- lm(a ~ b,data=data_cor)
  abline(lm_B)
}
dev.off()
#combines peptide, corr_age, and Corr_bitterness into one table  
data_filtered$Corr_age<-as.numeric(cor.results_A[,2])
data_filtered$Corr_bitter<-as.numeric(cor.results_B[,2])

####compares linear relationship of age and bitterness to peptide abundance####
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#age and bitterness are highly correlated as seen below
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
cor(cor_age,cor_bitter)
      #[1] 0.976
#age and bitterness are highly correlated as seen above(r = 0.975, p-value = 4e-12)



#####generates bitter peptide canadites list####
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#removes bulky extra columns
A_sample_list <-paste0("A.",samdf$sample_names)
data_filtered[A_sample_list] <-NULL#removes abundance  with abrevation "A" for data with out normilization or imputation
ADI_sample_list <-paste0("ADI.",samdf$sample_names)
data_filtered[ADI_sample_list] <-NULL#removes abundance origin with abrevation "ADI" for detected or imputed
ANI_sample_list <-paste0("ANI.",samdf$sample_names)
data_filtered[ANI_sample_list] <-NULL#removes abundance normalized imputed data with abrevation "ANI"
data_filtered["Sequence"] <-NULL
data_filtered["Protein"] <-NULL
data_filtered["Peptide.Groups.Peptide.Group.ID"] <-NULL
data_filtered["Modifications"] <-NULL
data_filtered["logpvalue"] <-NULL
data_filtered["logfoldchange"] <-NULL
data_filtered["Quan.Info"] <-NULL

#orders peptide based on rank average meam of stand diff mean and cor to bitter
data_filtered <- data_filtered[order(data_filtered$stand_diff_mean, decreasing = TRUE),] #oreders based on standard mean difference, largest value = 1
data_filtered$SMD_order <- order(data_filtered$stand_diff_mean, decreasing = TRUE) #saved the order from above to new column
data_filtered <-  data_filtered[order(data_filtered$Corr_bitter, decreasing = TRUE),] #oreders based on bitterness correlation, largest value = 1
data_filtered$Cor_order <- order(data_filtered$Corr_bitter, decreasing = TRUE) #saved the order from above to new column

#generates a new column called ROM for the rank order mean which is the average rank order of both standard mean difference and bitterness correlatio
data_filtered$ROM <-cbind((data_filtered$SMD_order + data_filtered$Cor_order)/2)

#re-orders data based on rank order mean and generates new column
data_filtered <-  data_filtered[order(data_filtered$ROM, decreasing = FALSE),]
data_filtered$ROM_order <- order(data_filtered$ROM, decreasing = FALSE)


#matches peptide list against literature reference database
colnames(lit_ref)[colnames(lit_ref) == 'Peptide.sequence'] <- 'combined_sequence' #renames so columns share names accross dataframes
nrow(lit_ref)
    #[1] 235, number of peptides in literature refernece database
lit_ref_short <- subset(lit_ref, select = c(combined_sequence,Groups,Bitterness.Threshold.value..umol.L., Mean.Bitterness.Intensity..1.15.)) #creates dataframe with only important information to merge with data_filtered

#generates a final dataframe which is the merged dataframe of filtered and short literature reference, ordered by top rand order mean first (most likely bitter peptide candidates at the top)
data_filtered_final <-left_join(data_filtered, lit_ref_short, by = "combined_sequence", copy = FALSE)
data_filtered_final <- data_filtered_final[,c("combined_sequence","Positions.in.Proteins","ROM_order","ROM","SMD_order","stand_diff_mean","Cor_order","Corr_bitter","diff_means","pvalue","mean_NB","mean_B","sd_NB","sd_B","pooled_sd",
"foldchange","Corr_age","Theo.MHplus.in.Da","Sequence.Length","Top.Apex.RT.in.min","Groups","Bitterness.Threshold.value..umol.L.","Mean.Bitterness.Intensity..1.15.")]

#generates dataframe of instrumental results with lit reference database
lit_ref_final <-subset(data_filtered_final, data_filtered_final$Groups >= 1)
nrow(lit_ref_final)
    #[1] 56, number of peptides in instrument data that are documented in the literature database



####Creates data frame for export####
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
write.csv(data_filtered_final, "data_filtered_final.csv")
write.csv(lit_ref_final, "lit_ref_final.csv")





# pruchased on the list YPFPGPIHN[S], YPFPGPIHNS,PFPGPIHNS
# RINKKIEKF	11.8 out of 15 on MBI scale in wrong peptie gourp

options(scipen = -1, digits = 2); View(data_filtered_final)


#examined data of incorectly labels top 20 pepetides
top_20<- read.csv("Top_20_List.csv", header = TRUE,)
grp_top_20=c(unlist(top_20$Sequence))
data_wrong_20<-data.frame(subset(data_filtered_final,data_filtered_final$combined_sequence%in%grp_top_20))
