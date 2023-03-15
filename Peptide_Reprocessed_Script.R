#Script for  identification of bitter peptides from Proteome Discover data
options( scipen = 1 )
options( digits = 5 )

####Install and load required packages##########################################
install.packages("gridExtra")
install.packages("gtools")
install.packages("Matrix")
install.packages("ggplot2")
install.packages("ggnewscale")
install.packages("dplyr")
install.packages("stringr")
install.packages("tidyverse")
install.packages("ggtext")
install.packages("ggrepel")

library(gtools)
library(dplyr)
library(matrixStats)
library(ggplot2) 
library(ggrepel)
library(tidyr)
library(Matrix)
library(gridExtra)
library(stringr)
library(ggnewscale)
library(tidyverse)
library(ggtext)
library(cowplot)

####Import raw data, creates important data frames############################## 
    # samdf = sample information,
    # lit_ref = literature reference database,
    # datarn = instrument data with new names
    # and renames columns
# #Imports literature reference database
lit_ref <- read.csv("Supplemental_Information_CSV_Database.csv", 
                    header = TRUE, 
                    stringsAsFactors = FALSE, 
                    fileEncoding = "latin1", 
                    comment.char = '#')	

#Imports peptide data list
data1 <- read.table("Kuhfeld_Bitter_Peptide_Reprocess-(1)_PeptideGroups.txt",
                    header = TRUE,
                    stringsAsFactors = FALSE,
                    fileEncoding = "latin1")

#imports file name list and extracts sample name, ID, and grouping

#input files line 14 and 15 column file name was manually edited from "_E_8_7_a_" & "_E_8_7_b_" to "_E_8_7a_" & "_E_8_7b_" to be consistent with the formatting for the threshold samples
sample_names_df <- read.table("Kuhfeld_Bitter_Peptide_Reprocess-(1)_InputFiles.txt",
                              header = TRUE,
                              stringsAsFactors = FALSE)
# this is a regular expression that extracts the sample name from the File.ID column. ie. "E_5.7"
sample_names_ext_name <- str_extract(sample_names_df$File.Name, "([A-Z]_[0-9]_[0-9][a-z]?)")
# this is a regular expression that replaces the "_" between numbers in sample name with a "."
sample_names_ext_rename <- str_replace(sample_names_ext_name, "(?<=[0-9])[_]", ".")
#creates new column with extracted renamed names
sample_names_df$sample_names <-sample_names_ext_rename
#creates list with file IDs
sample_names_ext_ID <-sample_names_df$File.ID

#imports sample overview to bring:
    # age,
    # Mean bitterness intensity = MBI, 
    # bitterness grouping = bitterness
    # categorical grouping = grouping
samdf_import <- read.csv("Sample_Overview.csv", 
                         header = TRUE, 
                         stringsAsFactors = FALSE)

#for loop that replaces the non-detect values in MBI with LOQ (0.5)
for (i in 1:nrow(samdf_import)) {
  if (samdf_import[i,'Mean.bitterness.value'] == "ND") {
    samdf_import[i,'Mean.bitterness.value'] = 0.5} # the LOQ of 0.5 is substituted for non-detected(ND) results to add in correlation analysis  
}

#renames columns and turns MBI to numeric
colnames(samdf_import)[colnames(samdf_import) == 'Mean.bitterness.value'] <- 'MBI'
samdf_import$MBI <- as.numeric(samdf_import$MBI)
colnames(samdf_import)[colnames(samdf_import) == 'Age'] <- 'age'
colnames(samdf_import)[colnames(samdf_import) == 'Bitterness.grouping'] <- 'grouping'
colnames(samdf_import)[colnames(samdf_import) == 'Categorical.grouping'] <- 'bitterness'
colnames(samdf_import)[colnames(samdf_import) == 'Sample.ID'] <- 'sample_names'

#creates a final version of samdf(sample data frame) which merges the sample overview and instrument input files and removes the extra columns, renames File.ID
sample_merg_df <- merge(samdf_import,
                        sample_names_df,
                        by.x = "sample_names",
                        all.x = TRUE, sort = FALSE)
samdf <- subset(sample_merg_df,
                select = c(sample_names, age, MBI, grouping, bitterness, File.ID))
colnames(samdf)[colnames(samdf) == 'File.ID'] <- 'sample_names_ext_ID'

#shortens names of columns in data1 (instrument results) dataframe
data1_names <-colnames(data1)
#replaces abundance normalized imputed data with abbreviation "ANI"
data1_names <- str_replace(data1_names, "Abundances.Normalized", "ANI")
#replaces abundance origin with abbreviation "ADI" for detected or imputed
data1_names <- str_replace(data1_names, "Abundances.Origin", "ADI")
#replaces abundance  with abbreviation "A" for data with out normalization or imputation
data1_names <- str_replace(data1_names, "Abundance", "A")
#removes sample and grouping information
data1_names <- str_replace(data1_names, ".Sample.....?.?.?.?.?.?", "")

#loop that replaces file id ie. "F1" with sample name ie. "T_0.2a"
data1_names_new <- list()
for (i in data1_names)  {
  temp_data1_name <-unlist(strsplit (x = i, split = "[/.]"))
  if(temp_data1_name[2] %in% sample_names_ext_ID)
    {
      temp_data2_name <- samdf$sample_names[samdf$sample_names_ext_ID==temp_data1_name[2]]
      temp_comb_name <- paste0(temp_data1_name[1],':',temp_data2_name)
      data1_names_new[[length(data1_names_new) + 1]] <- temp_comb_name
  }
  else
  {
    data1_names_new[[length(data1_names_new) + 1]] <- i
  }
}

#creates new dataframe datarn (instrument data re-named) and renames columns with shorter sample names
datarn <- data1
colnames(datarn) <- data1_names_new


####data crunching##############################################################
#counts the number of real Mass spec detection in the samples where "Det" = detected value and "Imp" = imputed and saves as column in df
datarn$count <- NA
for (i in 1:nrow(datarn)){
  datarn[i,'count'] <- (sum(c(datarn[i,grep("^ADI.",data1_names_new)] == "Det")))
}

####Q value Function####
#uses ney 1979 formula for determining hydrophobicity or Q-value of a peptide based on free energy transfer of amino acids which are the lines below assigning amino acids there energy transfer value
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

#loop that calculates Q-Value, 
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

#first group non-bitter peptides [1] "ANI:T_0.2a" "ANI:T_0.2b" "ANI:T_0.2c" "ANI:T_0.2d" "ANI:L_3.3"  "ANI:L_2.9"  "ANI:L_6.0" 
nb_sample_list <-paste0("ANI:",samdf$sample_names[samdf$bitterness == "Non-bitter"])
#second group bitter peptides "ANI:M_7.3"  "ANI:M_6.2"  "ANI:M_5.6"  "ANI:E_7.2"  "ANI:E_8.7a" "ANI:E_5.7"  "ANI:E_8.7b"
b_sample_list <-paste0("ANI:",samdf$sample_names[samdf$bitterness == "Bitter"])

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


####converts casein to Greek symbol and renames columns#########################
colnames(datarn)[colnames(datarn) == 'Master.Protein.Accessions'] <- 'Protein' #renames column for easier reference
#references Unicode for alpha, beta, and kappa Greek symbols
Beta<-"\U03B2"
BetaA1<- str_c(Beta,"A1")
BetaA2<- str_c(Beta, "A2")
Alphas1<-paste0("\U03B1","s1")
Alphas2<-paste0("\U03B1","s2")
Kappa<-paste("\U03BA")

#regular expressions and for loop to deal duplicate data due to beta casein P02666 genetic variant A1 & A2
datarn$Positions.in.Proteins <-gsub("; P02666A2[ ][//[].....?.?.?.?.?","",as.character(datarn$Positions.in.Proteins)) # regular expression to removes second instance of beta casein in position in protein column

#joins sequence and modification columns and checks for duplicates
#this combines the sequence ie. "APKHKEMPFPKYP" with modification "1xOxidation [M7]" to create new column combined_sequence ie"APKHKEMPFPKYP1xOxidation [M7]"
datarn$combined_sequence <-paste0(datarn$Sequence,datarn$Modifications)
  
#loop that labels any beta casein without overlap with the genetic variant site 67 as just beta casein (P02666) in both protein and position in protein columns.
#the loop extracts start and stop position and saves to seperate column
# Create new columns
datarn$start_position <- NA
datarn$stop_position <- NA
for (i in 1:nrow(datarn))  {
  temp_pnp <-unlist(strsplit (x=datarn[i,'Positions.in.Proteins'], split = " "))
  # Extract start and stop positions from temp_pnp[2]
  positions <- gsub("[^0-9-]+", "", temp_pnp[2])
  positions <- gsub("-", ",", positions)
  positions <- unlist(strsplit(positions, ","))
  
  # Save start and stop positions in the data frame
  datarn[i, "start_position"] <- as.numeric(positions[1])
  datarn[i, "stop_position"] <- as.numeric(positions[2])
  if (datarn[i,'Protein'] == "P02666A1; P02666A2")
  {
    datarn[i,'Positions.in.Proteins'] <-paste("P02666",temp_pnp[2])
    datarn[i,'Protein'] <- "P02666"
  }
}

#labels position in protein
#all beta casein with out a 67 position overlap
datarn$Positions.in.Proteins <-gsub("P02666",Beta,as.character(datarn$Positions.in.Proteins))
#beta casein with P67 position overlap labeled as A1
datarn$Positions.in.Proteins <-gsub("P02666A1",BetaA1,as.character(datarn$Positions.in.Proteins))
#beta casein with H67 position overlap labeled as A2
datarn$Positions.in.Proteins <-gsub("P02666A2",BetaA2,as.character(datarn$Positions.in.Proteins))
datarn$Positions.in.Proteins <-gsub("P02662",Alphas1,as.character(datarn$Positions.in.Proteins))
datarn$Positions.in.Proteins <-gsub("P02663",Alphas2,as.character(datarn$Positions.in.Proteins))
datarn$Positions.in.Proteins <-gsub("P02668",Kappa,as.character(datarn$Positions.in.Proteins))

#labels protein
#all beta casein with out a 67 position overlap
datarn$Protein[datarn$Protein == "P02666"] <- Beta
#beta casein with P67 position overlap labeled as A1
datarn$Protein[datarn$Protein == "P02666A1"] <- BetaA1
#beta casein with H67 position overlap labeled as A2
datarn$Protein[datarn$Protein == "P02666A2"] <- BetaA2
datarn$Protein[datarn$Protein == "P02662"] <- Alphas1
datarn$Protein[datarn$Protein == "P02663"] <- Alphas2
datarn$Protein[datarn$Protein == "P02668"] <- Kappa

#loop that labels the potion in protein with a * is a modification is present. This is needed later for graphs.
for (i in 1:nrow(datarn))  {
  if (datarn[i,"Modifications"] != ""){
    datarn[i,"Positions.in.Proteins"] <- paste0(datarn[i,"Positions.in.Proteins"],"*")    
  }
}

####Filters out data based on count, molecular weight and Q-value###############
#to reduce possible bitter peptide candidates to only those peptides with enough instrument data and literature based molecular weight and Q-value section criteria
#number of rows or original data set pre-filter
nrow(datarn)
      #[1] 2100

#filters for count >=4 (only includes peptides with at least four instrument responses) creates new data frame called datarnfc (re-named, filtered, count)
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


####Calculates linear correlations##############################################
cor.results_A <- c()
cor.results_B <- c()
sample_list_ANI <-paste0("ANI.",samdf$sample_names) #list of sample names used in correlation
#loop that calculates linear regression of instrument intensity to age and bitterness, generates plots for each peptide
pdf("Bar,box, and scatter plot of individual peptides.pdf",height=8,width=12)
par(mfrow=c(2,2))
for (i in 1:nrow(data_filtered)) {
  #creates internal dataframe, groups based on bitterness
  data_cor <- data.frame(b=samdf$MBI,
                         age=samdf$age,
                         grp=factor(samdf$grouping,
                                    levels=c("Threshold", "Low", "Moderate", "Extreme"))
                         ,a=c(unlist(data_filtered[i,sample_list_ANI])) )
  #correlates instrumental data to bitterness
  cor_B <- cor(data_cor[,c("a","b")])[2,1] 
  #saves data from above
  cor.results_B <-rbind(cor.results_B,c(data_filtered[i,1],cor_B)) 
  #correlates instrumental data to age
  cor_A <- cor(data_cor[,c("a","age")])[2,1] 
  #saves data from above
  cor.results_A <-rbind(cor.results_A,c(data_filtered[i,1],cor_A))
  #generate box plots based on bitterness grouping for each peptide
  boxplot(a~grp,
          data=data_cor,
          main=paste(data_filtered[i,'combined_sequence'],
                     "SMD:",
                     data_filtered[i,'stand_diff_mean']),
          xlab="bitterness grouping",
          ylab="Relative abundance") 
  
  #creates internal matrix for bar plot
  dat <- matrix(data_filtered[i,sample_list_ANI],1,14) 
  colnames(dat) <- sample_names_ext_rename 
  #generates bar plot showing relative abundance for each sample
  barplot(dat, 
          main="Relative abundance bar plot",
          xlab="Sample ID" ,
          ylab="Relative abundance",
          cex.names=.5) 
  
  #creates scatter plot with correlation age of relative abundance to age
  lm_A <- lm(a ~ age,data=data_cor)
  plot(samdf$age,data_filtered[i,sample_list_ANI],
       xlab="Cheese Age (years)",
       ylab="Relative abundance",
       main=paste("Relative abundance vs. cheese age",
                  "corr=",
                  round(cor_A,3)))
  abline(lm_A)
  
  #creates scatter plot with correlation age of relative abundance to bitterness
  lm_B <- lm(a ~ b,data=data_cor)
  plot(samdf$MBI,data_filtered[i,sample_list_ANI],
       xlab="Mean bitterness score",
       ylab="Relative abundance",
       main=paste("Relative abundance vs. cheese bitterness",
                  "corr=",round(cor_B,3)))
  abline(lm_B)
}
dev.off()

#combines peptide, corr_age, and Corr_bitterness into one table  
data_filtered$Corr_age<-as.numeric(cor.results_A[,2])
data_filtered$Corr_bitter<-as.numeric(cor.results_B[,2])


####compares linear relationship of age and bitterness to peptide abundance#####
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
#age and bitterness are  correlated as seen above(r = 0.975, p-value = 4e-12)


#####generates bitter peptide candidates list###################################
#removes bulky extra columns
A_sample_list <-paste0("A.",samdf$sample_names)
#removes abundance with abbreviation "A" for data with out normalization or imputation
data_filtered[A_sample_list] <-NULL
ADI_sample_list <-paste0("ADI.",samdf$sample_names)
#removes abundance origin with abbreviation "ADI" for detected or imputed
data_filtered[ADI_sample_list] <-NULL
ANI_sample_list <-paste0("ANI.",samdf$sample_names)
#removes unnesicary columns
data_filtered["Sequence"] <-NULL
data_filtered["Protein"] <-NULL
data_filtered["Peptide.Groups.Peptide.Group.ID"] <-NULL
data_filtered["Modifications"] <-NULL
data_filtered["Quan.Info"] <-NULL

#orders peptide based on rank average mean of stand diff mean and cor to bitter
#orders based on standard mean difference, largest value = 1
data_filtered <- data_filtered[order(data_filtered$stand_diff_mean, decreasing = TRUE),]
#saved the order from above to new column
data_filtered$SMD_order <- order(data_filtered$stand_diff_mean, decreasing = TRUE)
#orders based on bitterness correlation, largest value = 1
data_filtered <-  data_filtered[order(data_filtered$Corr_bitter, decreasing = TRUE),]
#saved the order from above to new column
data_filtered$Cor_order <- order(data_filtered$Corr_bitter, decreasing = TRUE)

#generates a new column called ROM for the rank order mean the average rank order of both standard mean difference and bitterness correlation
data_filtered$ROM <-cbind((data_filtered$SMD_order + data_filtered$Cor_order)/2)

#re-orders data based on rank order mean and generates new column
data_filtered <-  data_filtered[order(data_filtered$ROM, decreasing = FALSE),]
data_filtered$ROM_order <- order(data_filtered$ROM, decreasing = FALSE)


#matches peptide list against literature reference database
#renames so columns share names across dataframes
colnames(lit_ref)[colnames(lit_ref) == 'Peptide.sequence'] <- 'combined_sequence' 
nrow(lit_ref)
    #[1] 235, number of peptides in literature reference database

#creates dataframe with only important information to merge with data_filtered
lit_ref_short <- subset(lit_ref,
                        select = c(combined_sequence,
                                   Groups,Bitterness.Threshold.value..umol.L.,
                                   Mean.Bitterness.Intensity..1.15.,
                                   Literature.reference))

#generates a final data frame which is the merged data frame of filtered and short literature reference, 
    #Ordered by top rand order mean first (most likely bitter peptide candidates at the top)
data_filtered_final <-left_join(data_filtered,
                                lit_ref_short,
                                by = "combined_sequence", copy = FALSE)

data_filtered_final <- data_filtered_final[,c("combined_sequence",
                                              "Positions.in.Proteins",
                                              "ROM_order",
                                              "ROM",
                                              "SMD_order",
                                              "stand_diff_mean",
                                              "Cor_order",
                                              "Corr_bitter",
                                              "diff_means",
                                              "pvalue",
                                              "mean_NB",
                                              "mean_B",
                                              "sd_NB",
                                              "sd_B",
                                              "pooled_sd",
                                              "foldchange",
                                              "Corr_age",
                                              "Theo.MHplus.in.Da",
                                              "Sequence.Length",
                                              "Top.Apex.RT.in.min",
                                              "Groups",
                                              "Bitterness.Threshold.value..umol.L.",
                                              "Mean.Bitterness.Intensity..1.15.",
                                              "Literature.reference",
                                              "start_position",
                                              "stop_position")]

#generates dataframe of instrumental results with lit reference database
lit_ref_final <-subset(data_filtered_final, data_filtered_final$Groups >= 1)
nrow(lit_ref_final)
    #[1] 56, number of peptides in instrument data & the literature database


####Creates data frame for export###############################################
write.csv(data_filtered_final, "data_filtered_final.csv")
write.csv(lit_ref_final, "lit_ref_final.csv")

options(scipen = -1, digits = 2); View(data_filtered_final)


####short list of bitter peptide candidates#####################################
#selects top 12 bitter peptide candidates based on rank order mean

top_12 <- head(data_filtered, n=12)[,!(names(top_12) %in% c("Theo.MHplus.in.Da",
                                                            "Sequence.Length",
                                                            "Theo.MHplus.in.Da",
                                                            "Top.Apex.RT.in.min",
                                                            "count",
                                                            "Q_value",
                                                            "diff_means",
                                                            "sd_NB",
                                                            "sd_B",
                                                            "pooled_sd",
                                                            "stand_diff_mean",
                                                            "sd_NB",
                                                            "pvalue",
                                                            "SMD_order",
                                                            "Cor_order", 
                                                            "ROM",
                                                            "foldchange"))]


####Bitter peptide candidate categorical charts data crunching##################
#creates empty dataframe top_12_t for a box plot
columns = c("combined_sequence",
            "Positions.in.Proteins",
            "ANI",
            "Grp",
            "MBI")
top_12_t <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(top_12_t) = columns                        

counter <- 1
for (t12p in 1:nrow(top_12)){
  while (counter <=nrow(samdf)) {
    next_row <- nrow(top_12_t)+1
    top_12_t[next_row,"Positions.in.Proteins"] <- top_12[t12p,"Positions.in.Proteins"]    
    top_12_t[next_row,"combined_sequence"] <- top_12[t12p,"combined_sequence"]    
    top_12_t[next_row,"ANI"] <- top_12[t12p,paste0("ANI.",samdf[counter,'sample_names'])]
    top_12_t[next_row,"MBI"] <- samdf[counter,'MBI']
    grpID <- paste0("ANI:",samdf[counter,'sample_names'])
    #creates groups "T_L" and "M_E" for representing bitter and non-bitter groups
    if (grpID %in% nb_sample_list){
      top_12_t[next_row,"Grp"] <- "T_L"
    }
    if (grpID %in% b_sample_list){
      top_12_t[next_row,"Grp"] <- "M_E"
    }
    counter <- counter + 1
  }
  counter <- 1
}
top_12_t$Grp <- factor(top_12_t$Grp, levels=c('T_L','M_E'))


####Bitter peptide candidate box plot###########################################
#generates box plot
Box<-ggplot(data=top_12_t,
            aes(x=Positions.in.Proteins,
                y=log10(ANI),
                fill=factor(Grp)),) +
  geom_boxplot()+
  xlab(NULL)+
  xlab("Selected top 12 bitter peptide candidates") +
  ylab("Log10 (Relative abundance)") +
  labs(fill = "Categorical grouping") +
  scale_fill_manual(labels =c("Non-bitter","Bitter"),
                    values=c("#0096FF", "#FFA500")) +
  scale_y_continuous(breaks=c(4,5,6,7,8,9,10))+
#                     limits = c(2.75, 10))+
  geom_text(x=1.0,
            y=10,
            size=12,
            label="A",
            check_overlap = TRUE)+
  theme(
    legend.justification=c(1,0),
    legend.position = c(0.9975,0.01),
    legend.text = element_text(size=10),
    legend.title = element_text(size=12),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    axis.text.x = element_text(size=12,
                               colour="black",
                               angle = 90,
                               vjust = 0.5,
                               hjust=1),
    axis.text.y = element_text(size=12))+
    guides(fill = guide_legend(ncol = 2))    

png("Box_Plot.png",
    width = 900,
    height = 450)
grid.arrange(Box)
dev.off()


####Bitter peptide candidate linear regression##################################
LR <-ggplot(data=top_12_t,
            aes(x = MBI,
                y = log10(ANI))) +
  geom_point(aes(color=factor(Positions.in.Proteins)))+
  geom_smooth(method=lm,
              se=FALSE,
              aes(color=factor(Positions.in.Proteins))) +
  xlab("Selected top 12 bitter peptide candidates") +
  ylab("Log10 (Relative abundance)") +
  labs(color = "Peptides") +
  geom_text(x=0.5,
            y=10,
            size=12,
            label="A",
            check_overlap = TRUE)+
  scale_y_continuous(breaks=c(4,5,6,7,8,9,10))+
  scale_x_continuous(breaks=c(0,0.5,1,1.5,2,2.5,3,3.5))+
  theme(
    legend.justification=c(1,-0.02),
    legend.position = c(.998,0),
    legend.text = element_text(size=12),
    legend.title = element_text(size=12),
    axis.text.x = element_text(size=12))+
  guides(color = guide_legend(ncol = 5, title.position = "top"))

png("LR.png",
    width = 900,
    height = 450)
grid.arrange(LR)
dev.off()


####Bitter peptide candidate Volcano Plots Data Processing######################
#loop that creates separate protein column by stripping Positions.in.proteins
data_filtered$Protein<- NA
for (i in 1:nrow(data_filtered)){
  temp_PIP <- unlist(strsplit (x = data_filtered[i,"Positions.in.Proteins"], split = " "))
  if (temp_PIP[1] %in% protein_list){
    data_filtered[i,"Protein"] <- as.character(temp_PIP[1])
  }
  else {
    data_filtered[i,"Protein"] <- as.character("other")
  }
}

#creates list of casein proteins 
protein_list <- c("κ",
                  "αs1","αs2",
                   "β","βA1","βA2")

#adds protein column to dataframe and replaces all non-casein proteins with "other"
for (i in 1:nrow(data_filtered)){
  temp_PIP <- unlist(strsplit (x = data_filtered[i,"Positions.in.Proteins"], split = " "))
  if (temp_PIP[1] %in% protein_list){
    data_filtered[i,"Protein"] <- as.character(temp_PIP[1])
  }
  else {
    data_filtered[i,"Protein"] <- as.character("other")
  }
}

#creates data filter protein as factor for legend allignment
data_filtered$Protein<-factor(data_filtered$Protein, levels = c("κ", "αs1","αs2", "β","βA1","βA2","other"))


####Bitter peptide candidate fold change precentage calculations################

#list of top 12 peptides for manual reference 
grp_top_12 <- unique(top_12_t$Positions.in.Proteins)

#group the data_filtered data frame by protein
grouped_data <- group_by(data_filtered, Protein)

#calculate the total number of rows in the data_filtered data frame
total_rows <- nrow(data_filtered)

#calculate the total count of positive and negative logfoldchange for each unique protein
counts_pos_FC<- summarise(grouped_data, count_pos = sum(logfoldchange  > 0))
counts_neg_FC<- summarise(grouped_data, count_neg = sum(logfoldchange  < 0))
#add a column with the percentage of total number of rows
counts_pos_FC$percentage_pos<-NA
counts_neg_FC$percentage_neg<-NA
counts_neg_FC$pos_ratio<-NA
for (i in 1:nrow(counts_pos_FC)){
counts_pos_FC[i,"percentage_pos"] <- round(((counts_pos_FC[i,"count_pos"] / (counts_pos_FC[i,"count_pos"]+counts_neg_FC[i,"count_neg"])) * 100),1)
counts_neg_FC[i,"percentage_neg"] <- round(((counts_neg_FC[i,"count_neg"] / (counts_pos_FC[i,"count_pos"]+counts_neg_FC[i,"count_neg"])) * 100),1)
counts_neg_FC[i,"pos_ratio"] <- round(counts_pos_FC[i,"percentage_pos"] / counts_neg_FC[i,"percentage_neg"],1)
}
merged_protein_counts <- merge(counts_pos_FC, counts_neg_FC, by = "Protein")
      # Protein count_pos percentage_pos count_neg percentage_neg pos_ratio
      # 1   other         9           34.6        17           65.4       0.5
      # 2     αs1        21           12.7       145           87.3       0.1
      # 3     αs2        55           49.1        57           50.9       1.0
      # 4       β       220           51.6       206           48.4       1.1
      # 5     βA1        47           81.0        11           19.0       4.3
      # 6     βA2        37           75.5        12           24.5       3.1
      # 7       κ        20           57.1        15           42.9       1.3


####################Separated Volcano Plots#####################################
for (counter in 1:6){
Pro_highlighted <- protein_list[counter]
VCP <- ggplot(data=data_filtered,
       aes(x=logfoldchange,
           y=logpvalue)) +
  theme_bw()+
  xlab(NULL)+
  ylab(NULL)+
  xlim(-7.5, 7.5) +
  ylim(-0.1,4) +
  geom_hline(yintercept = 1.3,
             linetype='dotted',
             col = 'red')+
  geom_richtext(x=-8,
            y=4,
            size=4,
            label=paste0(abc_list[counter],") ",protein_list[counter],"-casein"),
            vjust = "inward", 
            hjust = "inward",
            fontface=2)+
  geom_point(aes(color=factor(Protein)),
             size = ifelse(1:nrow(data_filtered) <= 12, 3.5, 1.5),
             alpha = ifelse(data_filtered$Protein == Pro_highlighted, 1, 0))+
   scale_color_manual(
     values=c(β="#D82909",
              βA1="#f06625",
               βA2="yellow2",
               αs1="#baccdd",
               αs2="#517fab",
               κ="#0d0404",
               other="#7F7F7F"))+
  theme(legend.position="none")+
  geom_text_repel(
    data = subset(data_filtered,
                  data_filtered$Positions.in.Proteins %in% grp_top_12 & data_filtered$Protein %in% Pro_highlighted),
    aes(label = paste0(ROM_order,") ",Positions.in.Proteins)),
#    vjust = "outward", 
#    hjust = "inward",           
    direction = "both",
    color = "black",
    nudge_x =.1,
    nudge_y = .1,
    size = 4,
    box.padding = unit(.5, "lines"),
    point.padding = unit(.5, "lines"))

assign(str_c("VCP_", protein_list[counter]),VCP)
print(str_c("VCP_", protein_list[counter]))
}

#creates plot for extracting legend from of fold change for combined volcano plot
VCP_legend<-ggplot(data = subset(data_filtered, data_filtered$Protein %in% protein_list),
            aes(x=logfoldchange,
                y=logpvalue))+
  labs(color = "Fold change ratio") +
  theme(legend.text = element_text(size=12),
    legend.title = element_text(size=12, face=2))+
  ####################protein legend
  geom_point(aes(colour = factor(Protein)))+
  scale_color_manual(
    values=c(κ="#0d0404",
             αs1="#baccdd",
             αs2="#517fab",
             β="#950a11",
             βA1="#f06625",
             βA2="yellow2"),
    labels=c("κ = 1.3",
             "αs1 = 0.1",
             "αs2 = 1.0",
             "β = 1.1",
             "βA1 = 4.3",
             "βA2 = 1.3"))
    
# create a combined legend extracted from VCP_legend
legend <- get_legend(VCP_legend + theme(legend.position = "right"))

# create a combined plot of the six volcano plots
combined_VCP<-grid.arrange(VCP_κ,VCP_αs1,VCP_αs2,VCP_β,VCP_βA1,VCP_βA2, ncol=3, nrow=2, left=paste("-Log(p-value)"), bottom="Log2 Fold Change (relative abundance of threshold & low vs. moderate & extreme bitterness sample grouping)")

png("Volcano_Plot_Combined.png",
    width = 1000,
    height = 500)
grid.arrange(combined_VCP, right=legend)
dev.off()

  
####################Combined Volcano Plot#######################################
VCP<-ggplot(data=data_filtered,
       aes(x=logfoldchange,
           y=logpvalue,
           label = Positions.in.Proteins)) +
  theme_bw()+
  xlab("Log2 Fold Change (relative abundance of threshold & low vs. moderate & extreme bitterness sample grouping)") +
  ylab("-Log10(p-value)") +
  xlim(-12.5, 12.5) +
  labs(color = "Peptide's origin") +
  ylim(-0.1,4) +
  geom_point(aes(colour=Protein))+
  geom_hline(yintercept = 1.3,
             linetype='dotted',
             col = 'red')+
  theme(
    legend.justification=c(1,-0.02),
    legend.position = c(.998,0),
    legend.text = element_text(size=12),
    legend.title = element_text(size=12),
    axis.text.x = element_text(size=12),
    axis.title.y = element_text(size=12))+
  ####################protein legend
  geom_point(aes(colour = Protein),
             size = ifelse(1:nrow(data_filtered) <= 12, 3.5, 1.5))+
  scale_color_manual(
    values=c(β="#950a11",
             βA1="#f06625",
             βA2="yellow2",
             αs1="#baccdd",
             αs2="#517fab", ##517fab
             κ="#0d0404",
             other="#7F7F7F"))+
  #################Top_12 legend
  ggnewscale::new_scale_colour() + 
  labs(color = "Selected peptides") +
  geom_point(aes(colour = paste0(ROM_order,
                                 ") ",
                                 Positions.in.Proteins)),
             data = subset(data_filtered, 
                           data_filtered$Positions.in.Proteins %in% grp_top_12))+
  scale_color_manual(
    labels=c("1) β [60-65]",
             "2) κ [97-103]",
             "3) αs1 [180-187]",
             "4) βA2 [60-68]",
             "5) β [73-79]",
             "6) β [198-205]",
             "7) β [165-189]",
             "8) β [111-116]",
             "9) β [145-156]",
             "10) αs1 [181-190]",
             "11) βA1 [60-69]*",
             "12) βA1 [60-69]"),
    values=c("1) β [60-65]"="#950a11",
             "2) κ [97-103]"="#0d0404",
             "3) αs1 [180-187]"="#baccdd",
             "4) βA2 [60-68]"="yellow2",
             "5) β [73-79]"="#950a11",
             "6) β [198-205]"="#950a11",
             "7) β [165-189]"="#950a11",
             "8) β [111-116]"="#950a11",
             "9) β [145-156]"="#950a11",
             "10) αs1 [181-190]"="#baccdd",
             "11) βA1 [60-69]*"="#f06625",
             "12) βA1 [60-69]"="#f06625"))+
  theme(legend.key=element_rect(fill="white",
                                color = "white"))+
  guides(colour=guide_legend(override.aes=list(color="white"),
                             label.hjust	= 0,
                             label.position = "left"))+
  geom_text_repel(
    data = subset(data_filtered,
                  data_filtered$Positions.in.Proteins %in% grp_top_12),
    aes(label = ROM_order),
    direction = "both",
    nudge_x =.1,
    nudge_y = .1,
    size = 4,
    box.padding = unit(.5, "lines"),
    point.padding = unit(.5, "lines"))

png("Volcano_Plot.png",
    width = 900,
    height = 450)
grid.arrange(VCP)
dev.off()

########################beta heat map abundance#################################
#subests all beta-casein into seperate dataframe
data_filtered_beta <- data_filtered[grepl("\\b(β|βA1|βA2)\\b", data_filtered$Positions.in.Proteins), ]

#betaA1 fasta seq to list
beta_seq <- "RELEELNVPGEIVESLSSSEESITRINKKIEKFQSEEQQQTEDELQDKIHPFAQTQSLVYPFPGPIHNSLPQNIPPLTQTPVVVPPFLQPEVMGVSKVKEAMAPKHKEMPFPKYPVEPFTESQSLTLTDVENLHLPLPLLQSWMHQPHQPLPPTVMFPPQSVLSLSQSKVLPVPQKAVPYPQRDMPIQAFLLYQEPVLGPVRGPFPIIV"
beta_fasta_seq <- strsplit(beta_seq, "")[[1]]
length(beta_fasta_seq)
#[1] 209

# select columns that start with "ANI" and contain E,M,T,L
ani_E_cols <- grep("^ANI.*E.*", colnames(data_filtered_beta), value=TRUE)
ani_M_cols <- grep("^ANI.*M.*", colnames(data_filtered_beta), value=TRUE)
ani_L_cols <- grep("^ANI.*L.*", colnames(data_filtered_beta), value=TRUE)
ani_T_cols <- grep("^ANI.*T.*", colnames(data_filtered_beta), value=TRUE)

# calculate rowMeans for selected columns
beta_E_means <- rowMeans(data_filtered_beta[, ani_E_cols])
beta_M_means <- rowMeans(data_filtered_beta[, ani_M_cols])
beta_L_means <- rowMeans(data_filtered_beta[, ani_L_cols])
beta_T_means <- rowMeans(data_filtered_beta[, ani_T_cols])

#merges averaged E,T,L,M with beta subset df 
data_filtered_beta <- cbind(data_filtered_beta,beta_E_means,beta_M_means,beta_T_means,beta_L_means)

# create an empty multidimensional dataframe with dimensions defined by start_position, stop_position and number of rows in data_filtered
beta_hm_E<- array(NA, dim=c(length(beta_fasta_seq), nrow(data_filtered_beta)))
beta_hm_M<- array(NA, dim=c(length(beta_fasta_seq), nrow(data_filtered_beta)))
beta_hm_T<- array(NA, dim=c(length(beta_fasta_seq), nrow(data_filtered_beta)))
beta_hm_L<- array(NA, dim=c(length(beta_fasta_seq), nrow(data_filtered_beta)))
beta_hm_smd<- array(NA, dim=c(length(beta_fasta_seq), nrow(data_filtered_beta)))

# loop through each row in data_filtered
for(i in 1:nrow(data_filtered_beta)) {
  # get the relevant start_position, stop_position and mean_B values
  start_pos <- data_filtered_beta$start_position[i]
  stop_pos <- data_filtered_beta$stop_position[i]
  mean_E_val <- data_filtered_beta$beta_E_means[i]
  mean_T_val <- data_filtered_beta$beta_T_means[i]
  mean_M_val <- data_filtered_beta$beta_M_means[i]
  mean_L_val <- data_filtered_beta$beta_L_means[i]
  smd_val <- data_filtered_beta$stand_diff_mean[i]
  
  # assign mean_hm_val to the relevant positions in the beta_hm and rename the column
  # Create a new column in beta_hm for the amino acid
  beta_hm_E[start_pos:stop_pos,i] <- mean_E_val 
  beta_hm_T[start_pos:stop_pos,i] <- mean_T_val 
  beta_hm_M[start_pos:stop_pos,i] <- mean_M_val 
  beta_hm_L[start_pos:stop_pos,i] <- mean_L_val 
  beta_hm_smd[start_pos:stop_pos,i] <- smd_val 
  
}
#extracts peptides names as list
beta_pep_name <- as.list(c("index",
                           "AAs",
                           "count",
                           "avg_abs_E",
                           "avg_abs_M",
                           "avg_abs_L",
                           "avg_abs_T",
                           "avg_smd"))

# count non-NA values in each row
row_counts<- apply(beta_hm_E, 1, function(x) sum(!is.na(x)))

# compute row means, ignoring NA values
row_means_E <- rowMeans(beta_hm_E, na.rm = TRUE)
row_means_T <- rowMeans(beta_hm_T, na.rm = TRUE)
row_means_M <- rowMeans(beta_hm_M, na.rm = TRUE)
row_means_L <- rowMeans(beta_hm_L, na.rm = TRUE)
row_means_smd <-rowMeans(beta_hm_smd, na.rm = TRUE)

#generate index
index <- 1:length(beta_fasta)

#created merged fil for export
merged_mh <-cbind(index,
                  beta_fasta_seq,
                  row_counts,
                  row_means_E,
                  row_means_M,
                  row_means_L,
                  row_means_T,
                  row_means_smd) 
colnames(merged_mh) <- beta_pep_name

merged_mh <-t(merged_mh)
write.csv(merged_mh, "merged_dataframe_heatmap_beta_casein.csv")

########################as1 heat map abundance##################################
#subests all as1-casein into seperate dataframe
data_filtered_as1 <- data_filtered[grepl("\\b(αs1)\\b", data_filtered$Positions.in.Proteins), ]

#as1 fasta seq to list
as1_seq <- "RPKHPIKHQGLPQEVLNENLLRFFVAPFPEVFGKEKVNELSKDIGSESTEDQAMEDIKQMEAESISSSEEIVPNSVEQKHIQKEDVPSERYLGYLEQLLRLKKYKVPQLEIVPNSAEERLHSMKEGIHAQQKEPMIGVNQELAYFYPELFRQFYQLDAYPSGAWYYVPLGTQYTDAPSFSDIPNPIGSENSEKTTMPLW"
as1_fasta_seq <- strsplit(as1_seq, "")[[1]]
length(as1_fasta_seq)
    #[1] 199

# select columns that start with "ANI" and contain E,M,T,L
ani_E_cols <- grep("^ANI.*E.*", colnames(data_filtered_as1), value=TRUE)
ani_M_cols <- grep("^ANI.*M.*", colnames(data_filtered_as1), value=TRUE)
ani_L_cols <- grep("^ANI.*L.*", colnames(data_filtered_as1), value=TRUE)
ani_T_cols <- grep("^ANI.*T.*", colnames(data_filtered_as1), value=TRUE)

# calculate rowMeans for selected columns
as1_E_means <- rowMeans(data_filtered_as1[, ani_E_cols])
as1_M_means <- rowMeans(data_filtered_as1[, ani_M_cols])
as1_L_means <- rowMeans(data_filtered_as1[, ani_L_cols])
as1_T_means <- rowMeans(data_filtered_as1[, ani_T_cols])

#merges averaged E,T,L,M with as1 subset df 
data_filtered_as1 <- cbind(data_filtered_as1,as1_E_means,as1_M_means,as1_T_means,as1_L_means)

# create an empty multidimensional dataframe with dimensions defined by start_position, stop_position and number of rows in data_filtered
as1_hm_E<- array(NA, dim=c(length(as1_fasta_seq), nrow(data_filtered_as1)))
as1_hm_M<- array(NA, dim=c(length(as1_fasta_seq), nrow(data_filtered_as1)))
as1_hm_T<- array(NA, dim=c(length(as1_fasta_seq), nrow(data_filtered_as1)))
as1_hm_L<- array(NA, dim=c(length(as1_fasta_seq), nrow(data_filtered_as1)))
as1_hm_smd<- array(NA, dim=c(length(as1_fasta_seq), nrow(data_filtered_as1)))

# loop through each row in data_filtered
for(i in 1:nrow(data_filtered_as1)) {
  # get the relevant start_position, stop_position and mean_B values
  start_pos <- data_filtered_as1$start_position[i]
  stop_pos <- data_filtered_as1$stop_position[i]
  mean_E_val <- data_filtered_as1$as1_E_means[i]
  mean_T_val <- data_filtered_as1$as1_T_means[i]
  mean_M_val <- data_filtered_as1$as1_M_means[i]
  mean_L_val <- data_filtered_as1$as1_L_means[i]
  smd_val <- data_filtered_as1$stand_diff_mean[i]
  
  # assign mean_hm_val to the relevant positions in the as1_hm and rename the column
  # Create a new column in as1_hm for the amino acid
  as1_hm_E[start_pos:stop_pos,i] <- mean_E_val 
  as1_hm_T[start_pos:stop_pos,i] <- mean_T_val 
  as1_hm_M[start_pos:stop_pos,i] <- mean_M_val 
  as1_hm_L[start_pos:stop_pos,i] <- mean_L_val 
  as1_hm_smd[start_pos:stop_pos,i] <- smd_val 
  
}
#extracts peptides names as list
as1_pep_name <- as.list(c("index",
                           "AAs",
                           "count",
                           "avg_abs_E",
                           "avg_abs_M",
                           "avg_abs_L",
                           "avg_abs_T",
                           "avg_smd"))

# count non-NA values in each row
row_counts<- apply(as1_hm_E, 1, function(x) sum(!is.na(x)))

# compute row means, ignoring NA values
row_means_E <- rowMeans(as1_hm_E, na.rm = TRUE)
row_means_T <- rowMeans(as1_hm_T, na.rm = TRUE)
row_means_M <- rowMeans(as1_hm_M, na.rm = TRUE)
row_means_L <- rowMeans(as1_hm_L, na.rm = TRUE)
row_means_smd <-rowMeans(as1_hm_smd, na.rm = TRUE)

#generate index
index <- 1:length(as1_fasta_seq)

#created merged fil for export
merged_mh <-cbind(index,
                  as1_fasta_seq,
                  row_counts,
                  row_means_E,
                  row_means_M,
                  row_means_L,
                  row_means_T,
                  row_means_smd) 
colnames(merged_mh) <- as1_pep_name

merged_mh <-t(merged_mh)
write.csv(merged_mh, "merged_dataframe_heatmap_as1_casein.csv")


########################as2 heat map abundance##################################
#subests all as2-casein into seperate dataframe
data_filtered_as2 <- data_filtered[grepl("\\b(αs2)\\b", data_filtered$Positions.in.Proteins), ]

#as2 fasta seq to list
as2_seq <- "KNTMEHVSSSEESIISQETYKQEKNMAINPSKENLCSTFCKEVVRNANEEEYSIGSSSEESAEVATEEVKITVDDKHYQKALNEINQFYQKFPQYLQYLYQGPIVLNPWDQVKRNAVPITPTLNREQLSTSEENSKKTVDMESTEVFTKKTKLTEEEKNRLNFLKKISQRYQKFALPQYLKTVYQHQKAMKPWIQPKTKVIPYVRYL"
as2_fasta_seq <- strsplit(as2_seq, "")[[1]]
length(as2_fasta_seq)
      #[1] 207

# select columns that start with "ANI" and contain E,M,T,L
ani_E_cols <- grep("^ANI.*E.*", colnames(data_filtered_as2), value=TRUE)
ani_M_cols <- grep("^ANI.*M.*", colnames(data_filtered_as2), value=TRUE)
ani_L_cols <- grep("^ANI.*L.*", colnames(data_filtered_as2), value=TRUE)
ani_T_cols <- grep("^ANI.*T.*", colnames(data_filtered_as2), value=TRUE)

# calculate rowMeans for selected columns
as2_E_means <- rowMeans(data_filtered_as2[, ani_E_cols])
as2_M_means <- rowMeans(data_filtered_as2[, ani_M_cols])
as2_L_means <- rowMeans(data_filtered_as2[, ani_L_cols])
as2_T_means <- rowMeans(data_filtered_as2[, ani_T_cols])

#merges averaged E,T,L,M with as2 subset df 
data_filtered_as2 <- cbind(data_filtered_as2,as2_E_means,as2_M_means,as2_T_means,as2_L_means)


# create an empty multidimensional dataframe with dimensions defined by start_position, stop_position and number of rows in data_filtered
as2_hm_E<- array(NA, dim=c(length(as2_fasta_seq), nrow(data_filtered_as2)))
as2_hm_M<- array(NA, dim=c(length(as2_fasta_seq), nrow(data_filtered_as2)))
as2_hm_T<- array(NA, dim=c(length(as2_fasta_seq), nrow(data_filtered_as2)))
as2_hm_L<- array(NA, dim=c(length(as2_fasta_seq), nrow(data_filtered_as2)))
as2_hm_smd<- array(NA, dim=c(length(as2_fasta_seq), nrow(data_filtered_as2)))

# loop through each row in data_filtered
for(i in 1:nrow(data_filtered_as2)) {
  # get the relevant start_position, stop_position and mean_B values
  start_pos <- data_filtered_as2$start_position[i]
  stop_pos <- data_filtered_as2$stop_position[i]
  mean_E_val <- data_filtered_as2$as2_E_means[i]
  mean_T_val <- data_filtered_as2$as2_T_means[i]
  mean_M_val <- data_filtered_as2$as2_M_means[i]
  mean_L_val <- data_filtered_as2$as2_L_means[i]
  smd_val <- data_filtered_as2$stand_diff_mean[i]
  
  # assign mean_hm_val to the relevant positions in the as2_hm and rename the column
  # Create a new column in as2_hm for the amino acid
  as2_hm_E[start_pos:stop_pos,i] <- mean_E_val 
  as2_hm_T[start_pos:stop_pos,i] <- mean_T_val 
  as2_hm_M[start_pos:stop_pos,i] <- mean_M_val 
  as2_hm_L[start_pos:stop_pos,i] <- mean_L_val 
  as2_hm_smd[start_pos:stop_pos,i] <- smd_val 
  
}
#extracts peptides names as list
as2_pep_name <- as.list(c("index",
                          "AAs",
                          "count",
                          "avg_abs_E",
                          "avg_abs_M",
                          "avg_abs_L",
                          "avg_abs_T",
                          "avg_smd"))

# count non-NA values in each row
row_counts<- apply(as2_hm_E, 1, function(x) sum(!is.na(x)))

# compute row means, ignoring NA values
row_means_E <- rowMeans(as2_hm_E, na.rm = TRUE)
row_means_T <- rowMeans(as2_hm_T, na.rm = TRUE)
row_means_M <- rowMeans(as2_hm_M, na.rm = TRUE)
row_means_L <- rowMeans(as2_hm_L, na.rm = TRUE)
row_means_smd <-rowMeans(as2_hm_smd, na.rm = TRUE)

#generate index
index <- 1:length(as2_fasta_seq)

#created merged fil for export
merged_mh <-cbind(index,
                  as2_fasta_seq,
                  row_counts,
                  row_means_E,
                  row_means_M,
                  row_means_L,
                  row_means_T,
                  row_means_smd) 
colnames(merged_mh) <- as2_pep_name

merged_mh <-t(merged_mh)
write.csv(merged_mh, "merged_dataframe_heatmap_as2_casein.csv")


########################kapa heat map abundance#################################
#subests all kapa-casein into seperate dataframe
data_filtered_kapa <- data_filtered[grepl("\\b(κ)\\b", data_filtered$Positions.in.Proteins), ]

#kapaA1 fasta seq to list
kapa_seq <- "QEQNQEQPIRCEKDERFFSDKIAKYIPIQYVLSRYPSYGLNYYQQKPVALINNQFLPYPYYAKPAAVRSPAQILQWQVLSNTVPAKSCQAQPTTMARHPHPHLSFMAIPPKKNQDKTEIPTINTIASGEPTSTPTTEAVESTVATLEDSPEVIESPPEINTVQVTSTAV"
kapa_fasta_seq <- strsplit(kapa_seq, "")[[1]]
length(kapa_fasta_seq)
      #[1] 169

# select columns that start with "ANI" and contain E,M,T,L
ani_E_cols <- grep("^ANI.*E.*", colnames(data_filtered_kapa), value=TRUE)
ani_M_cols <- grep("^ANI.*M.*", colnames(data_filtered_kapa), value=TRUE)
ani_L_cols <- grep("^ANI.*L.*", colnames(data_filtered_kapa), value=TRUE)
ani_T_cols <- grep("^ANI.*T.*", colnames(data_filtered_kapa), value=TRUE)

# calculate rowMeans for selected columns
kapa_E_means <- rowMeans(data_filtered_kapa[, ani_E_cols])
kapa_M_means <- rowMeans(data_filtered_kapa[, ani_M_cols])
kapa_L_means <- rowMeans(data_filtered_kapa[, ani_L_cols])
kapa_T_means <- rowMeans(data_filtered_kapa[, ani_T_cols])

#merges averaged E,T,L,M with kapa subset df 
data_filtered_kapa <- cbind(data_filtered_kapa,kapa_E_means,kapa_M_means,kapa_T_means,kapa_L_means)


# create an empty multidimensional dataframe with dimensions defined by start_position, stop_position and number of rows in data_filtered
kapa_hm_E<- array(NA, dim=c(length(kapa_fasta_seq), nrow(data_filtered_kapa)))
kapa_hm_M<- array(NA, dim=c(length(kapa_fasta_seq), nrow(data_filtered_kapa)))
kapa_hm_T<- array(NA, dim=c(length(kapa_fasta_seq), nrow(data_filtered_kapa)))
kapa_hm_L<- array(NA, dim=c(length(kapa_fasta_seq), nrow(data_filtered_kapa)))
kapa_hm_smd<- array(NA, dim=c(length(kapa_fasta_seq), nrow(data_filtered_kapa)))

# loop through each row in data_filtered
for(i in 1:nrow(data_filtered_kapa)) {
  # get the relevant start_position, stop_position and mean_B values
  start_pos <- data_filtered_kapa$start_position[i]
  stop_pos <- data_filtered_kapa$stop_position[i]
  mean_E_val <- data_filtered_kapa$kapa_E_means[i]
  mean_T_val <- data_filtered_kapa$kapa_T_means[i]
  mean_M_val <- data_filtered_kapa$kapa_M_means[i]
  mean_L_val <- data_filtered_kapa$kapa_L_means[i]
  smd_val <- data_filtered_kapa$stand_diff_mean[i]
  
  # assign mean_hm_val to the relevant positions in the kapa_hm and rename the column
  # Create a new column in kapa_hm for the amino acid
  kapa_hm_E[start_pos:stop_pos,i] <- mean_E_val 
  kapa_hm_T[start_pos:stop_pos,i] <- mean_T_val 
  kapa_hm_M[start_pos:stop_pos,i] <- mean_M_val 
  kapa_hm_L[start_pos:stop_pos,i] <- mean_L_val 
  kapa_hm_smd[start_pos:stop_pos,i] <- smd_val 
  
}
#extracts peptides names as list
kapa_pep_name <- as.list(c("index",
                          "AAs",
                          "count",
                          "avg_abs_E",
                          "avg_abs_M",
                          "avg_abs_L",
                          "avg_abs_T",
                          "avg_smd"))

# count non-NA values in each row
row_counts<- apply(kapa_hm_E, 1, function(x) sum(!is.na(x)))

# compute row means, ignoring NA values
row_means_E <- rowMeans(kapa_hm_E, na.rm = TRUE)
row_means_T <- rowMeans(kapa_hm_T, na.rm = TRUE)
row_means_M <- rowMeans(kapa_hm_M, na.rm = TRUE)
row_means_L <- rowMeans(kapa_hm_L, na.rm = TRUE)
row_means_smd <-rowMeans(kapa_hm_smd, na.rm = TRUE)

#generate index
index <- 1:length(kapa_fasta_seq)

#created merged fil for export
merged_mh <-cbind(index,
                  kapa_fasta_seq,
                  row_counts,
                  row_means_E,
                  row_means_M,
                  row_means_L,
                  row_means_T,
                  row_means_smd) 
colnames(merged_mh) <- kapa_pep_name

merged_mh <-t(merged_mh)
write.csv(merged_mh, "merged_dataframe_heatmap_kapa_casein.csv")