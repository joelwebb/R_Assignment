getwd() #Wd was C:/Users/Joe/Desktop/R_Assignment
setwd("C:/Users/Joe/Desktop/R_Assignment") #Change wd to C:/Users/Joe/Desktop/R_Assignment
read.table("snp_position.txt", fill = TRUE) -> Snplist # create a variable titled Snplist for our data
read.table("fang_et_al_genotypes.txt", fill = TRUE) -> Fang_genotypes # create a variable titled Fang_genotypes
Snplist
Fang_genotypes
colnames(Fang_genotypes)


# file size
# # of columns
# # of lines 
attributes(Snplist)
names(Snplist)
dim(Snplist)
class(Snplist)
typeof(Snplist) # Type of object
attributes(Snplist) #additional arbitrary metadata
length(Snplist) # how many elements it contains

attributes(Fang_genotypes)
names(Fang_genotypes)
dim(Fang_genotypes)
class(Fang_genotypes)
typeof(Fang_genotypes) # Type of object
attributes(Fang_genotypes) #additional arbitrary metadata
length(Fang_genotypes) # how many elements it contains




library(data.table)
library("tibble")
# get column names
data.frame(colnames(Fang_genotypes)) -> Fang_column_Names
Fang_column_Names
Fang_new_Column_Names <- data.frame(Fang_genotypes[1,])

names(Fang_genotypes) <- Fang_new_Column_Names
Fang_genotypes

dim(Fang_column_Names)
dim(Fang_new_Column_Names)





#Add in name to all columns in fang et all genotypes

names(Fang_Genotypes)[1]<-paste("SNP_ID")
names(Fang_genotypes)[2]<-paste("Gene")
names(Fang_genotypes)[3]<-paste("Group")

snp_names <- data.frame(Fang_genotypes[,1])

#extract maize genes for ZMMIL, ZMMLR, and ZMMMR and teos genes for ZMPBAL, ZMPIL, ZMPJA


library(dplyr)
maize_genes <- filter(Fang_genotypes, Group=="ZMMIL" | Group=="ZMMLR" | Group=="ZMMR")

teos_genes <- filter(Fang_genotypes, Group=="ZMPBAL" | Group=="ZMPIL" | Group=="ZMPJA")


dim(teos_genes)
dim(maize_genes)
snp_names <- data.frame(Fang_genotypes[,1])
colnames(snp_names) 


headers <- data.frame(Transposed_Genotypes[3,])
headers
names(Transposed_Genotypes)=transposed_headers
head(Transposed_Genotypes)
t(headers) -> transposed_headers
transposed_headers
rownames(Transposed_Genotypes)=transposed_headers

Maize_genes_transposed_list <- rownames(Transposed_Genotypes=="Group"))]
Maize_genes_transposed_list

#Now we will manipulate the two files to `join` (these data sets so that we have both genotypes and positions formatted such that the first column is "SNP_ID", the second column is "Chromosome", the third column is "Position", and subsequent columns are genotype data from either maize or teosinte individuals.

Maize_genes_to_extract1 <- c('ZMMIL','ZMMLR', 'ZMMMR')

head(joined_genos)
Teos_genes_to_extract2 <- c('ZMPBA','ZMPIL','ZMPJA')

head(Snplist)
names(Snplist)[1]<-paste("SNP_ID")
names(Snplist)[3]<-paste("Chr")
names(Snplist)[4]<-paste("Position")
Snp_info <- data.frame(Snplist$SNP_ID, Snplist$Chr, Snplist$Position)
colnames(Snp_info)[1] = "Column_1"

as.data.frame(t(Fang_genotypes)) -> Transposed_Genotypes

teo_genes <- filter(Transposed_Genotypes, Group=="ZMPBA" | Group=="ZMPIL" | Group=="ZMPJA")
rownames(Transposed_Genotypes)[3]<- "Group"


Genotype_subset_teos <- subset(Transposed_Genotypes, grep(Teos_genes_to_extract2))


colnames(Transposed_Genotypes)[1] <- "Column_1"
merge(Snp_info, Transposed_Genotypes, by = "Column_1") -> joined_genos
colnames(joined_genos)[1] <- "SNP_ID"
colnames(joined_genos)[2] <- "Chromosome"
colnames(joined_genos)[3] <- 
  
  All_sample_ID_Names <- data.frame(Fang_genotypes$V1)

All_sample_ID_Names_added <- data.matrix(c("SNP_ID", "Chromosome", "BP_Position", All_sample_ID_Names))

#Extract out the genes from the joined file that we want for maize