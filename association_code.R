############################################################
# The code following is associated with the manuscript     #
# A Rosaceae Family-Level Approach to Identify Loci        # 
# Influencing Soluble Solids Content in Blackberry for     #
# DNA-Informed Breeding published in G3. The code will     #
# 1) Filter marker data, 2) Perform an interactive         #
# k-means clustering, 3) Model markers as per the two      #
# models presented in the manuscript, and 4) Output        #
# model results, PCA grouping results, and significance    #
# in accordence with the two models. Each step is descibed #
# prior to the code.                                       #
############################################################



##############################################################
# This block of code loads the marker table constructed      #
# from Python code as mydat and the SSC data for the samples #
# to be analyzed. The samples to be analyzed table is a 2    #
# table with the samples in column 1 and trait in column 2.  #
# Users should edit the allowed number of alleles per marker #
# allowed and the percent missing data per individual        #
# allowed. Defaults are allowing 2-4 alleles at a locus and  #
# 25% missing data.                                          #
##############################################################


mydat = read.table(choose.files(), header = TRUE, stringsAsFactors = FALSE)
samples_to_be_analyzed = read.table(choose.files(), header = TRUE, stringsAsFactors = FALSE)
alleles_per_marker_allowed = c(2, 3, 4)
percent_ind_missing_dat_allowed = 25

library(parallel)
cl = makeCluster(detectCores(), type = 'PSOCK')

######################################################
# Remove markers with a large abount of missing data #
######################################################

count_na <- function(x) sum(is.na(x))

missing_rate = parApply(cl = cl, mydat, MARGIN = 1, FUN = count_na)

missing_rate = (missing_rate/(length(names(mydat)) - 1))*100
missing_rate = (missing_rate < percent_ind_missing_dat_allowed)
mydat = mydat[missing_rate,]
rm(missing_rate)

#########################################
# Identify number of alleles per marker #
#########################################

allele_counter = function(x){
  alleles = unlist(strsplit(as.character(x), split = "/"))
  alleles = unique(alleles[2:length(alleles)])
  length(na.omit(alleles))
}


allele_count = parApply(cl = cl, mydat, MARGIN = 1, FUN = allele_counter)

mydat = cbind.data.frame(mydat, allele_count)

############################################################
# Filter markers by the number of observed alleles allowed #
############################################################

mydat = mydat[mydat$allele_count %in% alleles_per_marker_allowed,]
mydat = mydat[,1:(length(mydat)-1)]
rm(allele_count)

##################################################################
# This block of code construsts a presence absence allele matrix #
##################################################################

allele_matrix_contructor = function(x){
  alleles = unlist(strsplit(as.character(x), split = "/"))
  alleles = alleles[2:length(alleles)]
  alleles = unique(alleles)
  alleles = na.omit(alleles)
  marker_id2 = rep(x[1], length(alleles))
  temp_table = cbind.data.frame(marker_id2, alleles)
  temp_table[,1] = as.character(temp_table[,1])
  temp_table[,2] = as.character(temp_table[,2])
  
  allele_presence = list()
  
  count = 1
  
  while(count < (length(marker_id2)+1)){
    matrix_row = grepl(alleles[count], x)
    matrix_row = matrix_row[2:length(matrix_row)]
    matrix_row = sub("TRUE", "1", matrix_row)
    allele_presence[[count]] = sub("FALSE", "0", matrix_row)
    count = count + 1
  }
  
  allele_presence = do.call(rbind, allele_presence)
  
  temp_table = cbind.data.frame(temp_table, allele_presence)
  
  temp_table
}


allele_matrix = parApply(cl = cl, mydat, MARGIN = 1, FUN = allele_matrix_contructor)

allele_matrix = do.call(rbind, allele_matrix)
marker_id = paste(allele_matrix[,1], allele_matrix[,2], sep = "_")
allele_matrix = cbind.data.frame(marker_id, allele_matrix[,3:length(allele_matrix)])

colnames(allele_matrix) = colnames(mydat)

stopCluster(cl)

i = 2
while(i < length(colnames(allele_matrix))+1){
  allele_matrix[,i] = as.numeric(allele_matrix[,i])
  allele_matrix[,i] = sub("2", "0", allele_matrix[,i])
  allele_matrix[,i] = as.numeric(allele_matrix[,i])
  i = i +1
}
 

rownames(allele_matrix) = allele_matrix$marker_id
allele_matrix = allele_matrix[,2:length(colnames(allele_matrix))]
trans_all_mat = t(allele_matrix)


###############################################################
# This block of code performs a principal componant analysis  #
# and acompaning k-means clustering. k-means clustering is    #
# interactive requireing the input of the number of principal #
# componants to retain during analysis and the number of      #
# clusters based on a baysean information criterion (BIC).    #
# It is suggested to select the largest number of principal   #
# componants available and the number of clusters that gives  #
# the lowest BIC results of the first 4 principal componants  #
# and the clustering results will be written to files in the  #
# current working directory.                                  #
###############################################################

library(adegenet)

PCA_all <- prcomp(trans_all_mat, center = TRUE, scale. = FALSE)

nPCA <- c(1,2,3,4)

PCA <- PCA_all$x[,nPCA]
PCA = as.data.frame(PCA)

PCA_Var <- as.data.frame(cbind(seq(length(PCA_all$sdev)),(PCA_all$sdev^2)*100/sum(PCA_all$sdev^2)))
names(PCA_Var) <- c("PC", "Var")

clusters = find.clusters(trans_all_mat, center = TRUE, scale = FALSE)
Group = clusters$grp

PCA$Group = Group

write.table(PCA, file = "PCA_and_k_means_results.txt", quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

#######################################################################
# This block of code conducts an initial association analysis to      #
# identify markers to investigate further using a t-test. Users can   #
# change the significance level tested by changing the alpha variable.#
#######################################################################

samples_to_be_analyzed = samples_to_be_analyzed[order(samples_to_be_analyzed[,1]),]
SSC = samples_to_be_analyzed[,2]

clusterExport(cl, "SSC")

allele_identifier = function(x){
  library(car)
  alleles = unlist(strsplit(as.character(x), split = "/"))
  alleles = alleles[2:length(alleles)]
  alleles = unique(alleles)
  alleles = na.omit(alleles)
  marker_id2 = rep(x[1], length(alleles))
  
  temp_table = cbind.data.frame(marker_id2, alleles)
  temp_table[,1] = as.character(temp_table[,1])
  temp_table[,2] = as.character(temp_table[,2])
  
  missing_filter = is.na(x)
  geno_dat = x[!missing_filter]
  SSC_dat = SSC[!missing_filter[2:length(missing_filter)]]
  
  t_test_table = list()
  
  count = 1
  
  while(count < (length(temp_table[,1])+1)){
    allele_presence = grepl(temp_table[count, 2], geno_dat)
    allele_presence = allele_presence[2:length(allele_presence)]
    allele_presence = sub("TRUE", "Present", allele_presence)
    allele_presence = sub("FALSE", "Absent", allele_presence)
    
    if(length(grep("Present", allele_presence)) < 2){
      p_value = NA
    } else if (length(grep("Absent", allele_presence)) < 2){
      p_value = NA
    } else {
      allele_presence = as.factor(allele_presence)
      t.test_results = t.test(SSC_dat~allele_presence, alternative = "two.sided", var.equal = FALSE)
      p_value = t.test_results$p.value
    }
    
    t_test_table[[count]] = p_value
    count = count + 1
  }
  
  t_test_table = do.call(rbind, t_test_table)
  
  temp_table = cbind.data.frame(temp_table, t_test_table)
  
  temp_table
}


testing_table = parApply(cl = cl, mydat, MARGIN = 1, FUN = allele_identifier)

testing_table = do.call(rbind, testing_table)

colnames(testing_table) = c("marker_id", "allele", "t_test_p_value")
testing_table[,1] = as.character(testing_table[,1])
testing_table[,2] = as.character(testing_table[,2])
testing_table[,3] = as.character(testing_table[,3])


testing_table = na.omit(testing_table)

testing_table[,3] = as.numeric(testing_table[,3])

stopCluster(cl)

#####################################################################
# This block of code adjusts the p-values using B-H correction and  #
# filters the data base on the desired alpha level.                 #
#####################################################################

alpha = 0.05

adjustedvalues = p.adjust(testing_table$t_test_p_value, method = "BH")

testing_table = testing_table[adjustedvalues<=alpha,]

############################################################################
# The following code conducts the modeling for Allele x Group interactions #
############################################################################

testing_table[,1] = as.character(testing_table[,1])
testing_table[,2] = as.character(testing_table[,2])
testing_table[,3] = as.character(testing_table[,3])


#This block of code filters out the significant markers from the data

marker_id = mydat[,1]
mydat = mydat[,2:length(colnames(mydat))]
mydat = mydat[,order(colnames(mydat))]
mydat = cbind.data.frame(marker_id, mydat)

#This block of code establishes the factors for the model

samples_to_be_analyzed = samples_to_be_analyzed[order(samples_to_be_analyzed[,1]),]

Group = as.character(Group)

SSC_data = samples_to_be_analyzed[,2]

modeler = function(x){
  genotype = mydat[mydat[,1] %in% x[1], 2:length(colnames(mydat))]
  missing_filter = is.na(genotype)
  genotype = genotype[!missing_filter]
  genotype = grepl(x[2], genotype)
  genotype = sub("TRUE", "Present", genotype)
  genotype = sub("FALSE", "Absent", genotype)
  genotype = as.factor(genotype)
  
  SSC = SSC_data[!missing_filter]
  
  group = Group[!missing_filter]
  group = as.factor(group)
  
  aov.model = lm(SSC ~ group + genotype + group*genotype)
  results = anova(aov.model)
  
  P_group = results[1, 5]
  P_genotype = results[2, 5]
  P_group_x_genotype = results[3, 5]
  
  results = cbind.data.frame(P_group, P_genotype, P_group_x_genotype)
  
  results
}

results_table = apply(testing_table, MARGIN = 1, FUN = modeler)

results_table = do.call(rbind, results_table)

results_table = cbind.data.frame(testing_table, results_table)

############################################################################
# The following code conducts the modeling for Group nested in Allele      #
############################################################################

modeler = function(x){
  genotype = mydat[mydat[,1] %in% x[1], 2:length(colnames(mydat))]
  missing_filter = is.na(genotype)
  genotype = genotype[!missing_filter]
  genotype = grepl(x[2], genotype)
  genotype = sub("TRUE", "Present", genotype)
  genotype = sub("FALSE", "Absent", genotype)
  genotype = as.factor(genotype)
  
  SSC = SSC_data[!missing_filter]
  
  group = Group[!missing_filter]
  group = as.factor(group)
  
  aov.model = lm(SSC ~ group + group/genotype)
  results = anova(aov.model)
  
  P_group_nested_model = results[1, 5]
  P_genotype_in_group = results[2, 5]
  
  results = cbind.data.frame(P_group_nested_model, P_genotype_in_group)
  
  results
}

Group_nested_allele_results_table = apply(testing_table, MARGIN = 1, FUN = modeler)

Group_nested_allele_results_table = as.data.frame(t(do.call(rbind, results_table)))

results_table = cbind.data.frame(results_table, Group_nested_allele_results_table[,4:6])

write.table(results_table, file = "modeling_results.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

