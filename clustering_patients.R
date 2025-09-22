library(cluster)
library(dplyr)
library(pheatmap)

#predcition matrix of each SN datasets
#TAMBOOR predictions for each patient were coded as 1/-1/0 for metabolites predicted as 
#oversecreted/undersecreted/stable respectively in the patient. 
#you can use your own prediction matrix 

GSE20292 <- read.table('catMatrix_GSE20292.txt', header = TRUE)

GSE8397B <- read.table('catMatrix_GSE8397B.txt', header = TRUE)

GSE7621 <- read.table('catMatrix_GSE7621.txt', header = TRUE)

GSE20141 <- read.table('catMatrix_GSE20141.txt', header = TRUE)

GSE24378 <- read.table('catMatrix_GSE24378.txt', header = TRUE)

GSE26927 <- read.table('catMatrix_GSE26927.txt', header = TRUE)

GSE49036 <- read.table('catMatrix_GSE49036.txt', header = TRUE)

GSE114517 <- read.table('catMatrix_GSE114517.txt', header = TRUE)

colnames(GSE20292) <- paste("d1", colnames(GSE20292), sep = "_")
colnames(GSE20292) <- gsub("Cat_Matrix", "", colnames(GSE20292))

colnames(GSE8397B) <- paste("d2", colnames(GSE8397B), sep = "_")
colnames(GSE8397B) <- gsub("Cat_Matrix", "", colnames(GSE8397B))

colnames(GSE7621) <- paste("d3", colnames(GSE7621), sep = "_")
colnames(GSE7621) <- gsub("Cat_Matrix", "", colnames(GSE7621))

colnames(GSE20141) <- paste("d4", colnames(GSE20141), sep = "_")
colnames(GSE20141) <- gsub("Cat_Matrix", "", colnames(GSE20141))

colnames(GSE24378) <- paste("d5", colnames(GSE24378), sep = "_")
colnames(GSE24378) <- gsub("Cat_Matrix", "", colnames(GSE24378))

colnames(GSE26927) <- paste("d6", colnames(GSE26927), sep = "_")
colnames(GSE26927) <- gsub("Cat_Matrix", "", colnames(GSE26927))

colnames(GSE49036) <- paste("d7", colnames(GSE49036), sep = "_")
colnames(GSE49036) <- gsub("Cat_Matrix", "", colnames(GSE49036))

colnames(GSE114517) <- paste("d8", colnames(GSE114517), sep = "_")
colnames(GSE114517) <- gsub("Cat_Matrix", "", colnames(GSE114517))

mets <- read.table('mets.txt', header = FALSE) #metabolite names 
#comine each prediction matriz into a single matrix for SN datasets
xyz <- bind_cols(GSE20292,GSE8397B,GSE7621,GSE20141,GSE24378,GSE26927,GSE49036,GSE114517)

mets <- mets$V1
rownames(xyz) <- mets

#fiter metabolites that are dipeptides and tripeptides
#f_mets.txt includes filtered metabolite list
fmets <- read.table('f_mets.txt', header = FALSE)
fmets <- fmets$V1
xyz <- xyz[fmets,]

#filtered to exclude features/metabolites with more than 85% of 
#predictions being zero across 106 patients
# Count zeros in each row
zero_counts <- rowSums(xyz == 0)
# Filter rows 
filtered_data <- xyz[zero_counts <= 90, ]

#for anotate different datasets with different colors
#create annotation variable
d <- data.frame(t(filtered_data))
group <-  (substr(rownames(d), 1, 2))
annotation <- data.frame(Group = group)
rownames(annotation) <- rownames(d)

#clustering for prediction data
diss_matrix <- dist(d, method = "euclidean")
hc <- hclust(diss_matrix, method = "ward.D2")
plot(hc)
# Create and plot the heatmap
pheatmap(d, cluster_rows  = hc, show_rownames = F, show_colnames = F, cutree_rows = 3, treeheight_col = 0,annotation_row =annotation)


####################
#Feature selection 
####################

#NON-parametric Annova kruskal.test
###################################

d <- as.data.frame(t(filtered_data))
clusters <- as.factor(cutree(hc, k = 3))

#find significant genes with NON-parametric Annova
significant_features <- c()
for (feature in colnames(d)) {
  p_value <- kruskal.test(d[, feature] ~ clusters)$p.value
  if (p_value < 0.01) {
    significant_features <- c(significant_features, feature)
  }
}

#randomforest for feature selection
###################################
library(randomForest)

#rondomforest give error beacuse of some charcters so we gave names based on indexs
colnames(d)<-paste("n",seq(1:595), sep = "_")
# Set seed for reproducibility
set.seed(123)
rf <- randomForest(clusters ~ ., data = d, importance = TRUE, ntree=1000)
imp <- importance(rf)
# Sort by Mean Decrease Gini 
sorted_importance <- imp[order(imp[, "MeanDecreaseGini"], decreasing = TRUE), ]
#take their first 100 as important metabolites
imp1 <- sorted_importance[1:100,]
imp1_ind <- as.numeric(gsub("n_", "", rownames(imp1)))
d <- as.data.frame(t(filtered_data))
imp1_feat <- colnames(d)[imp1_ind]


#validation data (living brain (LB))data
#use prediction matrix of LB data
d2 <- read.table('catMatrix_Syn_sva_OrgFc.txt', header = TRUE)
mets <- read.table('mets.txt', header = FALSE)
xyz <- d2
mets <- mets$V1
rownames(xyz) <- mets
#fiter metabolites that are dipeptides and tripeptides
#f_mets.txt includes filtered metabolite list
fmets <- read.table('f_mets.txt', header = FALSE)
fmets <- fmets$V1
xyz <- xyz[fmets,]
filtered_data <- xyz
d2 <- as.data.frame(t(filtered_data))

# Select the important features of both SN and LB data
data1_selected <- d[, imp1_feat] 
data2_selected <- d2[, imp1_feat]


#Cluster assignment with KNN
#############################
library(class)


# Training data of our study: data1_selected (from SN data)
# Cluster labels : clusters of SN data
# New sample to classify: each sample of data2_selected

#Apply KNN
k <- 11  # Number of neighbors
predicted_cluster <- knn(train = data1_selected, test = data2_selected, cl = clusters, k = k)

table(predicted_cluster)

data1_clustered <- data1_selected %>%
  mutate(cluster = clusters)

data2_clustered <- data2_selected %>%
  mutate(cluster = predicted_cluster)

# Combine both datasets  for heatmap
combined_data <- rbind(data1_clustered, data2_clustered)

# Prepare data for heatmap of combined data
# Remove the cluster column for heatmap values
heatmap_data <- combined_data %>% select(-cluster)

# Create annotation for clusters to show them as different colors in heatmap
annotation <- data.frame(cluster = combined_data$cluster)
row.names(annotation) <- row.names(combined_data)
source <- c(rep("SN_data", 106), rep("LB_data", 81))
annotation <- cbind(annotation,source)

# Generate the heatmap
annotation_colors = list( source = c(SN_data='lightblue', LB_data='lightpink'),
                          cluster = c("1"="orange", "2"="turquoise","3"="lightgreen"))

#clustering for categorical data
diss_matrix2 <- dist(heatmap_data, method = "euclidean")
hc2 <- hclust(diss_matrix2, method = "ward.D2")

# Create and plot the heatmap
pheatmap(heatmap_data, cluster_rows  = hc2, show_rownames = F, show_colnames = F, cutree_rows = 3, treeheight_col = 0,
         annotation_row =annotation, annotation_colors=annotation_colors,
         border_color ="grey", cellwidth = 4, cellheight = 3)
