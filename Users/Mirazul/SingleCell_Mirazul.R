#####Md. Mirazul Islam
##03-17-2019
#Single_Cell_Class
##Final Project
############################
##STEP1: Data Download######
############################

#Data download link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056
##Downlod txt file
#Unzip the file.
###############

############################
##STEP2: Data Processing######
############################
data<-read.table("GSE72056_melanoma_single_cell_revised_v2.txt", header = T) #reading the data
#####Processing the matrix
cell.name<-colnames(data) ## Each single cell is marked by unique markers in column name.
##Different cy name represent different patient's IDs.
gene.name<-data[4:nrow(data),1] # All the gene names
data<- data[-c(1,3),] ##Seting up filter based on malignant status
data<-rbind(data, cell.name) 
data<-data[,-1]
colnames(data)<-data[1,] # make the column name as malignant values
data<-data[,colnames(data)=="2"] # Short columns based on column names, 2=malignant, 1=non-malignant, 0=unknown. 
#For Non-Malignant cells (Fig1, panel D), only need to change the value to maligname value to 1 duting column filtering.
colnames(data)<- data[nrow(data),]
data<-data[-nrow(data),]
data<- data[-1,]
data<- cbind(gene.name,data) # Final datamatrix with malignant cells only from 13 patients. 

#This is the final processed data matrix from 13 patients. 6 paritents' id is known, 7 patients' id
#is unknown. Need to upload in google drive. 


###############################################
##STEP3: Print datamatrix for 13 patients######
###############################################
write.table(data, file = "data.malig.txt", sep = "\t", row.names = T,col.names = T)

##Now need to process the data for only those 6 patients' shown in Fig1, panel c. 
#Let's identify those 6 patient's ID.

###############################################
##STEP4: Identify Cluster Info for 6 and then all (13) patients#
###############################################
#Extracting Cluster information. 

cluster<- colnames(data)
cluster<- as.data.frame(cluster)
cluster$gene<- substr(cluster$cluster, 1, 4)
colnames(cluster)<- c("gene.name","gene")
cluster<- as.data.frame(cluster)
cluster<- cluster[-1,]
#########
#From the paper, the 'cy..' marke is the patient ID. In figure 1m panel C, only 6 patients' cluster.
#cy80= Mel80
#cy78= Mel78
#cy79= Mel79
#cy81= Mel81
#cy84= Mel84
#cy88= Mel88
#All other cy.. = Unknown.
#summary(cluster$gene)
#cy53 Cy59 cy60 CY65 Cy71 CY75 cy78 cy79 cy80 Cy80 CY80 cy81 Cy81 cy82 cy84 CY84 cy88 CY88 CY89 
#16   54    9    4   54    3  120  468   50   40   35   70   63   32   11    3    1  116   98 
#cy94 CY94 
#6    4 

##Adding Cluster Numbers

cluster$group<- ifelse(cluster$gene == "cy80", 1.1,
                       ifelse(cluster$gene == "Cy80", 1.1,
                              ifelse(cluster$gene == "CY80", 1.1,
                       ifelse(cluster$gene == "cy78", 1.2,
                              ifelse(cluster$gene == "cy79", 1.3,
                                     ifelse(cluster$gene == "cy81", 1.4,
                                            ifelse(cluster$gene == "Cy81", 1.4,
                                                   ifelse(cluster$gene == "cy84", 1.5,
                                                          ifelse(cluster$gene == "CY84", 1.5,
                                                                 ifelse(cluster$gene == "cy88", 1.6,
                                                                        ifelse(cluster$gene == "CY88", 1.6, 1.7 )))))))))))

###############################################
##STEP5: Print ClusterID for 6/13 patients######
###############################################
##This is the ClusterID file for 13 patients, where 1.7 mean other 7 unknown patients' IDs. Need to upload in google drive.
write.table(cluster[,3], file = "clusterID.txt", sep = "\t", row.names = F,col.names = F)
cluster_6p<- cluster[cluster$group != '1.7', ] #Filtering for only 6 patients.
#This is clusterID for only 6 patients reported in panel c. 
write.table(cluster_6p[,3], file = "clusterID_6p.txt", sep = "\t", row.names = F,col.names = F)

###############################################
##STEP6: Print datamatrix for 6 patients######
###############################################
#Now we want to filter the datamatrix for only 6 patients. 
data.malig.6p<- data[ ,c(colnames(data) %in% cluster_6p$gene.name)]
gene.name<- data$gene.name
data.malig.6p<- cbind(gene.name, data.malig.6p)
#This is the final dataset for the paper's panel c representing 6 patients.
write.csv(data.malig.6p, file = "data.malig.6p.csv", sep = "\t", row.names = F,col.names = T)

###############################################
##STEP7: Print ClusterName for 6/13 patients######
###############################################
####
#Cluster name file.
Patient<- c(1:7)
ClusterID<- c('1.1', '1.2','1.3','1.4','1.5','1.6','1.7')
ClusterName<- c('Mel80',
                'Mel78',
                'Mel79',
                'Mel81',
                'Mel84',
                'Mel88',
                'Unknown')
ClusterNames_New<- cbind(Patient, ClusterID, ClusterName) # This ClusterNames is for 13 patients.
ClusterNames_New_6p<- ClusterNames_New[1:6,] # This ClusterNames is for 6 patients. 
##This are the ClusterName files. Need to upload in google drive.
write.csv(ClusterNames_New, file = "ClusterNames_New.csv", sep = "\t", row.names = F,col.names = T)# for 13 patieints
write.csv(ClusterNames_New_6p, file = "ClusterNames_New_6p.csv", sep = "\t", row.names = F,col.names = T) # for 6 patients


##All downstream analysis is done in Scanpy at Google Colab.
