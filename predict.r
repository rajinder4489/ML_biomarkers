#calculates different parameters for	whole list, only protien coding and only non-coding
#Removes highly correlated festures (cutoff=0.75)

#setwd("C:/Users/r.gupta/Desktop/livercancer/")

#loading libraries
library(caret)
library(MLeval)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
library(openxlsx)
library(cowplot)
library(openxlsx)
#library(pROC)

#output files 
stats_file <- "stats.txt"
write.table("Start", file=stats_file, append=FALSE)


######### Functions Start #############

#To preprocess the data; returns train and test dfs
preprocessing <- function(df_main, label)
{
	ret <- list()		#list to return, would contain 2 df: train and test

	df_t <- t(df_main)
	samples <- rownames(df_t)
	df_t <- cbind.data.frame(df_t, samples)
	df_t$samples <- as.character(df_t$samples)


	#Adding the labels using the label file
	df_l <- merge(df_t, label[,c(1,2)], by="samples")
	df_l$lbl <- factor(df_l$lbl)

	#Now remove the sample names; not required anymore
	rownames(df_l) <- df_l$samples
	df_l <- df_l[ , !names(df_l) %in% c("samples")]

	#splitting in train and test
	intraining <- createDataPartition(df_l$lbl, p=0.7, list=FALSE)
	ret$train <- df_l[intraining,]
	ret$test <- df_l[-intraining,]
	
	return(ret)
}

#Find highly correlated features (transcripts)
findhighlyCor <- function(dataframe)
{
	#calculate correlation matrix
	correlationMatrix <- cor(dataframe)	#transpose is done so that the correlaiton is calculated for the transcripts and not samples
	
	# find attributes that are highly correlated (ideally >0.75)
	highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75, names = TRUE)
	
	# print indexes of highly correlated attributes
#	print(paste("Highly correlated features (transcripts):", length(highlyCorrelated)))
	
	return(highlyCorrelated)
}

#Finding the important features (transcripts) for the model
findImp <- function(m)
{
	#estimate variable importance
	importance <- varImp(m, scale=FALSE)
	#summarize importance
#	print(importance)
	#plot importance
	plot(importance)
	return(importance)
}

### Prediction function
pred.function <- function(a, b, d)	# model, test, type: 1 for KNN, RF,NB and NNET; 2 for SVM
{
	if(d==1)
	{
		pred <- predict(a, b, "prob")
		pred$pr <- ifelse(pred[,1]>=0.5, 1, 2)
	}
	else
	{
		pred <- as.data.frame(predict(a, b)) #for using with SVM
		pred <- as.data.frame(cbind(pred, b[,"lbl"]))
		colnames(pred) <- c("pr", "lbl")
	}
	
	confusion_matrix = table(b[,"lbl"], pred$pr)
	return(confusion_matrix)
}
####

### All steps

allsteps <- function(dfs, f, excel_file)
{
#To remove the transcripts which are low expressed in train dataset
#	l <- dfs$train[["lbl"]]#extracting the lables
	train_fltr <- dfs$train[,c(c(colSums(dfs$train[, -which(names(dfs$train) %in% c("lbl"))])>f), TRUE)]		#also takes the factor column "lbl" somehow
	print(paste("After passing expression filter; Samples:", nrow(train_fltr), "Transcripts", ncol(train_fltr), sep=" "))
	
	write.table(paste("After passing expression filter; Samples:", nrow(train_fltr), "Transcripts", ncol(train_fltr), sep=" "), file = stats_file, append =TRUE)
	
#To find the highly correlated transcripts
	highlyCor <- findhighlyCor(train_fltr[, -which(names(train_fltr) %in% "lbl")])	#returns the names of the highly correlated features (transcripts)
#	highlyCor <- findhighlyCor(train_fltr)	#returns the names of the highly correlated features (transcripts)
	
	print(paste("Highly correlated features", length(highlyCor), sep=" "))
	write.table(paste("Highly correlated features (transcripts):", length(highlyCor)), file= stats_file, append = TRUE)
	
	removehc <- highlyCor[2:length(highlyCor)]	#1st is kept
	
#	if(itr <= 2)	#enter this for biomarkers 
#	{
#		train_subset <- train_fltr		#done this for biomarkers, did not remove the highly correlated features for them
#	}
#	else	#enter this for whole dataset
#	{
#		train_subset <- train_fltr[, -which(colnames(train_fltr) %in% removehc)]
#	}
	#remove all correlated features
	train_subset <- train_fltr[, -which(colnames(train_fltr) %in% removehc)]
	
	write.table(paste("Number of transcripts after removing highly correlated:", ncol(train_subset)), file= stats_file, append = TRUE)

#Running the model
	#Creating control
	fitControl <- trainControl(method = "repeatedcv", number = 5, repeats = 5, classProbs=T, savePredictions = T)

	#Running the model on train
	print("Running the ML models now")
	
	models <- list()		#empty list created
	runtime <- list()
	
	print("RF")
	runtime$rf_t <- system.time(models$r_f <- train(lbl ~., data=train_subset, method="rf", trControl=fitControl))
#	write.table(paste("RF", runtime$rf_t[3]), file= stats_file, append = TRUE)
	
	print("NB")
	runtime$nb_t <- system.time(models$nb <- train(lbl ~., data=train_subset, method="naive_bayes", trControl=fitControl))
#	write.table(paste("NB", runtime$nb_t[3]), file= stats_file, append = TRUE)
	
	print("KNN")
	runtime$knn_t <- system.time(models$k_nn <- train(lbl ~., data=train_subset, method="knn", trControl=fitControl))
#	write.table(paste("KNN", runtime$knn_t[3]), file= stats_file, append = TRUE)
	
	print("SVM")
	runtime$svm_t <- system.time(models$svm <- train(lbl ~ ., data=train_subset, method = 'svmLinear', preProcess = NULL, trainControl = fitControl, scale = FALSE, probability=T))
#	write.table(paste("SVM", runtime$svm_t[3]), file= stats_file, append = TRUE)
	
	print("NNET")
	runtime$nnet_t <- system.time(models$nnet <- train(lbl ~., data=train_subset, method="nnet", trControl=fitControl, MaxNWts=100000, verbose=FALSE))
#	write.table(paste("NNET", runtime$nnet_t[3]), file= stats_file, append = TRUE)
	
	conf_mat <- list()		#confusion matrix
	impfeatures <- list()	#important features
	
	rnames <- rownames(dfs$test)
	
	#creating a workbook
	wb <- createWorkbook()
	
	#Testing the model on test data
	for(j in 1:length(models))
	{
		nam <- names(models)[j]
		
		#predict on test
		if(nam != "svm")
		{
			pred <- predict(models[[j]], dfs$test, "prob")
			pred$pr <- ifelse(pred[,1]>=0.5, 1, 2)
		}
		else
		{
			pred <- as.data.frame(predict(models[[j]], dfs$test)) #for using with SVM
			pred <- as.data.frame(cbind(pred, dfs$test[,"lbl"]))
			colnames(pred) <- c("pr", "lbl")
		}
		
		pred$samples <- rnames
		towrite <- cbind(dfs$test[,"lbl"], pred[, c("samples", "pr")])
#		print(towrite)
		
		addWorksheet(wb, nam)
		writeData(wb, nam, towrite, rowNames = FALSE)
		
		#Creating confusion matrix
#		print(cbind(dfs$test[,"lbl"], pred$pr))
		conf_mat[[nam]] = table(dfs$test[,"lbl"], pred$pr)
#		print(paste("+", table(dfs$test[,"lbl"], pred$pr), sep = " "))
		#Finding the important features (transcripts)
		impfeatures[[nam]] <- findImp(models[[j]])
	}
	saveWorkbook(wb, file = excel_file, overwrite = TRUE)

	#Calculating AUC-ROC
	auc_roc <- evalm(list(models$r_f, models$nb, models$k_nn, models$nnet), gnames = c("RF", "NB", "KNN", "NNET"))		#for SVM it does not work
	
#	evalm(list(models$nnet), gnames = c("NNET"))
	allres <- list()		#all results, empty list
	
	allres$models <- models
	allres$conf_mat <- conf_mat
	allres$impfeatures <- impfeatures
	allres$auc_roc <- auc_roc
	allres$runtime <- runtime
	
	return(allres)
}

####

### Plot importance
importancegraphs <- function(lists, d, filename)
{
	#Importance graphs
	nams <- names(lists)
	titles <- c("Random Forest", "Naive Bayes", "KNN", "SVM", "NNET")		#known names for the algos used; used for printing on the graph
	tags <- LETTERS[1:length(nams)]
	plots <- list()		#empty list
	
	wb <- createWorkbook()

	for(z in 1:length(nams))
	{
	  imp <- lists[[nams[z]]][["importance"]]
		imp$id <- rownames(imp)
		imp_subset <- head(imp[order(imp[,1], decreasing = T),], 20)

		imp_name <- merge(imp_subset, d[, c("Transcript.stable.ID", "Transcript.name")], by.x="id", 
		                  by.y="Transcript.stable.ID", all.x=TRUE)
		
		#To replace the NA from the name, for the one which do not have a name value in Description file
		imp_name$Transcript.name <- ifelse(is.na(imp_name$Transcript.name), imp_name$id, imp_name$Transcript.name)
		
		melted <- melt(imp_name[,c("Transcript.name", colnames(imp_name)[2])])
		print(dim(melted))
		
		plots[[z]] <- ggplot(data=melted, aes(x=reorder(Transcript.name, value), y=value, fill=variable)) +
		geom_bar(stat="identity", position=position_dodge()) +
		coord_flip() +
		theme(axis.text.x = element_text(angle = 90), legend.position = "none") +
		labs(title = titles[z], x = "Transcripts", y= "Importance", tag=tags[z])
		
		
		addWorksheet(wb, nams[z])
		writeData(wb, nams[z], imp_name, rowNames = FALSE)
	}
	
	saveWorkbook(wb, file = filename, overwrite = TRUE)

	return(plots)
}
####

metric_graph <- function(df_g, ttl, ta)
{
  g <- list()
   
  level_order <- c("Senstivity", "Specificity", "MCC", "Informedness")
  
  g$graph <- ggplot(data=df_g, aes(x = factor(Metric, level = level_order), y = Value, fill = ML)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(legend.position = "right") +
  labs(title = ttl, x = "", y= "score", tag=ta) +
  scale_fill_manual("ML algo", values = c("RF" = "red", "NB" = "purple", "KNN" = "grey", "NNET" = "gold"))
  
  g_for_legend <- ggplot(data=df_g, aes(x = factor(Metric, level = level_order), y = Value, fill = ML)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme(legend.position = "right") +
  labs(title = ttl, x = "", y= "score", tag=ta) +
  scale_fill_manual("ML algo", values = c("RF" = "red", "NB" = "purple", "KNN" = "grey", "NNET" = "gold"))
  
  
  legends <- cowplot::get_legend(g_for_legend)
  
  p_legend <- ggplot(data=NULL, aes()) +
  geom_blank(mapping = NULL, data = NULL, stat = "identity",
  position = "identity", show.legend = TRUE, inherit.aes = TRUE) +
  labs(x = NULL, y = NULL, fill = NULL) +
  scale_fill_discrete(name = "", labels = legends)

  g$legends <- p_legend
  return(g)
}

performance_metrics <- function(dff, type, ta)	#df for data, type of transcript and tag
{
	df_all_metric <- c()
	
	for(i in 1:length(names(dff)))
	{
		df_all_metric <- rbind(df_all_metric, c("Senstivity", names(dff)[i], 
						dff[[names(dff)[i]]][1,1]))
		
		df_all_metric <- rbind(df_all_metric, c("Specificity", names(dff)[i], 
						dff[[names(dff)[i]]][2,1]))
		
		df_all_metric <- rbind(df_all_metric, c("MCC", names(dff)[i], 
						dff[[names(dff)[i]]][3,1]))
		
		df_all_metric <- rbind(df_all_metric, c("Informedness", names(dff)[i], 
						dff[[names(dff)[i]]][4,1]))
	}
	
	colnames(df_all_metric) <- c("Metric", "ML", "Value")
	df_all_metric <- as.data.frame(df_all_metric)
	
	write.table(df_all_metric, file=stats_file, sep="\t", append=TRUE)
	
	graph <- metric_graph(df_all_metric, type, ta)
	
	return(graph)
}

######### Functions End #############

### Main Starts ###
#rm(list = setdiff(ls(), lsf.str()))

#fltr <- c(0, 0, 0) #Filtering the expression; used in order for all, protein coding and non-coding
fltr <- c(10000, 10000, 10000) #Filtering the expression; used in order for all, protein coding and non-coding

set.seed(100)

files <- c("biomarkers_iso_norm_reads.txt", "biomarkers_iso_norm_fpkm-FuSe.txt",
			 "isoform_norm_reads.txt", "isoform_norm_fpkm-FuSe.txt")

name <- c("bm_iso_normreads", "bm_iso_normfpkmFuSe",
			 "all_iso_normreads", "all_iso_normfpkmFuSe")
name <-  paste(name, fltr[1], sep="_")

label <- read.delim(file = "labels.txt", sep="\t", stringsAsFactors = T)
desc <- read.delim(file = "biomart_desc.txt", sep="\t", stringsAsFactors = F)

#for(i in 1:length(files))
#for(i in 1:2)
for(i in 3:4)
{
#reading files
	file1 <- read.delim(file = files[i], sep="\t", stringsAsFactors = F)
	colnames(file1) <- gsub("Cypr0tex", "Cyprotex", colnames(file1))
	
	write.table(files[i], file=stats_file, sep="\t", append = TRUE)
	
	if(i <= 2)
	{
		#adding the transcript ids as the rownames
		rownames(file1) <- file1$id
		
		#removing the column: ids
		file1 <- file1[, -which(names(file1) %in% c("id"))]
	}
	
#Splitting in train and test and other preprocessing	
	dataframes <- preprocessing(file1, label)		#returns 2 dfs: train and test
	
#Running analyses for all transcripts
	print("====================================================================================")
	print("All transcripts")
	write.table("===================================================================================+", file=stats_file, sep="\t", append =TRUE)
	write.table("All Transcripts", file=stats_file, sep="\t", append =TRUE)
	
	all_excel <- paste(name[i], "_pred.xlsx", sep="")
	all_trans_res <- allsteps(dataframes, fltr[1], all_excel)
	
	saveRDS(all_trans_res, file = paste(name[i], "all_res.rds", sep="_"))
	
	write.table(paste("RF:", all_trans_res$runtime$rf_t[3], "NB:", all_trans_res$runtime$nb_t[3], "KNN:", all_trans_res$runtime$knn_t[3], "SVM:", all_trans_res$runtime$svm_t[3], "NNET:", all_trans_res$runtime$nnet_t[3], sep="\t"), file=stats_file, sep="\t", append = TRUE)
	
	#importance graph
	filename1 <- paste(name[i], "imp_all_biomarkers.xlsx", sep="_")
	filename2 <- paste(name[i], "imp_all_biomarkers.jpeg", sep="_")
	
	all_imp_graphs <- importancegraphs(all_trans_res$impfeatures, desc, filename1)
	
	all_g <- grid.arrange(all_imp_graphs[[1]], all_imp_graphs[[2]], all_imp_graphs[[3]], all_imp_graphs[[4]], all_imp_graphs[[5]], 
	                      top=textGrob("Importance - All transcripts", gp=gpar(fontsize=18,font=8)), 
	                      layout_matrix = rbind(c(1,2,3), c(4,5,NA)))
	
	ggsave(file=filename2, plot = all_g, device = NULL, path = NULL, scale = 1, width = 15, height = 10, units = "in", dpi = 200, limitsize = TRUE)
	
	gc()
	
#getting the protein coding and non coding transcripts
	pc <- desc[desc$Transcript.type == "protein_coding", c("Transcript.stable.ID", "Transcript.type")]	#protein coding
	nc <- desc[desc$Transcript.type != "protein_coding", c("Transcript.stable.ID", "Transcript.type")]	#non-coding
	
#Running analyses for protein coding
	dataframes_p <- list()
	dataframes_p$train <- dataframes$train[, -which(names(dataframes$train) %in% nc[,"Transcript.stable.ID"])]		#removing the non-coding ones
	dataframes_p$test <- dataframes$test[, -which(names(dataframes$train) %in% nc[,"Transcript.stable.ID"])]		#removing the non-coding ones
	print("Protein coding transcripts")

	p_excel <- paste(name[i], "_pred.xlsx", sep="")
	p_trans_res <- allsteps(dataframes_p, fltr[2], p_excel)		#protien coding transcripts
	
	saveRDS(p_trans_res, file = paste(name[i], "p_res.rds", sep="_"))
	
	write.table(paste("RF:", p_trans_res$runtime$rf_t[3], "NB:", p_trans_res$runtime$nb_t[3], "KNN:", p_trans_res$runtime$knn_t[3], "SVM:", p_trans_res$runtime$svm_t[3], "NNET:", p_trans_res$runtime$nnet_t[3], sep="\t"), file=stats_file, sep="\t", append = TRUE)
	
	#importance graph
	filename3 <- paste(name[i], "imp_prot_biomarkers.xlsx", sep="_")
	filename4 <- paste(name[i], "imp_prot_biomarkers.jpeg", sep="_")
	
	p_imp_graphs <- importancegraphs(p_trans_res$impfeatures, desc, filename3)
	
	p_g <- grid.arrange(p_imp_graphs[[1]], p_imp_graphs[[2]], p_imp_graphs[[3]], p_imp_graphs[[4]], p_imp_graphs[[5]], 
	                    top=textGrob("Importance - Protien coding transcripts", gp=gpar(fontsize=18,font=8)), 
	                    layout_matrix = rbind(c(1,2,3), c(4,5,NA)))
	
	ggsave(file=filename4, plot = p_g, device = NULL, path = NULL, scale = 1, width = 15, height = 10, units = "in", dpi = 200, limitsize = TRUE)
	gc()
	
#Running analyses for non-coding
	dataframes_n <- list()
	dataframes_n$train <- dataframes$train[, -which(names(dataframes$train) %in% pc[,"Transcript.stable.ID"])]		#removing the protein coding ones
	dataframes_n$test <- dataframes$test[, -which(names(dataframes$train) %in% pc[,"Transcript.stable.ID"])]		#removing the protein coding ones
	print("Non-coding transcripts")
	
	n_excel <- paste(name[i], "_pred.xlsx", sep="")
	n_trans_res <- allsteps(dataframes_n, fltr[3], n_excel)		#non-coding transcripts
	
	saveRDS(n_trans_res, file = paste(name[i], "n_res.rds", sep="_"))
	
	write.table(paste("RF:", n_trans_res$runtime$rf_t[3], "NB:", n_trans_res$runtime$nb_t[3], "KNN:", n_trans_res$runtime$knn_t[3], "SVM:", n_trans_res$runtime$svm_t[3], "NNET:", n_trans_res$runtime$nnet_t[3], sep="\t"), file=stats_file, sep="\t", append = TRUE)
	
	#importance graph
	filename5 <- paste(name[i], "imp_nc_biomarkers.xlsx", sep="_")
	filename6 <- paste(name[i], "imp_nc_biomarkers.jpeg", sep="_")
	
	n_imp_graphs <- importancegraphs(n_trans_res$impfeatures, desc, filename5)
	
	n_g <- grid.arrange(n_imp_graphs[[1]], n_imp_graphs[[2]], n_imp_graphs[[3]], n_imp_graphs[[4]], n_imp_graphs[[5]], 
	                    top=textGrob("Importance - Non-coding transcripts", gp=gpar(fontsize=18,font=8)), 
	                    layout_matrix = rbind(c(1,2,3), c(4,5,NA)))
	
	ggsave(file=filename6, plot = n_g, device = NULL, path = NULL, scale = 1, width = 15, height = 10, units = "in", dpi = 200, limitsize = TRUE)
	
	gc()

#Plotting graphs
	#ROC-AUC graphs
	print("Plots")
	#ROC
	a1 <- all_trans_res$auc_roc$roc + labs(title = "All transcripts", tag="A")
	a2 <- p_trans_res$auc_roc$roc + labs(title = "Protien coding transcripts", tag="B")
	a3 <- n_trans_res$auc_roc$roc + labs(title = "Non-coding transcripts", tag="C")
	
	g1 <- grid.arrange(a1, a2, a3, top=textGrob("ROC (Receiver operating characteristic)", gp=gpar(fontsize=18,font=8)), 
	                   layout_matrix = rbind(c(1,2), c(3,NA)))
	
	filename7 <- paste(name[i], "roc.jpeg", sep="_")
	ggsave(file=filename7, plot = g1, device = NULL, path = NULL, scale = 1, width = 12, height = 10, units = "in", dpi = 200, limitsize = TRUE)
	
	#
	a1 <- all_trans_res$auc_roc$proc + labs(title = "All transcripts", tag="A")
	a2 <- p_trans_res$auc_roc$proc + labs(title = "Protien coding transcripts", tag="B")
	a3 <- n_trans_res$auc_roc$proc + labs(title = "Non-coding transcripts", tag="C")
	
	g2 <- grid.arrange(a1, a2, a3, top=textGrob("Precision-Recall", gp=gpar(fontsize=18,font=8)), 
	                   layout_matrix = rbind(c(1,2), c(3,NA)))
	
	filename8 <- paste(name[i], "proc.jpeg", sep="_")
	ggsave(file=filename8, plot = g2, device = NULL, path = NULL, scale = 1, width = 12, height = 10, units = "in", dpi = 200, limitsize = TRUE)
	
	#
	a1 <- all_trans_res$auc_roc$prg + labs(title = "All transcripts", tag="A")
	a2 <- p_trans_res$auc_roc$prg + labs(title = "Protien coding transcripts", tag="B")
	a3 <- n_trans_res$auc_roc$prg + labs(title = "Non-coding transcripts", tag="C")
	
	g3 <- grid.arrange(a1, a2, a3, top=textGrob("Precision-Recall gain", gp=gpar(fontsize=18,font=8)), 
	                   layout_matrix = rbind(c(1,2), c(3,NA)))
	
	filename9 <- paste(name[i], "prg.jpeg", sep="_")
	ggsave(file=filename9, plot = g3, device = NULL, path = NULL, scale = 1, width = 12, height = 10, units = "in", dpi = 200, limitsize = TRUE)

	#Calibration curve
	a1 <- all_trans_res$auc_roc$cc + labs(title = "All transcripts", tag="A")
	a2 <- p_trans_res$auc_roc$cc + labs(title = "Protien coding transcripts", tag="B")
	a3 <- n_trans_res$auc_roc$cc + labs(title = "Non-coding transcripts", tag="C")
	
	g4 <- grid.arrange(a1, a2, a3, top=textGrob("Calibration Curve", gp=gpar(fontsize=18,font=8)), layout_matrix = rbind(c(1,2), c(3,NA)))
	
	filename10 <- paste(name[i], "cc.jpeg", sep="_")
	ggsave(file=filename10, plot = g4, device = NULL, path = NULL, scale = 1, width = 12, height = 10, units = "in", dpi = 200, limitsize = TRUE)
	
	
	b1 <- performance_metrics(all_trans_res$auc_roc$stdres, "All transcripts", "A")
	b2 <- performance_metrics(p_trans_res$auc_roc$stdres, "Protein coding transcripts", "B")
	b3 <- performance_metrics(n_trans_res$auc_roc$stdres, "Non-coding transcripts", "C")
	
	filename11 <- paste(name[i], "stats.jpeg", sep="_")
	
	g5 <- grid.arrange(b1$graph, b2$graph, b3$graph, layout_matrix = rbind(c(1,2), c(3,NA)))
	ggsave(file=filename11, plot = g5, device = NULL, path = NULL, scale = 1, width = 12, height = 10, units = "in", dpi = 200, limitsize = TRUE)
	
}


########## Main Ends #########
