### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ###  Libraries  ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
library(rgl)
library(ROCR)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# setting work directory 
dir_input_surf <- "/Users/edo/Documents/Progetti/COV-19/new_surf_screened/"
setwd(dir_input_surf)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
args <- commandArgs(trailingOnly = TRUE)
R_sphere_sm <- 6
plot_aus <- 0

  file1 <- as.character(args[1]) # args[1] contains "100"
  file2 <- as.character(args[2]) # args[2] contains "10"
  
  name_complex <- unlist(strsplit(file1, "_"))[3]
  
  print("********************")   
  print("The process is running...")
  
  surfProt1_aus <- read.table(file1)
  surfProt2_aus <- read.table(file2)
  surfProt1Coord <- surfProt1_aus[,1:3]
  surfProt2Coord <- surfProt2_aus[,1:3]
  
  surfProt1_aus <- surfProt1_aus[,c(1,2,3,4,5)]
  surfProt2_aus <- surfProt2_aus[,c(1,2,3,4,5)]
      
  if(ncol(surfProt2_aus) == 5 & ncol(surfProt1_aus) == 5){
        
    colnames(surfProt1_aus) <- c("x","y","z", "col","BS_Dist")
    colnames(surfProt2_aus) <- c("x","y","z", "col","BS_Dist")
       
    if(plot_aus == 1){
      rbPal <- colorRampPalette(c('red','gray','gray','gray','gray'))
      colaus <- rbPal(30)[as.numeric(cut(surfProt1_aus$col,breaks = 25))]
      plot3d(surfProt1_aus[,c(1,2,3)],col=colaus,box=F,axes=F, size=3, xlab = "", ylab = "", zlab = "", radius = 5)
    
      rbPal <- colorRampPalette(c('gray','gray','gray','gray','green'))
      colaus <- rbPal(30)[as.numeric(cut(surfProt2_aus$col,breaks = 25))]
      plot3d(surfProt2_aus[,c(1,2,3)],col=colaus,box=F,axes=F, size=3, xlab = "", ylab = "", zlab = "", radius = 5, add = T)
    } 
    
    # smooth protein 1
    surfProt1_aus$Col_sm <- 0
    distSurf <- as.matrix(dist(surfProt1Coord))
    for (i in 1:nrow(surfProt1Coord)){
      #print(i)
      mean_aus <- mean(surfProt1_aus$col[distSurf[i,] <= R_sphere_sm][!is.na(surfProt1_aus$col[distSurf[i,] <= R_sphere_sm])])
      surfProt1_aus[i,"Col_sm"] <- mean_aus
    }
    
    if(plot_aus == 1){
      rbPal <- colorRampPalette(c('red','gray','gray','gray','gray'))
      colaus <- rbPal(30)[as.numeric(cut(surfProt1_aus$Col_sm,breaks = 25))]
      plot3d(surfProt1_aus[,c(1,2,3)],col=colaus,box=F,axes=F, size=3, xlab = "", ylab = "", zlab = "", radius = 5)
    } 
    
    # smooth protein 2
    surfProt2_aus$Col_sm <- 0
    distSurf <- as.matrix(dist(surfProt2Coord))
    for (i in 1:nrow(surfProt2Coord)){
      #print(i)
      mean_aus <- mean(surfProt2_aus$col[distSurf[i,] <= R_sphere_sm][!is.na(surfProt2_aus$col[distSurf[i,] <= R_sphere_sm])])
      surfProt2_aus[i,"Col_sm"] <- mean_aus
    }
    
    if(plot_aus == 1){
      rbPal <- colorRampPalette(c('green','gray', 'gray', 'gray', 'gray'))
      colaus <- rbPal(30)[as.numeric(cut(surfProt2_aus$Col_sm,breaks = 25))]
      plot3d(surfProt2_aus[,c(1,2,3)],col=colaus,box=F,axes=F, size=3, xlab = "", ylab = "", zlab = "", radius = 5)
    }
    
    pred <- prediction(surfProt1_aus$Col_sm, surfProt1_aus$BS_Dist)
    perf <- performance(pred,"tpr","fpr")
    auc_ROCR <- performance(pred, measure = "auc")
    Ab_auc_ROCR <- 1-auc_ROCR@y.values[[1]]
    
    pred <- prediction(surfProt2_aus$Col_sm, surfProt2_aus$BS_Dist)
    perf <- performance(pred,"tpr","fpr")
    auc_ROCR <- performance(pred, measure = "auc")
    An_auc_ROCR <- 1-auc_ROCR@y.values[[1]]

    if(plot_aus == 1){
      plot(perf,colorize=TRUE)
      plot(perf,colorize=TRUE)
    }
    
    vet_info <- c(round(Ab_auc_ROCR,2), round(An_auc_ROCR,2), name_complex)
    
    print("********************")   
    print(paste("The complex analyzed is: ", vet_info[3], sep=""))   
    print(paste("The AUC of the first protein is: ", vet_info[1], sep=""))
    print(paste("The AUC of the second protein is: ", vet_info[2], sep=""))
    print("********************")   
    
    
  }   
      
    

