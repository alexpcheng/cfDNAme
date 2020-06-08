#!/usr/bin/env RScript
# Title: [ENTER TITLE]
# Authors: [ENTER AUTHORS]
# Brief description: [ENTER DESCRIPTION]

# Functions for tissue of origin measurement
group_by_celltype <- function(df){
  counter <- 1
  for (i in 1:ncol(df)){
    tissue_of_interest <- gsub("\\d+$", "", colnames(df)[[i]]) #tissue of the given column
    if (counter == 1){
      new_df <- data.frame(rowMeans(df[, grepl(tissue_of_interest, colnames(df)), drop=FALSE]))
      colnames(new_df)[counter] <- tissue_of_interest
      counter <- counter + 1
    }
    else{
      if (!(tissue_of_interest %in% colnames(new_df))){
        new_df <- cbind(new_df, rowMeans(df[, grepl(tissue_of_interest, colnames(df)), drop=FALSE]))
        colnames(new_df)[counter] <- tissue_of_interest
        counter <- counter + 1
      }
    }
  }
  return(new_df)
}
  
mixing_parameters_function<-function(b, A, other, sample_name, sum_to_one){
  #Accounting for missing tissues
  if(other){
    A$other<-1
    G<-diag(ncol(A))
    g<-rep(-1,ncol(A))
    g[ncol(A)]<-0
    G<-rbind(G, g)
    h<-rep(0,ncol(A))
    h[ncol(A)+1]<-(-1)
    E<-NULL
    f<-NULL
    sol<-lsei(A=A, B=b, G=G, H=h, E=E, F=f, type=2, fulloutput = T) #calls solve.QP from package quadprog
    sol.df<-data.frame(sol$X)
    sol.df$tissue<-rownames(sol.df)
    #sol.df<-sol.df[!rownames(sol.df) %in% "other",]
    if (sum_to_one){
      normalizing_factor<-sum(sol.df$sol.X)
      sol.df$sol.X<-sol.df$sol.X/normalizing_factor
    }
    rel_error<-(sol$solutionNorm/nrow(A))
    rel_error<-data.frame(rel_error, "rel_error")
    colnames(rel_error)<-c(sample_name, "tissue")
  }
  #Not Accounting for missing tissues
  if(!other){
    G<-diag(ncol(A))
    h<-rep(0,ncol(A))
    E<-(rep(1,ncol(A)))
    f<-1
    sol<-lsei(A=A, B=b, G=G, H=h, E=E, F=f, type=1, fulloutput = T)
    sol.df<-as.data.frame(sol$X)
    sol.df$tissue<-rownames(sol.df)
    if (sum_to_one){
      normalizing_factor<-sum(sol.df$sol.X)
      sol.df$sol.X<-sol.df$sol.X/normalizing_factor
    }
    rel_error<-(sol$solutionNorm/nrow(A))
    rel_error<-data.frame(rel_error, "rel_error")
    colnames(rel_error)<-c(sample_name, "tissue")
  }
  colnames(sol.df)[1]<-sample_name
  sol.df<-rbind(sol.df, rel_error)
  return(sol.df)
}


# Main function to call ----------------------------------------------------------------------------
mp_function<-function(reference_and_samples, other, nt, sum_to_one, group_tissues, sample_name){
  A<-as.data.frame(reference_and_samples[,4:(nt+3)])
  if (group_tissues){
    A<-group_by_celltype(A)
  }
  
  b<-reference_and_samples[,(nt+4):ncol(reference_and_samples)]
  percentages <- mixing_parameters_function(b, A, other, sample_name, sum_to_one)
  p <- data.frame(t(percentages))[1,]
  colnames(p)[ncol(p)]<-'rel_error'
  p$sample <- rownames(p)
  return(p)
}

