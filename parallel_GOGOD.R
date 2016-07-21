## GoGod.R
## Function for gene prioritization (finding causative genes for a phenotype)
## Serial version to be parallelized
##

library(GOSemSim)


# Input:
# filein is list of probes to prioritize
# list is the GO annotations List to compare with selected genes.
# ontol is the ontology that is used, either "CC", "MF", or "BP"
# measure is the semantic similarity measure used, e.g. "Lin"
goGoD <-
  function(filein,list,ontol,measure, file_choice) {
    
   # ptm <- proc.time()
    
    
    
    print("SELECT INPUT FILE ")
    
    
    
    
    
    # Read a data file, chosen by the user, as a table.
    #dati=read.table(file.choose(), header = TRUE, sep = "\t")
    dati=read.table("data.txt", header = TRUE, sep = "\t")
    # Converts the data to be in matrix form.
    dati1=as.matrix(dati)
    A=filein
    
    # b is the number of 
    b=length(A)
    
    # n1 is the index of the probe.
    n1=numeric(length=b)
    # Find indices of probes
    
    for(i in 1:length(A)){
      # assign to n1[i] the index of the probe A[i] in dati...
      # e.g. if A[i] = AM_10192, then n[i] = 1 since its the first probe in the data fi
      n1[i]=which(dati1==A[i])
    }
    # n2 is vector containing b 0's.
    n2=numeric(length=b)
    for(i in 1:length(n1)){
      # assign to n2[i] the GO annotations of probe n1[i], which are in 4th column of matrix.
      n2[i]=dati1[n1[i],4]
    }
    
    # Create a matrix out of the GO annotations
    # Each row contains GO terms.
    n3=as.matrix(n2)
    
    for (i in 1:dim(n3)[1]){
      # paste the GO annotations of each probe
      write(n3[i],paste0(i,"annotation_GO.txt", sep = ""))
    }
    # x is the vector of (1, 2, ..., length(n3))
    x=seq(1:length(n3))
    
    f=0
    for(i in 1:length(x)){
      # each element of f is a list of the annotations.
      f[i]=list(read.table(paste0(i,"annotation_GO.txt")))
    }
    
    # label the lists with their probe...
    # e.g. f[i] labeled with A[i], which is AM_10430
    names(f) = A;
    
    h=as.matrix(f[[1]])
    if (is.na(h[i])==TRUE){
      # is.na --> not available --> there aren't anny annotations for the probe.
      print(paste0("No annotations ", names(f)[i]))
    }
    score=numeric(length(n3))
    
    for(i in 1:length(x)){
      # h is f read in as a matrix.
      h = as.matrix(f[[i]])
      if (is.na(h[i])==TRUE){
        # is.na --> not available --> there aren't anny annotations for the probe.
        print(paste0("No annotations ", names(f)[i]))
      }
      
      # calculate semantic similarity between h, which a probe's annotations, and the list of terms
      # "list" (input into function)
      score[i] = mgoSim(h, list, ont=ontol, organism="human", measure=measure)
    }
    # label scores with probe names
    names(score)<-A
    
    # round scores to 3 digits
    score = round(score, digits = 3)
    write.table(score,paste(measure,"score.txt"),row.names=TRUE, col.names = FALSE)
    
    #    assign probe names to this variable, used in the x-axis of the histogram.
    names(score)->ll
    return (score)     
  }

parallel_GOGOD <- function(filein,list,ontol,measure, procs, file_choice){
  gene_probes <- filein
  GO_terms <- list
  # Set up parallelization by registering cluster/preparing processes.
  library(parallel)
  num_cores <- procs
  num_probes <- length(gene_probes)
  sublist_length <- ceiling(length(gene_probes)/procs)
  sublists <- split(gene_probes, rep(1:ceiling(num_probes/sublist_length), each=sublist_length)[1:num_probes])
  sublists <- unlist(sublists)
  library(foreach)
  library(doParallel)
  library(doMC)
  ptm<-proc.time()
  cl <- makePSOCKcluster(num_cores)
  registerDoParallel(cl)
  score <-foreach (i=1:length(sublists), .combine='c', .packages=c("GOSemSim"), .export=c("goGoD")) %dopar% {
    library(GO.db)
    goGoD(sublists[i], GO_terms, "BP", "Wang", file_choice)
  }
  stopCluster(cl)
  plotcolor=length(score)
  names(score)->ll
  # sort the scores, largest to smallest
  score = sort(score,decreasing=TRUE)
  write.table(score,paste(measure,"score.txt"),row.names=TRUE, col.names = FALSE)
  # Create histograms in serial
  barplot(score,  cex.axis = 0.7, cex.names=0.5,
          names.arg=ll, col=heat.colors(plotcolor), main=paste0("ONTOLOGY BASED PRIORITIZATION OF PROBES \n", measure, " method"), xlab="Probe", ylab="Probe-GOA List Similarity")
  pdf( paste("Graphic score", measure))
  barplot(score,  cex.axis = 0.7, cex.names=0.5,
          names.arg=ll, col=heat.colors(plotcolor), main=paste0("ONTOLOGY BASED PRIORITIZATION OF PROBES \n", measure, " method"), xlab="Probe", ylab="Probe-GOA List Similarity")
  dev.off();
  runtime = proc.time() - ptm
  r2=matrix(runtime)
  cat(r2[3, 1])
}

# 6 probes: "AM_10430","AM_10198","AM_10192", "AM_10229", "AM_10013", "AM_10339"
# 8 probes: "AM_10430","AM_10198","AM_10192", "AM_10229", "AM_10013", "AM_10339", "AM_10031", "AM_10179"
probes=c("AM_10430","AM_10198","AM_10192", "AM_10229" , "AM_10013", "AM_10339", "AM_10031", "AM_10179")
terms=c("GO:0015849","GO:0055085","GO:0005886","GO:0016021","GO:0015143","GO:0005886","GO:0016324","GO:0042493", "GO:0070330", "GO:0016098", "GO:0015347", "GO:0015711",
        "GO:0070062", "GO:0060316", "GO:0060315", "GO:0070062", "GO:0004364")
measure="Wang"

#time_to_run=100
# for (i in 1:2){
#   total_runtime = total_runtime + parallel_GOGOD(probes, terms, "BP", "Wang", 4, "data.txt")
# }
#print(total_runtime)
parallel_GOGOD(probes, terms, "BP", "Wang", 2, "data.txt")
