## Parallel gene annotation implementation in R
## Download annotation set and retrieve unique GO terms from it
## QuickGO database of EMBL-EBI used

## To do:
## - Decide how to represent genes for input (currently have vector of gene product ID's given)

library(parallel)
library(RCurl)

start_time <- proc.time()

## Return an annotation (vector of GO terms) for a protein.
DownloadAnnotation <- function(proteinID){
  ## Build up the URL string
  protein_URL <- paste("http://www.ebi.ac.uk/QuickGO/GAnnotation?protein=", proteinID, sep = "")
  protein_URL <- paste(protein_URL, "&format=tsv", sep = "")
  
  ## URL for one protein.
  geneURL <- getURL(protein_URL)
  ## Read in the data from the URL. 
  data <- read.csv(text = geneURL, header = TRUE, sep = "\t")
  ## Find and get the column that contains the GO IDs.
  ID_column_number <- which(names(data) == "GO.ID")
  ## Create sorted vector of terms, with duplicates removed via "unique" function
  ## [, ID_column_number] used because data is data frame, and the blank lets us get the whole column as a vector.
  annotation <- as.vector(unique(data[ , ID_column_number]))
  # Print unique terms
  #   for (term in annotation){
  #     print(term)
  #   }
  return (annotation)
}

## Return list of gene annotations
GetAnnotations <- function(genes){
  ## Download annotations of list of genes
  ## Parallel: Split up work of annotation downloading among processors
  annotations = mclapply(genes, DownloadAnnotation)
#   annotations = numeric(length(genes))
#   for (i in 1:length(genes)){
#     annotations[i]=DownloadAnnotation(genes[i])
#   }
  return annotations
}

## Sample annotation downloads
protein = "2950"
#proteins <- rep(proteins, times = 500)
#print(length(proteins))
#annots = GetAnnotations(proteins)
ann = DownloadAnnotation(protein)
print(ann)
## Timings
end_time <- (proc.time() - start_time)
print(end_time)