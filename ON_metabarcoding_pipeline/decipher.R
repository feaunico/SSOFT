library(DECIPHER)

# specify the path to the FASTA file (in quotes)
fas <- "racon.fasta"

# load the sequences from the file
seqs <- readDNAStringSet(fas) # or readRNAStringSet

# remove any gaps (if needed)
seqs <- RemoveGaps(seqs)

# load a training set object (trainingSet)
# see http://DECIPHER.codes/Downloads.html
load("UNITE_v2021_May2021.RData")

# classify the sequences
ids <- IdTaxa(seqs,
   trainingSet,
   strand="both", # or "top" if same as trainingSet
   threshold=60, # 60 (cautious) or 50 (sensible)
   processors=NULL) # use all available processors

# look at the results
print(ids)

assignment <- sapply(ids,
function(x)
paste(x$confidence,
collapse=";"))
write.table(assignment, "decipher_out", append = TRUE, sep = " ")
assignment <- sapply(ids,
function(x)
paste(x$taxon,
collapse=";"))
write.table(assignment, "decipher_out", append = TRUE, sep = " ")