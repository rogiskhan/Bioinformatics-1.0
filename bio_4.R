library("seqinr")
leprae <- read.fasta(file = "Q9CD83.fasta")
ulcerans <- read.fasta(file = "A0PQ23.fasta")
lepraeseq <- leprae[[1]]
ulceransseq <- ulcerans[[1]]

#OR... not working 

library("seqinr")
choosebank("swissprot")
query("leprae", "AC=Q9CD83")
lepraeseq <- getSequence(leprae$req[[1]])
query("ulcerans", "AC=A0PQ23")
ulceransseq <- getSequence(ulcerans$req[[1]])
closebank()

###############################################
#Find similarities
dotPlot(lepraeseq, ulceransseq)

###############################################
#ALLIGNMENT
library(Biostrings)
sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
###############

s1 <- "GAATTC"
s2 <- "GATTA"
globalAligns1s2 <- pairwiseAlignment(s1, s2, substitutionMatrix = sigma, gapOpening = -2,
                                       gapExtension = -8, scoreOnly = FALSE)
globalAligns1s2 # Print out the optimal alignment and its score
##################
data(BLOSUM50)
s3 <- "PAWHEAE"
s4 <- "HEAGAWGHEE"
globalAligns3s4 <- pairwiseAlignment(s3, s4, substitutionMatrix = "BLOSUM50", gapOpening = -2,
                                     gapExtension = -8, scoreOnly = FALSE)
globalAligns3s4 # Print out the optimal global alignment and its score


##########################
#convert list to string
lepraeseqstring <- c2s(lepraeseq)     # Make a string that contains the sequence in "lepraeseq"
ulceransseqstring <- c2s(ulceransseq) # Make a string that contains the sequence in "ulceransseq"
lepraeseqstring <- toupper(lepraeseqstring)
ulceransseqstring <- toupper(ulceransseqstring)
globalAlignLepraeUlcerans <- pairwiseAlignment(lepraeseqstring, ulceransseqstring,
                                               substitutionMatrix = BLOSUM50, gapOpening = -2, gapExtension = -8, scoreOnly = FALSE)
globalAlignLepraeUlcerans
#######################

printPairwiseAlignment <- function(alignment, chunksize=60, returnlist=FALSE)
{
  require(Biostrings)           # This function requires the Biostrings package
  seq1aln <- pattern(alignment) # Get the alignment for the first sequence
  seq2aln <- subject(alignment) # Get the alignment for the second sequence
  alnlen  <- nchar(seq1aln)     # Find the number of columns in the alignment
  starts  <- seq(1, alnlen, by=chunksize)
  n       <- length(starts)
  seq1alnresidues <- 0
  seq2alnresidues <- 0
  for (i in 1:n) {
    chunkseq1aln <- substring(seq1aln, starts[i], starts[i]+chunksize-1)
    chunkseq2aln <- substring(seq2aln, starts[i], starts[i]+chunksize-1)
    # Find out how many gaps there are in chunkseq1aln:
    gaps1 <- countPattern("-",chunkseq1aln) # countPattern() is from Biostrings package
    # Find out how many gaps there are in chunkseq2aln:
    gaps2 <- countPattern("-",chunkseq2aln) # countPattern() is from Biostrings package
    # Calculate how many residues of the first sequence we have printed so far in the alignment:
    seq1alnresidues <- seq1alnresidues + chunksize - gaps1
    # Calculate how many residues of the second sequence we have printed so far in the alignment:
    seq2alnresidues <- seq2alnresidues + chunksize - gaps2
    if (returnlist == 'FALSE')
    {
      print(paste(chunkseq1aln,seq1alnresidues))
      print(paste(chunkseq2aln,seq2alnresidues))
      print(paste(' '))
    }
  }
  if (returnlist == 'TRUE')
  {
    vector1 <- s2c(substring(seq1aln, 1, nchar(seq1aln)))
    vector2 <- s2c(substring(seq2aln, 1, nchar(seq2aln)))
    mylist <- list(vector1, vector2)
    return(mylist)
  }
}
printPairwiseAlignment(globalAlignLepraeUlcerans, 60)
########################## Smith-Waterman
localAlignLepraeUlcerans <- pairwiseAlignment(lepraeseqstring, ulceransseqstring,
                                              substitutionMatrix = BLOSUM50, gapOpening = -2, gapExtension = -8, scoreOnly = FALSE, type="local")
localAlignLepraeUlcerans # Print out the optimal local alignment and its score

printPairwiseAlignment(localAlignLepraeUlcerans, 60)
#########################################

generateSeqsWithMultinomialModel <- function(inputsequence, X)
{
  # Change the input sequence into a vector of letters
  require("seqinr") # This function requires the SeqinR package.
  inputsequencevector <- s2c(inputsequence)
  # Find the frequencies of the letters in the input sequence "inputsequencevector":
  mylength <- length(inputsequencevector)
  mytable <- table(inputsequencevector)
  # Find the names of the letters in the sequence
  letters <- rownames(mytable)
  numletters <- length(letters)
  probabilities <- numeric() # Make a vector to store the probabilities of letters
  for (i in 1:numletters)
  {
    letter <- letters[i]
    count <- mytable[[i]]
    probabilities[i] <- count/mylength
  }
  # Make X random sequences using the multinomial model with probabilities "probabilities"
  seqs <- numeric(X)
  for (j in 1:X)
  {
    seq <- sample(letters, mylength, rep=TRUE, prob=probabilities) # Sample with replacement
    seq <- c2s(seq)
    seqs[j] <- seq
  }
  # Return the vector of random sequences
  return(seqs)
}

randomseqs <- generateSeqsWithMultinomialModel('PAWHEAE',1000)
randomseqs[1:10]

s4 <- "HEAGAWGHEE"
pairwiseAlignment(s4, randomseqs[1], substitutionMatrix = "BLOSUM50", gapOpening = -2,
                  gapExtension = -8, scoreOnly = FALSE)
randomscores <- double(1000) # Create a numeric vector with 1000 elements
for (i in 1:1000)
{
  score <- pairwiseAlignment(s4, randomseqs[i], substitutionMatrix = "BLOSUM50",
                             gapOpening = -2, gapExtension = -8, scoreOnly = TRUE)
  randomscores[i] <- score
}

hist(randomscores, col="red") # Draw a red histogram
sum(randomscores >= -5)
