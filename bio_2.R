myvector1 <- c(10, 15, 22, 35, 43)
myvector2 <- c(3, 3.2, 3.9, 4.1, 5.2)
plot(myvector1, myvector2, xlab="myvector1", ylab="myvector2",type='b')
myfunction <- function(x) {return(20 + (x*x))}
myfunction
library("seqinr")
dengue <- read.fasta(file = "den1.fasta.txt")
dengueseq <- dengue[[1]]
starts <- seq(1, length(dengueseq)-2000, by = 2000)
starts
n <- length(starts)# Find the length of the vector "starts
chunkGCs <- numeric(n)
for (i in 1:n) {
  chunk <- dengueseq[starts[i]:(starts[i]+1999)]
  chunkGC <- GC(chunk)*100
  print (chunkGC)
  chunkGCs[i] <- chunkGC
}
plot(starts,chunkGCs,type="b",xlab="Nucleotide start position",ylab="GC content")


slidingwindowplot <- function(windowsize, inputseq)
{
  starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
  n <- length(starts)
  chunkGCs <- numeric(n)
  for (i in 1:n) {
    chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)]
    chunkGC <- GC(chunk)
    chunkGCs[i] <- chunkGC
  }
  plot(starts,chunkGCs,type="b",xlab="Nucleotide start position",ylab="GC content")
}

slidingwindowplot(1000, dengueseq)

#Draw a sliding window plot of GC content in the DEN-1 Dengue virus genome, using a window size of 200 nucleotides. Do you see any regions of unusual DNA content in the genome (eg. a high peak or low trough)?

slidingwindowplot(2000, dengueseq)

#Draw a sliding window plot of GC content in the genome sequence for the bacterium Mycobacterium leprae strain TN (accession NC_002677) using a window size of 20000 nucleotides. Do you see any regions of unusual DNA content in the genome (eg. a high peak or low trough)?
slidingwindowplot(200000, lepraeseq)

#Write a function to calculate the AT content of a DNA sequence (ie. the fraction of the nucleotides in the sequence that are As or Ts). What is the AT content of the Mycobacterium leprae TN genome?

AT_calculator <- function(inputseq)
{
  mytable <- count(inputseq, 1)
  len <- length(inputseq)
  A <- mytable[[1]]
  T <- mytable[[4]]
  myAT <- (A + T)/len
  return(myAT*100)
}
AT_calculator(lepraeseq)






