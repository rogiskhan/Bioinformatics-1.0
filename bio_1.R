##Read FASTA files

library("seqinr")
leprae <- read.fasta(file = "leprae.fasta")
lepraeseq <- leprae[[1]]

## length 
length(lepraeseq)

#How many of each of the four nucleotides 
#A, C, T and G, and any other symbols, are
#there in the Mycobacterium leprae TN genome sequence?

t <- table(lepraeseq)
print(t)

#What is the GC content of the Mycobacterium leprae TN genome sequence,
#when (i) all non-A/C/T/G nucleotides are included, 
#(ii) non-A/C/T/G nucleotides are discarded?

#(i)
percentage_with <- GC(lepraeseq)*100

#(ii)
percentage_without <- GC(lepraeseq, exact = FALSE)*100


