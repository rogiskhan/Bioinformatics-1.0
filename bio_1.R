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

#How many of each of the four nucleotides A, C, T and G are there in the complement of the Mycobacterium leprae TN genome sequence?

compleprae <- comp(lepraeseq)
t2 <- table(compleprae)
print(t2)

#How many occurrences of the DNA words CC, CG and GC occur in the Mycobacterium leprae TN genome sequence?

count(lepraeseq, 2)

#How many occurrences of the DNA words CC, CG and GC occur in the (i) first 1000 and (ii) last 1000 nucleotides of the Mycobacterium leprae TN genome sequence?
count(lepraeseq[1:1000],2)
count(lepraeseq[length(lepraeseq)-1000:length(lepraeseq)],2)

















