##Read FASTA files

library("seqinr")
leprae <- read.fasta(file = "leprae.fasta")
lepraeseq <- leprae[[1]]
length(lepraeseq)


