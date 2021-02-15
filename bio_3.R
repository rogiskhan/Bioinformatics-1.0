library("seqinr")
choosebank()
choosebank("genbank") # Specify that we want to search the 'genbank' ACNUC sub-database
choosebank("refseq") # Specify that we want to search the 'refseq' ACNUC sub-database
#then you use query() 
choosebank("genbank")
query("SchistosomamRNA", "SP=Schistosoma mansoni AND M=mrna")
closebank()
#example #2
choosebank("refseqViruses")
dengue1 <- query("Dengue1", "AC=NC_001477")
attributes(dengue1)
dengue1$nelem
dengue1$req
attr(dengue1, "names")

#retrieve sequence data
dengueseq <- getSequence(dengue1$req[[1]])

#retrieve annotation data
annots <- getAnnot(dengue1$req[[1]])

#retrieve data from articles
choosebank("genbank") # Specify that we want to search the 'genbank' ACNUC sub-database
query('naturepaper', 'R=Nature/460/352')



