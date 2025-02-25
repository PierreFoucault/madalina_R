parameters <- commandArgs(trailingOnly=TRUE)
fasta = readLines(parameters[1])
destino = parameters[2]
#fasta = readLines("/home/jesus/Escritorio/Madalina/Analisis_Madaline_All/Pseudomonas/Databases/LPSN/Pseudomonas_seqs.fasta")

ind_headers = grep(x=fasta,pattern = ">") #Vector containing the index of the headers

to_remove = c() #Vector that will contain the lines that will be removed after
for(i in 1:length(ind_headers)){ #As many loops as the number of headers
  n_of_fndgs = grep(x=fasta[ind_headers],pattern = fasta[ind_headers[i]]) #vector containing the index of the header being looped
  if(length(n_of_fndgs) == 2){ #If the header is found 2 times, we proced to:
     equal = identical(fasta[ind_headers[n_of_fndgs[1]]+1], fasta[ind_headers[n_of_fndgs[2]]+1]) #test if the sequences associated with those headers are also identicals
     if(equal){ #If so we:
       to_remove = c(to_remove, ind_headers[n_of_fndgs[2]]) #Add the index of the second header to be removed
     }
  }
}

#All values inside to_remove are repeated since they are found twice. We remove these lines (which are the headers of the repeated sequences) but also the next line which
#are the sequences itself. 
to_remove = unique(to_remove) #since all are duplicated, we sumplify it by keeping only one time each value.
fasta = fasta[-c(to_remove,to_remove+1)] #We remove the repeated headers and sequences so we only have it once in the fasta.

write.table(fasta,destino,quote = F,col.names=F, row.names=F, sep ="")
