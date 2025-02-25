#### Download all fasta for a Genus of interest from LPSN (https://lpsn.dsmz.de/) ####

# needed packages
library(rvest)
library(httr)

# get arguments for fetching function

#ask user for arguments
arguments <- commandArgs(trailingOnly=TRUE)

if(any(arguments=="-h") | any(arguments=="--help")){cat(" 1: Genus of interest (Ex: legionella or pseudomonas) \n 2: Ouptut path to a directory" ); stop("")}

genus <- arguments[1]
download_dir <- arguments[2]

# URL of the LPSN to scrape
lpsn_url <- 'https://lpsn.dsmz.de/' #url for the website
genus_url <- paste("https://lpsn.dsmz.de/search?word=",genus, sep="") #url for the Genus of interest

# Check if output filepath exists and create it if needed
if (!dir.exists(download_dir)) {
  dir.create(download_dir)
}

# Send a request to the LPSN website and list all possible urls for the Genus of interest
html_content <- read_html(genus_url) %>% html_nodes('a') %>% html_attr('href')

# Get a list of all species belonging to the Genus of interest
pattern = paste("/species/",genus,sep="")
html_all_species = html_content[grep(x=html_content,pattern = pattern)]
species.name = unlist(lapply(strsplit(x=html_all_species, split = "/"), function(x) x[length(x)]))

# Reset the genus_url as a list of all species urls for fasta dowloading
genus_url = paste(lpsn_url,html_all_species,sep="")

# Final fasta downloading loop
for(i in 1:length(genus_url)){  
  #Preparing the link
  html_content = read_html(genus_url[i]) %>% html_nodes('a') %>% html_attr('href')
 # %>% html_attr(html_nodes(.,"a"),'href')
  link = html_content[grep(x = html_content, pattern=(".fasta"))]
  if(length(link) == 0) {next}
  full_link = paste(lpsn_url,link, sep="")
  #Downloading the associated fasta
  file_name <- file.path(download_dir, basename(link))
  print(paste('Downloading', species.name[i], 'to', file_name))
  # Download the file
  GET(full_link, write_disk(file_name, overwrite = TRUE))
}
