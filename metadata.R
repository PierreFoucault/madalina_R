

metadata_16S<-read_delim("/home/jesus/Pierre/madalina/madalina_qiime2/metadata.txt",delim = "\t",show_col_types = F) %>%
  dplyr::mutate(to_remove=if_else(type %in% c("S66",'test'),"YES",
                                  if_else(location %in% c("CF",'CE','cf','contF'),"YES",
                                          if_else(week =="S12","YES","NO"))),
                month_number=if_else(week=="S1","2",
                              if_else(week %in% c("S2",'S3'),"3",
                                      if_else(week=="S4","4",
                                              if_else(week %in% c("S5",'S6','S7'),"5",
                                                      if_else(week %in% c("S8",'S9'),"6",
                                                              if_else(week %in% c("S10",'S11'),"7",NA)))))),
                month=if_else(week=="S1","Feb.",
                                     if_else(week %in% c("S2",'S3'),"Mar.",
                                             if_else(week=="S4","Apr.",
                                                     if_else(week %in% c("S5",'S6','S7'),"May",
                                                             if_else(week %in% c("S8",'S9'),"Jun.",
                                                                     if_else(week %in% c("S10",'S11'),"Jul.",NA)))))),
                code_location=if_else(location=="rivB","ABAM",
                                      if_else(location=="rivT","IREF",
                                              if_else(location=="souB","LHBRUT",
                                                      if_else(location %in% c("souT",'SouT'),"LHGREF",NA)))),
                UDI=if_else(code_location %in% c('ABAM','IREF'),"Marne",
                                                 if_else(code_location %in% c('LHBRUT','LHGREF'),"Vanne",NA)),
                status=if_else(code_location %in% c('ABAM','LHBRUT'),"pre-station",
                                            if_else(code_location %in% c('IREF','LHGREF'),"post-station",NA))) %>%
  write_tsv(.,"16S_data/metadata_complete.tsv")
                
View(metadata_16S)

tibble(metadata_16S$`sample-id`,colnames(madalina@otu_table))

dframe %>% 
  mutate(col1_dupe = !is.na(match(col1, col2))) ->
  dframe2
dframe2