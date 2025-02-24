#### Import packages ####

#core
library(BiocManager)
library(tidyverse)

#visual
library(tidytext)
library(ggh4x)
library(patchwork)
library(reshape2)
library(ggrepel)
library(ggConvexHull)
library(scales)
library(ggtext)

#analysis
library(phyloseq)
library(qiime2R)
library(vegan)

#stats
library(rstatix)
library(lmerTest)
library(lme4)

#### Palettes #####

palette_code<-c('#694065','#CB9EC6','#536102','#ABC705') 

palette_UDI<-c('#CB9EC6','#ABC705')

palette_easter<-c("#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2","#B3E2CD","#FDCDAC",
                  "#CBD5E8","#F4CAE4","#E6F5C9","#FFF2AE","#F1E2CC","#88F393","#2E6D74","#DCBCF3","#EFE0DC","#F18B6C",
                  "#B52239","#F5D2DB","#B7F7E6","#517BE6","#E0DBC1","#E19B8E","#BDCEC4","#99F8F3","#E45632","#D22A82",
                  "#70F4BB","#5D9BEF","#F391CC","#EEF4A3","#CAAB7D","#8AF3D4","#EAF7B9","#4263B2","#E0B333","#2CB550",
                  "#DA56EF","#EC97B3","#E0789F","#8CB8E2","#E89E34","#8D37A3","#71992A","#D8C5B2","#E75FA6","#EBE49A",
                  "#ADAC5A","#BDF6D3","#3BDCA8","#47D1C8","#F8EBD4","#F4F6D8","#B6A7BD","#D6F0F2","#3EA34E","#B8E2EC",
                  "#BE6AE7","#8B7664","#7161EC","#CACCA8","#F138F4","#8EAC73","#A899F0")

pie(rep(1,length(palette_UDI)), col=(palette_UDI))

#install.packages('')
#BiocManager::install(')
#devtools::install_github("")

# if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
# devtools::install_github("jbisanz/qiime2R")
#devtools::install_github("cmartin/ggConvexHull")
