#### Import data ####

## data from QIIME2
madalina <- qza_to_phyloseq(features = '16S_data/rar_filtered_table.qza',
                            tree='16S_data/rooted_tree.qza',
                            metadata = "16S_data/metadata_complete.tsv",
                            taxonomy = '16S_data/rar_filtered_taxonomy.qza')

## samples to remove (metadata_R script)
madalina<-madalina %>% subset_samples(to_remove !='YES')

## update metadata
madalina@sam_data<-
  read_delim("16S_data/metadata_complete.tsv","\t", escape_double = FALSE,
             trim_ws = TRUE,show_col_types = FALSE) %>%
  remove_rownames %>%
  column_to_rownames(var="sample-id") %>%
  subset(to_remove !='YES') %>%
  dplyr::mutate(across(c('week'), substr, 2, nchar(week))) %>%
  
  dplyr::mutate(.,
                samples=rownames(.),
                month=factor(month,levels=c("Feb.","Mar.","Apr.","May","Jun.","Jul.")),
                week=as.numeric(week),
                site=paste0(UDI," ",status),
                site_week=paste0(site," ",week),
                site=factor(site,levels=c("Marne pre-station","Marne post-station",
                                          "Vanne pre-station","Vanne post-station")),
                status=factor(status,levels=c("pre-station","post-station"))) %>%
  sample_data()

## set new ASV names for plots
ASV_names<-data.frame(features=rownames(madalina@otu_table),
                      ASV_ID=paste("ASV", 1:ntaxa(madalina), sep = "_"))

####| ####

#### Rar. curve ####

# get rarefaction threshold
raremax=paste('rar. at',sum(madalina@otu_table[,1]),'reads')

# Plot the rar. curve
Fig_rar_curve_16S<-
  rarecurve(otu_table(madalina) %>% as.data.frame() %>% t(), step=500, cex=0.5,tidy = T) %>%
  ggplot(.,aes(x=Sample,y=Species,group=Site))+
  geom_line(show.legend = F,color="black",linewidth=0.2)+theme_bw()+
  labs( y = "ASVs Richness", x = "Nb. of filtered reads\n(16S rRNA gene amplicon)")+
  theme(panel.grid = element_blank(),
        plot.title = element_blank(),
        axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        axis.ticks = element_blank())+
  #geom_vline(xintercept=raremax_B,color="blue")+
  #geom_text(x=2000, y=2000,label = raremax,check_overlap = T,size = 12,color="black")+
  annotate('text', x=3000, y=1800,label=raremax,size=10,hjust=0.5)+
  #scale_color_manual(values = c(palette_lake_3T,))+
  scale_x_continuous(expand=c(0,0))+ 
  scale_y_continuous(expand=c(0,NA))
Fig_rar_curve_16S

ggsave("/home/jesus/Pierre/madalina//Figures/rarecurve_16S.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)


####| ####

#### Barplot ####
#View(madalina@tax_table)

####____Phylum ####
madalina.barplot.df <-
  madalina %>% tax_glom(.,taxrank ="Phylum") %>%
  psmelt(.) %>% rename("ASV"="OTU") %>%
  # dplyr::mutate(.,
  #               Phylum_extended=ifelse(Phylum=="Proteobacteria",Class,Phylum)) %>%
  dplyr::group_by(site_week,Phylum,.drop = FALSE,.add = TRUE) %>%
  summarise(median_count=median(Abundance),.groups = "keep") %>%
  dplyr::group_by(site_week) %>%
  dplyr::mutate(median_abundance=
                  as.numeric(paste0((round(median_count/sum(median_count),4))))) %>%
  dplyr::mutate(.,
                Phylum_legend=ifelse(median_abundance<0.05,"Phyla < 5%",Phylum)) %>%
  dplyr::group_by(site_week,Phylum_legend,.drop = FALSE,.add = TRUE) %>%
  summarise(median_abundance=sum(median_abundance),.groups = "keep") %>%
  cbind(UDI=
          (madalina@sam_data$UDI[match(.$site_week,madalina@sam_data$site_week)]),
        site=
          (madalina@sam_data$site[match(.$site_week,madalina@sam_data$site_week)]),
        month=
          (madalina@sam_data$month[match(.$site_week,madalina@sam_data$site_week)]),
        week=
          (madalina@sam_data$week[match(.$site_week,madalina@sam_data$site_week)]),
        status=
          (madalina@sam_data$status[match(.$site_week,madalina@sam_data$site_week)]),
        .)
  

Phylum.barplot <-
  madalina.barplot.df %>%
  ggplot(.,aes(x=as.factor(week), y=(median_abundance),
               fill=fct_rev(fct_reorder(Phylum_legend,median_abundance))))+
  geom_col(width=0.8,color="black",size=0.2,show.legend = T)+
  theme_light()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color="black"),
        axis.title.y = element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.x=element_text(size=12,hjust=0.5),
                   axis.text.y=element_text(size=12),
                   axis.ticks =element_blank(),
                   legend.title = element_text(size=15,face = "bold"),
                   legend.position = "bottom",
                   legend.direction = "vertical",
                   legend.text = element_text(size=12))+
  scale_x_discrete(expand = c(0,0))+
  
  #palette_season<-c("#63849B","#EBAF47","#80A53F","#E36414")


  scale_y_continuous(expand = c(0,0),label=label_percent(accuracy = NULL,scale = 100))+
  coord_cartesian(xlim = c(0,12),ylim = c(0,1),clip = "off")+
  annotate("segment", x = 0.75, xend = 1.25, y = -0.16, yend = -0.16,color = "#63849B",linewidth=3)+
  annotate("text", x = 1.1,y = -0.27,colour = "black",label="Feb.",size=4,hjust=0.5)+
  annotate("segment",x = 1.75,xend = 3.25, y = -0.16, yend = -0.16,color = "#EBAF47",linewidth=3)+
  annotate("text", x = 2.6,y = -0.27,colour = "black",label="Mar.",size=4,hjust=0.5)+
  annotate("segment",x = 3.75,xend = 4.25, y = -0.16, yend = -0.16,color = "#EBAF47",linewidth=3)+
  annotate("text", x = 4.1,y = -0.27,colour = "black",label="Apr.",size=4,hjust=0.5)+
  annotate("segment",x = 4.75,xend = 7.25, y = -0.16, yend = -0.16,color = "#EBAF47",linewidth=3)+
  annotate("text", x = 6,y = -0.27,colour = "black",label="May",size=4,hjust=0.5)+
  annotate("segment",x = 7.75,xend = 9.25, y = -0.16, yend = -0.16,color = "#80A53F",linewidth=3)+
  annotate("text", x = 8.6,y = -0.27,colour = "black",label="Jun.",size=4,hjust=0.5)+
  annotate("segment",x = 9.75,xend = 11.25, y = -0.16, yend = -0.16,color = "#80A53F",linewidth=3)+
  annotate("text", x = 10.6,y = -0.27,colour = "black",label="Jul.",size=4,hjust=0.5)+
  scale_fill_manual(values=c("#FED9A6","#FFFFCC","#B3CDE3","#DECBE4","#CCEBC5",
                             "#CBD5E8","firebrick","#FDDAEC","#E5D8BD","violet","#F2F2F2",
                             "black"))+
  guides(fill=guide_legend(nrow=3,shape=22,size=5))+
  labs(fill="Prokaryotic Phyla")+
  facet_wrap2(~ site,scales = "fixed",
              strip = strip_themed(
                background_x = elem_list_rect(fill = palette_code,color = c("black","black","black","black")),
                text_x = elem_list_text(colour = "white",face = "bold",size=10)))

ggsave("/home/jesus/Pierre/madalina//Figures/barplot_16S.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####____Genus ####
madalina_genus.barplot.df <-
  madalina %>% #tax_glom(.,taxrank ="Genus") %>%
  psmelt(.) %>% rename("ASV"="OTU") %>%
  # dplyr::mutate(.,
  #               Phylum_extended=ifelse(Phylum=="Proteobacteria",Class,Phylum)) %>%
  dplyr::group_by(site_week,Genus,Species,ASV,.drop = FALSE,.add = TRUE) %>%
  #dplyr::group_by(site_week,Phylum,.drop = FALSE,.add = TRUE) %>%
  summarise(median_count=median(Abundance),.groups = "keep") %>%
  dplyr::group_by(site_week) %>%
  mutate(median_abundance=as.numeric(
    paste0((round(median_count/sum(median_count),4)))))

View(madalina_genus.barplot.df)

####______Legionella ####
legionella.barplot.df <-
  madalina_genus.barplot.df %>%
  subset(Genus == "Legionella") %>%
  cbind(ASV_ID=
           ASV_names$ASV_ID[match(.$ASV,ASV_names$features)]) %>%
  dplyr::mutate(.,
                Genus_legend=
                  ifelse(Species =="Legionella_pneumophila",paste0("L. pneumophila ",ASV_ID),"Other Legionella ASVs"),
                Genus_legend = coalesce(Genus_legend,"Other Legionella ASVs")) %>%
  dplyr::group_by(site_week,Genus_legend,.drop = FALSE,.add = TRUE) %>%
  summarise(median_abundance=sum(median_abundance),.groups = "keep") %>%
  cbind(UDI=
          madalina@sam_data$UDI[match(.$site_week,madalina@sam_data$site_week)],
        site=
          madalina@sam_data$site[match(.$site_week,madalina@sam_data$site_week)],
        month=
          madalina@sam_data$month[match(.$site_week,madalina@sam_data$site_week)],
        week=
          madalina@sam_data$week[match(.$site_week,madalina@sam_data$site_week)],
        status=
          madalina@sam_data$status[match(.$site_week,madalina@sam_data$site_week)],.) %>%
  dplyr::filter(median_abundance>0) %>%
  dplyr::mutate(.,
                Genus_legend=
                  ifelse(Genus_legend=="L. pneumophila ASV_16380","L. pneumophila ASV_1",
                         ifelse(Genus_legend=="L. pneumophila ASV_16378","L. pneumophila ASV_2",
                         Genus_legend)))
  

View(legionella.barplot.df)
legionella.barplot<-legionella.barplot.df %>%
  ggplot(.,aes(x=as.factor(week), y=(median_abundance),
               fill=Genus_legend))+
               #fill=fct_rev(fct_reorder(Phylum_legend,median_abundance))))+
  geom_col(width=0.8,color="black",linewidth=0.2,show.legend = T)+
  theme_light()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color="black"),
        axis.title.y = element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.x=element_text(size=12,hjust=0.5),
                   axis.text.y=element_text(size=12),
                   axis.ticks =element_blank(),
                   legend.title = element_text(size=15,face = "bold.italic"),
                   legend.position = "bottom",
                   legend.direction = "vertical",
                   legend.text =element_markdown(size=12))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0),label=label_percent())+
  scale_fill_manual(values=c("#C3D6C4","#87AE89","#0D5b11"),
                    labels = c("*L. pneumophila* ASV_1","*L. pneumophila* ASV_2","Other *Legionella* ASVs"))+
  guides(fill=guide_legend(ncol=1,shape=22,size=5))+
  labs(fill="Legionella ASVs")+
  coord_cartesian(xlim = c(0,12),ylim = c(0,0.04),clip = "off")+
  # annotate("rect", xmin = 0.75, xmax = 1.25, ymin = -0.13, ymax = -0.2,fill = "#63849B",color="black")+
  # annotate("text", x = 1.1,y = -0.27,colour = "black",label="Feb.",size=4,hjust=0.5)+
  # annotate("rect", xmin = 1.75, xmax = 3.25, ymin = -0.13, ymax = -0.2,fill = "#EBAF47",color="black")+
  # annotate("text", x = 2.7,y = -0.27,colour = "black",label="Mar.",size=4,hjust=0.5)+
  # annotate("rect", xmin = 3.75, xmax = 4.25, ymin = -0.13, ymax = -0.2,fill = "#EBAF47",color="black")+
  # annotate("text", x = 4.2,y = -0.27,colour = "black",label="Apr.",size=4,hjust=0.5)+
  # annotate("rect", xmin = 4.75, xmax = 7.25, ymin = -0.13, ymax = -0.2,fill = "#EBAF47",color="black")+
  # annotate("text", x = 6.1,y = -0.27,colour = "black",label="May",size=4,hjust=0.5)+
  # annotate("rect", xmin = 7.75, xmax = 9.25, ymin = -0.13, ymax = -0.2,fill = "#80A53F",color="black")+
  # annotate("text", x = 8.6,y = -0.27,colour = "black",label="Jun.",size=4,hjust=0.5)+
  # annotate("rect", xmin = 9.75, xmax = 11.25, ymin = -0.13, ymax = -0.2,fill = "#80A53F",color="black")+
  # annotate("text", x = 10.6,y = -0.27,colour = "black",label="Jul.",size=4,hjust=0.5)+
  annotate("segment", x = 0.75, xend = 1.25, y = -0.0045, yend = -0.0045,colour = "#63849B",linewidth=3)+
  #annotate("text", x = 1.1,y = -0.0045,colour = "black",label="Feb.",size=4,hjust=0.5)+
  annotate("segment", x = 1.75, xend = 3.25, y = -0.0045, yend = -0.0045,colour = "#EBAF47",linewidth=3)+
  #annotate("text", x = 2.7,y = -0.0088,colour = "black",label="Mar.",size=4,hjust=0.5)+
  annotate("segment", x = 3.75, xend = 4.25, y = -0.0045, yend = -0.0045,colour = "#EBAF47",linewidth=3)+
  #annotate("text", x = 4.2,y = -0.0088,colour = "black",label="Apr.",size=4,hjust=0.5)+
  annotate("segment", x = 4.75, xend = 7.25, y = -0.0045, yend = -0.0045,colour = "#EBAF47",linewidth=3)+
  #annotate("text", x = 6.1,y = -0.0088,colour = "black",label="May",size=4,hjust=0.5)+
  annotate("segment", x = 7.75, xend = 9.25, y = -0.0045, yend = -0.0045,colour = "#80A53F",linewidth=3)+
  #annotate("text", x = 8.6,y = -0.0088,colour = "black",label="Jun.",size=4,hjust=0.5)+
  annotate("segment", x = 9.75, xend = 11.25, y = -0.0045, yend = -0.0045,colour = "#80A53F",linewidth=3)+
  #annotate("text", x = 10.6,y = -0.0088,colour = "black",label="Jul.",size=4,hjust=0.5)+
  facet_wrap2( ~ site,scales = "fixed",nrow = 1, ncol = 4,
              strip = strip_themed(
                background_x = elem_list_rect(fill = palette_code,color = c("black","black","black","black")),
                text_x = elem_list_text(colour = "white",face = "bold",size=10)))


ggsave("/home/jesus/Pierre/madalina//Figures/barplot_legio.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

Phylum.barplot+BC_16S+legionella.barplot+
  plot_layout(design=design_diversity,guides = 'collect')+
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size=15,face="bold"),
        legend.position = "bottom",legend.direction = "vertical")

ggsave("/home/jesus/Pierre/madalina/Figures/diversity.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

#### Beta-Diversity ####

####_____Bray_Curtis ####
madalina.BC<-ordinate(madalina, "PCoA","bray") 

madalina_BC.df <- plot_ordination(madalina,madalina.BC,color="UDI",shape="status",label="month") %>%
  .[[1]]

ggplot(madalina_BC.df,aes(Axis.1,Axis.2))+
  # geom_convexhull(aes(group = season),
  #                 alpha=0.6,show.legend = F,size=0.6)+
  geom_point(size=5,shape=21,fill=site)+
             #aes(fill=UDI,shape=status))+
  theme_bw() +
  theme(#aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.title = element_text(size=15,face = "bold"),
        axis.text = element_text(size=10),
        axis.ticks = element_blank(),
        plot.title = element_text(size=20,face = "bold"),
        legend.title = element_text(size=15,face = "bold",hjust=0),
        legend.text = element_text(size=12),
        legend.position = "bottom",legend.direction = "vertical")+
  #scale_color_manual(values=palette_season)+
  scale_fill_manual(values=palette_UDI)+
  scale_shape_manual(values = c(21,24))+
  #scale_x_reverse()+
  #cale_y_reverse()+
  guides(fill=guide_legend(override.aes=aes(alpha=1,shape=21,size=6),color="black",nrow=1),
         shape=guide_legend(override.aes=aes(size=5),nrow=1))+
  labs(fill="UDI",shape='Status',
       title="Procaryotic community structure",
       x=paste0("Axis 1 [",round(madalina.BC$values$Relative_eig[1],3)*100,"%]"),
       y=paste0("Axis 2 [",round(madalina.BC$values$Relative_eig[2],3)*100,"%]"))

uniquemonth <- c("2","3","4","5","6","7")
labelmonth<- unlist(lapply(uniquemonth, utf8ToInt))

ggplot(madalina_BC.df,aes(Axis.1,Axis.2,color=site,fill=site,shape = month))+
  # geom_convexhull(aes(group = season),
  #                 alpha=0.6,show.legend = F,size=0.6)+
  geom_point(size=7)+
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.title = element_text(size=15,face = "bold"),
        axis.text = element_text(size=10),
        axis.ticks = element_blank(),
        plot.title = element_text(size=20,face = "bold"),
        legend.title = element_text(size=15,face = "bold",hjust=0),
        legend.text = element_text(size=12),
        legend.position = "bottom",legend.direction = "vertical")+
  scale_fill_manual(values=palette_code)+
  scale_color_manual(values=palette_code)+
  scale_shape_manual(values=labelmonth)+
  scale_x_reverse()+
  guides(fill=guide_legend(override.aes=aes(alpha=1,shape=21,size=6),color="black",nrow=2),
         color=FALSE,shape=guide_legend(override.aes=aes(size=6),nrow=2))+
  labs(fill="UDI",shape="Month",
       title="Procaryotic community structure",
       x=paste0("Axis 1 [",round(madalina.BC$values$Relative_eig[1],3)*100,"%]"),
       y=paste0("Axis 2 [",round(madalina.BC$values$Relative_eig[2],3)*100,"%]"))

BC_16S<-
  ggplot(madalina_BC.df,aes(Axis.1,Axis.2,
                            fill=site,#shape=status,
                            label=month_number))+
  geom_point(size=5,shape=21)+
  geom_text_repel(show.legend = F, max.overlaps = 30,size=5,
                  min.segment.length = Inf,force_pull = 10)+
  theme_light() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color="black"),
        axis.title = element_text(size=15,face = "bold"),
        axis.text = element_text(size=10),
        axis.ticks = element_blank(),
        legend.title = element_text(size=15,face = "bold",hjust=0),
        legend.text = element_text(size=12),
        legend.position = "bottom",legend.direction = "vertical")+
  scale_fill_manual(values=palette_code)+
  #scale_shape_manual(values = c(21,24))+
  scale_x_reverse()+
  guides(fill=guide_legend(override.aes=aes(alpha=1,shape=21,size=6),
                           color="black",nrow=2),
         shape=guide_legend(override.aes=aes(size=5),nrow=1))+
  labs(fill="UDI",shape='Status',
       x=paste0("Axis 1 [",round(madalina.BC$values$Relative_eig[1],3)*100,"%]"),
       y=paste0("Axis 2 [",round(madalina.BC$values$Relative_eig[2],3)*100,"%]"))

ggsave("/home/jesus/Pierre/madalina/Figures/BC_16S.pdf",units = "in",dpi = "retina",width = 13.4,height = 9.8)

####_____Jaccard ####

####____________####