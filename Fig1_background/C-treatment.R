library(readr)
library(dplyr)
library(ggplot2)
library(swimplot)
library(ggsci)
library(pals)
library(readxl)
library(ggpubr)
library(RColorBrewer)
library(scales)
library(ggpubr)
library(patchwork)
setwd("~/CASCADEpaper/paper/Fig1_background")
# version 1 ---------------------------------------------------------------
treatment_history <- read_excel("treatment_history.xlsx")
colnames(treatment_history) <- c("ID", "Treatment", 'Duration',"order")
temp <- treatment_history %>% group_by(ID) %>% mutate(time = cumsum(Duration))
temp$start <- 0
temp$start <- temp$time - temp$Duration
temp <- temp[-30,]
temp <- as.data.frame(temp)
temp$ID <- gsub("00", "", temp$ID)
temp$phenotype <- ifelse(temp$ID %in% c("CA27", "CA58"), "Mixed", ifelse(temp$ID %in% c("CA46", "CA90"), "NE", "AD"))
temp$ID <- factor(temp$ID, levels = rev(unique(temp$ID[order(temp$phenotype)])))



temp$treatment_ordered <- factor(paste(temp$Treatment, temp$ID), levels = paste(temp$Treatment, temp$ID)[order(temp$time)])


n <- 18
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))


colors <- setNames(object = kelly(18), unique(temp$Treatment))

temp$tile_start <- temp$start + temp$Duration/2 #START AT TILE CENTER


[1] "#7FC97F" "#BEAED4" "#FDC086" "#FFFF99" "#386CB0" "#F0027F" "#BF5B17" "#666666" "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E"
[14] "#E6AB02" "#A6761D" "#666666" "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A"
[27] "#FFFF99" "#B15928" "#FBB4AE" "#B3CDE3" "#CCEBC5" "#DECBE4" "#FED9A6" "#FFFFCC" "#E5D8BD" "#FDDAEC" "#F2F2F2" "#B3E2CD" "#FDCDAC"
[40] "#CBD5E8" "#F4CAE4" "#E6F5C9" "#FFF2AE" "#F1E2CC" "#CCCCCC" "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#A65628"
[53] "#F781BF" "#999999" "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494" "#B3B3B3" "#8DD3C7" "#FFFFB3" "#BEBADA"
[66] "#FB8072" "#80B1D3" "#FDB462" "#B3DE69" "#FCCDE5" "#D9D9D9" "#BC80BD" "#CCEBC5" "#FFED6F"
# after JL new comments on treatment histroy ------------------------------
library(cowplot)
library(tools)
treatment_history <- read_excel("Cascade_Prostate_Clinical Data _Dec23_JC.xlsx", sheet = 5)

points <- read_excel("Cascade_Prostate_Clinical Data _Dec23_JC.xlsx", sheet = 6)

deaths <- data.frame(
  stringsAsFactors = FALSE,
                          ID = c("CA27",
                                 "CA34","CA35","CA43","CA46","CA58","CA76",
                                 "CA83","CA90"),
                       values = c(68.5,60,
                                 175.0666667,71.7,7.833333333,117.7333333,
                                 41.3,12.06666667,64.4),
                   treatment = c("Death",
                                 "Death","Death","Death","Death","Death",
                                 "Death","Death","Death")
          )
deaths$values <- deaths$values+4
points <- rbind(points, deaths)
points$values[1] <- 5

temp <- merge(treatment_history, points, by = "ID", all = T)
temp$height <- 0.8* temp$height
temp$Treatment <- factor(temp$Treatment, levels = unique(temp$Treatment[order(temp$order, temp$`Treatment class`)]))
temp$ID <- factor(temp$ID, levels = rev(sort(unique(temp$ID))))
cols <- structure(c("black", "#FFFF99","#FDBF6F","#FF7F00","#A65628","plum1","#984EA3", "lightblue","#386CB0","steelblue1","cadetblue","#E7298A","#A6D854"), 
                  names = c("LHRH", "Bicalutamide","Nilutamide","Stilbesterol","Aminoglutatemide and cortisone acetate","Abiraterone","Enzalutamine","Docetaxel","Cabazitaxel","Carboplatin and etoposide","Ketaconazole","Strontium"))

temp <- temp[order(temp$order, temp$ID),]

head(temp)
temp$`Treatment class` <- Hmisc::capitalize(temp$`Treatment class`)

p <- ggplot(temp)+
  geom_tile(aes(x = tile_start, y = ID, fill = Treatment, width = Duration, height = height),colour = "black")+  
  scale_fill_manual(values = cols) +
  geom_point( aes(x = values, y= ID, shape = treatment, color = treatment), size = 7, stroke = 1)+
  scale_shape_manual(values = c("\u006c" ,"\u2715", "\u26a1", "\u26a1"), breaks = unique(temp$treatment)[c(2, 3, 1, 4)])+
  scale_color_manual(values = c("black", "black", "gold2", "gold2"), breaks = unique(temp$treatment)[c(2, 3, 1, 4)])+
  xlab("Duration (months)") + 
  ylab("Patient ID") + 
  facet_grid(phenotype~., scales = "free", space = "free") + 
  theme_classic(base_size = 20) +
  theme(legend.position = "none")
p

grDevices::cairo_pdf("treatment_history.pdf", width = 11, height = 4.5)
p
dev.off()

temp$legends <- temp$`Treatment class`
temp$legends[temp$legends %in% c("Ketaconazole","Strontium")] <- "Other"

legends <- lapply(split(temp, temp$legends),function(df){
  g <- ggplot(df)+
    geom_tile(aes(x = tile_start, y = ID, fill = Treatment, width = Duration, height = height),colour = "black")+  
    scale_fill_manual(values = cols, name = unique(df$legends)) +
    xlab("Duration (months)") + 
    ylab("Patient ID") + 
    facet_grid(phenotype~., scales = "free", space = "free") + 
    theme_classic(base_size = 20)
  leg <- get_legend(g)
  # return(as_ggplot(leg) + theme(plot.margin=grid::unit(c(1, 1, 1, 1), "cm"),plot.background = element_rect(colour = "black", fill=NA, linewidth =1)))
  return(leg)

})

f <- ggplot(temp)+
  geom_point( aes(x = values, y= ID, shape = treatment, color = treatment), size = 7, stroke = 1)+
  scale_shape_manual(values = c("\u26a1" ,"\u2715", "\u26a1", "\u26a1"), breaks = unique(temp$treatment)[c(1, 3, 2, 4)], limits = unique(temp$treatment)[c(1, 3)], name = "Point treatment")+
  scale_color_manual(values = c( "gold2", "black", "gold2", "black"), breaks = unique(temp$treatment)[c(1, 3, 2, 4)], limits = unique(temp$treatment)[c(1, 3)], name = "Point treatment")+
  theme_classic(base_size = 20)

f2 <- ggplot(temp)+
  geom_point( aes(x = values, y= ID, shape = treatment, color = treatment), size = 7, stroke = 1)+
  scale_shape_manual(values = c("\u26a1" ,"\u2715", "\u006c" , "\u26a1"), breaks = unique(temp$treatment)[c(1, 3, 2, 4)], limits = unique(temp$treatment)[2], name = "Death")+
  scale_color_manual(values = c( "gold2", "black", "black", "black"), breaks = unique(temp$treatment)[c(1, 3, 2, 4)], limits = unique(temp$treatment)[2], name = "Death")+
  theme_classic(base_size = 20)

# legs <- legends[[3]] + legends[[1]] + legends[[2]] + legends[[4]] + legends[[5]] + f + f2 +guide_area() +plot_layout(guides = "collect")
# as_ggplot(get_legend(legs))

legends <- legends[c(3, 1, 2, 4, 5)]
legends[[6]] <- get_legend(f)
legends[[7]] <-get_legend(f2)
grDevices::cairo_pdf("legends/legends.pdf", height = 15, width = 8)
plot_grid(plotlist = legends, align = "v",axis = "tl", ncol = 1,nrow = 7, greedy = F, byrow = F, rel_heights = c(1, 3, 1, 2), rel_widths = c(2, 1))
dev.off()

legs <- plot_grid(plotlist = legends, align = "v",axis = "tl", ncol = 2,nrow = 4, greedy = F, byrow = F, rel_heights = c(1, 3, 1, 2), rel_widths = c(2, 1))
ggdraw(legs)
plot_grid(p, legs, rel_widths = c(8, 3), nrow = 1, greedy = F)

for(i in seq_along(legends)){
  temp2 <- legends[[i]]
  pdf(paste0("legends/", i, ".pdf"))
  plot(temp2)
  dev.off()
  print(i)
}
