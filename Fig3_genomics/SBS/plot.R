library(patchwork)
setwd("~/CASCADEpaper/paper/Fig3_genomics/SBS")
chord_df$sample <- paste(chord_df$id, chord_df$site, sep = "_") 
chord_df%>% select(c("sample", "p_hrd"))
df <- merge(plt_df, chord_df%>% select(c("sample", "p_hrd")), by = "sample")
lesion <- df$site
lesion[grep("brain", lesion)] <- "brain"
lesion[grep("dura", lesion)] <- "dura"
lesion[grep("prostate", lesion)] <- "prostate"
lesion[grep("liver", lesion)] <- "liver"
lesion[grep("fat", lesion)] <- "fat"
lesion[grep("lymph|hilar", lesion)] <- "LN"
lesion[grep("rib|vertebra", lesion)] <- "bone"
lesion[grep("lung", lesion)] <- "lung"
lesion[grep("abdomen", lesion)] <- "abdomen"
lesion[grep("bladder", lesion)] <- "bladder"
length(lesion)
lesion <- paste(lesion, sapply(strsplit(df$sample, split = "_"), function(x) tail(x, 1)))
lesion[df$sample == "CA0035_left_paraaortic_lymph_node_track_1A"] <- "LN 1"
lesion <- lesion[df$sample != "CA0035_left_paraaortic_lymph_node_track_1B"]
df$id <- gsub("00", "", df$id)
saveRDS(df, "plotdf.Rds")
df <- df[df$sample != "CA0035_left_paraaortic_lymph_node_track_1B",]

f1 <- ggplot(df, aes(x = sample, y = p_hrd, fill = id)) + geom_col() + ylim(0,1) + geom_hline(yintercept = 0.5) + 
  theme_bw(base_size = 18)+ 
  ggh4x::facet_nested(.~pathology+id, scales = "free", space = "free")+
  theme(axis.text.x = element_blank(), legend.position = "none", 
        axis.ticks.x = element_blank(),strip.background = element_rect(fill ="white"), axis.title.x = element_blank(),
        panel.spacing = unit(1, "mm") ) 
f1
g <- ggplot(df, aes(x = sample, y = count, fill = sig)) + geom_col() +  theme_bw(base_size = 18)+ 
  facet_grid(.~pathology+id, scales = "free", space = "free") + 
  scale_fill_manual(values = c(paletteer_d("ggsci::springfield_simpsons") %>% as.character(), "lightgrey")) + 
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90),   strip.background = element_blank(),strip.text.x = element_blank(),
        panel.spacing = unit(1, "mm") )+
  scale_x_discrete(breaks = df$sample, labels=lesion )
g

p <- f1/g + plot_layout(axes = "collect", guides = "collect", heights = c(1, 3))
p
ggsave("sbs.pdf", width = 7.5, height = 8)
png(paste0("~/CASCADEpaper/paper/Fig3/SBS/png/sbs.png"), height = 700, width = 500)
p
dev.off()
saveRDS(p, "~/CASCADEpaper/paper/Fig3/SBS/plot_obj/sbs.Rds")