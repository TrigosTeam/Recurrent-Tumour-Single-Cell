normal3 <- readRDS("~/normal_cell_annotation/subtype/final_normalsrt.Rds")
normal3$cell_type <- ifelse(normal3$cell_anno%in% c("T cells","B cells","Plasma cells",  "Macrophages"), "immune cells", 
                            ifelse(normal3$cell_anno%in% c("Fibroblasts","Adipocytes","Pericytes"), "stromal cells", "others"))
normal3$pathology <- ifelse(normal3$patient %in% c("CA0090", "CA0046"), "NE", ifelse(normal3$patient %in% c("CA0027", "CA0058"), "Mixed", "AD"))

df <- normal3@meta.data %>% group_by(patient,site,sample, cell_anno) %>% summarise(n = n()) %>% mutate (freq = n/sum(n))

df$pathology <- ifelse(df$patient %in% c("CA0090", "CA0046"), "NE", ifelse(df$patient %in% c("CA0027", "CA0058"), "Mixed", "AD"))
df$patient <- gsub("00", "", df$patient)
df$labs <- paste0(signif(df$freq*100, 2), "%")
df$labs[df$freq < 0.03] <- NA

sites <- df$sample
sites[grep("brain", sites)] <- "brain"
sites[grep("dura", sites)] <- "dura"
sites[grep("prostate", sites)] <- "prostate"
sites[grep("liver", sites)] <- "liver"
sites[grep("fat", sites)] <- "fat"
sites[grep("lymph|hilar", sites)] <- "LN"
sites[grep("rib|vertebra", sites)] <- "bone"
sites[grep("lung", sites)] <- "lung"
sites[grep("abdomen", sites)] <- "abdomen"
sites[grep("bladder", sites)] <- "bladder"
df$site <- sites



df2 <- meta %>% group_by(sample) %>% summarise(avg1 = mean(Module1),
                                                 avg2 = mean(Module2),
                                                 avg3 = mean(Module3),
                                                 avg4 = mean(Module4),
                                                 avg5 = mean(Module5),
                                                 avg6 = mean(Module6)) %>% 
  mutate(tile1 = ordered(ntile(avg1, 5)),
         tile2 =  ordered(ntile(avg2, 5)),
         tile3 =  ordered(ntile(avg3, 5)),
         tile4 =  ordered(ntile(avg4, 5)),
         tile5 =  ordered(ntile(avg5, 5)),
         tile6 =  ordered(ntile(avg6, 5)))

temp1 <- reshape2::melt(df2[,1:7], id.vars = "sample", variable.name = "avg", value.name = "exp")
temp1 <- temp1 %>% group_by(avg) %>% mutate(tile = ordered(ntile(exp, 5)))
temp1$avg <- as.character(temp1$avg)
for (i in 1:6){
  temp1$avg[temp1$avg == paste0("avg", i)] <- new_name[i]
}
colnames(temp1) <- c("sample", "module", "exp", "tile")
temp1 <- merge(df, temp1, by = "sample", all.x = T, all.y = F)

temp  <- merge(df, df2, by = "sample", all.x = T, all.y = F)

plotdf <- data.frame()
pdf("pdf/normalcell_propVSbinedmodule_n=5.pdf", width = 10, height = 5)
for (n in levels(normal3$cell_anno)){
  subdf <-temp[temp$cell_anno == n,]
  for (m in 1:6){
    test <- DescTools::JonckheereTerpstraTest(subdf$freq, g = subdf[, paste0("tile", m )], alternative = "increasing", nperm = 5000)
    test2 <- DescTools::JonckheereTerpstraTest(subdf$freq, g = subdf[, paste0("tile", m )], alternative = "decreasing", nperm = 5000)
  f <- ggplot(subdf, aes(x = subdf[, paste0("tile", m )], y = freq)) + geom_violin()+
    geom_boxplot(outlier.shape = NA, width = 0.1) + geom_jitter(aes(colour = patient)) + 
    scale_color_manual(values  = patient_cols)+
    labs(x = paste("bins of", new_name[m] , "Module mean expression"), y = paste(n, "proportion"))+theme_minimal(base_size = 16)
  
  g <- ggplot(subdf, aes(x = subdf[, paste0("avg", m)], y = freq ,colour = patient)) + geom_point() +  
    geom_smooth(method=lm, color="black", se=FALSE)+ 
    scale_color_manual(values =patient_cols)+
    labs(x = paste( new_name[m], "intensity"), y = paste(n, "proportion"))+theme_minimal(base_size = 16)
  cor = cor.test(subdf$freq[!is.na(subdf[, paste0("avg", m)])], subdf[, paste0("avg", m)][!is.na(subdf[, paste0("avg", m)])], method = "spearman")
  print(f+g+plot_annotation(title = paste("correlation between", new_name[m], "Module expression and", n, "proportion:\nrho = ", round(cor$estimate,3), "p-value", round(cor$p.value,3)),
                            subtitle =  paste0("JonckheereTerpstraTest increasing p-value:", test$p.value, " decreasing p-value:",test2$p.value))+
          plot_layout(axes = "collect", guides = "collect"))
  if (any(c(test$p.value, test2$p.value, cor$p.value)<0.05)){
    plotdf <- rbind(plotdf, data.frame(subdf[, 1:7], module = new_name[m], 
                                       exp = subdf[, paste0("avg", m)], tile = subdf[, paste0("tile", m )], 
                                       JTincP = test$p.value, 
                                       JTdecP = test2$p.value,
                                       Cor = cor$estimate,
                                       CorP = cor$p.value))
  }
  
  }
  print (n)
}
dev.off()

saveRDS(plotdf, "sigModulenormal.Rds")

# df2 <- meta %>% group_by(sample) %>% summarise(avg1 = mean(Module1),
#                                                avg2 = mean(Module2),
#                                                avg3 = mean(Module3),
#                                                avg4 = mean(Module4),
#                                                avg5 = mean(Module5),
#                                                avg6 = mean(Module6)) %>% 
#   mutate(tile1 = ordered(ntile(avg1, 3)),
#          tile2 =  ordered(ntile(avg2, 3)),
#          tile3 =  ordered(ntile(avg3, 3)),
#          tile4 =  ordered(ntile(avg4, 3)),
#          tile5 =  ordered(ntile(avg4, 3)),
#          tile6 =  ordered(ntile(avg6, 3)))
# 
# temp  <- merge(df, df2, by = "sample", all.x = T, all.y = F)
# 
# plotdf <- data.frame()
# pdf("pdf/normalcell_propVSbinedmodule_n=3.pdf", width = 10, height = 5)
# for (n in levels(normal3$cell_anno)){
#   subdf <-temp[temp$cell_anno == n,]
#   for (m in 1:6){
#     test <- DescTools::JonckheereTerpstraTest(subdf$freq, g = subdf[, paste0("tile", m )], alternative = "increasing", nperm = 5000)
#     test2 <- DescTools::JonckheereTerpstraTest(subdf$freq, g = subdf[, paste0("tile", m )], alternative = "decreasing", nperm = 5000)
#     f <- ggplot(subdf, aes(x = subdf[, paste0("tile", m )], y = freq)) + geom_violin()+
#       geom_boxplot(outlier.shape = NA, width = 0.1) + geom_jitter(aes(colour = patient)) + 
#       scale_color_manual(values  = patient_cols)+
#       labs(x = paste("bins of", new_name[m] , "Module mean expression"), y = paste(n, "proportion"))+theme_minimal(base_size = 16)
#     
#     g <- ggplot(subdf, aes(x = subdf[, paste0("avg", m)], y = freq ,colour = patient)) + geom_point() +  
#       geom_smooth(method=lm, color="black", se=FALSE)+ 
#       scale_color_manual(values =patient_cols)+
#       labs(x = paste( new_name[m], "intensity"), y = paste(n, "proportion"))+theme_minimal(base_size = 16)
#     cor = cor.test(subdf$freq[!is.na(subdf[, paste0("avg", m)])], subdf[, paste0("avg", m)][!is.na(subdf[, paste0("avg", m)])], method = "spearman")
#     print(f+g+plot_annotation(title = paste("correlation between", new_name[m], "Module expression and", n, "proportion:\nrho = ", round(cor$estimate,3), "p-value", round(cor$p.value,3)),
#                               subtitle =  paste0("JonckheereTerpstraTest increasing p-value:", test$p.value, " decreasing p-value:",test2$p.value))+
#             plot_layout(axes = "collect", guides = "collect"))
#     if (any(c(test$p.value, test2$p.value, cor$p.value)<0.05)){
#       plotdf <- rbind(plotdf, data.frame(subdf[, 1:7], module = new_name[m], 
#                                          exp = subdf[, paste0("avg", m)], tile = subdf[, paste0("tile", m )], 
#                                          JTincP = test$p.value, 
#                                          JTdecP = test2$p.value,
#                                          Cor = cor$estimate,
#                                          CorP = cor$p.value))
#     }
#     
#   }
#   print (n)
# }
# dev.off()
# saveRDS(plotdf, "sigModulenormaln=3.Rds")

sigModulenormal%>% select(c("cell_anno", "module")) %>% distinct()
head(sigModulenormal)

sigModulenormal <- readRDS("~/Fig4-5_archetype/sigModulenormal.Rds")
sigModulenormal$module[sigModulenormal$module == "Hypoxia"] <- "Glycolysis"
sigModulenormal$module <- factor(sigModulenormal$module, levels = new_name)
# ggplot(temp1 %>% filter(cell_anno %in% sigModulenormal$cell_anno), aes(x = tile, y = freq,fill = module)) + geom_violin()+
#   geom_boxplot(outlier.shape = NA, width = 0.1) + geom_jitter(aes(colour = patient)) + 
#   scale_color_manual(values  = patient_cols)+
#   scale_fill_manual(values = module_cols2)+
#   theme_minimal(base_size = 16)+
#   facet_wrap(.~cell_anno)
#   
ggplot(sigModulenormal, aes(x = tile, y = freq,fill = module)) + geom_violin()+
  geom_boxplot(outlier.shape = NA, width = 0.1) + geom_jitter(aes(colour = patient)) + 
  scale_color_manual(values  = patient_cols)+
  scale_fill_manual(values = module_cols2)+
  theme_bw()+
  ggh4x::facet_nested_wrap(.~cell_anno+module, scale = "free", nrow = 1)+
  theme(strip.background = element_rect(fill = "white"), 
        axis.title.x = element_text(size = 18),
        strip.text = element_text(size = 16, face = "bold"), 
        axis.title.y = element_text(size = 18),
        legend.position = "bottom", 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16))+
  guides(fill = guide_legend(nrow = 1),
         color = guide_legend(nrow = 1))+
  labs(y = "Proportion", x = "Expression Levels of Modules", fill = "Module", color = "Patient")
ggsave("pdf/JTtest.pdf", width = 18, height = 3.5)
