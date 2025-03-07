library(ggplot2)
library(pals)
df <- normal_clean5@meta.data %>% group_by(sample, cell_id) %>% summarise(n = n()) %>% mutate (freq = n/sum(n))
df$patient <- substr(df$sample, 1, 6)
df$cell_type <- ifelse(df$cell_id%in% c("T cell","B cell","Plasma cell",  "Naive cell", "Macrophage", "Erythroblast"), "immune cells", "stromal cells")
df$cell_id <- factor(df$cell_id, levels = c("T cell","B cell","Plasma cell",  "Naive cell", "Macrophage", "Erythroblast", "Epithelial cell", "Adipocyte","Chondrocyte", "Smooth muscle", "Fibroblast", "Endothelial cell", "Hepatocyte", "Neuron"))


ggplot(df, aes(x = sample, y = freq, fill = cell_id)) + geom_bar(stat = "identity") + 
  facet_grid(.~patient, scales = "free", space = "free") + 
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) + 
  scale_fill_manual(values = brewer.set1(14))

df <- normal_clean5@meta.data %>% group_by(site, cell_id) %>% summarise(n = n()) %>% mutate (freq = n/sum(n))
df$cell_type <- ifelse(df$cell_id%in% c("T cell","B cell","Plasma cell",  "Naive cell", "Macrophage", "Erythroblast"), "immune cells", "stromal cells")
df$cell_id <- factor(df$cell_id, levels = c("T cell","B cell","Plasma cell",  "Naive cell", "Macrophage", "Erythroblast", "Epithelial cell", "Adipocyte","Chondrocyte", "Smooth muscle", "Fibroblast", "Endothelial cell", "Hepatocyte", "Neuron"))

ggplot(df, aes(x = site, y = freq, fill = cell_id)) + geom_bar(stat = "identity") + 
  facet_grid(cell_type~., scales = "free", space = "free") + 
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) + 
  scale_fill_manual(values = brewer.set1(14))


normal_clean5$cell_type <- ifelse(normal_clean5$cell_id%in% c("T cell","B cell","Plasma cell",  "Naive cell", "Macrophage", "Erythroblast"), "immune cells", "stromal cells")
df <- normal_clean5@meta.data %>% group_by(site, cell_type) %>% summarise(n = n()) %>% mutate (freq = n/sum(n))

ggplot(df, aes(x = site, y = freq, fill = cell_type)) + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) + 
  scale_fill_manual(values = brewer.set2(5))


df <- normal_clean5@meta.data %>% group_by(site, cell_id) %>% summarise(n = n()) %>% mutate (freq = n/sum(n))
df$cell_type <- ifelse(df$cell_id%in% c("T cell","B cell","Plasma cell",  "Naive cell", "Macrophage", "Erythroblast"), "immune cells", "stromal cells")
df$cell_id <- factor(df$cell_id, levels = c("T cell","B cell","Plasma cell",  "Naive cell", "Macrophage", "Erythroblast", "Epithelial cell", "Adipocyte","Chondrocyte", "Smooth muscle", "Fibroblast", "Endothelial cell", "Hepatocyte", "Neuron"))

ggplot(df, aes(x = site, y = freq, fill = cell_id)) + geom_bar(stat = "identity") + 
  facet_grid(cell_type~., scales = "free", space = "free") + 
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) + 
  scale_fill_manual(values = brewer.set1(14))

df <- normal_clean5@meta.data %>% group_by(sample, cell_id) %>% summarise(n = n()) %>% mutate (freq = n/sum(n))
df$patient <- substr(df$sample, 1, 6)
df$pathology <- ifelse(df$patient %in% c("CA0090", "CA0046"), "NE", ifelse(df$patient %in% c("CA0027", "CA0058"), "Mixing", "AD"))
df$cell_type <- ifelse(df$cell_id%in% c("T cell","B cell","Plasma cell",  "Naive cell", "Macrophage", "Erythroblast"), "immune cells", "stromal cells")
df$cell_id <- factor(df$cell_id, levels = c("T cell","B cell","Plasma cell",  "Naive cell", "Macrophage", "Erythroblast", "Epithelial cell", "Adipocyte","Chondrocyte", "Smooth muscle", "Fibroblast", "Endothelial cell", "Hepatocyte", "Neuron"))

ggplot(df, aes(x = sample, y = freq, fill = cell_id)) + geom_bar(stat = "identity") + 
  facet_grid(.~pathology+patient, scales = "free", space = "free") + 
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) + 
  scale_fill_manual(values = brewer.set1(14))


normal_clean5$pathology <- ifelse(normal_clean5$patient %in% c("CA0090", "CA0046"), "NE", ifelse(normal_clean5$patient %in% c("CA0027", "CA0058"), "Mixing", "AD"))
df <- normal_clean5@meta.data %>% group_by(pathology, cell_id) %>% summarise(n = n()) %>% mutate (freq = n/sum(n))
df$cell_id <- factor(df$cell_id, levels = c("T cell","B cell","Plasma cell",  "Naive cell", "Macrophage", "Erythroblast", "Epithelial cell", "Adipocyte","Chondrocyte", "Smooth muscle", "Fibroblast", "Endothelial cell", "Hepatocyte", "Neuron"))
ggplot(df, aes(x = pathology, y = freq, fill = cell_id)) + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) + 
  scale_fill_manual(values = brewer.set1(14))

df <- normal_clean5@meta.data %>% group_by(pathology, cell_type) %>% summarise(n = n()) %>% mutate (freq = n/sum(n))
ggplot(df, aes(x = pathology, y = freq, fill = cell_type)) + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) + 
  scale_fill_manual(values = brewer.set2(5))


df <- normal_clean5@meta.data %>% group_by(patient,site, cell_id) %>% summarise(n = n()) %>% mutate (freq = n/sum(n))
ggplot(df, aes(x = site, y = freq, fill = cell_id)) + geom_bar(stat = "identity") + 
  facet_grid(.~ patient, scales = "free", space = "free") + 
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90)) + 
  scale_fill_manual(values =  brewer.set1(14))
