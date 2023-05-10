if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Mm.eg.db")



library(tidyverse)
library("readxl")
res <- read_excel("C:/data/CORBIN/Graber/new/B.vs.A_adult and older.xlsx")
res

res_e <- read_excel("C:/data/CORBIN/Graber/new/B.vs.control_E vs. A.xlsx")

res2 <- res %>% 
  dplyr::select(gene_id, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene_id) %>% 
  summarize(stat=mean(stat))
res2

library(fgsea)
ranks <- deframe(res2)
head(ranks, 20)
pathways.hallmark <- gmtPathways("C:/data/CORBIN/Graber/new713/MousePath_All_gmt-Format.gmt")
pathways.hallmark %>% 
  head() %>% 
  lapply(head)

fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

older <- pathways.hallmark %>% 
  enframe("pathway", "gene_id") %>% 
  unnest(cols = c(gene_id)) %>% 
  
  inner_join(res, by="gene_id")

elderly <- pathways.hallmark %>% 
  enframe("pathway", "gene_id") %>% 
  unnest(cols = c(gene_id)) %>% 
  
  inner_join(res_e, by="gene_id")

older1 <- older[!duplicated(older$gene_id), ]
elderly1 <- elderly[!duplicated(older$gene_id), ]
older2 <- older %>% distinct(gene_id, .keep_all = TRUE)
install.packages('xlsx')
library(xlsx)
write.xlsx(older1, "C:/data/CORBIN/Graber/new713/GSEA_adult_older_result.xlsx")
write.xlsx(elderly1, "C:/data/CORBIN/Graber/new713/GSEA_adult_elderly_result.xlsx")

