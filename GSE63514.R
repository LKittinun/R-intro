library(tibble)
library(dplyr)
library(purrr)
library(limma)
library(affy)
library(AnnotationDbi)
library(hgu133plus2.db)
library(GEOquery)

library(rvest)
gse<- getGEO("GSE63514", destdir = "D:/")
gse <- getGEO(filename="D:/GSE63514_series_matrix.txt.gz")
pData(gse) |> write.csv("GSE63514_meta.csv")



GSE63514_files <- list.files(path = "Resource/GSE63514_RAW", full.names = TRUE)
GSE63514_raw <- ReadAffy(filenames = GSE63514_files)
GSE63514 <- rma(GSE63514_raw)
?select
GSE63514_genes <- AnnotationDbi::select(hgu133plus2.db, keys = keys(hgu133plus2.db), columns=c("SYMBOL", "ENTREZID", "GENENAME")) |> 
  summarize(across(c(SYMBOL, ENTREZID, GENENAME), ~paste(.x, collapse = ", ")), .by = PROBEID)

write.csv(GSE63514_genes, "Resource/hgu133plus2_genenames.csv")
write.csv(exprs(GSE63514), "GSE63514_norm.csv")

pGroup <- pData(GSE63514) |>
  mutate(group = gsub(".*(Normal|CIN1|CIN2|CIN3|Cancer).*", "\\1" , rownames(pData(GSE63514)) )) |> 
  pull(group)

library(limma)

design_sep_CIN <- model.matrix(~pGroup+0)
design_sep_CIN_contrast <- apply(combn(colnames(design_sep_CIN),2)[2:1,], 2, paste, collapse = "-")
design_sep_CIN_contrast <- makeContrasts(contrasts = design_sep_CIN_contrast, levels = colnames(design_sep_CIN))

sep_CIN_fit <- lmFit(GSE63514, design_sep_CIN)
sep_CIN_fit <- contrasts.fit(sep_CIN_fit, design_sep_CIN_contrast)
sep_CIN_fit <- eBayes(sep_CIN_fit)

sep_CIN_sig <- map(colnames(sep_CIN_fit$coef), 
    safely(\(x) topTable(sep_CIN_fit, coef = x, number = Inf, p = 0.05) %>%
             rownames_to_column("PROBEID") %>%
             inner_join(GSE63514_genes) %>%
             relocate(c("SYMBOL", "ENTREZID", "GENENAME"), .before = "logFC")) ) |> 
  set_names(colnames(sep_CIN_fit$coef)) |> 
  janitor::clean_names()

####


