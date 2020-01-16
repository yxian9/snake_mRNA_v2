
library(ballgown)
library(magrittr)
args <- commandArgs(trailingOnly = TRUE)

sampleids = args[1:length(args)-1]
output = args[length(args)]
print(sampleids)
## generate the  expression table for each sample
bg <-  ballgown(samples = sampleids , meas='FPKM')

print('process sample data')

print(sampleNames(bg))
gene_expression = ballgown::gexpr(bg)
gene_table = tibble::as.tibble(gene_expression)
gene_table$geneid = rownames(gene_expression)

geneID = tibble::tibble(genename = ballgown::geneNames(bg) ,
                        geneid =ballgown::geneIDs(bg) ) %>%
  dplyr::filter(! genename =='.' ) %>% dplyr::distinct()

gene_table <- dplyr::inner_join(gene_table ,geneID )
readr::write_csv(gene_table, path = output)
