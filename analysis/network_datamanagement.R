library(biomaRt)
library(gtools)

# load data and functions
source("functions.R")
load("data/cvn.rda")
data <- readRDS("data/KiKmedataset.rds")


# Data management -------------------------------------------------------------------
# Core graph contains only edges that are shared by all graphs
core <- CVN::find_core_graph(cvn)[[1]]

# Extract edges that are shared by all graphs
el_core <- make_edge_list(core)

# Extract edges that are shared between graphs
el_between <- reduce_my_adj(cvn)
# Convert CVN results to igraph objects and extract edge lists (length = 36)
tmp <- cvn2igraph(cvn)[[1]]
tmp  <- lapply(tmp, as_edgelist)

# list of edges for each of the 9 graphs
el_graphs <- list()
for(i in 1:9){
  el_graphs[[i]] <- paste(tmp[[i]][,1], tmp[[i]][,2], sep = "-")
}


# Make dataset to track edge inclusion across graphs
posedges <- gtools::combinations(cvn$p, 2, 
                                 v = 1:cvn$p, 
                                 set = TRUE, 
                                 repeats.allowed = FALSE)

el_data <- data.frame(posedges = paste(posedges[,1], posedges[,2], sep = "-")) 
# add info whether edge is in core graph
el_data <- el_data %>%
  mutate(core = if_else(posedges %in% el_core, 1, 0))


# Update el_data with shared edges information
for (i in seq_along(el_between)){                    
  varname <- names(el_between)[i]
  el_data <- el_data %>%
    mutate(!!varname := if_else(posedges %in% el_between[[i]], 1, 0))
}

# Update el_data with edge inclusion information across different graphs
# el_data: 18145 x 57
el_data <- el_data %>%
  mutate(N00 = if_else(posedges %in% el_graphs[[1]], 1, 0),
         N01 = if_else(posedges %in% el_graphs[[2]], 1, 0),
         N02 = if_else(posedges %in% el_graphs[[3]], 1, 0),
         N10 = if_else(posedges %in% el_graphs[[4]], 1, 0),
         N11 = if_else(posedges %in% el_graphs[[5]], 1, 0),
         N12 = if_else(posedges %in% el_graphs[[6]], 1, 0),
         N20 = if_else(posedges %in% el_graphs[[7]], 1, 0),
         N21 = if_else(posedges %in% el_graphs[[8]], 1, 0),
         N22 = if_else(posedges %in% el_graphs[[9]], 1, 0),
         D00 = if_else(posedges %in% el_graphs[[1]], 1, 0),
         D01 = if_else(posedges %in% el_graphs[[4]], 1, 0),
         D02 = if_else(posedges %in% el_graphs[[7]], 1, 0),
         D10 = if_else(posedges %in% el_graphs[[2]], 1, 0),
         D11 = if_else(posedges %in% el_graphs[[5]], 1, 0),
         D12 = if_else(posedges %in% el_graphs[[8]], 1, 0),
         D20 = if_else(posedges %in% el_graphs[[3]], 1, 0),
         D21 = if_else(posedges %in% el_graphs[[6]], 1, 0),
         D22 = if_else(posedges %in% el_graphs[[9]], 1, 0),
         none = if_else(posedges %in% purrr::simplify(el_graphs), 0, 1))


# add gene names from ensemble
genes <- gtools::combinations(length(data$gene_labels), 2, 
                              v = data$gene_labels, 
                              set = TRUE, repeats.allowed = FALSE)

genes <- data.frame(genes, edge = paste(genes[,1], genes[,2], sep = "-")) 
names(genes) <- c("gene1", "gene2", "genes")
el_data <- data.frame(el_data, genes)


# Edge differences between dosage 0 Gy and 2 Gy in all 3 groups
# Calculate edge differences between dosage 0 Gy and 2 Gy in all 3 groups
dosis <- el_data %>% 
  filter(N00 == 1 | N01 == 1 | N02 == 1 | N10 == 1 | N11 == 1 | N12 == 1 | N20 == 1 | 
         N21 == 1 | N22 == 1) %>%
  mutate(
    N0_02 = if_else(N00 != N02, 1, 0),
    N1_02 = if_else(N10 != N12, 1, 0),
    N2_02 = if_else(N20 != N22, 1, 0),
    N0_p02 = if_else(N00 == 0 & N02 == 1, 1, 0),
    N0_m02 = if_else(N00 == 1 & N02 == 0, 1, 0),
    N1_p02 = if_else(N10 == 0 & N12 == 1, 1, 0),
    N1_m02 = if_else(N10 == 1 & N12 == 0, 1, 0),
    N2_p02 = if_else(N20 == 0 & N22 == 1, 1, 0),
    N2_m02 = if_else(N20 == 1 & N22 == 0, 1, 0),
    D0_p02 = if_else(N00 == 0 & N20 == 1, 1, 0),
    D0_m02 = if_else(N00 == 1 & N20 == 0, 1, 0),
    D1_p02 = if_else(N01 == 0 & N21 == 1, 1, 0),
    D1_m02 = if_else(N01 == 1 & N21 == 0, 1, 0),
    D2_p02 = if_else(N02 == 0 & N22 == 1, 1, 0),
    D2_m02 = if_else(N02 == 1 & N22 == 0, 1, 0),
    N0 = rowSums(across(c(N00, N01, N02), ~ . == 1)) == 3,
    N1 = rowSums(across(c(N10, N11, N12), ~ . == 1)) == 3,
    N2 = rowSums(across(c(N20, N21, N22), ~ . == 1)) == 3) 


# Connect to the Ensembl database to check for TP53-related edges
library(biomaRt)
ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Define unique genes from dosage data
gene_dosis <- unique(c(dosis$gene1, dosis$gene2))
gene_id <- data$gene_labels

# Get gene information from Ensembl
gene_dosis_info <- getBM(
  attributes = c('ensembl_gene_id', 'entrezgene_id', 'external_gene_name', 'chromosome_name', 'description'),
  filters = "ensembl_gene_id",
  values = gene_id,
  mart = ensembl
)
genes <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name'),
  filters = "ensembl_gene_id",
  values = gene_id,
  mart = ensembl
)
# data correction
genes[genes[,2]=="",] <- c("novel transcript chr9", "novel transcript chr20")


# Merge gene information with dosage data
dosis2 <- merge(dosis, 
                subset(gene_dosis_info, select = -c(entrezgene_id)), 
                by.x = "gene1",
                by.y = "ensembl_gene_id",
                all.x = TRUE)
dosis2 <- merge(dosis2, 
                subset(gene_dosis_info, select = -c(entrezgene_id)), 
                by.x = "gene2",
                by.y = "ensembl_gene_id",
                all.x = TRUE)

# Clean and reorder columns
dosis <- dosis2 %>% 
  dplyr::select(-c(gene1, gene2, posedges)) %>%
  dplyr::rename(
    gene1 = external_gene_name.x,
    chr1 = chromosome_name.x,
    gene2 = external_gene_name.y,
    chr2 = chromosome_name.y,
    des1 = description.x,
    des2 = description.y
  ) %>%
  mutate(des1 = stringr::str_to_title(
    stringr::str_sub(des1, 1, stringr::str_locate(des1, "\\[")[, 1] - 2)),
    des2 = stringr::str_to_title(
      stringr::str_sub(des2, 1, stringr::str_locate(des2, "\\[")[, 1] - 2))) %>%
  relocate(gene1, gene2, chr1, chr2, .before = N00) %>% 
  arrange(desc(N0_02), desc(N1_02), desc(N2_02)) 



# Remove temporary variables
rm(dosis2, gene_dosis_info, gene_dosis, gene_id, tmp,
   ensembl, posedges, el_between, el_core, el_graphs)
