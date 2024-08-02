# **************************************************************************
#
# CVN analysis
#
# **************************************************************************

# load data, define weight matrix and tuning parameters ---------------- 
data <- readRDS("data/KiKmedataset.rds")

# define function to create weight matrix that output a 9x9 weight matrix
W_radition <- function(a,b,g) {
  matrix(c(
    0, a, 0, g, 0, 0, 0, 0, 0, 
    a, 0, b, 0, g, 0, 0, 0, 0, 
    0, b, 0, 0, 0, g, 0, 0, 0, 
    g, 0, 0, 0, a, 0, g, 0, 0, 
    0, g, 0, a, 0, b, 0, g, 0, 
    0, 0, g, 0, b, 0, 0, 0, g, 
    0, 0, 0, g, 0, 0, 0, a, 0, 
    0, 0, 0, 0, g, 0, a, 0, b, 
    0, 0, 0, 0, 0, g, 0, b, 0
  ), 
  ncol = 9)}
a <- 975 / 1000



# Apply the CVN ---------------------------------------------------------------------
  cvn <- CVN::CVN(data$X, 
                  W = W_radition(a = a, b = 0.025, g = 0.5), 
                  gamma1 = 5e-5, 
                  gamma2 = 5e-6, 
                  eps = 1e-3, maxiter = 10000, 
                  n_cores = 11, verbose = FALSE)
save(cvn, file = "data/cvn.rda")


# plot graph ------------------------------------------------------------------------

# 1) plot cvn object using the plot function of the CVN package
pl_cvn <- plot(cvn)

# Convert visNetwork plots to gridExtra grobs
for (g in seq_along(pl_cvn$plots[[1]])){
  plot_file_name <- paste0('results/kikme-graph_', g)
  
  # save plot as html file
  htmlwidgets::saveWidget(pl_cvn$plots[[1]][g][[1]], 
                          paste0(plot_file_name, ".html"), 
                          selfcontained = TRUE)
  
  # capture html as png image
  img <- webshot::webshot(paste0(plot_file_name, ".html"), 
                          file = paste0(plot_file_name, ".png"),
                          cliprect = c(50, 260, 475, 475)
                          )
  
}

# 2) plot hamming distance
# look at it
hd <- CVN::hamming_distance(cvn)

# plot it
pl_hd <- CVN::plot_hamming_distances_cvn(cvn)

# Prepare data for plotting HD with different color scheme (ggplot2)
distance_matrix <- matrix(pl_hd$distances[[1]], 9, 9)
hd_data <- reshape2::melt(distance_matrix)
hd_data$label1 <- labels[hd_data[,1]]
hd_data$label2 <- labels[hd_data[,2]]

# Plot Hamming distance
ggplot(data = hd_data, aes(x = label1, y = label2, fill = value)) + 
  geom_tile() + 
  #ggtitle("Weights") +  
  xlab("") + 
  ylab("") + 
  theme(axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_text()
        ) +
  scale_y_discrete(limits = rev(unique(hd_data$label2))) +
  scale_x_discrete(limits = unique(hd_data$label1)) +
  viridis::scale_fill_viridis(
    option = "cividis", 
    direction = -1, 
    name = 'Hamming\ndistance') + # Use viridis color scale
  geom_text(aes(label = value), color = 'white', size = 3)  # Add values to cells

ggsave(filename = "results/heatmap.pdf", width = 16, height = 10, unit = "cm")




# Network characteristics -----------------------------------------------------------
source('analysis/network_datamanagement.R')

# convert cvn result to igraph object
icvn <- cvn2igraph(cvn)[[1]]
# name nodes
for(i in 1:9){
  V(icvn[[i]])$name <- genes$external_gene_name
}
# add network descriptives to igraph object
icvn <- lapply(icvn, make_my_graph)
# create graph representations and calculate graph descriptives 
gd <- lapply(icvn, graph_descriptives)
gd <- do.call(cbind, gd)
names(gd) <- c(paste0(names(gd), rep(1:9, each = 2)))
# update graph descriptives with gene names
 for (i in 1:9){
   gd[,i*2] <- V(icvn2[[i]])$name[gd[,i*2]]
 }
# Convert graph descriptives to a data frame
gd <- data.frame(characteristics = rownames(gd), gd)
rownames(gd) <- NULL
 
# Transpose data to create a table "Characteristics of the 9 estimated graphs" 
# in the text
 tab_data <- gd[c(7,8,1,5,10,9),c(1,2,4,6,8,10,12,14,16,18)]
 rownames(tab_data) <- tab_data[,1]
 tab_data <- t(tab_data[,2:10])
 rownames(tab_data) <- c(paste('N0', c('0 Gy', '0.05 Gy', '2 Gy'), sep = ', '), 
                         paste('N1', c('0 Gy', '0.05 Gy', '2 Gy'), sep = ', '), 
                         paste('N2', c('0 Gy', '0.05 Gy', '2 Gy'), sep = ', '))
 
# save as latex table
latex_table <- xtable::xtable(tab_data)
xtable::print.xtable(latex_table, 
                     file = 'results/table_networkcharacteristics.tex', 
                     booktabs = TRUE, 
                     include.rownames = TRUE)

 

# tp53 ------------------------------------------------------------------------------
# just keep edges of node TP53 and rearange data set
tp53 <- dosis %>% 
         filter(gene1 ==  "TP53" | gene2 ==  "TP53") %>%
         dplyr::select(gene1, gene2, des1, des2, chr1, chr2, core, N00:N22) %>%
         mutate(Gene = if_else(gene1 == "TP53", gene2, gene1),
                Description = if_else(des1 == "Tumor Protein P53", des2, des1),
                Chr = if_else(gene1 == "TP53", chr2, chr1),
                total = rowSums(across(starts_with("N")))) %>%
         relocate(gene1, gene2, chr1, chr2, core, Gene, Description, Chr, total, 
                  .before = N00) %>%
         arrange(desc(core), desc(N00))
# correct chromosome of gene SESN2
tp53$Chr[tp53$Gene == "SESN2"] <- 1

# columns used in manuscript
tp53_out <- tp53[,c(8,10,9,11,12:20)] %>%
  arrange(N00, N01, N02, N10, N11, N12, N20, N21, N22, Chr) %>%
  mutate(comment = case_when(total == 9 ~ 'in all graphs',
                             N00 == 0 & N01 == 0 & N02 == 0 & N10 == 0 & N11 == 0 & N12 == 0 ~ 'N2',
                             N20 == 0 & N21 == 0 & N22 == 0 & N10 == 0 & N11 == 0 & N12 == 0 ~ 'control group',
                             N00 == 0 & N01 == 0 & N10 == 0 & N11 == 0 & N20 == 0 & N21 == 0 ~ '2 Gy',
                             Gene == 'TAF3' ~ 'in N0 and N1 with <= 0.05 Gy',
                             Gene %in% c('XRCC6', 'FBXO22') ~ 'in cancer group(s) and for N0 with 2 Gy',
                             N00 == 1 & N01 == 1 & N02 == 1 & total %in% c(5:7) ~ '0-2 Gy in control group, unclear pattern for cancer groups')) %>%
  dplyr::select(Gene, Chr, Description, comment, N00, N01, N02, N10, N11, N12, N20, N21, N22)


# make latex table
latex_table_tp53 <- xtable::xtable(tp53_out)
xtable::digits(latex_table_tp53)[6:14] <- 0

# Save results
xtable::print.xtable(latex_table_tp53, 
                     file = 'results/table_tp53.tex', 
                     booktabs = TRUE, include.rownames = FALSE)
