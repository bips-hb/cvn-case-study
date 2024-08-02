# **************************************************************************
#
# Preprocess the KiKme dataset
#
# **************************************************************************
# load the raw data provided by KiKME project (not allowed to share)
# 
stop('We are not allowed to share the data because of privacy reasons.')
load(file = "data/genetic_data_response_to_ionizing_radiation.RData")

# Patient 122 has only 2 (not 3) experiments and will be removed
meta_sens <- meta[-which(meta$Patient == 122),]

# meta contains the data on the individual patients, 
# which is something we do not keep. The only thing 
# we have is the ID & patient number that fall into the 9 
# different groups, i.e., {0 Gy, 0.05 Gy, 2 Gy} x {0, 1, 2}

# filter for the 9 different groups
none_2tumors <- meta %>% filter(Dose == "0 Gy", !is.na(ICCC.3.2.))
none_1tumors <- meta %>% filter(Dose == "0 Gy", !is.na(ICCC.3.1.), is.na(ICCC.3.2.))
none_0tumors <- meta %>% filter(Dose == "0 Gy", is.na(ICCC.3.1.))

low_2tumors <- meta %>% filter(Dose == "0.05 Gy", !is.na(ICCC.3.2.))
low_1tumors <- meta %>% filter(Dose == "0.05 Gy", !is.na(ICCC.3.1.), is.na(ICCC.3.2.))
low_0tumors <- meta %>% filter(Dose == "0.05 Gy", is.na(ICCC.3.1.))

high_2tumors <- meta %>% filter(Dose == "2 Gy", !is.na(ICCC.3.2.))
high_1tumors <- meta %>% filter(Dose == "2 Gy", !is.na(ICCC.3.1.), is.na(ICCC.3.2.))
high_0tumors <- meta %>% filter(Dose == "2 Gy", is.na(ICCC.3.1.))

# provide labels
obs_labels <- meta$ID
gene_labels <- rownames(df2)

# genetic information
X <- list(
  t(df2[,none_0tumors$ID]),
  t(df2[,low_0tumors$ID]),
  t(df2[,high_0tumors$ID]),
  t(df2[,none_1tumors$ID]),
  t(df2[,low_1tumors$ID]),
  t(df2[,high_1tumors$ID]),
  t(df2[,none_2tumors$ID]),
  t(df2[,low_2tumors$ID]),
  t(df2[,high_2tumors$ID])
)

# create meta-data
meta <- data.frame(
  index = 1:9, 
  n_tumors = c(0, 0, 0, 1, 1, 1, 2, 2, 2), 
  radiation_Gy = c(0, .05, 2, 0, .05, 2, 0, .05, 2), 
  radiation = c("none", "low", "high", "none", "low", "high", "none", "low", "high")
)

# build data set
data <- list(
  X = X, 
  meta = meta, 
  gene_labels = gene_labels,
  obs_labels = obs_labels
)


# Save ------------------------------------------------------------------------------
saveRDS(data, "data/KiKmedataset.rds")
labels <- paste(data$meta$n_tumors, data$meta$radiation, sep = "-")
saveRDS(labels, "data/labels.rds")

# remove temporary files
rm(list = ls())
