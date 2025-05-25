
library(Seurat) 
library(data.table)
library(arrow)

panc_integrated = readRDS("./rawdata/integrated_data.rds")
table(panc_integrated$cell_type,panc_integrated$disease_status)


# ------------------------------------------------------------------------------
# 1. beta cells
# ------------------------------------------------------------------------------

# scale data
panc_beta = subset(panc_integrated,cell_type == "beta")
panc_beta = NormalizeData(panc_beta, normalization.method = "LogNormalize", scale.factor = 10000)
panc_beta = FindVariableFeatures(panc_beta, selection.method = "vst", nfeatures = 2000)
panc_beta = ScaleData(panc_beta,features = rownames(panc_beta))

beta_df = as.data.frame(t(panc_beta@assays$RNA$scale.data))
beta_df$disease = panc_beta$disease_status
beta_df$donor = panc_beta$donor
beta_df$cell_id = rownames(panc_beta@meta.data)
write_parquet(beta_df, "./rawdata/beta_scale.parquet")



# ------------------------------------------------------------------------------
# 2. alpha cells
# ------------------------------------------------------------------------------

# scale data
panc_alpha = subset(panc_integrated,cell_type == "alpha")
panc_alpha = NormalizeData(panc_alpha, normalization.method = "LogNormalize", scale.factor = 10000)
panc_alpha = FindVariableFeatures(panc_alpha, selection.method = "vst", nfeatures = 2000)
panc_alpha = ScaleData(panc_alpha,features = rownames(panc_alpha))

alpha_df = as.data.frame(t(panc_alpha@assays$RNA$scale.data))
alpha_df$disease = panc_alpha$disease_status
alpha_df$donor = panc_alpha$donor
alpha_df$cell_id = rownames(panc_alpha@meta.data)
write_parquet(alpha_df, "./rawdata/alpha_scale.parquet")

# ------------------------------------------------------------------------------
# 3. delta cells
# ------------------------------------------------------------------------------

# scale data
panc_delta = subset(panc_integrated,cell_type == "delta")
panc_delta = NormalizeData(panc_delta, normalization.method = "LogNormalize", scale.factor = 10000)
panc_delta = FindVariableFeatures(panc_delta, selection.method = "vst", nfeatures = 2000)
panc_delta = ScaleData(panc_delta,features = rownames(panc_delta))

delta_df = as.data.frame(t(panc_delta@assays$RNA$scale.data))
delta_df$disease = panc_delta$disease_status
delta_df$donor = panc_delta$donor
delta_df$cell_id = rownames(panc_delta@meta.data)
write_parquet(delta_df, "./rawdata/delta_scale.parquet")

# ------------------------------------------------------------------------------
# 4. gamma cells
# ------------------------------------------------------------------------------

# scale data
panc_gamma = subset(panc_integrated,cell_type == "PP/gamma")
panc_gamma = NormalizeData(panc_gamma, normalization.method = "LogNormalize", scale.factor = 10000)
panc_gamma = FindVariableFeatures(panc_gamma, selection.method = "vst", nfeatures = 2000)
panc_gamma = ScaleData(panc_gamma,features = rownames(panc_gamma))

gamma_df = as.data.frame(t(panc_gamma@assays$RNA$scale.data))
gamma_df$disease = panc_gamma$disease_status
gamma_df$donor = panc_gamma$donor
gamma_df$cell_id = rownames(panc_gamma@meta.data)
write_parquet(gamma_df, "./rawdata/gamma_scale.parquet")


# ------------------------------------------------------------------------------
# 5. epsilon cells
# ------------------------------------------------------------------------------

# scale data
panc_epsilon = subset(panc_integrated,cell_type == "epsilon")
panc_epsilon = NormalizeData(panc_epsilon, normalization.method = "LogNormalize", scale.factor = 10000)
panc_epsilon = FindVariableFeatures(panc_epsilon, selection.method = "vst", nfeatures = 2000)
panc_epsilon = ScaleData(panc_epsilon,features = rownames(panc_epsilon))

epsilon_df = as.data.frame(t(panc_epsilon@assays$RNA$scale.data))
epsilon_df$disease = panc_epsilon$disease_status
epsilon_df$donor = panc_epsilon$donor
epsilon_df$cell_id = rownames(panc_epsilon@meta.data)
write_parquet(epsilon_df, "./rawdata/epsilon_scale.parquet")

