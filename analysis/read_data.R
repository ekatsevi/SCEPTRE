####################################################
#
# Read in the CRISPR screen raw data
#
####################################################

# base directory should already be in the workspace when this script is called
# base_dir = "/home/ekatsevi/project-files/gene-enhancer/"

# cell "phenotype" names (phenotypes include including guide counts and confounders)
pd_colnames = read_tsv(sprintf("%s/data/processed/phenoData_selected_colnames.tsv", base_dir), 
                       col_names = FALSE, col_types = "c") %>% pull()
num_phenotypes = length(pd_colnames)
num_grnas = 6598

# cell names
cell_names = read_tsv(sprintf("%s/data/raw/CRISPR/GSE120861_at_scale_screen.cells.txt", base_dir), 
                      col_names = FALSE, col_types = "c") %>% pull()
num_cells = length(cell_names)

# gene names
gene_ids = read_tsv(sprintf("%s/data/raw/CRISPR/GSE120861_at_scale_screen.genes.txt", base_dir), 
                    col_names = FALSE, col_types = "c") %>% pull()
num_genes = length(gene_ids)

# gene expressions (as a file-backed matrix)
expression_FBM = FBM(num_cells, num_genes, type = "integer",
                     backingfile = sprintf("%s/data/processed/expression_matrix", base_dir), 
                     create_bk = FALSE)

# gene expressions (as a file-backed matrix)
# expression_transpose_FBM = FBM(num_genes, num_cells, type = "integer",
#                      backingfile = sprintf("%s/data/processed/expression_matrix_transpose", base_dir),
#                      create_bk = FALSE)


# "phenotypes" for each cell, including guide counts and confounders (as a file-backed matrix)
phenodata_FBM = FBM(num_cells, num_phenotypes, 
                    backingfile = sprintf("%s/data/processed/phenoData_selected", base_dir), 
                    create_bk = FALSE)

# phenodata_transpose_FBM = FBM(num_phenotypes, num_cells,
#                     backingfile = sprintf("%s/data/processed/phenoData_selected_transpose", base_dir), 
#                     create_bk = FALSE)

# reference cells used by Gasperini et al for computational purposes
reference_cells = readRDS(sprintf("%s/data/raw/CRISPR/GSE120861_50k_reference_cells.rds", base_dir))

# # dispersion coefficients alpha_0 and alpha_1 from equation (6) in DESeq2 paper
disp_coeffs = as.numeric(read_tsv(sprintf("%s/data/processed/dispCoefficients.tsv", base_dir), col_types = "dd"))
# 
# ## dispersion table (mean expressions and dispersion estimates for each gene)
disp_table = read_tsv(sprintf("%s/data/processed/dispTable.tsv", base_dir), col_types = "cdd")

# extract confounders
confounder_names = c("guide_count", "prep_batch", "percent.mito", "total_umis", "umi_count", "read_count")
confounders_full = phenodata_FBM[,match(confounder_names, pd_colnames)]
colnames(confounders_full) = confounder_names
confounders_full = as.data.frame(confounders_full)
confounders_full$gene_count = read_tsv(sprintf("%s/data/processed/gene_counts.tsv", base_dir), col_types = "i") %>% pull()

reference_cells_bool = cell_names %in% reference_cells

# extract size factors
Size_Factors_full = phenodata_FBM[,match("Size_Factor", pd_colnames)]
confounders_full$Size_Factor = Size_Factors_full