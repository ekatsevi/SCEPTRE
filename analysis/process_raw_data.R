#######################################################################
### read in the expression matrix and write to a file-backed matrix ###
#######################################################################

# raw expression matrix file
raw_expression_file = sprintf("%s/data/raw/CRISPR/GSE120861_at_scale_screen.exprs.mtx", base_dir)

# header contains number of genes and cells
header = as.numeric(read.table(raw_expression_file, nrows=1, skip = 1, header=FALSE, fill=TRUE))
num_genes = header[1]
num_cells = header[2]

# set up a file-backed matrix for the transpose, if necessary
expression_transpose_filename = sprintf("%s/data/processed/expression_matrix_transpose", base_dir)
if(!file.exists(sprintf("%s.bk", expression_transpose_filename))){
  # initialize file-backed matrix
  expression_FBM_transpose = FBM(num_genes, num_cells, type = "integer",
                                 backingfile = expression_transpose_filename, 
                                 create_bk = TRUE)
  
  # specify chunk size based on available RAM
  chunk_size = 1000000
  
  # callback function to write to FBM
  write_to_FBM = function(expression_data_chunk, pos){
    cells_in_chunk = expression_data_chunk %>% pull(cell) %>% unique() %>% sort()
    print(sprintf("Processing cells %d to %d...", head(cells_in_chunk,1), tail(cells_in_chunk,1)))
    for(cell in cells_in_chunk){
      expressions_per_cell = expression_data_chunk %>% filter(cell == !!cell)
      expression_FBM_transpose[expressions_per_cell %>% pull(gene),cell] = expressions_per_cell %>% pull(expression)
    }  
  }
  
  # read raw matrix in chunks and write to FBM
  read_delim_chunked(raw_expression_file, delim = " ",
                     SideEffectChunkCallback$new(write_to_FBM), 
                     skip = 2,
                     chunk_size = chunk_size, 
                     col_names = c("gene", "cell", "expression"), 
                     col_types = "iii",
                     progress = FALSE)  
}

# transpose matrix so that columns correspond to genes, if necessary
expression_filename = sprintf("%s/data/processed/expression_matrix", base_dir)
if(!file.exists(sprintf("%s.bk", expression_filename))){
  expression_FBM = big_transpose(expression_FBM_transpose, backingfile = expression_filename)
}

##########################################################################
### read in the gRNA and confounder data and write to a file-backed matrix
##########################################################################

if(!file.exists(sprintf("%s/data/processed/phenoData_selected.bk", base_dir))){
  # read in the cell data set object (takes a few minutes)
  cds = readRDS(sprintf("%s/data/raw/CRISPR/GSE120861_at_scale_screen.cds.rds", base_dir))
  
  # extract phenoData (207324 cells by 6616 features, including 18 covariates and 
  #                    indicators for the 6598 guide RNAs)
  pd_tbl = as_tibble(pData(cds))
  
  # remove a few columns of phenoData that are strings (annoying because they cannot be saved
  # as a file-backed matrix) and recast three string columns we might need (batch effects) 
  # as integer dummy variables
  pd_tbl_selected = pd_tbl %>% 
    select(-c(sample, cell, gene, all_gene, barcode, sample_directory, ko_barcode_file, id)) %>% 
    mutate(prep_batch = factor(prep_batch, labels = 1:length(unique(prep_batch))), 
           within_batch_chip = factor(within_batch_chip, labels = 1:length(unique(within_batch_chip))), 
           within_chip_lane = factor(within_chip_lane, labels = 1:length(unique(within_chip_lane))))
  
  # create file backed matrices out of phenoData
  write_tsv(pd_tbl_selected, sprintf("%s/data/processed/phenoData_selected.tsv", base_dir))
  pd_FBM = big_read(sprintf("%s/data/processed/phenoData_selected.tsv", base_dir), 
                    select = 1:6608,
                    backingfile = sprintf("%s/data/processed/phenoData_selected", base_dir))
  
  # save coefficients of fitted mean-dispersion relationship
  write_tsv(attr(cds@dispFitInfo$blind$disp_func, "coefficients"), 
            sprintf("%s/data/processed/dispCoefficients.tsv", base_dir))
  
  # save the raw dispersions
  write_tsv(cds@dispFitInfo$blind$disp_table, 
            sprintf("%s/data/processed/disp_table.tsv", base_dir))
  
  # write the feature names to file
  write_tsv(as_tibble(colnames(pd_tbl_selected)), 
            sprintf("%s/data/processed/phenoData_selected_colnames.tsv", base_dir), col_names = FALSE) 
}

# calculate number of expressed genes per cell
if(!file.exists(sprintf("%s/data/processed/gene_counts.tsv", base_dir))){
  gene_counts = num_genes - unlist(big_apply(expression_transpose_FBM, 
                                             function(X, ind)(colSums(X[,ind] == 0)), 
                                             ncores = 4))
  write_tsv(tibble(gene_counts), sprintf("%s/data/processed/gene_counts.tsv", base_dir))  
}