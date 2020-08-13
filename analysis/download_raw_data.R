####################################################################
# 
# Download raw data from the web
# 
# Note: The GeneHancer database is proprietary and therefore
# must be accessed piecemeal via https://genealacart.genecards.org/.
# 
####################################################################

######################
# Download CRISPR data 
######################

# path to store raw data
raw_data_dir = sprintf("%s/data/raw")

# URL of data
remote = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE120861&format=file&file="

# Gasperini et al results
all_deg_results_filename = "GSE120861_all_deg_results.at_scale.txt"

# names of genes
genes_filename = "GSE120861_at_scale_screen.genes.txt"

# names of cells
cells_filename = "GSE120861_at_scale_screen.cells.txt"

# "reference cells" used by Gasperini et al for computational purposes
reference_cells_filename = "GSE120861_50k_reference_cells.rds"

# all (gRNA, gene) pairs 
gRNAgroup_pair_table_filename = "GSE120861_gene_gRNAgroup_pair_table.at_scale.txt"

# list of gRNA groups used
gRNA_groups_filename = "GSE120861_grna_groups.at_scale.txt"

# Monocle Cell Data Set object with all data 
cds_filename = "GSE120861_at_scale_screen.cds.rds"

# Expression data
expression_filename = "GSE120861_at_scale_screen.exprs.mtx"

# list of files to download
filenames = c( 
              all_deg_results_filename,
              genes_filename,
              cells_filename,
              reference_cells_filename,
              cds_filename,
              expression_filename,
              gRNAgroup_pair_table_filename,
              gRNA_groups_filename
              )

# download files if not already present
for(filename in filenames){
  if(!file.exists(sprintf("%s/%s", raw_data_dir, filename))){
    cat(sprintf("Downloading %s...\n", filename))
    download.file(sprintf("%s%s.gz", remote, filename), sprintf("%s/CRISPR/%s.gz", raw_data_dir, filename))
    gunzip(sprintf("%s/%s.gz", raw_data_dir, filename))
  }
}

# Download supplementary Table S2 from Cell
supplementary_table_file = "https://www.cell.com/cms/10.1016/j.cell.2018.11.029/attachment/7319ccb0-a8c0-45f3-8203-26b9159b0102/mmc2.xlsx"
download.file(supplementary_table_file, sprintf("%s/data/raw/CRISPR/Gasperini_TableS2.xlsx", base_dir))

####################
# Download HI-C data
####################

# URL of data
remote = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63525&format=file&file="
TADs_filename = "GSE63525_K562_Arrowhead_domainlist.txt"
contact_matrices_dirname = "GSE63525_K562_intrachromosomal_contact_matrices"

# download TADs file, if necessary
if(!file.exists(sprintf("%s/HIC/%s", raw_data_dir, TADs_filename))){
  cat(sprintf("Downloading %s...\n", TADs_filename))
  download.file(sprintf("%s%s.gz", remote, TADs_filename), 
                sprintf("%s/HIC/%s.gz", raw_data_dir, TADs_filename))
  gunzip(sprintf("%s/%s.gz", raw_data_dir, TADs_filename))
}

# download contact matrices, if necessary
if(!dir.exists(sprintf("%s/HIC/%s", raw_data_dir, contact_matrices_dirname))){
  cat(sprintf("Downloading %s...\n", contact_matrices_dirname))
  download.file(sprintf("%s%s.tar.gz", remote, contact_matrices_dirname), 
                sprintf("%s/HIC/%s.tar.gz", raw_data_dir, contact_matrices_dirname))
  gunzip(sprintf("%s/%s.tar.gz", raw_data_dir, contact_matrices_dirname))
  untar(sprintf("%s/%s.tar", raw_data_dir, contact_matrices_dirname))
}

########################
# Download ChIP-seq data
########################

# Available courtesy of Shendure Lab at 
# https://drive.google.com/drive/folders/177djZEEPV-udBkdqOtdjjV061O7s8dKj;
# original data can be downloaded at https://www.encodeproject.org/.

##########################
# Download GeneHancer data
##########################

# Available online at https://www.genecards.org/