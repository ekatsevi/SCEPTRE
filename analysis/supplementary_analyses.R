#########################################################
#
# Carry out all the analyses of the of the association
# results, including the ChIP-seq and HIC enrichment
# analyses.
#
#########################################################

# Read in the association results
original_results = read_tsv(sprintf("%s/results/original_results.tsv", base_dir),
                            col_types = "cclcciiilccciicdddddd")
resampling_results = read_tsv(sprintf("%s/results/resampling_results.tsv", base_dir),
                              col_types = "cclcciiilccciiccddddddd")
likelihood_results = read_tsv(sprintf("%s/results/likelihood_results.tsv", base_dir),
                              col_types = "cccciiilccciiccccddd")
source("analysis/read_data.R")

### Confounding analysis, per gRNA
conditional_effects_filename = sprintf("%s/results/gRNA_confounding.tsv", base_dir)
if(!file.exists(conditional_effects_filename)){
  all_grnas = original_results %>% pull(grna_group) %>% unique()
  conditional_effects = sapply(all_grnas, 
                               function(grna_group){
                                 df = confounders_full
                                 df$grna = phenodata_FBM[,match(grna_group, pd_colnames)]
                                 summary(glm(grna ~ guide_count + total_umis, data = df))$coefficients["total_umis", "t value"]
                               })
  tibble(grna_group = all_grnas, conditional_effect = conditional_effects) %>% 
    write_tsv(conditional_effects_filename)
}

# ChIP-seq enrichment analysis
if(!file.exists(sprintf("%s/results/TF_paired_enhancer_fractions.tsv", base_dir))|
   !file.exists(sprintf("%s/results/TF_enrichments.tsv", base_dir))) {
  # read chipseq data
  important_TFs = c("H3K27ac", "EP300", "BRD4", "GATA2", "TAL1", "TBL1XR1", "DPF2", "RNF2")
  num_TFs = length(important_TFs)
  chipseq_data = vector("list", num_TFs)
  names(chipseq_data) = important_TFs
  for(TF in important_TFs){
    filename = sprintf("%s/data/raw/ChIP-seq/%s.bed", base_dir, TF)
    data = suppressMessages(read_tsv(filename, col_names = FALSE))
    if(ncol(data) %in% c(8,10)){
      data = data[,c(1,2,3,7)]
    } else if(ncol(data) %in% c(5,6)){
      data = data[,c(1,2,3,5)]
    } else{
      cat(sprintf("Skipping %s because score not provided.\n", TF))
      next
    }
    names(data) = c("chrom", "chromStart", "chromEnd", "signalValue")
    data$TF = TF
    chipseq_data[[TF]] = data
  }
  chipseq_data = do.call("rbind", chipseq_data)
  
  # extract which enhancers are paired to genes in original and new analyses
  df_cand_enhancers = original_results %>% 
    filter(site_type == "DHS") %>%
    select(chr, target_site, target_site.start, target_site.stop, rejected) %>%
    group_by(chr, target_site, target_site.start, target_site.stop) %>% 
    summarise(rejected_old = any(rejected)) %>%
    ungroup() %>%
    inner_join(resampling_results %>% 
                 filter(method == "conditional_randomization", site_type == "DHS") %>%
                 select(target_site, rejected) %>%
                 group_by(target_site) %>%
                 summarise(rejected_new = any(rejected)),
               by = "target_site") %>%
    unique()
  
  # GRanges object for chipseq data
  gr_chipseq <- GRanges(
    seqnames = chipseq_data$chrom,
    ranges = IRanges(start = chipseq_data$chromStart, end = chipseq_data$chromEnd),
    score = chipseq_data$signalValue,
    TF = chipseq_data$TF)
  
  # GRanges object for candidate enhancers
  gr_cand = GRanges(
    seqnames = df_cand_enhancers$chr,
    ranges = IRanges(start = df_cand_enhancers$target_site.start, end = df_cand_enhancers$target_site.stop),
    rejected_old = df_cand_enhancers$rejected_old,
    rejected_new = df_cand_enhancers$rejected_new,
    target_site = df_cand_enhancers$target_site)
  
  # Split chipseq data into quintiles
  gr_chipseq_quintiles = gr_chipseq %>% 
    subsetByOverlaps(gr_cand) %>% 
    group_by(TF) %>% 
    mutate(quintile = 1+floor((5-1e-10)*percent_rank(score))) %>% 
    ungroup() 
  
  # Compute which quintile each candidate enhancer falls into
  assign_quintiles = function(TF){
    join_overlap_left(gr_cand, gr_chipseq_quintiles %>% 
                        filter(TF == !!TF)) %>% 
      unique() %>% 
      mutate(TF = !!TF, 
             quintile = ifelse(is.na(quintile), 0, quintile))
  }
  TFs = chipseq_data %>% pull(TF) %>% unique()
  gr_cand_overlaps = do.call("c", lapply(TFs, assign_quintiles))
  
  # compute paired fractions in each quintile
  paired_fractions = gr_cand_overlaps %>% 
    group_by(TF, quintile) %>% 
    summarise(rejected_old = mean(rejected_old), 
              rejected_new = mean(rejected_new)) %>%
    as_tibble()
  write_tsv(paired_fractions, path = sprintf("%s/results/TF_paired_enhancer_fractions.tsv", base_dir))
  
  # compute odds ratios for old and new methods
  old_enrichments = sapply(important_TFs, function(TF){
    enrichment = gr_cand_overlaps %>% 
      filter(TF == !!TF, quintile %in% c(0,5)) %>% 
      as_tibble() %>% 
      select(rejected_old, quintile) %>% 
      table() %>% 
      fisher.test()
    enrichment$estimate
  })
  new_enrichments = sapply(important_TFs, function(TF){
    enrichment = gr_cand_overlaps %>% 
      filter(TF == !!TF, quintile %in% c(0,5)) %>% 
      as_tibble() %>% 
      select(rejected_new, quintile) %>% 
      table() %>% 
      fisher.test()
    enrichment$estimate
  })
  
  TF_enrichments = tibble(TF = important_TFs, old_enrichments, new_enrichments) %>% 
    gather(method, enrichment, -TF) %>% 
    mutate(method = factor(method, 
                           levels = c("old_enrichments", "new_enrichments"), 
                           labels = c("Original", "Proposed")))
  write_tsv(TF_enrichments, path = sprintf("%s/results/TF_enrichments.tsv", base_dir))
}

# HI-C enrichment analysis
if(!file.exists(sprintf("%s/results/rejected_pairs_HIC.tsv", base_dir))){
  domains = read_tsv(sprintf("%s/data/raw/HIC/GSE63525_K562_Arrowhead_domainlist.txt", base_dir), 
                     col_types = "ciiciicddddd") %>% 
    mutate(chr1 = sprintf("chr%s", chr1), chr2 = sprintf("chr%s", chr2))
  
  all_pairs = original_results %>% 
    filter(site_type == "DHS", quality_rank_grna == "top_two") %>%
    select(chr, gene_id, target_gene.start, target_gene.stop, TSS, 
           target_site, target_site.start, target_site.stop, rejected) %>%
    rename(rejected_old = rejected) %>%
    left_join(resampling_results %>% 
                filter(method == "conditional_randomization", 
                       site_type == "DHS", 
                       quality_rank_grna == "top_two") %>%
                select(gene_id,  target_site, rejected) %>%
                rename(rejected_new = rejected),
              by = c("gene_id", "target_site"))
  
  all_enhancers = all_pairs %>% select(target_site, chr, target_site.start, target_site.stop) %>% unique()
  all_genes = all_pairs %>% select(gene_id, chr, target_gene.start, target_gene.stop, TSS) %>% unique()
  
  gr_enhancers = GRanges(
    seqnames = all_enhancers$chr,
    ranges = IRanges(start = all_enhancers$target_site.start, end = all_enhancers$target_site.stop),
    target_site = all_enhancers$target_site)
  
  gr_genes = GRanges(
    seqnames = all_genes$chr,
    ranges = IRanges(start = all_genes$target_gene.start, end = all_genes$target_gene.stop),
    gene_id = all_genes$gene_id)
  
  gr_domains <- GRanges(
    seqnames = domains$chr1,
    ranges = IRanges(start = domains$x1, end = domains$x2),
    domain_id = 1:nrow(domains))
  
  gr_genes = gr_genes %>% join_overlap_left(gr_domains)
  gr_enhancers = gr_enhancers %>% join_overlap_left(gr_domains)
  
  rejected_pairs = all_pairs %>% filter(rejected_old | rejected_new)
  num_rejected_pairs = nrow(rejected_pairs)
  TAD_left = integer(num_rejected_pairs)
  TAD_right = integer(num_rejected_pairs)
  for(pair_idx in 1:num_rejected_pairs){
    print(pair_idx)
    overlapping_domains = intersect(gr_enhancers %>% 
                                      filter(target_site == rejected_pairs$target_site[pair_idx]) %>% 
                                      as_tibble() %>% 
                                      pull(domain_id),
                                    gr_genes %>% 
                                      filter(gene_id == rejected_pairs$gene_id[pair_idx]) %>% 
                                      as_tibble() %>% 
                                      pull(domain_id))
    merged_domain = gr_domains %>% filter(domain_id %in% overlapping_domains) %>% GenomicRanges::reduce() %>% as_tibble()
    if(nrow(merged_domain) > 0){
      TAD_left[pair_idx] = merged_domain$start
      TAD_right[pair_idx] = merged_domain$end    
    } else{
      TAD_left[pair_idx] = NA
      TAD_right[pair_idx] = NA
    }
  }
  rejected_pairs$TAD_left = TAD_left
  rejected_pairs$TAD_right = TAD_right
  
  quality = "MAPQG0"
  resolution = 5000
  resolution_name = "5kb"
  chrs = rejected_pairs %>% filter(!is.na(TAD_left)) %>% pull(chr) %>% unique() %>% sort()
  rejected_pairs_chr_list = vector("list", length(chrs))
  names(rejected_pairs_chr_list) = chrs
  
  for(chr in chrs){
    cat(sprintf("Working on %s...\n", chr))
    cat(sprintf("Reading HI-C data...\n"))
    observed = read_tsv(sprintf("%s/data/raw/HIC/GSE63525_K562_intrachromosomal_contact_matrices/K562/%s_resolution_intrachromosomal/%s/%s/%s_%s.RAWobserved", 
                                base_dir,resolution_name,chr,quality,chr,resolution_name), 
                        col_names = c("Start1", "Start2", "count"), 
                        col_types = "iid")
    
    KRnorm = read_tsv(sprintf("%s/data/raw/HIC/GSE63525_K562_intrachromosomal_contact_matrices/K562/%s_resolution_intrachromosomal/%s/%s/%s_%s.KRnorm", 
                              base_dir,resolution_name,chr,quality,chr,resolution_name), 
                      col_names = "normalization", col_types = "d") %>%
      pull()
    
    rejected_pairs_chr = rejected_pairs %>% 
      filter(chr == !!chr, !is.na(TAD_left)) %>% 
      mutate(enhancer = 0.5*(target_site.start + target_site.stop)) %>%
      select(enhancer, TSS, TAD_left, TAD_right, gene_id, target_site, rejected_old, rejected_new) %>%
      mutate_at(c("enhancer", "TSS", "TAD_left", "TAD_right"), ~floor(./resolution)+1)
    
    observed_normalized = observed %>% 
      mutate(Start1 = Start1/resolution+1, 
             Start2 = Start2/resolution+1,
             score = count/(KRnorm[Start1]*KRnorm[Start2])) %>%
      select(Start1, Start2, score)
    
    num_rejected_pairs = nrow(rejected_pairs_chr)
    score_ranks = numeric(num_rejected_pairs)
    
    for(idx in 1:num_rejected_pairs){
      cat(sprintf("Working on rejected pair %d out of %d...\n", idx, num_rejected_pairs))
      enhancer = rejected_pairs_chr$enhancer[idx]
      TSS = rejected_pairs_chr$TSS[idx]
      TAD_left = rejected_pairs_chr$TAD_left[idx]
      TAD_right = rejected_pairs_chr$TAD_right[idx]
      pair_left = min(enhancer, TSS)
      pair_right = max(enhancer, TSS)
      score_ranks[idx] = observed_normalized %>% 
        filter(Start1 >= TAD_left, Start2 <= TAD_right, Start2 - Start1 == pair_right - pair_left) %>%
        mutate(score_rank = percent_rank(score)) %>% 
        filter(Start1 == pair_left, Start2 == pair_right) %>%
        pull(score_rank)  
    } 
    rejected_pairs_chr$score_rank = score_ranks
    rejected_pairs_chr_list[[chr]] = rejected_pairs_chr
    rm(observed)
    gc()
  }
  
  rejected_pairs_chr = do.call("rbind", rejected_pairs_chr_list)
  
  rejected_pairs = rejected_pairs %>% 
    left_join(rejected_pairs_chr %>% 
                select(gene_id, target_site, score_rank), 
              by = c("gene_id", "target_site"))
  write_tsv(rejected_pairs, sprintf("%s/results/rejected_pairs_HIC.tsv", base_dir))
}

# matching with GeneHancer database
if(!file.exists(sprintf("%s/results/GeneHancer_annotated_rejections.tsv", base_dir))){
  gene_hancer_genes = 
    rbind(read_excel(sprintf("%s/data/raw/GeneHancer/GeneHancer.xlsx", base_dir), sheet = "Genomics"),
          read_excel(sprintf("%s/data/raw/GeneHancer/GeneHancer_2.xlsx", base_dir), sheet = "Genomics")) %>%
    group_by(InputTerm) %>% 
    mutate(TSS = ifelse(`Strand (GRCh37/hg19)` == "Plus", 
                        strsplit(strsplit(`Locations (GRCh37/hg19)`, split = ":")[[1]][1], split = "\\(")[[1]][2], 
                        strsplit(strsplit(`Locations (GRCh37/hg19)`, split = ":")[[1]][2], split = "\\)")[[1]][1])) %>% 
    ungroup() %>% 
    select(InputTerm, Symbol, `Chromosome (GRCh37/hg19)`, TSS)
  
  gene_hancer_aliases = rbind(read_excel(sprintf("%s/data/raw/GeneHancer/GeneHancer.xlsx", base_dir), sheet = "Aliases"),
                              read_excel(sprintf("%s/data/raw/GeneHancer/GeneHancer_2.xlsx", base_dir), sheet = "Aliases"))
  
  gene_alias_list = lapply(1:nrow(gene_hancer_aliases), function(idx)(c(gene_hancer_aliases$Symbol[idx], 
                                                                        strsplit(gene_hancer_aliases$Alias[idx], split = "\\|\\|")[[1]])))
  names(gene_alias_list) = gene_hancer_aliases$Symbol
  
  gene_hancer_enhancers = rbind(read_excel(sprintf("%s/data/raw/GeneHancer/GeneHancer.xlsx", base_dir), sheet = "GeneHancers"),
                                read_excel(sprintf("%s/data/raw/GeneHancer/GeneHancer_2.xlsx", base_dir), sheet = "GeneHancers")) %>%
    left_join(gene_hancer_genes, by = c("InputTerm", "Symbol")) %>%
    mutate(location = as.numeric(TSS) + as.numeric(TSSdistance)) %>%
    rename(gene_id = InputTerm, gene_short_name = Symbol, chrom = `Chromosome (GRCh37/hg19)`) %>%
    mutate(chrom = sprintf("chr%s", chrom))
  
  rejections = resampling_results %>% 
    filter(site_type == "DHS", quality_rank_grna == "top_two", method == "conditional_randomization", rejected) %>%
    select(gene_id, gene_short_name, grna_group, chr, target_site.start, target_site.stop) 
  
  get_GHID = function(grna_group_chr, gene_short_name, grna_group_start, grna_group_end){
    matching = gene_hancer_enhancers %>% 
      group_by(gene_short_name) %>%
      filter(chrom == grna_group_chr, 
             !!gene_short_name %in% gene_alias_list[[unique(gene_short_name)]],
             location > grna_group_start - 2000, 
             location < grna_group_end + 2000) %>%
      ungroup()
    if(nrow(matching) == 0){
      return(NA)
    } else{
      matching %>% arrange(abs(location - grna_group_start)) %>% head(1) %>% pull(GHIDs)
    }
  }
  
  annotated_rejections = rejections %>% 
    group_by(gene_short_name, grna_group) %>% 
    mutate(GHID = get_GHID(unique(chr), unique(gene_short_name), 
                           unique(target_site.start), unique(target_site.stop))) %>%
    ungroup() %>% 
    left_join(gene_hancer_enhancers %>% select(GHIDs, Score, GHtype, Sources, NumGenesAway, location) %>% 
                dplyr::rename(GHID = GHIDs, GH_Score = Score), 
              by = c("GHID")) %>%
    select(chr, gene_short_name, grna_group, GHID, GH_Score, GHtype, Sources, NumGenesAway, target_site.start, target_site.stop, location) 
  
  write_tsv(annotated_rejections, sprintf("%s/results/GeneHancer_annotated_rejections.tsv", base_dir))  
}