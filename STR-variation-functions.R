library('dplyr')

standard_chroms = paste('chr', c(1:22, 'X', 'Y'), sep = '')

clean_STR_calls = function(STRs_annotated) {
  # Remove non-canonical chromosomes
  STRs_annotated = subset(STRs_annotated, chrom %in% standard_chroms)
  # Set missing data to NA
  STRs_annotated[STRs_annotated == ''] = NA
  STRs_annotated[STRs_annotated == 'nan'] = NA
  #STRs_annotated[is.na(STRs_annotated$closest_TSS_transcript_name), 'closest_TSS_distance'] = NaN
  #STRs_annotated$feature[is.na(STRs_annotated$feature)] = 'intergenic'
  STRs_annotated$locuscoverage_log = STRs_annotated$total_assigned_log
  return(STRs_annotated)
}

# get_annotations = function(all_STR_calls) {
#   # Assign a unique locus ID
#   all_STR_calls$locus = paste(all_STR_calls$chrom, all_STR_calls$start, all_STR_calls$end, sep = '_')
#   # Extract annotations only
#   STR_annotations = subset(all_STR_calls, select=-c(4,7:16))
#   STR_annotations = unique(STR_annotations)
#   # Check that the number of rows matches the number of unique locus IDs
#   stopifnot(dim(STR_annotations)[1] == length(unique(all_STR_calls$locus)))
#   return(STR_annotations)
# }

split_dot = function(s) {
  strsplit(s,'.',fixed = T)[[1]][1]
}

get_annotations = function(annotation_file, gnomad_genes_file, disease_genes_file, omim_genes_file) {
  # Parse and clean up annotation data
  all_annotations = read.csv(annotation_file,
                             sep='\t', stringsAsFactors = F)
  # Remove non-canonical chromosomes
  all_annotations = subset(all_annotations, chrom %in% standard_chroms)
  
  all_annotations[all_annotations == ''] = NA # Set empty strings to NA
  all_annotations$locus = paste(all_annotations$chrom, all_annotations$start, all_annotations$end, sep = '_')
  all_annotations$gene_id = sapply(all_annotations$gene_id, split_dot)
  all_annotations$repeatunitlen = nchar(all_annotations$repeatunit)
  all_annotations$known_pathogenic = !is.na(all_annotations$pathogenic)
  
  # Remove loci with repeat unit length greater than 6 bp
  all_annotations = all_annotations[all_annotations$repeatunitlen <= 6,]
  
  # Annotate the two HD loci
  HD_indices = all_annotations$pathogenic == "HD_HTT" & !is.na(all_annotations$pathogenic)
  all_annotations$pathogenic[HD_indices] = paste(all_annotations$pathogenic[HD_indices], all_annotations$repeatunit[HD_indices])
  
  # Set loci with no gene annotation to intergenic, except those on non-standard chromosomes which should be set to NA
  standard_chroms = paste('chr', c(1:22, 'X', 'Y'), sep = '')
  all_annotations$feature[all_annotations$feature == 'nan' & all_annotations$chrom %in% standard_chroms] = 'intergenic'
  all_annotations$feature[all_annotations$feature == 'nan'] = NA
  
  # Parse and clean up gnomad genes
  gnomad_genes = read.csv(gnomad_genes_file, sep='\t', stringsAsFactors = F)
  #keep_columns_gnomad = c('gene',	'transcript', 'chromosome',	'start_position',	'end_position', 'gene_id', 'pLI', 'max_af')
  keep_columns_gnomad = c('gene_id', 'pLI', 'oe_lof', 'oe_lof_upper')
  gnomad_genes = gnomad_genes[,keep_columns_gnomad]
  
  # Parse and clean up disease genes
  disease_genes_table = read.csv(disease_genes_file,
                                 sep='\t', stringsAsFactors = F)
  disease_genes_table = disease_genes_table[disease_genes_table$diseaseType == 'disease',]
  disease_genes = unique(disease_genes_table$geneSymbol)
  
  # Parse and clean up omim genes
  omim_genes_table = read.table(omim_genes_file, sep='\t', stringsAsFactors = F,
                                col.names = c('mim_number', 'mim_type', 'entrez_gene_id', 'hgnc_symbol', 'ensembl_gene_id'))
  omim_genes_table[omim_genes_table == ''] = NA # Replace empty strings with NA
  omim_ensembl_gene_ids = unique(omim_genes_table$ensembl_gene_id)
  omim_ensembl_gene_ids = omim_ensembl_gene_ids[!is.na(omim_ensembl_gene_ids)]
  omim_genes = unique(omim_genes_table$hgnc_symbol)
  omim_genes = omim_genes[!is.na(omim_genes)]
  
  # add disease genes and gnomad pLI, max_af to annotations
  all_annotations$disease_gene = all_annotations$gene_name %in% disease_genes
  all_annotations$in_omim = all_annotations$gene_id %in% omim_ensembl_gene_ids | all_annotations$gene_name %in% omim_genes
  all_annotations = merge(all_annotations, gnomad_genes, all.x = T)
  # Set missing gene names to NA (otherwise will be False)
  all_annotations$disease_gene[is.na(all_annotations$gene_name)] = NA
  all_annotations$in_omim[is.na(all_annotations$gene_id) & is.na(all_annotations$gene_name)] = NA
  
  # Calculate oe_lof_upper and pLI deciles over all genes that overlap STRs
  all_annotations = within(all_annotations, 
    oe_decile <- as.integer(cut(oe_lof_upper, quantile(oe_lof_upper, seq(0, 1, 0.1), na.rm = T), include.lowest=TRUE))
    )
  all_annotations = within(all_annotations, 
    pLI_decile <- as.integer(cut(pLI, quantile(pLI, seq(0, 1, 0.1), na.rm = T), include.lowest=TRUE))
  )
  
  # Check that annotation loci are unique (so that they can be used as unique identifiers)
  stopifnot(length(all_annotations$locus) == length(unique(all_annotations$locus)))
  
  return(all_annotations)
}

summary_stats = function(all_STR_calls, all_annotations) {
  # Use hubers, not huber because huber is just using MAD to estimate scale then using that to estimate mu
  # Note locuscoverage_log is already coverage-normalised
  all_STR_calls %>% group_by(locus) %>%  
    do(data.frame(hubers(.$locuscoverage_log))) ->
    hubers_by_locus
  setnames(hubers_by_locus, old=c("mu","s"), new=c("mu_locuscoverage_log", "s_locuscoverage_log"))
  
  # Calculate number of significant calls per locus at various thresholds
  all_STR_calls %>% 
    group_by(locus) %>% 
    summarise(signif_0_01 = sum(p_adj < 0.01), signif_0_05 = sum(p_adj < 0.05), 
              n_zeros = sum(locuscoverage == 0), 
              total = length(p_adj),
              mean_locuscoverage_log = mean(locuscoverage_log),
              var_locuscoverage_log = var(locuscoverage_log)) ->
    STR_calls_locus
  
  # Add annotations to the summary stats
  STR_calls_locus = merge(STR_calls_locus, hubers_by_locus)
  STR_calls_locus = merge(STR_calls_locus, all_annotations, all.x = T)
  
  return(STR_calls_locus)
}

mark_subsets = function(all_STR_calls){
  all_STR_calls$project = 'other'
  all_STR_calls$project[startsWith(all_STR_calls$sample, 'RGP_')] = 'RGP'
  
  # Mark samples previously identified as highly variable
  hyper_samples = read.table('high_var_samples.txt', stringsAsFactors = F)[[1]]
  all_STR_calls$hyper = all_STR_calls$sample %in% hyper_samples
  
  return(all_STR_calls)
}

read_annotate_calls = function(all_STR_calls_csv, all_annotations) {
  all_STR_calls = read.csv(all_STR_calls_csv, sep='\t', stringsAsFactors = F)
  
  all_STR_calls = mark_subsets(all_STR_calls)
  all_STR_calls = clean_STR_calls(all_STR_calls) # Set missing values to NA
  all_STR_calls$locus = paste(all_STR_calls$chrom, all_STR_calls$start, all_STR_calls$end, sep = '_')
  
  # Annotate all calls
  all_STR_calls = merge(all_STR_calls, all_annotations, all.x = T)
  
  return(all_STR_calls)
}

### How many STR loci show evidence of expansion in our sample?
expansion_stats_locus = function(all_annotated_STRs_locus_calls) {
  # Define evidence of expansion as having at least one read reported for that locus in STRetch:
  all_annotated_STRs_locus_calls$has_counts = !is.na(all_annotated_STRs_locus_calls$total)
  # loci with any counts
  nonzero_loci_n = table(all_annotated_STRs_locus_calls$has_counts)['TRUE']
  nonzero_loci_percent = table(all_annotated_STRs_locus_calls$has_counts)['TRUE']/length(all_annotated_STRs_locus_calls$has_counts) * 100
  
  # Has at least one significant call from STRetch
  # sig at 0.05
  signif_0_05_loci_n = table(all_annotated_STRs_locus_calls$signif_0_05 > 0)['TRUE']
  signif_0_05_loci_percent = table(all_annotated_STRs_locus_calls$signif_0_05 > 0) ['TRUE']/length(all_annotated_STRs_locus_calls$signif_0_05) * 100
  
  # sig at 0.01
  signif_0_01_loci_n = table(all_annotated_STRs_locus_calls$signif_0_01 > 0)['TRUE']
  signif_0_01_loci_percent = table(all_annotated_STRs_locus_calls$signif_0_01 > 0) ['TRUE']/length(all_annotated_STRs_locus_calls$signif_0_01) * 100

  return(setNames(c(nonzero_loci_n, nonzero_loci_percent, 
                    signif_0_05_loci_n, signif_0_05_loci_percent, signif_0_01_loci_n, signif_0_01_loci_percent), 
                  c('nonzero_loci_n','nonzero_loci_percent',
                    'signif_0_05_loci_n','signif_0_05_loci_percent','signif_0_01_loci_n','signif_0_01_loci_percent')))
    
}

expansion_stats_individual = function(all_STR_calls) {
  ### What number/proportion of individuals have STR expansions?
  all_STR_calls %>% 
    group_by(sample, project) %>% 
    dplyr::summarise(signif_0_01 = sum(p_adj < 0.01), signif_0_05 = sum(p_adj < 0.05), 
                     non_zero = sum(locuscoverage > 0), genomecov = mean(genomecov)
    ) -> all_STR_calls_sample
  return(summary(all_STR_calls_sample))
}
