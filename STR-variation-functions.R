library('dplyr')

# Set missing data to NA
clean_annotations = function(STRs_annotated) {
  STRs_annotated[STRs_annotated == ''] = NA
  STRs_annotated[STRs_annotated == 'nan'] = NaN
  STRs_annotated[is.na(STRs_annotated$closest_TSS_transcript_name), 'closest_TSS_distance'] = NaN
  STRs_annotated$feature[is.na(STRs_annotated$feature)] = 'intergenic'
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
  all_annotations[all_annotations == ''] = NA # Set empty strings to NA
  all_annotations$locus = paste(all_annotations$chrom, all_annotations$start, all_annotations$end, sep = '_')
  all_annotations$gene_id = sapply(all_annotations$gene_id, split_dot)
  all_annotations$repeatunitlen = nchar(all_annotations$repeatunit)
  all_annotations$known_pathogenic = !is.na(all_annotations$pathogenic)
  
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
  keep_columns_gnomad = c('gene_id', 'pLI', 'oe_lof', 'max_af')
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
  
  # Check that annotation loci are unique (so that they can be used as unique identifiers)
  stopifnot(length(all_annotations$locus) == length(unique(all_annotations$locus)))
  
  return(all_annotations)
}



