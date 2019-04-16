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

get_annotations = function(annotation_file, gnomad_genes_file, disease_genes_file) {
  # Parse and clean up annotation data
  all_annotations = read.csv(annotation_file,
                             sep='\t', stringsAsFactors = F)
  all_annotations$locus = paste(all_annotations$chrom, all_annotations$start, all_annotations$end, sep = '_')
  all_annotations$gene_id = sapply(all_annotations$gene_id, split_dot)
  all_annotations$repeatunitlen = nchar(all_annotations$repeatunit)
  # Set loci with no gene annotation to intergenic, except those on non-standard chromosomes which should be set to NA
  standard_chroms = paste('chr', c(1:22, 'X', 'Y'), sep = '')
  all_annotations$feature[all_annotations$feature == 'nan' & all_annotations$chrom %in% standard_chroms] = 'intergenic'
  all_annotations$feature[all_annotations$feature == 'nan'] = NA
  
  # Parse and clean up gnomad genes
  gnomad_genes = read.csv(gnomad_genes_file, sep='\t', stringsAsFactors = F)
  #keep_columns_gnomad = c('gene',	'transcript', 'chromosome',	'start_position',	'end_position', 'gene_id', 'pLI', 'max_af')
  keep_columns_gnomad = c('gene_id', 'pLI', 'max_af')
  gnomad_genes = gnomad_genes[,keep_columns_gnomad]
  
  # Parse and clean up disease genes
  disease_genes_table = read.csv(disease_genes_file,
                                 sep='\t', stringsAsFactors = F)
  disease_genes_table = disease_genes_table[disease_genes_table$diseaseType == 'disease',]
  disease_genes = unique(disease_genes_table$geneSymbol)
  
  # add disease genes and gnomad pLI, max_af to annotations
  all_annotations$disease_gene = all_annotations$gene_name %in% disease_genes
  all_annotations = merge(all_annotations, gnomad_genes, all.x = T)
  
  # Check that annotation loci are unique (so that they can be used as unique identifiers)
  stopifnot(length(all_annotations$locus) == length(unique(all_annotations$locus)))
  
  return(all_annotations)
}



