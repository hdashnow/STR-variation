library('dplyr')
library('ggplot2')
library("ggrepel")

source('STR-variation-functions.R')

### Settings
data_dir = '/group/bioi1/harrietd/STR_variation/working_dir/old/'
#annotation_file = '/group/bioi1/harrietd/STR_variation/working_dir2/RGP_149_1.STRs.annotated.tsv'

# Annotation files
annotation_file = '/group/bioi1/harrietd/STR_variation/working_dir/hg19.simpleRepeat_period1-6_dedup.sorted.annotated.tsv'
gnomad_genes_file = '/group/bioi1/harrietd/STR_variation/annotation_data/gnomad.v2.1.1.lof_metrics.by_gene.txt'
disease_genes_file = '/group/bioi1/harrietd/STR_variation/annotation_data/curated_gene_disease_associations.tsv'

all_annotations = get_annotations(annotation_file, gnomad_genes_file, disease_genes_file)

# Should be 362 samples total across all datasets - need to combine them
### Parse/clean up data
#controls_143.STRs = read.csv(paste0(data_dir,"PCRfreeWGS_143_STRetch_controls.STRs.annotated.tsv"), sep='\t')
RGP_219.STRs = read.csv(paste0(data_dir,"RGP_219.STRs.annotated.tsv"), sep='\t', stringsAsFactors = F)

all_STR_calls = RGP_219.STRs
all_STR_calls = clean_annotations(all_STR_calls) # Set missing values to NA
all_STR_calls$locus = paste(all_STR_calls$chrom, all_STR_calls$start, all_STR_calls$end, sep = '_')

# # Create a version of all_STR_calls with zeros rows removed
# all_STR_calls_nonzero = all_STR_calls[all_STR_calls$locuscoverage > 0,]
# # Extract a smaller subset for plot development
# nrows = dim(all_STR_calls_nonzero)[1]
# STR_subset = all_STR_calls_nonzero[sample(nrows, 10000),]
# 
# RGP_single_sample_STRs = read.csv(annotation_file,
#                                   sep='\t', stringsAsFactors = F)
# RGP_STR_annotations = get_annotations(RGP_single_sample_STRs)

# Remove annotations and replace with new annotations XXX code not working
#all_STR_calls = all_STR_calls[1:16]
#all_STR_calls = merge(all_STR_calls, RGP_STR_annotations, all.x = T)

### Calculate summary stats

# Annotate the two HD loci
HD_indices = all_annotations$pathogenic == "HD_HTT" & !is.na(all_annotations$pathogenic)
all_annotations$pathogenic[HD_indices] = paste(all_annotations$pathogenic[HD_indices], all_annotations$repeatunit[HD_indices])

# Calculate number of significant calls per locus at various thresholds
all_STR_calls %>% 
  group_by(locus) %>% 
  summarise(signif_0_01 = sum(p_adj < 0.01), signif_0_05 = sum(p_adj < 0.05), 
            n_zeros = sum(locuscoverage == 0), 
            mad_locuscoverage_log = mad(locuscoverage_log), # Note locuscoverage_log is already coverage-normalised
            total = length(p_adj)) ->
  STR_calls_locus
# Add annotations to the summary stats
STR_calls_locus = merge(STR_calls_locus, all_annotations, all.x = T)

# Calculate MAD on data with zero-read results removed
all_STR_calls_nonzero %>% 
  group_by(locus) %>% 
  summarise(mad_locuscoverage_log_nonzero = mad(locuscoverage_log)) ->
  STR_calls_nonzero_locus
# Add annotations to the summary stats
STR_calls_nonzero_locus = merge(STR_calls_nonzero_locus, all_annotations, all.x = T)
STR_calls_nonzero_locus$repeatunitlen = nchar(STR_calls_nonzero_locus$repeatunit)


# Plot feature

ggplot(data = STR_subset, # all_STR_calls_nonzero, 
       aes(x = feature, y = locuscoverage_log)) + 
  geom_jitter() +
  geom_violin(colour = 'red') +
  facet_wrap(~repeatunitlen)
ggsave('plots/feature_locuscovlog_subset.jpg')
  
ggplot(data = STR_calls_locus, aes(x = feature, y = signif_0_01)) + 
  geom_jitter() +
  geom_violin(colour = 'red') +
  facet_wrap(~repeatunitlen)
ggsave('plots/feature_signif_all.jpg')

ggplot(data = STR_calls_locus, aes(x = feature, y = n_zeros)) + 
  geom_jitter() +
  geom_violin(colour = 'red') +
  facet_wrap(~repeatunitlen)
ggsave('plots/feature_nzeros_all.jpg')

ggplot(data = STR_calls_locus, aes(x = feature, y = mad_locuscoverage_log)) + 
  geom_jitter() +
  geom_violin(colour = 'red') +
  facet_wrap(~repeatunitlen)
ggsave('plots/feature_madlocuscov_all.jpg')

ggplot(data = STR_calls_nonzero_locus, aes(x = feature, y = mad_locuscoverage_log_nonzero)) + 
  geom_jitter() +
  geom_violin(colour = 'red') +
  facet_wrap(~repeatunitlen)
ggsave('plots/feature_madlocuscov-nonzero_all.jpg')

# Plot TSS

min_tss = min(STR_subset$closest_TSS_distance, na.rm = T)

# ggplot(data = all_STR_calls_nonzero, 
#        aes(x = closest_TSS_distance, y = locuscoverage_log)) + 
#   geom_vline(xintercept = 0, colour = 'blue') +
#   geom_point()
# ggsave('plots/tss_locuscovlog_all.jpg')

ggplot(data = STR_subset, # all_STR_calls_nonzero, 
       aes(x = closest_TSS_distance, y = locuscoverage_log)) + 
  geom_vline(xintercept = 0, colour = 'blue') +
  geom_point()
ggsave('plots/tss_locuscovlog_subset.jpg')

ggplot(data = STR_subset, # all_STR_calls_nonzero, 
       aes(x = closest_TSS_distance, y = locuscoverage_log)) + 
  geom_vline(xintercept = 0, colour = 'blue') +
  geom_point() +
  xlim(min_tss, abs(min_tss))
ggsave('plots/tss_locuscovlog_subset_zoom.jpg')

ggplot(data = STR_subset, # all_STR_calls_nonzero, 
       aes(x = closest_TSS_distance, y = locuscoverage_log)) + 
  geom_vline(xintercept = 0, colour = 'blue') +
  geom_point() +
  facet_wrap(~repeatunitlen, scales = 'free_x')
ggsave('plots/tss_locuscovlog_subset_facet.jpg')

ggplot(data = STR_calls_locus, 
       aes(x = closest_TSS_distance, y = signif_0_01)) + 
  geom_vline(xintercept = 0, colour = 'blue') +
  geom_point()
ggsave('plots/tss_signif_all.jpg')

ggplot(data = STR_calls_locus, 
       aes(x = closest_TSS_distance, y = signif_0_01)) + 
  geom_vline(xintercept = 0, colour = 'blue') +
  geom_point() +
  xlim(min_tss, abs(min_tss))
ggsave('plots/tss_signif_all_zoom.jpg')

ggplot(data = STR_calls_locus, 
       aes(x = closest_TSS_distance, y = signif_0_01)) + 
  geom_vline(xintercept = 0, colour = 'blue') +
  geom_point() + 
  facet_wrap(~repeatunitlen) +
  xlim(min_tss, abs(min_tss))
ggsave('plots/tss_signif_all_zoom_facet.jpg')

ggplot(data = STR_calls_nonzero_locus, 
       aes(x = closest_TSS_distance, y = mad_locuscoverage_log_nonzero)) + 
  geom_vline(xintercept = 0, colour = 'blue') +
  geom_point() + 
  facet_wrap(~repeatunitlen) +
  xlim(min_tss, abs(min_tss))
ggsave('plots/tss_madlocuscov-nonzero_zoom_facet.jpg')


# Plot pathogenic

ggplot(data = subset(STR_calls_locus, !is.na(pathogenic)), 
       aes(y = signif_0_01, x = feature, label=pathogenic, colour=factor(repeatunitlen))) + 
  geom_point(size = 2) +
  geom_text_repel() +
  labs(subtitle = 'N=219 samples',
       y = 'number of significant calls per locus at p<0.01')
ggsave('plots/pathogenic_signif.jpg')

ggplot(data = subset(STR_calls_locus, !is.na(pathogenic)), 
       aes(y = mad_locuscoverage_log, x = feature, label=pathogenic, colour=factor(repeatunitlen))) + 
  geom_point(size = 2) +
  geom_text_repel() +
  labs(subtitle = 'N=219 samples')
ggsave('plots/pathogenic_mad.jpg')

ggplot(data = subset(STR_calls_nonzero_locus, !is.na(pathogenic)), 
       aes(y = mad_locuscoverage_log_nonzero, x = feature, label=pathogenic, colour=factor(repeatunitlen))) + 
  geom_point(size = 2) +
  geom_text_repel() +
  labs(subtitle = 'N=219 samples')
ggsave('plots/pathogenic_mad-nonzero.jpg')

STR_calls_locus$known_pathogenic = !is.na(STR_calls_locus$pathogenic)
STR_calls_nonzero_locus$known_pathogenic = !is.na(STR_calls_nonzero_locus$pathogenic)

ggplot(data = STR_calls_locus, 
       aes(y = mad_locuscoverage_log, x = known_pathogenic)) + 
  geom_violin() +
  labs(subtitle = 'N=219 samples')
ggsave('plots/pathogenic_mad_violin.jpg')

ggplot(data = STR_calls_nonzero_locus, 
       aes(y = mad_locuscoverage_log_nonzero, x = known_pathogenic)) + 
  geom_violin() +
  labs(subtitle = 'N=219 samples')
ggsave('plots/pathogenic_mad-nonzero_violin.jpg')








ggplot(data = STR_calls_locus, aes(y = signif_0_01, x = pathogenic)) + 
  geom_boxplot() +
  labs(subtitle = 'N=219 samples',
       x = "overlaps known pathogenic STR locus", 
       y = 'number of significant calls per locus at p<0.01')

ggplot(data = subset(STR_calls_locus, !is.na(pathogenic)), aes(y = signif_0_01, x = feature, label=pathogenic)) + 
  geom_point() +
  geom_text_repel() +
  labs(subtitle = 'N=219 samples',
       y = 'number of significant calls per locus at p<0.01')



ggplot(data = STR_calls_locus, aes(y = signif_0_01, x = closest_TSS_distance,
                                   group = closest_TSS_distance)) + 
  geom_boxplot() +
  labs(subtitle = 'N=219 samples', 
       x = "distance to closest TSS", 
       y = 'number of significant calls per locus at p<0.01')



### Plots
# ggplot(data = RGP_subset, aes(x = closest_TSS_distance)) + 
#   geom_density() +
#   xlim(-100000, 100000)
# 
# ggplot(data = RGP_subset, aes(y = outlier, x = closest_TSS_distance)) + 
#   geom_boxplot(aes(group = cut_width(closest_TSS_distance, 20000))) +
#   xlim(-100000, 100000)
# 
# ggplot(data = subset(RGP_subset, outlier > 1), aes(y = outlier, x = closest_TSS_distance)) + 
#   geom_point() +
#   xlim(-100000, 100000)
# 
# ggplot(data = RGP_subset, aes(y = outlier, x = abs(closest_TSS_distance))) + 
#   geom_point() +
#   xlim(0, 700000)
# 
# 
# ggplot(data = RGP_subset, aes(y = outlier, x = feature)) + 
#   geom_boxplot() 

ggplot(data = STR_calls_locus, aes(x = signif_0_01, fill = feature)) +
  geom_bar() 

#zoom in
ggplot(data = subset(STR_calls_locus, signif_0_01 > 10), aes(x = signif_0_01, fill = feature)) +
  geom_bar() 


subset(STR_calls_locus, signif_0_01 > 0 & feature == 'CDS')[,c(2,3,4,5,10,11,13,17,18,19)]

ggplot(data = STR_calls_locus, aes(x = signif_0_01)) +
  facet_wrap(~feature) +
  geom_bar() 

ggplot(data = subset(STR_calls_locus, feature = 'CDS'), aes(x = signif_0_01)) +
  geom_bar() 


