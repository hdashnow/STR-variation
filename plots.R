library('ggplot2'); theme_set(theme_classic(base_size = 16))
library('scales')
library('RColorBrewer')
source('prepare_data.R')

#library('plyr') <- having both of these loaded at once causes problems because they both contain summarise
library('dplyr')


# Look at some pathogenic loci
RGP_sig_pathogenic_loci = all_STR_calls[!is.na(all_STR_calls$pathogenic) & all_STR_calls$project == 'RGP' & all_STR_calls$p_adj < 0.05,]


### Some general stats about STRs in the genome

all_annotated_STRs_locus_calls = all_annotated_STRs_locus_calls[all_annotated_STRs_locus_calls$repeatunitlen <= 6,]

# Number of STRs in the STRetch annotation file
dim(all_annotated_STRs_locus_calls)[1]
standard_chroms = paste('chr', c(1:22, 'X', 'Y'), sep = '')
all_annotated_STRs_locus_calls = subset(all_annotated_STRs_locus_calls, chrom %in% standard_chroms)
# Set chromosome plotting order
all_annotated_STRs_locus_calls$chrom = factor(all_annotated_STRs_locus_calls$chrom, levels = standard_chroms)
dim(all_annotated_STRs_locus_calls)[1]

hg19.genome = read.table('/group/bioi1/harrietd/git/STRetch/reference-data/hg19.STRdecoys.sorted.fasta.genome', 
                         col.names = c('chrom', 'size_bp'))
hg19.genome = hg19.genome[hg19.genome$chrom %in% standard_chroms,]
hg19.genome$chrom = factor(hg19.genome$chrom, levels = standard_chroms)

# Plot number of bp in each chromosome
ggplot(data = hg19.genome, aes(x = chrom, y = size_bp)) + geom_bar(stat = "identity")

# plot number of STRs annotated in each chomosome
ggplot(data = all_annotated_STRs_locus_calls,
       aes(x = chrom)) +
  geom_bar(aes(fill = factor(repeatunitlen)))

# plot number of STRs with each repeat unit length
ggplot(data = all_annotated_STRs_locus_calls,
       aes(x = repeatunitlen)) +
  geom_bar() + scale_x_continuous(breaks = 1:6)

# Counts/percent of STRs overlapping CDS/exons
feature_table = table(all_annotated_STRs_locus_calls$feature)
sum(feature_table[c('CDS','exon')])
sum(feature_table[c('CDS','exon')])/sum(feature_table) * 100

# number of STRs in various gene reigons
ggplot(data = all_annotated_STRs_locus_calls,
       aes(x = feature, fill=feature)) +
  geom_bar() +
  facet_wrap(~repeatunitlen) 

# first make a dataframe with frequencies
df <- as.data.frame(with(all_annotated_STRs_locus_calls, table(feature, repeatunitlen)))
# next: compute percentages per group
df <- plyr::ddply(df, plyr::.(repeatunitlen), transform, p = Freq/sum(Freq))
# and plot
ggplot(df, aes(feature, p, fill=feature)) +
  geom_bar(stat = 'identity') + 
  facet_wrap(~repeatunitlen) +
  scale_y_continuous(labels=percent_format()) +
  labs(y = 'percentage of STR loci with that repeat unit length overlapping that feature')

# Counts for the above plots:
table(all_annotated_STRs_locus_calls$feature, all_annotated_STRs_locus_calls$repeatunitlen)

### How many STR loci show evidence of expansion in our sample?

# Define evidence of expansion as having at least one read reported for that locus in STRetch:
all_annotated_STRs_locus_calls$has_counts = !is.na(all_annotated_STRs_locus_calls$total)
table(all_annotated_STRs_locus_calls$has_counts)
table(all_annotated_STRs_locus_calls$has_counts)['TRUE']/length(all_annotated_STRs_locus_calls$has_counts) * 100

# Has at least one significant call from STRetch
table(all_annotated_STRs_locus_calls$signif_0_05 > 0)
table(all_annotated_STRs_locus_calls$signif_0_05 > 0) ['TRUE']/length(all_annotated_STRs_locus_calls$signif_0_05) * 100
table(all_annotated_STRs_locus_calls$signif_0_01 > 0)
table(all_annotated_STRs_locus_calls$signif_0_01 > 0) ['TRUE']/length(all_annotated_STRs_locus_calls$signif_0_01) * 100

### What number/proportion of individuals have STR expansions?
all_STR_calls %>% 
  group_by(sample, project) %>% 
  dplyr::summarise(signif_0_01 = sum(p_adj < 0.01), signif_0_05 = sum(p_adj < 0.05), 
            non_zero = sum(locuscoverage > 0), genomecov = mean(genomecov)
            ) -> all_STR_calls_sample
summary(all_STR_calls_sample)

ggplot(data = all_STR_calls_sample, aes(non_zero, fill = project)) + geom_histogram() + 
  labs(x='Number of STR loci with any reads assigned', y='Number of samples')
ggplot(data = all_STR_calls_sample, aes(signif_0_05, fill = project)) + geom_histogram() + 
  labs(x='Number of STR loci significant a p<0.05', y='Number of samples')
ggplot(data = all_STR_calls_sample, aes(signif_0_01, fill = project)) + geom_histogram() + 
  labs(x='Number of STR loci significant a p<0.01', y='Number of samples')

sequenced_together = read.table('/group/bioi1/harrietd/STR_variation/G84373_samples.txt', stringsAsFactors = F)[[1]]

# A subset of the control samples have a surprisingly large number of expanded STRs
high_var_samples = all_STR_calls_sample[all_STR_calls_sample$non_zero > 10000,]$sample
low_var_samples = all_STR_calls_sample[all_STR_calls_sample$non_zero <= 10000,]$sample

# Write the high and low var sample IDs to file
write.table(high_var_samples,
            file = 'high_var_samples.txt', row.names = F, col.names = F, quote = F)
write.table(low_var_samples,
            file = 'low_var_samples.txt', row.names = F, col.names = F, quote = F)



all_STR_calls_sample[all_STR_calls_sample$non_zero < 10000 & all_STR_calls_sample$project == 'controls',]$sample

tmp_df = cbind(high_var_samples, sequenced_together)

ggplot(data = all_STR_calls_sample, aes(y = non_zero, x = genomecov, colour = project)) + geom_point() + 
  labs(y='Number of STR loci with any reads assigned', y='Sequencing depth')
ggplot(data = all_STR_calls_sample, aes(y = signif_0_01, x = genomecov, colour = project)) + geom_point()

all_STR_calls_sample$seq_depth[all_STR_calls_sample$genomecov < 35] = '<35'
all_STR_calls_sample$seq_depth[all_STR_calls_sample$genomecov >= 35] = '35-40'
all_STR_calls_sample$seq_depth[all_STR_calls_sample$genomecov > 40] = '>40'
# Set plotting order
all_STR_calls_sample$seq_depth = factor(all_STR_calls_sample$seq_depth, levels = c('<35', '35-40', '>40'))

ggplot(data = all_STR_calls_sample, aes(genomecov, fill = seq_depth)) + geom_histogram()
ggplot(data = all_STR_calls_sample, aes(non_zero, fill = seq_depth)) + geom_histogram()


### What number/proportion of loci in feature categories have STR expansions?
all_annotated_STRs_locus_calls %>% 
  group_by(feature) %>% 
  dplyr::summarise(any_signif_0_01 = sum(signif_0_01 > 0, na.rm = T), 
                   any_signif_0_05 = sum(signif_0_05 > 0, na.rm = T), 
                   num_non_zero = sum(has_counts, na.rm = T), total = length(has_counts)
  ) -> all_annotated_STRs_locus_calls_feature


# Are STRs overlapping certain features more/less likely to be significant/expanded?

ggplot(data = all_annotated_STRs_locus_calls_feature, 
       aes(x = feature, y = any_signif_0_01/total)) + 
  geom_bar(stat = 'identity') +
  scale_y_continuous(labels=percent_format()) 

ggplot(data = all_annotated_STRs_locus_calls_feature, 
       aes(x = feature, y = any_signif_0_05/total)) + 
  geom_bar(stat = 'identity') +
  scale_y_continuous(labels=percent_format()) 

ggplot(data = all_annotated_STRs_locus_calls_feature, 
       aes(x = feature, y = num_non_zero/total)) + 
  geom_bar(stat = 'identity') +
  scale_y_continuous(labels=percent_format()) 

# Gene set test
# x = number of significant calls within the category
# m = number of significant genes overall
# n = number of non-siginficant genes overall
# k = number of loci in the category

all_annotated_STRs_locus_calls_feature$num_zero = all_annotated_STRs_locus_calls_feature$total - all_annotated_STRs_locus_calls_feature$num_non_zero
all_annotated_STRs_locus_calls_feature$not_signif_0_01 = all_annotated_STRs_locus_calls_feature$total - all_annotated_STRs_locus_calls_feature$any_signif_0_01
all_annotated_STRs_locus_calls_feature$not_signif_0_05 = all_annotated_STRs_locus_calls_feature$total - all_annotated_STRs_locus_calls_feature$any_signif_0_05

# Are there more significant genes in introns?
gene_set_test = function(df, feature, col_name){
  x = df[[col_name]][df$feature == feature]
  m = sum(df[[col_name]])
  n = sum(df$total) - m
  k = df$total[df$feature == feature]
  
  phyper(x,m,n,k)
}
for (col_name in c("any_signif_0_01", "any_signif_0_05", "num_non_zero","num_zero") ){
  print(paste0('test column: ', col_name))
  for (feature in all_annotated_STRs_locus_calls_feature$feature){
    print(feature)
    print(gene_set_test(all_annotated_STRs_locus_calls_feature, feature, col_name))
  }
}






# Disease gene and constraint measures

ggplot(data = all_annotations,
       aes(x = pLI, y = oe_lof, colour = disease_gene)) +
  geom_point() +
  scale_y_log10() +
  geom_vline(xintercept = 0.9, colour = 'blue') +
  geom_hline(yintercept = 0.35, colour = 'blue')

ggplot(data = all_annotations,
       aes(x = oe_lof, y = oe_lof_upper, colour = disease_gene)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, colour = 'blue') +
  xlim(NA,2) 




# Check that my measures of disease genes are relating to each other in a sensible way.
# In this dataset the 0.35 oe_lof value is in the 4th decile
ggplot(data = all_annotations %>% group_by(oe_decile) %>% summarise(
  disgenet= sum(disease_gene, na.rm = T)/length(disease_gene),
  omim = sum(in_omim, na.rm = T)/length(in_omim)
) %>% melt('oe_decile', variable.name = "database", value.name = "prop_disease_genes")) +
  geom_point(aes(x = oe_decile, y = prop_disease_genes, colour = database))

# pLI doesn't work as well (at least with these gene databases and this particular subet of genes)
ggplot(data = all_annotations %>% group_by(pLI_decile) %>% summarise(
  disgenet= sum(disease_gene, na.rm = T)/length(disease_gene),
  omim = sum(in_omim, na.rm = T)/length(in_omim)
) %>% melt('pLI_decile', variable.name = "database", value.name = "prop_disease_genes")) +
  geom_point(aes(x = pLI_decile, y = prop_disease_genes, colour = database))

# Will use oe_lof_upper or oe_decile from now on


# Plot mean~variance

ggplot(data = STR_calls_locus, 
       aes(x = mu_locuscoverage_log, y = s_locuscoverage_log)) +
  geom_point(aes(colour = n_zeros/total)) #+
#facet_wrap(~repeatunitlen, scales = 'free') +
#geom_abline(intercept = -0.4056, slope = 0.3435, colour = 'grey') +
#geom_abline(intercept = -1.3, slope = 1, colour = 'red') +
#geom_smooth(method = 'lm')
ggsave('plots/locuscov_mean_var.jpg')

lm(s_locuscoverage_log ~ mu_locuscoverage_log, STR_calls_locus)



ggplot(data = STR_calls_locus,
       aes(x = mu_locuscoverage_log, y = s_locuscoverage_log)) +
  geom_point(aes(colour = known_pathogenic), alpha = .1)

ggplot(data = STR_calls_locus,
       aes(x = mu_locuscoverage_log, y = s_locuscoverage_log)) +
  geom_point(aes(colour = disease_gene), alpha = .3)

ggplot(data = subset(STR_calls_locus, signif_0_01 > 0 | (n_zeros/total < 0.90 & mu_locuscoverage_log > 2)),
       aes(x = mu_locuscoverage_log, y = s_locuscoverage_log)) +
  geom_point(aes(colour = signif_0_01))

ggplot(data = subset(STR_calls_locus, n_zeros/total < 0.90 & mu_locuscoverage_log > 2.2),
       aes(x = mu_locuscoverage_log, y = s_locuscoverage_log)) +
  geom_point(aes(colour = signif_0_01))



ggplot(data = STR_calls_locus, 
       aes(x = mu_locuscoverage_log, y = s_locuscoverage_log)) +
  geom_point(colour = 'grey') +
  geom_point(data = subset(STR_calls_locus, known_pathogenic),
             aes(colour = factor(signif_0_01))
             ) + 
  ggsave('plots/locuscov_mean_var_disease.jpg')

ggplot(data = STR_calls_locus, 
       aes(x = mu_locuscoverage_log, y = s_locuscoverage_log)) +
  geom_point(colour = 'grey') +
  geom_point(data = subset(STR_calls_locus, known_pathogenic),
             aes(colour = factor(signif_0_01))
  ) +
  geom_text_repel(data = subset(STR_calls_locus, known_pathogenic & s_locuscoverage_log > 0.14),
                  aes(label=pathogenic)
  ) +
  xlim(NA, 1.9) + 
  ylim(NA, 0.6) + 
  ggsave('plots/locuscov_mean_var_disease_zoom.jpg')

ggplot(data = STR_calls_locus, 
       aes(x = mu_locuscoverage_log, y = s_locuscoverage_log, colour = feature)) +
  geom_point(colour = 'grey') +
  geom_point(data = subset(STR_calls_locus, feature == 'exon')) +
  geom_point(data = subset(STR_calls_locus, feature == 'CDS')) + 
  ggsave('plots/locuscov_mean_var_exons.jpg')
             
   #+
  xlim(NA, 1.9) + 
  ylim(NA, 0.6)


# Plot observed/expected lof

ggplot(data = STR_calls_locus,
       aes(x = oe_lof_upper, colour = in_omim)) +
  geom_density() +
  geom_vline(xintercept = 0.35, colour = 'blue')

ggplot(data = STR_calls_locus,
       aes(x = oe_lof_upper, colour = disease_gene)) +
  geom_density() +
  geom_vline(xintercept = 0.35, colour = 'blue')


# STR_calls_locus %>% group_by(oe_decile) %>% summarise(
#   v = median(mu_locuscoverage_log, na.rm = T)) -> tmp
# tmp = tmp[complete.cases(tmp),]
ggplot(data = subset(STR_calls_locus, feature %in% c('exon', 'CDS')),
       aes(x = oe_decile, y = s_locuscoverage_log, group = oe_decile)) +
  geom_violin() + ylim(NA, 0.15)
#ggsave('plots/oeLOF_signif')

ggplot(data = STR_calls_locus,
       aes(x = oe_lof_upper, y = s_locuscoverage_log, colour = in_omim)) +
  geom_point() +
  xlim(NA, 2.5) +
  geom_vline(xintercept = 0.35, colour = 'blue')

ggplot(data = STR_calls_locus,
       aes(x = oe_lof_upper, y = s_locuscoverage_log, colour = disease)) +
  geom_point() +
  xlim(NA, 2.5) +
  geom_vline(xintercept = 0.35, colour = 'blue')






# Calculate variance (s_locuscoverage_log) deciles over all loci ***may need to increase number of bins***
s_quantiles = quantile(STR_calls_locus$s_locuscoverage_log, seq(0, 1, 0.1), na.rm = T)
STR_calls_locus = within(STR_calls_locus, 
  s_decile <- as.integer(cut(s_locuscoverage_log, quantile(s_locuscoverage_log, seq(0, 1, 0.1), na.rm = T), include.lowest=TRUE))
)

STR_loci_feature = table(STR_calls_locus$feature, useNA = 'always')
# expected counts for each feature:
STR_loci_feature/10

STR_calls_locus$feature = factor(STR_calls_locus$feature)
STR_calls_locus %>% group_by(s_decile, feature) %>% summarise(obs = length(feature)) %>% 
  tidyr::complete(s_decile, feature) -> STR_calls_locus_s_decile_counts
STR_calls_locus_s_decile_counts[is.na(STR_calls_locus_s_decile_counts)] = 0
STR_calls_locus_s_decile_counts$exp = rep(STR_loci_feature/10, 10)


ggplot(data = subset(STR_calls_locus_s_decile_counts, feature %in% c('CDS', 'exon', 'intron', 'intergenic')), 
       aes(y = obs/exp, x = s_decile, colour = feature)) +
  geom_hline(yintercept = 1, colour = 'grey', linetype = "longdash") +
  geom_line() +
  geom_point() +
  scale_color_manual(values=brewer.pal(4, 'Set2')) +
  labs(x = 'decile of variance', 
       y = 'obs/exp loci in this cateogry') +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(breaks = 0:6) +
  ggsave('plots/obs_exp_feature_deciles.jpg')

ggplot(data = subset(STR_calls_locus_s_decile_counts, feature %in% c('CDS', 'exon', 'intron', 'intergenic')), 
       aes(y = obs/exp, x = s_decile, colour = feature)) +
  geom_hline(yintercept = 1, colour = 'grey', linetype = "longdash") +
  geom_line() +
  geom_point() +
  scale_color_manual(values=brewer.pal(4, 'Set2')) +
  labs(x = 'decile of variance', 
       y = 'obs/exp loci in this cateogry') +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(breaks = 0:6, limits = c(NA,2)) +
  ggsave('plots/obs_exp_feature_deciles_zoom.jpg')

# ggplot(data = subset(STR_calls_locus_s_decile_counts, feature %in% c('exon', 'intron', 'intergenic')), 
#        aes(y = obs/exp, x = s_decile, colour = feature)) +
#   geom_hline(yintercept = 1, colour = 'grey') +
#   geom_line() +
#   geom_point() + 
#   scale_color_manual(values=brewer.pal(4, 'Set2')[2:4]) +
#   labs(x = 'decile of variance', 
#        y = 'obs/exp loci in this cateogry') +
#   scale_x_continuous(breaks = 1:10)




#XXX up to here



# Plot feature

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

ggplot(data = STR_calls_locus, aes(x = feature, y = s_locuscoverage_log)) +
  geom_jitter() +
  geom_violin(colour = 'red') +
  facet_wrap(~repeatunitlen)
ggsave('plots/feature_slocuscov_all.jpg')

ggplot(data = STR_calls_locus, aes(x = feature, y = mu_locuscoverage_log)) +
  geom_jitter() +
  geom_violin(colour = 'red') +
  facet_wrap(~repeatunitlen)
ggsave('plots/feature_mulocuscov_all.jpg')

# Plot TSS

min_tss = min(STR_subset$closest_TSS_distance, na.rm = T)


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



# Plot pathogenic

ggplot(data = subset(STR_calls_locus, !is.na(pathogenic)),
       aes(y = signif_0_01, x = feature, label=pathogenic, colour=factor(repeatunitlen))) +
  geom_point(size = 2) +
  geom_text_repel() +
  labs(subtitle = 'N=219 samples',
       y = 'number of significant calls per locus at p<0.01')
ggsave('plots/pathogenic_signif.jpg')

ggplot(data = subset(STR_calls_locus, !is.na(pathogenic)),
       aes(y = s_locuscoverage_log, x = feature, label=pathogenic, colour=factor(repeatunitlen))) +
  geom_point(size = 2) +
  geom_text_repel() +
  labs(subtitle = 'N=219 samples')
ggsave('plots/pathogenic_s.jpg')


ggplot(data = STR_calls_locus,
       aes(y = s_locuscoverage_log, x = known_pathogenic)) +
  geom_violin() +
  labs(subtitle = 'N=219 samples')
ggsave('plots/pathogenic_s_violin.jpg')






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


#subset(STR_calls_locus, signif_0_01 > 0 & feature == 'CDS')[,c(2,3,4,5,10,11,13,17,18,19)]

ggplot(data = STR_calls_locus, aes(x = signif_0_01)) +
  facet_wrap(~feature) +
  geom_bar()

ggplot(data = subset(STR_calls_locus, feature = 'CDS'), aes(x = signif_0_01)) +
  geom_bar()
