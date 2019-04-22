library('ggplot2')
library('limma')
#library('tidyr')

# Set ggplot theme for all plots
theme_set(theme_classic())

# Data cleaning and rearrangement takes a while, 
# so done in a separate script and saved R objects to disk.
#source('prepare_metadata.R')
load('STR_data.RData', verbose = T)

table(all_genomecov$project)

# Check concordance of coverage estimates for the RGP data
RGP_meancov = RGP_sample_info[c('collaborator_participant_id','mean_coverage')]
names(RGP_meancov)[1] = 'sample'
RGP_compare_cov = merge(RGP_genomecov, RGP_meancov)

ggplot(data = RGP_compare_cov, aes(y=genomecov, x=mean_coverage)) + 
  geom_abline(intercept = 0, slope = 1, colour = 'blue') + 
  geom_point() + 
  labs(y='median coverage estimated by STRetch/mosdepth', 
       x='mean coverage estimated by GATK/Picard')

# Plot genome coverage for all samples
ggplot(data = all_genomecov, aes(x=genomecov, fill = project)) +
  geom_histogram()
ggplot(data = all_genomecov, aes(x=genomecov, colour = project)) +
  geom_density()

# Check for duplicate samples
data_wide_trans = t(data_wide)
duplicated(data_wide_trans)

# Remove homopolymers
loci_1bp = locus_data$locus[locus_data$repeatunit %in% c('A','C')]
loci_2_6bp = locus_data$locus[!(locus_data$repeatunit %in% c('A','C'))]
#data_wide = data_wide[loci_2_6bp,]

# Get MDS plotting data
mds_plot = plotMDS(data_wide, plot = FALSE)
mds_df = data.frame(mds_plot$cmdscale.out)
# Categorise samples by project and coverage
all_samples = names(data_wide)
RGP_samples = unique(RGP_data$sample)
#mds_df$project = 'other'
#mds_df$project[all_samples %in% RGP_samples] = 'RGP'
#rownames(all_genomecov) = all_genomecov$sample
mds_df$sample = rownames(mds_df)
mds_df = merge(mds_df, all_genomecov)

ggplot(data=mds_df, aes(x=X1, y=X2, colour = project)) + 
  geom_point() +
  labs(x='dim 1', y='dim 2')

ggplot(data=mds_df, aes(x=X1, y=X2, colour = genomecov)) + 
  geom_point() +
  labs(x='dim 1', y='dim 2')

# ggplot(data=mds_df, aes(x=X1, y=X2, shape = project, colour = genomecov)) +
#   geom_point() +
#   labs(x='dim 1', y='dim 2') +
#   scale_colour_gradient(low = "grey", high = "dark blue")

ggplot(data=mds_df[mds_df$project == 'RGP',], 
       aes(x=X1, y=X2, colour = genomecov)) + 
  geom_point() +
  labs(x='dim 1', y='dim 2', title = 'RGP only') +
  scale_colour_gradient(low = "grey", high = "red")

ggplot(data=mds_df[mds_df$project == 'other',], 
       aes(x=X1, y=X2, colour = genomecov)) + 
  geom_point() +
  labs(x='dim 1', y='dim 2', title = 'non-RGP only') +
  scale_colour_gradient(low = "grey", high = "red")
