
# source("scripts/pre_figure_data_wrangling.R")

################################################################################
######################## Figures for Manuscript (selected) #####################
################################################################################

## Figure 1
## integrative figure: methods + results (rcc heatmap)


## Figure 2(a)
wv_contribution = ggplot(spectra_gene_0.5, aes(gene_func, fill = wv_range)) +
  geom_bar(position = "fill") +
  scale_fill_manual( values = c(vis="#007f4e",nir="#ff6361",swir="#2c4875"),
                     labels = c(vis="VIS", nir="NIR", swir="SWIR")
                     ) +
  xlab(NULL) +
  ylab("Percentage") +
  scale_x_discrete(labels = 
                     c("ABA" = "ABA-related gene", 
                       "drought" = "drought-related gene", 
                       "pathogen" = "pathogen-related gene", 
                       "photosynthesis" = "photosynthesis-related gene",
                       "pigment" = "pigment-related gene",
                       "water" = "water-related gene")) +
  guides(fill=guide_legend(title = "Wave Range" )) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(size= 20),
    axis.text.x = element_text(
      angle = 45, vjust = 1, hjust = 1,
      size = 20),
  )

pdf(file = file.path(outdir, paste0("wv_contribution_", SPECIES, ".pdf")),width = 8, height = 4.5)
plot(wv_contribution)
graphics.off()

## Figure 2(b): spectra loading
f_abs_spectra_loading = ggplot(spectra_loading_corr2, aes(wv_range, abs_loading, colour = wv_range))+
  geom_violin() + geom_point(position = "jitter", size =0.05) + 
  scale_color_manual(values = c(vis="#007f4e",nir="#ff6361",swir="#2c4875"),
                     labels = c(vis="VIS", nir="NIR", swir="SWIR")
                     ) +
  xlab(NULL) +
  ylab("Magnitude of Loading") +
  labs(colour = "Wave Band Range") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  scale_x_discrete(labels = c("VIS", "NIR", "SWIR")) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size= 20),
    axis.text.y = element_text(size= 20),
    axis.text.x = element_text(size = 20)
  )

pdf(file = file.path(outdir, paste0("f_abs_spectra_loading_", SPECIES, "_2.pdf")),width = 6, height = 3)
plot(f_abs_spectra_loading)
graphics.off()

## Figure S1: gene loading
f_gene_loading = ggplot(gene_loading_func, aes(gene_func, abs_loading, colour = gene_func)) +
  geom_violin() + geom_point(position = "jitter", size = 0.05) + 
  scale_color_manual(values = c(ABA = "#ff6361", drought = "#58508d", pathogen = "#c7522a", 
                                photo = "#adc178", pigment = "#007f4e", water = "#43b0f1"),
                     labels = c("ABA-related gene", "Drought-related gene","Pathogen-related gene",
                                "Photosynthesis-related gene", "Pigment-related gene", "Water-related gene")) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  labs(colour = "Gene Functional Categories") +
  xlab("Gene Functional Category") +
  ylab("Magnitude of Loading") +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  scale_x_discrete(labels = c("ABA", "Drought","Pathogen","Photosynthesis", "Pigment", "Water")) +
  theme_minimal()

pdf(file = file.path(outdir, paste0("f_gene_loading_", SPECIES, "_2.pdf")), width = 10, height = 3)
plot(f_gene_loading)
graphics.off()

## Figure S2 (a) cor > 0.3
###### histogram

#f_cor_brks_low_all = ggplot(df_cor_brks_low, aes(x = cor_perc)) +
#  geom_histogram(binwidth = (1/50)) +
#  theme_bw()

#pdf(file.path(outdir, paste0("0.3_hist_cor_all_",SPECIES,".pdf")), width = 10, height = 5)
#plot(f_cor_brks_low_all)
#dev.off()


f_cor_brks_low = ggplot(df_cor_brks_flt_low, aes(x = cor_perc)) +
  geom_histogram(binwidth = (1/50)) + 
  xlab("Percent of Canonical Correlation Above Threshold (0.3)") +
  ylab("Gene Counts") + 
  theme_bw()

pdf(file.path(outdir, paste0("0.3_hist_cor_",SPECIES,".pdf")), width = 10, height = 5)
plot(f_cor_brks_low)
dev.off()


## Figure S2 (b) > 0.7
## spectra ~ histogram

f_spectra_cor_strong = ggplot(spectr_cor_strong, aes(wv, XY_correlation, colour = gene_func))+
  geom_point(alpha = 1, size = 0.3) +
  scale_x_continuous(limits = c(400, 2400))+
  scale_y_continuous(limits = c(-1, 1)) +
  scale_color_manual(values = c(ABA = "#ff6361", photosynthesis = "#adc178", water = "#43b0f1"),
                     labels = c("ABA-related gene", "photosynthesis-related gene", "water-related gene")) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  labs(colour = "Gene Functional Categories") +
  xlab("Wave Length (nm)") +
  ylab("Canonical Correlation") + 
  theme_minimal() +
  theme(
    axis.title.x = element_text(size= 20),
    axis.title.y = element_text(size= 20),
    axis.text.x = element_text(size= 20),
    axis.text.y = element_text(size = 20)
  )

pdf(file.path(outdir, 
              paste0("0.7_spectra_cor_strong_",SPECIES,".pdf")), width = 12, height = 4)
plot(f_spectra_cor_strong)
dev.off() 

# one_point = spectr_cor_strong %>% group_by(gene) %>% slice(c(200,500))

# f_spectra_cor_strong2 = ggplot(spectr_cor_strong, aes(wv, XY_correlation, colour = gene_func))+
#   geom_point(alpha = 0.7, size = 0.1) + 
#   geom_label(data = one_point, aes(label = gene), vjust = -1, size = 2) + 
#   labs(colour = "Gene Functional Annotation Group") +
#   xlab("Wave Length") +
#   theme_minimal()
  
# pdf(file.path(outdir, 
#                 paste0("0.7_spectra_cor_strong2_",SPECIES,".pdf")), width = 10)
# plot(f_spectra_cor_strong2)
# dev.off()


# f_spectra_cor_strong = ggplot(spectr_cor_strong, aes(wv, XY_correlation, colour = gene))+
#   geom_point(alpha = 0.7, size = 0.1) +
#   scale_color_manual(values = c( ACSA_27016 = "#62cff4",
#                                  MSTRG.4642.1 = "#057dcd", 
#                                  ACSA_21008 = "#809bce", 
#                                  MSTRG.3500.4 = "#29f69f", 
#                                  ACSA_24367 = "#b5ea8c",
#                                  ACSA_21177 = "#007f4e", 
#                                  ACSA_00472 = "#ffadad", 
#                                  ACSA_35147 ="#ff6361")) +
#   guides(color = guide_legend(override.aes = list(size = 6))) +
#   labs(colour = "Gene Functional Annotation Group") +
#   theme_minimal()