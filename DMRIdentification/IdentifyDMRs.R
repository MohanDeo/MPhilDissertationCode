library(bsseq)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(Homo.sapiens)
library(clusterProfiler)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(gridExtra)
library(rGREAT)

#-------------------METHDIST and LINEAR MODEL-----------------------------------------------------------------------------#

prague_class_df = readRDS()
# Get identified DMRs
setwd()
dss = readRDS()

dmr_coords = dss$DMR

plot_adder = function (plot) {
  rel_size = 5.5
  plot + theme(legend.position = "right", text = element_text(size = rel(rel_size)),
               legend.text = element_text(size = rel(rel_size)), plot.title = element_text(size = rel(rel_size))) +
    guides(color = guide_legend(override.aes = list(linewidth = rel(rel_size), size = rel(rel_size))),
           shape = guide_legend(override.aes = list(size = rel(rel_size))))
}

clean_dmr_coords = function() {
  number_of_coords_to_select = round(0.1 * nrow(dmr_coords))
  most_hypermethylated_dmrs = dmr_coords[order(-dmr_coords$areaStat)[1:number_of_coords_to_select],]
  most_hypomethylated_dmrs = dmr_coords[order(dmr_coords$areaStat)[1:number_of_coords_to_select],]

  dmr_coords = rbind(most_hypermethylated_dmrs, most_hypomethylated_dmrs)
  return(dmr_coords)
}

# We have decided not to clean the dmr coords
#dmr_coords = clean_dmr_coords()

filter_dmr_coords = function(prague_and_meth_data, return_p_values = F, plot_ecdf = F) {
  # First we need to separate the large and small values
  large_meth_data_original = prague_and_meth_data[which(prague_and_meth_data$PragueClass == "Large"),]$MethLevels
  small_meth_data_original = prague_and_meth_data[which(prague_and_meth_data$PragueClass == "Small"),]$MethLevels

  # Now we need to reshape the methylation data
  number_of_dmrs = nrow(dmr_coords)
  large_meth_data = matrix(large_meth_data_original,
                           ncol = number_of_dmrs,
                           nrow = length(
                             unique(prague_and_meth_data$UniqueID[which(prague_and_meth_data$PragueClass == "Large")])
                           ),
                           byrow = T)
  small_meth_data = matrix(small_meth_data_original,
                           ncol = number_of_dmrs,
                           nrow = length(
                             unique(prague_and_meth_data$UniqueID[which(prague_and_meth_data$PragueClass == "Small")])
                           ),
                           byrow = T)

  # Now we need to carry out a KS test column-wise on each matrix
  ks_p_values = numeric(number_of_dmrs)
  wilcox_p_values = numeric(number_of_dmrs)
  t_test_p_values = numeric(number_of_dmrs)
  plot_list = list()
  for (i in 1:ncol(large_meth_data)) {
    large_vec = rep("Large", length(large_meth_data[, i]))
    small_vec = rep("Small", length(small_meth_data[, i]))
    temp_df = data.frame(meth = c(large_meth_data[, i], small_meth_data[, i]), PragueClass = c(large_vec, small_vec))
    wilcox_p_values[i] = wilcox.test(large_meth_data[, i], small_meth_data[, i])$p.value
    t_test_p_values[i] = t.test(large_meth_data[, i], small_meth_data[, i])$p.value
    ks_p_values[i] = ks.test(large_meth_data[, i], small_meth_data[, i])$p.value
    # make an ECDF plot
    if (plot_ecdf == T) {
      plot_list[[i]] = ggplot(data = temp_df, aes(x = meth, color = PragueClass)) + stat_ecdf()
      if (i == 1) {
        print(ggplot(data = temp_df, aes(x = meth, color = PragueClass)) + stat_ecdf())
        print(ggplot(data = temp_df, aes(x = meth, color = PragueClass)) + geom_histogram(fill = "white",
                                                                                          alpha = 0.5, position = "identity"))
        print(qqplot(large_meth_data[, i], small_meth_data[, i]))
        svg("ecdf_plot.svg")
        levels(temp_df$PragueClass) = c("Low", "High")
        ecdf_plot_1 = plot_adder(ggplot(data = temp_df, aes(x = meth, color = PragueClass)) + stat_ecdf(linewidth = 1.0) +
                                   labs(y="Cumulative proportion", x = "Methylation level", color = "Risk group") + 
                                   annotate("text", x = 0.3, y = 0.15, label = paste("p-value:", signif(ks_p_values[i], digits = 4)), 
                                            size = 5.5))
        print(ecdf_plot_1)
        dev.off()
      }
      if (i == 100) {
        print(ggplot(data = temp_df, aes(x = meth, color = PragueClass)) + stat_ecdf())
      }

      if (i == 1000) {
        print(ggplot(data = temp_df, aes(x = meth, color = PragueClass)) + stat_ecdf())
      }
    }
  }
  # Only keep the DMRs which have a p-value less than 0.05 from the Wilcox test
  number_of_samples = length(prague_and_meth_data$UniqueID)
  dmr_regions_to_keep = which(ks_p_values < 0.05)
  prague_and_meth_data_indicies_to_keep = rep(dmr_regions_to_keep, number_of_samples) + (rep((1:number_of_samples - 1) * number_of_dmrs,
                                                                                             each = length(dmr_regions_to_keep)))
  prague_and_meth_data = prague_and_meth_data[prague_and_meth_data_indicies_to_keep,]
  prague_and_meth_data = prague_and_meth_data[complete.cases(prague_and_meth_data),]

  if (return_p_values == T) {
    return(list(ks_p_values, wilcox_p_values, t_test_p_values, prague_and_meth_data))
  }
  return(prague_and_meth_data)
}

# We have now narrowed the smoothed objects down - we only have one which we think fits the best
smoothed_object = readRDS("")

# Reset working directory
setwd()

create_dmr_granges = function(return_strings = F) {
  coordinate_table = tibble()
  for (i in 1:nrow(dmr_coords)) {
    current_dmr = dmr_coords[i,]
    current_dmr = paste0(current_dmr$chr, ":", current_dmr$start, "-", current_dmr$end)
    coordinate_table = rbind(coordinate_table, tibble(id = current_dmr))
  }
  if (return_strings) {
    return(coordinate_table$id)
  }

  coordinate_table = coordinate_table %>%
    mutate(chr = stringr::str_split(id, "[:-]", simplify = TRUE)[, 1],
           start = stringr::str_split(id, "[:-]", simplify = TRUE)[, 2],
           end = stringr::str_split(id, "[:-]", simplify = TRUE)[, 3]) %>%
    makeGRangesFromDataFrame()

  return(coordinate_table)
}

dmr_granges = create_dmr_granges()

methylation_values = getMeth(smoothed_object, dmr_granges, type = "smooth", what = "perRegion") #[[1]][, 1:ncol(smoothed_object)]
colnames(methylation_values) = gsub(".*\\.(.*?)_.*", "\\1", colnames(smoothed_object))

# We then need to join the smoothed object with the prague class df
make_matrix = function() {
  meth_df = as.data.frame(t(methylation_values))
  #colnames(meth_df) = create_dmr_granges(T)
  meth_df = cbind(prague_class_df, meth_df)
  return(meth_df)
}

prague_meth_df = make_matrix()

# We then want to visualise the methylation values across each length group (Small, Medium, Large)
make_methylation_level_distribution_plots_per_region = function(prague_meth_df, cleaned_p_values) {
  reshaped_methylation_values = matrix(prague_meth_df$MethLevels, nrow = length(unique(prague_meth_df$UniqueID)),
                                       ncol = nrow(prague_meth_df) / length((unique(prague_meth_df$UniqueID))),
                                       byrow = T)

  df_to_plot = as.data.frame(reshaped_methylation_values)
  colnames(df_to_plot) = paste0("DMR", 1:ncol(df_to_plot))

  df_to_plot$PragueClass = factor(prague_class_df$prague_class_order)
  df_to_plot$PragueClass[which(df_to_plot$PragueClass == "Medium")] = "Small"

  p_values_as_asterisk = cleaned_p_values
  p_values_as_asterisk = ifelse(p_values_as_asterisk < 0.05, "*", "")

  plot_list = list()
  # The last column in the PragueClass which we don't want to plot
  levels(df_to_plot$PragueClass) = c("Low", "Medium", "High")
  for (i in 1:(ncol(df_to_plot) - 1)) {
    # grob <- grobTree(textGrob(round(cleaned_p_values[i], digits=2), x=0.1,  y=0.95, hjust=0,
    # gp=gpar(col="black", fontsize=13, fontface="italic")))
    plot_list[[i]] = ggplot(df_to_plot, aes(x = PragueClass, y = df_to_plot[, i], fill = PragueClass)) +
      geom_boxplot() +
      xlab("Risk group") +
      ylab(paste("Methylation level (%) in", dmr_granges[i])) +
      labs(fill = 'Risk group') +
      annotate("text", x = 2, y = 1, label = paste("p-value:", signif(cleaned_p_values[i], digits = 4)))
    #print(plot_list[[i]])
  }

  for (i in 1:(ncol(df_to_plot) - 1)) {
    print(plot_list[[i]])
    if (i==1) {
      plot_to_save = plot_list[[i]]
      plot_to_save = plot_adder(plot_to_save)
        svg(paste0("methylationDistributionPerRegion_", i, ".svg"))
        print(plot_to_save)
        dev.off()
    }

    if (i==2) {
      plot_to_save = plot_list[[i]]
      plot_to_save = plot_adder(plot_to_save)
        svg(paste0("methylationDistributionPerRegion_", i, ".svg"))
        print(plot_to_save)
        dev.off()
    }
  }

  #Use ggplot to make a histogram of the cleaned p-values
  bw <- 2 * IQR(cleaned_p_values) / length(cleaned_p_values)^(1/3)
  svg("CleanedDmrPValuesDistribution.svg")
  p_value_hist = ggplot() + aes(cleaned_p_values)+ geom_histogram(binwidth = bw) + xlab("P-value") +
    ylab("Count") + scale_y_continuous(breaks = 0:10, labels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
  p_value_hist = plot_adder(p_value_hist)
  print(p_value_hist)
  dev.off()

  pdf("", height = 40, width = 40)
  #gridExtra::grid.arrange(grobs = plot_list, ncol = 2)
  grid.arrange(grobs = plot_list, ncol = 2)
  #do.call(grid.arrange,plot_list)
  dev.off()
}


make_methylation_level_distribution_plots = function(two_prague_classes, non_filtered_df_return = F, make_overall_plots = F, make_per_region_plots = F, plot_pre_filtering=F) {
  number_of_repeats_for_each_level = nrow(methylation_values)
  prague_class_column = which(colnames(prague_meth_df) == "prague_class_order")

  # The methylation value data frame is joined on the right hand side of the prague class df, so if we only wa
  methylation_values_flattened = c(methylation_values)

  methylation_values_prague_class_df = data.frame(UniqueID = prague_class_df$uniqueID,
                                                  PragueClass = factor(prague_class_df$prague_class_order))
  methylation_values_prague_class_df$'(M)' = prague_class_df$'(M)'
  meth_df = as.data.frame(t(methylation_values))
  colnames(meth_df) = paste0("DMR", 1:ncol(meth_df))

  methylation_values_prague_class_df = cbind(methylation_values_prague_class_df, meth_df)
  saveRDS(methylation_values_prague_class_df, file = "DMRMatrix_output.rds")

  # This matrix is just in an easier format to plot
  methylation_values_prague_class_df_to_plot = methylation_values_prague_class_df[, c(
    which(colnames(methylation_values_prague_class_df) == "PragueClass"),
    grep("DMR", colnames(methylation_values_prague_class_df))
  )]

  methylation_values_prague_class_df_to_plot = data.frame(MethLevels = methylation_values_flattened)
  methylation_values_prague_class_df_to_plot$PragueClass = factor(rep(prague_meth_df[, prague_class_column],
                                                                      each = number_of_repeats_for_each_level))
  methylation_values_prague_class_df_to_plot$'(M)' = rep(prague_meth_df$'(M)', each = number_of_repeats_for_each_level)
  methylation_values_prague_class_df_to_plot$UniqueID = rep(prague_meth_df$uniqueID,
                                                            each = number_of_repeats_for_each_level)
  methylation_values_prague_class_df_to_plot$Case = rep(prague_meth_df$Case, each = number_of_repeats_for_each_level)
  # Make boxplot with violin plot
  # ggsave("MethylationDistributionPlots.png", dpi = 300)

  if (two_prague_classes) {
    methylation_values_prague_class_df_to_plot$PragueClass[which(methylation_values_prague_class_df_to_plot$PragueClass == "Medium")] = "Small"
  }

  if (non_filtered_df_return) {
    return(methylation_values_prague_class_df_to_plot)
  }

  if (plot_pre_filtering & make_overall_plots) {
    filter_dmr_coords(methylation_values_prague_class_df_to_plot, F, T)
    wilcox_test_res = pairwise.wilcox.test(methylation_values_prague_class_df_to_plot$MethLevels, methylation_values_prague_class_df_to_plot$PragueClass, p.adjust.method = "BH")
    wilcox_p_value = wilcox_test_res$p.value
    # The third p-value is NULL as the it corresponds to the comparison between Large and Large
    wilcox_p_value = wilcox_p_value[-3]

    p_value_labels = c("Low and medium", "Low and high", "Medium and high")
    levels(methylation_values_prague_class_df_to_plot$PragueClass) = c("Low", "Medium", "High")
    svg("MethylationBoxPlotDistributionsPreFilter.svg")
    box_plot = ggplot(methylation_values_prague_class_df_to_plot, aes(x = PragueClass, y = MethLevels, fill = PragueClass)) +
      geom_boxplot() +
      xlab("Risk group") +
      ylab("Methylation level") +
      labs(fill = 'Risk group') +
      annotate(geom = "text", x = 2.5, y = c(1.2, 1.15, 1.1), label = paste(p_value_labels, "p-value:", signif(wilcox_p_value, digits = 4)), size = 4.5)
    box_plot = plot_adder(box_plot)
    print(box_plot)
    dev.off()

    return(wilcox_test_res)
  }

  p_values = filter_dmr_coords(methylation_values_prague_class_df_to_plot, T)[[1]]
  meth_df = meth_df[, which(p_values < 0.05)]
  print(which(p_values < 0.05))
  methylation_values_prague_class_df = data.frame(UniqueID = prague_class_df$uniqueID,
                                                  PragueClass = factor(prague_class_df$prague_class_order))
  if (two_prague_classes) {
    methylation_values_prague_class_df$PragueClass[which(methylation_values_prague_class_df$PragueClass == "Medium")] = "Small"
  }
  methylation_values_prague_class_df$'(M)' = prague_class_df$'(M)'
  methylation_values_prague_class_df = cbind(methylation_values_prague_class_df, meth_df)
  print("NUMBER OF COLS")
  print(ncol(methylation_values_prague_class_df))
  print(head(methylation_values_prague_class_df))

  saveRDS(methylation_values_prague_class_df, file = "CleanedDMRMatrix_output.rds")

  methylation_values_prague_class_df_to_plot = filter_dmr_coords(methylation_values_prague_class_df_to_plot)
  p_values = p_values[which(p_values < 0.05)]


  if (make_overall_plots == T) {
    png("MethylationDistributionPlotsViolin.png", width = 10, height = 10, units = "in", res = 300)
    violin_box_plot = ggplot(methylation_values_prague_class_df_to_plot, aes(x = PragueClass, y = MethLevels, fill = PragueClass)) +
      geom_violin(scale = "width") +
      geom_boxplot(width = 0.1, color = "grey", alpha = 0.2) +
      xlab("Prague classification") +
      ylab("Methylation level (%)") +
      labs(fill = 'Prague classification')
    print(violin_box_plot)
    dev.off()

    wilcox_test_res = pairwise.wilcox.test(methylation_values_prague_class_df_to_plot$MethLevels, methylation_values_prague_class_df_to_plot$PragueClass, p.adjust.method = "BH")
    wilcox_p_value = wilcox_test_res$p.value

    svg("MethylationBoxPlotDistributions.svg")
    box_plot = ggplot(methylation_values_prague_class_df_to_plot, aes(x = PragueClass, y = MethLevels, fill = PragueClass)) +
      geom_boxplot() +
      xlab("Prague classification") +
      ylab("Methylation level (%)") +
      labs(fill = 'Prague classification') +
      annotate(geom = "text", x = 2, y = 1.0, label = paste("p-value:", signif(wilcox_p_value, digits = 4))) +
      scale_y_continuous(limits = c(0, 1))
    print(box_plot)
    dev.off()

  }

  if (make_per_region_plots == T) {
    make_methylation_level_distribution_plots_per_region(methylation_values_prague_class_df_to_plot, p_values)
  }
  return(methylation_values_prague_class_df_to_plot)
}

filter_granges = function(granges) {
  prague_meth_data = make_methylation_level_distribution_plots(T, T)
  make_methylation_level_distribution_plots(F, F, T, F, T)
  dmr_p_values = filter_dmr_coords(prague_meth_data, T)

  new_dmr_granges = granges[which(dmr_p_values[[1]] < 0.05)]
  saveRDS(new_dmr_granges, file = "CleanedDmrGranges.rds")
  #colnames(new_dmr_granges) = c("Chromosome", "Start", "End", "Width", "Strand")

  # Save new granges as a .csv file - don't save the strand
  # write.csv(new_dmr_granges[, 1:4], file = "CleanedDmrGranges.csv")
  return(new_dmr_granges)
}

dmr_granges = filter_granges(dmr_granges)

make_scatter_plot_all_dmrs = function(meth_values_subset_1, meth_values_subset_2) {
  # If length of one vector is larger than the other, then take a random sample of the larger one - the smallest vectors
  # are still large enough that a representative value can be taken
  if (length(meth_values_subset_1) < length(meth_values_subset_2)) {
    meth_values_subset_2 = sample(meth_values_subset_2, length(meth_values_subset_1))
  }

  if (length(meth_values_subset_1) > length(meth_values_subset_2)) {
    meth_values_subset_1 = sample(meth_values_subset_1, length(meth_values_subset_2))
  }

  temp_df = data.frame(x = meth_values_subset_1, y = meth_values_subset_2)
  ggscatter(temp_df, x = "x", y = "y",
            add = "reg.line", conf.int = TRUE,
            cor.coef = TRUE, cor.method = "spearman",
            xlab = "Small Prague class methylation levels (%)", ylab = "Prague class methylation levels (%)")

}

make_scatter_plot_per_dmr = function(meth_values_subset_1, meth_values_subset_2, meth_values_subset_1_class = "Small",
                                     meth_values_subset_2_class = "Large") {
  # Similar approach here as in `make_scatter_plot_all_dmrs` - but just take all possible values of the smaller meth
  # values subset
  number_of_meth_value_dmrs = nrow(methylation_values)
  number_of_repeats_1 = length(meth_values_subset_1) / number_of_meth_value_dmrs
  number_of_repeats_2 = length(meth_values_subset_2) / number_of_meth_value_dmrs

  if (number_of_repeats_1 < number_of_repeats_2) {
    meth_values_subset_2 = meth_values_subset_2[1:number_of_repeats_1 * number_of_meth_value_dmrs]
  }
  if (number_of_repeats_1 > number_of_repeats_2) {
    meth_values_subset_1 = meth_values_subset_1[1:number_of_repeats_2 * number_of_meth_value_dmrs]
  }

  temp_df = data.frame(x = meth_values_subset_1, y = meth_values_subset_2)

  ggscatter(temp_df, x = "x", y = "y",
            add = "reg.line", conf.int = TRUE,
            cor.coef = TRUE, cor.method = "spearman",
            xlab = paste(meth_values_subset_1_class, "Prague class methylation levels (%)"),
            ylab = paste(meth_values_subset_2_class, "Prague class methylation levels (%)"))

  # ggplot(temp_df, aes(x = x, y = y)) +
  # geom_density_2d_filled()

}

test_spearman_correlation = function(two_methylation_classes = F) {
  methylation_prague_class_df = make_methylation_level_distribution_plots(two_methylation_classes)

  # if (two_methylation_classes) {
  #   methylation_prague_class_df$PragueClass[which(methylation_prague_class_df$PragueClass == "Medium")] = "Small"
  # }

  small_methylation_levels = methylation_prague_class_df$MethLevels[which(methylation_prague_class_df$PragueClass == "Small")]
  medium_methylation_levels = NA
  if (!two_methylation_classes) {
    medium_methylation_levels = methylation_prague_class_df$MethLevels[which(methylation_prague_class_df$PragueClass == "Medium")]
  }
  large_methylation_levels = methylation_prague_class_df$MethLevels[which(methylation_prague_class_df$PragueClass == "Large")]

  png("SmallPragueClassMethLevelsHist.png", width = 10, height = 10, units = "in", res = 300)
  meth_dist = gghistogram(data.frame(MethLevels = small_methylation_levels), x = "MethLevels") +
    xlab("Methylation levels (%)") +
    ylab("Count")
  print(meth_dist)
  dev.off()

  if (!two_methylation_classes) {
    png("MediumPragueClassMethLevelsHist.png", width = 10, height = 10, units = "in", res = 300)
    meth_dist = gghistogram(data.frame(MethLevels = medium_methylation_levels), x = "MethLevels") +
      xlab("Methylation levels (%)") +
      ylab("Count")
    print(meth_dist)
    dev.off()
  }

  png("LargePragueClassMethLevelsHist.png", width = 10, height = 10, units = "in", res = 300)
  meth_dist = gghistogram(data.frame(MethLevels = large_methylation_levels), x = "MethLevels") +
    xlab("Methylation levels (%)") +
    ylab("Count")
  print(meth_dist)
  dev.off()

  png("OverlaidMethDesnsityPlots.png", width = 10, height = 10, units = "in", res = 300)
  meth_dist = ggplot(methylation_prague_class_df, aes(x = MethLevels, fill = PragueClass)) +
    geom_density(position = "identity", alpha = 0.4) +
    labs(x = "Methylation levels (%)", y = "Density") +
    scale_fill_discrete(name = "Prague Classification")
  print(meth_dist)
  dev.off()


  #make_scatter_plot_per_dmr(small_methylation_levels, large_methylation_levels)
}

test_spearman_correlation(T)


linear_model_prague_length_methylation_level = function(two_methylation_classes = F) {
  methylation_prague_class_df = make_methylation_level_distribution_plots(two_methylation_classes, F, two_methylation_classes, two_methylation_classes)

  # if (two_methylation_classes) {
  #   methylation_prague_class_df$PragueClass[which(methylation_prague_class_df$PragueClass == "Medium")] = "Small"
  # }

  ggplot(methylation_prague_class_df, aes(x = `(M)`, y = MethLevels)) + geom_point()

  prague_class_conf_int_plot = ggplot(methylation_prague_class_df, aes(x = `(M)`, y = MethLevels, color = factor(PragueClass))) +
    geom_point() +
    ylab("Methylation level (%)") +
    xlab("Prague (M) Length (cm)") +
    theme_bw() +
    geom_smooth(method = "lm") +
    scale_fill_discrete(name = "Prague Classification")

  print(prague_class_conf_int_plot)

  linear_model = lm(MethLevels ~ factor(PragueClass) - 1, data = methylation_prague_class_df)
  print(anova(linear_model))

  print(kruskal.test(MethLevels ~ PragueClass, data = methylation_prague_class_df))
  print(pairwise.wilcox.test(methylation_prague_class_df$MethLevels, methylation_prague_class_df$PragueClass, p.adjust.method = "BH"))
}

linear_model_prague_length_methylation_level(F)
linear_model_prague_length_methylation_level(T)

#-------------------------------------GO TERMS-------------------------------------------------------------------------#
get_entrez_ids_from_granges = function(granges_obj) {
  genes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
  overlaps = findOverlaps(genes, granges_obj)
  entrez_ids = elementMetadata(genes)[queryHits(overlaps), "gene_id"]
  print(head(entrez_ids))
  return(as.integer(unlist(entrez_ids)))
}

get_go_terms_from_granges = function(granges_obj) {
  enriched_go_terms = enrichGO(gene = get_entrez_ids_from_granges(granges_obj),
                               OrgDb = org.Hs.eg.db,
                               ont = "MF",
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.2,
                               qvalueCutoff = 0.1,
                               readable = TRUE)

  #pdf("MfAcyclicGraph.pdf", width=10, height = 10)
  acyclic_graph = goplot(enriched_go_terms)
  print(acyclic_graph)
  #dev.off()

  enriched_go_df_mf = data.frame(Description = enriched_go_terms$Description,
                                 p.adjust = enriched_go_terms$p.adjust)

  enriched_go_terms = enrichGO(gene = get_entrez_ids_from_granges(granges_obj),
                               OrgDb = org.Hs.eg.db,
                               ont = "CC",
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.99,
                               qvalueCutoff = 0.99,
                               readable = TRUE)

  #pdf("CcAcyclicGraph.pdf", width=10, height = 10)
  acyclic_graph = goplot(enriched_go_terms)
  print(acyclic_graph)
  #dev.off()

  enriched_go_df_cc = data.frame(Description = enriched_go_terms$Description,
                                 p.adjust = enriched_go_terms$p.adjust)

  enriched_go_terms = enrichGO(gene = get_entrez_ids_from_granges(granges_obj),
                               OrgDb = org.Hs.eg.db,
                               ont = "BP",
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.3,
                               qvalueCutoff = 0.2,
                               readable = TRUE)

  #pdf("BpAcyclicGraph.pdf", width=10, height = 10)#, paper = "a4")
  acyclic_graph = goplot(enriched_go_terms)
  print(acyclic_graph)
  #dev.off()

  enriched_go_df_bp = data.frame(Description = enriched_go_terms$Description,
                                 p.adjust = enriched_go_terms$p.adjust)

  return(list(mf = enriched_go_df_mf, cc = enriched_go_df_cc, bp = enriched_go_df_bp))
}

enriched_go_term_tables = get_go_terms_from_granges(dmr_granges)

print("MF: ")
enriched_go_term_tables[1]
print("CC: ")
enriched_go_term_tables[2]
print("BP: ")
enriched_go_term_tables[3]

coerce_rgreat_annotated_regions_to_csv = function (annotated_granges) {
  annotated_granges_df = as.data.frame(annotated_granges)
  # The annotated_genes column is a list, so each item needs to be unpacked
  annotated_granges_df$annotated_genes = sapply(annotated_granges_df$annotated_genes, paste, collapse = ",")
  annotated_granges_df$dist_to_TSS = sapply(annotated_granges_df$dist_to_TSS, paste, collapse = ",")

  colnames(annotated_granges_df) = c("Chromosome", "Start", "End", "Width", "Strand", "Annotated genes", "Distance to TSS")
  # We don't include the strand information when saving
  write.csv(annotated_granges_df[, c(1:4, 6:7)], file = "rGreatGeneAnnotations.csv", row.names = F)
}

create_region_summary_table_and_plot_graphs = function() {
  great_object = great(dmr_granges, "msigdb:C4:CGN",
                       "TxDb.Hsapiens.UCSC.hg38.knownGene")
  svg("GeneAssociationsGREAT.svg")
  gene_assoc_graph = plotRegionGeneAssociations(great_object)
  print(gene_assoc_graph)
  dev.off()
  table = getRegionGeneAssociations(great_object)
  go_terms_great = getEnrichmentTable(great_object)
  coerce_rgreat_annotated_regions_to_csv(table)
  return(list(regions = table, go = go_terms_great))
}

great_analysis = create_region_summary_table_and_plot_graphs()

print(great_analysis$regions)
print(great_analysis$go)

