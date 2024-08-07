library("data.table")
library(dplyr)
library(ggplot2)

#------------------HELPER FUNCTIONS-----------------------------------------#
'%!in%' <- function(x,y)!('%in%'(x,y))

create_elbow_plot = function(unsorted_data, y_label, threshold_value) {
  # sorted_data = sort(unsorted_data, decreasing = T)
  # 
  # plot(1:length(sorted_data), sorted_data, ylab = y_label, xlab="Sorted index")
  # abline(h=threshold_value)
  
  # Sort the data in decreasing order
  sorted_data = sort(unsorted_data, decreasing = TRUE)
  
  # Create a data frame for ggplot
  df = data.frame(
    index = 1:length(sorted_data),
    value = sorted_data
  )
  
  # Create the ggplot
  p <- ggplot(df, aes(x = index, y = value)) +
    geom_point() +
    labs(y = y_label, x = "Sorted index") +
    geom_hline(yintercept = threshold_value, linetype = "dashed", color = "red",
               linewidth = 1.5) +
    theme_minimal()
  
  return(plot_adder(p))
  
}

integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

plot_adder = function (plot) {
  rel_size = 5.5
  plot + theme(legend.position = "right", text = element_text(size = rel(rel_size)),
               legend.text = element_text(size = rel(rel_size)), plot.title = element_text(size = rel(rel_size)),
               strip.text.y = element_text(size=rel(rel_size)), strip.text.x = element_text(size=rel(rel_size))) +
    guides(color = guide_legend(override.aes = list(linewidth = rel(rel_size), size = rel(rel_size))),
           shape = guide_legend(override.aes = list(size = rel(rel_size))))
}

#-----------------------READ IN DATA-------------------------------------------#
hsmetrics_summary <- function(l, dataBatch = dataBatch){
  # Databatch is the name of the dataBatch you use
  hsmetrics.list <- lapply(l, function(x) { fread(cmd=paste0("grep '## METRICS CLASS	picard.analysis.directed.HsMetrics' -A 2 ", x)) })
  PicardCollectHSmetrics.summary <- do.call(rbind,
                                            c(hsmetrics.list, list(fill = TRUE)))
  print(head(PicardCollectHSmetrics.summary))
  PicardCollectHSmetrics.summary <- PicardCollectHSmetrics.summary %>%
    dplyr::mutate(sampleID = sub("^(S\\d+).*", "\\1", sub(".*/HSmetrics/([^/]+)*", "\\1", l)),
                  plateID = dataBatch
    )
  invisible(PicardCollectHSmetrics.summary)
}

# We need to be where the data is
setwd()

# The samples are split into different directories - we don't want the SLX-23587
# or the CYTRD0004_20240104_FFPE directory though and the directory itself
directory_names = list.dirs(getwd())
directory_names = directory_names[-grep("SLX-23587", directory_names)]
directory_names = directory_names[-grep("CYTRD0004_20240104_FFPE", 
                                        directory_names)]

directory_names = directory_names[-which(directory_names==getwd())]

concatenated_tables = c()

# We need to make a column so that we can join them later down the line

for(directory in directory_names) {
  filenames = list.files(
    directory, pattern="_raw"
  )
  plate_number = which(directory_names==directory)
  colnames1=c()
  colnames2=c()
  setwd(directory)
  sample_table = hsmetrics_summary(filenames, directory)
  # Some columns need to be filled as not all files have the same number of columns
  if (directory==""){
    colnames1 = c(colnames1, colnames(sample_table))
  }
  if (directory==""){
    colnames2 = colnames(sample_table)
    # print(colnames1[which(colnames1%!in%colnames2)])
  }
  sample_table$uniqueID = paste0("P", plate_number, "_", sample_table$sampleID)
  concatenated_tables = rbind(concatenated_tables, sample_table, fill=T)
}

# The demographics data in a different directory
setwd()
demographics_data = readRDS("ALLPLATES_CaseDemographics.Rds")

#---------------------------SET SCREE THRESHOLDS--------------------------------#
ten_x_cov_threshold = 0.85
dup_rate_threshold = 0.2
off_target_threshold = 0.15
at_dropout_threshold = 17
gc_droput_threshold = 3

thresholds = c(ten_x_cov_threshold, dup_rate_threshold, off_target_threshold, 
               at_dropout_threshold, gc_droput_threshold)

# cov10X = Threshold of desired 10X coverage percentage, default = 0.5
percent_coverages = c(concatenated_tables$PCT_TARGET_BASES_10X,
                      concatenated_tables$PCT_TARGET_BASES_20X,
                      concatenated_tables$PCT_TARGET_BASES_30X)
coverage_labels = c(rep("10X", length(concatenated_tables$PCT_TARGET_BASES_10X)))#,
                    # rep("20X", length(concatenated_tables$PCT_TARGET_BASES_20X)),
                    # rep("30X", length(concatenated_tables$PCT_TARGET_BASES_30X)))
percent_coverages = data.frame(Coverage = coverage_labels, Value = percent_coverages)


#-------------------------------MAKE SEQUENCING PLOTS--------------------------#
setwd()
svg("percent_coverage_box.svg")
plot_adder(ggplot(data = percent_coverages, aes(x=Coverage, y=Value)) +
  geom_boxplot() + 
  geom_jitter(color="black", size=0.4, alpha=0.9) + ylab("Bases with given coverage (%)"))
dev.off()

svg("percent_coverage_jitter.svg")
create_elbow_plot(concatenated_tables$PCT_TARGET_BASES_10X, "Coverage (%)",
                  ten_x_cov_threshold)
dev.off()

# dupRate = Threshold of desidered duplication rate, default = 0.4
duplication_rates = data.frame(dup_rates = concatenated_tables$PCT_EXC_DUPE)
svg("duplication_rates_box.svg")
plot_adder(ggplot(data = duplication_rates, aes(y=dup_rates, x = 0)) +
  geom_boxplot() + 
  geom_jitter(color="black", size=0.4, alpha=0.9) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + labs(y="Duplication rate (%)"))
dev.off()

svg("duplication_rates_jitter.svg")
create_elbow_plot(duplication_rates[,1], "Duplication rate (%)", dup_rate_threshold)
dev.off()

# offTarget = Threshold of desidered off-target percentage, default = 0.2
off_target = data.frame(off_tgt = concatenated_tables$PCT_EXC_OFF_TARGET)
svg("off_target_box.svg")
plot_adder(ggplot(data = off_target, aes(y=off_tgt, x = 0)) +
  geom_boxplot() + 
  geom_jitter(color="black", size=0.4, alpha=0.9) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + labs(y="Off target (%)"))
dev.off()

svg("off_target_jitter.svg")
create_elbow_plot(off_target[,1], "Off target (%)", off_target_threshold)
dev.off()

# ATdropout = Threshold of A/T dropout proportion allowed, default = NA
at_dropouts = data.frame(at_drop = concatenated_tables$AT_DROPOUT)
svg("at_dropouts_hist.svg")
plot_adder(ggplot(at_dropouts, aes(x = at_drop)) +
  geom_histogram() +
  labs(x = "A/T dropouts (%)", y = "Frequency"))
dev.off()

svg("at_dropouts_box.svg")
plot_adder(ggplot(data = at_dropouts, aes(y = at_drop, x = 0)) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = "A/T dropouts (%)"))
dev.off()

svg("at_dropouts_jitter.svg")
create_elbow_plot(at_dropouts[, 1], "A/T dropouts (%)", at_dropout_threshold)
dev.off()

# GCdropout = Threshold of G/C dropout proportion allowed, default = NA
gc_dropouts = data.frame(gc_drop = concatenated_tables$GC_DROPOUT)
svg("gc_dropouts_hist.svg")
plot_adder(ggplot(gc_dropouts, aes(x = gc_drop)) +
  geom_histogram() +
  labs(x = "G/C dropouts (%)", y = "Frequency"))
dev.off()

svg("gc_dropouts_box.svg")
plot_adder(ggplot(data = gc_dropouts, aes(y = gc_drop, x = 0)) +
  geom_boxplot() +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = "G/C dropouts (%)"))
dev.off()

svg("gc_dropouts_jitter.svg")
create_elbow_plot(gc_dropouts[, 1], "G/C dropouts (%)", gc_droput_threshold)
dev.off()


#-----------------CLEAN SEQUENCING DATA----------------------------------------#
get_outlier_indicies = function(data, threshold, column_of_interest = 1, 
                                minimum_threshold=F) {
  # All input data will be a numeric data frame
  # Operator needs to change here
  if(minimum_threshold==T){
    return(which(data[,column_of_interest] < threshold))
  }
  return(which(data[,column_of_interest] > threshold))
}

IQR.outliers <- function(x) {
  if(any(is.na(x)))
    stop("x is missing values")
  if(!is.numeric(x))
    stop("x is not numeric")
  Q3<-quantile(x,0.75)
  Q1<-quantile(x,0.25)
  IQR<-(Q3-Q1)
  left<- (Q1-(1.5*IQR))
  right<- (Q3+(1.5*IQR))
  print(left)
  print(right)
  c(which(x <left), which(x>right))
}

get_cleaned_data_indicies = function() {
  # For some metrics eg. coverage, we want to threshold at a minimum value
  # and in others, we want to threshold at a maximum value
  all_dataframes = list(data.frame(ten_x_cov=percent_coverages[,2][which(percent_coverages[,1]=='10X')]), duplication_rates, off_target,
                        at_dropouts, gc_dropouts)
  cleaned_indicies = rownames(concatenated_tables)
  for(i in 1:length(all_dataframes)) {
    outlier_indicies = NA
    if (colnames(all_dataframes[[i]])=="ten_x_cov"){
      # We want to threshold for a minimum value for percentage coverage
      outlier_indicies = get_outlier_indicies(all_dataframes[[i]],
                                              thresholds[i], minimum_threshold=T)
    } else {
      # outlier_indicies = get_outlier_indicies(all_dataframes[[i]], thresholds[i])
      outlier_indicies = IQR.outliers(all_dataframes[[i]][,1])
    }
    # outlier_indicies = IQR.outliers(all_dataframes[[i]][,1])
    cleaned_indicies = cleaned_indicies[cleaned_indicies%!in%outlier_indicies]#cleaned_indicies[-which(cleaned_indicies%in%outlier_indicies)]
  }
  return(as.numeric(cleaned_indicies))
}

concatenated_tables_cleaned = concatenated_tables[get_cleaned_data_indicies(), ]

case_demo_data = demographics_data$CaseDemoComplete

collate_demographics_data_from_all_plates = function (demographics_split_by_plate) {
  plate_names = rownames(summary(demographics_split_by_plate))
  concatenated_demographics_data = c()
  
  for (plate_name in plate_names) {
    # We need to make the plate numbers again, so we can annotate each sample individually
    plate_number = which(plate_names == plate_name)
    plate_data = demographics_split_by_plate[[plate_name]]
    plate_data$uniqueID = paste0("P", plate_number, "_", plate_data$sample)
    concatenated_demographics_data = rbind(concatenated_demographics_data, plate_data)
  }
  # Now we want to apply the same cleaning to demographics data
  concatenated_demographics_data_indicies_to_keep = which(concatenated_demographics_data$uniqueID%in%concatenated_tables_cleaned$uniqueID)
  concatenated_demographics_data = concatenated_demographics_data[concatenated_demographics_data_indicies_to_keep, ]
  return(concatenated_demographics_data)
}

case_demo_data = collate_demographics_data_from_all_plates(case_demo_data)



#-----------------------GROUP BY PRAGUE LENGTH----------------------------------#
group_samples_by_prague_classification = function(demographics){
  # We don't want to include samples with NA as their prague length
  demographics = demographics[!is.na(demographics$`(M)`), ]
  prague_lengths = demographics$`(M)`
  
  small_lower_limit = 0
  small_upper_limit = 3
  medium_lower_limit = 4
  medium_upper_limit = 6
  large_lower_limit = 7
  large_upper_limit = 100000000
  
  small_length_indicies = prague_lengths >= small_lower_limit & prague_lengths <= small_upper_limit
  medium_length_indicies = prague_lengths >= medium_lower_limit & prague_lengths <= medium_upper_limit
  large_length_indicies = prague_lengths >= large_lower_limit & prague_lengths <= large_upper_limit
  
  small_length_table =  demographics[small_length_indicies, ]
  medium_length_table = demographics[medium_length_indicies,  ]
  large_length_table = demographics[large_length_indicies, ]
  
  return(list(small_prague = small_length_table, medium_prague = medium_length_table, large_prague = large_length_table))
}



#---------------------PRAGUE CLASS VISUALISATION-------------------------------#
prepare_prague_classification_tables_for_data_visualisation = function (demographics) {
  # We don't want to include samples with NA as their prague length
  demographics = demographics[!is.na(demographics$`(M)`), ]
  prague_lengths = demographics$`(M)`
  
  small_lower_limit = 0
  small_upper_limit = 3
  medium_lower_limit = 4
  medium_upper_limit = 6
  large_lower_limit = 7
  large_upper_limit = 100000000
  
  small_length_indicies = prague_lengths >= small_lower_limit & prague_lengths <= small_upper_limit
  medium_length_indicies = prague_lengths >= medium_lower_limit & prague_lengths <= medium_upper_limit
  large_length_indicies = prague_lengths >= large_lower_limit & prague_lengths <= large_upper_limit
  
  classifications = numeric(length(prague_lengths))
  classifications[small_length_indicies] = 'Small'
  classifications[medium_length_indicies] = 'Medium'
  classifications[large_length_indicies] = 'Large'
  classifications = as.factor(classifications)
  
  demographics = as.data.frame(demographics)
  demographics$prague_class = classifications
  return(demographics)
}

get_number_of_atypia_and_p53_positive = function (prague_classification_data_for_visualisation, p53_overlap=F) {
  indicies = which(prague_classification_data_for_visualisation$Atypia == "Pos")
  if (p53_overlap == T) {
    indicies = which(prague_classification_data_for_visualisation$Atypia == "Pos" & prague_classification_data_for_visualisation$p53 == "Pos")
  }
  reduced_df = prague_classification_data_for_visualisation[indicies, ]
  return(reduced_df)
}

produce_visualisations_for_prague_classification_tables = function (demographics) {
  # We will produce a visualisation of:
  # 1. The number of samples in each prague classification
  # 2. The distrubtion of ages in each prague classification
  # 3. The Numbers of Male and Female in each prague classification
  # 4. The number of TFF3 positive samples in each prague classification
  # 5. The number of P53 positive samples in each prague classification
  # 6. The number of atypia positive samples in each prague classification
  # 7. The number of atypia and p53 positive samples in each prague classification
  
  prague_classification_data_for_visualisation = prepare_prague_classification_tables_for_data_visualisation(demographics)
  # We want the order of the facets to be small, medium, large
  prague_classification_data_for_visualisation$prague_class_order = factor(prague_classification_data_for_visualisation$prague_class, levels = c("Small", "Medium", "Large"))
  
  # For TFF3 and p53, "Yes" is the same as positive and "No" is the same as negative, "" is the same as "Not reported"
  tff3_indicies_to_change = which(prague_classification_data_for_visualisation$TFF3 == "Yes")
  prague_classification_data_for_visualisation$TFF3[tff3_indicies_to_change] = "Pos"
  tff3_indicies_to_change = which(prague_classification_data_for_visualisation$TFF3 == "No")
  prague_classification_data_for_visualisation$TFF3[tff3_indicies_to_change] = "Neg"
  tff3_indicies_to_change = which(prague_classification_data_for_visualisation$TFF3 == "")
  prague_classification_data_for_visualisation$TFF3[tff3_indicies_to_change] = "Not Reported"
  
  p53_indicies_to_change = which(prague_classification_data_for_visualisation$p53 == "Yes")
  prague_classification_data_for_visualisation$p53[p53_indicies_to_change] = "Pos"
  p53_indicies_to_change = which(prague_classification_data_for_visualisation$p53 == "No")
  prague_classification_data_for_visualisation$p53[p53_indicies_to_change] = "Neg"
  p53_indicies_to_change = which(prague_classification_data_for_visualisation$p53 == "")
  prague_classification_data_for_visualisation$p53[p53_indicies_to_change] = "Not Reported"
  
  
  # 1. The number of samples in each prague classification - bar chart
  bar_prague_classification_data_for_visualisation = prague_classification_data_for_visualisation
  print(levels(bar_prague_classification_data_for_visualisation$prague_class_order))
  levels(bar_prague_classification_data_for_visualisation$prague_class_order) = c("Low", "Medium", "High")
  svg("number_of_samples.svg")
  number_of_samples_plot = ggplot(data = bar_prague_classification_data_for_visualisation, aes(x = 1)) +
    geom_bar(fill = "steelblue") +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.3, size = 0)+ #3.5) +
    labs(x = "Prague classification", y = "Number of samples") +
    # Order the x axis by prague classification (Small, Medium, Large)
    scale_x_discrete(limits = c(0, 1, 2)) + facet_wrap(~prague_class_order) + theme(axis.title.x=element_blank(),
                                                                                    axis.text.x=element_blank(),
                                                                                    axis.ticks.x=element_blank())
  number_of_samples_plot = plot_adder(number_of_samples_plot)
  
  print(number_of_samples_plot)
  dev.off()
  
  # 2. The distrubtion of ages in each prague classification - histogram
  svg("distribution_of_ages.svg")
  distribution_of_ages_plot = ggplot(data = prague_classification_data_for_visualisation, aes(x = Age)) +
    geom_histogram(fill = "steelblue") +
    labs(x = "Age", y = "Number of samples") +
    facet_wrap(~prague_class_order)
  
  print(distribution_of_ages_plot)
  dev.off()
  
  # 3. The Numbers of Male and Female in each prague classification - bar chart
  svg("numbers_of_males_and_females.svg")
  numbers_of_males_and_females_plot = ggplot(data = prague_classification_data_for_visualisation, aes(x = 1, fill = Sex)) +
    geom_bar(position = "dodge") +
    geom_text(stat = "count", aes(label = after_stat(count)), position = position_dodge(width = 0.9), vjust = -0.3, size = 0)+ #3.5) +
    labs(x = "Prague classification", y = "Number of samples") +
    scale_fill_discrete(name = "Sex", labels = c("Unclassified", "Female", "Male"), ) + facet_wrap(~prague_class_order) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  print(numbers_of_males_and_females_plot)
  dev.off()
  
  # 4. The number of TFF3 positive samples in each prague classification - bar chart
  svg("number_of_tff3_positive_samples.svg")
  tff3_positive_plot = ggplot(data = prague_classification_data_for_visualisation, aes(x = 1, fill = TFF3)) +
    geom_bar(position = "dodge") +
    geom_text(stat = "count", aes(label = after_stat(count)), position = position_dodge(width = 0.9), vjust = -0.3, size = 0)+ #3.5) +
    labs(x = "Prague classification", y = "Number of samples") + facet_wrap(~prague_class_order) + theme(
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()
    )
  
  print(tff3_positive_plot)
  dev.off()
  
  # 5. The number of P53 positive samples in each prague classification - bar chart
  svg("number_of_p53_positive_samples.svg")
  p53_positive_plot = ggplot(data = prague_classification_data_for_visualisation, aes(x = 1, fill = p53)) +
    geom_bar(position = "dodge") +
    geom_text(stat = "count", aes(label = after_stat(count)), position = position_dodge(width = 0.9), vjust = -0.3, size = 0)+ #3.5) +
    labs(x = "Prague classification", y = "Number of p53 positive samples") + facet_wrap(~prague_class_order) + theme(
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()
    )
  
  print(p53_positive_plot)
  dev.off()
  
  # 6. The number of atypia positive samples in each prague classification
  prague_df_for_atypia_plot = get_number_of_atypia_and_p53_positive(prague_classification_data_for_visualisation)
  svg("number_of_atypia_positive_samples.svg")
  atypia_plot = ggplot(data = prague_df_for_atypia_plot, aes(x = 1)) +
    geom_bar(position = "dodge") +
    geom_text(stat = "count", aes(label = after_stat(count)), position = position_dodge(width = 0.9), vjust = -0.3, size = 0)+ #3.5) +
    labs(x = "Prague classification", y = "Number of atypia positive samples") + facet_wrap(~prague_class_order) + scale_x_discrete(limits = c(0, 1, 2)) + facet_wrap(~prague_class_order) + theme(axis.title.x=element_blank(),
                                                                                                                                                                                   axis.text.x=element_blank(),
                                                                                                                                                                                   axis.ticks.x=element_blank())
  print(atypia_plot)
  dev.off()
  
  # 7. The number of atypia and p53 positive samples in each prague classification
  prague_df_for_atypia_and_p53_plot = get_number_of_atypia_and_p53_positive(prague_classification_data_for_visualisation, T)
  svg("number_of_atypia_and_p53_positive_samples.svg")
  atypia_and_p53_plot = ggplot(data = prague_df_for_atypia_and_p53_plot, aes(x = 1)) +
    geom_bar(position = "dodge") +
    geom_text(stat = "count", aes(label = after_stat(count)), position = position_dodge(width = 0.9), vjust = -0.3, size = 0)+ #3.5) +
    labs(x = "Prague classification", y = "Number of atypia and p53 positive samples") + facet_wrap(~prague_class_order) + scale_x_discrete(limits = c(0, 1, 2)) + facet_wrap(~prague_class_order) + theme(axis.title.x=element_blank(),
                                                                                                                                                                                   axis.text.x=element_blank(),
                                                                                                                                                                                   axis.ticks.x=element_blank()) +
    scale_y_continuous(breaks = integer_breaks())
  print(atypia_and_p53_plot)
  dev.off()
  
  return(prague_classification_data_for_visualisation)
  
  # We could make a UMAP of the results and then cluster by Louvain, and label points based on their classification
}


#-------------SUMMARY TABLES---------------------------------------------------#
calculate_mean_by_factor <- function(data, column_name, factor_name) {
  # Ensure the column names are provided as character strings
  column_name <- as.character(substitute(column_name))
  factor_name <- as.character(substitute(factor_name))
  
  # Calculate mean by factor, ignoring NA values
  result <- data %>%
    group_by(!!sym(factor_name)) %>%
    summarise(mean_value = mean(!!sym(column_name), na.rm = TRUE))
  
  return(as.data.frame(result))
}

create_summary_table_for_visualisations = function (prague_classification_data_for_visualisation, outliers=F) {
  number_of_samples_in_each_prague_classification = table(prague_classification_data_for_visualisation$prague_class_order)
  mean_age_in_each_prague_classification = calculate_mean_by_factor(prague_classification_data_for_visualisation, Age, prague_class_order)
  number_of_males_and_females_in_each_prague_classification = table(
    prague_classification_data_for_visualisation$prague_class_order,
    prague_classification_data_for_visualisation$Sex
  )
  
  number_of_TFF3_positive_samples_in_each_prague_classification= table(
    prague_classification_data_for_visualisation$prague_class_order,
    prague_classification_data_for_visualisation$TFF3
  )
  
  number_of_P53_positive_samples_in_each_prague_classification = table(
    prague_classification_data_for_visualisation$prague_class_order,
    prague_classification_data_for_visualisation$p53
  )
  
  number_of_atypia_positive_samples_in_each_prague_classification = table(
    get_number_of_atypia_and_p53_positive(prague_classification_data_for_visualisation)$prague_class_order
  )
  
  number_of_atypia_and_p53_positive_samples_in_each_prague_classification = table(
    get_number_of_atypia_and_p53_positive(prague_classification_data_for_visualisation, T)$prague_class_order
  )
  
  # Transpose all tables and then join them by row, as they have a different number of columns
  number_of_samples_in_each_prague_classification = t(number_of_samples_in_each_prague_classification)
  mean_age_in_each_prague_classification = t(mean_age_in_each_prague_classification)
  number_of_males_and_females_in_each_prague_classification = t(number_of_males_and_females_in_each_prague_classification)
  number_of_TFF3_positive_samples_in_each_prague_classification = t(number_of_TFF3_positive_samples_in_each_prague_classification)
  number_of_P53_positive_samples_in_each_prague_classification = t(number_of_P53_positive_samples_in_each_prague_classification)
  
  summary_table = rbind(number_of_samples_in_each_prague_classification,
                        mean_age_in_each_prague_classification,
                        number_of_males_and_females_in_each_prague_classification,
                        number_of_atypia_positive_samples_in_each_prague_classification,
                        number_of_atypia_and_p53_positive_samples_in_each_prague_classification,
                        number_of_TFF3_positive_samples_in_each_prague_classification,
                        number_of_P53_positive_samples_in_each_prague_classification
  )
  
  # Drop the second row as it just contains "Small", "Medium", "Large" as an artefact from calculate_mean_by_factor
  summary_table = summary_table[-2,]
  
  # Format the new second row (mean ages) to 4 significant figures
  summary_table[2,] = as.numeric(summary_table[2,])
  summary_table[2,] = format(as.numeric(summary_table[2,]), scientific = FALSE, digits = 4)
  # The rownames should be named ready for figures
  rownames(summary_table)[1] = "Samples"
  rownames(summary_table)[2] = "Mean age"
  rownames(summary_table)[3] = "Unreported genders"
  rownames(summary_table)[4] = "Males"
  rownames(summary_table)[5] = "Females"
  rownames(summary_table)[6] = "Atypia pos"
  rownames(summary_table)[7] = "Atypia and p53 pos"
  tff3_new_rownames = paste("TFF3", tolower(rownames(number_of_TFF3_positive_samples_in_each_prague_classification)))
  tff3_start_index = 8
  rownames(summary_table)[tff3_start_index:(tff3_start_index + length(tff3_new_rownames) - 1)] = tff3_new_rownames
  p53_start_row = nrow(number_of_TFF3_positive_samples_in_each_prague_classification) + tff3_start_index
  p53_new_rownames = paste("p53", tolower(rownames(number_of_P53_positive_samples_in_each_prague_classification)))
  rownames(summary_table)[p53_start_row:(length(p53_new_rownames) + p53_start_row - 1)] = p53_new_rownames
  
  summary_df = data.frame(summary_table)
  
  # We only want the rows that are reported, or are positive
  summary_df = summary_df[-c(3, 8, 9, 10, 12, 13, 14),]

  # Save df as a .csv file
  setwd()
  write.csv(summary_df, file = paste0("summary", outliers, ".csv"))
  
  return(summary_df)
}

#---------CREATE VIS AND SUMMARY------------------------------------------------#
prague_class_df = produce_visualisations_for_prague_classification_tables(case_demo_data)
create_summary_table_for_visualisations(prague_class_df)

#-----------------------------CREATE OUTLIER TABLE-----------------------------#
create_outliers_table = function() {
  outliers = which(concatenated_tables$uniqueID%!in%concatenated_tables_cleaned$uniqueID)
  
  demo_outliers = case_demo_data[outliers,]
  demo_outliers_visualtisation_table = prepare_prague_classification_tables_for_data_visualisation(demo_outliers)
  # We need this factor for create_summary_table_for_visualisations() to work
  demo_outliers_visualtisation_table$prague_class_order = factor(demo_outliers_visualtisation_table$prague_class, levels = c("Small", "Medium", "Large"))
  
  return(create_summary_table_for_visualisations(demo_outliers_visualtisation_table, outliers=T))
}

create_outliers_table()


#------findNumberOfRefluxAndBarrettsInOutliers--------------------------------#
find_numbers_of_barrets_and_reflux_in_outliers = function () {
  outlier_indicies = which(case_demo_data$Case %!in% prague_class_df$Case)
  outlier_data = case_demo_data[outlier_indicies,]
  number_of_reflux_outliers = length(which(outlier_data$Pathway == "Reflux"))
  number_of_barrett_outliers = length(which(outlier_data$Pathway == "Barretts"))
  
  return(c(number_of_reflux_outliers, number_of_barrett_outliers))
}

number_of_relfux_barrets_outliers = find_numbers_of_barrets_and_reflux_in_outliers()

paste0("Number of reflux outliers from demographics cleaning: ", number_of_relfux_barrets_outliers[1])
paste0("Number of barrets outliers from demographics cleaning: ", number_of_relfux_barrets_outliers[2])

















