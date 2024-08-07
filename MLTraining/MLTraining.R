library(ggplot2)
library(dplyr)
library(caret)
library(randomForest)
library(Metrics)
library(smotefamily)
library(umap)
library(pROC)
library(geomtextpath)
library(PRROC)

set.seed(511)

setwd()

# Set this for plots - make labels bigger for paper
plot_adder = function (plot) {
  rel_size = 5.5
  plot + theme(legend.position = "right", text = element_text(size = rel(rel_size)),
               legend.text = element_text(size = rel(rel_size)), plot.title = element_text(size = rel(rel_size))) +
    guides(color = guide_legend(override.aes = list(linewidth = rel(rel_size), size = rel(rel_size))),
           shape = guide_legend(override.aes = list(size = rel(rel_size))))
}

# Read in prague_meth_df to plot - it is just the DMRs that are significant
clean_df = readRDS()
clean_df_orig = readRDS()
dmr_granges = readRDS()

# Make this df so we can pass it to Caret to apply SMOTE on it within its function
clean_df_orig_preprocessed = clean_df_orig[, c(2:ncol(clean_df))]
clean_df_orig_preprocessed$PragueClass = as.integer(clean_df_orig_preprocessed$PragueClass)
clean_df_orig_preprocessed$PragueClass = factor(clean_df_orig_preprocessed$PragueClass)
# clean_df_orig_preprocessed[, 3:ncol(clean_df_orig_preprocessed)] = predict(preProcess(
#   clean_df_orig_preprocessed[, 3:ncol(clean_df_orig_preprocessed)], method = c("center", "scale")),
#                                                                            clean_df_orig_preprocessed[, 3:ncol(clean_df_orig_preprocessed)]
# )

clean_df_fully_pre_processed = clean_df_orig_preprocessed[, c(1, 3:ncol(clean_df_orig_preprocessed))]
# Need to make "X3" the positive class
clean_df_fully_pre_processed$PragueClass = make.names(clean_df_fully_pre_processed$PragueClass)
clean_df_fully_pre_processed$PragueClass = as.factor(clean_df_fully_pre_processed$PragueClass)
clean_df_fully_pre_processed$PragueClass = relevel(clean_df_fully_pre_processed$PragueClass, ref = "X3")
index = createDataPartition(clean_df_fully_pre_processed$PragueClass, p = 0.7, list = FALSE)
training_df = clean_df_fully_pre_processed[index,]
test_df = clean_df_fully_pre_processed[-index,]
training_test_split = 0.8

'%!in%' <- function(x, y)!('%in%'(x, y))

# Maybe use technique to sort out class imbalance:
plot_smote_umap = function(smote_output, k_value = NA, c_value = NA, eps_value = NA) {
  # Combine original and synthetic data
  original_data <- rbind(smote_output$orig_N, smote_output$orig_P)
  synthetic_data <- smote_output$syn_data

  # The column labels do not match between the two data frames, so we must change them before rbind()
  colnames(original_data)[1] <- "(M)"
  colnames(original_data)[ncol(original_data)] <- "PragueClass"

  # The first column also needs to be moved to the last column
  original_data <- original_data[, c(ncol(original_data), 1:(ncol(original_data) - 1))]

  # Create a label column to differentiate original and synthetic data
  original_data$Label <- "Original"
  synthetic_data$Label <- "Synthetic"

  # Combine the data frames
  combined_data <- rbind(original_data, synthetic_data)

  # Perform UMAP
  umap_result <- umap(combined_data[, -c(1, ncol(combined_data))])

  # Add UMAP coordinates to the data frame
  combined_data$UMAP1 <- umap_result$layout[, 1]
  combined_data$UMAP2 <- umap_result$layout[, 2]

  prague_labels <- c("X1" = "Low", "X3" = "High")

  # Plot using ggplot2
  smote_umap = ggplot(combined_data, aes(x = UMAP1, y = UMAP2, color = factor(PragueClass, labels = prague_labels), shape = Label)) +
    #geom_point(size = 3, alpha = 0.7) +
    geom_point(size = 3, fill = NA, stroke = 0.5) +
    scale_shape_manual(values = c(1, 2)) +
    labs(color = "Risk group",
         shape = "Data Type") +
    theme_minimal()

  smote_umap  = plot_adder(smote_umap)

  # Create a title for the plot displaying the values of k and c, just k, or just c, or just epsilon when applicable
  if (is.na(k_value) && is.na(c_value)) {
    return(smote_umap)
  } else if (!is.na(k_value) && is.na(c_value)) {
    return(smote_umap + ggtitle(paste0("k = ", k_value)))
  } else if (is.na(k_value) && !is.na(c_value)) {
    return(smote_umap + ggtitle(paste0("c = ", c_value)))
  } else if (!is.na(eps_value)) {
    return(smote_umap + ggtitle(paste0("epsilon = ", eps_value)))
  } else {
    return(smote_umap + ggtitle(paste0("k = ", k_value, ", c = ", c_value)))
  }

  return(smote_umap)
}

# 50 or more for k is too high - Requesting more near neighbors than data points
smote_on_clean_df = function(without_M = F, k_values = c(5, 10, 20, 40)) {
  # smote_df = clean_df[, c(2, 4:ncol(clean_df))]
  smote_df = clean_df[, c(2:ncol(clean_df))]

  if (without_M) {
    smote_df = clean_df[, c(2, 4:ncol(clean_df))]
  }

  smote_df$PragueClass = as.integer(smote_df$PragueClass)
  smote_df$PragueClass = factor(smote_df$PragueClass)

  for (k in k_values) {
    smote_output = SMOTE(smote_df[, -1], smote_df[, 1], K = k)
    # The algorithm makes the target variable a character
    smote_output$data$class = as.integer(smote_output$data$class)
    smote_output$data$class = factor(smote_output$data$class)
    smote_output$syn_data$class = as.integer(smote_output$syn_data$class)
    smote_output$syn_data$class = factor(smote_output$syn_data$class)
    # Reorder the data frame as the algorithm moves the target variable to the end
    num_cols = ncol(smote_output$data)
    new_order = c(num_cols, 1:(num_cols - 1))

    smote_output$data = smote_output$data[, new_order]
    smote_output$syn_data = smote_output$syn_data[, new_order]
    colnames(smote_output$data) = c("PragueClass", colnames(smote_output$data)[-1])
    colnames(smote_output$syn_data) = c("PragueClass", colnames(smote_output$syn_data)[-1])

    svg(paste0("smote_", k, ".svg"))
    print(plot_smote_umap(smote_output, k_value = k))
    dev.off()
  }


  return(smote_output$data)

  #return(smote(clean_df[, c(2:ncol(clean_df))], "PragueClass"))
}

borderline_smote_on_clean_df = function(without_M = F, k_value = 5, c_value = 5) {
  bsmote_df = clean_df[, c(2:ncol(clean_df))]

  if (without_M) {
    bsmote_df = clean_df[, c(2, 4:ncol(clean_df))]
  }

  bsmote_df$PragueClass = as.integer(bsmote_df$PragueClass)
  bsmote_df$PragueClass = as.integer(bsmote_df$PragueClass)

  # Doesn't work: K cannot exceed the size of DANGER
  bsmote_output = BLSMOTE(bsmote_df[, -1], bsmote_df[, 1], K = k_value, C = c_value)

  return(bsmote_output)
}

db_smote_on_clean_df = function(without_M = F) {
  dbsmote_df = clean_df[, c(2:ncol(clean_df))]

  if (without_M) {
    dbsmote_df = clean_df[, c(2, 4:ncol(clean_df))]
  }

  dbsmote_df$PragueClass = as.integer(dbsmote_df$PragueClass)
  dbsmote_df$PragueClass = factor(dbsmote_df$PragueClass)

  # Doesn't work: K cannot exceed the size of DANGER
  dbsmote_output = DBSMOTE(dbsmote_df[, -1], dbsmote_df[, 1])

  dbsmote_output$data$class = as.integer(dbsmote_output$data$class)
  dbsmote_output$syn_data$class = as.integer(dbsmote_output$syn_data$class)
  dbsmote_output$data$class = factor(dbsmote_output$data$class)
  dbsmote_output$syn_data$class = factor(dbsmote_output$syn_data$class)

  # Reorder the data frame as the algorithm moves the target variable to the end
  num_cols = ncol(dbsmote_output$data)
  new_order = c(num_cols, 1:(num_cols - 1))

  dbsmote_output$data = dbsmote_output$data[, new_order]
  dbsmote_output$syn_data = dbsmote_output$syn_data[, new_order]
  colnames(dbsmote_output$data) = c("PragueClass", colnames(dbsmote_output$data)[-1])
  colnames(dbsmote_output$syn_data) = c("PragueClass", colnames(dbsmote_output$syn_data)[-1])

  svg(paste0("dbsmote.svg"))
  print(plot_smote_umap(dbsmote_output))
  dev.off()


  return(dbsmote_output)
}

safe_level_smote_on_clean_df = function(without_M = F, k_values = c(5, 10, 20, 40), c_values = c(5, 10, 50, 100)) {
  slsmote_df = clean_df[, c(2:ncol(clean_df))]

  if (without_M) {
    slsmote_df = clean_df[, c(2, 4:ncol(clean_df))]
  }

  slsmote_df$PragueClass = as.integer(slsmote_df$PragueClass)
  slsmote_df$PragueClass = factor(slsmote_df$PragueClass)

  for (k in k_values) {
    for (c in c_values) {
      slsmote_output = SLS(slsmote_df[, -1], slsmote_df[, 1], K = k, C = c)

      # The algorithm makes the target variable a character
      slsmote_output$data$class = as.integer(slsmote_output$data$class)
      slsmote_output$data$class = factor(slsmote_output$data$class)
      slsmote_output$syn_data$class = as.integer(slsmote_output$syn_data$class)
      slsmote_output$syn_data$class = factor(slsmote_output$syn_data$class)
      # Reorder the data frame as the algorithm moves the target variable to the end
      num_cols = ncol(slsmote_output$data)
      new_order = c(num_cols, 1:(num_cols - 1))

      slsmote_output$data = slsmote_output$data[, new_order]
      slsmote_output$syn_data = slsmote_output$syn_data[, new_order]
      colnames(slsmote_output$data) = c("PragueClass", colnames(slsmote_output$data)[-1])
      colnames(slsmote_output$syn_data) = c("PragueClass", colnames(slsmote_output$syn_data)[-1])

      svg(paste0("slsmote_", k, "_", c, ".svg"))
      print(plot_smote_umap(slsmote_output, k_value = k, c_value = c))
      dev.off()
    }
  }
  # slsmote_output = SLS(slsmote_df[, -1], slsmote_df[, 1], K=20)
  #
  # return(slsmote_output)
}

make_roc_curve = function(model) {
  #Substitute in or out the rev depending on what class you want to be interpreted as +ve - first given to levels param
  # is interpreted as the +ve one - do I want to optimise for large or small class
  roc_curve = roc(response = model$pred$obs, predictor = as.numeric(model$pred$pred),
                  levels = rev(levels(model$pred$obs)))
  auc_value = auc(roc_curve)
  auc_roc_curve = plot(roc_curve, col = "blue", lwd = 2) +
    text(x = 0.1, y = 0.1, labels = sprintf("AUC = %.2f", auc_value), col = "red")


  # Make roc for test as well


  print(auc_roc_curve)
}

make_roc_curve_for_test = function(model, test_data) {
  predictions = predict(model, newdata = test_data, type = "prob")

  # Compute the ROC curve and AUC
  roc_curve <- roc(response = test_data$PragueClass, predictor = predictions$`3`) #,
  # Uncomment this if you want to change it to the +ve class - will need to pass train data as an arg
  #levels = rev(levels(train_data$PragueClass)))
  auc_value <- auc(roc_curve)

  # Plot the ROC curve
  roc_plot <- ggroc(roc_curve, aes = c("colour" = "blue"), lwd = 2) +
    theme_minimal() +
    annotate("text", x = 0.1, y = 0.1, label = sprintf("AUC = %.2f", auc_value), color = "red", size = 5)

  print(roc_plot)

  return(roc_plot)
}

# Then we want to go through each region and fit a linear model to it
fit_linear_model_per_region = function() {
  plot_list = list()
  # The first three columns are UniqueID, PragueClass, (M)
  for (i in 4:ncol(clean_df)) {
    plot_list[[i - 3]] = ggplot(clean_df, aes(x = `(M)`, y = clean_df[, i], color = factor(PragueClass))) +
      geom_point() +
      ylab(paste("Methylation level (%) in", dmr_granges[i - 3])) +
      xlab("Prague (M) Length (cm)") +
      theme_bw() +
      geom_smooth(method = "lm") +
      scale_fill_discrete(name = "Prague Classification")
  }

  for (i in 1:(ncol(clean_df) - 3)) {
    print(plot_list[[i]])
    if (i == 1) {
      #svg(paste0("methylationDistributionPerRegion_", i, ".svg"))
      #print(plot_list[[i]])
      #dev.off()
    }

    if (i == 2) {
      #svg(paste0("methylationDistributionPerRegion_", i, ".svg"))
      #print(plot_list[[i]])
      #dev.off()
    }
  }
}

# fit_linear_model_per_region()
clean_df = smote_on_clean_df()

# For SMOTE and safe level SMOTE, k was set to 5, 10, 20 and 40 and in safe level SMOTE c was set to 5, 10, 50 and 100
# Therefore, we need a separate list for each combination of parameters

# Keeping this one to make original functions work
smote_list_for_caret = list(
  name = "Smote from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = SMOTE(x, y)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)

smote_list_for_caret_5 = list(
  name = "Smote from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = SMOTE(x, y)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)

smote_list_for_caret_10 = list(
  name = "Smote from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = SMOTE(x, y)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)

smote_list_for_caret_20 = list(
  name = "Smote from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = SMOTE(x, y)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)

smote_list_for_caret_40 = list(
  name = "Smote from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = SMOTE(x, y)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)


db_smote_list_for_caret = list(
  name = "DBSmote from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = DBSMOTE(as.data.frame(x), y)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    # We need to return a list with elements of the same name
    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)

# d = training_df
# index1 = createDataPartition(d$PragueClass, p = 0.8, list = FALSE)
# print(levels(d$PragueClass))
# levels(d$PragueClass) = make.names(levels(d$PragueClass))
# print(levels(d$PragueClass))
# d = d[index1, ]
# n = d[, 2:ncol(d)]
# m = d[, 1]
#
# gr = function(x, y) {
#     # This function must take arguments named x and y
#     library(smotefamily)
#
#     print(levels(y))
#     x = as.data.frame(x)
#     rownames(x) = names(y)
#
#     dat = SLS(as.data.frame(x), y, K = 5, C = 5)$data
#     # dat$class = as.integer(dat$class)
#     dat$class = factor(dat$class)
#     # We need to relevel so that X3 is the positive class
#     dat$class = relevel(dat$class, ref = "X3")
#
#     # We need to return a list with elements of the same name
#     list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
#          y = dat$class)
#   }
#
# e = gr(n, m)

# Now do the same for SLS, varying k and c as in the list below
safe_level_smote_list_for_caret_k5_c5 = list(
  name = "Safe level SMOTE from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = SLS(as.data.frame(x), y, K = 5, C = 5)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    # We need to return a list with elements of the same name
    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)

safe_level_smote_list_for_caret_k5_c10 = list(
  name = "Safe level SMOTE from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = SLS(as.data.frame(x), y, K = 5, C = 10)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    # We need to return a list with elements of the same name
    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)

safe_level_smote_list_for_caret_k5_c50 = list(
  name = "Safe level SMOTE from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = SLS(as.data.frame(x), y, K = 5, C = 50)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    # We need to return a list with elements of the same name
    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)

safe_level_smote_list_for_caret_k5_c100 = list(
  name = "Safe level SMOTE from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = SLS(as.data.frame(x), y, K = 5, C = 100)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    # We need to return a list with elements of the same name
    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)

safe_level_smote_list_for_caret_k10_c5 = list(
  name = "Safe level SMOTE from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = SLS(as.data.frame(x), y, K = 10, C = 5)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    # We need to return a list with elements of the same name
    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)

safe_level_smote_list_for_caret_k10_c10 = list(
  name = "Safe level SMOTE from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = SLS(as.data.frame(x), y, K = 10, C = 10)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    # We need to return a list with elements of the same name
    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)

safe_level_smote_list_for_caret_k10_c50 = list(
  name = "Safe level SMOTE from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = SLS(as.data.frame(x), y, K = 10, C = 50)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    # We need to return a list with elements of the same name
    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)

safe_level_smote_list_for_caret_k10_c100 = list(
  name = "Safe level SMOTE from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = SLS(as.data.frame(x), y, K = 10, C = 100)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    # We need to return a list with elements of the same name
    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)

safe_level_smote_list_for_caret_k20_c5 = list(
  name = "Safe level SMOTE from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = SLS(as.data.frame(x), y, K = 20, C = 5)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    # We need to return a list with elements of the same name
    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)

safe_level_smote_list_for_caret_k20_c10 = list(
  name = "Safe level SMOTE from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = SLS(as.data.frame(x), y, K = 20, C = 10)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    # We need to return a list with elements of the same name
    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)

safe_level_smote_list_for_caret_k20_c50 = list(
  name = "Safe level SMOTE from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = SLS(as.data.frame(x), y, K = 20, C = 50)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    # We need to return a list with elements of the same name
    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)

safe_level_smote_list_for_caret_k20_c100 = list(
  name = "Safe level SMOTE from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = SLS(as.data.frame(x), y, K = 20, C = 100)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    # We need to return a list with elements of the same name
    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)

safe_level_smote_list_for_caret_k25_c5 = list(
  name = "Safe level SMOTE from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = SLS(as.data.frame(x), y, K = 25, C = 5)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    # We need to return a list with elements of the same name
    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)

safe_level_smote_list_for_caret_k25_c10 = list(
  name = "Safe level SMOTE from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = SLS(as.data.frame(x), y, K = 25, C = 10)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    # We need to return a list with elements of the same name
    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)

safe_level_smote_list_for_caret_k25_c50 = list(
  name = "Safe level SMOTE from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = SLS(as.data.frame(x), y, K = 25, C = 50)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    # We need to return a list with elements of the same name
    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)

safe_level_smote_list_for_caret_k25_c100 = list(
  name = "Safe level SMOTE from smotefamily",
  func = function(x, y) {
    # This function must take arguments named x and y
    library(smotefamily)

    x = as.data.frame(x)
    # rownames(x) = names(y)

    dat = SLS(as.data.frame(x), y, K = 25, C = 100)$data
    dat$class = factor(dat$class)
    dat$class = relevel(dat$class, ref = "X3")

    # We need to return a list with elements of the same name
    list(x = dat[, !grepl("class", colnames(dat), fixed = TRUE)],
         y = dat$class)
  },
  first = TRUE
)

total_smote_list = list(
  # SMOTE #
  smote_list_for_caret_5,
  smote_list_for_caret_10,
  smote_list_for_caret_20,
  smote_list_for_caret_40,
  # DBSMOTE #
  db_smote_list_for_caret,
  # SLS #
  safe_level_smote_list_for_caret_k5_c5,
  safe_level_smote_list_for_caret_k5_c10,
  safe_level_smote_list_for_caret_k5_c50,
  safe_level_smote_list_for_caret_k5_c100,
  safe_level_smote_list_for_caret_k10_c5,
  safe_level_smote_list_for_caret_k10_c10,
  safe_level_smote_list_for_caret_k10_c50,
  safe_level_smote_list_for_caret_k10_c100,
  safe_level_smote_list_for_caret_k20_c5,
  safe_level_smote_list_for_caret_k20_c10,
  safe_level_smote_list_for_caret_k20_c50,
  safe_level_smote_list_for_caret_k20_c100,
  safe_level_smote_list_for_caret_k25_c5,
  safe_level_smote_list_for_caret_k25_c10,
  safe_level_smote_list_for_caret_k25_c50,
  safe_level_smote_list_for_caret_k25_c100
)

names(total_smote_list) = c(paste0("SMOTE(", c(5, 10, 20, 40), ") "), "DBSMOTE", paste0("SLS(", rep(c(5, 10, 20, 25), each = 4), ", ", c(5, 10, 50, 100), ")"))
#names(total_smote_list) = c(paste0("SLS(", rep(c(5), each = 4), ", ", c(5, 10, 50, 100), ")"))

make_roc_curve_multiple_models = function(models, test_data, top_three = T) {
  number_of_models = length(models)
  roc_curve_results = list()
  auc_values = numeric(number_of_models)
  levels(test_data$PragueClass) = make.names(levels(test_data$PragueClass))

  for (i in 1:number_of_models) {
    model = models[[i]]
    predictions = predict(model, newdata = test_data, type = "prob")
    roc_curve = roc(response = test_data$PragueClass, predictor = predictions$`X3`) #,
    # Uncomment this if you want to change it to the +ve class - will need to pass train data as an arg
    #levels = rev(levels(train_data$PragueClass)))

    auc_value = auc(roc_curve)
    roc_coords = coords(roc_curve)
    # roc_curve$model = names(models)[i]
    # roc_curve_results[[i]] = roc_curve

    roc_curve_results[[i]] = data.frame(
      specificity = roc_curve$specificities,
      sensitivity = roc_curve$sensitivities,
      model = names(models)[i]
    )

    auc_values[i] = auc_value
  }

  if (!top_three) {
    roc_curves = do.call(rbind, roc_curve_results)

    # Plot the ROC curve
    roc_auc_plot = ggplot(roc_curves, aes(1 - specificity, sensitivity, color = factor(model))) +
      geom_line(size = 1) +
      geom_abline(linetype = "dashed") +
      scale_x_continuous("1-Specificity") +
      scale_y_continuous("Sensitvity") +
      labs(color = "Algorithm") +
      theme(legend.position = "right")

    return(list(plots = roc_auc_plot,
                aucs = auc_values,
                models = models))
  }
  # Select the three best auc values and make a roc plot of them
  top_three_indicies = order(auc_values, decreasing = T)[1:3]

  top_three_aucs = auc_values[top_three_indicies]
  top_three_roc_curves = do.call(rbind, roc_curve_results[top_three_indicies])

  # Plot the ROC curve
  top_three_roc_auc_plot = ggplot(top_three_roc_curves, aes(1 - specificity, sensitivity, color = factor(model))) +
    geom_line(size = 1) +
    geom_abline(linetype = "dashed") + #, slope = 1, intercept = 0, ) +
    # scale_color_brewer(palette = "Set1", guide = "none", direction = -1) +
    scale_x_continuous("1-Specificity") + #, labels = scales::percent) +
    scale_y_continuous("Sensitivity") + #, labels = scales::percent) +
    #coord_equal(expand = FALSE) +
    labs(title = paste0(names(models)[top_three_indicies][1], " - AUC = ",
                        signif(top_three_aucs[1], 3), ", ",
                        names(models)[top_three_indicies][2], " - AUC = ",
                        signif(top_three_aucs[2], 3), ", ",
                        names(models)[top_three_indicies][3], " - AUC = ",
                        signif(top_three_aucs[3], 3)),
         color = "Model") +
    theme(legend.position = "right")

  return(list(plots = top_three_roc_auc_plot,
              aucs = top_three_aucs,
              models = models[top_three_indicies]))
}

make_training_test_roc_curve_one_model = function(model, test_data) {
  roc_curve_results = list()
  auc_values = numeric(2)
  # levels(test_data$PragueClass) = make.names(levels(test_data$PragueClass))

  # Get the training ROC curve
  train_predictions = model$pred$`X3`
  roc_curve_train = roc(response = model$pred$obs, predictor = train_predictions)

  train_auc_value = auc(roc_curve_train)

  roc_curve_results[[1]] = data.frame(
    specificity = roc_curve_train$specificities,
    sensitivity = roc_curve_train$sensitivities,
    model = "Training"
  )

  auc_values[1] = train_auc_value

  # Get the test ROC curve
  test_predictions = predict(model, newdata = test_data, type = "prob")
  roc_curve_test = roc(response = test_data$PragueClass, predictor = test_predictions$`X3`)

  test_auc_value = auc(roc_curve_test)
  # roc_coords = coords(roc_curve_test)

  roc_curve_results[[2]] = data.frame(
    specificity = roc_curve_test$specificities,
    sensitivity = roc_curve_test$sensitivities,
    model = "Test"
  )

  auc_values[2] = test_auc_value


  roc_curves = do.call(rbind, roc_curve_results)

  # Plot the ROC curve
  roc_auc_plot = ggplot(roc_curves, aes(1 - specificity, sensitivity, color = factor(model))) +
    geom_line(size = 1) +
    geom_abline(linetype = "dashed") +
    scale_x_continuous("False Positive Rate") +
    scale_y_continuous("True Positive Rate") +
    labs(color = "Model") +
    theme(legend.position = "right")

  return(roc_auc_plot)
}


make_precision_recall_curve_multiple_models = function(models, test_data) {
  number_of_models = length(models)
  precision_recall_curve_results = list()
  auc_results = numeric(number_of_models)
  for (i in 1:number_of_models) {
    model = models[[i]]
    predictions = predict(model, newdata = test_data, type = "prob")
    pr_curve = pr.curve(scores.class0 = predictions$`X3`, scores.class1 = predictions$`X1`, curve = T)
    auc_results[i] = pr_curve$auc.integral

    pr_curve = as.data.frame(pr_curve$curve)
    pr_curve$model = names(models)[i]
    precision_recall_curve_results[[i]] = pr_curve
  }

  top_three_indicies = order(auc_results, decreasing = T)[1:3]
  top_three_aucs = auc_results[top_three_indicies]
  top_three_pr_curves = do.call(rbind, precision_recall_curve_results[top_three_indicies])

  # Plot the ROC curve
  top_three_pr_auc_plot = ggplot(top_three_pr_curves, aes(V1, V2, color = factor(model))) +
    geom_line(size = 1) +
    scale_x_continuous("Recall") +
    scale_y_continuous("Precision") +
    labs(title = paste0(names(models)[top_three_indicies][1], " - AUC = ",
                        signif(top_three_aucs[1], 3), ", ",
                        names(models)[top_three_indicies][2], " - AUC = ",
                        signif(top_three_aucs[2], 3), ", ",
                        names(models)[top_three_indicies][3], " - AUC = ",
                        signif(top_three_aucs[3], 3)),
         color = "Model") +
    theme(legend.position = "right")


  return(list(plots = top_three_pr_auc_plot,
              aucs = top_three_aucs,
              models = models[top_three_indicies]))
}


# Then we want to look into random forests
random_forest_classification = function(smote_in_caret = TRUE) {
  x_train = training_df[, 2:ncol(training_df)]
  y_train = training_df[, 1]

  metric <- 'accuracy'
  customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
  customRF$parameters <- data.frame(parameter = c("maxnodes", "ntree"), class = rep("numeric", 2), label = c("maxnodes", "ntree"))
  customRF$grid <- function(x, y, len = NULL, search = "grid") { }

  customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    randomForest(x, y, maxnodes = param$maxnodes, ntree = param$ntree, ...)
  }

  customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
    predict(modelFit, newdata)

  customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
    predict(modelFit, newdata, type = "prob")

  customRF$sort <- function(x) x[order(x[, 1]),]
  customRF$levels <- function(x) x$classes
  # Set grid search parameters
  control = NA
  if (smote_in_caret) {
    control <- trainControl(method = "repeatedcv", number = 10, repeats = 3, sampling = smote_list_for_caret,
                            p = training_test_split, savePredictions = T)
  } else {
    control <- trainControl(method = "repeatedcv", number = 10, repeats = 3,
                            p = training_test_split, savePredictions = T)
  }

  # Outline the grid of parameters
  # The accuracy improves around 50 nodes so sample more around it 
  tunegrid <- expand.grid(.maxnodes = c(5, 30, 40, 50, 75, 85, 90, 100, 500), .ntree = c(10, 100, 500, 1000))


  rf_gridsearch <- train(x = x_train, y = y_train, method = customRF, metric = metric, tuneGrid = tunegrid,
                         trControl = control)

  random_forest_error_plots = plot(rf_gridsearch)
  print(random_forest_error_plots)

  feature_importance_plot = varImpPlot(rf_gridsearch$finalModel, main = 'Feature importance')
  print(feature_importance_plot)
  roc_curve = make_roc_curve_for_test(rf_gridsearch, test_df)
  return(rf_gridsearch)
}

#random_forest_model = random_forest_classification()

random_forest_classification_all = function() {
  # levels(training_df$PragueClass) = make.names(levels(training_df$PragueClass))
  x_train = training_df[, 2:ncol(training_df)]
  y_train = training_df[, 1]

  metric <- 'ROC'
  customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
  customRF$parameters <- data.frame(parameter = c("maxnodes", "ntree"), class = rep("numeric", 2), label = c("maxnodes", "ntree"))
  customRF$grid <- function(x, y, len = NULL, search = "grid") { }

  customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    randomForest(x, y, maxnodes = param$maxnodes, ntree = param$ntree, ...)
  }

  customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
    predict(modelFit, newdata)

  customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
    predict(modelFit, newdata, type = "prob")

  customRF$sort <- function(x) x[order(x[, 1]),]
  customRF$levels <- function(x) x$classes
  # Set grid search parameters
  # The accuracy improves around 50 nodes so sample more around it
  tunegrid <- expand.grid(.maxnodes = c(5, 30, 40, 50, 75, 85, 90, 100, 500), .ntree = c(10, 100, 500, 1000))

  number_of_smoted_datasets = length(total_smote_list)
  models = list()
  for (i in 1:number_of_smoted_datasets) {
    # models[[i]] = train(x = x_train, y = y_train, method = customRF, metric = metric, tuneGrid = tunegrid,
    #                     trControl = trainControl(method = "repeatedcv", number = 10, repeats = 3, sampling = total_smote_list[[i]],
    #                                              p = training_test_split, savePredictions = T))

    models[[i]] = train(x = x_train, y = y_train, method = customRF, metric = metric, tuneGrid = tunegrid,
                        trControl = trainControl(method = "repeatedcv", number = 5, repeats = 3, p = training_test_split,
                                                 savePredictions = T, sampling = total_smote_list[[i]],
                                                 classProbs = T, summaryFunction = twoClassSummary))
  }

  names(models) = names(total_smote_list)
  top_three_models = make_roc_curve_multiple_models(models, test_df)

  svg(paste0("random_forest_classification_all.svg"))
  print(top_three_models$plots)
  dev.off()
  print(top_three_models$models[[1]]$finalModel)

  top_three_models_pr = make_precision_recall_curve_multiple_models(models, test_df)

  svg("rf_all_pr.svg")
  print(top_three_models_pr$plots)
  dev.off()

  return(list(roc = top_three_models, pr = top_three_models_pr))
}

rf_top_three = random_forest_classification_all()

# K means clustering
k_means_clustering = function() {
  model = train(PragueClass ~ ., data = training_df, method = "knn",
                tuneGrid = expand.grid(k = 2),
                trControl = trainControl(method = "repeatedcv", number = 5, repeats = 3, p = training_test_split,
                                         savePredictions = T, sampling = smote_list_for_caret))

  model_db_smote = train(PragueClass ~ ., data = training_df, method = "knn",
                         tuneGrid = expand.grid(k = 2),
                         trControl = trainControl(method = "repeatedcv", number = 5, repeats = 3, p = training_test_split,
                                                  savePredictions = T, sampling = db_smote_list_for_caret))

  roc_curve = make_roc_curve_for_test(model, test_df)
  roc_curve_db = make_roc_curve_for_test(model_db_smote, test_df)

  return(model)
}

#k_means_model = k_means_clustering()

k_means_all = function() {
  number_of_smoted_datasets = length(total_smote_list)
  models = list()
  levels(training_df$PragueClass) = make.names(levels(training_df$PragueClass))
  for (i in 1:number_of_smoted_datasets) {
    model = train(PragueClass ~ ., data = training_df, method = "knn",
                  tuneGrid = expand.grid(k = 2), preProc = c("center", "scale"),
                  trControl = trainControl(method = "repeatedcv", number = 5, repeats = 3, p = training_test_split,
                                           savePredictions = T, sampling = total_smote_list[[i]],
                                           classProbs = T, summaryFunction = twoClassSummary), metric = "ROC")

    models[[i]] = model
  }
  names(models) = names(total_smote_list)
  top_three_models = make_roc_curve_multiple_models(models, test_df)

  svg("k_means_all.svg")
  print(top_three_models$plots)
  dev.off()

  print(top_three_models$models[[1]]$finalModel)

  top_three_models_pr = make_precision_recall_curve_multiple_models(models, test_df)

  svg("k_means_all_pr.svg")
  print(top_three_models_pr$plots)
  dev.off()

  return(list(roc = top_three_models, pr = top_three_models_pr))
}

km_top_three = k_means_all()

self_organising_map = function() {
  model = train(PragueClass ~ ., data = training_df, method = "xyf",
                trControl = trainControl(method = "repeatedcv", number = 5, repeats = 3, p = training_test_split,
                                         savePredictions = T, sampling = smote_list_for_caret))
  roc_curve = make_roc_curve_for_test(model, test_df)
  return(model)
}

#self_organising_map_model = self_organising_map()

self_organising_map_all = function() {
  number_of_smoted_datasets = length(total_smote_list)
  models = list()
  for (i in 1:number_of_smoted_datasets) {
    print(names(total_smote_list)[i])
    model = train(PragueClass ~ ., data = training_df, method = "xyf",
                  trControl = trainControl(method = "repeatedcv", number = 5, repeats = 3, p = training_test_split,
                                           savePredictions = T, sampling = total_smote_list[[i]]))
    models[[i]] = model
  }
  names(models) = names(total_smote_list)
  top_three_models = make_roc_curve_multiple_models(models, test_df)

  svg("self_organising_map_all.svg")
  print(top_three_models$plots)
  dev.off()

  print(top_three_models$models[[1]]$finalModel)
  return(top_three_models)
}

#som_top_three = self_organising_map_all()

logisitc_regression = function() {
  model = train(PragueClass ~ ., data = training_df, method = "glm",
                trControl = trainControl(method = "repeatedcv", number = 5, repeats = 3, p = training_test_split,
                                         savePredictions = T, sampling = smote_list_for_caret))
  roc_curve = make_roc_curve_for_test(model, test_df)
  return(model)
}

logistic_regression_all = function() {
  number_of_smoted_datasets = length(total_smote_list)
  models = list()
  levels(training_df$PragueClass) = make.names(levels(training_df$PragueClass))
  for (i in 1:number_of_smoted_datasets) {
    model = train(PragueClass ~ ., data = training_df, method = "glm", preProc = c("center", "scale"),
                  trControl = trainControl(method = "repeatedcv", number = 5, repeats = 3, p = training_test_split,
                                           savePredictions = T, sampling = total_smote_list[[i]],
                                           classProbs = T, summaryFunction = twoClassSummary), metric = "ROC")
    models[[i]] = model
  }
  names(models) = names(total_smote_list)
  top_three_models = make_roc_curve_multiple_models(models, test_df)

  svg("logistic_regression_all.svg")
  print(top_three_models$plots)
  dev.off()

  print(top_three_models$models[[1]]$finalModel)

  top_three_models_pr = make_precision_recall_curve_multiple_models(models, test_df)

  svg("lr_all_pr.svg")
  print(top_three_models_pr$plots)
  dev.off()

  return(list(roc = top_three_models, pr = top_three_models_pr))
}

#logistic_regression_model = logisitc_regression()
lr_top_three = logistic_regression_all()

svm_clustering = function(return_all_models = F, smote_in_caret = T) {
  # Do SVM
  model_linear = train(PragueClass ~ ., data = training_df, method = "svmLinear",
                       tuneGrid = expand.grid(C = c(0.001, 0.01, 0.1, 1, 10, 100, 1000)),
                       trControl = trainControl(method = "repeatedcv", number = 5, repeats = 3, p = training_test_split,
                                                savePredictions = T, sampling = smote_list_for_caret))
  # Non-linear SVMs choose the optimal parameter for you
  model_radial = train(PragueClass ~ ., data = training_df, method = "svmRadial",
                       trControl = trainControl(method = "repeatedcv", number = 5, repeats = 3, p = training_test_split,
                                                savePredictions = T, sampling = smote_list_for_caret))
  model_poly = train(PragueClass ~ ., data = training_df, method = "svmPoly",
                     trControl = trainControl(method = "repeatedcv", number = 5, repeats = 3, p = training_test_split,
                                              savePredictions = T, sampling = smote_list_for_caret))

  if (return_all_models) {
    # make_roc_curve(model_linear)
    # make_roc_curve(model_radial)
    # make_roc_curve(model_poly)
    return(list(linear_svm = model_linear, radial_svm = model_radial, poly_svm = model_poly))
  }

  # Now determine the best model
  model = NA
  if (max(model_linear$results$Accuracy) > max(model_radial$results$Accuracy) & max(model_linear$results$Accuracy) > max(model_poly$results$Accuracy)) {
    model = model_linear
  } else if (max(model_radial$results$Accuracy) > max(model_poly$results$Accuracy)) {
    model = model_radial
  } else {
    model = model_poly
  }

  #roc_curve = make_roc_curve_for_test(model, test_df)
  return(model)
}

# svm_model = svm_clustering(T)


svm_clustering_all = function() {
  number_of_smoted_datasets = length(total_smote_list)
  models = list()
  # levels(training_df$PragueClass) = make.names(levels(training_df$PragueClass))
  for (i in 1:number_of_smoted_datasets) {
    model_linear = train(PragueClass ~ ., data = training_df, method = "svmLinear", preProc = c("center", "scale"),
                         trControl = trainControl(method = "repeatedcv", number = 5, repeats = 3, p = training_test_split,
                                                  savePredictions = T, sampling = total_smote_list[[i]],
                                                  classProbs = T, summaryFunction = twoClassSummary), metric = "ROC")
    # Non-linear SVMs choose the optimal parameter for you
    model_radial = train(PragueClass ~ ., data = training_df, method = "svmRadial", preProc = c("center", "scale"),
                         trControl = trainControl(method = "repeatedcv", number = 5, repeats = 3, p = training_test_split,
                                                  savePredictions = T, sampling = total_smote_list[[i]],
                                                  classProbs = T, summaryFunction = twoClassSummary), metric = "ROC")
    model_poly = train(PragueClass ~ ., data = training_df, method = "svmPoly", preProc = c("center", "scale"),
                       trControl = trainControl(method = "repeatedcv", number = 5, repeats = 3, p = training_test_split,
                                                savePredictions = T, sampling = total_smote_list[[i]],
                                                classProbs = T, summaryFunction = twoClassSummary), metric = "ROC")

    model = NA
    if (max(model_linear$results$ROC) > max(model_radial$results$ROC) & max(model_linear$results$ROC) > max(model_poly$results$ROC)) {
      model = model_linear
    } else if (max(model_radial$results$ROC) > max(model_poly$results$ROC)) {
      model = model_radial
    } else {
      model = model_poly
    }

    models[[i]] = model
  }
  names(models) = names(total_smote_list)
  return(models)
  # top_three_models = make_roc_curve_multiple_models(models, test_df)
  # print(top_three_models$plots)
  # print(top_three_models$models[[1]]$finalModel)
  # return(top_three_models)
}

pca_neural_network_clustering = function() {
  model = train(PragueClass ~ ., data = training_df, method = "pcaNNet",
                tuneGrid = expand.grid(size = c(1, 2, 3), decay = c(0.1, 0.01, 0.001)),
                trControl = trainControl(method = "repeatedcv", number = 5, repeats = 3, p = training_test_split,
                                         savePredictions = T, sampling = smote_list_for_caret))
  roc_curve = make_roc_curve_for_test(model, test_df)
  return(model)
}

#pca_nnet_model = pca_neural_network_clustering()

pca_neural_network_clustering_all = function() {
  number_of_smoted_datasets = length(total_smote_list)
  models = list()
  for (i in 1:number_of_smoted_datasets) {
    model = train(PragueClass ~ ., data = training_df, method = "pcaNNet",
                  tuneGrid = expand.grid(size = c(1, 2, 3), decay = c(0.1, 0.01, 0.001)),
                  trControl = trainControl(method = "repeatedcv", number = 5, repeats = 3, p = training_test_split,
                                           savePredictions = T, sampling = total_smote_list[[i]]))
    models[[i]] = model
  }
  names(models) = names(total_smote_list)
  top_three_models = make_roc_curve_multiple_models(models, test_df)

  svg("pca_neural_network_clustering_all.svg")
  print(top_three_models$plots)
  dev.off()

  print(top_three_models$models[[1]]$finalModel)
  return(top_three_models)
}

# pca_top_three = pca_neural_network_clustering_all()

mlp_neural_network_clustering = function() {
  model = train(PragueClass ~ ., data = training_df, method = "mlpML",
                trControl = trainControl(method = "repeatedcv", number = 5, repeats = 3, p = training_test_split,
                                         savePredictions = T, sampling = smote_list_for_caret))
  roc_curve = make_roc_curve_for_test(model, test_df)
  return(model)
}

#mlp_nnet_model = mlp_neural_network_clustering()

mlp_neural_network_clustering_all = function() {
  number_of_smoted_datasets = length(total_smote_list)
  models = list()
  levels(training_df$PragueClass) = make.names(levels(training_df$PragueClass))
  for (i in 1:number_of_smoted_datasets) {
    model = train(PragueClass ~ ., data = training_df, method = "mlpML", preProc = c("center", "scale"),
                  trControl = trainControl(method = "repeatedcv", number = 5, repeats = 3, p = training_test_split,
                                           savePredictions = T, sampling = total_smote_list[[i]],
                                           classProbs = T, summaryFunction = twoClassSummary), metric = "ROC")
    models[[i]] = model
  }
  names(models) = names(total_smote_list)
  top_three_models = make_roc_curve_multiple_models(models, test_df)

  svg("mlp_neural_network_clustering_all.svg")
  print(top_three_models$plots)
  dev.off()

  print(top_three_models$models[[1]]$finalModel)

  top_three_models_pr = make_precision_recall_curve_multiple_models(models, test_df)

  svg("mlp_all_pr.svg")
  print(top_three_models_pr$plots)
  dev.off()

  return(list(roc = top_three_models, pr = top_three_models_pr))
}

mlp_top_three = mlp_neural_network_clustering_all()


xgb_all = function(bootstrap_method = F) {
  number_of_smoted_datasets = length(total_smote_list)
  models = list()
  levels(training_df$PragueClass) = make.names(levels(training_df$PragueClass))
  rownames(training_df) = NULL
  for (i in 1:number_of_smoted_datasets) {
    model = NA
    if (bootstrap_method) {
      model = train(PragueClass ~ ., data = training_df, method = "xgbTree",
                    trControl = trainControl(p = training_test_split,
                                             savePredictions = T, sampling = total_smote_list[[i]],
                                             classProbs = T, summaryFunction = twoClassSummary), metric = "ROC")
    } else {
      model = train(PragueClass ~ ., data = training_df, method = "xgbTree",
                    trControl = trainControl(method = "repeatedcv", number = 5, repeats = 3, p = training_test_split,
                                             savePredictions = T, sampling = total_smote_list[[i]],
                                             classProbs = T, summaryFunction = twoClassSummary), metric = "ROC")
    }
    models[[i]] = model
  }
  names(models) = names(total_smote_list)
  top_three_models = make_roc_curve_multiple_models(models, test_df)

  if (bootstrap_method) {
    svg("xgb_all_bootstrap.svg")
    print(top_three_models$plots)
    dev.off()
  } else {
    svg("xgb_all.svg")
    print(top_three_models$plots)
    dev.off()
  }

  print(top_three_models$models[[1]]$finalModel)


  top_three_models_pr = make_precision_recall_curve_multiple_models(models, test_df)

  if (bootstrap_method) {
    svg("xgb_all_pr_bootstrap.svg")
    print(top_three_models_pr$plots)
    dev.off()
  } else {
    svg("xgb_all_pr.svg")
    print(top_three_models_pr$plots)
    dev.off()
  }

  return(list(roc = top_three_models, pr = top_three_models_pr))
}

xgb_top_three = xgb_all()
xgp_top_three_boot = xgb_all(T)

make_roc_curve_for_best_models = function(test_df) {

  rf_best_model = rf_top_three$roc$models[[1]]
  knn_best_model = km_top_three$roc$models[[1]]
  lr_best_model = lr_top_three$roc$models[[1]]
  mlp_best_model = mlp_top_three$roc$models[[1]]
  xgb_best_model = xgb_top_three$roc$models[[1]]

  models = list(rf_best_model, knn_best_model, lr_best_model, mlp_best_model, xgb_best_model)
  names(models) = c("RF", "KM", "LR", "MLP", "XGB")

  roc_curve = make_roc_curve_multiple_models(models, test_df, F)
  roc_curve$plots = plot_adder(roc_curve$plots)

  svg("all_roc.svg")
  print(roc_curve$plots)
  dev.off()

  auc_df = data.frame(Model = names(models), AUC = c(
    signif(rf_top_three$roc$auc[1], 3),
    signif(km_top_three$roc$auc[1], 3),
    signif(lr_top_three$roc$auc[1], 3),
    signif(mlp_top_three$roc$auc[1], 3),
    signif(xgb_top_three$roc$auc[1], 3)
  ))

  write.csv(auc_df, "allRocAuc.csv", row.names = F)

  return(roc_curve)
}

make_roc_curve_with_decision_threshold = function() {
  pr_xg = predict(xgb_top_three$roc$models[[1]], newdata = test_df, type = "prob")
  rc_xgb = roc(response = test_df$PragueClass, predictor = pr_xg$`X3`)
  coords_xg = coords(rc_xgb, ret = c('sensitivity', 'specificity', 'threshold'), x = "all")
  # There is an entry with 0.80 for sens and 0.714 for spec - 31st res
  opt_dt = coords_xg[20,]


  roc_auc_plot = ggplot(coords_xg, aes(x = 1 - specificity, y = sensitivity)) +
    geom_line(size = 1, colour = "purple") +
    geom_abline(linetype = "dashed") +
    # annotate("point", x = 1 - opt_dt['specificity'], y = opt_dt['sensitivity'], color = 'black', size = 3) +
    # annotate("text", x = 1 - opt_dt['specificity'], y = opt_dt['sensitivity'],
    #          label = sprintf("Threshold: %.2f\nSens: %.2f\nSpec: %.2f",
    #                          opt_dt['threshold'], opt_dt['sensitivity'], opt_dt['specificity']),
    #          vjust = -1, hjust = 1, color = 'black') +
    annotate("point", x = (1 - opt_dt$specificity), y = opt_dt$sensitivity, colour = "black") +
    annotate("text", x = (1 - opt_dt$specificity - 0.05), y = (opt_dt$sensitivity),
             label = paste0("(", signif(opt_dt$threshold, 3), ", ",
                            signif(opt_dt$specificity, 3), ", ", "\n", signif(opt_dt$sensitivity, 3), ")"),
             colour = "black", vjust = -0.5) +
    scale_x_continuous("False Positive Rate") +
    scale_y_continuous("True Positive Rate") +
    theme_minimal() +
    theme(legend.position = "right")

  svg("decisionThresholdRoc.svg")
  print(roc_auc_plot)
  dev.off()
}

make_roc_curve_with_decision_threshold_with_test = function(model, test_data) {
  roc_curve_results = list()
  auc_values = numeric(2)
  # levels(test_data$PragueClass) = make.names(levels(test_data$PragueClass))

  # Get the training ROC curve
  train_predictions = model$pred$`X3`
  roc_curve_train = roc(response = model$pred$obs, predictor = train_predictions)

  train_auc_value = auc(roc_curve_train)

  roc_curve_results[[1]] = data.frame(
    specificity = roc_curve_train$specificities,
    sensitivity = roc_curve_train$sensitivities,
    model = "Training"
  )

  auc_values[1] = train_auc_value

  # Get the test ROC curve
  test_predictions = predict(model, newdata = test_data, type = "prob")
  roc_curve_test = roc(response = test_data$PragueClass, predictor = test_predictions$`X3`)

  test_auc_value = auc(roc_curve_test)
  coords_xg = coords(roc_curve_test, , ret = c('sensitivity', 'specificity', 'threshold'), x = "all")

  roc_curve_results[[2]] = data.frame(
    specificity = roc_curve_test$specificities,
    sensitivity = roc_curve_test$sensitivities,
    model = "Test"
  )

  # There is an entry with 0.80 for sens and 0.714 for spec - 31st res
  opt_dt = coords_xg[20,]

  auc_values[2] = test_auc_value


  roc_curves = do.call(rbind, roc_curve_results)

  # Plot the ROC curve - row 20001
  spec = 0.6000000
  sens = 0.6114638
  thresh = 0.2272270
  roc_auc_plot = ggplot(roc_curves, aes(1 - specificity, sensitivity, color = factor(model))) +
    geom_line(size = 1) +
    geom_abline(linetype = "dashed") +
    annotate("point", x = (1 - spec), y = sens, colour = "black") +
    annotate("text", x = (1 - spec - 0.05), y = (sens),
             label = paste0("(", signif(0.2272270, 3), ", ",
                            signif(spec, 3), ", ", "\n", signif(sens, 3), ")"),
             colour = "black", vjust = -0.5) +
    scale_x_continuous("1-Specificity") +
    scale_y_continuous("Sensitivity") +
    labs(color = "Model") +
    theme(legend.position = "right") +
    #add AUC values in title
    ggtitle(paste0("Test AUC: ", signif(test_auc_value, 3), " Train AUC: ", signif(train_auc_value, 3)))

  svg("decisionThresholdRocWithTest.svg")
  print(roc_auc_plot)
  dev.off()
  return(roc_auc_plot)
}

make_test_rocs_multiple_models = function(models, test_data) {
  number_of_models = length(models)
  roc_curve_results = list()
  auc_values = numeric(number_of_models)
  levels(test_data$PragueClass) = make.names(levels(test_data$PragueClass))

  for (i in 1:number_of_models) {
    model = models[[i]]
    predictions = model$pred$`X3`
    roc_curve = roc(response = model$pred$obs, predictor = predictions) #,
    # Uncomment this if you want to change it to the +ve class - will need to pass train data as an arg
    #levels = rev(levels(train_data$PragueClass)))

    auc_value = auc(roc_curve)
    roc_coords = coords(roc_curve)

    roc_curve_results[[i]] = data.frame(
      specificity = roc_curve$specificities,
      sensitivity = roc_curve$sensitivities,
      model = names(models)[i]
    )
    auc_values[i] = auc_value
  }

  roc_curves = do.call(rbind, roc_curve_results)
  roc_auc_plot = ggplot(roc_curves, aes(1 - specificity, sensitivity, color = factor(model))) +
    geom_line(size = 1) +
    geom_abline(linetype = "dashed") +
    scale_x_continuous("1-Specificity") +
    scale_y_continuous("Sensitvity") +
    labs(color = "Algorithm")

  roc_auc_plot = plot_adder(roc_auc_plot)

  return(list(plots = roc_auc_plot,
              aucs = auc_values,
              models = models))
}

make_test_roc_curves_for_best_models = function() {
  rf_best_model = rf_top_three$roc$models[[1]]
  knn_best_model = km_top_three$roc$models[[1]]
  lr_best_model = lr_top_three$roc$models[[1]]
  mlp_best_model = mlp_top_three$roc$models[[1]]
  xgb_best_model = xgb_top_three$roc$models[[1]]

  models = list(rf_best_model, knn_best_model, lr_best_model, mlp_best_model, xgb_best_model)
  names(models) = c("RF", "KM", "LR", "MLP", "XGB")

  roc_curve = make_test_rocs_multiple_models(models, test_df)

  svg("all_roc_test.svg")
  print(roc_curve$plots)
  dev.off()


  write.csv(roc_curve$aucs, file = "all_roc_test.csv")

  return(roc_curve)
}
