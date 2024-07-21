effect_size_calculation <- function(raw_data,Group,index,control_group){
  require(tidyverse)
  require(metafor)
  

  group_col <- sym(Group)
  index_col <- sym(index)
  

  if(!control_group %in% raw_data[[Group]]){
    stop(paste("Data does not contain specified control group:", control_group))
  }
  

  control_data <- raw_data %>%
    filter(!!group_col == control_group)
  
  CMean <- mean(control_data[[index]], na.rm = TRUE)
  CSD <- sd(control_data[[index]], na.rm = TRUE)
  CN <- nrow(control_data)
  
  summary_data <- raw_data %>%
    filter(!!group_col != control_group) %>%  
    group_by(!!group_col) %>%
    summarize(
      TMean = mean(!!index_col, na.rm = TRUE),
      TSD = sd(!!index_col, na.rm = TRUE),
      TN = n()  
    ) %>%
    mutate(
      CMean = CMean,
      CSD = CSD,
      CN = CN
    )
  
  if(nrow(summary_data) == 0){
    stop("No data available after filtering specified control group.")
  }
  
  
  d2 <- escalc(measure = "ROM",
               data = summary_data, 
               m1i = TMean, 
               sd1i = TSD, 
               n1i = TN, 
               m2i = CMean, 
               sd2i = CSD, 
               n2i = CN)
  
  
  if(nrow(d2) == 0){
    stop("Effect size calculation resulted in no data.")
  }
  
  
 
  r1 <- rma(yi, vi, data = d2, method = "REML")
  summary(r1)
  
  
  kk1 <- summary(r1)
  
  
  names(d2)[1] <- "Group"
  subgroups <- unique(d2$Group)  
  subgroup_results <- list()  
  
  for (subgroup in subgroups) {
    subgroup_data <- subset(d2, Group == subgroup)  
    subgroup_result <- rma(yi, vi, data = subgroup_data, method = "REML")  # 计算当前亚组的效应值
    subgroup_results[[subgroup]] <- subgroup_result  
  }
  
  
  for (subgroup in subgroups) {
    print(subgroup_results[[subgroup]])
  }
  
  
  
  combined_data1 <- data.frame(
    estimate = kk1$b,
    se = kk1$se,
    pval = kk1$pval
  )
  
  
  combined_subgroup_data <- list()
  
  for (subgroup in subgroups) {
    subgroup_result <- subgroup_results[[subgroup]]
    combined_subgroup_data[[subgroup]] <- data.frame(
      estimate = subgroup_result$b,
      se = subgroup_result$se,
      pval = subgroup_result$pval
    )
  }
  
  
  combined_data_subgroup <- do.call(rbind, combined_subgroup_data)
  
 
  combined_data_subgroup$Subgroup <- rownames(combined_data_subgroup)
  
  
  data_df <- combined_data_subgroup
  
  data_df$sig <- "ns"
  
 
  for (i in 1:nrow(data_df)) {
    if (data_df$pval[i] < 0.001) {
      data_df$sig[i] <- "***"
    } else if (data_df$pval[i] < 0.01) {
      data_df$sig[i] <- "**"
    } else if (data_df$pval[i] < 0.05) {
      data_df$sig[i] <- "*"
    }
  }
  return(data_df)
}