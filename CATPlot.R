## Create Concordance at the Top Plot from two lists of microbial biomarkers

## Packages -------------
suppressPackageStartupMessages({
  library(ffpe)
  library(dplyr)
  library(tidyr)
  library(tidyverse)
  library(stats)
})

list1 <- readRDS("~/EOCRC/locrc_res.rds")
list2 <- readRDS("~/EOCRC/eocrc_res.rds")

format <- function(list) {
  ## Need two ranked lists of microbial signatures in the order from most significant
    ## to least significant based on adjusted p-values from the output of 
    ## differential analyses
  ranked_list <- 
    list[order(list$padj, decreasing =
                          FALSE),]
  
  ## Convert to data frame and drop NA for padj
  ranked_list_df <- as.data.frame(ranked_list) %>% drop_na(padj)
  
  return(ranked_list_df)

}

CAT_main <- function() {
  ranked_list1 <- format(list1)
  ranked_list2 <- format(list2)
  ## Creating 95% CI for distribution of concordance values between two ranked 
    ## lists of microbial signatures sampled 100 times
  set.seed(1)
  concordance_repeat <- replicate(100, CATplot(sample(rownames(ranked_list1), 
                                                      replace = FALSE), 
                                               rownames(ranked_list2), 
                                               make.plot = FALSE)$concordance)
  
  quants = apply(concordance_repeat, 1, function(r) quantile(r, c(0.025, 0.975)))
  
  cat_ci <- CATplot(rownames(ranked_list1), rownames(ranked_list2), make.plot = TRUE)
  lines(x = cat_ci$rank, y = quants["2.5%",], lty = 3)
  lines(x = cat_ci$rank, y = quants["97.5%",], lty = 3)
  legend("bottomright",lty=1:3,legend=c("Actual concordance",
                                        "Concordance expected by chance", 
                                        "95% Confidence Intervals"))

  return(cat_ci)
}

CAT_main()
