library(tidyverse)
setwd("C:/Users/Michelle/OneDrive/MPhil Genomic Medicine/Dissertation Work/CKD")

#FISHER'S TEST FOR CENSUS COMPARISON
census.data <- read_tsv("edited.Census comparison.tsv")
census.data <- select(census.data, -6, -7)
census.data <- drop_na(census.data)

perform_fishers_test <- function(Total, GenPop, Ethnicity) {
  # Create a contingency table for the group
  contingency_table <- matrix(c(Total, GenPop,
                                sum(census.data$Total) - Total,
                                sum(census.data$GenPop) - GenPop),
                              nrow = 2, byrow = TRUE)
  rownames(contingency_table) <- c("Total", "GenPop")
  colnames(contingency_table) <- c(Ethnicity, paste("Not", Ethnicity))
  # Perform Fisher's exact test
  fisher_test <- fisher.test(contingency_table)
  
  # Return results
  return(list(group = Ethnicity, p_value = fisher_test$p.value, odds_ratio = fisher_test$estimate))
}

census.fisher <- lapply(1:nrow(census.data), function(i) {
  perform_fishers_test(census.data$Total[i], census.data$GenPop[i], census.data$Ethnicity[i])
})

census.fisher <- do.call(rbind, lapply(census.fisher, as.data.frame))
census.fisher <- do.call(data.frame, census.fisher)
census.fisher$adjusted_p_value <- p.adjust(census.fisher$p_value, method = "bonferroni") #not used given the small sample sizes in this patient cohort


#FISHER'S TESTS FOR GREEN(DIAGNOSTIC) GENES ON EACH PANEL
coded.data <- read_tsv("CKD Tables for R (sre) - colour codes.tsv")
green.test <- select(coded.data, 1:4)
green.test <- green.test %>%
  mutate(Not.Green.reg = Total - Green.reg) %>%
  mutate(Not.Green.super = Total - Green.super)

#general fisher's test for regular panel
group_pairs <- combn(green.test$Ethnicity, 2, simplify = FALSE)
perform_fisher_test <- function(pair, green.test) {
  group1 <- pair[1]
  group2 <- pair[2]

  green.reg_group1 <- green.test %>% filter(Ethnicity == group1) %>% select(Green.reg) %>% pull()
  not.green.reg_group1 <- green.test %>% filter(Ethnicity == group1) %>% select(Not.Green.reg) %>% pull()
  green.reg_group2 <- green.test %>% filter(Ethnicity == group2) %>% select(Green.reg) %>% pull()
  not.green.reg_group2 <- green.test %>% filter(Ethnicity == group2) %>% select(Not.Green.reg) %>% pull()
 
  contingency_table <- matrix(c(
    green.reg_group1, not.green.reg_group1,
    green.reg_group2, not.green.reg_group2
  ), nrow = 2)
  
  #Perform Fisher's exact test
  test_result <- fisher.test(contingency_table)
  
  return(tibble(
    Group1 = group1,
    Group2 = group2,
    FisherPValue = test_result$p.value,
    OddsRatio = test_result$estimate
  ))
}
# Apply Fisher's exact test to all pairwise combinations
green.fisher <- group_pairs %>%
  map_df(~ perform_fisher_test(.x, green.test))

green.fisher <- green.fisher %>%
  mutate(AdjustedP = p.adjust(FisherPValue, method = "BH"))  #not used given the small sample sizes in this patient cohort

x <- filter(green.fisher, FisherPValue < 0.05)

#general fisher's test for super panel
group_pairs <- combn(green.test$Ethnicity, 2, simplify = FALSE)
perform_fisher_test <- function(pair, green.test) {
  group1 <- pair[1]
  group2 <- pair[2]
  
  green.super_group1 <- green.test %>% filter(Ethnicity == group1) %>% select(Green.super) %>% pull()
  not.green.super_group1 <- green.test %>% filter(Ethnicity == group1) %>% select(Not.Green.super) %>% pull()
  green.super_group2 <- green.test %>% filter(Ethnicity == group2) %>% select(Green.super) %>% pull()
  not.green.super_group2 <- green.test %>% filter(Ethnicity == group2) %>% select(Not.Green.super) %>% pull()
  
  contingency_table <- matrix(c(
    green.super_group1, not.green.super_group1,
    green.super_group2, not.green.super_group2
  ), nrow = 2)
  
  #Perform Fisher's exact test
  test_result <- fisher.test(contingency_table)
  
  return(tibble(
    Group1 = group1,
    Group2 = group2,
    FisherPValue = test_result$p.value,
    OddsRatio = test_result$estimate
  ))
}
# Apply Fisher's exact test to all pairwise combinations
green.fisher.super <- group_pairs %>%
  map_df(~ perform_fisher_test(.x, green.test))

green.fisher.super <- green.fisher.super %>%
  mutate(AdjustedP = p.adjust(FisherPValue, method = "BH")) #not used given the small sample sizes in this patient cohort

y <- filter(green.fisher.super, FisherPValue < 0.05)

#fisher's exact test compared to White: British
Non.White <- filter(green.test, Ethnicity != "White: British")
White <- filter(green.test, Ethnicity == "White: British") ##contingency tables


fisher_regular <- green.test %>%
  filter(Ethnicity != "White: British") %>%
  rowwise() %>%
  mutate(
    test_result = list(fisher.test(matrix(c(
      White$Green.reg, White$Not.Green.reg,
      Green.reg, Not.Green.reg
    ), nrow = 2))),
    p.value = test_result$p.value,
    odds.ratio = (White$Green.reg * Not.Green.reg) /
      (White$Not.Green.reg * Green.reg)
  ) %>%
  select(Ethnicity, p.value, odds.ratio) #no need to correct p value for multiple comparisons cause of small sample size

##rmr that odds ratio > 1 is the no. of times it's more likely that the White: British group has a green gene variant and vice versa

fisher_super <- green.test %>%
  filter(Ethnicity != "White: British") %>%
  rowwise() %>%
  mutate(
    test_result = list(fisher.test(matrix(c(
      White$Green.super, White$Not.Green.super,
      Green.super, Not.Green.super
    ), nrow = 2))),
    p.value = test_result$p.value,
    odds.ratio = (White$Green.super * Not.Green.super) /
      (White$Not.Green.super * Green.super)
  ) %>%
  select(Ethnicity, p.value, odds.ratio) #no need to correct p value for multiple comparisons cause of small sample size

##zero odds ratios are simply because there are no people without green gene variants in these groups