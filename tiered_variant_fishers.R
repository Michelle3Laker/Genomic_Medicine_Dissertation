library(tidyverse)
setwd("C:/Users/Michelle/OneDrive/MPhil Genomic Medicine/Dissertation Work/CKD")

#READ IN DATA
tiered.data <- read_tsv("CKD Tables for R (sre) - tiered variants.tsv")
tiered.data <- drop_na(tiered.data)

#Tier 1 general
tier1 <- select (tiered.data, 1:3)
tier1 <- tier1 %>%
  mutate(Not.Tier1 = Total - Tier1)

group_pairs <- combn(tier1$Ethnicity, 2, simplify = FALSE)
perform_fisher_test <- function(pair, tier1) {
  group1 <- pair[1]
  group2 <- pair[2]
  
  tier1_group1 <- tier1 %>% filter(Ethnicity == group1) %>% select(Tier1) %>% pull()
  not.tier1_group1 <- tier1 %>% filter(Ethnicity == group1) %>% select(Not.Tier1) %>% pull()
  tier1_group2 <- tier1 %>% filter(Ethnicity == group2) %>% select(Tier1) %>% pull()
  not.tier1_group2 <- tier1 %>% filter(Ethnicity == group2) %>% select(Not.Tier1) %>% pull()
  
  contingency_table <- matrix(c(
    tier1_group1, not.tier1_group1,
    tier1_group2, not.tier1_group2
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
tier1.fisher <- group_pairs %>%
  map_df(~ perform_fisher_test(.x, tier1))
x <- filter(tier1.fisher, FisherPValue < 0.05)

tier1.fisher <- tier1.fisher %>%
  mutate(AdjustedP = p.adjust(FisherPValue, method = "BH")) #not used given the small sample sizes in this patient cohort

#Tier1 fisher's exact test compared to White: British
Non.White <- filter(tier1, Ethnicity != "White: British")
White <- filter(tier1, Ethnicity == "White: British") ##contingency tables


fisher_tier1 <- tier1 %>%
  filter(Ethnicity != "White: British") %>%
  rowwise() %>%
  mutate(
    test_result = list(fisher.test(matrix(c(
      White$Tier1, White$Not.Tier1,
      Tier1, Not.Tier1
    ), nrow = 2))),
    p.value = test_result$p.value,
    odds.ratio = (White$Tier1 * Not.Tier1) /
      (White$Not.Tier1 * Tier1)
  ) %>%
  select(Ethnicity, p.value, odds.ratio) #no need for p value corrections cause of small sample size


#Tier 2 general
tier2 <- select (tiered.data, 1, 2, 4)
tier2 <- tier2 %>%
  mutate(Not.Tier2 = Total - Tier2)

group_pairs <- combn(tier2$Ethnicity, 2, simplify = FALSE)
perform_fisher_test <- function(pair, tier2) {
  group1 <- pair[1]
  group2 <- pair[2]
  
  tier2_group1 <- tier2 %>% filter(Ethnicity == group1) %>% select(Tier2) %>% pull()
  not.tier2_group1 <- tier2 %>% filter(Ethnicity == group1) %>% select(Not.Tier2) %>% pull()
  tier2_group2 <- tier2 %>% filter(Ethnicity == group2) %>% select(Tier2) %>% pull()
  not.tier2_group2 <- tier2 %>% filter(Ethnicity == group2) %>% select(Not.Tier2) %>% pull()
  
  contingency_table <- matrix(c(
    tier2_group1, not.tier2_group1,
    tier2_group2, not.tier2_group2
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
tier2.fisher <- group_pairs %>%
  map_df(~ perform_fisher_test(.x, tier2))
y <- filter(tier2.fisher, FisherPValue < 0.05)

tier2.fisher <- tier2.fisher %>%
  mutate(AdjustedP = p.adjust(FisherPValue, method = "BH")) #not used given the small sample sizes in this patient cohort

#Tier2 fisher's exact test compared to White: British
Non.White <- filter(tier2, Ethnicity != "White: British")
White <- filter(tier2, Ethnicity == "White: British") ##contingency tables


fisher_tier2 <- tier2 %>%
  filter(Ethnicity != "White: British") %>%
  rowwise() %>%
  mutate(
    test_result = list(fisher.test(matrix(c(
      White$Tier2, White$Not.Tier2,
      Tier2, Not.Tier2
    ), nrow = 2))),
    p.value = test_result$p.value,
    odds.ratio = (White$Tier2 * Not.Tier2) /
      (White$Not.Tier2 * Tier2)
  ) %>%
  select(Ethnicity, p.value, odds.ratio) #no need for p value corrections cause of small sample size


#Tier 1&2 for Genetic Diagnosis
tierGD <- select (tiered.data, 1, 2, 5)
tierGD <- tierGD %>%
  mutate(Not.TierGD = Total - TierGD)

group_pairs <- combn(tierGD$Ethnicity, 2, simplify = FALSE)
perform_fisher_test <- function(pair, tierGD) {
  group1 <- pair[1]
  group2 <- pair[2]
  
  tierGD_group1 <- tierGD %>% filter(Ethnicity == group1) %>% select(TierGD) %>% pull()
  not.tierGD_group1 <- tierGD %>% filter(Ethnicity == group1) %>% select(Not.TierGD) %>% pull()
  tierGD_group2 <- tierGD %>% filter(Ethnicity == group2) %>% select(TierGD) %>% pull()
  not.tierGD_group2 <- tierGD %>% filter(Ethnicity == group2) %>% select(Not.TierGD) %>% pull()
  
  contingency_table <- matrix(c(
    tierGD_group1, not.tierGD_group1,
    tierGD_group2, not.tierGD_group2
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
tierGD.fisher <- group_pairs %>%
  map_df(~ perform_fisher_test(.x, tierGD))
w <- filter(tierGD.fisher, FisherPValue < 0.05)

tierGDfisher <- tierGD.fisher %>%
  mutate(AdjustedP = p.adjust(FisherPValue, method = "BH")) #not used given the small sample sizes in this patient cohort

#TierGD fisher's exact test compared to White: British
Non.White <- filter(tierGD, Ethnicity != "White: British")
White <- filter(tierGD, Ethnicity == "White: British") ##contingency tables

fisher_tierGD <- tierGD %>%
  filter(Ethnicity != "White: British") %>%
  rowwise() %>%
  mutate(
    test_result = list(fisher.test(matrix(c(
      White$TierGD, White$Not.TierGD,
      TierGD, Not.TierGD
    ), nrow = 2))),
    p.value = test_result$p.value,
    odds.ratio = (White$TierGD * Not.TierGD) /
      (White$Not.TierGD * TierGD)
  ) %>%
  select(Ethnicity, p.value, odds.ratio) #no need for p value corrections cause of small sample size


#Tier 3 general
tier3 <- select (tiered.data, 1, 2, 6)
tier3 <- tier3 %>%
  mutate(Not.Tier3 = Total - Tier3)

group_pairs <- combn(tier3$Ethnicity, 2, simplify = FALSE)
perform_fisher_test <- function(pair, tier3) {
  group1 <- pair[1]
  group2 <- pair[2]
  
  tier3_group1 <- tier3 %>% filter(Ethnicity == group1) %>% select(Tier3) %>% pull()
  not.tier3_group1 <- tier3 %>% filter(Ethnicity == group1) %>% select(Not.Tier3) %>% pull()
  tier3_group2 <- tier3 %>% filter(Ethnicity == group2) %>% select(Tier3) %>% pull()
  not.tier3_group2 <- tier3 %>% filter(Ethnicity == group2) %>% select(Not.Tier3) %>% pull()
  
  contingency_table <- matrix(c(
    tier3_group1, not.tier3_group1,
    tier3_group2, not.tier3_group2
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
tier3.fisher <- group_pairs %>%
  map_df(~ perform_fisher_test(.x, tier3))
z <- filter(tier3.fisher, FisherPValue < 0.05)

tier3.fisher <- tier3.fisher %>%
  mutate(AdjustedP = p.adjust(FisherPValue, method = "BH")) #not used given the small sample sizes in this patient cohort

#Tier3 fisher's exact test compared to White: British
Non.White <- filter(tier3, Ethnicity != "White: British")
White <- filter(tier3, Ethnicity == "White: British") ##contingency tables

fisher_tier3 <- tier3 %>%
  filter(Ethnicity != "White: British") %>%
  rowwise() %>%
  mutate(
    test_result = list(fisher.test(matrix(c(
      White$Tier3, White$Not.Tier3,
      Tier3, Not.Tier3
    ), nrow = 2))),
    p.value = test_result$p.value,
    odds.ratio = (White$Tier3 * Not.Tier3) /
      (White$Not.Tier3 * Tier3)
  ) %>%
  select(Ethnicity, p.value, odds.ratio) #no need for p value corrections cause of small sample size
