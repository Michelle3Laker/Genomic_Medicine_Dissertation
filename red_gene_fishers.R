#FISHER'S TESTS FOR GREEN(DIAGNOSTIC) GENES ON EACH PANEL
coded.data <- read_tsv("CKD Tables for R (sre) - colour codes.tsv")
red.test <- select(coded.data, 1, 2, 7, 8)
red.test <- red.test %>%
  mutate(Not.Red.reg = Total - Red.reg) %>%
  mutate(Not.Red.super = Total - Red.super)

#general fisher's test for regular panel
group_pairs <- combn(red.test$Ethnicity, 2, simplify = FALSE)
perform_fisher_test <- function(pair, red.test) {
  group1 <- pair[1]
  group2 <- pair[2]
  
  red.reg_group1 <- red.test %>% filter(Ethnicity == group1) %>% select(Red.reg) %>% pull()
  not.red.reg_group1 <- red.test %>% filter(Ethnicity == group1) %>% select(Not.Red.reg) %>% pull()
  red.reg_group2 <- red.test %>% filter(Ethnicity == group2) %>% select(Red.reg) %>% pull()
  not.red.reg_group2 <- red.test %>% filter(Ethnicity == group2) %>% select(Not.Red.reg) %>% pull()
  
  contingency_table <- matrix(c(
    red.reg_group1, not.red.reg_group1,
    red.reg_group2, not.red.reg_group2
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
red.fisher <- group_pairs %>%
  map_df(~ perform_fisher_test(.x, red.test))

red.fisher <- red.fisher %>%
  mutate(AdjustedP = p.adjust(FisherPValue, method = "BH")) #not used given the small sample sizes in this patient cohort

x <- filter(red.fisher, FisherPValue < 0.05)

#general fisher's test for super panel
group_pairs <- combn(red.test$Ethnicity, 2, simplify = FALSE)
perform_fisher_test <- function(pair, red.test) {
  group1 <- pair[1]
  group2 <- pair[2]
  
  red.super_group1 <- red.test %>% filter(Ethnicity == group1) %>% select(Red.super) %>% pull()
  not.red.super_group1 <- red.test %>% filter(Ethnicity == group1) %>% select(Not.Red.super) %>% pull()
  red.super_group2 <- red.test %>% filter(Ethnicity == group2) %>% select(Red.super) %>% pull()
  not.red.super_group2 <- red.test %>% filter(Ethnicity == group2) %>% select(Not.Red.super) %>% pull()
  
  contingency_table <- matrix(c(
    red.super_group1, not.red.super_group1,
    red.super_group2, not.red.super_group2
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
red.fisher.super <- group_pairs %>%
  map_df(~ perform_fisher_test(.x, red.test))

red.fisher.super <- red.fisher.super %>%
  mutate(AdjustedP = p.adjust(FisherPValue, method = "BH")) #not used given the small sample sizes in this patient cohort

y <- filter(red.fisher.super, FisherPValue < 0.05)

#fisher's exact test compared to White: British
Non.White <- filter(red.test, Ethnicity != "White: British")
White <- filter(red.test, Ethnicity == "White: British") ##contingency tables


fisher_regular <- red.test %>%
  filter(Ethnicity != "White: British") %>%
  rowwise() %>%
  mutate(
    test_result = list(fisher.test(matrix(c(
      White$Red.reg, White$Not.Red.reg,
      Red.reg, Not.Red.reg
    ), nrow = 2))),
    p.value = test_result$p.value,
    odds.ratio = (White$Red.reg * Not.Red.reg) /
      (White$Not.Red.reg * Red.reg)
  ) %>%
  select(Ethnicity, p.value, odds.ratio) #no need for p value corrections cause of small sample size

##Rmr that odds ratio is the no. of times it's more likely that the White: British group has a red gene and vice versa.

fisher_super <- red.test %>%
  filter(Ethnicity != "White: British") %>%
  rowwise() %>%
  mutate(
    test_result = list(fisher.test(matrix(c(
      White$Red.super, White$Not.Red.super,
      Red.super, Not.Red.super
    ), nrow = 2))),
    p.value = test_result$p.value,
    odds.ratio = (White$Red.super * Not.Red.super) /
      (White$Not.Red.super * Red.super)
  ) %>%
  select(Ethnicity, p.value, odds.ratio)

##zero odds ratios are simply because there are no people without red gene variants in these groups