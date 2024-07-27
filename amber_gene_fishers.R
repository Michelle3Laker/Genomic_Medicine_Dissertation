#FISHER'S TESTS FOR GREEN(DIAGNOSTIC) GENES ON EACH PANEL
coded.data <- read_tsv("CKD Tables for R (sre) - colour codes.tsv")
amber.test <- select(coded.data, 1, 2, 5, 6)
amber.test <- amber.test %>%
  mutate(Not.Amber.reg = Total - Amber.reg) %>%
  mutate(Not.Amber.super = Total - Amber.super)

#general fisher's test for regular panel
group_pairs <- combn(amber.test$Ethnicity, 2, simplify = FALSE)
perform_fisher_test <- function(pair, amber.test) {
  group1 <- pair[1]
  group2 <- pair[2]
  
  amber.reg_group1 <- amber.test %>% filter(Ethnicity == group1) %>% select(Amber.reg) %>% pull()
  not.amber.reg_group1 <- amber.test %>% filter(Ethnicity == group1) %>% select(Not.Amber.reg) %>% pull()
  amber.reg_group2 <- amber.test %>% filter(Ethnicity == group2) %>% select(Amber.reg) %>% pull()
  not.amber.reg_group2 <- amber.test %>% filter(Ethnicity == group2) %>% select(Not.Amber.reg) %>% pull()
  
  contingency_table <- matrix(c(
    amber.reg_group1, not.amber.reg_group1,
    amber.reg_group2, not.amber.reg_group2
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
amber.fisher <- group_pairs %>%
  map_df(~ perform_fisher_test(.x, amber.test))
amber.fisher <- amber.fisher %>%
  mutate(AdjustedP = p.adjust(FisherPValue, method = "BH")) 
##what is actually the best correction method? that isn't super conservative

x <- filter(amber.fisher, FisherPValue < 0.05)

#general fisher's test for super panel
group_pairs <- combn(amber.test$Ethnicity, 2, simplify = FALSE)
perform_fisher_test <- function(pair, amber.test) {
  group1 <- pair[1]
  group2 <- pair[2]
  
  amber.super_group1 <- amber.test %>% filter(Ethnicity == group1) %>% select(Amber.super) %>% pull()
  not.amber.super_group1 <- amber.test %>% filter(Ethnicity == group1) %>% select(Not.Amber.super) %>% pull()
  amber.super_group2 <- amber.test %>% filter(Ethnicity == group2) %>% select(Amber.super) %>% pull()
  not.amber.super_group2 <- amber.test %>% filter(Ethnicity == group2) %>% select(Not.Amber.super) %>% pull()
  
  contingency_table <- matrix(c(
    amber.super_group1, not.amber.super_group1,
    amber.super_group2, not.amber.super_group2
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
amber.fisher.super <- group_pairs %>%
  map_df(~ perform_fisher_test(.x, amber.test))
amber.fisher.super <- amber.fisher.super %>%
  mutate(AdjustedP = p.adjust(FisherPValue, method = "BH")) 
##what is actually the best correction method? that isn't super conservative

y <- filter(amber.fisher.super, FisherPValue < 0.05)

#fisher's exact test compared to White: British
Non.White <- filter(amber.test, Ethnicity != "White: British")
White <- filter(amber.test, Ethnicity == "White: British") ##contingency tables


fisher_regular <- amber.test %>%
  filter(Ethnicity != "White: British") %>%
  rowwise() %>%
  mutate(
    test_result = list(fisher.test(matrix(c(
      White$Amber.reg, White$Not.Amber.reg,
      Amber.reg, Not.Amber.reg
    ), nrow = 2))),
    p.value = test_result$p.value,
    odds.ratio = (White$Amber.reg * Not.Amber.reg) /
      (White$Not.Amber.reg * Amber.reg)
  ) %>%
  select(Ethnicity, p.value, odds.ratio) #no need for p value corrections cause of small sample size

##rmr that odds ratio > 1 is the no. of times it's more likely that 
##the White ethnic group has a green gene and vice versa

fisher_super <- amber.test %>%
  filter(Ethnicity != "White: British") %>%
  rowwise() %>%
  mutate(
    test_result = list(fisher.test(matrix(c(
      White$Amber.super, White$Not.Amber.super,
      Amber.super, Not.Amber.super
    ), nrow = 2))),
    p.value = test_result$p.value,
    odds.ratio = (White$Amber.super * Not.Amber.super) /
      (White$Not.Amber.super * Amber.super)
  ) %>%
  select(Ethnicity, p.value, odds.ratio)

##zero odds ratios are simply because there are no people without green genes
##in these groups