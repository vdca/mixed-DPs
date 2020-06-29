
#-------------------------------------------------------------------------------
# Global
#-------------------------------------------------------------------------------

rm(list=ls())     # remove previous objects from workspace

library(tidyverse)

# statistical models
library(lmerTest)
library(blme)
library(broom)

# show GAS examples and proportions in table form
library(flextable)
library(officer)
std_border <- fp_border(color="gray")
thick_border <- fp_border(color="black", width = 2)

#-------------------------------------------------------------------------------
# Helper functions
#-------------------------------------------------------------------------------

# combi_props() takes n grouping variables and computes proportion of x for each combination
combi_props <- function(dx, depvar, ...) {
  dx %>% 
    count(..., {{depvar}}, .drop = T) %>% 
    group_by(...) %>% 
    mutate(p = n/sum(n)) %>% 
    ungroup() %>% 
    mutate(proportion = round(p*100, 0),
           proportion = paste0(proportion, '%'),
           p = NULL)
}

# print tables for visual inspection
flx <- function(x, ...) {
  x %>% 
    flextable() %>% 
    bold(i = 1, part = "header", bold = T) %>% 
    autofit(add_w = 0, add_h = 0.01) %>% 
    fontsize(size = 10, part = 'all')
}

# print fixed effects of mixed model
print_model <- function(mx) {
  mx %>% 
    tidy(effects = 'fixed') %>% 
    mutate(signif = gtools::stars.pval(p.value),
           odds.ratio = exp(estimate),
           z.value = statistic,
           statistic = NULL,
           p.value = format.pval(p.value, digits = 2),
           term = str_replace_all(term, 'es.gender', 'es.gender_'),
           term = str_replace_all(term, 'noun.end', 'noun.end_'),
           term = str_replace_all(term, 'root.end', 'root.end_'),
           term = str_replace_all(term, 'L1_S1', 'L1S_yes'),
           term = str_replace_all(term, 'S_area1', 'S_area_yes')) %>% 
    mutate_if(is.numeric, ~ round(., digits = 2)) %>% 
    select(term, estimate, std.error, odds.ratio, z.value, p.value, signif) %>% 
    flx()
}

#-------------------------------------------------------------------------------
# Read data
#-------------------------------------------------------------------------------

# read responses and background
responses <- read_tsv('../data/responses.tsv')
bg_edb <- read_tsv('../data/bg_edb.tsv')

# merge responses and background info,
# and convert some variables to factors.
d <- responses %>%
  left_join(bg_edb) %>% 
  mutate_at(vars(noun.end, root.end, detGender, L1_S, L1_B, S_area), factor)

# create subset of data for illustration purposes
d_illustration <- d %>% 
  mutate(gloss = paste0(det, ' ', noun, ' ‘', es, '.', toupper(es.gender), ', ', en, '’'),
         gloss2 = paste0(noun, ' ‘', es, '.', toupper(es.gender), ', ', en, '’')) %>% 
  filter(lemma != 'komun') %>% # filter out komun, just because spa-equivalent is multiword
  select(example = gloss, det, noun, es, en, det.gender = detGender,
         es.gender, noun.end, root.end, L1_S, L1_B, S_area, gloss2) %>% 
  mutate(L1_S = recode(L1_S, '0' = 'no', '1' = 'yes')) %>% 
  rename(L1S = L1_S)

#-------------------------------------------------------------------------------
# Table A.1: proportion of F/M determiner by es.gender and noun.end
#-------------------------------------------------------------------------------

# proportion of: F/M det by es.gender * noun.end
d_an <- d_illustration %>% 
  combi_props(det.gender, es.gender, noun.end)

# examples for: F/M det by es.gender * noun.end
d_an <- d_illustration %>% 
  arrange(L1S) %>% 
  select(es.gender, noun.end, det.gender, example) %>% 
  group_by(es.gender, noun.end, det.gender) %>%
  slice(1) %>% 
  ungroup() %>% 
  right_join(d_an) %>% 
  select(-example, example)

# print table: F/M det by es.gender * noun.end
d_an %>% 
  flextable() %>% 
  bold(i = 1, part = "header", bold = T) %>%
  merge_h_range(i = 1, j1 = 4, j2 = 4, part = 'header') %>% 
  align_text_col('center') %>% 
  align(i = 1, j = 4, align = 'center', part = 'header') %>% 
  align(j = 6, align = 'left', part = 'body') %>% 
  hline(i = c(2,4,6,8,10,11,13), border = std_border) %>%
  hline(i = c(8), border = thick_border) %>% 
  autofit(add_w = -0.05, add_h = 0.01) %>% 
  fontsize(size = 10, part = 'all')

#-------------------------------------------------------------------------------
# Table A.2: proportion of F/M determiner by es.gender and root.end (only nouns ending in -a)
#-------------------------------------------------------------------------------

# proportion of: F/M det by es.gender * root.end (only -a noun-end)
d_anr <- d_illustration %>%
  filter(noun.end == 'a') %>% 
  combi_props(det.gender, es.gender, noun.end, root.end)

# examples for: F/M det by es.gender * root.end (only -a noun-end)
d_anr <- d_illustration %>% 
  arrange(L1S) %>% 
  select(es.gender, noun.end, root.end, det.gender, example) %>% 
  group_by(es.gender, noun.end, root.end, det.gender) %>%
  slice(1) %>% 
  right_join(d_anr) %>% 
  select(-example, example)

# print table: F/M det by es.gender * root.end (only -a noun-end)
d_anr %>% 
  flextable() %>% 
  bold(i = 1, part = "header", bold = T) %>%
  merge_h_range(i = 1, j1 = 4, j2 = 4, part = 'header') %>% 
  align_text_col('center') %>% 
  align(i = 1, j = 4, align = 'center', part = 'header') %>% 
  align(j = 7, align = 'left', part = 'body') %>% 
  hline(i = c(2,4,6,8,10,12,14), border = std_border) %>%
  hline(i = c(8), border = thick_border) %>% 
  autofit(add_w = -0.05, add_h = 0.01) %>% 
  fontsize(size = 10, part = 'all')

#-------------------------------------------------------------------------------
# Table A.3: proportion of F/M determiner by es.gender and noun.end, separated by L1S
#-------------------------------------------------------------------------------

# proportion of: F/M det by es.gender * noun.end * L1S
d_lan <- d_illustration %>%
  combi_props(det.gender, L1S, es.gender, noun.end)

# examples for: F/M det by es.gender * noun.end * L1S
d_lan <- d_illustration %>% 
  arrange(L1S) %>% 
  select(L1S, es.gender, noun.end, det.gender, example) %>% 
  group_by(L1S, es.gender, noun.end, det.gender) %>%
  slice(1) %>% 
  right_join(d_lan) %>% 
  select(-example, example)

# print table: F/M det by es.gender * noun.end * L1S
d_lan %>% 
  arrange(L1S) %>% 
  select(L1S, es.gender, noun.end, det.gender, example) %>% 
  group_by(L1S, es.gender, noun.end, det.gender) %>%
  slice(1) %>% 
  right_join(d_lan) %>% 
  select(-example, example) %>% 
  as_grouped_data(c("L1S")) %>% 
  as_flextable() %>% 
  bold(j = 1, i = ~ !is.na(L1S), bold = T, part = "body") %>% 
  bg(i = c(1,16),  bg = "#D3D3D3", part = "body") %>% 
  bold(i = 1, part = "header", bold = T) %>%
  merge_h_range(i = 1, j1 = 4, j2 = 4, part = 'header') %>% 
  align_text_col('center') %>% 
  align(i = 1, j = 4, align = 'center', part = 'header') %>% 
  align(j = 6, align = 'left', part = 'body') %>% 
  hline(i = c(3,4,6,8,10,11,13,18,20,22,24,26,27,29), border = std_border) %>%
  hline(i = c(8,24), border = thick_border) %>% 
  autofit(add_w = -0.05, add_h = 0.01) %>% 
  fontsize(size = 10, part = 'all')

#-------------------------------------------------------------------------------
# Appendix B: mixed-effects logistic models
#-------------------------------------------------------------------------------

#' ## Table B.1: Model 1
#' 
#' det.gender ~ es.gender + noun.end + root.end + (1|participant)

# only cues
phon.glm1 <- d %>%
  glmer(detGender ~ es.gender + noun.end + root.end +
          (1|participant), ., family = 'binomial',
        control=glmerControl(optimizer="bobyqa"))

phon.glm1 %>% 
  print_model %>% 
  hline(i = c(1,2,5), border = std_border)

#' ## Table B.2: Model 2
#' 
#' det.gender ~ es.gender + noun.end + root.end + L1S + S_area +
#' L1S:es.gender + S_area:es.gender + L1S:noun.end + (1|participant)

# full model
blme.m1 <- bglmer(detGender ~ es.gender + noun.end + root.end +
                    L1_S:es.gender + S_area:es.gender + L1_S:noun.end +
                    L1_S + S_area + (1|participant),
                  data=d,
                  family=binomial,
                  fixef.prior = normal,
                  control=glmerControl(optimizer='bobyqa'))

blme.m1 %>%
  print_model %>% 
  hline(i = c(1,2,5,8,12), border = std_border)

