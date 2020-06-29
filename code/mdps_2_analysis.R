
#-------------------------------------------------------------------------------
# Global
#-------------------------------------------------------------------------------

rm(list=ls())     # remove previous objects from workspace

library(tidyverse)

# statistical models
library(lmerTest)
library(blme)

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

#-------------------------------------------------------------------------------
# Model 1: only cues (without participant background)
#-------------------------------------------------------------------------------

# analysis 1:
# det.gender = es.gender + noun.end + root.end

phon.glm1 <- d %>%
  glmer(detGender ~ es.gender + noun.end + root.end +
          (1|participant), ., family = 'binomial',
        control=glmerControl(optimizer="bobyqa"))

summary(phon.glm1)
performance::r2(phon.glm1)

# likelihood ratio tests:
# test each cue by comparing full model to model with relevant cue removed.

# analysis 2: does phon.end improve the model?
# --> it does

phon.glm2 <- d %>%
  glmer(detGender ~ es.gender + root.end +
          (1|participant), ., family = 'binomial',
        control=glmerControl(optimizer="bobyqa"))

anova(phon.glm1, phon.glm2)

# analysis 3: does root.end improve the model?
# --> it does

phon.glm3 <- d %>%
  glmer(detGender ~ es.gender + noun.end +
          (1|participant), ., family = 'binomial',
        control=glmerControl(optimizer="bobyqa"))

anova(phon.glm1, phon.glm3)

# analysis 4: does es.gender improve the model?
# --> it does

phon.glm4 <- d %>%
  glmer(detGender ~ noun.end + root.end +
          (1|participant), ., family = 'binomial',
        control=glmerControl(optimizer="bobyqa"))

anova(phon.glm1, phon.glm4)

#-------------------------------------------------------------------------------
# Model 2: cues + participant background (L1 and sociolinguistic area)
#-------------------------------------------------------------------------------

# run stepwise model selection on a saturated model.
# warning: takes some time to run (ca. 5 minutes);
# read saved version instead:
detgen.full <- readRDS('../data/detgenfull_model.rds')

# detgen.full <- StatisticalModels::GLMERSelect(
#   modelData = d,
#   responseVar = "detGender",
#   fitFamily = "binomial",
#   fixedFactors = c('es.gender', 'noun.end', 'root.end', 'L1_S', 'L1_B', 'S_area'),
#   fixedInteractions = c('L1_S:es.gender', 'L1_B:es.gender', 'S_area:es.gender',
#                         'L1_S:noun.end', 'L1_B:noun.end', 'S_area:noun.end',
#                         'L1_S:root.end', 'L1_B:root.end', 'S_area:root.end'),
#   randomStruct = "(1|participant)",
#   verbose = T)

# detgen.full <- readRDS('../data/detgenfull_model.rds')

# the final model step-wise selected includes:
# 
# 1. all three cues (es.gender, noun.end, root.end)
# 2. L1_S and S_area, although not statistically significant;
#    (L1_B not selected by model, simultaneous bilinguals seem to pattern with L1_S.)
# 3. interactions between L1_S/S_area and es.gender (more analogical if L1_S or S_area)
#                         L1_S and noun.end (less phonological if L1_S)

detgen.full$final.call
detgen.full$model %>% summary()

# this model has a 'complete separation' problem, which causes convergence issues.
#   see: https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#penalizationhandling-complete-separation
#   see: https://stats.stackexchange.com/search?q=complete+separation+mixed
# the complete separation stems from the interaction between L1_S and noun.end:
# all nouns with final -i produced by L1_S=0 participants have M determiner.
d %>% 
  filter(noun.end == 'i') %>% 
  count(L1_S, detGender, .drop = F)

# re-run the selected model.
# blme pakcage, bglmer() function, bayesian glmer
# solves the 'complete separation' issue (see above)
blme.m1 <- blme::bglmer(detGender ~ es.gender + noun.end + root.end +
                          L1_S:es.gender + S_area:es.gender + L1_S:noun.end +
                          L1_S + S_area + (1|participant),
                        data=d,
                        family=binomial,
                        fixef.prior = normal,
                        control=glmerControl(optimizer='bobyqa'))
summary(blme.m1)
performance::r2(blme.m1)

