
#-------------------------------------------------------------------------------
# Global
#-------------------------------------------------------------------------------

rm(list=ls())     # remove previous objects from workspace

library(tidyverse)

#-------------------------------------------------------------------------------
# Load response data
#-------------------------------------------------------------------------------

responses <- read_tsv('../data/responses.tsv')

#-------------------------------------------------------------------------------
# Load participants' background data
#-------------------------------------------------------------------------------

bg <- read_tsv('../data/participant_background.tsv')

#-------------------------------------------------------------------------------
# Add sociolinguistically dominant language for participant's (1) age group and (2) town
#-------------------------------------------------------------------------------

# Get knowledge-of-Basque data (by town and age group) from Euskararen Datu Basea:
# http://www.soziolinguistika.eus/edb/index.php?erakus=orriak&zo=106

edb <- read_csv('../data/edb_gaitasuna_adina.csv', skip = 1)

# clean data.
# description of variables:
# Euskaldunak = B_speakers = bilingual Basque+Spanish speakers
# Erdaldunak = S_speakers = monolingual Spanish speakers
# S_area = area where Spanish is the sociolinguistically dominant language;
#          i.e. monolingual Spanish speakers constitute more than half of the population.
edb_age <- edb %>% 
  mutate(area = str_replace(Lurraldea, ' \\(.{1,2}\\)', ''),
         area = str_replace_all(area, 'Ã±', 'ñ')) %>%
  separate(adina, c('age.range', 'temp'), sep = ' ') %>% 
  mutate(population = Erdaldunak + Euskaldunak + `Ia Euskaldunak`,
         B_speakers = (Euskaldunak) / population,
         S_speakers = (population-Euskaldunak) / population,
         S_area = if_else(S_speakers >= .5, 1, 0),
         area = if_else(area == 'Iruñea', 'Iruñerria', area),
         area = if_else(area == 'Ezkio-Itsaso', 'Ezkio', area)) %>%
  select(area, age.range, B_speakers, S_speakers, S_area, age.range)

# assign participant age to age.range.
# merge EDB data with participant data.
# all participants within 20--50 range.
summary(bg$age)

# hence: if <= 34 --> age.range = 15-34
#        else     --> age.range = 35-64
bg_edb <- bg %>% 
  mutate(age.range = if_else(age <= 34, '15-34', '35-64')) %>% 
  left_join(edb_age)

# write new dataframe to disk
# write_tsv(bg_edb, '../data/bg_edb.tsv')

