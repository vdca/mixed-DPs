
#-------------------------------------------------------------------------------
# Global
#-------------------------------------------------------------------------------

rm(list=ls())     # remove previous objects from workspace

library(tidyverse)

# statistical models
library(lmerTest)
library(blme)
library(broom)

# function takes a string like 'a/a' and returns a simplified version ('a')
# if the character before and after the slash are the same.
simplify_repetition <- function(x) if_else(str_detect(x, '(.)/\\1'), str_sub(x, 1, 1), x)

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
  mutate_at(vars(noun.end, root.end, detGender, L1_S, L1_B, S_area), factor) %>% 
  arrange(new_id, trial)

#-------------------------------------------------------------------------------
# Nested sets of data grouped by evidence: cue-based GAS
#-------------------------------------------------------------------------------

# for each participant, get:
# n.all = total number of mixed DPs produced
# d.GAS, n.GAS, p.GAS = all DPs **satisfying** the relevant GAS, with their n and proportion
# d.notGAS, n.notGAS, p.notGAS = all DPs **violating** the relevant GAS, with their n and proportion
d_nested <- d %>% 
  group_by(across(c(new_id, participant, role:S_area))) %>% 
  nest() %>% 
  ungroup() %>% 
  rename(d.all = data) %>% 
  mutate(n.all = map_int(d.all, nrow),
         # analog
         d.analog = map(d.all, ~ filter(.x, analogical==1)),
         n.analog = map_int(d.analog, nrow),
         p.analog = n.analog/n.all,
         # not-analog
         d.notanalog = map(d.all, ~ filter(.x, analogical==0)),
         n.notanalog = map_int(d.notanalog, nrow),
         p.notanalog = n.notanalog/n.all,
         # noun
         d.noun = map(d.all, ~ filter(.x, phon_noun==1)),
         n.noun = map_int(d.noun, nrow),
         p.noun = n.noun/n.all,
         # not-noun
         d.notnoun = map(d.all, ~ filter(.x, phon_noun==0)),
         n.notnoun = map_int(d.notnoun, nrow),
         p.notnoun = n.notnoun/n.all,
         # root
         d.root = map(d.all, ~ filter(.x, phon_lemma==1)),
         n.root = map_int(d.root, nrow),
         p.root = n.root/n.all,
         # not-root
         d.notroot = map(d.all, ~ filter(.x, phon_lemma==0)),
         n.notroot = map_int(d.notroot, nrow),
         p.notroot = n.notroot/n.all)

#-------------------------------------------------------------------------------
# Nested sets of data grouped by evidence: default-gender GAS
#-------------------------------------------------------------------------------

# for each participant, get:
# - DPs where noun has evidence for feminine gender (d.femevid)
# - n and proportion of these d.femevid where the produced det is actually feminine (femevid.f)
# - DPs with no evidence for feminine gender (d.nofemevid)
# - n and proportion of these d.nofemevid where the det is feminine (thus contradicting all 3 fem-cues)
# same for masculine gender.

d_defaultf <- d_nested %>% 
  select(new_id, participant, d.all, n.all) %>% 
  mutate(d.femevid = map(d.all, ~ filter(.x, (noun.end == 'a') |
                                           (root.end == 'a') |
                                           (es.gender == 'f'))),
         n.femevid = map_int(d.femevid, nrow),
         n.femevid.f = map_dbl(d.femevid, ~ sum(.x$default_f)),
         p.femevid.f = n.femevid.f/n.femevid,
         d.nofemevid = map(d.all, ~ filter(.x, (noun.end != 'a'),
                                           (root.end != 'a'),
                                           (es.gender != 'f'))),
         n.nofemevid = map_int(d.nofemevid, nrow),
         n.nofemevid.f = map_dbl(d.nofemevid, ~ sum(.x$default_f)),
         p.nofemevid.f = n.nofemevid.f/n.nofemevid) %>% 
  mutate(fem_where_evidence = paste0(n.femevid.f, '/', n.femevid),
         ev_freq = round(p.femevid.f, 2),
         fem_where_noevidence = paste0(n.nofemevid.f, '/', n.nofemevid),
         noev_freq = round(p.nofemevid.f, 2),
         p.explained = (n.femevid.f + n.nofemevid.f) / n.all,
         p.explained = round(p.explained, 2))

d_defaultm <- d_nested %>% 
  select(new_id, participant, d.all, n.all) %>% 
  mutate(d.mascevid = map(d.all, ~ filter(.x, (noun.end %in% c('o', 'i')) |
                                            (root.end %in% c('o', 'i')) |
                                            (es.gender == 'm'))),
         n.mascevid = map_int(d.mascevid, nrow),
         n.mascevid.m = map_dbl(d.mascevid, ~ sum(.x$default_m)),
         p.mascevid.m = n.mascevid.m/n.mascevid,
         d.nomascevid = map(d.all, ~ filter(.x, (!noun.end %in% c('o', 'i')),
                                            (!root.end %in% c('o', 'i')),
                                            (!es.gender == 'm'))),
         n.nomascevid = map_int(d.nomascevid, nrow),
         n.nomascevid.m = map_dbl(d.nomascevid, ~ sum(.x$default_m)),
         p.nomascevid.m = n.nomascevid.m/n.nomascevid) %>% 
  mutate(masc_where_evidence = paste0(n.mascevid.m, '/', n.mascevid),
         ev_freq = round(p.mascevid.m, 2),
         masc_where_noevidence = paste0(n.nomascevid.m, '/', n.nomascevid),
         noev_freq = round(p.nomascevid.m, 2),
         p.explained = (n.mascevid.m + n.nomascevid.m) / n.all,
         p.explained = round(p.explained, 2))

#-------------------------------------------------------------------------------
# Variance explained by each combined-evidence GAS candidate (only cue-based)
#-------------------------------------------------------------------------------

# for each participant:
# - which DPs are **not** explained by A(nalogy)?
# - out of the non-analogical: which are explained by N(oun.end)?
# - after considering the A and N GAS, which DPs are left unexplained?

# A>N (Analogical GAS > Noun.end GAS)
d_AN <- d_nested %>% 
  select(new_id, participant, contains('all'), contains('analog')) %>% 
  mutate(d.explained = map(d.notanalog, ~ filter(.x, phon_noun==1)),
         n.explained = map_int(d.explained, nrow),
         p.explained = n.explained/n.notanalog,
         p.AN = (n.analog+n.explained)/n.all,
         d.notexplained = map(d.notanalog, ~ filter(.x, phon_noun==0)))

# A>R (Analogical GAS > Root.end GAS)
d_AR <- d_nested %>% 
  select(new_id, participant, contains('all'), contains('analog')) %>% 
  mutate(d.explained = map(d.notanalog, ~ filter(.x, phon_lemma==1)),
         n.explained = map_int(d.explained, nrow),
         p.explained = n.explained/n.notanalog,
         p.AR = (n.analog+n.explained)/n.all,
         d.notexplained = map(d.notanalog, ~ filter(.x, phon_lemma==0)))

# N>R (Noun.end GAS > Root.end GAS)
d_NR <- d_nested %>% 
  select(new_id, participant, contains('all'), contains('noun')) %>% 
  mutate(d.explained = map(d.notnoun, ~ filter(.x, phon_lemma==1)),
         n.explained = map_int(d.explained, nrow),
         p.explained = n.explained/n.notnoun,
         p.NR = (n.noun+n.explained)/n.all,
         d.notexplained = map(d.notnoun, ~ filter(.x, phon_lemma==0)))

#-------------------------------------------------------------------------------
# Determine proportions explained by each GAS
#-------------------------------------------------------------------------------

# For each pair of combined-evidence GAS, determine whether
# a>b or b>a is preferred, for each participant.
gas_pref <- d_nested %>% 
  mutate(pref.AN = if_else(n.analog>n.noun, 'A', 'N'),
         pref.AN = if_else(n.analog==n.noun, 'A/N', pref.AN),
         pref.AR = if_else(n.analog>n.root, 'A', 'R'),
         pref.AR = if_else(n.analog==n.root, 'A/R', pref.AR),
         pref.NR = if_else(n.noun>n.root, 'N', 'R'),
         pref.NR = if_else(n.noun==n.root, 'N/R', pref.NR)) %>% 
  select(new_id, participant, starts_with('pref')) %>% 
  gather('candidate', 'GAS_1', starts_with('pref')) %>% 
  mutate(candidate = str_replace(candidate, 'pref.', '')) %>% 
  separate(candidate, c('candi1', 'candi2'), sep = 1, remove=F) %>% 
  mutate(GAS_2 = if_else(GAS_1==candi1, candi2, candi1),
         GAS_2 = if_else(str_detect(GAS_1, '/'), GAS_1, GAS_2),
         GAS = paste(GAS_1, GAS_2, sep='>')) %>% 
  select(-candi1, -candi2)

# join the data on how much each GAS (cue-based and defaults) explains for each participant
gas_summary <- select(d_AN, new_id, participant, p.AN) %>%
  left_join(select(d_AR, new_id, participant, p.AR)) %>% 
  left_join(select(d_NR, new_id, participant, p.NR)) %>%
  left_join(select(d_defaultf, new_id, participant, Fdef = p.explained)) %>%
  left_join(select(d_defaultm, new_id, participant, Mdef = p.explained)) %>% 
  gather('candidate', 'explained', p.AN:Mdef) %>% 
  mutate(candidate = str_replace(candidate, 'p.', '')) %>% 
  left_join(gas_pref) %>% 
  mutate(GAS = if_else(is.na(GAS), candidate, GAS)) %>%
  mutate(GAS_1 = if_else(is.na(GAS_1), candidate, GAS_1)) %>% 
  arrange(new_id, desc(explained))

#-------------------------------------------------------------------------------
# Get winner GAS and make tie situations explicit
#-------------------------------------------------------------------------------

# get the winner GAS for each participant (i.e. the GAS with highest explanatory power).
# calculate also whether there is a tie between two top-GAS (i.e. more than one winner)
gas_top <- gas_summary %>% 
  group_by(new_id) %>% 
  mutate(max_explained = max(explained),
         is_winner = explained == max_explained) %>% 
  filter(is_winner == T) %>%
  nest() %>% 
  mutate(n_winners = map_int(data, nrow)) %>% 
  unnest(data) %>% 
  ungroup() %>% 
  select(-is_winner, -max_explained) %>% 
  arrange(new_id, explained)

# 7 participants show a tie between two winner-GAS:
# for 2 participants (7, 21):
#       the same pair of strategies can be reversed (a>b, b>a);
#       it cannot be decided which of the two is primary or secondary.
# for 4 participants (1, 9, 22, 27):
#       the primary strategy is clear, but the secondary is not (a>b, a>c).
# for 1 participant (9): the primary is not clear, but the secondary is (a>c, b>c)
gas_top %>% 
  filter((n_winners > 1) |
         str_detect(GAS, '/')) %>% 
  group_by(new_id)

# re-encode winning GAS
gas_with_ties <- gas_top %>% 
  group_by(new_id, n_winners, explained) %>% 
  nest() %>% 
  mutate(topGAS_1 = map_chr(data, ~ paste(.$GAS_1, collapse = '/')),
         topGAS_2 = map_chr(data, ~ paste(.$GAS_2, collapse = '/')),
         topGAS_1 = simplify_repetition(topGAS_1),
         topGAS_2 = simplify_repetition(topGAS_2),
         topGAS_2 = str_replace_all(topGAS_2, 'NA', '--'),
         topGAS = paste(topGAS_1, topGAS_2, sep = '>')) %>% 
  select(-data) %>%
  ungroup()

# show participants with ties
gas_with_ties %>% 
  filter(str_detect(topGAS, '/'))

# merge background data
gas_bg <- left_join(gas_with_ties, bg_edb)

# subset participants who show no ties at all
no_ties <- gas_bg %>% 
  filter(!str_detect(topGAS, '/'))

#-------------------------------------------------------------------------------
# Summarise winner GAS (only clear cases, with no ties)
#-------------------------------------------------------------------------------

# filter out columns with tie symbol (= "/")
filter_tie <- function(dx, filter_var) dx %>% filter(!str_detect({{filter_var}}, '/'))

# calculate proportion in count data; apply filter_tie()
filter_prop <- function(dx, filter_var, ...) {
  dx %>% 
    filter_tie({{filter_var}}) %>% 
    group_by(...) %>% 
    mutate(total = sum(n),
           p = n/total*100) %>% 
    ungroup()
}

# winner GAS (no ties)
count(gas_bg, topGAS, sort = T) %>% filter_prop(topGAS)
count(gas_bg, topGAS_1, sort = T) %>% filter_prop(topGAS_1)
count(gas_bg, topGAS_2, sort = T) %>% filter_prop(topGAS_2)

# winner GAS sorted by L1_S (wide format).

# total n of participants with no tie at the combined GAS level is:
# 23: 11 with L1_S=0, 12 with L1_S=1.
# most frequent combined GAS for L1_S=0 is N>R (7/11)
# most frequent combined GAS for L1_S=1 is A>N and A>R (8/12; 4 each)
count(gas_bg, topGAS, L1_S, sort = T) %>% 
  filter_tie(topGAS) %>% 
  spread(L1_S, n)

# total n of participants with no tie at the primary GAS level is:
# 27: 11 with L1_S=0, 16 with L1_S=1.
# most frequent primary GAS for L1_S=0 is N  (8/11)
# most frequent primary GAS for L1_S=1 is A (11/16) -> 8/11 L1=S, 3/5 L1=SB
count(gas_bg, topGAS_1, L1_S, sort = T) %>% 
  filter_tie(topGAS_1) %>% 
  spread(L1_S, n)

count(gas_bg, topGAS_1, L1, sort = T) %>% 
  filter_tie(topGAS_1) %>% 
  spread(L1, n)

# total n of participants with no tie at the secondary GAS level is:
# 24: 12 with L1_S=0, 12 with L1_S=1.
# most frequent secondary GAS for L1_S=0 is R (7/12)
# most frequent secondary GAS for L1_S=1 is N (5/12)
count(gas_bg, topGAS_2, L1_S, sort = T) %>% 
  filter_tie(topGAS_2) %>% 
  spread(L1_S, n)

#-------------------------------------------------------------------------------
# merge winner-GAS information with explained DPs and exception DPs
#-------------------------------------------------------------------------------

winner_summary <- left_join(gas_top, gas_with_ties)

# gas_full:
# unify the explained and unexplained responses
# of both winner and non-winner GAS candidates
AN_full <- winner_summary %>% 
  left_join(d_AN) %>% 
  mutate(dataset = 'AN')
AR_full <- winner_summary %>% 
  left_join(d_AR) %>% 
  mutate(dataset = 'AR')
NR_full <- winner_summary %>% 
  left_join(d_NR) %>% 
  mutate(dataset = 'NR')

# merge all three datasets.
# *fix* the only participant with default_gender evidence has the winner GAS manually encoded
gas_full <- bind_rows(AN_full, AR_full, NR_full) %>% 
  mutate(winner = if_else(dataset == candidate, 1, 0),
         n.notexplained = map_int(d.notexplained, nrow),
         winner = if_else((participant=='B5') & (dataset=='AN'), 1, winner),
         dataset = if_else((participant=='B5') & (dataset=='AN'), 'defaultf', dataset),
         n.notexplained = if_else((participant=='B5') & (dataset=='AN'), as.integer(0), n.notexplained)) %>%
  arrange(desc(winner))

# subset only the winner strategies, and add participant bg data.
# gas_win dataframe contains 35 rows (while there are only 30 participants),
# that's because 5 participants show a tie between 2 different sets of candidate GAS (e.g. AN and AR)
gas_win <- gas_full %>%
  filter(winner == 1) %>% 
  left_join(bg_edb) %>% 
  arrange(new_id)

#-------------------------------------------------------------------------------
# helper functions to explore GAS
#-------------------------------------------------------------------------------

# default round() to digits = 2
rnd <- function(x, digits = 2) {round(x, digits)}

# get list of participants with a given GAS.
# by default, subset only participants with no ties; to include ties, change nwinners to 2
fun_participants <- function(gasx, nwinners = 1) {
  gas_win %>%
    filter(GAS == gasx,
           n_winners <= nwinners) %>% 
    select(new_id, participant, GAS, topGAS, p_explained=explained, n_unexplained=n.notexplained, L1, area) %>% 
    mutate(p_explained = rnd(p_explained)) %>% 
    arrange(-p_explained)
}

# get list of exceptions for a given winner GAS.
# by default, subset only participants with no ties; to include ties, change nwinners to 2
fun_exceptions <- function(gasx, nwinners = 1) {
  gas_win %>%
    filter(GAS == gasx,
           n_winners <= nwinners) %>% 
    select(new_id, participant, d.notexplained) %>% 
    unnest(d.notexplained) %>% 
    select(new_id, participant, det, noun, es, es.gender, en, fluidgen)
}

# mainly to show example DPs when looking for default GAS
fun_defaults <- function(x) {
  x %>%
    arrange(det) %>% 
    select(new_id, participant, det, noun, es, es.gender, en, fluidgen)
}

#-------------------------------------------------------------------------------
# Browse winner GAS and exceptions
#-------------------------------------------------------------------------------

# some examples:

fun_participants('A>N')

fun_participants('A>R')
fun_exceptions('A>R')

fun_participants('N>R')
fun_exceptions('N>R')

#-------------------------------------------------------------------------------
# Default feminine
#-------------------------------------------------------------------------------

# which participants show strongest evidence for default-F GAS?
d_defaultf %>%
  arrange(desc(p.explained)) %>% 
  select(new_id, participant, fem_where_evidence:p.explained)

# The participant with the strongest evidence for a default-f strategy is B5.
# In all cases (27/27) where there is some motivation to use a feminine det
# (i.e. lemma ends in -a, or produced noun ends in -a, or translation is f)
# B5 produces a feminine det.
# Besides, and more crucially, B5 also produces 2 mixed DPs where there is no
# independent motivation to produce a fem det and, still,
# B5 produces a fem det in both these cases:
# *la eraztun* ('el.MASC anillo', 'the ring'),
# *la soineko* ('el.MASC vestido', 'the dress')
d_defaultf %>%
  filter(participant == 'B5') %>% 
  unnest(d.nofemevid) %>% 
  fun_defaults()

# Compare these results to those of the second participant
# with the highest frequency of fem DET in cases
# where there is evidence for fem gender: A5.
# There are 5 exceptions to a putative default f.
# 3 (out of 29) items for which fem evidence is available are produced with a masc DET.
d_defaultf %>%
  filter(participant == 'A5') %>%
  unnest(d.femevid) %>% 
  filter(detGender == 'm') %>% 
  fun_defaults()

# A5 also produces 3 additional items with no fem evidence.
# 1/3 does conform to a putative default-f (la soineko);
# but the other two are produced with a masc DET.
d_defaultf %>%
  filter(participant == 'A5') %>%
  unnest(d.nofemevid) %>% 
  fun_defaults()

# These 5 exceptions could be explained as following the R (root-final-phoneme) GAS:
# since both the final phoneme /i/ (begi-a, eguzki-a)
# and /n/ (babarrun, kirten, eraztun)
# have been shown to have a masc preference.

# The mixed DPs produced by A5 are better explained by the N>A combined GAS,
# which explains 97% of the productions (i.e. only one exception):

gas_win %>% filter(participant == 'A5')

d_AN %>%
  filter(participant == 'A5') %>%
  unnest(d.notexplained) %>%
  fun_defaults()

#-------------------------------------------------------------------------------
# Default masculine
#-------------------------------------------------------------------------------

# In order to test whether a speaker uses a default gender (e.g. m),
# we select the cases where there is no independent motivation to use the gender in question,
# and check whether a putative default (m) is consistently used there.
# 
# The participant for which a default-m would explain more data is A9 (93%).

# which participants show strongest evidence for default-M GAS?
d_defaultm %>%
  arrange(desc(p.mascevid.m)) %>% 
  select(participant, masc_where_evidence:p.explained)

# In all cases (21/21) where there is some motivation to use a masculine det
# (i.e. lemma or noun ends in -o or -i, or translation is m)
# A9 produced a masculine det.
# A9 also produced 7 mixed DPs where there is no
# independent motivation to produce a masc det.

d_defaultm %>%
  filter(participant == 'A9') %>% 
  unnest(d.nomascevid) %>% 
  fun_defaults()

# In 71% of these cases (5/7), A9 produces a masc det (el: kipula, txanpon, babarrun, belaun, bekain).
# The two instances where A9 produces a fem det are *la kipula* and *la soka*.
# An ad-hoc explanation can be provided for these two exceptions.
# Out of a total of 28 mixed DPs,
# these two are the only ones for which A9 provided a self-correction.
# Namely, A9 produced *la kipula* as a rectification for *el kipula*
# ("el kipula, la kipula más bien mejor"),
# and A9 produced *la soka* as a rectification for a slip-of-the-tongue *la soga*,
# which is a non-mixed Spanish DP which includes the noun *soga*, a Spanish cognate
# with Basque *soka* ("la soga, la soka perdón, jajaja").
# 
# Side note: in Spanish there exists the word *soga* 'rope', and this is is a
# transparent cognate with Basque *soka* 'rope'. Still, we included this stimulus in our
# design because it is much more likely the image would elicit the Spanish word *cuerda* 'rope'
# (which is not cognate with Basque *soka*).
# Based on a picture naming experiment conducted with 100 speakers of Spanish
# (Duñabeitia et al. 2017 Multipic...), the same 'rope' picture we used was named
# *cuerda* by 91 participants, and *soga* by only 3 participants.
# 
# Alternatively, we can analyse A9's production as following an N>A combined GAS.
# 23/28 items can be explained by the N (noun-final-phoneme) GAS.
# 4 out of the 5 remaining items can be explained by analogical GAS.
# The single remaining un-explained item is *el kipula*,
# which is one of the two cases where A9 self-corrects
# ("el kipula, la kipula más bien mejor").

d %>% 
  filter(participant == 'A9',
         phon_noun == 0) %>%
  fun_defaults()

# Based on the previous data, it is unclear whether A9
# is following a default-m strategy or a combined N>A strategy.
# The experimental materials include 8 items which are likely to elicit
# Basque lemmas ending in -a.

