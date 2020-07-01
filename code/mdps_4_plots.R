
#-------------------------------------------------------------------------------
# Global
#-------------------------------------------------------------------------------

rm(list=ls())     # remove previous objects from workspace

library(tidyverse)
library(ggrepel)
library(ggmosaic)

# general settings for plots
theme_set(theme_bw() + theme(strip.background = element_blank()))
figdir <- "../plots/"
figscale <- 1
fw <- 5 # figure width
fh <- 5 # figure height

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
# Statistics per participant (byparticipant)
#-------------------------------------------------------------------------------

# For each participant: which proportion of the total produced DPs can be explained by each GAS?
# (1) analogical, (2) phon_lemma, (3) phon_noun
byparticipant <- d %>% 
  select(participant, analogical, phon_lemma, phon_noun) %>% 
  nest(data = c(analogical, phon_lemma, phon_noun)) %>% 
  mutate(n = map_int(data, nrow),
         freq_analogical = map_dbl(data, ~ sum(.x$analogical)),
         freq_phonlemma = map_dbl(data, ~ sum(.x$phon_lemma)),
         freq_phonnoun = map_dbl(data, ~ sum(.x$phon_noun)),
         data = NULL) %>% 
  mutate_at(vars(matches('freq_')), list(~ ./n)) %>% 
  left_join(bg_edb) %>% 
  mutate_at(vars(L1_S, S_area), factor)

#-------------------------------------------------------------------------------
# Figure 1: analogy by L1_S and area
#-------------------------------------------------------------------------------

# more readable L1_S labels
facetL1 <- c('0' = 'Spanish = L2',
             '1' = 'Spanish = L1')

#' ### Frequency of analogy ~ L1_S + % of S_speakers in population (add symbol for bilinguals)
set.seed(3)
byparticipant %>%
  ggplot() +
  aes(x = S_speakers*100, y = freq_analogical*100) +
  geom_label_repel(aes(fill = S_area, label = consecu.id.bil),
                   colour = "white", fontface = "bold", box.padding = 0) +
  scale_color_brewer(palette="Dark2") +
  scale_fill_brewer(palette="Dark2",
                    labels = c('B-dominant area', 'S-dominant area')) +
  facet_wrap(~ L1_S,
             labeller = labeller(L1_S = facetL1)) +
  geom_vline(xintercept = 50, linetype = 2, colour = 'black') +
  xlab("% of S monolinguals in participant's area") + ylab('% of analogical DPs') +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=.5),
        strip.background = element_blank(),
        legend.position = c(.12, .87),
        axis.title = element_text(face="bold", size = rel(1.1)),
        strip.text = element_text(face="bold", size = rel(1.1)),
        legend.title = element_text(face="bold", size = rel(1.1)),
        legend.text = element_text(size = rel(1.1)),
        axis.text = element_text(size = rel(1.2))) +
  guides(fill = guide_legend(title = 'Sociolinguistic area'))

# figfile <- paste0(figdir, 'analogy_area_biling.pdf')
# ggsave(figfile, scale = figscale*1.3, width = fw*1.5, height = fh, device=cairo_pdf)

#-------------------------------------------------------------------------------
# Figure 2: mosaic for det.gender by noun.end and L1_S
#-------------------------------------------------------------------------------

# rather empty theme for mosaic plot
theme_opts<-list(theme(panel.grid.minor = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.background = element_blank(),
                       panel.border = element_blank(),
                       plot.background = element_blank(),
                       axis.line = element_blank(),
                       axis.line.x = element_blank(),
                       axis.line.y = element_blank(),
                       axis.ticks = element_blank(),
                       plot.title = element_blank(),
                       legend.position = 'none',
                       strip.background = element_blank(),
                       axis.title = element_text(face="bold", size = rel(1.1)),
                       strip.text = element_text(face="bold", size = rel(1.1)),
                       axis.text = element_text(size = rel(1.2)),
                       axis.text.x = element_text(margin = margin(-15, 0, 0, 0)),
                       axis.text.y = element_text(margin = margin(0, -5, 0, 0)),
                       panel.spacing = unit(-.5, "lines")))

# manually define colours for data categories in mosaic.
# first F, then M, then gray
mosaic_palette <- c('#1b9e77', '#d95f02', '#d3d3d3')
# first F, then M
mosaic_palette <- c('#e41a1c', '#377eb8')
# highlight significant associations using alpha
mosaicalpha <- c(rep(.7, 8), rep(.4, 6), rep(.7, 2))

# base mosaic with no labels
g <- d %>% 
  mutate(noun.end = paste0('-', noun.end)) %>%
  ggplot() +
  geom_mosaic(aes(x = product(detGender, noun.end), fill = detGender),
              offset = .01, na.rm = T, alpha = mosaicalpha) +
  facet_grid(~ L1_S, labeller = labeller(L1_S = facetL1)) +
  scale_fill_manual(values=mosaic_palette) +
  theme_opts +
  labs(x="noun-final phoneme", y='determiner gender')

# get coordinates for mosaic squares (for correct label positioning)
temp <- ggplot_build(g)$data %>%
  as.data.frame %>%
  mutate(prop = as.character(round(ymax - ymin, 3)),
         x.position = (xmax + xmin) / 2,
         y.position = (ymax + ymin) / 2,
         L1_S = if_else(PANEL == 1, 0, 1)) %>% 
  separate(label, c('detGender', 'noun.end')) %>% 
  mutate_at(vars(L1_S, detGender, noun.end), factor)

# generate labels for each mosaic square
temp2 <- d %>% 
  count(L1_S, noun.end, detGender, .drop = F) %>%
  group_by(L1_S, noun.end) %>% 
  mutate(p = (n/sum(n)*100) %>% round(0)) %>%
  left_join(temp) %>% 
  mutate(labeltext = paste0(p, '%\nn=', n),
         labeltext = if_else((L1_S == 0) & (noun.end == 'i') & (detGender == 'f'),
                             NA_character_, labeltext))

# combine mosaic plot and descriptive labels
g + geom_label(data = temp2,
               aes(x = x.position, y = y.position, label = labeltext),
               inherit.aes = F,
               size = rel(2.5))

# figfile <- paste0(figdir, 'noun_mosaic.pdf')
# ggsave(figfile, scale = figscale*1.3, width = fw*1.5, height = fh, device=cairo_pdf)

