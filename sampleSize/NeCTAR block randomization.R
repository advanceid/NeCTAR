rm(list = ls())
library(tidyverse)
library(magrittr)
#### test (block randomization with varying block size) ####
sampleSize <- 520           # pre-specify sample size
blockSizeList <- c(2, 4, 6) # pre-specify possible block sizes
pattern <- 1:20             # assign a code to different patterns
treatment <- c('colistin', 'saline')  # treatment groups

patternTreatment <- expand.grid(  
  # List out all combinations of pattern, treatment and block sizes
  # Prepare more than enough blocks to randomize, eradicate redundant ones later
  pattern = pattern, treatment = treatment, 
  blocki = blockSizeList, blockID = 1:sampleSize)

patternTreatment %<>% 
  arrange(pattern, blockID) %>% group_by(pattern, blockID) %>% 
  mutate(blockSize = sample(blocki, replace = T, prob = 1/blocki),  
         blockSize = first(blockSize)) %>%
  filter(blocki <= blockSize) %>% 
  # randomize block size per block, 
  # set prob so patients have equal chance to be allocated to block of the 3 sizes
  # drop records above block size, before randomize treatment
  group_by(pattern) %>% mutate(sequence = row_number()) %>%
  # each record within a pattern has a unique sequence
  group_by(pattern, blockID) %>%
  mutate(allocation = sample(treatment)) %>%  
  # randomize treatment within blocks
  arrange(pattern, blockID, sequence) %>% 
  dplyr::select(pattern, sequence, blockID, blockSize, allocation)
# write.csv(patternTreatment, 'blockList.csv', row.names = F)


#### Check ####
# plus-minus versus number of subj by pattern
check1 <- patternTreatment %>% group_by(pattern) %>%
  mutate(pos = ifelse(allocation == 'colistin', 1, -1),
    cum = cumsum(pos), group = paste0('pattern', pattern)) %>%
  group_by(pattern, blockID) %>% mutate(
    i = row_number(), cum2 = abs(nth(cum, 2)), cum4 = abs(nth(cum, 4)))
check1 %<>% group_by(pattern, blockID) %>%
  mutate(maxCum = max(abs(cum)), cum24 = (cum2 == 2)&(cum4 == 2), cum24 = max(cum24)) %>%
  mutate(guess1 = (maxCum == 3)&(i %in% 4:6),
         guess2 = cum24&(i %in% 5:6),
         guess = guess1|guess2)
guessRate <- check1 %>% group_by(pattern) %>%
  summarise(rate = sum(guess)/n())
for(i in 1:nrow(guessRate)){
  print(
    paste0('For pattern ', guessRate$pattern[i], 
           ', allocation of ', round(100*guessRate$rate[i], 1), 
           ' in every 100 patients can be determined early')
  )
}
100*mean(guessRate$rate)

# distribution of block size
table(patternTreatment$blockSize)
table(patternTreatment$blockSize)[1]/blockSizeList[1]
table(patternTreatment$blockSize)[2]/blockSizeList[2]
table(patternTreatment$blockSize)[3]/blockSizeList[3]

# even distribution of treatment combination (in sequence)
# block size = 2
check2a <- patternTreatment %>% filter(blockSize == 2) %>%
  group_by(pattern, blockID) %>%
  mutate(rank = row_number(), rank = paste0('comb', rank)) %>%
  pivot_wider(values_from = 'allocation', names_from = rank,
              id_cols = c('pattern', 'blockID')) %>%
  group_by(comb1, comb2) %>% summarise(n = n())
hist(check2a$n)

# block size = 4
check2b <- patternTreatment %>% filter(blockSize == 4) %>%
  group_by(pattern, blockID) %>%
  mutate(rank = row_number(), rank = paste0('comb', rank)) %>%
  pivot_wider(values_from = 'allocation', names_from = rank,
              id_cols = c('pattern', 'blockID')) %>%
  group_by(comb1, comb2, comb3, comb4) %>% summarise(n = n())
hist(check2b$n)

check2c <- patternTreatment %>% filter(blockSize == 6) %>%
  group_by(pattern, blockID) %>%
  mutate(rank = row_number(), rank = paste0('comb', rank)) %>%
  pivot_wider(values_from = 'allocation', names_from = rank,
              id_cols = c('pattern', 'blockID')) %>%
  group_by(comb1, comb2, comb3, comb4, comb5, comb6) %>% summarise(n = n())
hist(check2c$n)

# treatment groups appear at equal sequence
check3 <- patternTreatment %>% group_by(pattern, allocation) %>%
  summarise(mean = mean(sequence), 
            sd = sd(sequence)) %>%
  pivot_wider(id_cols = pattern, names_from = 'allocation', values_from = c('mean', 'sd')) %>%
  mutate(meanDif = mean_colistin - mean_saline,
         sdDif = sd_colistin - sd_saline)
hist(check3$meanDif)
hist(check3$sdDif)