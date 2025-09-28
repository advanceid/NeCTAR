#### Sample size calculation ####
# Preparation ----
rm(list = ls())
# setwd('D:/github')
setwd('NeCTAR/sampleSize')



# Library ----
# install.packages("brms", repos = "https://cloud.r-project.org")
library(brms)
library(doParallel)
library(magrittr)
library(tidyverse)



# Different scenarios ----
# N <- (1:20)*50
# p1 <- c(3:10)/20
# delta <- c(1:15)/100

N <- c(3:7)*100 + 50
p1 <- c(1:5)*0.1
delta <- c(2:6)/40
n_sampling <- 500   # number of sampling in estimation



# Function ----
recruitSim <- function(N = N[1], p1 = p1[1], delta = delta[1]){
  # assume 10 blocks, uneven distribution
  n_block <- 10
  comp <- runif(n_block, min = 0, max = 1)
  comp <- comp/sum(comp)
  block <- sample(1:n_block, N, replace = T, prob = comp)
  # block <- sample(1:10)
  # block <- rbinom(N, 10, 2:8/10) + 1
  block <- data.frame(block = block) 
  block %<>% arrange(block) %>%
    group_by(block) %>%
    summarise(count = n())
  
  recruit <- blockList %>% left_join(block, by = 'block') %>%
    filter(!is.na(count), sequence <= count) %>%
    mutate(probE = ifelse(allocation == 1,
                          p1 - delta,
                          p1)) %>%
    rowwise %>%
    mutate(outcome = rbinom(1, 1, probE))
  return(recruit)
}

func <- function(index = index, scenario = scenario){
  for(i in 1:nrow(scenario)){
    dfSim <- recruitSim(N = scenario$N[i],
                        p1 = scenario$p1[i],
                        delta = scenario$delta[i])
    
    fit <- brm(data = dfSim,
               outcome ~ allocation,
               family  = bernoulli(link = "logit"),
               prior = prior(normal(0, 10), class = b),
               chains = 4,
               cores = 4)

    coef <- fixef(fit)
    scenario[i, c('est', 's.e.', 'QL', 'QU')] <- coef[2,]
    scenario$index <- index
  }
  return(scenario)
}



# Load block randomization list ----
blockList <- read.csv('blockList.csv') %>%
  rename(block = pattern) %>%
  mutate(allocation = recode(allocation, 
                             'colistin' = 1,
                             'saline' = 0))



# Loop ----
scenario <- expand.grid(N = N, p1 = p1, delta = delta)
scenario <- scenario %>% filter(p1 == 0.4)
index <- split(1:n_sampling, 1:n_sampling)

cl <- makeCluster(50)
clusterExport(cl, c('recruitSim', 'blockList'))
clusterEvalQ(cl, {
  library(dplyr)
  library(magrittr)
  library(brms)})
results <- parLapply(cl, index, func, scenario = scenario)
stopCluster(cl)

combined <- do.call(rbind, results)
write.csv(combined, 'result/repeat500.csv', row.names = F)

q()



# Result summary ----
combined <- read.csv('result/repeat1.csv')



# Summarise ----
stat <- combined %>% group_by(N, p1, delta) %>%
  summarize(count = n(), mean = mean(est), 
            pass = sum(QU < 0)) %>%
  mutate(beta = pass/count,
         group = as.character(p1),
         delta = as.character(delta))



# plot ----
pic <- stat %>% #filter(delta == 0.1) %>%
  ggplot(data = ., aes(x = N, y = beta, group = p1)) +
  # geom_line(aes(col = group))
  geom_smooth(aes(col = group)) + 
  facet_wrap(vars(delta), nrow = 2)
pic
  



  