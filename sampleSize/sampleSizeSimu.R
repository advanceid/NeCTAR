#### Sample size calculation ####
# Preparation ----
rm(list = ls())
# setwd('D:/NUS Dropbox/Xiangyuan Huang/github/')
setwd('NeCTAR/sampleSize')



# Library ----
# install.packages("brms", repos = "https://cloud.r-project.org")
library(brms)
library(doParallel)
library(magrittr)
library(tidyverse)



# Different scenarios ----
N <- 400 + (1:8)*25
p1 <- c(0.3, 0.35, 0.4)
delta <- c(0.05, 0.06, 0.08, 0.1)
n_sampling <- 1000   # number of sampling in estimation



# Function ----
recruitSim <- function(N = N[1], p1 = p1[1], delta = delta[1]){
  # assume 10 blocks, uneven distribution
  n_block <- 10
  comp <- runif(n_block, min = 0, max = 1)
  comp <- comp/sum(comp)
  block <- sample(1:n_block, N, replace = T, prob = comp)
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
    
    # fit <- brm(data = dfSim,
    #            outcome ~ allocation,
    #            family  = bernoulli(link = "logit"),
    #            prior = prior(normal(0, 10), class = b),
    #            chains = 4, cores = 4) 
    #            # file = "model_cache", save_model = "model.stan")
    fit <- update(fit0, newdata = dfSim)
    rm(list = 'dfSim')

    coef <- fixef(fit, probs = c(0.025, 0.9, 0.95, 0.975))
    scenario[i, c('est', 's.e.', 'q025', 'q900', 'q950', 'q975')] <- coef[2,]
    scenario$index <- index
    rm(list = c('fit'))
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
# scenario <- scenario[1:3, ]
index <- split(1:n_sampling, 1:n_sampling)

dfSim <- data.frame(outcome = c(0,1), allocation = c(0,1))
fit0 <- brm(data = dfSim,
           outcome ~ allocation,
           family  = bernoulli(link = "logit"),
           prior = prior(normal(0, 10), class = b),
           chains = 0) 

cl <- makeCluster(10)
clusterExport(cl, c('recruitSim', 'blockList', 'fit0'))
clusterEvalQ(cl, {
  library(dplyr)
  library(magrittr)
  library(brms)})

objs <- ls()
sapply(objs, function(x) object.size(get(x)))

results <- parLapply(cl, index, func, scenario = scenario)

objs <- ls()
sapply(objs, function(x) object.size(get(x)))

stopCluster(cl)

combined <- do.call(rbind, results)
write.csv(combined, 'result/repeat.csv', row.names = F)

q()



# Result summary ----
result1 <- read.csv('result/repeat100a.csv')
result2 <- read.csv('result/repeat100b.csv') %>%
  mutate(index = index + 100)
# result <- read.csv('result/repeat.csv')

result1 <- read.csv('result/repeat0930.csv')
result2 <- read.csv('result/repeat.csv')

result <- rbind(result1, result2)



# Summarise ----
stat <- result %>% group_by(N, p1, delta) %>%
  summarize(count = n(), mean = mean(est), 
            pass = sum(q950 < 0)) %>%
  mutate(beta = pass/count,
         group = as.character(p1),
         delta = as.character(delta))



# plot ----
pic <- stat %>% #filter(delta == 0.1) %>%
  ggplot(data = ., aes(x = N, y = beta, group = p1)) +
  geom_point() +
  # geom_line(aes(col = group))
  geom_smooth(aes(col = group), se = F) + 
  facet_wrap(vars(delta), nrow = 2)
pic

ggsave('power.tiff', dpi = 300, width = 12, height = 8)
  



  