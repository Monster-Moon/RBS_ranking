rm(list = ls())
gc(reset = T)

setwd('/Users/moon/Documents/CBS_model/code/simulation_code/ex3')
source('../model.R')

load('ex3_data_list_500.Rdata')

model_list = vector('list', length(x_data_list))
for(i in 1:100)
{
  model_list[[i]] = btm_fit(x = x_data_list[[i]],
                            sim_n = 10,
                            eta = 0.05,
                            max_iter = 1000,
                            lambda1 = 0,
                            lambda2 = 0,
                            eps = 1e-4)
  cat(i, '\n')
}

save(model_list, file = 'ex3_btm_n10.Rdata')


for(i in 1:100)
{
  model_list[[i]] = btm_fit(x = x_data_list[[i]],
                            sim_n = 50,
                            eta = 0.05,
                            max_iter = 1000,
                            lambda1 = 0,
                            lambda2 = 0,
                            eps = 1e-4)
  cat(i, '\n')
}

save(model_list, file = 'ex3_btm_n50.Rdata')


for(i in 1:100)
{
  model_list[[i]] = btm_fit(x = x_data_list[[i]],
                            sim_n = 100,
                            eta = 0.05,
                            max_iter = 1000,
                            lambda1 = 0,
                            lambda2 = 0,
                            eps = 1e-4)
  cat(i, '\n')
}

save(model_list, file = 'ex3_btm_n100.Rdata')

