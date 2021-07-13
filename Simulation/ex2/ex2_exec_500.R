rm(list = ls())
gc(reset = T)

setwd('/Users/moon/Documents/CBS_model/code/simulation_code/ex2')
source('../model.R')

load('ex2_data_list_500.Rdata')

lambda_candidate = seq(0, 0.2, 0.05)

model_list = vector('list', length(x_data_list))
for(i in 1:100)
{
  tmp_model_list = vector('list', 5)
  for(j in 1:5)
  {
    tmp_model_list[[j]] = cbs_fit(x = x_data_list[[i]],
                                  sim_n = 10,
                                  eta = 0.05,
                                  max_iter = 1000,
                                  lambda1 = lambda_candidate[j],
                                  lambda2 = 0,
                                  eps = 1e-4)
  }
  model_list[[i]] = tmp_model_list
  cat(i, '\n')
}

save(model_list, file = 'ex2_500_n10.Rdata')

model_list = vector('list', length(x_data_list))
for(i in 1:100)
{
  tmp_model_list = vector('list', 5)
  for(j in 1:5)
  {
    tmp_model_list[[j]] = cbs_fit(x = x_data_list[[i]],
                                  sim_n = 50,
                                  eta = 0.05,
                                  max_iter = 1000,
                                  lambda1 = lambda_candidate[j],
                                  lambda2 = 0,
                                  eps = 1e-4)
  }
  model_list[[i]] = tmp_model_list
  cat(i, '\n')
}

save(model_list, file = 'ex2_500_n50.Rdata')


model_list = vector('list', length(x_data_list))
for(i in 1:100)
{
  tmp_model_list = vector('list', 5)
  for(j in 1:5)
  {
    tmp_model_list[[j]] = cbs_fit(x = x_data_list[[i]],
                                  sim_n = 100,
                                  eta = 0.05,
                                  max_iter = 1000,
                                  lambda1 = lambda_candidate[j],
                                  lambda2 = 0,
                                  eps = 1e-4)
  }
  model_list[[i]] = tmp_model_list
  cat(i, '\n')
}

save(model_list, file = 'ex2_500_n100.Rdata')





