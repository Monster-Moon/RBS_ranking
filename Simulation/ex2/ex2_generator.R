rm(list = ls())
gc(reset = T)

setwd('D:\\Paper\\cbs_model\\code\\simulation_code')
source('model.R')

n_items = 5
n_items_comb = n_items * (n_items - 1) / 2

A_mat = A_mat_fun(n_items = n_items) ## design matrix

gamma_star = c(0, -0.1, 0.05, -0.05, 0.1, 0.05, 0.05, -0.1, -0.1, 0.1)
gamma_ratio = 1
beta_star = c(-1, -0.8, 0.5, 0.7, 0.6, gamma_star * gamma_ratio)

theta_star = A_mat %*% beta_star
theta_star = theta_star[1:n_items_comb]

test_rank_data = permutations(n_items)
test_data = x_data_generates(data = test_rank_data)

n = 500L
sample_prob = exp(test_data %*% theta_star) / sum(exp(test_data %*% theta_star))

x_data_list = vector('list', 100)
for(i in 1:100)
{
  sample_inx = sample(1:nrow(test_data), size = n, replace = T, prob = sample_prob)
  sample_inx = factor(sample_inx, levels = 1:nrow(test_data))
  sample_inx_prop = table(sample_inx) / sum(table(sample_inx))
  
  x_data_list[[i]] = test_rank_data[sample_inx, ]
  cat(i, '\n')
}

save(x_data_list, file = 'ex2_data_list_500.Rdata')


