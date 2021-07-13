rm(list = ls())
gc(reset = T)

setwd('/Users/moon/Documents/CBS_model/code/real_data/soc_code')

require(BradleyTerry2)
require(dplyr)

#### soc mtpuz ####
source('..//model.R')

readLines('../soc/ED-00025-00000004.soc', n = 6)
mtpuz_data = read.table('../soc/ED-00025-00000004.soc', skip = 6, sep = ',')

head(mtpuz_data)
mtpuz_data_list = lapply(1:nrow(mtpuz_data), function(i)
  mtpuz_data[rep(i, mtpuz_data[i, 1]), -1])

mtpuz_data = mtpuz_data_list %>% do.call('rbind', .) %>% as.matrix()
rownames(mtpuz_data) = colnames(mtpuz_data) = NULL
head(mtpuz_data)

barplot(table(apply(mtpuz_data, 1, paste0, collapse = '')))


#### data generates ####
data_list = vector('list', 100)
for(i in 1:100)
{
  set.seed(i)
  split_sample = sample(1:3, size = nrow(mtpuz_data), prob=c(0.7,0.15,0.15), replace = TRUE)
  train = mtpuz_data[split_sample == 1, ] %>% as.matrix()
  valid = mtpuz_data[split_sample == 2,] %>% as.matrix()
  test =  mtpuz_data[split_sample == 3,] %>% as.matrix()
  data_list[[i]] = list(train = train, valid = valid, test = test)
  cat(i, '\n')
}
save(data_list, file = 'r2_data_list.Rdata')

#### model fitting ####
setwd('/Users/moon/Documents/CBS_model/code/real_data/soc_code')
load('r2_data_list.Rdata')
lambda_candidate = c(0, 0.3, 0.5, 0.7, 1)

model_list = vector('list', 100)
for(i in 1:100)
{
  tmp_model_list = vector('list', length(lambda_candidate))
  for(j in 1:length(lambda_candidate))
  {
    tmp_fit = cbs_fit(x = data_list[[i]]$train, sim_n = 100, eta = 0.005, 
                      max_iter = 2000, lambda1 = lambda_candidate[j], lambda2 = 0)    
    tmp_model_list[[j]] = tmp_fit
  }
  model_list[[i]] = tmp_model_list
  cat(i, '\n')
}

save(model_list, file = 'r2_model_list.Rdata')

#### btm model fitting ####
source('..//model.R')
load('r2_data_list.Rdata')

n_items = ncol(data_list[[1]]$train)
n_items_comb = n_items * (n_items - 1) / 2

btm_model_list = vector('list', 100)
for(i in 1:100)
{
  btm_model_list[[i]] = btm_fit(x = data_list[[i]]$train,
                                sim_n = 100,
                                eta = 0.005,
                                max_iter = 2000,
                                lambda1 = 0,
                                lambda2 = 0,
                                eps = 1e-4)
  cat(i, '\n') 
}
save(btm_model_list, file = 'r2_btm_list.Rdata')




#### bt model fitting #####
source('..//model.R')
load('r2_data_list.Rdata')

n_items = ncol(data_list[[1]]$train)
n_items_comb = n_items * (n_items - 1) / 2

A_mat = A_mat_fun(n_items)
tmp_bt_data = expand.grid(1:n_items, 1:n_items)
tmp_bt_data = tmp_bt_data[tmp_bt_data[, 1] < tmp_bt_data[, 2], ] %>% arrange(Var1)

bt_model_list = vector('list', 100)
for(i in 1:100)
{
  bt_data = cbind(tmp_bt_data, t(apply(x_data_generates(data_list[[i]]$train), 2, table))[,2:1])
  bt_data$Var1 = factor(bt_data$Var1, levels = 1:n_items)
  bt_data$Var2 = factor(bt_data$Var2, levels = 1:n_items)
  colnames(bt_data) = c('home_team', 'away_team', 'home_win', 'away_win')
  bt_model_list[[i]] = BTm(cbind(home_win, away_win), home_team, away_team, data = bt_data, id = "team")
  cat(i, '\n') 
}
save(bt_model_list, file = 'r2_bt_list.Rdata')

#### configure ####
rm(list = ls())
gc(reset = T)
setwd('/Users/moon/Documents/CBS_model/code/real_data/soc_code')
source('..//model.R')
load('r2_data_list.Rdata')
load('r2_model_list.Rdata')
load('r2_btm_list.Rdata')
load('r2_bt_list.Rdata')

n_items = ncol(data_list[[1]]$train)
n_items_comb = n_items * (n_items - 1) / 2
permu_items = permutations(n_items)
all_data = x_data_generates(permu_items)


train_kl_mat = train_tv_mat = train_hd_mat = matrix(0, nrow = 100, ncol = 5)
for(i in 1:100)
{
  train_g = x_data_generates(data_list[[i]]$train)
  train_g_p = apply(train_g, 1, paste0, collapse = '')
  train_order = order(train_g_p)
  
  train_g_p_df = data.frame(table(train_g_p))
  train_g_p_df$prop = train_g_p_df$Freq / sum(train_g_p_df$Freq)
  
  for(j in 1:5)
  {
    prob_est = exp(all_data %*% model_list[[i]][[j]]$theta_hat) / 
      sum(exp(all_data %*% model_list[[i]][[j]]$theta_hat))
    prob_est_data = data.frame(train_g_p = apply(all_data, 1, paste0, collapse = ''), prob_est)
    
    tmp_train_g_p_df = train_g_p_df %>% left_join(prob_est_data, by = 'train_g_p')
    
    train_kl_mat[i, j] = -sum(tmp_train_g_p_df$prop * log(tmp_train_g_p_df$prob_est))
    train_tv_mat[i, j] = sum(abs(tmp_train_g_p_df$prop - tmp_train_g_p_df$prob_est))
    train_hd_mat[i, j] = sum((sqrt(tmp_train_g_p_df$prop) - sqrt(tmp_train_g_p_df$prob_est))^2)
  }
}

valid_kl_mat = valid_tv_mat = valid_hd_mat = matrix(0, nrow = 100, ncol = 5)
for(i in 1:100)
{
  valid_g = x_data_generates(data_list[[i]]$valid)
  valid_g_p = apply(valid_g, 1, paste0, collapse = '')
  valid_order = order(valid_g_p)
  
  valid_g_p_df = data.frame(table(valid_g_p))
  valid_g_p_df$prop = valid_g_p_df$Freq / sum(valid_g_p_df$Freq)
  
  for(j in 1:5)
  {
    prob_est = exp(all_data %*% model_list[[i]][[j]]$theta_hat) / 
      sum(exp(all_data %*% model_list[[i]][[j]]$theta_hat))
    prob_est_data = data.frame(valid_g_p = apply(all_data, 1, paste0, collapse = ''), prob_est)
    
    tmp_valid_g_p_df = valid_g_p_df %>% left_join(prob_est_data, by = 'valid_g_p')
    
    valid_kl_mat[i, j] = -sum(tmp_valid_g_p_df$prop * log(tmp_valid_g_p_df$prob_est))
    valid_tv_mat[i, j] = sum(abs(tmp_valid_g_p_df$prop - tmp_valid_g_p_df$prob_est))
    valid_hd_mat[i, j] = sum((sqrt(tmp_valid_g_p_df$prop) - sqrt(tmp_valid_g_p_df$prob_est))^2)
  }
}

test_kl_mat = test_tv_mat = test_hd_mat = matrix(0, nrow = 100, ncol = 5)
for(i in 1:100)
{
  test_g = x_data_generates(data_list[[i]]$test)
  test_g_p = apply(test_g, 1, paste0, collapse = '')
  test_order = order(test_g_p)
  
  test_g_p_df = data.frame(table(test_g_p))
  test_g_p_df$prop = test_g_p_df$Freq / sum(test_g_p_df$Freq)
  
  for(j in 1:5)
  {
    prob_est = exp(all_data %*% model_list[[i]][[j]]$theta_hat) / 
      sum(exp(all_data %*% model_list[[i]][[j]]$theta_hat))
    prob_est_data = data.frame(test_g_p = apply(all_data, 1, paste0, collapse = ''), prob_est)
    
    tmp_test_g_p_df = test_g_p_df %>% left_join(prob_est_data, by = 'test_g_p')
    
    test_kl_mat[i, j] = -sum(tmp_test_g_p_df$prop * log(tmp_test_g_p_df$prob_est))
    test_tv_mat[i, j] = sum(abs(tmp_test_g_p_df$prop - tmp_test_g_p_df$prob_est))
    test_hd_mat[i, j] = sum((sqrt(tmp_test_g_p_df$prop) - sqrt(tmp_test_g_p_df$prob_est))^2)
  }
}

A_mat = A_mat_fun(n_items)
train_bt_kl_vec = train_bt_tv_vec = train_bt_hd_vec = numeric(100)
for(i in 1:100)
{
  train_g = x_data_generates(data_list[[i]]$train)
  train_g_p = apply(train_g, 1, paste0, collapse = '')
  train_order = order(train_g_p)
  
  train_g_p_df = data.frame(table(train_g_p))
  train_g_p_df$prop = train_g_p_df$Freq / sum(train_g_p_df$Freq)
  
  tmp_theta_bt_vec = A_mat[1:n_items_comb, 1:n_items] %*% c(0, bt_model_list[[i]]$coefficients)
  prob_est_bt = exp(all_data %*% tmp_theta_bt_vec) / sum(exp(all_data %*% tmp_theta_bt_vec))
  
  prob_est_bt_data = data.frame(train_g_p = apply(all_data, 1, paste0, collapse = ''), prob_est_bt)
  train_g_p_df = train_g_p_df %>% left_join(prob_est_bt_data, by = 'train_g_p')
  
  train_bt_kl_vec[i] = -sum(train_g_p_df$prop * log(train_g_p_df$prob_est_bt))
  train_bt_tv_vec[i] = sum(abs(train_g_p_df$prop - train_g_p_df$prob_est_bt))
  train_bt_hd_vec[i] = sum((sqrt(train_g_p_df$prop) - sqrt(train_g_p_df$prob_est_bt))^2)
}

test_bt_kl_vec = test_bt_tv_vec = test_bt_hd_vec = numeric(100)
for(i in 1:100)
{
  test_g = x_data_generates(data_list[[i]]$test)
  test_g_p = apply(test_g, 1, paste0, collapse = '')
  test_order = order(test_g_p)
  
  test_g_p_df = data.frame(table(test_g_p))
  test_g_p_df$prop = test_g_p_df$Freq / sum(test_g_p_df$Freq)
  
  tmp_theta_bt_vec = A_mat[1:n_items_comb, 1:n_items] %*% c(0, bt_model_list[[i]]$coefficients)
  prob_est_bt = exp(all_data %*% tmp_theta_bt_vec) / sum(exp(all_data %*% tmp_theta_bt_vec))
  
  prob_est_bt_data = data.frame(test_g_p = apply(all_data, 1, paste0, collapse = ''), prob_est_bt)
  test_g_p_df = test_g_p_df %>% left_join(prob_est_bt_data, by = 'test_g_p')
  
  test_bt_kl_vec[i] = -sum(test_g_p_df$prop * log(test_g_p_df$prob_est_bt))
  test_bt_tv_vec[i] = sum(abs(test_g_p_df$prop - test_g_p_df$prob_est_bt))
  test_bt_hd_vec[i] = sum((sqrt(test_g_p_df$prop) - sqrt(test_g_p_df$prob_est_bt))^2)
}

train_btm_kl_vec = train_btm_tv_vec = train_btm_hd_vec = numeric(100)
for(i in 1:100)
{
  train_g = x_data_generates(data_list[[i]]$train)
  train_g_p = apply(train_g, 1, paste0, collapse = '')
  train_order = order(train_g_p)
  
  train_g_p_df = data.frame(table(train_g_p))
  train_g_p_df$prop = train_g_p_df$Freq / sum(train_g_p_df$Freq)
  
  prob_est_btm = exp(all_data %*% btm_model_list[[i]]$theta_hat) / 
    sum(exp(all_data %*% btm_model_list[[i]]$theta_hat))
  prob_est_btm_data = data.frame(train_g_p = apply(all_data, 1, paste0, collapse = ''), prob_est_btm)
  train_g_p_df = train_g_p_df %>% left_join(prob_est_btm_data, by = 'train_g_p')
  
  train_btm_kl_vec[i] = -sum(train_g_p_df$prop * log(train_g_p_df$prob_est_btm))
  train_btm_tv_vec[i] = sum(abs(train_g_p_df$prop - train_g_p_df$prob_est_btm))
  train_btm_hd_vec[i] = sum((sqrt(train_g_p_df$prop) - sqrt(train_g_p_df$prob_est_btm))^2)
}

test_btm_kl_vec = test_btm_tv_vec = test_btm_hd_vec = numeric(100)
for(i in 1:100)
{
  test_g = x_data_generates(data_list[[i]]$test)
  test_g_p = apply(test_g, 1, paste0, collapse = '')
  test_order = order(test_g_p)
  
  test_g_p_df = data.frame(table(test_g_p))
  test_g_p_df$prop = test_g_p_df$Freq / sum(test_g_p_df$Freq)
  
  prob_est_btm = exp(all_data %*% btm_model_list[[i]]$theta_hat) / 
    sum(exp(all_data %*% btm_model_list[[i]]$theta_hat))
  prob_est_btm_data = data.frame(test_g_p = apply(all_data, 1, paste0, collapse = ''), prob_est_btm)
  test_g_p_df = test_g_p_df %>% left_join(prob_est_btm_data, by = 'test_g_p')
  
  test_btm_kl_vec[i] = -sum(test_g_p_df$prop * log(test_g_p_df$prob_est_btm))
  test_btm_tv_vec[i] = sum(abs(test_g_p_df$prop - test_g_p_df$prob_est_btm))
  test_btm_hd_vec[i] = sum((sqrt(test_g_p_df$prop) - sqrt(test_g_p_df$prob_est_btm))^2)
}

train = list(train_hd_mat = train_hd_mat,
             train_kl_mat = train_kl_mat,
             train_tv_mat = train_tv_mat,
             train_btm_tv_vec = train_btm_tv_vec,
             train_btm_kl_vec = train_btm_kl_vec,
             train_btm_hd_vec = train_btm_hd_vec,
             train_bt_tv_vec = train_bt_tv_vec,
             train_bt_kl_vec = train_bt_kl_vec,
             train_bt_hd_vec = train_bt_hd_vec)

test = list(test_hd_mat = test_hd_mat,
            test_kl_mat = test_kl_mat,
            test_tv_mat = test_tv_mat,
            test_btm_tv_vec = test_btm_tv_vec,
            test_btm_kl_vec = test_btm_kl_vec,
            test_btm_hd_vec = test_btm_hd_vec,
            test_bt_tv_vec = test_bt_tv_vec,
            test_bt_kl_vec = test_bt_kl_vec,
            test_bt_hd_vec = test_bt_hd_vec)

valid = list(valid_hd_mat = valid_hd_mat,
             valid_kl_mat = valid_kl_mat,
             valid_tv_mat = valid_tv_mat)

r2_result = list(train = train, test = test, valid = valid)
save(r2_result, file = 'r2_result.Rdata')




