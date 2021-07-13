rm(list = ls())
gc(reset = T)

if(!require(dplyr)) install.packages('dplyr'); require(dplyr)

setwd('/Users/moon/Documents/CBS_model/code/simulation_code/ex1')
source('../model.R')

n_items = 5
n_items_comb = n_items * (n_items - 1) / 2

A_mat = A_mat_fun(n_items = n_items) ## design matrix

gamma_star = c(0, -0.2, 0, 0, 0.2, 0, 0, 0, 0, 0)
gamma_ratio = 1
beta_star = c(0.7, 0.5, 0, -0.4, -0.8, gamma_star * gamma_ratio)

theta_star = A_mat %*% beta_star
theta_star = theta_star[1:n_items_comb]

permu_items = permutations(5)
all_data = x_data_generates(data = permu_items)

prob_star = exp(all_data %*% theta_star) / sum(exp(all_data %*% theta_star))

########### BT #############
load('ex1_500_bt.Rdata')

bt_theta_list = lapply(bt_model_list, 
                       function(x) A_mat[1:n_items_comb, 1:n_items] %*% c(0, x$coefficients))
prob_bt_est = lapply(bt_theta_list, function(x)
  exp(all_data %*% x) / sum(exp(all_data %*% x)))

kl_bt_vec = lapply(prob_bt_est, function(x) -sum(prob_star * log(x))) %>% unlist()
tv_bt_vec = lapply(prob_bt_est, function(x) sum(abs(x - prob_star))) %>% unlist()
hd_bt_vec = lapply(prob_bt_est, function(x) sum((sqrt(prob_star) - sqrt(x))^2)) %>% unlist()
mean(kl_bt_vec)
mean(tv_bt_vec)
mean(hd_bt_vec)

######### BTM ##############
load('ex1_btm_n10.Rdata')

prob_btm_est = lapply(model_list, function(x)
  exp(all_data %*% x$theta_hat) / sum(exp(all_data %*% x$theta_hat)))

kl_btm_10_vec = lapply(prob_btm_est, function(x) -sum(prob_star * log(x))) %>% unlist()
tv_btm_10_vec = lapply(prob_btm_est, function(x) sum(abs(x - prob_star))) %>% unlist()
hd_btm_10_vec = lapply(prob_btm_est, function(x) sum((sqrt(prob_star) - sqrt(x))^2)) %>% unlist()

load('ex1_btm_n50.Rdata')

prob_btm_est = lapply(model_list, function(x)
  exp(all_data %*% x$theta_hat) / sum(exp(all_data %*% x$theta_hat)))

kl_btm_50_vec = lapply(prob_btm_est, function(x) -sum(prob_star * log(x))) %>% unlist()
tv_btm_50_vec = lapply(prob_btm_est, function(x) sum(abs(x - prob_star))) %>% unlist()
hd_btm_50_vec = lapply(prob_btm_est, function(x) sum((sqrt(prob_star) - sqrt(x))^2)) %>% unlist()

load('ex1_btm_n100.Rdata')

prob_btm_est = lapply(model_list, function(x)
  exp(all_data %*% x$theta_hat) / sum(exp(all_data %*% x$theta_hat)))

kl_btm_100_vec = lapply(prob_btm_est, function(x) -sum(prob_star * log(x))) %>% unlist()
tv_btm_100_vec = lapply(prob_btm_est, function(x) sum(abs(x - prob_star))) %>% unlist()
hd_btm_100_vec = lapply(prob_btm_est, function(x) sum((sqrt(prob_star) - sqrt(x))^2)) %>% unlist()



###### CBS 10 ##########
load('ex1_500_n10.Rdata')

prob_est = lapply(model_list, function(y)
  lapply(y, function(x) exp(all_data %*% x$theta_hat) / sum(exp(all_data %*% x$theta_hat))))

kl_list = lapply(prob_est, function(y)
  lapply(y, function(x) -sum(prob_star * log(x))))

kl_mat = lapply(kl_list, unlist) %>% do.call('rbind', .)
boxplot(kl_mat)

tv_list = lapply(prob_est, function(y)
  lapply(y, function(x) sum(abs(x - prob_star))))

tv_mat = lapply(tv_list, unlist) %>% do.call('rbind', .)

hd_list = lapply(prob_est, function(y)
  lapply(y, function(x) sum((sqrt(prob_star) - sqrt(x))^2)))

hd_mat = lapply(hd_list, unlist) %>% do.call('rbind', .)

kl_return_mat_10 = cbind(kl_mat[, 1], kl_mat[, 5])
tv_return_mat_10 = cbind(tv_mat[, 1], tv_mat[, 5])
hd_return_mat_10 = cbind(hd_mat[, 1], hd_mat[, 5])

colMeans(kl_return_mat_10)
colMeans(tv_return_mat_10)
colMeans(hd_return_mat_10)


###### CBS 50 ##########
load('ex1_500_n50.Rdata')

prob_est = lapply(model_list, function(y)
  lapply(y, function(x) exp(all_data %*% x$theta_hat) / sum(exp(all_data %*% x$theta_hat))))

kl_list = lapply(prob_est, function(y)
  lapply(y, function(x) -sum(prob_star * log(x))))

kl_mat = lapply(kl_list, unlist) %>% do.call('rbind', .)
boxplot(kl_mat)

tv_list = lapply(prob_est, function(y)
  lapply(y, function(x) sum(abs(x - prob_star))))

tv_mat = lapply(tv_list, unlist) %>% do.call('rbind', .)

hd_list = lapply(prob_est, function(y)
  lapply(y, function(x) sum((sqrt(prob_star) - sqrt(x))^2)))

hd_mat = lapply(hd_list, unlist) %>% do.call('rbind', .)

kl_return_mat_50 = cbind(kl_mat[, 1], kl_mat[, 5])
tv_return_mat_50 = cbind(tv_mat[, 1], tv_mat[, 5])
hd_return_mat_50 = cbind(hd_mat[, 1], hd_mat[, 5])

colMeans(kl_return_mat_50)
colMeans(tv_return_mat_50)
colMeans(hd_return_mat_50)


###### CBS 100 ##########
load('ex1_500_n100.Rdata')

prob_est = lapply(model_list, function(y)
  lapply(y, function(x) exp(all_data %*% x$theta_hat) / sum(exp(all_data %*% x$theta_hat))))


kl_list = lapply(prob_est, function(y)
  lapply(y, function(x) -sum(prob_star * log(x))))

kl_mat = lapply(kl_list, unlist) %>% do.call('rbind', .)
boxplot(kl_mat)

tv_list = lapply(prob_est, function(y)
  lapply(y, function(x) sum(abs(x - prob_star))))

tv_mat = lapply(tv_list, unlist) %>% do.call('rbind', .)

hd_list = lapply(prob_est, function(y)
  lapply(y, function(x) sum((sqrt(prob_star) - sqrt(x))^2)))

hd_mat = lapply(hd_list, unlist) %>% do.call('rbind', .)

kl_return_mat_100 = cbind(kl_mat[, 1], kl_mat[, 5])
tv_return_mat_100 = cbind(tv_mat[, 1], tv_mat[, 5])
hd_return_mat_100 = cbind(hd_mat[, 1], hd_mat[, 5])

colMeans(kl_return_mat_100)
colMeans(tv_return_mat_100)
colMeans(hd_return_mat_100)




#### summary ####
kl_return_mat = cbind(kl_return_mat_10, kl_return_mat_50, kl_return_mat_100, 
                      kl_btm_10_vec, kl_btm_50_vec, kl_btm_100_vec, kl_bt_vec)
tv_return_mat = cbind(tv_return_mat_10, tv_return_mat_50, tv_return_mat_100, 
                      tv_btm_10_vec, tv_btm_50_vec, tv_btm_100_vec, tv_bt_vec)
hd_return_mat = cbind(hd_return_mat_10, hd_return_mat_50, hd_return_mat_100, 
                      hd_btm_10_vec, hd_btm_50_vec, hd_btm_100_vec, hd_bt_vec)

# kl_return_mat = cbind(kl_return_mat_10, kl_return_mat_50, kl_return_mat_100, kl_return_mat_300, kl_bt_vec)
# tv_return_mat = cbind(tv_return_mat_10, tv_return_mat_50, tv_return_mat_100, tv_return_mat_300, tv_bt_vec)
# hd_return_mat = cbind(hd_return_mat_10, hd_return_mat_50, hd_return_mat_100, hd_return_mat_300, hd_bt_vec)

colnames(kl_return_mat) = colnames(tv_return_mat) = colnames(hd_return_mat) =
  c(paste(rep(c('BS', 'CBS'), 3), rep(c(10, 50, 100), each = 2)), paste('BTM', c(10, 50, 100)), 'BT')

# colnames(kl_return_mat) = colnames(tv_return_mat) = colnames(hd_return_mat) =
#   c(paste(rep(c('BS', 'CBS'), 4), rep(c(10, 50, 100, 300), each = 2)), 'BT')

head(kl_return_mat)

boxplot(kl_return_mat)
boxplot(tv_return_mat)
boxplot(hd_return_mat)

colMeans(kl_return_mat)
colMeans(tv_return_mat)
colMeans(hd_return_mat)

mean(kl_bt_vec)
mean(tv_bt_vec)


sim1_result = list(kl_return_mat, tv_return_mat, hd_return_mat)
save(sim1_result, file = '../sim1_result.Rdata')














