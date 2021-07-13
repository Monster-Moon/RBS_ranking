
rm(list = ls())
gc(reset = T)

setwd('/Users/moon/Documents/CBS_model/code/real_data/soc_code')

list.files()

load('r1_result.Rdata')
load('r2_result.Rdata')
load('r3_result.Rdata')
load('r4_result.Rdata')
load('r5_result.Rdata')


mean_fun = function(x)
{
  if(is.matrix(x))
  {
    return(colMeans(x))
  }else
  {
    return(mean(x))
  }
}

sd_fun = function(x)
{
  if(is.matrix(x))
  {
    return(apply(x, 2, sd))
  }else
  {
    return(sd(x))
  }
}

### r1
r1_mean_list = lapply(r1_result, function(y) lapply(y, mean_fun))
r1_sd_list = lapply(r1_result, function(y) lapply(y, sd_fun))

r1_mean_list$test

lapply(r1_mean_list$valid[1:3], which.min) ## 4 :: KL 기준
r1_mean_list$test[1:3]
tmp = lapply(r1_mean_list$test[1:3], function(x) x[c(1, 2)])
tmp_sd = lapply(r1_sd_list$test[1:3], function(x) x[c(1, 2)])

r1_vec = c('Mtdot', 
           paste0(round(c(tmp$test_kl_mat, r1_mean_list$test$test_btm_kl_vec, r1_mean_list$test$test_bt_kl_vec), 2), '(',
                  round(c(tmp_sd$test_kl_mat, r1_sd_list$test$test_btm_kl_vec, r1_sd_list$test$test_bt_kl_vec), 2), ')'),
           paste0(round(c(tmp$test_hd_mat, r1_mean_list$test$test_btm_hd_vec, r1_mean_list$test$test_bt_hd_vec), 2), '(',
                  round(c(tmp_sd$test_hd_mat, r1_mean_list$test$test_btm_hd_vec, r1_sd_list$test$test_bt_hd_vec), 2), ')'))

r1_vec

#### r2
r2_mean_list = lapply(r2_result, function(y) lapply(y, mean_fun))
r2_sd_list = lapply(r2_result, function(y) lapply(y, sd_fun))

lapply(r2_mean_list$valid[1:3], which.min) ## 5 :: KL 기준

tmp = lapply(r2_mean_list$test[1:3], function(x) x[c(1, 5)])
tmp_sd = lapply(r2_sd_list$test[1:3], function(x) x[c(1, 5)])

r2_vec = c('Mtpuz', 
           paste0(round(c(tmp$test_kl_mat, r2_mean_list$test$test_btm_kl_vec, r2_mean_list$test$test_bt_kl_vec), 2), '(',
                  round(c(tmp_sd$test_kl_mat, r2_sd_list$test$test_btm_kl_vec, r2_sd_list$test$test_bt_kl_vec), 2), ')'),
           paste0(round(c(tmp$test_hd_mat, r2_mean_list$test$test_btm_hd_vec, r2_mean_list$test$test_bt_hd_vec), 2), '(',
                  round(c(tmp_sd$test_hd_mat, r2_mean_list$test$test_btm_hd_vec, r2_sd_list$test$test_bt_hd_vec), 2), ')'))
#### r3
r3_mean_list = lapply(r3_result, function(y) lapply(y, mean_fun))
r3_sd_list = lapply(r3_result, function(y) lapply(y, sd_fun))

lapply(r3_mean_list$valid[1:3], which.min) ## 4 :: KL 기준

tmp = lapply(r3_mean_list$test[1:3], function(x) x[c(1, 3)])
tmp_sd = lapply(r3_sd_list$test[1:3], function(x) x[c(1, 3)])

r3_vec = c('APA1', 
           paste0(round(c(tmp$test_kl_mat, r3_mean_list$test$test_btm_kl_vec, r3_mean_list$test$test_bt_kl_vec), 2), '(',
                  round(c(tmp_sd$test_kl_mat, r3_sd_list$test$test_btm_kl_vec, r3_sd_list$test$test_bt_kl_vec), 2), ')'),
           paste0(round(c(tmp$test_hd_mat, r3_mean_list$test$test_btm_hd_vec, r3_mean_list$test$test_bt_hd_vec), 2), '(',
                  round(c(tmp_sd$test_hd_mat, r3_mean_list$test$test_btm_hd_vec, r3_sd_list$test$test_bt_hd_vec), 2), ')'))
r3_vec

#### r4
r4_mean_list = lapply(r4_result, function(y) lapply(y, mean_fun))
r4_sd_list = lapply(r4_result, function(y) lapply(y, sd_fun))

lapply(r4_mean_list$valid[1:3], which.min) ## 4 :: KL 기준

tmp = lapply(r4_mean_list$test[1:3], function(x) x[c(1, 5)])
tmp_sd = lapply(r4_sd_list$test[1:3], function(x) x[c(1, 5)])

r4_vec = c('APA2', 
           paste0(round(c(tmp$test_kl_mat, r4_mean_list$test$test_btm_kl_vec, r4_mean_list$test$test_bt_kl_vec), 2), '(',
                  round(c(tmp_sd$test_kl_mat, r4_sd_list$test$test_btm_kl_vec, r4_sd_list$test$test_bt_kl_vec), 2), ')'),
           paste0(round(c(tmp$test_hd_mat, r4_mean_list$test$test_btm_hd_vec, r4_mean_list$test$test_bt_hd_vec), 2), '(',
                  round(c(tmp_sd$test_hd_mat, r4_mean_list$test$test_btm_hd_vec, r4_sd_list$test$test_bt_hd_vec), 2), ')'))
##### r5
r5_mean_list = lapply(r5_result, function(y) lapply(y, mean_fun))
r5_sd_list = lapply(r5_result, function(y) lapply(y, sd_fun))

lapply(r5_mean_list$valid[1:3], which.min) ## 5 :: KL 기준

tmp = lapply(r5_mean_list$test[1:3], function(x) x[c(1, 3)])
tmp_sd = lapply(r5_sd_list$test[1:3], function(x) x[c(1, 3)])

r5_vec = c('APA3', 
           paste0(round(c(tmp$test_kl_mat, r5_mean_list$test$test_btm_kl_vec, r5_mean_list$test$test_bt_kl_vec), 2), '(',
                  round(c(tmp_sd$test_kl_mat, r5_sd_list$test$test_btm_kl_vec, r5_sd_list$test$test_bt_kl_vec), 2), ')'),
           paste0(round(c(tmp$test_hd_mat, r5_mean_list$test$test_btm_hd_vec, r5_mean_list$test$test_bt_hd_vec), 2), '(',
                  round(c(tmp_sd$test_hd_mat, r5_mean_list$test$test_btm_hd_vec, r5_sd_list$test$test_bt_hd_vec), 2), ')'))
r5_vec

result_mat = rbind(r1_vec, r2_vec, r3_vec, r4_vec, r5_vec)
# result_mat = rbind(r1_vec, r2_vec, r5_vec)

colnames(result_mat) = c('Dataset', rep(c('BS', 'CBS', 'BTM', 'BT'), 2))

result_mat = cbind(rep(result_mat[, 1], 2),
      rep(c('KL', 'HD'), each = 5),
      rbind(result_mat[,2:5], 
            result_mat[,c(6:9)]))

xtable::xtable(result_mat) %>% print(., include.rownames = F)




result_mat

r1_result_mat$test[4:6]
colMeans(r1_result$train$train_hd_mat)
mean(r1_result$train$train_bt_hd_vec)




r2_result_mat = lapply(r2_result, function(y)
  lapply(y, mean_fun))

lapply(r2_result_mat$valid[1:3], which.min) ## 4
lapply(r2_result_mat$test[1:3], which.min) 
lapply(r2_result_mat$test[1:3], function(x) x[c(1, 4)])


