rm(list = ls())
gc(reset = T)

setwd('D:\\Paper\\cbs_model\\code\\simulation_code\\ex4')
source('..\\model.R')

if(!require(BradleyTerry2)) install.packages('BradleyTerry2'); require(BradleyTerry2)
load('ex4_data_list_500.Rdata')

n_items = 3
n_items_comb = n_items * (n_items - 1) / 2

A_mat = A_mat_fun(n_items)
tmp_bt_data = expand.grid(1:n_items, 1:n_items)
tmp_bt_data = tmp_bt_data[tmp_bt_data[, 1] < tmp_bt_data[, 2], ] %>% arrange(Var1)

bt_data_list = vector('list', 100)
bt_model_list = vector('list', 100)

for(i in 1:100)
{
  bt_data = cbind(tmp_bt_data, t(apply(x_data_generates(x_data_list[[i]]), 2, table))[,2:1])
  bt_data$Var1 = factor(bt_data$Var1, levels = 1:n_items)
  bt_data$Var2 = factor(bt_data$Var2, levels = 1:n_items)
  colnames(bt_data) = c('home_team', 'away_team', 'home_win', 'away_win')
  bt_data_list[[i]] = bt_data
  bt_model_list[[i]] = BTm(cbind(home_win, away_win), home_team, away_team, data = bt_data, id = "team")
  cat(i, '\n')
}

save(bt_model_list, file = 'ex4_500_bt.Rdata')
