rm(list = ls())
gc(reset = T)

source('..//model.R')
n_items = 3
n_items_comb = n_items * (n_items - 1) / 2

gamma_star = c(0, -0.5, 0.5)
gamma_ratio = 1
beta_star = c(0.1, -0.15, 0.05, gamma_star * gamma_ratio)

A_mat = A_mat_fun(n_items)
theta_star = A_mat %*% beta_star
theta_star = theta_star[1:n_items_comb]

permu_items = permutations(n_items)
all_data = x_data_generates(data = permu_items)

prob_star = exp(all_data %*% theta_star) / sum(exp(all_data %*% theta_star))
prob_mat = matrix(0, nrow = n_items, ncol = n_items)

tt = lapply(1:ncol(all_data), function(i) c(sum(all_data[, i] * prob_star))) %>% unlist()
l = 1
for(i in 1:(ncol(prob_mat)-1))
{
  for(j in (i+1):ncol(prob_mat))
  {
    prob_mat[i, j] = tt[l]
    l = l + 1
  }
}
prob_mat
l = 1
for(j in 1:(ncol(prob_mat)-1))
{
  for(i in (j+1):nrow(prob_mat))
  {
    prob_mat[i, j] = 1 - tt[l]
    l = l + 1
  }
}

comb_mat = t(combn(1:3, 2))
count_vec = numeric(nrow(comb_mat))
for(i in 1:nrow(comb_mat))
{
  count_vec[i] = xor((prob_mat[comb_mat[i, 1], comb_mat[i, 2]] > prob_mat[comb_mat[i, 2], comb_mat[i, 1]]),
      (sum(prob_mat[, comb_mat[i, 2]]) > sum(prob_mat[, comb_mat[i, 1]])))
}

count_vec
comb_mat[count_vec == 1, ]
prob_mat[2, 3] >= prob_mat[3, 2]


n = 500L
sample_prob = exp(all_data %*% theta_star) / sum(exp(all_data %*% theta_star))

x_data_list = vector('list', 100)
for(i in 1:100)
{
  sample_inx = sample(1:nrow(all_data), size = n, replace = T, prob = sample_prob)
  sample_inx = factor(sample_inx, levels = 1:nrow(all_data))
  sample_inx_prop = table(sample_inx) / sum(table(sample_inx))
  
  x_data_list[[i]] = permu_items[sample_inx, ]
  cat(i, '\n')
}

save(x_data_list, file = 'ex4_data_list_500.Rdata')




count
i = 1
j = 2


xor(F, F)
xor(T, F)
xor(F, T)
xor(T, T)


prob_mat[3, 4] > prob_mat[4, 3]
sum(prob_mat[, 4]) > sum(prob_mat[, 3])
colSums(prob_mat)
prob_mat



upper.tri(prob_mat)[1]
lapply(1:ncol(all_data), function(i) c(sum(all_data[, i] * prob_star))) %>% unlist()

lapply(1:ncol(all_data), function(i) c(sum(all_data[, i] * prob_star), 1- sum(all_data[, i] * prob_star)))

head(all_data)

x = model_list[[1]][[1]]

prob_tmp = exp(all_data %*% x$theta_mat[421, ]) / sum(exp(all_data %*% x$theta_mat[421, ]))
-sum(prob_star * log(prob_tmp))


colSums(all_data * matrix(prob_tmp, nrow = nrow(all_data), ncol = ncol(all_data)))
sum(all_data[, 1] * prob_tmp)

colSums(all_data) / nrow(all_data)
prob_tmp

