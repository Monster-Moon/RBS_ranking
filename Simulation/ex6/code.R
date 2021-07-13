rm(list = ls())
gc(reset = T)

setwd('/Users/moon/Documents/CBS_model/code/simulation_code/ex6')
source('../model.R')

n_items = 5
n_items_comb = n_items * (n_items - 1) / 2

A_mat = A_mat_fun(n_items = n_items) ## design matrix
            # 12,  13,   14,   15,    23,   24,   25,    34,   35,  45
gamma_star = c(0, 0.1, 0.1, -0.1, -0.1, 0.1, 0.1, -0.2, -0.1, 0.1)
gamma_ratio = 5
beta_star = c(-0.01, -0.01, 0, 0.01, 0.01, gamma_star * gamma_ratio)

constraints_check_fun(beta_star, A_mat) %>% lapply(., sum)

theta_star = A_mat %*% beta_star
theta_star = theta_star[1:n_items_comb]

permu_items = permutations(n_items)
all_data = x_data_generates(data = permu_items)

prob_star = exp(all_data %*% theta_star) / sum(exp(all_data %*% theta_star))
prob_mat = matrix(0, nrow = n_items, ncol = n_items)

tmp_prob = lapply(1:ncol(all_data), function(i) c(sum(all_data[, i] * prob_star))) %>% unlist()
l = 1
for(i in 1:(ncol(prob_mat)-1))
{
  for(j in (i+1):ncol(prob_mat))
  {
    prob_mat[i, j] = tmp_prob[l]
    l = l + 1
  }
}

l = 1
for(j in 1:(ncol(prob_mat)-1))
{
  for(i in (j+1):nrow(prob_mat))
  {
    prob_mat[i, j] = 1 - tmp_prob[l]
    l = l + 1
  }
}

prob_mat
comb_mat = t(combn(1:n_items, 2))
z = numeric(nrow(comb_mat))
v = matrix(0, nrow = nrow(comb_mat), ncol = n_items)
for(i in 1:nrow(comb_mat))
{
  z[i] = sign(prob_mat[comb_mat[i, 2], comb_mat[i, 1]] - prob_mat[comb_mat[i, 1], comb_mat[i, 2]])
  v[i, ] = prob_mat[ , comb_mat[i, 1]] - prob_mat[, comb_mat[i, 2]]  ## v_12
}

z = factor(z, levels = c(-1, 1))
svm_model = e1071::svm(z ~ v, type = 'C-classification', kernel = 'linear')

all(z == predict(svm_model, decision.values = T))
which(z != predict(svm_model, decision.values = T))

prob_mat

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

save(x_data_list, file = 'ex6_data_list_500.Rdata')




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

