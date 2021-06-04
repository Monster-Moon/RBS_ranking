######################################################
###################### CBS model #####################
######################################################

if(!require(PerMallows)) install.packages('PerMallows'); require(PerMallows);
if(!require(dplyr)) install.packages('dplyr'); require(dplyr);
if(!require(BradleyTerry2)) install.packages('BradleyTerry2'); require(BradleyTerry2)

########### CBS model fitting code ############## 
#### rank <-> inx form
x_vec_generates = function(i, data) ## change one ranking sample into index form
{
  n_items = ncol(data)
  order_data = order(data[i, ])
  x_vec = lapply(1:(n_items - 1), 
                 function(x) ifelse(order_data[-c(1:x)] - order_data[x] > 0, 1, 0)) %>% 
    unlist()
  return(x_vec)
}

x_data_generates = function(data) ## change ranking samples into (1->2), (1->3), (1->4), .... , (4 -> 5) data form
{
  return(lapply(1:nrow(data), x_vec_generates, data = data) %>%
           do.call('rbind', .) %>% as.matrix())
}

permutations = function(n) ## makes permutation matrix
{
  if(n==1){
    return(matrix(1))
  } else {
    sp = permutations(n-1)
    p = nrow(sp)
    A = matrix(nrow=n*p,ncol=n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] = cbind(i,sp+(sp>=i))
    }
    return(A)
  }
}

A_mat_fun = function(n_items) ## A :: A %*% beta = theta 
{
  n_items_comb = n_items * (n_items - 1) / 2
  A_mat_left_upper_list = lapply(1:(n_items - 1), function(j) 
    lapply(1:(n_items - j), function(i) c(rep(0, j - 1), 
                                          1, 
                                          rep(0, i-1), 
                                          -1, 
                                          rep(0, n_items - i - j))))
  
  A_mat_left_upper = lapply(A_mat_left_upper_list, function(x) do.call('rbind', x)) %>% do.call('rbind', .)
  A_mat_right_upper = diag(1, n_items_comb)
  A_mat_upper = cbind(A_mat_left_upper, A_mat_right_upper)
  
  A_mat_left_lower = rbind(matrix(rep(1, n_items), nrow = 1),
                           matrix(0, nrow = n_items - 1, ncol = n_items))

  # A_mat_left_lower = rbind(matrix(c(1, rep(0, n_items - 1)), nrow = 1),
  #                          matrix(0, nrow = n_items - 1, ncol = n_items))

  
  # A_mat_right_lower = rbind(matrix(rep(0, n_items_comb), nrow = 1),
  #                           do.call('cbind',
  #                                   lapply((n_items - 1):1,
  #                                          function(i) rbind(matrix(0, nrow = n_items - i - 1, ncol = i), diag(c(1, rep(0, i-1)), i)))))
  
  A_mat_right_lower = rbind(matrix(rep(0, n_items_comb), nrow = 1),
                            do.call('cbind',
                                    lapply((n_items - 1):1,
                                           function(i) rbind(matrix(0, nrow = n_items - i - 1, ncol = i), diag(1, i)))))
  A_mat_lower = cbind(A_mat_left_lower, A_mat_right_lower)                    
  A_mat = rbind(A_mat_upper, A_mat_lower)
  return(A_mat)
}

transitivity_test_one_fun = function(x, test_data) ## transivity test for one data
{
  tmp_x = matrix(x, nrow = nrow(test_data), ncol = ncol(test_data), byrow = T) == test_data
  tmp_x_inx = apply(tmp_x, 1, all)
  return(list(test_result = any(tmp_x_inx), test_inx = which(tmp_x_inx)))
}

transitivity_test_fun = function(data, test_data, test_rank_data) ## transivity test for whole data
{
  test_result = apply(data, 1, transitivity_test_one_fun, test_data = test_data)
  tmp_inx = lapply(test_result, function(x) x$test_result) %>% unlist()
  tmp_rank_inx = lapply(test_result, function(x) x$test_inx) %>% unlist()
  return(list(data = data[tmp_inx, ], rank_data = test_rank_data[tmp_rank_inx, ]))
}

btm_fit = function(x, sim_n, eta, max_iter, lambda1, lambda2, eps = 1e-5) ## subgradient
{
  x_data = x_data_generates(x) ## rank form -> inx form
  lmm_fit = lmm(data = x) ## sigma0, phi estimates
  sim_sigma0 = lmm_fit$mode
  
  n_items = ncol(x)
  n_items_comb = n_items * (n_items - 1) / 2
  
  sim_sigma0_inx = combn(sim_sigma0, 2) %>% apply(., 2, function(x) paste(sort(x), collapse = ''))
  ordered_inx = combn(1:n_items, 2) %>% apply(., 2, function(x) paste(sort(x), collapse = ''))
  
  change_inx = match(ordered_inx, sim_sigma0_inx)
  sign_inx = ifelse(x_data_generates(data = matrix(sim_sigma0, nrow = 1)) == matrix(1, ncol = n_items_comb, nrow = 1),
                    1, -1) %>% as.numeric()
  change_mat = diag(length(change_inx))[change_inx, ] * sign_inx
  
  # x_ordered = apply(x, 1, function(i) sim_sigma0[i]) %>% t()
  x_ordered = apply(x, 1, function(y) match(y, sim_sigma0)) %>% t()
  x_data_ordered = x_data_generates(x_ordered)
  
  sim_phi = exp(-lmm_fit$theta)
  A_mat = A_mat_fun(n_items = n_items) ## design matrix
  A_solve_mat = solve(A_mat)
  
  #### Optimization ####
  theta_vec_0 = rep(lmm_fit$theta, n_items_comb)
  # Z_theta0 = lapply(1:(n_items-1), function(x) sum(sim_phi^c(0:x))) %>% do.call('prod', .)
  # Z_theta0 = Z_theta0 / (sim_phi^n_items_comb)
  
  # H = matrix(0, nrow = n_items_comb - 1, ncol = n_items_comb)
  # H[row(H) == col(H)] = 1
  # H[row(H) + 1 == col(H)] = -1
  # HtH = t(H) %*% H
  
  # set.seed(1)
  # tmp_bt_data = expand.grid(1:n_items, 1:n_items)
  # tmp_bt_data = tmp_bt_data[tmp_bt_data[, 1] < tmp_bt_data[, 2], ] %>% arrange(Var1)
  # 
  # bt_data = cbind(tmp_bt_data, t(apply(x_data_generates(x_ordered), 2, function(x) table(factor(x, levels = 1:0)))))
  # bt_data$Var1 = factor(bt_data$Var1, levels = 1:n_items)
  # bt_data$Var2 = factor(bt_data$Var2, levels = 1:n_items)
  # colnames(bt_data) = c('home_team', 'away_team', 'home_win', 'away_win')
  # bt_model = BTm(cbind(home_win, away_win), home_team, away_team, data = bt_data, id = "team")
  # theta_vec = A_mat[1:n_items_comb, 1:n_items] %*% (c(0, bt_model$coefficients) - mean(c(0, bt_model$coefficients)))
  # 

  theta_vec = rnorm(n_items_comb)

  log_likelihood_vec = numeric(max_iter)
  theta_mat = matrix(0, ncol = length(theta_vec), nrow = max_iter)
  
  D_mat = cbind(matrix(0, n_items_comb, n_items), diag(1, n_items_comb, n_items_comb)) 
  D_theta_mat = rbind(diag(1, n_items_comb, n_items_comb), matrix(0, n_items, n_items_comb))
  
  S_mat = D_mat %*% A_solve_mat %*% D_theta_mat %*% change_mat
  theta_update_fixed_term = -colMeans(x_data_ordered)
  
  grad_norm_vec = numeric(max_iter)
  log_likelihood_vec = numeric(max_iter)
  for(iter in 1:max_iter)
  {
    sim_rank_data = rim(n = sim_n, sigma0 = 1:n_items, phi = sim_phi) ## random sample generates
    x_sim_data = x_data_generates(sim_rank_data)
    
    theta_mat[iter, ] = theta_vec
    h_l = exp(x_sim_data %*% (theta_vec - theta_vec_0))
    xh_l = x_sim_data * matrix(h_l, nrow = nrow(x_sim_data), ncol = ncol(x_sim_data))
    
    # Z_theta_approx = colMeans(h_l) * Z_theta0
    # log_likelihood_vec[iter] = colSums(x_data %*% theta_vec) - log(Z_theta_approx))
    agrad = theta_update_fixed_term + 
      colSums(xh_l) / colSums(h_l) + 
      lambda1 * t(S_mat) %*% sign(S_mat %*% theta_vec) 
    # lambda2 * HtH %*% theta_vec
    
    theta_vec = theta_vec - eta * agrad
    beta_vec = A_solve_mat %*% c(theta_vec, rep(0, n_items))
    beta_vec[(n_items+1):length(beta_vec)] = 0
    theta_vec = (A_mat %*% beta_vec)[1:n_items_comb]
    
    grad_norm_vec[iter] = norm(agrad, '2')
    if(grad_norm_vec[iter] <= eps) break;
    # cat(iter, '\n')
  }
  
  theta_hat = colMeans(theta_mat[round(0.7 * max_iter):max_iter, ])
  theta_hat_return = change_mat %*% theta_hat
  beta_hat = solve(A_mat) %*% c(theta_hat_return, rep(0, n_items))
  
  df_mat = A_mat[-c(1:(n_items_comb)), ]
  df_mat[, abs(beta_hat) <= 1e-3] = 0
  df = sum(pmax(0, rowSums(df_mat) - 1))
  
  return(list(mm = lmm_fit,
              iter = iter,
              grad_norm_vec = grad_norm_vec[1:iter],
              # log_likelihood_vec = log_likelihood_vec[1:iter],
              change_mat = change_mat,
              # Z_theta0 = Z_theta0,
              # Z_theta_approx = Z_theta_approx,
              theta_mat = theta_mat[1:iter, ],
              theta_hat = theta_hat_return,
              beta_hat = beta_hat,
              df = df))
}

cbs_fit = function(x, sim_n, eta, max_iter, lambda1, lambda2, eps = 1e-5) ## subgradient
{
  x_data = x_data_generates(x) ## rank form -> inx form
  lmm_fit = lmm(data = x) ## sigma0, phi estimates
  sim_sigma0 = lmm_fit$mode
  
  n_items = ncol(x)
  n_items_comb = n_items * (n_items - 1) / 2
  
  sim_sigma0_inx = combn(sim_sigma0, 2) %>% apply(., 2, function(x) paste(sort(x), collapse = ''))
  ordered_inx = combn(1:n_items, 2) %>% apply(., 2, function(x) paste(sort(x), collapse = ''))
  
  change_inx = match(ordered_inx, sim_sigma0_inx)
  sign_inx = ifelse(x_data_generates(data = matrix(sim_sigma0, nrow = 1)) == matrix(1, ncol = n_items_comb, nrow = 1),
                    1, -1) %>% as.numeric()
  change_mat = diag(length(change_inx))[change_inx, ] * sign_inx
  
  # x_ordered = apply(x, 1, function(i) sim_sigma0[i]) %>% t()
  x_ordered = apply(x, 1, function(y) match(y, sim_sigma0)) %>% t()
  x_data_ordered = x_data_generates(x_ordered)
  
  sim_phi = exp(-lmm_fit$theta)
  A_mat = A_mat_fun(n_items = n_items) ## design matrix
  
  #### Optimization ####
  theta_vec_0 = rep(lmm_fit$theta, n_items_comb)
  # Z_theta0 = lapply(1:(n_items-1), function(x) sum(sim_phi^c(0:x))) %>% do.call('prod', .)
  # Z_theta0 = Z_theta0 / (sim_phi^n_items_comb)

  # H = matrix(0, nrow = n_items_comb - 1, ncol = n_items_comb)
  # H[row(H) == col(H)] = 1
  # H[row(H) + 1 == col(H)] = -1
  # HtH = t(H) %*% H
  
  # set.seed(1)
  # tmp_bt_data = expand.grid(1:n_items, 1:n_items)
  # tmp_bt_data = tmp_bt_data[tmp_bt_data[, 1] < tmp_bt_data[, 2], ] %>% arrange(Var1)
  # 
  # bt_data = cbind(tmp_bt_data, t(apply(x_data_generates(x_ordered), 2, function(x) table(factor(x, levels = 1:0)))))
  # bt_data$Var1 = factor(bt_data$Var1, levels = 1:n_items)
  # bt_data$Var2 = factor(bt_data$Var2, levels = 1:n_items)
  # colnames(bt_data) = c('home_team', 'away_team', 'home_win', 'away_win')
  # bt_model = BTm(cbind(home_win, away_win), home_team, away_team, data = bt_data, id = "team")
  # theta_vec = A_mat[1:n_items_comb, 1:n_items] %*% (c(0, bt_model$coefficients) - mean(c(0, bt_model$coefficients)))
  # 
  
  theta_vec =  rnorm(n_items_comb)
  log_likelihood_vec = numeric(max_iter)
  theta_mat = matrix(0, ncol = length(theta_vec), nrow = max_iter)
  
  D_mat = cbind(matrix(0, n_items_comb, n_items), diag(1, n_items_comb, n_items_comb)) 
  D_theta_mat = rbind(diag(1, n_items_comb, n_items_comb), matrix(0, n_items, n_items_comb))
  
  S_mat = D_mat %*% solve(A_mat) %*% D_theta_mat %*% change_mat
  theta_update_fixed_term = -colMeans(x_data_ordered)
  
  grad_norm_vec = numeric(max_iter)
  log_likelihood_vec = numeric(max_iter)
  for(iter in 1:max_iter)
  {
    sim_rank_data = rim(n = sim_n, sigma0 = 1:n_items, phi = sim_phi) ## random sample generates
    x_sim_data = x_data_generates(sim_rank_data)
    
    theta_mat[iter, ] = theta_vec
    h_l = exp(x_sim_data %*% (theta_vec - theta_vec_0))
    xh_l = x_sim_data * matrix(h_l, nrow = nrow(x_sim_data), ncol = ncol(x_sim_data))

    # Z_theta_approx = colMeans(h_l) * Z_theta0
    # log_likelihood_vec[iter] = colSums(x_data %*% theta_vec) - log(Z_theta_approx))
    agrad = theta_update_fixed_term + 
      colSums(xh_l) / colSums(h_l) + 
      lambda1 * t(S_mat) %*% sign(S_mat %*% theta_vec) 
      # lambda2 * HtH %*% theta_vec
    
    theta_vec = theta_vec - eta * agrad

    grad_norm_vec[iter] = norm(agrad, '2')
    if(grad_norm_vec[iter] <= eps) break;
    # cat(iter, '\n')
  }
  
  theta_hat = colMeans(theta_mat[round(0.7 * max_iter):max_iter, ])
  theta_hat_return = change_mat %*% theta_hat
  beta_hat = solve(A_mat) %*% c(theta_hat_return, rep(0, n_items))
  
  df_mat = A_mat[-c(1:(n_items_comb)), ]
  df_mat[, abs(beta_hat) <= 1e-3] = 0
  df = sum(pmax(0, rowSums(df_mat) - 1))
  
  return(list(mm = lmm_fit,
              iter = iter,
              grad_norm_vec = grad_norm_vec[1:iter],
              # log_likelihood_vec = log_likelihood_vec[1:iter],
              change_mat = change_mat,
              # Z_theta0 = Z_theta0,
              # Z_theta_approx = Z_theta_approx,
              theta_mat = theta_mat[1:iter, ],
              theta_hat = theta_hat_return,
              beta_hat = beta_hat,
              df = df))
}

cbs_fit_CD = function(x, eta, max_iter, lambda1, lambda2, eps = 1e-5, transition_n = 5) ## subgradient
{
  x_data = x_data_generates(x) ## rank form -> inx form
  lmm_fit = lmm(data = x) ## sigma0, phi estimates
  
  sim_sigma0 = lmm_fit$mode
  sim_phi = exp(-lmm_fit$theta)
  
  n_items = ncol(x)
  n_items_comb = n_items * (n_items - 1) / 2
  
  theta_vec_0 = rep(lmm_fit$theta, n_items_comb)
  
  A_mat = A_mat_fun(n_items = n_items) ## design matrix
  
  sim_sigma0_inx = combn(sim_sigma0, 2) %>% apply(., 2, function(x) paste(sort(x), collapse = ''))
  ordered_inx = combn(1:n_items, 2) %>% apply(., 2, function(x) paste(sort(x), collapse = ''))
  
  change_inx = match(ordered_inx, sim_sigma0_inx)
  sign_inx = ifelse(x_data_generates(data = matrix(sim_sigma0, nrow = 1)) == matrix(1, ncol = n_items_comb, nrow = 1),
                    1, -1) %>% as.numeric()
  change_mat = diag(length(change_inx))[change_inx, ] * sign_inx
  
  # x_ordered = apply(x, 1, function(i) sim_sigma0[i]) %>% t()
  x_ordered = apply(x, 1, function(y) match(y, sim_sigma0)) %>% t()
  x_data_ordered = x_data_generates(x_ordered)
  theta_update_fixed_term = -colMeans(x_data_ordered)
  
  sim_phi = exp(-lmm_fit$theta)
  
  #### Optimization ####
  H = matrix(0, nrow = n_items_comb - 1, ncol = n_items_comb)
  H[row(H) == col(H)] = 1
  H[row(H) + 1 == col(H)] = -1
  HtH = t(H) %*% H
  
  # set.seed(1)
  theta_vec =  rnorm(n_items_comb)
  
  log_likelihood_vec = numeric(max_iter)
  theta_mat = matrix(0, ncol = length(theta_vec), nrow = max_iter)
  
  D_mat = cbind(matrix(0, n_items_comb, n_items), diag(1, n_items_comb, n_items_comb)) 
  D_theta_mat = rbind(diag(1, n_items_comb, n_items_comb), matrix(0, n_items, n_items_comb))
  
  S_mat = D_mat %*% solve(A_mat) %*% D_theta_mat %*% change_mat
  
  grad_norm_vec = numeric(max_iter)
  log_likelihood_vec = numeric(max_iter)
  x_nrow = nrow(x_data_ordered)
  for(iter in 1:max_iter)
  {
    theta_mat[iter, ] = theta_vec
    transition_data = x_data_ordered
    for(k in 1:transition_n)
    {
      sim_rank_data = rim(n = x_nrow, sigma0 = 1:n_items, phi = sim_phi)
      x_sim_data = x_data_generates(sim_rank_data)

      sam_prob_vec = pmin(exp((x_sim_data - transition_data) %*% (theta_vec - theta_vec_0)), 1)
      binom_vec = rbinom(n = x_nrow, size = 1, prob = as.numeric(sam_prob_vec))
      binom_mat = matrix(binom_vec, nrow = x_nrow, ncol = n_items_comb)
      transition_data = binom_mat * x_sim_data + (1 - binom_mat) * transition_data
    }

    # for(j in 1:nrow(transition_data))
    # {
    #   for(k in 1:transition_n)
    #   {
    #     sim_rank_data = rim(n = 1, sigma0 = 1:n_items, phi = sim_phi) ## random sample generates
    #     x_sim_data = x_data_generates(sim_rank_data)
    #     sam_prob = min(exp((x_sim_data - transition_data[j, ]) %*% (theta_vec - theta_vec_0)), 1)
    #     if(sample(0:1, size = 1, prob = c(1 - sam_prob, sam_prob)) == 1)
    #     {
    #       transition_data[j, ] = x_sim_data
    #     }
    #   }
    # }
    # apply(x_data_ordered == transition_data, 1, all) %>% summary()
    
    log_likelihood_vec[iter] = colSums(x_data_ordered %*% theta_vec)
    
    agrad = theta_update_fixed_term + colMeans(transition_data) +
      lambda1 * t(S_mat) %*% sign(S_mat %*% theta_vec) + 
      lambda2 * HtH %*% theta_vec
    
    theta_vec = theta_vec - eta * agrad
    grad_norm_vec[iter] = norm(agrad, '2')
    if(grad_norm_vec[iter] <= eps) break;
    # cat(iter, '\n')
  }
  
  theta_hat = theta_vec
  theta_hat_return = change_mat %*% theta_hat
  beta_hat = solve(A_mat) %*% c(theta_hat_return, rep(0, n_items))
  return(list(mm = lmm_fit,
              iter = iter,
              theta_mat = theta_mat[1:iter, ],
              chage_mat = change_mat,
              theta_hat = theta_hat_return,
              beta_hat = beta_hat))
}

#### generates data from mallows model with RIM ####
sample_position = function(x, phi) ## generates random position
{
  phi_vec = phi^c(x:1)
  return(sample(1:x, size = 1, prob = phi_vec / sum(phi_vec), replace = T))
}

insertion_position = function(i, insertion_vec, sigma0) ## generates insertion vector 
{
  tmp_vec = numeric(i)
  tmp_vec[insertion_vec[i]] = sigma0[i]
  return(tmp_vec)
}

rim_onesample = function(sigma0, phi) ## generates one rank sample
{
  insertion_vec = lapply(1:length(sigma0), sample_position, phi = phi) %>% unlist()
  insertion_position_list = lapply(length(insertion_vec):1, FUN = insertion_position,
                                   insertion_vec = insertion_vec, sigma0 = sigma0)
  
  rank_sample = insertion_position_list[[1]]
  for(i in 2:length(insertion_position_list))
  {
    rank_sample[rank_sample == 0] = insertion_position_list[[i]]
  }
  return(rank_sample)
}

rim = function(n, sigma0, phi) ## generates n rank sample
{
  return(t(replicate(n, rim_onesample(sigma0 = sigma0, phi = phi))))
}

constraints_check_fun = function(x, A_mat)
{
  non_zero_inx = apply(A_mat[-c(1:n_items_comb), ], 1, function(x) which(x != 0))
  return(lapply(non_zero_inx, function(x) beta_star[x]))
}

# cbs_fit = function(x, sim_n, eta, max_iter, lambda, eps = 1e-5)
# {
#   x_data = x_data_generates(x) ## rank form -> inx form
#   
#   lmm_fit = lmm(data = x) ## sigma0, phi estimates
#   sim_sigma0 = lmm_fit$mode
#   # sim_sigma0 = c(1, 2, 3, 4)
#   sim_phi = exp(-lmm_fit$theta)
#   
#   # sim_rank_data = rmm(sim_n, sigma0 = lmm_fit$mode, theta = lmm_fit$mode)
#   sim_rank_data = rim(n = sim_n, sigma0 = sim_sigma0, phi = sim_phi) ## random sample generates
#   x_sim_data = x_data_generates(sim_rank_data)
#   
#   n_items = ncol(x)
#   n_items_comb = n_items * (n_items - 1) / 2
#   
#   A_mat = A_mat_fun(n_items = n_items) ## design matrix
#   
#   #### Optimization ####
#   theta_vec_0 = rep(-log(sim_phi), n_items_comb)
#   Z_theta0 = lapply(1:(n_items-1), function(x) sum(sim_phi^c(0:x))) %>% do.call('prod', .)
#   
#   # H = matrix(0, nrow = n_items_comb - 1, ncol = n_items_comb)
#   # H[row(H) == col(H)] = 1
#   # H[row(H) + 1 == col(H)] = -1
#   # HtH = t(H) %*% H
#   
#   grad_norm_vec = numeric(max_iter)
#   log_likelihood_vec = numeric(max_iter)
#   
#   set.seed(1)
#   theta_vec = rnorm(n_items_comb)
#   theta_mat = matrix(0, ncol = length(theta_vec), nrow = max_iter)
#   
#   for(iter in 1:max_iter)
#   {
#     theta_mat[iter, ] = theta_vec
#     h_l = exp(x_sim_data %*% (theta_vec - theta_vec_0))
#     xh_l = x_sim_data * matrix(h_l, nrow = nrow(x_sim_data), ncol = ncol(x_sim_data))
#     
#     Z_theta_approx = colMeans(h_l) * Z_theta0
#     log_likelihood_vec[iter] = colSums(x_data %*% theta_vec - log(Z_theta_approx))
#     agrad = -colMeans(x_data) + colSums(xh_l) / colSums(h_l) # + lambda * HtH %*% theta_vec
#     theta_vec = theta_vec - eta * agrad
#     
#     grad_norm_vec[iter] = norm(agrad, '2')
#     if(grad_norm_vec[iter] <= eps) break;
#   }
#   
#   beta_vec = solve(A_mat) %*% c(theta_vec, rep(0, n_items))
#   return(list(mm = lmm_fit,
#               iter = iter, 
#               grad_norm_vec = grad_norm_vec[1:iter],
#               log_likelihood_vec = log_likelihood_vec[1:iter],
#               Z_theta0 = Z_theta0,
#               Z_theta_approx = Z_theta_approx,
#               theta_mat = theta_mat[1:iter, ],
#               theta_vec = theta_vec,
#               beta_vec = beta_vec))
# }
# cbs_fit_lasso = function(x, sim_n, eta, eta_z, admm_iter, max_iter, lambda1, lambda2, rho, eps = 1e-6) ## ADMM
# {
#   x_data = x_data_generates(x) ## rank form -> inx form
#   lmm_fit = lmm(data = x) ## sigma0, phi estimates
#   sim_sigma0 = lmm_fit$mode
#   
#   n_items = ncol(x)
#   n_items_comb = n_items * (n_items - 1) / 2
#   
#   sim_sigma0_inx = combn(sim_sigma0, 2) %>% apply(., 2, function(x) paste(sort(x), collapse = ''))
#   ordered_inx = combn(1:n_items, 2) %>% apply(., 2, function(x) paste(sort(x), collapse = ''))
#   
#   change_inx = match(ordered_inx, sim_sigma0_inx)
#   sign_inx = ifelse(x_data_generates(data = matrix(sim_sigma0, nrow = 1)) == matrix(1, ncol = n_items_comb, nrow = 1),
#                     1, -1) %>% as.numeric()
#   change_mat = diag(length(change_inx))[change_inx, ] * sign_inx
#   
#   x_ordered = apply(x, 1, function(i) sim_sigma0[i]) %>% t()
#   x_data_ordered = x_data_generates(x_ordered)
#   
#   sim_phi = exp(-lmm_fit$theta)
#   sim_rank_data = rim(n = sim_n, sigma0 = 1:n_items, phi = sim_phi) ## random sample generates
#   
#   x_sim_data = x_data_generates(sim_rank_data)
#   A_mat = A_mat_fun(n_items = n_items) ## design matrix
#   
#   #### Optimization ####
#   theta_vec_0 = rep(lmm_fit$theta, n_items_comb)
#   Z_theta0 = lapply(1:(n_items-1), function(x) sum(sim_phi^c(0:x))) %>% do.call('prod', .)
#   Z_theta0 = Z_theta0 / (sim_phi^n_items_comb)
#   
#   
#   # H = matrix(0, nrow = n_items_comb - 1, ncol = n_items_comb)
#   # H[row(H) == col(H)] = 1
#   # H[row(H) + 1 == col(H)] = -1
#   # HtH = t(H) %*% H
#   
#   # grad_norm_vec = numeric(max_iter)
#   log_likelihood_vec = numeric(max_iter)
#   
#   # set.seed(1)
#   theta_vec =  rnorm(n_items_comb)
#   z_vec = u_vec = rnorm(n_items + n_items_comb)
#   
#   theta_mat = matrix(0, ncol = length(theta_vec), nrow = max_iter)
#   z_mat = matrix(0, ncol = length(z_vec), nrow = max_iter)
#   u_mat = matrix(0, ncol = length(u_vec), nrow = max_iter)
#   
#   D_mat = cbind(matrix(0, n_items_comb, n_items), diag(1, n_items_comb, n_items_comb)) 
#   D_theta_mat = rbind(diag(1, n_items_comb, n_items_comb), matrix(0, n_items, n_items_comb))
#   
#   S_mat = solve(A_mat) %*% D_theta_mat %*% change_mat
#   theta_update_fixed_term = -colMeans(x_data_ordered) ## 확률 형태로 계산 변경
#   
#   lambda2_vec = c(rep(0, n_items), rep(lambda2, n_items_comb))
#   
#   for(iter in 1:admm_iter)
#   {
#     theta_mat[iter, ] = theta_vec
#     z_mat[iter, ] = z_vec
#     u_mat[iter, ] = u_vec
#     
#     ## theta update
#     for(inner_iter in 1:max_iter) ## 
#     {
#       h_l = exp(x_sim_data %*% (theta_vec - theta_vec_0))
#       xh_l = x_sim_data * drop(h_l)
#       # matrix(h_l, nrow = nrow(x_sim_data), ncol = ncol(x_sim_data))
#       
#       Z_theta_approx = colMeans(h_l) * Z_theta0
#       agrad =  theta_update_fixed_term + colSums(xh_l) / colSums(h_l) + 
#         # lambda1 * HtH %*% theta_vec +
#         - rho * t(S_mat) %*% (z_vec - S_mat %*% theta_vec + u_vec)
#       theta_vec = theta_vec - eta * agrad
#       if(norm(agrad, '2') <= 1e-8) break;
#     }
#     
#     ## z update
#     z_vec = S_mat %*% theta_vec - u_vec
#     z_vec = pmax(0, (z_vec - lambda2_vec / rho)) - 
#       pmax(0, (-z_vec - lambda2_vec / rho)) ## prox
#     
#     ## u update
#     u_vec = u_vec + z_vec - S_mat %*% theta_vec ## u update
#     
#     ## log likelihood
#     log_likelihood_vec[iter] = colSums(x_data_ordered %*% theta_vec - log(Z_theta_approx))
#     
#     primal_resi = z_vec - S_mat %*% theta_vec
#     dual_resi = t(S_mat) %*% (z_vec - z_mat[iter, ])
#     
#     primal_resi_norm = norm(primal_resi, '2')
#     dual_resi_norm = norm(dual_resi, '2')
#     # if(primal_resi_norm >= dual_resi_norm * 10) rho = 2 * rho
#     # if(dual_resi_norm >= primal_resi_norm * 10) rho = rho / 2
#     # 
#     # 
#     # cat(iter, '\n')
#     # cat('primal residual: ', primal_resi_norm, '\n')
#     # cat('dual residual: ', dual_resi_norm, '\n')
#     # cat('s_mat & theta_vec: ', S_mat %*% theta_vec, '\n')
#     # cat('z_vec :', z_vec, '\n')
#     if(primal_resi_norm < eps & dual_resi_norm < eps) break
#   }
#   
#   theta_vec_return = change_mat %*% theta_vec
#   
#   df_mat = A_mat[-c(1:(n_items_comb)), ]
#   df_mat[, z_vec == 0] = 0
#   df = sum(pmax(0, rowSums(df_mat) - 1))
#   
#   aic = 2 * df - 2 * log_likelihood_vec[iter] ## threshold 값 확인 필요
#   bic = nrow(x) * df - 2 * log_likelihood_vec[iter] ## beta vec의 제약조건이 threshold에 따라서 달라진다 (문제점)
#   
#   return(list(iter = iter, 
#               theta_mat = theta_mat[1:iter, ],
#               u_mat = u_mat[1:iter, ],
#               z_mat = z_mat[1:iter, ],
#               log_likelihood_vec = log_likelihood_vec[1:iter],
#               aic = aic,
#               bic = bic,
#               Z_theta0 = Z_theta0,
#               Z_theta_approx = Z_theta_approx,
#               theta_vec = theta_vec,
#               u_vec = u_vec,
#               z_vec = z_vec,
#               df = df
#   ))
# }
