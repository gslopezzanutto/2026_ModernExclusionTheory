f_df = function(n){
  data = as.data.frame(matrix(NA, nrow = n, ncol = 21)) # Data frame containing the results from one iteration
  names(data) = c("n", "condition", "n_n_1", "N_star", "MVN_p-value", "ND_p-value", "FD_p-value", "c_ij", "var_c_ij", "r_i", "var_r_i", "eta_i", "var_eta_i", "ND", "mean_ND", "var_ND", "FD", "mean_FD", "var_FD", "cov_NFD", "cor_NFD")
  return(data)
}
transA = function(A = A, mu_ip = mu_ip){
  Ap = A
  for(i in 1:length(A[1,])){
    for(j in 1:length(A[1,])){
      Ap[i,j] = A[i,j]*mu_ip[j]/mu_ip[i]/A[j,j]
    }
  }
  Ap = as.matrix(Ap)
  return(Ap)
}
LV_model = function(Time, N_t, Pars){
  # Lotka-Volterra model
  with(as.list(c(N_t, Pars)),{
    dN_dt = (N_t*(mu_i-A%*%N_t))*(1/(1+exp(-10**7*(N_t-10**(-6)))))
    return(list(c(dN_dt)))
  })
}
get_LV_sim = function(State = runif(n, min = 0, max = 1), mu_i, A, Time = seq(0, 1000000, 1), thold) {
  ode(y = State, parms = list(A = A, mu_i = mu_i),
      times = Time, func = LV_model, method = "bdf", atol = 1e-8, rtol = 1e-6) |>
    as.data.frame() |>
    as_tibble() |>
    filter(time == max(time)) |>
    select(!time) |>
    pivot_longer(cols = everything(), names_to = "species", values_to = "density") |>
    mutate(species = as.integer(species))
}
get_survivors_indices = function(community_table) {
  community_table |>
    filter(density > thold) |>
    pull(species)
}
get_nonsurvivors_indices = function(community_table) {
  community_table |>
    filter(density < thold) |>
    pull(species)
}
get_solve_densities = function(n, mu_i, A, survivors_indices) {
  solve_densities = rep(0, n)
  A_surv = A[survivors_indices, survivors_indices]
  mu_i_surv = mu_i[survivors_indices]
  solve_densities[survivors_indices] = solve(A_surv, mu_i_surv)
  return(solve_densities)
}
get_condit = function(n, solve_densities, survivors_indices, nonsurvivors_indices){
  condition = rep(0, n)
  condition[survivors_indices] = "resident"
  condition[nonsurvivors_indices] = "excluded"
  return(condition)
}
get_inv_test = function(solve_densities, mu_i, A){
  mu_i - A%*%solve_densities
}
modify_condit_inv = function(inv_test, condition, nonsurvivors_indices){
  if(inv == TRUE){
    condition[nonsurvivors_indices] = ifelse(inv_test[nonsurvivors_indices] > 0, "invader", "excluded")
  } else {
    condition[nonsurvivors_indices] = "excluded"
  }
  return(condition)
}
get_N_i_star_surv = function(survivors_indices, mu_i, A){
  survs = length(survivors_indices)
  A_surv = A[survivors_indices, survivors_indices]
  mu_i_surv = mu_i[survivors_indices]
  
  N_i_star_surv = matrix(0, survs, survs) # Calculate n-1 equilibrium abundances
  for (i in 1:survs){
    ind = (1:survs)[1:survs != i]
    N_i_star_surv[i,ind] = solve(A_surv[ind,ind], mu_i_surv[ind])
    if (any(N_i_star_surv[i,] < 0)) N_i_star_surv[i,] = NaN
  }
  return(N_i_star_surv)
}
get_sub_equi = function(solve_densities, mu_i, A, survivors_indices, nonsurvivors_indices){
  sub_equi = matrix(0, length(mu_i), length(mu_i))
  N_i_star_surv = get_N_i_star_surv(survivors_indices, mu_i, A)
  
  sub_equi[survivors_indices, survivors_indices] = N_i_star_surv # columns with surviving species are filled with n-1 communities
  sub_equi[nonsurvivors_indices,] = rep(solve_densities, each = length(nonsurvivors_indices)) # columns with excluded species are filled with the equilibrium abundances
  
  return(sub_equi)
}
modify_condit_sub = function(sub_equi, condition){
  for (i in 1:length(sub_equi[1,])){
    if (any(is.nan(sub_equi[i,]))) condition[i] = "no_n_1"
  }
  return(condition)
}
get_r_i = function(mu_i, A, sub_equi){
  r_i = mu_i - einsum::einsum("ij, ij -> i", A, sub_equi)
}
get_c_ij = function(A, n){
  c_ij = matrix(0, n, n) # create c_ij matrix
  for (i in 1:n){
    for (j in 1:n){
      c_ij[i,j] = sqrt(abs(A[j,j]*A[i,j]/A[i,i]/A[j,i]))
    }
  }
  return(c_ij)
}
get_eta_i = function(c_ij, mu_i, A, sub_equi){
  eta_i = mu_i - diag(A)*einsum::einsum("ij, ij -> i", c_ij, sub_equi)
}
get_ND = function(mu_i, r_i, eta_i){
  ND = 1 - ((mu_i-r_i)/(mu_i-eta_i))
}
get_FD = function(mu_i, eta_i){
  FD = 1 - mu_i/(mu_i-eta_i)
}
CommNFD = function(A = A, mu_i = mu_i, thold = thold){
  n = length(A[1,])
  data = f_df(n)
  
  ode_densities = get_LV_sim(State = runif(n, min = 0, max = 1), mu_i, A, Time = seq(0, 1000000, 1), thold)
  equi = ode_densities$density
  survivors_indices = get_survivors_indices(ode_densities)
  nonsurvivors_indices = get_nonsurvivors_indices(ode_densities)
  
  survs = length(survivors_indices)
  
  solve_densities = get_solve_densities(n, mu_i, A, survivors_indices)
  condition = get_condit(n, solve_densities, survivors_indices, nonsurvivors_indices)
  if (identical(which(equi > thold), which(solve_densities > 0)) == FALSE || any(solve_densities < 0) == TRUE){
    print("Mismatch between ode and solve densities. Try increasing time sequence to ode.")
  }
  inv_test = get_inv_test(solve_densities, mu_i, A)
  condition = modify_condit_inv(inv_test, condition, nonsurvivors_indices)
  
  sub_equi = get_sub_equi(solve_densities, mu_i, A, survivors_indices, nonsurvivors_indices)
  condition = modify_condit_sub(sub_equi, condition)
  
  r_i = get_r_i(mu_i, A, sub_equi)
  c_ij = get_c_ij(A, n)
  eta_i = get_eta_i(c_ij, mu_i, A, sub_equi)
  ND = get_ND(mu_i, r_i, eta_i)
  FD = get_FD(mu_i, eta_i)
  
  df_comm_mvfit = data.frame(ND, FD)
  if(length(which(df_comm_mvfit$ND > 0)) > 7){
    mvn_result_hz = mvn(data = df_comm_mvfit, mvn_test = "hz")
    data[,"MVN_p-value"] = as.numeric(mvn_result_hz$multivariate_normality[3])
    data[,"ND_p-value"] = as.numeric(mvn_result_hz$univariate_normality[1,4])
    data[,"FD_p-value"] = as.numeric(mvn_result_hz$univariate_normality[2,4])
  } else {
    print("Not enough species to perform multivariate analysis. At least 8 species with calculated niche and fitness differences necessary")
  }
  
  
  data[,"n"] = n
  data[,"N_star"] = solve_densities
  data[,"condition"] = condition
  data[,"n_n_1"] = sum(data$condition == "no_n_1")/n
  
  data[,"c_ij"] = mean(c_ij[row(c_ij) != col(c_ij)])
  data[,"var_c_ij"] = var(c(c_ij), na.rm = TRUE)
  data[,"r_i"] = r_i
  data[,"var_r_i"] = var(r_i, na.rm = TRUE)
  data[,"eta_i"] = eta_i
  data[,"var_eta_i"] = var(eta_i, na.rm = TRUE)
  
  data[,"ND"] = ND
  data[,"mean_ND"] = mean(ND, na.rm = TRUE)
  data[,"var_ND"] = var(ND, na.rm = TRUE)
  data[,"FD"] = FD
  data[,"mean_FD"] = mean(FD, na.rm = TRUE)
  data[,"var_FD"] = var(FD, na.rm = TRUE)
  data[,"cov_NFD"] = cov(ND, FD, use = "complete.obs")
  data[,"cor_NFD"] = cor(ND, FD, use = "complete.obs")
  
  
  return(data)
}

AIC_test = function(data, mean_vec, cov_mat, k){
  like = c()
  for(i in 1:length(data[,1])){
    cte = 1/sqrt((2*pi)**2*det(cov_mat))
    x_min_mu = as.numeric(c(data[i,] - mean_vec))
    p_xi = cte*exp(-0.5*t(x_min_mu)%*%solve(cov_mat)%*%x_min_mu)
    like = append(like, p_xi)
  }
  p_data = sum(log(like))
  res = 2*k-2*p_data
}

comp_lambda = function(A = A){
  diag(A) = NA
  lambda = mean(A, na.rm = TRUE)
}
make_ellipse = function(mean_vec, cov_mat, level = 0.95, n_points = 500) {
  theta = seq(0, 2 * pi, length.out = n_points)
  circle = cbind(cos(theta), sin(theta))
  chisq_val = qchisq(level, df = 2)
  chol_mat = chol(cov_mat)
  ellipse_pts = circle %*% chol_mat * sqrt(chisq_val)
  ellipse_pts = sweep(ellipse_pts, 2, mean_vec, "+")
  data.frame(x = ellipse_pts[, 1], y = ellipse_pts[, 2])
}