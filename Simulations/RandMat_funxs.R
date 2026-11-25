f_OneMatrix = function(n){
  data = as.data.frame(matrix(NA, nrow = n, ncol = 36)) # Data frame containing the results from one iteration
  names(data) = c("matrix", "n", "sigma", "lambda", "rho", "condition", "proportion", "N_star", "no_n_1", "c_ij", "var_c_ij", "r_i", "var_r_i", "eta_i", "var_eta_i", "ND", "mean_ND", "var_ND", "FD", "mean_FD", "var_FD", "cov_NFD", "cor_NFD", "pvalue", "theo_c_ij", "theo_var_c_ij", "theo_r_i", "theo_var_r_i", "theo_eta_i", "theo_var_eta_i", "theo_ND", "theo_FD", "theo_var_ND", "theo_var_FD", "theo_cov_NFD", "theo_cor_NFD")
  return(data)
}
create_matrix = function(n, sig, rho) {
  # Creates a matrix with random entries, entries have variance sig and correlation rho
  # Returns itera number of matrices with shape (n,n)
  # Entries follow normal distribution N(0, sig^2)
  # Covariance between X[i,j] and X[j,i] is sig^2*rho
  
  # X = array(diag(rnorm(n, mean = 0, sd = sig)), dim = c(n, n,itera))
  X = matrix(NaN, n, n)
  
  cov_matrix = sig**2*matrix(c(1, rho, rho, 1), nrow = 2)
  
  random_pairs = MASS::mvrnorm(n = as.integer((n-1)*(n)/2), mu = rep(0, 2), Sigma = cov_matrix)
  
  indices = which(upper.tri(matrix(0, n, n), diag = FALSE), arr.ind = TRUE)
  
  X[upper.tri(X)] = random_pairs[,2]
  X = t(X)
  X[upper.tri(X)] = random_pairs[,1]
  diag(X) = rnorm(n, mean = 0, sd = sig)
  return(X)
}
LV_model = function(Time, N_t, Pars){
  # Lotka-Volterra model
  with(as.list(c(N_t, Pars)),{
    dN_dt = (N_t*(mu_i-A%*%N_t))*(1/(1+exp(-10**7*(N_t-10**(-6)))))
    return(list(c(dN_dt)))
  })
}
get_LV_sim = function(State = runif(n, min = 0, max = 1), mu_i, A, Time = seq(0, 10000, 2), thold) {
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
  condition[nonsurvivors_indices] = "excluded"
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
calc_NFDs = function(n, sig, rho, lambda, mu_i = rep(1, n), Time = seq(0, 10000, 2), thold = 1e-5){
  Id = diag(n) # identity matrix
  J = matrix(1, n, n) # 1 matrix array
  
  X = create_matrix(n, sig, rho)
  A = Id + lambda*J - X
  A = A/(1+lambda)
  
  lambda = lambda/(1+lambda)
  
  ode_densities = get_LV_sim(State = runif(n, min = 0, max = 1), mu_i, A, Time = Time, thold)
  equi = ode_densities$density
  survivors_indices = get_survivors_indices(ode_densities)
  nonsurvivors_indices = get_nonsurvivors_indices(ode_densities)
  
  survs = length(survivors_indices)
  
  solve_densities = get_solve_densities(n, mu_i, A, survivors_indices)
  condition = get_condit(n, solve_densities, survivors_indices, nonsurvivors_indices)
  
  inv_test = get_inv_test(solve_densities, mu_i, A)
  condition = modify_condit_inv(inv_test, condition, nonsurvivors_indices)
  
  sub_equi = get_sub_equi(solve_densities, mu_i, A, survivors_indices, nonsurvivors_indices)
  condition = modify_condit_sub(sub_equi, condition)
  
  r_i = get_r_i(mu_i, A, sub_equi)
  c_ij = get_c_ij(A, n)
  eta_i = get_eta_i(c_ij, mu_i, A, sub_equi)
  ND = get_ND(mu_i, r_i, eta_i)
  FD = get_FD(mu_i, eta_i)
  
  return(vals = list(A = A,
                     solve_densities = solve_densities,
                     equi = equi,
                     survivors_indices = survivors_indices,
                     nonsurvivors_indices = nonsurvivors_indices,
                     inv_test = inv_test,
                     ND = ND,
                     FD = FD,
                     c_ij = c_ij,
                     r_i = r_i,
                     eta_i = eta_i,
                     condition = condition))
}
get_NFD = function(n_itera, lambda_itera, sig_itera, rho_itera, thold = thold, inv = FALSE){
  repeat{
    n = n_itera[replicate]
    mu_i = rep(1, n) # intrinsic growth-rate of the species
    lambda = lambda_itera[replicate]
    sig = sig_itera[replicate]
    rho = rho_itera[replicate]
    gamma = lambda/(1+lambda*n)
    M = replicate
    print(replicate)
    
    data = f_OneMatrix(n)
    
    vals = calc_NFDs(n, sig, rho, lambda, mu_i = rep(1, n), Time = seq(0, 10000, 2), thold = 1e-5)
    
    df_comm_mvfit = data.frame("ND" = vals$ND, "FD" = vals$FD)
    if(length(df_comm_mvfit$ND) < 8){next}
    mvn_result_hz = mvn(data = df_comm_mvfit, mvn_test = "hz") # Henze-Zirkler
    
    if (identical(which(vals$equi > thold), which(vals$solve_densities > 0)) == FALSE || any(vals$solve_densities < 0) == TRUE){
      next
    }
    
    data[,"matrix"] = M
    data[,"n"] = n
    data[,"sigma"] = sig
    data[,"lambda"] = lambda
    data[,"rho"] = rho
    
    data[,"condition"] = vals$condition
    data[,"N_star"] = vals$solve_densities
    
    data[,"c_ij"] = mean(vals$c_ij[row(vals$c_ij) != col(vals$c_ij)])
    data[,"var_c_ij"] = var(c(vals$c_ij), na.rm = TRUE)
    data[,"r_i"] = mean(vals$r_i, na.rm = TRUE)
    data[,"var_r_i"] = var(vals$r_i, na.rm = TRUE)
    data[,"eta_i"] = mean(vals$eta_i, na.rm = TRUE)
    data[,"var_eta_i"] = var(vals$eta_i, na.rm = TRUE)
    
    data[,"ND"] = vals$ND
    data[,"mean_ND"] = mean(vals$ND, na.rm = TRUE)
    data[,"var_ND"] = var(vals$ND, na.rm = TRUE)
    data[,"FD"] = vals$FD
    data[,"mean_FD"] = mean(vals$FD, na.rm = TRUE)
    data[,"var_FD"] = var(vals$FD, na.rm = TRUE)
    data[,"cov_NFD"] = cov(vals$ND, vals$FD, use = "complete.obs")
    data[,"cor_NFD"] = cor(vals$ND, vals$FD, use = "complete.obs")
    
    data[,"pvalue"] = mvn_result_hz$multivariate_normality[,3]
    
    data[data$condition == "resident","proportion"] = (sum(data$condition == "resident")+sum(data$condition == "no_n_1"))/n # Proportion of residents and species with no n-1 community, thus, proportion of extant species in the community
    data[data$condition == "excluded","proportion"] = sum(data$condition == "excluded")/n # Proportion of excluded species that cannot reinvade
    
    survs = length(vals$solve_densities > 0)
    
    data[,"theo_c_ij"] = 1
    data[,"theo_var_c_ij"] = (1/(2*lambda**2) + 1/(2*(1+lambda)**2) - rho/(2*lambda**2))*sig**2
    data[,"theo_r_i"] = (1-lambda)/(1-lambda+n*lambda)
    data[,"theo_var_r_i"] = ((sig**2*n)/(1+lambda*n)**2)*(1+n*sig**2)
    data[,"theo_eta_i"] = (1-survs)/(1+survs*lambda)
    data[,"theo_var_eta_i"] = ((n*sig**2)/(1+lambda*n)**2)*((n/4)+((1+lambda)**2*(1-rho))/(2*lambda**2))
    data[,"theo_ND"] = 1-lambda
    data[,"theo_FD"] = ((survs-1)/survs)*(1-lambda)
    data[,"theo_var_ND"] = (sig**2/(survs*(1+lambda)**2))*(((1+rho)/2)+((survs*lambda**2)/(4*(1+lambda)**2)))
    data[,"theo_var_FD"] = (sig**2/(survs*(1+lambda)**2))*(((1-rho)/2)+((survs*lambda**2)/(4*(1+lambda)**2)))*((1+lambda*n)**2/(n**2+lambda**2))
    data[,"theo_cov_NFD"] = sig**2/n*(1-lambda)**2*1/4*((1-lambda)*lambda+n*lambda**2)
    data[,"theo_cor_NFD"] = 1/sqrt(((1+rho)/2*4/n*1/lambda**2+1)*((1-rho)/2*4/n*1/lambda**2+1))
    
    NFDs = rbind(data, NFDs)
    
    replicate = replicate + 1
    
    if (replicate > itera) break
  }
}

get_Rsquare = function(data = NFDs, y, x){
  m = lm(y ~ x,
         data = data,
         method = "qr")
  return(summary(m)$r.squared)
}