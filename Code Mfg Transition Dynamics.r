
# generating synthetic data
set.seed(123)

n=1000
k=650

# Gamma realizations by Transformation Method
# alpha-----> positive integer
# beta------> positive number

alpha=2
beta=0.1
Gam=c()
for(i in 1:k){
  u1=runif(alpha)
  Gam[i]=(1/beta)*sum(-log(u1))
}
hist(Gam)

# Exponential realizations by Inverse method
lam=0.5
u2=runif(n-k)
Exp=-(1/lam)*log(1-u2)
hist(Exp)

# Joining the two sets of synthetic data

x=c(Gam,Exp)

x=(x-min(x))/(max(x)-min(x))  #normalizing or scaling the data

hist(x,breaks=30)

  

# Single Component MH
SCMH=function(nsim){

start=Sys.time()
# initial parameters of alpha,beta,lam priors
a=3 ; b=1
a1=a; a2=a
b1=b; b2=b
lam1=1

# current values of the parameters
alpha_cur = 2
beta_cur=1
lam_cur=0.5
k_cur=floor(n/2)

alpha_chain=numeric(nsim)
beta_chain=numeric(nsim)
lam_chain=numeric(nsim)
k_chain=numeric(nsim)

#log-posterior calculation
log_posterior=function(alpha,beta,lam,k){
  if(alpha<=0 || beta<=0 || lam<=0 || k<=0 || k>=n){
    return(-Inf)
  } #return negative inf in order to prevent invalid values
  
  log_prior_alpha= a1*log(b1) - log(gamma(a1)) + (a1-1)*log(alpha) - b1*alpha
  log_prior_beta= a2*log(b2) - log(gamma(a2)) + (a2-1)*log(beta) - b2*beta
  log_prior_lam= log(lam1) - (lam1*lam)
  log_prior_k= -log(n-1) #this quantity will not used in rho calculation
  
  #likelihood calculations
  x_gamma=x[1:k]    #first k points (Gamma)
  x_exp=x[(k+1):n]  #remaining n-k points (Exp)
  
  f_gamma=(beta^alpha)*(x_gamma^(alpha-1))*exp(-beta*x_gamma)/gamma(alpha)
  f_exp = lam*exp(-lam*x_exp)
  
  log_L=sum(log(f_gamma*k/n)) + sum(log(f_exp*(n-k)/n))
  
  log_posterior=log_L+log_prior_alpha+log_prior_beta+log_prior_lam + log_prior_k
  
  return(log_posterior)
}

#Single component MH update

for(j in 1:nsim){
  k_new=sample(1:(n-1),1)
  alpha_new=abs(rnorm(1,mean=alpha_cur,sd=0.5))
  beta_new=abs(rnorm(1,mean=beta_cur,sd=0.5))
  lam_new=abs(rnorm(1,mean=lam_cur,sd=0.5))
  
  log_rho=log_posterior(alpha_new,beta_new,lam_new,k_new)-log_posterior(alpha_cur,beta_cur,lam_cur,k_cur)
  
  if(is.nan(log_rho)|is.na(log_rho)){next}
  
  #accept or reject the new value for each component
  rho=min(exp(log_rho),1)
  
  u=runif(1)
  if(u<rho){
    alpha_cur=alpha_new
    beta_cur=beta_new
    lam_cur=lam_new
    k_cur=k_new
  }
  
  alpha_chain[j]=alpha_cur
  beta_chain[j]=beta_cur
  lam_chain[j]=lam_cur
  k_chain[j]=k_cur
  
}

#estimating mean of each component
alpha_mean=mean(alpha_chain)
beta_mean=mean(beta_chain)
lam_mean=mean(lam_chain)
k_mean=mean(k_chain)

end=Sys.time()
runtime=end-start

return(list("mean of alpha"=alpha_mean,"mean of beta"=beta_mean,"mean of lam"=lam_mean,"mean of k"=k_mean,"runtime"=runtime))
}

SCMH(500)
SCMH(1000)
SCMH(5000)
SCMH(10000)




# Independent MH
INDMH=function(nsim){

start=Sys.time()
# initializing parameters of alpha,beta,lam
a=3 ; b=1
a1=a; a2=a; a3=a
b1=b; b2=b; b3=b

# log-likelihood calculation
log_L=function(x,alpha,beta,lam,k){
  if(alpha<=0 || beta<=0 || lam<=0 || k<=0 || k>=n){
    return(-Inf)} #return negative inf in order to prevent invalid values
  gamma_L=sum((alpha-1)*log(x[1:k])- x[1:k]*beta -log(gamma(alpha))+ alpha*log(beta))
  exp_L=sum(log(lam) - lam*x[(k+1):n])
  return(gamma_L+exp_L)
}

# log-prior calculation
log_prior=function(alpha,beta,lam,k){
  log_prior_alpha=(a1-1)*log(alpha) - alpha*b1 - log(gamma(a1)) + a1*log(b1)
  log_prior_beta=(a2-1)*log(beta) - beta*b2 - log(gamma(a2)) + a2*log(b2)
  log_prior_lam=(a3-1)*log(lam) - lam*b3 - log(gamma(a3)) + a3*log(b3)
  log_prior_k=log(1/n)
  return(log_prior_alpha+ log_prior_beta + log_prior_lam + log_prior_k)
}

# custom gamma proposal generation function
generate_gamma=function(alpha,beta){
  u=runif(alpha)
  return((1/beta)*sum(-log(u)))
}

# parameters for independent MH algorithm
alpha=1; beta=1; lam=1; k=n/2

# store output
output=matrix(NA,nsim,4)

for(i in 1:nsim){
  # generate from proposal distributions
  alpha_proposal=generate_gamma(1,1)
  beta_proposal=generate_gamma(1,1)
  lam_proposal=generate_gamma(1,1)
  k_proposal=sample(1:n,1)
  
  # compute log(rho)
  log_rho=log_L(x,alpha_proposal,beta_proposal,lam_proposal,k_proposal) 
          + log_prior(alpha_proposal,beta_proposal,lam_proposal,k_proposal) 
          - log_L(x,alpha,beta,lam,k) - log_prior(alpha,beta,lam,k)
  
  if(!is.na(log_rho) && log(runif(1)) < log_rho){
    alpha=alpha_proposal
    beta=beta_proposal
    lam=lam_proposal
    k=k_proposal
  }
  output[i,]=c(alpha,beta,lam,k)
}

alpha_mean=mean(alpha)
beta_mean=mean(beta)
lam_mean=mean(lam)
k_mean=mean(k)

end=Sys.time()
runtime=end-start

return(list("mean of alpha"=alpha_mean,"mean of beta"=beta_mean,"mean of lam"=lam_mean,"mean of k"=k_mean,"runtime"=runtime))
}

INDMH(500)
INDMH(1000)
INDMH(5000)
INDMH(10000)




# Table output

# True values
true_alpha <- 2
true_beta <- 0.1
true_lambda <- 0.5
true_k <- 650

# Data: SCMH and INDMH results
results <- data.frame(
  Method = c("SCMH", "INDMH", "SCMH", "INDMH", "SCMH", "INDMH", "SCMH", "INDMH"),
  Iterations = c(500, 500, 1000, 1000, 5000, 5000, 10000, 10000),
  Alpha = c(0.8361905, 0.5049638, 0.2573558, 0.2078914, 0.1357002, 1.717153, 0.716214, 0.1613593),
  Beta = c(0.8683006, 0.9283596, 1.594198, 0.7464388, 0.7058617, 1.152869, 1.02531, 0.799877),
  Lambda = c(0.8359284, 0.4304201, 1.23521, 2.090191, 0.8530741, 2.892499, 1.398508, 0.02134403),
  K = c(785.688, 916, 640.073, 394, 674.5192, 22, 633.9616, 725),
  Runtime = c(0.2039139, 0.05020595, 0.3548911, 0.08906603, 1.764668, 0.412102, 3.33606, 0.8283839)
)

# Calculate relative errors
results$RE_Alpha <- abs(results$Alpha - true_alpha) / true_alpha * 100
results$RE_Beta <- abs(results$Beta - true_beta) / true_beta * 100
results$RE_Lambda <- abs(results$Lambda - true_lambda) / true_lambda * 100
results$RE_K <- abs(results$K - true_k) / true_k * 100

# Select and format columns for the table
table_data <- results[, c("Method", "Iterations", "RE_Alpha", "RE_Beta", "RE_Lambda", "RE_K", "Runtime")]

# Output the table using kable to render it as a markdown table
knitr::kable(table_data, format = "markdown", digits = 2, 
             col.names = c("Method", "Iterations", "RE_Alpha (%)", "RE_Beta (%)", "RE_Lambda (%)", "RE_K (%)", "Runtime (s)"),
             caption = "Relative Error and Runtime for Single Component MH and Independent MH Methods")







