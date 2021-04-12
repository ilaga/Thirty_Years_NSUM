library(NSUM)
library(LaplacesDemon) ## For inverse chi-square distribution
library(rjags) ## For Teo models
library(ggplot2)

## Load data
data(Curitiba)

## Simulate data using NSUM package and basic R functions
# y = as.matrix(read.csv("Curitiba_y.csv", header = F))
# x_cov = as.matrix(read.csv("Curitiba_x_int_cov.csv", header = F))
# x_cov_teo = x_cov[,-1]
y = with(Curitiba, nsum.simulate(500, known, unknown, N, model = "degree",
                             mu, sigma))$y
x_cov_teo = cbind(rnorm(500), rbinom(500, 1, 0.3), rbinom(500, 1, 0.8),
                  rbinom(500, 1, 0.5))
x_cov_teo = scale(x_cov_teo, center = TRUE, scale = FALSE)
x_cov_teo = as.data.frame(x_cov_teo)
names(x_cov_teo) = c("age", "born.in.curitiba", "employed", "gender")
N = Curitiba$N
known = Curitiba$known
n.known = length(known)
N.i = nrow(y)
N.k = ncol(y)




# Killworth 1998b ---------------------------------------------------------
killworth.est = killworth.start(y, known, N)




# Zheng 2006 --------------------------------------------------------------

## Gibbs sampler
N.mc = 10000


prevalences = Curitiba$known / Curitiba$N
pg1.ind = 1:20
Pg1 = sum(prevalences[pg1.ind])

## Declare parameters
alphas = matrix(NA, nrow = N.mc, ncol = N.i)
betas = matrix(NA, nrow = N.mc, ncol = N.k)
omegas = matrix(NA, nrow = N.mc, ncol = N.k)
mu.alpha = mu.beta = sigma.sq.alpha = sigma.sq.beta = rep(NA, N.mc)
C1 = C2 = C = NA

alphas[1,] = 0
betas[1,] = 0
omegas[1,] = 1.5
mu.alpha[1] = mu.beta[1] = sigma.sq.alpha[1] = sigma.sq.beta[1] = 1


for(ind in 2:N.mc){
  ## Step 1
  for(i in 1:N.i){
    alpha.prop = alphas[ind - 1,i] + rnorm(1, 0, 0.4)
    zeta.prop = exp(alpha.prop + betas[ind - 1,]) / (omegas[ind - 1,] - 1)
    zeta.old = exp(alphas[ind - 1, i] + betas[ind - 1,]) / (omegas[ind - 1,] - 1)
    sum1 = sum(lgamma(y[i,] + zeta.prop) - lgamma(zeta.prop) - zeta.prop * log(omegas[ind - 1,]) +
                 dnorm(alpha.prop, mu.alpha[ind - 1], sqrt(sigma.sq.alpha[ind - 1])))
    sum2 = sum(lgamma(y[i,] + zeta.old) - lgamma(zeta.old) - zeta.old * log(omegas[ind - 1,]) +
                 dnorm(alphas[ind - 1,i], mu.alpha[ind - 1], sqrt(sigma.sq.alpha[ind - 1])))
    prob.acc = exp(sum1 - sum2)
    
    if(prob.acc > runif(1)){
      alphas[ind, i] = alpha.prop
    }else{
      alphas[ind, i] = alphas[ind - 1, i]
    }
  }
  
  ## Step 2
  for(k in 1:N.k){
    beta.prop = betas[ind - 1,k] + rnorm(1, 0, 0.05)
    zeta.prop = exp(alphas[ind, ] + beta.prop) / (omegas[ind - 1,k] - 1)
    zeta.old = exp(alphas[ind, ] + betas[ind - 1,k]) / (omegas[ind - 1,k] - 1)
    sum1 = sum(lgamma(y[,k] + zeta.prop) - lgamma(zeta.prop) - zeta.prop * log(omegas[ind - 1,k]) +
                 dnorm(beta.prop, mu.beta[ind - 1], sqrt(sigma.sq.beta[ind - 1])))
    sum2 = sum(lgamma(y[,k] + zeta.old) - lgamma(zeta.old) - zeta.old * log(omegas[ind - 1,k]) +
                 dnorm(betas[ind - 1,k], mu.beta[ind - 1], sqrt(sigma.sq.beta[ind - 1])))
    prob.acc = exp(sum1 - sum2)
    
    if(prob.acc > runif(1)){
      betas[ind, k] = beta.prop
    }else{
      betas[ind, k] = betas[ind - 1, k]
    }
  }
  
  ## Step 3
  mu.alpha.hat = mean(alphas[ind,])
  mu.alpha[ind] = rnorm(1, mu.alpha.hat, sqrt(sigma.sq.alpha[ind - 1] / 2))
  
  ## Step 4
  sigma.alpha.hat = mean((alphas[ind,] - mu.alpha[ind])^2)
  sigma.sq.alpha[ind] = rinvchisq(1, N.i - 1, sigma.alpha.hat)
  
  ## Step 5
  mu.beta.hat = mean(betas[ind,])
  mu.beta[ind] = rnorm(1, mu.beta.hat, sqrt(sigma.sq.beta[ind - 1] / 2))
  
  ## Step 6
  ## Sampling is poor if using inverse-chi squared distribution
  ## Instead put weak prior of sigma.sq.beta
  # sigma.beta.hat = mean((betas[ind,] - mu.beta[ind])^2)
  # sigma.sq.beta[ind] = rinvchisq(1, N.k - 1, sigma.beta.hat)
  sigma.sq.beta[ind] = 10
  
  
  ## Step 7
  for(k in 1:N.k){
    omega.prop = omegas[ind - 1,k] + rnorm(1, 0, 0.1)
    if(omega.prop > 1){
      zeta.prop = exp(alphas[ind, ] + betas[ind,k]) / (omega.prop - 1)
      zeta.old = exp(alphas[ind, ] + betas[ind,k]) / (omegas[ind - 1,k] - 1)
      sum1 = sum(lgamma(y[,k] + zeta.prop) - lgamma(zeta.prop) - zeta.prop * log(omega.prop) +
                   y[,k] * log((omega.prop - 1)/omega.prop))
      sum2 = sum(lgamma(y[,k] + zeta.old) - lgamma(zeta.old) - zeta.old * log(omegas[ind - 1,k]) +
                   y[,k] * log((omegas[ind - 1, k] - 1)/omegas[ind - 1,k]))
      prob.acc = exp(sum1 - sum2)
      
      if(prob.acc > runif(1)){
        omegas[ind, k] = omega.prop
      }else{
        omegas[ind, k] = omegas[ind - 1, k]
      }
    }else{
      omegas[ind,k] = omegas[ind - 1, k]
    }
  }
  
  ## Step 8
  C1 = log(sum(exp(betas[ind, pg1.ind]) / Pg1))
  C = C1
  alphas[ind,] = alphas[ind,] + C
  mu.alpha[ind] = mu.alpha[ind] + C
  betas[ind,] = betas[ind,] - C
  mu.beta[ind] = mu.beta[ind] - C
  
  if(ind %% 1 == 0){
    print(ind)
  }
}

## Burn-in and thin
alphas = alphas[-c(1:1000),]
betas = betas[-c(1:1000),]
alphas = alphas[seq(1, nrow(alphas), by = 10),]
betas = betas[seq(1, nrow(betas), by = 10),]



# Maltiel 2015 ------------------------------------------------------------

## Degree model
maltiel.randeg.est = nsum.mcmc(y, known, N, model = "degree", iterations = 5000, burnin = 2000,
                        size = 2500)
maltiel.randeg.degree.est = rowMeans(maltiel.randeg.est$d.values)
maltiel.randeg.unknown.est = mean(maltiel.randeg.est$NK.values[1,])

## Barrier Effect Model
maltiel.barrier.est = nsum.mcmc(y, known, N, model = "barrier", iterations = 5000, burnin = 2000,
                        size = 2500)
maltiel.barrier.degree.est = rowMeans(maltiel.barrier.est$d.values)
maltiel.barrier.unknown.est = mean(maltiel.barrier.est$NK.values[1,])



# Habecker 2015 -----------------------------------------------------------

habecker.MoS.degree = rep(0, N.i)
for(i in 1:N.i){
  habecker.MoS.degree[i] = mean(y[i,1:n.known] / Curitiba$known) * N
}

habecker.MoS.unknown = mean(y[,N.k] / habecker.MoS.degree) * N




# Teo 2019 ----------------------------------------------------------------

## indexu is unknown subpopulation indices
## indexk is known subpopulation indices
## Sk is known subpopulation sizes

indexu = 21
indexk = 1:20


## No covariates for model1
model1 = 'model {

for(i in 1:N)

{
  
  for(k in 1:Ku)
  
  {
  
  nu[i,k] ~ dpois(lambda*alpha[i]*Su[k])
  
  }
  
  for(k in 1:Kk)
  
  {
  
  nk[i,k] ~ dpois(lambda*alpha[i]*Sk[k])
  
  }
  
  alpha[i]~dlnorm(0,tau)
  
}

for(k in 1:Ku)

{
  
  Su[k]~dunif(0,2500000)
  
}

for(k in 1:Kk)

{
  
  Sk[k]~dunif(0,2500000)
  
}

lambda ~ dunif(0,10)

tau ~ dunif(0,10)

}
'

runModel = function(indexu,indexk,NITERATION)
{
  dataset=list(
    N=dim(y)[1],
    Kk=length(indexk),
    nk=y[,indexk],
    Ku=length(indexu),
    nu=as.matrix(y[,indexu], ncol = length(indexu)),
    Sk=Curitiba$known[indexk],
    Su=rep(NA,length(indexu)))
  
  initialisation=list(lambda=0.1)
  jagmod=jags.model(textConnection(model1),data=dataset,inits=initialisation,n.chains=2)
  update(jagmod, n.iter=5000, progress.bar="text")
  posterior = coda.samples(jagmod, c("alpha","lambda","tau","Su"),n.iter=NITERATION,progress.bar="text",thin=10)
  results = list(indexk=indexk,indexu=indexu,dataset = dataset,posterior=posterior)
  return(results)
}


teo.basic.res = runModel(indexu, indexk, 10000)
teo.basic.post = teo.basic.res$posterior

teo.basic.su = c(teo.basic.post[[1]][,1], teo.basic.post[[2]][,1])
teo.basic.alpha = rbind(teo.basic.post[[1]][,2:501], teo.basic.post[[2]][,2:501])
teo.basic.lambda = c(teo.basic.post[[1]][,502], teo.basic.post[[2]][,502])
teo.basic.tau = c(teo.basic.post[[1]][,503], teo.basic.post[[2]][,503])

## Calculate network size
teo.basic.network.size = rep(NA, N.i)
for(i in 1:N.i){
  teo.basic.network.size[i] = mean(teo.basic.lambda * teo.basic.alpha[,i])*Curitiba$N
}






# Use respondent demographic covariates -----------------------------------

model4 = 'model {
  
for(i in 1:N)

{
  
  for(k in 1:Ku)
  
  {
  
  nu[i,k] ~ dpois(lambda*alpha[i]*exp(b1u[k]*age[i])*exp(b2u[k]*born[i])*exp(b3u[k]*employed[i])*exp(b4u[k]*gender[i])*Su[k])
  
  }
  
  for(k in 1:Kk)
  
  {
  
  nk[i,k] ~ dpois(lambda*alpha[i]*exp(b1k[k]*age[i])*exp(b2k[k]*born[i])*exp(b3k[k]*employed[i])*exp(b4k[k]*gender[i])*Sk[k])
  
  }
  
  alpha[i]~dlnorm(0,tau)
  
}

for(k in 1:Ku)

{
  
  Su[k]~dunif(0,2500000)
  
}

for(k in 1:Kk)

{
  
  Sk[k]~dunif(0,2500000)
  b1k[k]~dnorm(0,(1/sigmab)^2)
  b2k[k]~dnorm(0,(1/sigmab)^2)
  b3k[k]~dnorm(0,(1/sigmab)^2)
  b4k[k]~dnorm(0,(1/sigmab)^2)
}

for(k in 1:Ku)

{
  b1u[k]~dnorm(0,(1/sigmab)^2)
  b2u[k]~dnorm(0,(1/sigmab)^2)
  b3u[k]~dnorm(0,(1/sigmab)^2)
  b4u[k]~dnorm(0,(1/sigmab)^2)
}

lambda ~ dunif(0,10)

tau ~ dunif(0,10)

sigma ~ dgamma(1,0.01)
sigmab ~ dgamma(1,0.01)
}
'


runModel.4 = function(indexu,indexk,NITERATION)
{
  dataset=list(
    N=dim(y)[1],
    Kk=length(indexk),
    nk=y[,indexk],
    Ku=length(indexu),
    nu=as.matrix(y[,indexu], ncol = length(indexu)),
    Sk=Curitiba$known[indexk],
    Su=rep(NA,length(indexu)),
    age=x_cov_teo$age,
    born=x_cov_teo$born.in.curitiba,
    employed=x_cov_teo$employed,
    gender=x_cov_teo$gender)
  
  initialisation=list(lambda=0.1)
  jagmod=jags.model(textConnection(model4),data=dataset,inits=initialisation,n.chains=2)
  update(jagmod, n.iter=5000, progress.bar="text")
  posterior = coda.samples(jagmod, c("alpha","lambda",'sigma','sigmab',"tau","Su","b1u","b2u","b3u",'b4u',"b1k","b2k",'b3k','b4k'),n.iter=NITERATION,progress.bar="text",thin=10)
  results = list(indexk=indexk,indexu=indexu,dataset = dataset,posterior=posterior)
  return(results)
}


teo.barrier.est = runModel.4(indexu, indexk, 10000)
teo.barrier.post = teo.barrier.est$posterior

teo.barrier.su = c(teo.barrier.post[[1]][,1], teo.barrier.post[[2]][,1])
teo.barrier.alpha = rbind(teo.barrier.post[[1]][,2:501], teo.barrier.post[[2]][,2:501])
teo.barrier.bk = rbind(teo.barrier.post[[1]][,502:585], teo.barrier.post[[2]][,502:585])
teo.barrier.lambda = c(teo.barrier.post[[1]][,586], teo.barrier.post[[2]][,586])
teo.barrier.sigma = c(teo.barrier.post[[1]][,587], teo.barrier.post[[2]][,587])
teo.barrier.sigmab = c(teo.barrier.post[[1]][,588], teo.barrier.post[[2]][,588])
teo.barrier.tau = c(teo.barrier.post[[1]][,589], teo.barrier.post[[2]][,589])

## Calculate network size
teo.barrier.network.size = rep(NA, N.i)
for(i in 1:N.i){
  teo.barrier.network.size[i] = mean(teo.barrier.lambda * teo.barrier.alpha[,i])*Curitiba$N
}






# Save data ---------------------------------------------------------------
save.image("Curitiba_Analysis_Results.RData")





# Degree estimates --------------------------------------------------------

killworth.degree.df = data.frame(degree = killworth.est$d.start, model = "Killworth MLE")
habecker.degree.df = data.frame(degree = habecker.MoS.degree, model = "Habecker MoS")
maltiel.randeg.degree.df = data.frame(degree = maltiel.randeg.degree.est, model = "Maltiel Random Degree")
maltiel.barrier.degree.df = data.frame(degree = maltiel.barrier.degree.est, model = "Maltiel Barrier Effects")
zheng.degree.df = data.frame(degree = colMeans(exp(alphas)), model = "Zheng Overdispersed")
teo.basic.degree.df = data.frame(degree = teo.basic.network.size, model = "Teo Basic")
teo.barrier.degree.df = data.frame(degree = teo.barrier.network.size, model = "Teo Barrier Effects")
degree.df = data.frame(rbind(killworth.degree.df, habecker.degree.df,
                             maltiel.randeg.degree.df, maltiel.barrier.degree.df,
                             zheng.degree.df,
                             teo.basic.degree.df, teo.barrier.degree.df))
degree.df$model = factor(degree.df$model, levels = c("Killworth MLE",
                                                     "Zheng Overdispersed",
                                                     "Maltiel Random Degree",
                                                     "Maltiel Barrier Effects",
                                                     "Habecker MoS",
                                                     "Teo Basic",
                                                     "Teo Barrier Effects"))


ggplot(degree.df) +
  geom_density(aes(degree, color = model, group = model), 
               key_glyph = "path", lwd = 1, bw = "bcv") +
  xlim(0, 1000) +
  xlab("Estimated Degree") + ylab("Density") + theme_gray(base_size = 17) +
  guides(color = guide_legend(title = "Model"))


# Unknown subpopulation size estimates ------------------------------------

killworth.est$NK.start
habecker.MoS.unknown
mean(maltiel.unknown.est)
mean(exp(betas[,21]) * N)
mean(teo.basic.su)
mean(teo.barrier.su)

unknown.df = data.frame(est = c(killworth.est$NK.start,
                   habecker.MoS.unknown,
                   mean(maltiel.randeg.est$NK.values),
                   mean(maltiel.barrier.est$NK.values),
                   mean(exp(betas[,21]) * N),
                   mean(teo.basic.su),
                   mean(teo.barrier.su),
                   0.063 * N, 0.07 * N, 0.055 * N),
           lower = c(killworth.est$NK.start - qnorm(0.975) * sqrt(N * killworth.est$NK.start / sum(killworth.est$d.start)),
                     habecker.MoS.unknown - qnorm(0.975) * sqrt(N * habecker.MoS.unknown / N * mean(habecker.MoS.degree)),
                     quantile(maltiel.randeg.est$NK.values, probs = c(0.025)),
                     quantile(maltiel.barrier.est$NK.values, probs = c(0.025)),
                     quantile(exp(betas[,21]) * N, probs = c(0.025)),
                     quantile(teo.basic.su, probs = c(0.025)),
                     quantile(teo.barrier.su, probs = c(0.025)),
                     0.045 * N, 0.055 * N, 0.048 * N),
           upper = c(killworth.est$NK.start + qnorm(0.975) * sqrt(N * killworth.est$NK.start / sum(killworth.est$d.start)),
                     habecker.MoS.unknown + qnorm(0.975) * sqrt(N * habecker.MoS.unknown / N * mean(habecker.MoS.degree)),
                     quantile(maltiel.randeg.est$NK.values, probs = c(0.975)),
                     quantile(maltiel.barrier.est$NK.values, probs = c(0.975)),
                     quantile(exp(betas[,21]) * N, probs = c(0.975)),
                     quantile(teo.basic.su, probs = c(0.975)),
                     quantile(teo.barrier.su, probs = c(0.975)),
                     0.08 * N, 0.085 * N, 0.062 * N),
           model = c("Killworth MLE",
                     "Habecker MoS",
                     "Maltiel Random Degree",
                     "Maltiel Barrier Effects",
                     "Zheng Overdispersed",
                     "Teo Basic",
                     "Teo Barrier Effects",
                     "Feehan Generalized", "Maltiel Transmission", "Maltiel Combined"),
           coded = c(rep("No", 7), rep("Yes", 3)))
unknown.df$model = factor(unknown.df$model, levels = c("Maltiel Combined", "Maltiel Transmission", "Feehan Generalized",
                                                       "Teo Barrier Effects",
                                                       "Teo Basic",
                                                       "Habecker MoS",
                                                       "Maltiel Barrier Effects",
                                                       "Maltiel Random Degree",
                                                       "Zheng Overdispersed",
                                                       "Killworth MLE"))
unknown.df$coded = factor(unknown.df$coded, levels = c("No", "Yes"))


unknown.df$est = unknown.df$est / N
unknown.df$lower = unknown.df$lower / N
unknown.df$upper = unknown.df$upper / N

ggplot(unknown.df, group = model, aes(x = model, y = est, col = factor(coded))) +
  geom_point(aes(x = model, y = est)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, lwd = 1) +
  coord_flip() +
  ylab("Estimated Prevalence in the General Population") +
  xlab("\n") +
  theme_grey(base_size = 17) +
  labs(color = "From Literature") +
  scale_color_manual(labels = c("No", "Yes"), values = c(1, 2))