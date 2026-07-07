## Data generation function for Joint Multi-Layer LTA simulation

sim_JointLTA_data <- function(N, L, K, Time, dk, R, br, pi, tau, item_p, beta, type, J = NULL){

  ### N: number of individuals
  ### L: number of latent classes
  ### K: number of items
  ### Time: number of time points
  ### dk: number of responses for each item
  ### R: number of sequence clusters
  ### br: sequence cluster membership probabilities
  ### pi: item-response latent class membership probabilities at time=1
  ### tau: transition matrix
  ### item_p: item-response probabilities of the K items
  ### beta: outcome model coefficients
  ### type: type of outcome variable, "binary" or "multinomial"
  ### J: number of categories for multinomial outcome variable

  ## simulate sequence cluster memberships for each individual
  lambda <- sample(1:R, size = N, replace = TRUE, prob = br)
  print(prop.table(table(lambda)))

  ## simulate item-response latent class memberships over time
  eta_over_time <- matrix(nrow = N, ncol = Time)
  for(i in 1:N){
    eta_over_time[i,1] <- sample(1:L, size = 1, prob = pi[lambda[i],])
  }
  for(t in 2:Time){
    for(i in 1:N){
      eta_over_time[i,t] <- sample(1:L, size = 1,
                                   prob = tau[lambda[i],eta_over_time[i,t-1],])
    }
  }

  eta_over_time <- as.data.frame(eta_over_time)
  eta_over_time$lambda <- lambda
  colnames(eta_over_time)[1:Time] <- sapply(1:Time, FUN = function(t) paste0("Visit", t))

  ## check latent class memberships distribution at time1
  print(eta_over_time %>%
          count(lambda, Visit1) %>%
          group_by(lambda) %>%
          mutate(prop = n/sum(n)) %>%
          arrange(lambda, Visit1))
  
  ## simulate item-responses over time
  item <- array(dim = c(Time, N, K))
  for(t in 1:Time){
    for(i in 1:N){
      for(k in 1:K){
        item[t,i,k] <- sample(1:dk[k], size = 1, replace = TRUE,
                              prob = item_p[eta_over_time[i,t],k,1:dk[k]])
      }
    }
  }

  ## randomly pick a few responses (0-10) for each individual visit to be missing
  for(t in 1:Time){
    for(i in 1:N){
      item[t,i,sample(1:K, sample(0:10, size = 1), replace = FALSE)] <- NA
    }
  }

  ## randomly pick a few visits (0-4) for each individual to be missing
  having_visits <- matrix(1, nrow = N, ncol = Time)
  for(i in 1:N){
    having_visits[i,sample(1:Time, size = sample(0:4, size = 1,
                                                 prob = c(0.4, 0.2, 0.2, 0.1, 0.1)),
                           replace = FALSE)] <- 0
  }

  ## set the item-responses to be NAs for missing visits
  for(i in 1:N){
    item[which(having_visits[i,]==0),i,] <- NA
  }

  ## simulate outcomes
  if(type == "binary"){
    ## simulate binary outcomes based on the latent class memberships at time=10 and sequence cluster memberships
    logit <- beta[1]+beta[2]*(eta_over_time$lambda==2)+
      beta[3]*(eta_over_time$lambda==3)+beta[4]*(eta_over_time$Visit10==2)+
      beta[5]*(eta_over_time$Visit10==3)+beta[6]*(eta_over_time$Visit10==4)
    prob <- exp(logit)/(1+exp(logit))
    outcome <- rbinom(n = N, size = 1, prob = prob)
  } else {
    ## simulate multinomial outcomes based on the latent class memberships at time=10 and sequence cluster memberships
    outcome <- c()
    for(i in 1:N){
      odds_i <- c()
      for(j in 2:J){
        odds_ij <- exp(beta[[j-1]][1]+beta[[j-1]][2]*(eta_over_time[i,Time]==2)+beta[[j-1]][3]*(eta_over_time[i,Time]==3)+
          beta[[j-1]][4]*(eta_over_time[i,Time]==4)+beta[[j-1]][5]*(lambda[i]==2)+beta[[j-1]][6]*(lambda[i]==3))
        odds_i <- c(odds_i, odds_ij)
      }
      pi_i <- c()
      pi_i0 <- 1/(1+sum(odds_i))
      pi_i <- c(pi_i, pi_i0)
      for(j in 2:J){
        pi_i <- c(pi_i, pi_i0*odds_i[j-1])
      }
      outcome <- c(outcome, sample(1:J, size = 1, prob = pi_i))
    }
  }

  return(list(eta_over_time = eta_over_time, item = item, outcome = outcome, having_visits = having_visits))
}
