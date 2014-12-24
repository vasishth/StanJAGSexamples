univ_stan <- 'data {
    int<lower=1> N;
    real logRTresidual[N];			// outcome
    real f1[N];					// predictor
    real f2[N];					// predictor
    real f1_X_f2[N];				// interaction
    int<lower=1> I;   				// number of subjects
    int<lower=1> J;				// number of items
    int<lower=1,upper=I> subject_id[N];	// subjects
    int<lower=1,upper=J> item_id[N];	// items
    vector[4] mu_prior; 			//vector of zeros passed in from R

}

transformed data{
   real ZERO; // like #define ZERO 0 in C/C++
   ZERO <- 0.0;
}

parameters{
    real Intercept;
    real beta_f1;
    real beta_f2;
    real beta_f1_X_f2;
    vector[4] u[I];							// random intercept for subjs & 3 slopes for f1, f2, & f12
    vector[4] w[J];							// random intercept for items & 3 slopes for f1, f2, & f12
    real<lower=0> sigma_e;						// residual sd
    vector<lower=0>[4] sigma_u;						// subject sds 
    vector<lower=0>[4] sigma_w;						// item sds
    corr_matrix[4] Omega_subj;					// corr matrix for subjs
    corr_matrix[4] Omega_item;					// corr matrix for items
}

transformed parameters {

	cov_matrix[4] Sigma_subj;
	cov_matrix[4] Sigma_item;

    	Sigma_subj <- diag_matrix(sigma_u) * Omega_subj * diag_matrix(sigma_u);
	Sigma_item <- diag_matrix(sigma_w) * Omega_item * diag_matrix(sigma_w);
}
model{
    	real mu[N];
	matrix[4,4] L;
	matrix[4,4] DL;
	matrix[4,4] M;
	matrix[4,4] DM;

   	 // Priors
	Intercept ~ normal(0,10);
    	beta_f1 ~ normal( 0 , 10 );
	beta_f1 ~ normal( 0 , 10);
    	beta_f1_X_f2 ~ normal( 0 , 10 );

    	sigma_e ~ uniform( 0 , 10 );			// residuals
    	sigma_u ~ uniform( 0 , 10 );			// subj sd
    	sigma_w ~ uniform( 0 , 10 );			// item sd

    	Omega_subj ~ lkj_corr(2.0);
    	Omega_item ~ lkj_corr(2.0);

	L <- cholesky_decompose(Omega_subj);
	M <- cholesky_decompose(Omega_item);

	for (m in 1:4)
		for (n in 1:m)
			DL[m,n] <- L[m,n] * sigma_u[m];
	for (m in 1:4)
		for (n in (m+1):4)
			DL[m,n] <- ZERO;

	
	for (m in 1:4)
		for (n in 1:m)
			DM[m,n] <- M[m,n] * sigma_w[m];
	for (m in 1:4)
		for (n in (m+1):4)
			DM[m,n] <- ZERO;

    // Varying effects
    for ( i in 1:I ) u[i] ~ multi_normal_cholesky( mu_prior, DL );
    for ( j in 1:J ) w[j] ~ multi_normal_cholesky( mu_prior, DM );

    // Fixed effects
    for ( i in 1:N ) {
        mu[i] <- Intercept + beta_f1*f1[i] + beta_f2*f2[i] + beta_f1_X_f2*f1_X_f2[i] + u[subject_id[i],1] + 
			u[subject_id[i],2]*f1[i] + u[subject_id[i],3]*f2[i] + u[subject_id[i],4]*f1_X_f2[i] + 
			w[item_id[i],1] + w[item_id[i],2]*f1[i] + w[item_id[i],3]*f2[i] + w[item_id[i],4]*f1_X_f2[i];


    }
    logRTresidual ~ normal( mu , sigma_e );
}'

## load preprocessed data Expt 1:
load("univ3.base.RData")
## data: r0, r3, r4, r7, r8, r9


######## Prep data - region 9 - similar for other regions ##############

r9 <- subset(univ,region==9 & logRTresidual!="NA")

head(r9)

summary(r9)

r9$f1[r9$cond=="complex_complex"] <- 1
r9$f1[r9$cond=="complex_simple"] <- 1
r9$f1[r9$cond=="simple_complex"] <- 0
r9$f1[r9$cond=="simple_simple"] <- 0

r9$f2[r9$cond=="complex_complex"] <- 1
r9$f2[r9$cond=="complex_simple"] <- 0
r9$f2[r9$cond=="simple_complex"] <- 1
r9$f2[r9$cond=="simple_simple"] <- 0



univ_dat <- list(mu_prior=c(0,0,0,0),
  subject_id=as.integer(factor(r9$subj)),
  item_id=as.integer(factor(r9$item)),
  logRTresidual = r9$logRTresidual, 
  f1 = ifelse(r9$f1%in%c("0"),-1,1),
  f2 = ifelse(r9$f2%in%c("0"),-1,1),
  f1_X_f2 = ifelse(r9$f1%in%c("0"),-1,1) * ifelse(r9$f2%in%c("0"),-1,1),  
  N = nrow(r9),
  I = length(unique(r9$subject.index)),
  J = length(unique(r9$item))
)


##### Compile model

library(rstan)



mod <- stan(model_code = univ_stan,
            data = univ_dat,
            iter = 1000, warmup=500, chains = 4,
            pars=c("beta_f1","beta_f2","beta_f1_X_f2",
                   "sigma_e","sigma_u",
                   "sigma_w"))

print(mod)

MCMCsamp<-as.matrix(mod)

## check post. distribution:
hist(MCMCsamp[,1])
hist(MCMCsamp[,2])
hist(MCMCsamp[,3])

mean(MCMCsamp[,1]<0)
mean(MCMCsamp[,2]<0)
mean(MCMCsamp[,3]<0)






