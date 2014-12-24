color.stan <- 'data {
    int<lower=1> N;
    real logRTresidual[N];			// outcome
    real f1[N];					// predictor
    real f2[N];					// predictor
    real f1_X_f2[N];				// interaction
    
    int<lower=1> I;   				// number of subjects
    int<lower=1> J;				// number of items
    int<lower=1,upper=I> subject_id[N];	// subjects
    int<lower=1,upper=J> item_id[N];	// items
    vector[4] subj_prior; 			//vector of zeros passed in from R
    vector[4] item_prior; 			//vector of zeros passed in from R
}

transformed data{
   real ZERO; // like #define ZERO 0 in C/C++   ZERO <- 0.0;
}

parameters{
    real Intercept;
    real beta_f1;
    real beta_f2;
    real beta_f1_X_f2;

    real<lower=0> sigma_e;						// residual sd
    vector[4] u[I];							// random intercept for subjs & 3 slopes for f1, f2, & f12
    vector<lower=0>[4] sigma_u;						// subject sds 
    corr_matrix[4] Omega_subj;					// corr matrix for subjs
    vector[4] w[J];							// random intercept for items & 3 slopes for f1, f2, & f12
    vector<lower=0>[4] sigma_w;						// item sds
    corr_matrix[4] Omega_item;					// corr matrix for items
}




transformed parameters {

	cov_matrix[4] Sigma_subj;	cov_matrix[4] Sigma_item;

    	Sigma_subj <- diag_matrix(sigma_u) * Omega_subj * diag_matrix(sigma_u);
	Sigma_item <- diag_matrix(sigma_w) * Omega_item * diag_matrix(sigma_w);}
model{
    	real mu[N];
	matrix[4,4] L;	matrix[4,4] DL;
	matrix[4,4] M;	matrix[4,4] DM;

   	 // Priors
	Intercept ~ normal(0,10);
    	beta_f1 ~ normal( 0 , 10 );
	beta_f1 ~ normal( 0 , 10 );
    	beta_f1_X_f2 ~ normal( 0 , 10 );

    	sigma_e ~ uniform( 0 , 10 );			// residuals
    	sigma_u ~ uniform( 0 , 10 );			// subj sd
    	sigma_w ~ uniform( 0 , 10 );			// item sd

    	Omega_subj ~ lkj_corr(4.0);
    	Omega_item ~ lkj_corr(4.0);

	L <- cholesky_decompose(Omega_subj);
	M <- cholesky_decompose(Omega_item);

	// create var-cov matrixes for subjects
	for (m in 1:4)		for (n in 1:m)			DL[m,n] <- L[m,n] * sigma_u[m];	for (m in 1:4)		for (n in (m+1):4)			DL[m,n] <- ZERO;

		for (m in 1:4)		for (n in 1:m)			DM[m,n] <- M[m,n] * sigma_w[m];	for (m in 1:4)		for (n in (m+1):4)			DM[m,n] <- ZERO;

    // Varying effects
    for ( i in 1:I ) u[i] ~ multi_normal_cholesky( subj_prior, DL );
    for ( j in 1:J ) w[j] ~ multi_normal_cholesky( item_prior, DM );

    // Fixed effects
    for ( i in 1:N ) {
        mu[i] <- 	Intercept + 
			beta_f1*f1[i] + 
			beta_f2*f2[i] + 
			beta_f1_X_f2*f1_X_f2[i] +
 
			u[subject_id[i],1] + 
			u[subject_id[i],2]*f1[i] + 
			u[subject_id[i],3]*f2[i] + 
			u[subject_id[i],4]*f1_X_f2[i] + 
			w[item_id[i],1] + 
			w[item_id[i],2]*f1[i] + 
			w[item_id[i],3]*f2[i] + 
			w[item_id[i],4]*f1_X_f2[i];


    }
    logRTresidual ~ normal( mu , sigma_e );
}'

spr <- transform(spr, subject.index=as.integer(as.factor(spr$subj)))

r5 <- subset(spr,region==5 & logRTresidual!="NA")


color.dat <- list(subj_prior=c(0,0,0,0),
  item_prior=c(0,0,0,0),  subject_id=sort(as.integer(factor(r5$subject.index))),
  item_id=sort(as.integer(factor(r5$item))),  logRTresidual = r5$logRTresidual,   f1 = ifelse(r5$f1%in%c("0"),-1,1),
  f2 = ifelse(r5$f2%in%c("0"),-1,1),
  f1_X_f2 = ifelse(r5$f1%in%c("0"),-1,1) * ifelse(r5$f2%in%c("0"),-1,1),    N = nrow(r5),  I = length(unique(r5$subject.index)),
  J = length(unique(r5$item))
)

r5.sm <- stan_model(model_code="color.stan",model_name="color")
r5.sf1 <- sampling(r5.sm, data = color.dat,chain_id=1,chains=1,warmup=2500,iter=5000,thin=1,seed=7654321)
r5.sf2 <- sampling(r5.sm, data = color.dat,chain_id=2,chains=1,warmup=2500,iter=5000,thin=1,seed=7654321)
r5.sf3 <- sampling(r5.sm, data = color.dat,chain_id=3,chains=1,warmup=2500,iter=5000,thin=1,seed=7654321)
r5.sf4 <- sampling(r5.sm, data = color.dat,chain_id=4,chains=1,warmup=2500,iter=5000,thin=1,seed=7654321)

------------------------------------------------

r13 <- subset(spr,region==13 & logRTresidual!="NA")



color.dat <- list(subj_prior=c(0,0,0,0),
  item_prior=c(0,0,0,0),  subject_id=sort(as.integer(factor(r13$subject.index))),
  item_id=sort(as.integer(factor(r13$item))),  logRTresidual = r13$logRTresidual,   f1 = ifelse(r13$f1%in%c("0"),-1,1),
  f2 = ifelse(r13$f2%in%c("0"),-1,1),
  f1_X_f2 = ifelse(r13$f1%in%c("0"),-1,1) * ifelse(r13$f2%in%c("0"),-1,1),    N = nrow(r13),  I = length(unique(r13$subject.index)),
  J = length(unique(r13$item))
)

r13.sm <- stan_model(model_code="color.stan",model_name="color")

r13.sf1 <- sampling(r13.sm, data = color.dat,chain_id=1,chains=1,warmup=2500,iter=5000,thin=1,seed=7654321)
r13.sf2 <- sampling(r13.sm, data = color.dat,chain_id=2,chains=1,warmup=2500,iter=5000,thin=1,seed=7654321)
r13.sf3 <- sampling(r13.sm, data = color.dat,chain_id=3,chains=1,warmup=2500,iter=5000,thin=1,seed=7654321)
r13.sf4 <- sampling(r13.sm, data = color.dat,chain_id=4,chains=1,warmup=2500,iter=5000,thin=1,seed=7654321)

