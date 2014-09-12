data {
    int<lower=0> N;
    real lrt[N];                     //outcome
    real dist[N];                     //predictor
    real rctype[N];                     //predictor
    real interaction[N];                     //predictor
    int<lower=1> I;                 //number of subjects
    int<lower=1> K;                 //number of items
    int<lower=1, upper=I> subj[N];    //subject id
    int<lower=1, upper=K> item[N];    //item id
    vector[4] mu_prior;             //vector of zeros passed in from R
}
transformed data {
    real ZERO;                      // like #define ZERO 0 in C/C++
    ZERO <- 0.0;
}
parameters {
    vector[4] beta;                 // intercept and slope
    vector[4] u[I];                 // random intercept and slopes subj
    vector[4] w[K];
    real<lower=0> sigma_e;          // residual sd
    vector<lower=0>[4] sigma_u;     // subj sd
    vector<lower=0>[4] sigma_w;     // item sd
    corr_matrix[4] Omega_u;           // correlation matrix for random intercepts and slopes subj
    corr_matrix[4] Omega_w;           // correlation matrix for random intercepts and slopes item
}
transformed parameters {
    matrix[4,4] D_u;
    matrix[4,4] D_w;
    D_u <- diag_matrix(sigma_u);
    D_w <- diag_matrix(sigma_w);
}
model {
    matrix[4,4] L_u;
    matrix[4,4] DL_u;
    matrix[4,4] L_w;
    matrix[4,4] DL_w;
    real mu[N]; // mu for likelihood
    //priors:
    beta ~ normal(0,10);
    sigma_e ~ normal(0,10);
    sigma_u ~ normal(0,10);
    sigma_w ~ normal(0,10);
    Omega_u ~ lkj_corr(4.0);
    Omega_w ~ lkj_corr(4.0);
    L_u <- cholesky_decompose(Omega_u);
    L_w <- cholesky_decompose(Omega_w);
    for (m in 1:4) {
        for (n in 1:m) {
            DL_u[m,n] <- L_u[m,n] * sigma_u[m];
        }
    }
    for (m in 1:4){
        for (n in (m+1):4){
            DL_u[m,n] <- ZERO;
        }
    }
    for (m in 1:4){
        for (n in 1:m){
            DL_w[m,n] <- L_w[m,n] * sigma_w[m];
        }
    }
    for (m in 1:4){
        for (n in (m+1):4){
            DL_w[m,n] <- ZERO;
        }
    }
    for (i in 1:I)                  // loop for subj random effects
        u[i] ~ multi_normal_cholesky(mu_prior, DL_u);
    for (k in 1:K)                  // loop for item random effects
        w[k] ~ multi_normal_cholesky(mu_prior, DL_w);    
    for (n in 1:N) {
        mu[n] <- beta[1] + beta[2]*dist[n] + beta[3]*rctype[n] + beta[4]*interaction[n] + u[subj[n], 1] + u[subj[n], 2]*dist[n] + u[subj[n], 3]*rctype[n] + u[subj[n], 4]*interaction[n] + w[item[n], 1] + w[item[n], 2]*dist[n] + w[item[n], 3]*rctype[n] + w[item[n], 4]*interaction[n];
    }
    lrt ~ normal(mu,sigma_e);        // likelihood
}
generated quantities {
    cov_matrix[4] Sigma_u;
    cov_matrix[4] Sigma_w;
    Sigma_u <- D_u * Omega_u * D_u;
    Sigma_w <- D_w * Omega_w * D_w;
}

