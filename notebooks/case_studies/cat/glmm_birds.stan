
data{
    int<lower=1> N;
    int<lower=1> N_id;
    int notes[N];
    int cat[N];
    int id[N];
}

parameters{
    vector<lower=0>[N_id] alpha;
    vector<lower=0>[N_id] beta;
    real<lower=0> alpha_bar;
    real<lower=0> beta_bar;
}

model{
    vector[N] lambda;
    beta_bar ~ exponential( 0.1 );
    alpha_bar ~ exponential( 0.1 );
    for ( j in 1:N_id ) beta[j] ~ exponential( 1.0/beta_bar );
    for ( j in 1:N_id ) alpha[j] ~ exponential( 1.0/alpha_bar );
    for ( i in 1:N ) {
        lambda[i] = (1 - cat[i]) * alpha[id[i]] + cat[i] * beta[id[i]];
    }
    notes ~ poisson( lambda );
}

generated quantities{
    vector[N] lambda;
    for ( i in 1:N ) {
        lambda[i] = (1 - cat[i]) * alpha[id[i]] + cat[i] * beta[id[i]];
    }

}