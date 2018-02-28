
data{
    int N;
    int notes[N];
    int cat[N];
}

parameters{
    real<lower=0,upper=1> kappa;
    real<lower=0> beta;
    real<lower=0> alpha;
}

model{
    beta ~ exponential( 0.1 );
    alpha ~ exponential( 0.1 );
    kappa ~ beta(4,4);
    for ( i in 1:N ) {
        if ( cat[i]==-1 ) { // cat missing
            target += log_mix( kappa ,
                    poisson_lpmf( notes[i] | beta ),
                    poisson_lpmf( notes[i] | alpha )
                );
        } else { // cat not missing
            cat[i] ~ bernoulli(kappa);
            notes[i] ~ poisson( (1-cat[i])*alpha + cat[i]*beta );
        }
    }
}

generated quantities{
    vector[N] cat_impute;
    for ( i in 1:N ) {
        real logPxy;
        real logPy;
        if ( cat[i]==-1 ) {
            // need P(cat=1|notes)
            // P(cat=1|notes) = P(cat=1,notes)/P(notes)
            // P(cat=1,notes) = P(cat=1)P(notes|cat=1 or rate = beta)
            // P(notes) = P(cat==1)P(notes|cat==1) + P(cat==0)P(notes|cat==0)
            logPxy = log(kappa) + poisson_lpmf(notes[i]|beta);
            logPy = log_mix( kappa ,
                    poisson_lpmf( notes[i] | beta ),
                    poisson_lpmf( notes[i] | alpha ) );
            cat_impute[i] = exp( logPxy - logPy );
        } else {
            cat_impute[i] = cat[i];
        }
    }
}