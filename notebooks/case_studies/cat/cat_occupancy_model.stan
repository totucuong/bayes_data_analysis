
data {
    int N;
    int notes[N];
    int cat[N];
}

parameters {
    real<lower=0> beta;
    real<lower=0> alpha;
    real<lower=0,upper=1> kappa;    // prob cat present
    real<lower=0,upper=1> delta;    // prob of detecting cat
}

model {
    beta ~ exponential( 0.1 );
    alpha ~ exponential( 0.1 );
    kappa ~ beta(4,4);
    delta ~ beta(4,4);
    for ( i in 1:N ) {
        if ( cat[i]==1 )
            // cat present and detected
            target += log(kappa) + log(delta) + poisson_lpmf( notes[i] | beta );
        if ( cat[i]==0 ) {
            // cat not observed, but cannot be sure not there
            // marginalize over unknown cat state:
            // (1) cat present and not detected
            // (2) cat absent
            target += log_sum_exp(
                    log(kappa) + log1m(delta) + poisson_lpmf( notes[i] | beta ),
                    log1m(kappa) + poisson_lpmf( notes[i] | alpha ) );
        }
    }
}