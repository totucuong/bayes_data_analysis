
data{
    int<lower=1> N;
    int notes[N];
    int cat[N];
}

parameters{
    real<lower=0> alpha; //rate of singing when cat is absent
    real<lower=0> beta;  // rate of singing when cate is present
}

model{
    vector[N] lambda;
    beta ~ exponential( 0.1 );
    alpha ~ exponential( 0.1 );
    for ( i in 1:N ) {
        lambda[i] = (1 - cat[i]) * alpha + cat[i] * beta;
    }
    notes ~ poisson( lambda );
}
generated quantities{
    vector[N] lambda;
    for ( i in 1:N ) {
        lambda[i] = (1 - cat[i]) * alpha + cat[i] * beta;
    }

}