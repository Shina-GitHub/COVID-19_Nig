functions {
  real[] dz_dt(real t, real[] z, real[] theta,
             real[] x_r, int[] x_i) {
      real s = z[1];
      real e1 = z[2];
      real e2 = z[3];
      real i1 = z[4];
      real i2 = z[5];
      real i3 = z[6];
      real r =  z[7];
      real d =  z[8];
      real c =  z[9];
      real v =  z[10];
      
      real beta = theta[1];
      real sigma1 = 0.3125;
      real sigma2 = 0.5;
      real rho1 = 0.5;
      real rho2 = theta[2];
      real rho3 = theta[3];
      real gamma1 = 0.5;
      real gamma2 = 0.176;
      real gamma3 = 0.5;
      real k1 = theta[4];
      //real k2 = theta[5];
      real n1  = 198000000;
      
      real ds_dt = - beta*(i1+2*i2+3*i3)*s/n1;
      real de1_dt = beta*(i1+2*i2+3*i3)*s/n1 - sigma1*e1;
      real de2_dt = sigma1*e1 - sigma2*e2;
      real di1_dt = rho1*sigma2*e2 - gamma1*i1;
      real di2_dt = k1*c + gamma1*i1-gamma2*i2;
      real di3_dt = (1-rho1)*sigma2*e2 - gamma3*i3;
      real dr_dt = rho2*gamma2*i2 + gamma3*i3;
      real dd_dt = (1-rho2)*gamma2*i2;
      real dc_dt = k1*c;
      real dv_dt = rho3*(k1*c + gamma1*i1);
      return { ds_dt,de1_dt,de2_dt,di1_dt,di2_dt,di3_dt,dr_dt,dd_dt,dc_dt,dv_dt };
    }
}

data {
  int N;           // num measurements
  real ts[N];                 // measurement times > 0
  real y_me_init[10];             // initial measured population
  real y_me[N, 3];    // measured population at measurement times [c,v,d]
}

parameters {
  real <lower=0,upper=1> theta[4];   // theta = {alpha, beta, gamma, delta}
  //real z_init[10];  // initial population
  real <lower=0> sigma[3];   // error scale
}

transformed parameters {
  real z[N, 10]
    = integrate_ode_rk45(dz_dt, y_me_init, 0, ts, theta,
                         rep_array(0.0, 0), rep_array(0, 0),
                         1e-6, 1e-5, 1e3);
}

model {
  //theta[{1, 3}] ~ normal(1, 0.5);
  //theta[{2, 4}] ~ normal(0.05, 0.05);
  //sigma ~ normal(0,1);
  //sigma ~ cauchy(0,1);
  //sigma[{1}] ~ normal(0,1);
 // sigma[{2,3}] ~ normal(0,10);
  //sigma ~ lognormal(-1, 1);
  //z_init ~ lognormal(log(10), 1);
  for (k in 1:3) {
    y_me[ , k] ~ normal(z[, (k+7)], sigma[k]);
  }
}
