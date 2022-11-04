//
// This Stan program defines the copula model for competing risks.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//


// defined functions for calculating likelihood functions

functions {
// baseline hazard functions as piecewise constant functions
   real lambda0U(real u, real[] LU, real[] CPU) {
     int i;
     int num;
     i=1;
     num=1;
     while (i < (num_elements(CPU)+1)){
       if (u>CPU[num_elements(CPU)-i+1]) {
         num=num_elements(CPU)-i+1;
         i=num_elements(CPU)+100;
       }else {
         i=i+1;
       }
     }
     return(LU[num]);
   }



// marginal suvival functions
  real SU(real u, vector x, vector beta1, real[] LU, real[] CPU) {
    int i;
    int num;
    real CH;
    real su;
    i=1;
    num=1;
    CH=0;
    while (i < num_elements(CPU)){
      if (u > CPU[i+1]){
        CH=CH+(CPU[i+1]-CPU[i])*LU[i];
        i=i+1;
        num=i;
      }else {
        i=num_elements(CPU)+100;
      }
    }
    CH=CH+(u-CPU[num])*LU[num];
    su=exp(-CH)^(exp(sum(beta1 .* x)));
    if(su>0.001){
      return(su);
    }else{
      return(0.001);
    }
  }



// marginal CDFs
  real FU(real u, vector x, vector beta1, real[] LU, real[] CPU) {
    return 1-SU(u, x, beta1, LU, CPU);
  }

  
// marginal pdfs
  real fU(real u, vector x, vector beta1, real[] LU, real[] CPU) {
    return lambda0U(u, LU, CPU)*exp(sum(beta1 .* x))*SU(u, x, beta1, LU, CPU);
  }

// joint survival function
  real SUV(real u, real v, vector x, vector beta1, vector beta2, real theta, real[] LU, real[] LV, real[] CPU, real[] CPV) {
    return (SU(u, x, beta1, LU, CPU)^(-theta)+SU(v, x, beta2, LV, CPV)^(-theta)-1)^(-1/theta);
  }




// sub-density functions
  real fT1(real t, vector x, vector beta1, vector beta2, real theta, real[] LU, real[] LV, real[] CPU, real[] CPV) {
    return (SU(t,x, beta1, LU, CPU)^(-theta)+SU(t,x, beta2, LV, CPV)^(-theta)-1)^(-1/theta-1)*SU(t,x, beta1, LU, CPU)^(-theta-1)*fU(t, x, beta1, LU, CPU);
  }


  real fT2(real t, vector x, vector beta1, vector beta2, real theta, real[] LU, real[] LV, real[] CPU, real[] CPV) {
    return (SU(t,x, beta1, LU, CPU)^(-theta)+SU(t,x, beta2, LV, CPV)^(-theta)-1)^(-1/theta-1)*SU(t,x, beta2, LV, CPV)^(-theta-1)*fU(t, x, beta2, LV, CPV);
  }




  real FT(real t, vector x, vector beta1, vector beta2, real theta, real[] LU, real[] LV, real[] CPU, real[] CPV, real[] GQ_ti, real[] GQ_wi){
    int nx=num_elements(GQ_ti);
    real FT3=0;
    for (i in 1:nx){
      FT3=FT3+(t/2)*GQ_wi[i]*fT1(((t/2)*GQ_ti[i]+(t/2)),  x,  beta1,  beta2,  theta, LU, LV, CPU, CPV)+(t/2)*GQ_wi[i]*fT2(((t/2)*GQ_ti[i]+(t/2)),  x,  beta1,  beta2,  theta, LU, LV, CPU, CPV);
    }
    if(FT3<=0.001){FT3=0.001;}
    if(FT3>=0.999){FT3=0.999;}
    return FT3;
  }


///////////////////////////////


// individual likelihood
  real Li_T(real t, vector x, real D, vector beta1, vector beta2, real theta, real[] LU, real[] LV, real[] CPU, real[] CPV) {
    if (D==1){
      return fT1(t, x, beta1, beta2, theta, LU, LV, CPU, CPV);
    }
    else if (D==2){
      return fT2(t, x, beta1, beta2, theta, LU, LV, CPU, CPV);
    }else {
      return SUV(t, t, x, beta1, beta2, theta, LU, LV, CPU, CPV);
    }
  }

// individual log likelihood
  real li_T(real t, vector x, real D, vector beta1, vector beta2, real theta, real[] LU, real[] LV, real[] CPU, real[] CPV) {
    return log(Li_T(t, x, D, beta1, beta2, theta, LU, LV, CPU, CPV));
  }




// longitudinal contributions:

  real qt_loss(real u, real qt){
    if(u<0){
      return u*(qt-1);
    }else{
      return u*qt;
    }
  }

//  real fY(real y, vector w, vector z, vector betaT, real qt, real rrho, vector b){
 //   return (qt*(1-qt)/rrho)*exp(-qt_loss(((y-sum(betaT .* w)-sum(z .* b))/rrho),qt));
//  }


  real FY(real y, vector w, vector z, vector betaT, real qt, real rrho, vector b){
    real fyy;
    if(y>(sum(betaT .* w)+sum(z .* b))){
      fyy=1-(1-qt)*exp(-qt*(y-sum(betaT .* w)-sum(z .* b))/rrho);
    }else{
      fyy=qt*exp(-(qt-1)*(y-sum(betaT .* w)-sum(z .* b))/rrho);
    }
    if(fyy<=0.001){fyy=0.001;}
    if(fyy>=0.999){fyy=0.999;}
    return fyy;
  }


  //real rhoj(real t, real rho, real ty){
    //real rhoj;
    //rhoj=(rho/fabs(rho))*(fabs(rho)^(t-ty));
    //if (fabs(rhoj)>0.99){
      //return (rho/fabs(rho))*0.99;
    //}else{
      //return rhoj;
    //}
  //}

  real c_Rho(int K, vector y, matrix w, matrix z, matrix betaT, vector qt, vector rrho, matrix B, matrix Rho, real QFt){
    real QFy;
    real CRho;
    row_vector[K+1] QF;
    vector[K+1] TQF;
    matrix[K+1,K+1] I=diag_matrix(rep_vector(1.0,K+1));

    int i=1;
    QF[i]=QFt;
    while (i < K+1){
       QFy=inv_Phi(FY(y[i], to_vector(w[i,]), to_vector(z[i,]), to_vector(betaT[i,]), qt[i], rrho[i], to_vector(B[i,])));
       i=i+1;
       QF[i]=QFy;
     }

    TQF=QF';
    
    CRho=1/sqrt(determinant(Rho))*exp(-0.5*QF*(inverse_spd(Rho)-I)*TQF);
    if(CRho<=1e-100){
      return 1e-100;
    }else{
      return CRho;
    }
  }


  real li_Y(real y, vector w, vector z, vector betaT, real qt, real rrho, vector b){
    return log(qt*(1-qt)/rrho)+(-qt_loss(((y-sum(betaT .* w)-sum(z .* b))/rrho),qt));
  }




}



data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> N1;
  int<lower=0> N2;
  int<lower=0> N_L;
  int<lower=0> N_Y;
  int<lower=0> Ny;
  int<lower=0> Ni[N];
  int<lower=0> N_GQ;
  int<lower=0> p1;
  int<lower=0> p2;
  int<lower=0> p3;
  real ID[N];
  real T[N];
  real D[N];
  matrix[N,p1] X;
  real T1[N1];
  real T2[N2];
  //matrix[N, N_L] B;
  //matrix[N, N_L] IB;
  //matrix[N, N_L] M;
  //matrix[N, N_L] IM;
  //int Ni[N];
  real ID_Y[N_Y];
  real Ty[N_Y];
  real TI[N_Y];
  matrix[N_Y,p1] Xi;
  matrix[N_Y,p2] W;
  matrix[N_Y,p3] Z;
  real Y[N_Y];
  real Qt[K];
  real CPU[N_L];
  real CPV[N_L];
  real GQ_ti[N_GQ];
  real GQ_wi[N_GQ];
}


transformed data {
  matrix[p3,p3] identity;
  vector[p3] a;
  identity=diag_matrix(rep_vector(1.0,p3));
  a=rep_vector(0,p3);
}


parameters {
  vector[p1] beta1;
  vector[p1] beta2;
  real<lower=-1> theta;
  real LU[N_L];
  real LV[N_L];
  cholesky_factor_corr[K+1] L_omega;
  matrix[K, p2] betaT;
  vector<lower=0>[K] rrho;
  cov_matrix[p3] sig[K];
  matrix[K*N, p3] b;
}

transformed parameters {
  real Ktau;
  matrix[K+1,K+1] Rho;
  Ktau = theta/(theta+2);
  Rho = L_omega*L_omega';
}


model {
  // priors
  int J0=0;
  beta1 ~ normal(0,100);
  beta2 ~ normal(0,100);
  theta ~ normal(0,100);
  LU ~ gamma(0.2,0.4);
  LV ~ gamma(0.2,0.4);
  L_omega ~ lkj_corr_cholesky(2);

  for (i in 1:K){
    betaT[i] ~ normal(0,100);
    rrho[i] ~ normal(1,5);
    sig[i] ~ inv_wishart(p3+2, identity);
    for (j in 1:N){
      b[((i-1)*N+j)] ~ multi_normal(a,sig[i]);
    }
  }
  

  // log likelihood
  for (i in 1:N){
    int Ji=Ni[i];
    real id=ID[i];
    real t=T[i];
    real d=D[i];
    vector[p1] x=to_vector(X[i,]);
    
    real QFt=inv_Phi(FT(t, x, beta1, beta2, theta, LU, LV, CPU, CPV, GQ_ti, GQ_wi));


    target += li_T(t, x, d, to_vector(beta1), to_vector(beta2), theta, LU, LV, CPU, CPV);
    

      for (j in (J0+1):(J0+Ji)){
        real liY=0;
        vector[K] y;
        matrix[K,p2] w;
        matrix[K,p3] z;
        matrix[K,p3] bi;
        
        for (k in 1:K){
          y[k]=Y[((k-1)*Ny+j)];
          w[k,] = W[((k-1)*Ny+j), ];
          z[k,] = Z[((k-1)*Ny+j), ];
          bi[k,] = b[((k-1)*N+i), ];
          liY += li_Y(Y[((k-1)*Ny+j)],to_vector(W[((k-1)*Ny+j), ]), to_vector(Z[((k-1)*Ny+j), ]), to_vector(betaT[k,]), Qt[k], rrho[k], to_vector(b[((k-1)*N+i), ]));
        }
        
        target += log(c_Rho(K, y, w, z, betaT, to_vector(Qt), rrho, bi, Rho, QFt))+liY;
      
      }
    J0=J0+Ji;
    
 }

}


generated quantities{
   vector[K+2] DEV;
   real DEV_J;
   int J00=0;
   DEV = rep_vector(0,K+2);

   for (i in 1:N){
    real id=ID[i];
    real t=T[i];
    real d=D[i];
    vector[p1] x=to_vector(X[i,]);
    int Ji=Ni[i];
    real QFt=inv_Phi(FT(t, x, beta1, beta2, theta, LU, LV, CPU, CPV, GQ_ti, GQ_wi));
    DEV[1] += -2*li_T(t, x, d, to_vector(beta1), to_vector(beta2), theta, LU, LV, CPU, CPV);


    for (j in (J00+1):(J00+Ji)){
      vector[K] y;
      matrix[K,p2] w;
      matrix[K,p3] z;
      matrix[K,p3] bi;
      
      for (k in 1:K){
          y[k]=Y[((k-1)*Ny+j)];
          w[k,] = W[((k-1)*Ny+j), ];
          z[k,] = Z[((k-1)*Ny+j), ];
          bi[k,] = b[((k-1)*N+i), ];
          DEV[2+k] += -2*li_Y(Y[((k-1)*Ny+j)],to_vector(W[((k-1)*Ny+j), ]), to_vector(Z[((k-1)*Ny+j), ]), to_vector(betaT[k,]), Qt[k], rrho[k], to_vector(b[((k-1)*N+i), ]));
      }
      DEV[2] += -2*log(c_Rho(K, y, w, z, betaT, to_vector(Qt), rrho, bi, Rho, QFt));
    }
    J00=J00+Ji;


   }
   
   DEV_J = sum(DEV);

}

