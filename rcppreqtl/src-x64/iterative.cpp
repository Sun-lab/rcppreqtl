#include <iostream>
#include <cmath>
#include <Rmath.h>
#include "meta.h"
#include "problem.h"
#include "solver/bfgssolver.h"
#include "solver/lbfgsbsolver.h"
//#include "solver/neldermeadsolver.h"
#include "brent.h"

#include "iterative.h"

unsigned long seed[3] = {4247, 3582, 3172};

// to use this library just use the namespace "cppoptlib"
namespace cppoptlib {

//negative loglik
//define negative beta-binomial (ASE) assuming we know AS counts, mean proportion for each sample
//asn - paternal + maternal allele specific counts
//asnp - paternal allele specific counts
double ASE_lik(Vector<double> asn, Vector<double> asnp, Vector<double> prob, double lphi, int n){
  int k;
  double out = 0.0, out2  = 0.0;
  double phi = exp(lphi);
  for(int i=0; i<n; i++){
    out2 = 0.0;
    if(asnp[i] > 0){
      for(k=0; k<asnp[i]; k++) out2 -= log(prob[i] + k*phi);
    }
    
    if(asnp[i] < asn[i]){
      for(k=0; k<(asn[i]-asnp[i]); k++) out2 -= log(1-prob[i] + k*phi);
    }
    
    for(k=0; k<asn[i]; k++){
      out2 += log(1+k*phi);
    }
    out += out2;
  }        
  return out;
}


double ASE_lik2(Vector<double> asn, Vector<double> asnp, Vector<double> prob, double lphi, int n){
  int k;
  double alpha, beta;
  double out = 0.0, out2  = 0.0;
  double phi = exp(lphi);
  for(int i=0; i<n; i++){
    out2 = 0.0;

    alpha = prob[i]/phi;
    beta  = (1.0 - prob[i])/phi;    
    out2 = R::lbeta(asnp[i] + alpha, asn[i] - asnp[i] + beta) - R::lbeta(alpha, beta);
    out += out2;
  }        
  return -out;
}

//negative loglik
//define negative binomial TREC assuming we know total read counts, log mean value for each sample
double TREC_lik(Vector<double> y, Vector<double> yf, Vector<double> lmu, double liphi, int n){
//liphi = log(phi^-1) or -log(phi)
  double out = 0; 
  double iphi = exp(liphi);
  for(int i=0; i<n; i++){
    out -= lgamma(y[i]+iphi)-yf[i] - (iphi+y[i])*log(iphi+exp(lmu[i])) + y[i]*lmu[i];
  }    
  out -= n*(iphi*liphi-lgamma(iphi));
  if(!std::isfinite(out))out = MLIK;
  return out;
}

Vector<double> wt(Vector<double> lmu, int n, double iphi){
  Vector<double> out(n);
  double ldet = 0;
  for(int i=0;i<n;i++){
    out[i] = lmu[i] - log1p(iphi*exp(lmu[i]));
    ldet += out[i];
  }
//  std::cout << ldet << " ";
  ldet /= n;
  for(int i=0;i<n;i++){
    out[i] = exp(out[i]-ldet);
  }
  return out;
}
Vector<double> vexp(Vector<double> v, int n){
  Vector<double> out(n);
  for(int i=0;i<n;i++){
    out[i] = exp(v[i]);
  }
  return out;
}
Vector<double> vlog1p(Vector<double> v, int n){
  Vector<double> out(n);
  for(int i=0;i<n;i++){
    out[i] = log1p(v[i]);
  }
  return out;
}
}

extern "C" {

//this function employs 1-dimensional Brent optimizer for overdispersion
//and iteratively reweighting for mean structure
// general procedure: 
//(by default first estimate mean structure using linear model for log counts, can turn it off)
// tune betas at a given initial value of overdispersion
// start iteratively reweighting (get new optimal overdispersion keeping linear mean structure constant)-(get new linear mean structure keeping overdispersion constant)

double max_trec1irw0(int* RnInd, int* Rngenes, int* Rthp, double* Rtrc, int* RnPar, double* Rinits, double* RXmatr, int* Restinit) {
    int maxiter = 100,maxiter2=50,iter,iter2,tryiter;
    double eps = 1e-10;
//    double gold = (1+sqrt(5))/2; //maybe instead of dividing will think about doing goldern rule in future?
    const int np = RnPar[0];//number of param
    const int n = RnInd[0];//number of individuals
    const int ngenes = Rngenes[0];//number of genes, 1 gene for now    
    const int estinit = Restinit[0];
    const int nnlin = 1;//number of nonlinear par
    const int nb = np-nnlin;

    //vector with parameters to be optimized using numerical optimizer (the rest are betas - to be IRW in the code)
    cppoptlib::Vector<double> yy(1);

    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> trc(Rtrc,n,1);    
    Eigen::MatrixXd k(n,1);

    //expected counts for each individual on log scale - tmn and without log scale etmn
    cppoptlib::Vector<double> tmn(n);//log mean of total counts
    cppoptlib::Vector<double> etmn(n); //meat of total read counts
    
    Eigen::Map<Eigen::Matrix<int,Eigen::Dynamic,1>> thp(Rthp,n,1);    
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>> X(RXmatr,n,nb);
    Eigen::DiagonalMatrix<double,Eigen::Dynamic> Wt(n);

    cppoptlib::Vector<double> inits(np);for(int i=0;i<np;i++)inits[i] = Rinits[i];   //all inits
    cppoptlib::Vector<double> betas(nb);for(int i=nnlin;i<np;i++)betas[i-nnlin] = Rinits[i]; //linear betas
    cppoptlib::Vector<double> betas0(nb);betas0=betas; //prev iter betas
    cppoptlib::Vector<double> trybetas(nb);trybetas=betas; //for line search
    cppoptlib::Vector<double> delta(nb);delta.setZero(); //prev iter betas
    cppoptlib::Vector<double> nlin(1);for(int i=0;i<1;i++)nlin[i] = Rinits[i];       //nonlinear part
    cppoptlib::Vector<double> nlinf(2);nlinf.setZero();nlinf[0]=exp(nlin[0]);       //nonlinear part
    double trynlin = nlin[0]; //for line search
    cppoptlib::Matrix<double> pre;
    cppoptlib::TRECirw0 <double> f_irw;

    double a = LPHI; //lower boundary for log(phi^-1)
    double b = UPHI; //upper boundary for log(phi^-1)
    double t = 1e-15;
    double res, x = 0;
    double ll1=1e6,ll0=1e6,ll2=1e6,tryll=1e6,ll_old=1e6;

    yy[0] = nlin[0];
    int flag = 1, flag0 = 1, flag1 = 1;
//    for(int j=0;j<N;j++){ //this will be when we will add multiple genes per run
      yy = nlin;
      //generally I use the attributes to pass additional arguments - they are to be defined in problem.h
      f_irw.setNind(n);
      f_irw.setNbet(nb);    
      f_irw.setTRC(trc);
      f_irw.setTRCF(trc);    
      f_irw.setTHP(thp);
      f_irw.setTMN(trc);
      f_irw.setETMN(trc);//do we use it?
      f_irw.setX(X);
      f_irw.setBetas(betas);
      f_irw.setPartmu();

      //estimate betas with simple lm
      if(estinit){
        k = cppoptlib::vlog1p(trc,n);
        betas = (X.transpose()*X).inverse()*(X.transpose()*k);        
        f_irw.setBetas(betas);
        f_irw.setPartmu();        
      }
      //do IRW of these initials betas assuming initial overdispersion
      iter = 0;
      ll1 = f_irw(yy);
      do{
        ll0 = ll1;
        tmn = f_irw.updTMNb0(nlinf);
        //weights are mean/(1.0+mean*phi^-1)
        Wt.diagonal() = cppoptlib::wt(tmn,n,exp(yy[0]));
        etmn = cppoptlib::vexp(tmn,n);

        k = (trc-etmn).cwiseQuotient(etmn);
        pre = (X.transpose()*Wt*X).inverse()*(X.transpose()*Wt);
        delta = pre*k;        
        betas0 = betas;      
        betas = betas + delta;
        f_irw.setBetas(betas);
        f_irw.setPartmu();
        ll1 = f_irw(yy);

        //get reasonable initial values of betas assuming not large overdispersion 
        tryiter = 0;
        flag = 1;
        do{
          trybetas = betas0;
          delta /= 2;
          trybetas = trybetas + delta;
          f_irw.setBetas(trybetas);
          f_irw.setPartmu();
          tryll = f_irw(yy);
          if(tryll>ll1){            
            f_irw.setBetas(betas);
            f_irw.setPartmu();
            flag = 0;
          }else{
            betas = trybetas;
            ll1 = tryll;
          }
          tryiter++;
        }while((tryiter<maxiter2) && (flag==1));
        iter++;
      }while( (iter< maxiter2) && (ll0-ll1)>eps);
      
      iter = 0;
      yy = nlin;

      //now update between estimating non-linear and IRW of linear
      do{
        ll0 = ll1;
        x = yy[0];
        res = local_min (a, b, t, f_irw, x );
        yy[0] = x;

        if((f_irw(yy)-ll0)>eps){
          //if after optimization we get worse result get back to old overdispersion
          yy = nlin;
        }else{
          //otherwise update overdispersion and update mean structure
          nlin = yy;
          nlinf[0] = yy[0];
        }
        
        tmn = f_irw.updTMNb0(nlinf);
        
        //for current overdispersion
        //weights are mean/(1.0+mean*phi^-1)
        Wt.diagonal() = cppoptlib::wt(tmn,n,exp(yy[0]));
        etmn = cppoptlib::vexp(tmn,n);

        k = (trc-etmn).cwiseQuotient(etmn);        
        delta = (X.transpose()*Wt*X).inverse()*(X.transpose()*Wt)*k;        
        betas0 = betas;      
        betas = betas + delta;
        f_irw.setBetas(betas);
        f_irw.setPartmu();
        ll1 = f_irw(yy);

        //get betas for current overdispersion 
        tryiter = 0;
        flag = 1;
        do{
          trybetas = betas0;
          delta /= 2;
          trybetas = trybetas + delta;
          f_irw.setBetas(trybetas);
          f_irw.setPartmu();
          tryll = f_irw(yy);
          if(ll1 - tryll < eps){            
            f_irw.setBetas(betas);
            f_irw.setPartmu();
            flag = 0;
          }else{
            betas = trybetas;
            ll1 = tryll;
          }
          tryiter++;
        }while((delta.norm()>1)&&(tryiter<maxiter2) && (flag==1));

        iter++;
        flag = (int)(((ll0-ll1)>eps)&&(iter<maxiter));        
      }while(flag==1);      
      //if the last step failed to optimize return to previous values
      if(ll0<ll1){
        betas=betas0;
        f_irw.setBetas(betas);
        f_irw.setPartmu();      
      }        
      ll1 = f_irw(yy);
//    }          
    for(int i=nnlin;i<np;i++)Rinits[i]=betas[i-nnlin];
    Rinits[0] = yy[0];
    Rinits[np] = ll1;

    return ll1;
}


//this function employs step-wise 1-dimensional BRENT for poo effect and overdispersion
//and iteratively reweighting for mean structure
// general procedure: 
//(by default first estimate mean structure using linear model for log counts, can turn it off)
// tune betas at a given initial value of overdispersion
// start iteratively reweighting (get new optimal overdispersion keeping linear mean structure constant)-(get new linear mean structure keeping overdispersion constant)
double max_trec1irw_b0(int* RnInd, int* Rngenes, int* Rthp, double* Rtrc, int* RnPar, double* Rinits, double* RXmatr, int* Restinit) {
    int maxiter = 100,maxiter2=50,iter,iter2,tryiter;
    int eps = 1e-8;

    double tol = 1e-15;
    double lbnd_phi = LPHI; //lower boundary for log(phi^-1)
    double ubnd_phi = UPHI; //upper boundary for log(phi^-1)
    double lbnd_b = LBND2; //lower boundary for genetic effect
    double ubnd_b = UBND2; //upper boundary for genetic effect
    double res1,res2;
    
    const int n = RnInd[0];//number of individuals
    const int np = RnPar[0];//number of param
    const int ngenes = Rngenes[0];//number of genes, 1 gene for now    
    const int estinit = Restinit[0];
    const int nnlin = 2;//number of nonlinear par
    const int nb = np-nnlin;
    cppoptlib::Vector<double> opt_liphi(1);
    cppoptlib::Vector<double> opt_b0(1);    
    cppoptlib::Vector<double> nlin_liphi(1);
    cppoptlib::Vector<double> nlin_b0(1);

    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> trc(Rtrc,n,1);    
    Eigen::MatrixXd k(n,1);
    //expected counts for each individual on log scale - tmn and without log scale etmn
    cppoptlib::Vector<double> tmn(n);  
    cppoptlib::Vector<double> etmn(n);  
        
    Eigen::Map<Eigen::Matrix<int,Eigen::Dynamic,1>> thp(Rthp,n,1);    
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>> X(RXmatr,n,nb);
    Eigen::DiagonalMatrix<double,Eigen::Dynamic> Wt(n);

    cppoptlib::Vector<double> inits(np);for(int i=0;i<np;i++)inits[i] = Rinits[i];
    cppoptlib::Vector<double> betas(nb);for(int i=nnlin;i<np;i++)betas[i-nnlin] = Rinits[i];
    nlin_liphi[0] = inits[0];
    nlin_b0[0] = inits[1];
    
    cppoptlib::Vector<double> betas0(nb);betas0 = betas;
    cppoptlib::Vector<double> trybetas(nb);trybetas = betas;    
    cppoptlib::Vector<double> delta(nb);delta.setZero(); //prev iter betas
    

    cppoptlib::Matrix<double> pre;    

    double ll1=1e6,ll0=1e6,ll2=1e6,tryll=1e6,ll_old=1e6;
    int flag = 1, flag0 = 1, flag1 = 1;

    //do optim using minimizer
    cppoptlib::TRECirw_b0 <double> f_irw;
      opt_liphi = nlin_liphi;
      opt_b0 = nlin_b0;
      
      f_irw.setNind(n);
      f_irw.setNbet(nb);    
      f_irw.setTRC(trc);
      f_irw.setTRCF(trc);    
      f_irw.setTHP(thp);
      f_irw.setTMN(trc);
      f_irw.setETMN(trc);
      f_irw.setBetas(betas);
      f_irw.setX(X);
      f_irw.setPartmu();
      f_irw.setLiphi(nlin_liphi[0]);

      if(estinit){
        k = cppoptlib::vlog1p(trc,n);
        betas = (X.transpose()*X).inverse()*(X.transpose()*k);        
        f_irw.setBetas(betas);
        f_irw.setPartmu();        
      }
      iter = 0;
      f_irw.setDir(2);      //dir = 1 for liphi
      ll1 = f_irw(nlin_b0);
      do{
        ll0 = ll1;
        f_irw.setDir(2);
        tmn = f_irw.updTMNb0(nlin_b0);
        //weights are mean/(1.0+mean*phi^-1)
        Wt.diagonal() = cppoptlib::wt(tmn,n,exp(nlin_liphi[0]));//weights are 1/(1+a mu)
        etmn = cppoptlib::vexp(tmn,n);

        k = (trc-etmn).cwiseQuotient(etmn);
        delta = (X.transpose()*Wt*X).inverse()*(X.transpose()*Wt)*k;        
        betas0 = betas;      
        betas = betas + delta;
        f_irw.setBetas(betas);
        f_irw.setPartmu();
        f_irw.setDir(2);
        ll1 = f_irw(nlin_b0);

        //get reasonable initial values of betas assuming initial large overdispersion 
        tryiter = 0;
        flag = 1;
        do{
          trybetas = betas0;
          delta /= 2;
          trybetas = trybetas + delta;
          f_irw.setBetas(trybetas);
          f_irw.setPartmu();
          tryll = f_irw(nlin_b0);
          if(tryll>ll1){            
            f_irw.setBetas(betas);
            f_irw.setPartmu();
            flag = 0;
          }else{
            betas = trybetas;
            ll1 = tryll;
          }
          tryiter++;
        }while((tryiter<maxiter2) && (flag==1));
        iter++;
      }while( (iter< maxiter2) && (ll0-ll1)>eps);

      iter = 0;
      opt_b0 = nlin_b0;
      opt_liphi = nlin_liphi;      
      //now update non-linear - recalculate betas iteratively
      do{
        ll0 = ll1;
        //dir = 1 for iphi and dir = 2 for b0
        f_irw.setDir(2);        
        res1 = local_min (nlin_b0[0]+lbnd_b, nlin_b0[0]+ubnd_b, tol, f_irw, opt_b0[0] );
        if(res1 < ll0){
          f_irw.setB0(opt_b0[0]);
        }else{
          opt_b0=nlin_b0;
        }
        f_irw.setDir(1);        
        res2 = local_min (lbnd_phi, ubnd_phi, tol, f_irw, opt_liphi[0] );
        if(res2 > res1)opt_liphi=nlin_liphi;

        if((f_irw(opt_liphi)-ll0)>eps){
          //if after optimization we get worse result get back to old overdispersion
          opt_liphi = nlin_liphi;
          opt_b0 = nlin_b0;
        }else{
          //otherwise update overdispersion and update mean structure
          nlin_liphi = opt_liphi;
          nlin_b0 = opt_b0;
        }
        f_irw.setLiphi(nlin_liphi[0]);
        f_irw.setB0(nlin_b0[0]);        

        tmn = f_irw.updTMNb0(nlin_b0);
        //for current overdispersion
        //weights are mean/(1.0+mean*phi^-1)
        Wt.diagonal() = cppoptlib::wt(tmn,n,exp(nlin_liphi[0]));
        etmn = cppoptlib::vexp(tmn,n);
        k = (trc-etmn).cwiseQuotient(etmn);        
        delta = (X.transpose()*Wt*X).inverse()*(X.transpose()*Wt)*k;        
        betas0 = betas;      
        betas = betas + delta;
        f_irw.setBetas(betas);
        f_irw.setPartmu();
        f_irw.setLiphi(nlin_liphi[0]);
        f_irw.setB0(nlin_b0[0]);
        f_irw.setDir(2);
        ll1 = f_irw(nlin_b0);

        //get betas for current overdispersion 
        tryiter = 0;
        flag = 1;
        do{
          trybetas = betas0;
          delta /= 2;
          trybetas = trybetas + delta;
          f_irw.setBetas(trybetas);
          f_irw.setPartmu();
          f_irw.setDir(2);
          tryll = f_irw(nlin_b0);
          if(ll1 - tryll < eps){            
            f_irw.setBetas(betas);
            f_irw.setPartmu();
            flag = 0;
          }else{
            betas = trybetas;
            ll1 = tryll;
          }
          tryiter++;
        }while((delta.norm()>1)&&(tryiter<maxiter2) && (flag==1));
        iter++;
        flag = (int)(((ll0-ll1)>eps)&&(iter<maxiter));        
      }while(flag==1);
      if(ll0<ll1){
        betas=betas0;
        f_irw.setBetas(betas);
        f_irw.setPartmu();      
      }        
      f_irw.setDir(2);
      ll1 = f_irw(nlin_b0);
//    }    
    Rinits[0]=nlin_liphi[0];
    Rinits[1]=nlin_b0[0];    
    for(int i=nnlin;i<np;i++)Rinits[i]=betas[i-nnlin];
    Rinits[np] = ll1;
    
    return ll1;
}

//this function employs 2-dimensional L-BFGS-B optimizer for overdispersion
//and iteratively reweighting for mean structure
// general procedure: 
//(by default first estimate mean structure using linear model for log counts, can turn it off)
// tune betas at a given initial value of overdispersion
// start iteratively reweighting (get new optimal overdispersion keeping linear mean structure constant)-(get new linear mean structure keeping overdispersion constant)
double max_trec1irw_b1(int* RnInd, int* Rngenes, int* Rthp, double* Rtrc, int* RnPar, double* Rinits, double* RXmatr, int* Restinit, double* Rb0){
    int maxiter = 100,maxiter2=50,iter,iter2,tryiter;
    int eps = 1e-8;

    double tol = 1e-20;
    double lbnd_phi = LPHI; //lower boundary for log(phi^-1)
    double ubnd_phi = UPHI; //upper boundary for log(phi^-1)
    double lbnd_b = LBND2; //lower boundary for genetic effect
    double ubnd_b = UBND2; //upper boundary for genetic effect
    double res1,res2;
    
    const int n = RnInd[0];//number of individuals
    const int np = RnPar[0];//number of param
    const int ngenes = Rngenes[0];//number of genes, 1 gene for now
    const int estinit = Restinit[0];
    const int nnlin = 2;//number of nonlinear par
    const int nb = np-nnlin;
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>> b0(Rb0,ngenes,1);
    cppoptlib::Vector<double> opt_liphi(1);
    cppoptlib::Vector<double> opt_b1(1);    
    cppoptlib::Vector<double> nlin_liphi(1);
    cppoptlib::Vector<double> nlin_b1(1);

    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> trc(Rtrc,n,1);    
    Eigen::MatrixXd k(n,1);
    //expected counts for each individual on log scale - tmn and without log scale etmn
    cppoptlib::Vector<double> tmn(n);  
    cppoptlib::Vector<double> etmn(n);  
        
    Eigen::Map<Eigen::Matrix<int,Eigen::Dynamic,1>> thp(Rthp,n,1);    
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>> X(RXmatr,n,nb);
    Eigen::DiagonalMatrix<double,Eigen::Dynamic> Wt(n);

    cppoptlib::Vector<double> inits(np);for(int i=0;i<np;i++)inits[i] = Rinits[i];
    cppoptlib::Vector<double> betas(nb);for(int i=nnlin;i<np;i++)betas[i-nnlin] = Rinits[i];
    nlin_liphi[0] = inits[0];
    nlin_b1[0] = inits[1];
    
    cppoptlib::Vector<double> betas0(nb);betas0 = betas;
    cppoptlib::Vector<double> trybetas(nb);trybetas = betas;    
    cppoptlib::Vector<double> delta(nb);delta.setZero(); //prev iter betas
    

    cppoptlib::Matrix<double> pre;    

    double ll1=1e6,ll0=1e6,ll2=1e6,tryll=1e6,ll_old=1e6;
    int flag = 1, flag0 = 1, flag1 = 1;

    //do optim using minimizer
    int genej = 0;
    cppoptlib::TRECirw_b1 <double> f_irw;
      opt_liphi = nlin_liphi;
      opt_b1 = nlin_b1;
      
      f_irw.setNind(n);
      f_irw.setNbet(nb);    
      f_irw.setTRC(trc);
      f_irw.setTRCF(trc);    
      f_irw.setTHP(thp);
      f_irw.setTMN(trc);
      f_irw.setETMN(trc);
      f_irw.setBetas(betas);
      f_irw.setX(X);
      f_irw.setPartmu();
      f_irw.setLiphi(nlin_liphi[0]);
      f_irw.setB1(nlin_b1[0]);
      f_irw.setB0((double)b0(genej,0));    

      if(estinit){
        k = cppoptlib::vlog1p(trc,n);
        betas = (X.transpose()*X).inverse()*(X.transpose()*k);        
        f_irw.setBetas(betas);
        f_irw.setPartmu();        
      }
      iter = 0;
      f_irw.setDir(2);      //dir = 1 for liphi
      ll1 = f_irw(nlin_b1);
      do{
        ll0 = ll1;
        f_irw.setDir(2);
        tmn = f_irw.updTMNb1(nlin_b1);
        //weights are mean/(1.0+mean*phi^-1)
        Wt.diagonal() = cppoptlib::wt(tmn,n,exp(nlin_liphi[0]));//weights are 1/(1+a mu)
        etmn = cppoptlib::vexp(tmn,n);
        k = (trc-etmn).cwiseQuotient(etmn);
        delta = (X.transpose()*Wt*X).inverse()*(X.transpose()*Wt)*k;        
        betas0 = betas;      
        betas = betas + delta;
        f_irw.setBetas(betas);
        f_irw.setPartmu();
        f_irw.setDir(2);
        ll1 = f_irw(nlin_b1);

        //get reasonable initial values of betas assuming initial large overdispersion 
        tryiter = 0;
        flag = 1;
        do{
          trybetas = betas0;
          delta /= 2;
          trybetas = trybetas + delta;
          f_irw.setBetas(trybetas);
          f_irw.setPartmu();
          tryll = f_irw(nlin_b1);
          if(tryll>ll1){            
            f_irw.setBetas(betas);
            f_irw.setPartmu();
            flag = 0;
          }else{
            betas = trybetas;
            ll1 = tryll;
          }
          tryiter++;
        }while((tryiter<maxiter2) && (flag==1));
        iter++;
      }while( (iter< maxiter2) && (ll0-ll1)>eps);
      f_irw.setBetas(betas);
      f_irw.setPartmu();
      iter = 0;
      opt_b1 = nlin_b1;
      opt_liphi = nlin_liphi;      
      //now update non-linear - recalculate betas iteratively
      do{
        ll0 = ll1;
        f_irw.setDir(2);        
        res1 = local_min (lbnd_b, ubnd_b, tol, f_irw, opt_b1[0] );
        if(res1 < ll0){
          f_irw.setB1(opt_b1[0]);
        }else{
          opt_b1=nlin_b1;
        }        
//         glomin ( double a, double b, double c, double m, double e, double t,  func_base& f, double &x );
        f_irw.setDir(1);        
        res2 = local_min (lbnd_phi, ubnd_phi, tol, f_irw, opt_liphi[0] );
        if(res2 < res1){
          f_irw.setLiphi(opt_liphi[0]);
        }else{
          opt_liphi=nlin_liphi;
        }

        if((f_irw(opt_liphi)-ll0)>eps){
          //if after optimization we get worse result get back to old overdispersion
          opt_liphi = nlin_liphi;
          opt_b1 = nlin_b1;
        }else{
          //otherwise update overdispersion and b1
          nlin_liphi = opt_liphi;
          nlin_b1 = opt_b1;
        }
        tmn = f_irw.updTMNb1(nlin_b1);
        //for current overdispersion
        //weights are mean/(1.0+mean*phi^-1)
        Wt.diagonal() = cppoptlib::wt(tmn,n,exp(nlin_liphi[0]));
        etmn = cppoptlib::vexp(tmn,n);
        k = (trc-etmn).cwiseQuotient(etmn);        
        delta = (X.transpose()*Wt*X).inverse()*(X.transpose()*Wt)*k;        
        betas0 = betas;      
        betas = betas + delta;
        f_irw.setBetas(betas);
        f_irw.setPartmu();
        f_irw.setDir(2);
        ll1 = f_irw(nlin_b1);

        //get betas for current overdispersion 
        tryiter = 0;
        flag = 1;
        do{
          trybetas = betas0;
          delta /= 2;
          trybetas = trybetas + delta;
          f_irw.setBetas(trybetas);
          f_irw.setPartmu();
          f_irw.setDir(2);
          tryll = f_irw(nlin_b1);
          if(ll1 - tryll < eps){            
            f_irw.setBetas(betas);
            f_irw.setPartmu();
            flag = 0;
          }else{
            betas = trybetas;
            ll1 = tryll;
          }
          tryiter++;
        }while((delta.norm()>1)&&(tryiter<maxiter2) && (flag==1));
        iter++;
        flag = (int)(((ll0-ll1)>eps)&&(iter<maxiter));        
      }while(flag==1);
      if(ll0<ll1){
        betas=betas0;
        f_irw.setBetas(betas);
        f_irw.setPartmu();      
      }        
      f_irw.setDir(2);
      ll1 = f_irw(nlin_b1);
//    }    
    Rinits[0]=nlin_liphi[0];
    Rinits[1]=nlin_b1[0];    
    for(int i=nnlin;i<np;i++)Rinits[i]=betas[i-nnlin];
    Rinits[np] = ll1;
    
    return ll1;
}


//full TRECASE, fitting b0 - additive
//as.integer(nInd), as.integer(nGenes), as.integer(thp),as.double(trc),as.double(asn),as.double(asnp),
//  as.integer(nPar),optval=as.double(inits),as.double(Xmatr), as.integer(est.inits)

double max_trecase1irw_b0(int* RnInd, int* RnAse, int* Rngenes, int* Rthp, double* Rtrc, double* Rasn, double* Rasnp, 
                          int* RnPar, double* Rinits, double* RXmatr, int* Restinit) {
    int maxiter = 100,maxiter2=50,iter,iter2,tryiter;
    int eps = 1e-8;
    double tol = 1e-20;
    double lbnd_phi = LPHI; //lower boundary for log(phi^-1)
    double ubnd_phi = UPHI; //upper boundary for log(phi^-1)
    double lbnd_b = LBND2; //lower boundary for genetic effect
    double ubnd_b = UBND2; //upper boundary for genetic effect
    double res1,res2,res3;
    
    const int n = RnInd[0];//number of individuals
    const int nase = RnAse[0];
    const int np = RnPar[0];//number of param
    const int ngenes = Rngenes[0];//number of genes, 1 gene for now
        
    const int estinit = Restinit[0];
    const int nnlin = 3;//number of nonlinear par
    const int nb = np-nnlin;
    cppoptlib::Vector<double> opt_liphi(1);
    cppoptlib::Vector<double> opt_lphi(1);
    cppoptlib::Vector<double> opt_b0(1);    
    cppoptlib::Vector<double> nlin_liphi(1);
    cppoptlib::Vector<double> nlin_lphi(1);
    cppoptlib::Vector<double> nlin_b0(1);
    
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> trc(Rtrc,n,1);    
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> asn(Rasn,nase,1);    
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> asnp(Rasnp,nase,1);    

    Eigen::MatrixXd k(n,1);
    //expected counts for each individual on log scale - tmn and without log scale etmn
    cppoptlib::Vector<double> tmn(n);  
    cppoptlib::Vector<double> etmn(n);  
        
    Eigen::Map<Eigen::Matrix<int,Eigen::Dynamic,1>> thp(Rthp,nase,1);    
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>> X(RXmatr,n,nb);
    Eigen::DiagonalMatrix<double,Eigen::Dynamic> Wt(n);

    cppoptlib::Vector<double> inits(np);for(int i=0;i<np;i++)inits[i] = Rinits[i];
    cppoptlib::Vector<double> betas(nb);for(int i=nnlin;i<np;i++)betas[i-nnlin] = Rinits[i];
    nlin_liphi[0] = inits[0];
    nlin_lphi[0] = inits[1];
    nlin_b0[0] = inits[2];
    
    cppoptlib::Vector<double> betas0(nb);betas0 = betas;
    cppoptlib::Vector<double> trybetas(nb);trybetas = betas;    
    cppoptlib::Vector<double> delta(nb);delta.setZero(); //prev iter betas
    

    cppoptlib::Matrix<double> pre;    

    double ll1=1e6,ll0=1e6,ll2=1e6,tryll=1e6,ll_old=1e6;
    int flag = 1, flag0 = 1, flag1 = 1;

    //do optim using minimizer
    cppoptlib::TRECASEirw_b0 <double> f_irw;
      opt_liphi = nlin_liphi;
      opt_lphi = nlin_lphi;
      opt_b0 = nlin_b0;
      
      f_irw.setNind(n);
      f_irw.setNase(nase);
      f_irw.setNbet(nb);    
      f_irw.setTRC(trc);
      f_irw.setAS(asn);
      f_irw.setASp(asnp);      
      f_irw.setTRCF(trc);    
      f_irw.setTHP(thp);
      f_irw.setTMN(trc);
      f_irw.setETMN(trc);
      f_irw.setBetas(betas);
      f_irw.setX(X);
      f_irw.setPartmu();
      f_irw.setProb(asn);
      f_irw.setLiphi(nlin_liphi[0]);
      f_irw.setLphi(nlin_lphi[0]);

      if(estinit){
        k = cppoptlib::vlog1p(trc,n);
        betas = (X.transpose()*X).inverse()*(X.transpose()*k);        
        f_irw.setBetas(betas);
        f_irw.setPartmu();        
      }
      iter = 0;
      f_irw.setDir(3);      //dir = 1 for liphi, dir = 3 for lphi
      ll1 = f_irw(nlin_b0);
      do{
        ll0 = ll1;
        f_irw.setDir(3);
        tmn = f_irw.updTMNb0(nlin_b0);
        //weights are mean/(1.0+mean*phi^-1)
        Wt.diagonal() = cppoptlib::wt(tmn,n,exp(nlin_liphi[0]));//weights are 1/(1+a mu)
        etmn = cppoptlib::vexp(tmn,n);

        k = (trc-etmn).cwiseQuotient(etmn);
        delta = (X.transpose()*Wt*X).inverse()*(X.transpose()*Wt)*k;        
        betas0 = betas;      
        betas = betas + delta;
        f_irw.setBetas(betas);
        f_irw.setPartmu();
        f_irw.setDir(3);
        ll1 = f_irw(nlin_b0);

        //get reasonable initial values of betas assuming initial large overdispersion 
        tryiter = 0;
        flag = 1;
        
        do{
          trybetas = betas0;
          delta /= 2;
          trybetas = trybetas + delta;
          f_irw.setBetas(trybetas);
          f_irw.setPartmu();
          tryll = f_irw(nlin_b0);
          if(tryll>ll1){            
            f_irw.setBetas(betas);
            f_irw.setPartmu();
            flag = 0;
          }else{
            betas = trybetas;
            ll1 = tryll;
          }
          tryiter++;
        }while((tryiter<maxiter2) && (flag==1));
        iter++;
      }while( (iter< maxiter2) && (ll0-ll1)>eps);

      iter = 0;
      opt_b0 = nlin_b0;
      opt_liphi = nlin_liphi;      
      //now update non-linear - recalculate betas iteratively
      do{
        ll0 = ll1;
        //dir = 1 for iphi and dir = 2 for b0
        f_irw.setDir(3);        
        res1 = local_min (nlin_b0[0]+lbnd_b, nlin_b0[0]+ubnd_b, tol, f_irw, opt_b0[0] );
        if(res1 < ll0){
          f_irw.setB0(opt_b0[0]);
        }else{
          opt_b0=nlin_b0;
        }


        f_irw.setDir(1);        
        lbnd_phi = LPHI + opt_liphi[0];
        if(lbnd_phi<LPHIU)lbnd_phi=LPHIU;
        ubnd_phi = UPHI + opt_liphi[0];
        if(ubnd_phi>UPHIU)ubnd_phi=UPHIU;
        res2 = local_min (lbnd_phi, ubnd_phi, tol, f_irw, opt_liphi[0] );
        if(res2 > res1)opt_liphi=nlin_liphi;

        f_irw.setDir(2);        
        lbnd_phi = LPHI + opt_lphi[0];
        if(lbnd_phi<LPHIU)lbnd_phi=LPHIU;
        ubnd_phi = UPHI + opt_lphi[0];
        if(ubnd_phi>UPHIU)ubnd_phi=UPHIU;
        res3 = local_min (lbnd_phi, ubnd_phi, tol, f_irw, opt_lphi[0] );
        if(res3 > res2)opt_lphi=nlin_lphi;

        if((f_irw(opt_lphi)-ll0)>eps){
          //if after optimization we get worse result get back to old overdispersion
          opt_liphi = nlin_liphi;
          opt_lphi = nlin_lphi;
          opt_b0 = nlin_b0;
        }else{
          //otherwise update overdispersion and update mean structure
          nlin_liphi = opt_liphi;
          nlin_lphi = opt_lphi;
          nlin_b0 = opt_b0;
        }
        f_irw.setLiphi(nlin_liphi[0]);
        f_irw.setLphi(nlin_lphi[0]);
        f_irw.setB0(nlin_b0[0]);        

        tmn = f_irw.updTMNb0(nlin_b0);
        //for current overdispersion
        //weights are mean/(1.0+mean*phi^-1)
        Wt.diagonal() = cppoptlib::wt(tmn,n,exp(nlin_liphi[0]));
        etmn = cppoptlib::vexp(tmn,n);
        k = (trc-etmn).cwiseQuotient(etmn);        
        delta = (X.transpose()*Wt*X).inverse()*(X.transpose()*Wt)*k;        
        betas0 = betas;      
        betas = betas + delta;
        f_irw.setBetas(betas);
        f_irw.setPartmu();
        f_irw.setLiphi(nlin_liphi[0]);
        f_irw.setB0(nlin_b0[0]);
        f_irw.setDir(3);
        ll1 = f_irw(nlin_b0);

        //get betas for current overdispersion 
        tryiter = 0;
        flag = 1;
        do{
          trybetas = betas0;
          delta /= 2;
          trybetas = trybetas + delta;
          f_irw.setBetas(trybetas);
          f_irw.setPartmu();
          f_irw.setDir(3);
          tryll = f_irw(nlin_b0);
          if(ll1 - tryll < eps){            
            f_irw.setBetas(betas);
            f_irw.setPartmu();
            flag = 0;
          }else{
            betas = trybetas;
            ll1 = tryll;
          }
          tryiter++;
        }while((delta.norm()>1)&&(tryiter<maxiter2) && (flag==1));
        iter++;
        flag = (int)(((ll0-ll1)>eps)&&(iter<maxiter));        
      }while(flag==1);
      if(ll0<ll1){
        betas=betas0;
        f_irw.setBetas(betas);
        f_irw.setPartmu();      
      }        
      f_irw.setDir(3);
      ll1 = f_irw(nlin_b0);
//    }    
    Rinits[0]=nlin_liphi[0];
    Rinits[1]=nlin_lphi[0];    
    Rinits[2]=nlin_b0[0];    
    for(int i=nnlin;i<np;i++)Rinits[i]=betas[i-nnlin];
    Rinits[np] = ll1;

    return ll1;
}


//full TRECASE, fitting b1 - poo
//as.integer(nInd), as.integer(thp),as.double(trc),as.double(asn),as.double(asnp),
//  as.integer(nPar),optval=as.double(inits),as.double(Xmatr), as.integer(est.inits), as.integer(b0)
double max_trecase1irw_b1(int* RnInd, int* RnAse, int* Rngenes, int* Rthp, double* Rtrc, double* Rasn, double* Rasnp, 
                          int* RnPar, double* Rinits, double* RXmatr, int* Restinit,  double* Rb0) {
    int maxiter = 100,maxiter2=50,iter,iter2,tryiter;
    int eps = 1e-8;
    double tol = 1e-20;
    double lbnd_phi = LPHI; //lower boundary for log(phi^-1)
    double ubnd_phi = UPHI; //upper boundary for log(phi^-1)
    double lbnd_b = LBND2; //lower boundary for genetic effect
    double ubnd_b = UBND2; //upper boundary for genetic effect
    double res1,res2,res3;
    
    const int n = RnInd[0];//number of individuals
    const int nase = RnAse[0];
    const int np = RnPar[0];//number of param
    const int ngenes = Rngenes[0];//number of genes, 1 gene for now    
    const int estinit = Restinit[0];
    const int nnlin = 3;//number of nonlinear par
    const int nb = np-nnlin;
    cppoptlib::Vector<double> opt_liphi(1);
    cppoptlib::Vector<double> opt_lphi(1);
    cppoptlib::Vector<double> opt_b1(1);    
    cppoptlib::Vector<double> nlin_liphi(1);
    cppoptlib::Vector<double> nlin_lphi(1);
    cppoptlib::Vector<double> nlin_b1(1);

    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>> b0(Rb0,ngenes,1);    
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> trc(Rtrc,n,1);    
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> asn(Rasn,nase,1);    
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> asnp(Rasnp,nase,1);    

    Eigen::MatrixXd k(n,1);
    //expected counts for each individual on log scale - tmn and without log scale etmn
    cppoptlib::Vector<double> tmn(n);  
    cppoptlib::Vector<double> etmn(n);  
        
    Eigen::Map<Eigen::Matrix<int,Eigen::Dynamic,1>> thp(Rthp,nase,1);    
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>> X(RXmatr,n,nb);
    Eigen::DiagonalMatrix<double,Eigen::Dynamic> Wt(n);

    cppoptlib::Vector<double> inits(np);for(int i=0;i<np;i++)inits[i] = Rinits[i];
    cppoptlib::Vector<double> betas(nb);for(int i=nnlin;i<np;i++)betas[i-nnlin] = Rinits[i];
    nlin_liphi[0] = inits[0];
    nlin_lphi[0] = inits[1];
    nlin_b1[0] = inits[2];
    
    cppoptlib::Vector<double> betas0(nb);betas0 = betas;
    cppoptlib::Vector<double> trybetas(nb);trybetas = betas;    
    cppoptlib::Vector<double> delta(nb);delta.setZero(); //prev iter betas
    

    cppoptlib::Matrix<double> pre;    

    double ll1=1e6,ll0=1e6,ll2=1e6,tryll=1e6,ll_old=1e6;
    int flag = 1, flag0 = 1, flag1 = 1;

    int genej = 0;
    //do optim using minimizer
    cppoptlib::TRECASEirw_b1 <double> f_irw;
      opt_liphi = nlin_liphi;
      opt_lphi = nlin_lphi;
      opt_b1 = nlin_b1;
      
      f_irw.setNind(n);
      f_irw.setNase(nase);
      f_irw.setNbet(nb);    
      f_irw.setTRC(trc);
      f_irw.setAS(asn);
      f_irw.setASp(asnp);      
      f_irw.setTRCF(trc);    
      f_irw.setTHP(thp);
      f_irw.setTMN(trc);
      f_irw.setETMN(trc);
      f_irw.setBetas(betas);
      f_irw.setX(X);
      f_irw.setPartmu();
      f_irw.setProb(asn);
      f_irw.setLiphi(nlin_liphi[0]);
      f_irw.setLphi(nlin_lphi[0]);
      f_irw.setB1(nlin_b1[0]);
      f_irw.setB0((double)b0(genej,0));          

      if(estinit){
        k = cppoptlib::vlog1p(trc,n);
        betas = (X.transpose()*X).inverse()*(X.transpose()*k);        
        f_irw.setBetas(betas);
        f_irw.setPartmu();        
      }
      iter = 0;
      f_irw.setDir(4);      //dir = 1 for liphi, dir = 3 for lphi
      ll1 = f_irw(nlin_b1);
      do{
        ll0 = ll1;
        f_irw.setDir(4);
        tmn = f_irw.updTMNb1(nlin_b1);
        //weights are mean/(1.0+mean*phi^-1)
        Wt.diagonal() = cppoptlib::wt(tmn,n,exp(nlin_liphi[0]));//weights are 1/(1+a mu)
        etmn = cppoptlib::vexp(tmn,n);

        k = (trc-etmn).cwiseQuotient(etmn);
        delta = (X.transpose()*Wt*X).inverse()*(X.transpose()*Wt)*k;        
        betas0 = betas;      
        betas = betas + delta;
        f_irw.setBetas(betas);
        f_irw.setPartmu();
        f_irw.setDir(4);
        ll1 = f_irw(nlin_b1);

        //get reasonable initial values of betas assuming initial large overdispersion 
        tryiter = 0;
        flag = 1;
        
        do{
          trybetas = betas0;
          delta /= 2;
          trybetas = trybetas + delta;
          f_irw.setBetas(trybetas);
          f_irw.setPartmu();
          tryll = f_irw(nlin_b1);
          if(tryll>ll1){            
            f_irw.setBetas(betas);
            f_irw.setPartmu();
            flag = 0;
          }else{
            betas = trybetas;
            ll1 = tryll;
          }
          tryiter++;
        }while((tryiter<maxiter2) && (flag==1));
        iter++;
      }while( (iter< maxiter2) && (ll0-ll1)>eps);

      iter = 0;
      opt_b1 = nlin_b1;
      opt_liphi = nlin_liphi;      
      //now update non-linear - recalculate betas iteratively
      do{
        ll0 = ll1;
        //dir = 1 for iphi and dir = 2 for b1
        f_irw.setDir(4);        
        res1 = local_min (nlin_b1[0]+lbnd_b, nlin_b1[0]+ubnd_b, tol, f_irw, opt_b1[0] );
        if(res1 < ll0){
          f_irw.setB1(opt_b1[0]);
        }else{
          opt_b1=nlin_b1;
        }

        f_irw.setDir(1);        
        lbnd_phi = LPHI + opt_liphi[0];
        if(lbnd_phi<LPHIU)lbnd_phi=LPHIU;
        ubnd_phi = UPHI + opt_liphi[0];
        if(ubnd_phi>UPHIU)ubnd_phi=UPHIU;
        res2 = local_min (lbnd_phi, ubnd_phi, tol, f_irw, opt_liphi[0] );
        if(res2 > res1)opt_liphi=nlin_liphi;

        f_irw.setDir(2);        
        lbnd_phi = LPHI + opt_lphi[0];
        if(lbnd_phi<LPHIU)lbnd_phi=LPHIU;
        ubnd_phi = UPHI + opt_lphi[0];
        if(ubnd_phi>UPHIU)ubnd_phi=UPHIU;
        res3 = local_min (lbnd_phi, ubnd_phi, tol, f_irw, opt_lphi[0] );
        if(res3 > res2)opt_lphi=nlin_lphi;

        if((f_irw(opt_lphi)-ll0)>eps){
          //if after optimization we get worse result get back to old overdispersion
          opt_liphi = nlin_liphi;
          opt_lphi = nlin_lphi;
          opt_b1 = nlin_b1;
        }else{
          //otherwise update overdispersion and update mean structure
          nlin_liphi = opt_liphi;
          nlin_lphi = opt_lphi;
          nlin_b1 = opt_b1;
        }
        f_irw.setLiphi(nlin_liphi[0]);
        f_irw.setLphi(nlin_lphi[0]);
        f_irw.setB1(nlin_b1[0]);        

        tmn = f_irw.updTMNb1(nlin_b1);
        //for current overdispersion
        //weights are mean/(1.0+mean*phi^-1)
        Wt.diagonal() = cppoptlib::wt(tmn,n,exp(nlin_liphi[0]));
        etmn = cppoptlib::vexp(tmn,n);
        k = (trc-etmn).cwiseQuotient(etmn);        
        delta = (X.transpose()*Wt*X).inverse()*(X.transpose()*Wt)*k;        
        betas0 = betas;      
        betas = betas + delta;
        f_irw.setBetas(betas);
        f_irw.setPartmu();
        f_irw.setLiphi(nlin_liphi[0]);
        f_irw.setB1(nlin_b1[0]);
        f_irw.setDir(4);
        ll1 = f_irw(nlin_b1);

        //get betas for current overdispersion 
        tryiter = 0;
        flag = 1;
        do{
          trybetas = betas0;
          delta /= 2;
          trybetas = trybetas + delta;
          f_irw.setBetas(trybetas);
          f_irw.setPartmu();
          f_irw.setDir(4);
          tryll = f_irw(nlin_b1);
          if(ll1 - tryll < eps){            
            f_irw.setBetas(betas);
            f_irw.setPartmu();
            flag = 0;
          }else{
            betas = trybetas;
            ll1 = tryll;
          }
          tryiter++;
        }while((delta.norm()>1)&&(tryiter<maxiter2) && (flag==1));
        iter++;
        flag = (int)(((ll0-ll1)>eps)&&(iter<maxiter));        
      }while(flag==1);
      if(ll0<ll1){
        betas=betas0;
        f_irw.setBetas(betas);
        f_irw.setPartmu();      
      }        
      f_irw.setDir(4);
      ll1 = f_irw(nlin_b1);
//    }    
    Rinits[0]=nlin_liphi[0];
    Rinits[1]=nlin_lphi[0];    
    Rinits[2]=nlin_b1[0];    
    for(int i=nnlin;i<np;i++)Rinits[i]=betas[i-nnlin];
    Rinits[np] = ll1;
    
    return ll1;
}

//full TRECASE, fitting b1 - poo
//as.integer(nInd), as.integer(thp),as.double(trc),as.double(asn),as.double(asnp),
//  as.integer(nPar),optval=as.double(inits),as.double(Xmatr), as.integer(est.inits), as.integer(b0)
double max_trecase1irw_b01(int* RnInd, int* RnAse, int* Rngenes, int* Rthp, double* Rtrc, double* Rasn, double* Rasnp, 
                          int* RnPar, double* Rinits, double* RXmatr, int* Restinit) {
    int maxiter = 100,maxiter2=50,iter,iter2,tryiter;
    int eps = 1e-8;
    double tol = 1e-20;
    double lbnd_phi = LPHI; //lower boundary for log(phi^-1)
    double ubnd_phi = UPHI; //upper boundary for log(phi^-1)
    double lbnd_b = LBND2; //lower boundary for genetic effect
    double ubnd_b = UBND2; //upper boundary for genetic effect
    double res1,res2,res3;
    
    const int n = RnInd[0];//number of individuals
    const int nase = RnAse[0];
    const int np = RnPar[0];//number of param
    const int ngenes = Rngenes[0];//number of genes, 1 gene for now    
    const int estinit = Restinit[0];
    const int nnlin = 4;//number of nonlinear par
    const int nb = np-nnlin;
    cppoptlib::Vector<double> opt_liphi(1);
    cppoptlib::Vector<double> opt_lphi(1);
    cppoptlib::Vector<double> opt_b0(1);    
    cppoptlib::Vector<double> opt_b1(1);    
    cppoptlib::Vector<double> nlin_liphi(1);
    cppoptlib::Vector<double> nlin_lphi(1);
    cppoptlib::Vector<double> nlin_b0(1);
    cppoptlib::Vector<double> nlin_b1(1);

    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> trc(Rtrc,n,1);    
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> asn(Rasn,nase,1);    
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>> asnp(Rasnp,nase,1);    

    Eigen::MatrixXd k(n,1);
    //expected counts for each individual on log scale - tmn and without log scale etmn
    cppoptlib::Vector<double> tmn(n);  
    cppoptlib::Vector<double> etmn(n);  
        
    Eigen::Map<Eigen::Matrix<int,Eigen::Dynamic,1>> thp(Rthp,nase,1);    
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>> X(RXmatr,n,nb);
    Eigen::DiagonalMatrix<double,Eigen::Dynamic> Wt(n);

    cppoptlib::Vector<double> inits(np);for(int i=0;i<np;i++)inits[i] = Rinits[i];
    cppoptlib::Vector<double> betas(nb);for(int i=nnlin;i<np;i++)betas[i-nnlin] = Rinits[i];
    nlin_liphi[0] = inits[0];
    nlin_lphi[0] = inits[1];
    nlin_b0[0] = inits[2];
    nlin_b1[0] = inits[3];
    
    cppoptlib::Vector<double> betas0(nb);betas0 = betas;
    cppoptlib::Vector<double> trybetas(nb);trybetas = betas;    
    cppoptlib::Vector<double> delta(nb);delta.setZero(); //prev iter betas
    

    cppoptlib::Matrix<double> pre;    

    double ll1=1e6,ll0=1e6,ll2=1e6,tryll=1e6,ll_old=1e6;
    int flag = 1, flag0 = 1, flag1 = 1;
    double curlbnd,curubnd;
    int genej = 0;
    //do optim using minimizer
    cppoptlib::TRECASEirw_b01 <double> f_irw;
      opt_liphi = nlin_liphi;
      opt_lphi = nlin_lphi;
      opt_b0 = nlin_b0;
      opt_b1 = nlin_b1;
      
      f_irw.setNind(n);
      f_irw.setNase(nase);
      f_irw.setNbet(nb);    
      f_irw.setTRC(trc);
      f_irw.setAS(asn);
      f_irw.setASp(asnp);      
      f_irw.setTRCF(trc);    
      f_irw.setTHP(thp);
      f_irw.setTMN(trc);
      f_irw.setETMN(trc);
      f_irw.setBetas(betas);
      f_irw.setX(X);
      f_irw.setPartmu();
      f_irw.setProb(asn);
      f_irw.setLiphi(nlin_liphi[0]);
      f_irw.setLphi(nlin_lphi[0]);
      f_irw.setB0(nlin_b0[0]);
      f_irw.setB1(nlin_b1[0]);
      cppoptlib::BfgsSolver<double> solver;
      int rotflag = 1;
      cppoptlib::Vector<double> nlin_b0_b1(2);
      cppoptlib::Vector<double> opt_b0_b1(2);
      nlin_b0_b1[0] = nlin_b0[0];
      nlin_b0_b1[1] = nlin_b1[0];
      opt_b0_b1 = nlin_b0_b1;
      
      f_irw.setDir(4);      //dir = 1 for liphi, dir = 2 for lphi
      ll1 = f_irw(nlin_b1);

      if(estinit){
        k = cppoptlib::vlog1p(trc,n);
        betas = (X.transpose()*X).inverse()*(X.transpose()*k);        
        f_irw.setBetas(betas);
        f_irw.setPartmu();        
      }
      iter = 0;

      f_irw.setDir(4);      //dir = 1 for liphi, dir = 2 for lphi
      ll1 = f_irw(nlin_b1);
      do{
        ll0 = ll1;
        f_irw.setDir(4);
        tmn = f_irw.updTMNb1(nlin_b1);
        //weights are mean/(1.0+mean*phi^-1)
        Wt.diagonal() = cppoptlib::wt(tmn,n,exp(nlin_liphi[0]));//weights are 1/(1+a mu)
        etmn = cppoptlib::vexp(tmn,n);

        k = (trc-etmn).cwiseQuotient(etmn);
        delta = (X.transpose()*Wt*X).inverse()*(X.transpose()*Wt)*k;        
        betas0 = betas;      
        betas = betas + delta;
        f_irw.setBetas(betas);
        f_irw.setPartmu();
        f_irw.setDir(4);
        ll1 = f_irw(nlin_b1);

        //get reasonable initial values of betas assuming initial large overdispersion 
        tryiter = 0;
        flag = 1;
        
        do{
          trybetas = betas0;
          delta /= 2;
          trybetas = trybetas + delta;
          f_irw.setBetas(trybetas);
          f_irw.setPartmu();
          tryll = f_irw(nlin_b1);
          if(tryll>ll1){            
            f_irw.setBetas(betas);
            f_irw.setPartmu();
            flag = 0;
          }else{
            betas = trybetas;
            ll1 = tryll;
          }
          tryiter++;
        }while((tryiter<maxiter2) && (flag==1));
        iter++;
      }while( (iter< maxiter2) && (ll0-ll1)>eps);

      iter = 0;
      opt_b1 = nlin_b1;
      opt_liphi = nlin_liphi;      
      //now update non-linear - recalculate betas iteratively
      do{
        ll0 = ll1;
//        std::cout << "init ll:" << ll0 << " ";
        //dir = 1 for iphi and dir = 2 for phi, dir = 3 for b0, dir = 4 for b1
        if(rotflag){
          f_irw.setDir(5);
//        std::cout << opt_b0_b1.transpose() << " " <<f_irw(opt_b0_b1) << " ";
          solver.minimize(f_irw, opt_b0_b1);
          ll1 = f_irw(opt_b0_b1);
//          std::cout << ll1 << " at " << opt_b0_b1.transpose() << " ";
//          rotflag=0;          
          nlin_b0_b1=opt_b0_b1;
          nlin_b0[0] = nlin_b0_b1[0];          
          nlin_b1[0] = nlin_b0_b1[1];          
          opt_b0[0] = nlin_b0_b1[0];          
          opt_b1[0] = nlin_b0_b1[1];
          res1 = ll1;
        }else{
          f_irw.setDir(3);        
          res3 = local_min (nlin_b0[0]+lbnd_b, nlin_b0[0]+ubnd_b, tol, f_irw, opt_b0[0] );
          if(res1 < ll0){
            f_irw.setB0(opt_b0[0]);
          }else{
            opt_b0=nlin_b0;
          }
          f_irw.setDir(4);        
          res1 = local_min (nlin_b1[0]+lbnd_b, nlin_b1[0]+ubnd_b, tol, f_irw, opt_b1[0] );
          if(res1 < ll0){
            f_irw.setB1(opt_b1[0]);
          }else{
            opt_b1=nlin_b1;
          }
//          rotflag=1;
          nlin_b0_b1[0] = nlin_b0[0];          
          nlin_b0_b1[0] = nlin_b0[0];
        }
        f_irw.setDir(1);        
//        curlbnd = max(opt_liphi[0]+lbnd_phi,
        lbnd_phi = LPHI + opt_liphi[0];
        if(lbnd_phi<LPHIU)lbnd_phi=LPHIU;
        ubnd_phi = UPHI + opt_liphi[0];
        if(ubnd_phi<UPHIU)ubnd_phi=UPHIU;
        res2 = local_min (lbnd_phi, ubnd_phi, tol, f_irw, opt_liphi[0] );
//        std::cout << "liphi: " << res2 << " at " << opt_liphi[0] << " ";
        if(res2 > res1)opt_liphi=nlin_liphi;

        f_irw.setDir(2);        
        lbnd_phi = LPHI + opt_lphi[0];
        if(lbnd_phi<LPHIU)lbnd_phi=LPHIU;
        ubnd_phi = UPHI + opt_lphi[0];
        if(ubnd_phi>UPHIU)ubnd_phi=UPHIU;
        res3 = local_min (lbnd_phi, ubnd_phi, tol, f_irw, opt_lphi[0] );
//        std::cout << "lphi: " << res3 << " at " << opt_lphi[0] << std:endl;
        if(res3 > res1)opt_lphi=nlin_lphi;

        if((f_irw(opt_lphi)-ll0)>eps){
          //if after optimization we get worse result get back to old overdispersion
          opt_liphi = nlin_liphi;
          opt_lphi = nlin_lphi;
          opt_b0 = nlin_b0;
          opt_b1 = nlin_b1;
        }else{
          //otherwise update overdispersion and update mean structure
          nlin_liphi = opt_liphi;
          nlin_lphi = opt_lphi;
          nlin_b0 = opt_b0;
          nlin_b1 = opt_b1;
        }
        f_irw.setLiphi(nlin_liphi[0]);
        f_irw.setLphi(nlin_lphi[0]);
        f_irw.setB0(nlin_b0[0]);        
        f_irw.setB1(nlin_b1[0]);        

        f_irw.setDir(4);
        tmn = f_irw.updTMNb1(nlin_b1);
        //for current overdispersion
        //weights are mean/(1.0+mean*phi^-1)
        Wt.diagonal() = cppoptlib::wt(tmn,n,exp(nlin_liphi[0]));
        etmn = cppoptlib::vexp(tmn,n);
        k = (trc-etmn).cwiseQuotient(etmn);        
        delta = (X.transpose()*Wt*X).inverse()*(X.transpose()*Wt)*k;        
        betas0 = betas;      
        betas = betas + delta;
        f_irw.setBetas(betas);
        f_irw.setPartmu();
        f_irw.setLiphi(nlin_liphi[0]);
        f_irw.setLphi(nlin_lphi[0]);
        f_irw.setB0(nlin_b0[0]);
        f_irw.setB1(nlin_b1[0]);
        f_irw.setDir(4);
        ll1 = f_irw(nlin_b1);

        //get betas for current overdispersion 
        tryiter = 0;
        flag = 1;
        do{
          trybetas = betas0;
          delta /= 2;
          trybetas = trybetas + delta;
          f_irw.setBetas(trybetas);
          f_irw.setPartmu();
          f_irw.setDir(4);
          tryll = f_irw(nlin_b1);
          if(ll1 - tryll < eps){            
            f_irw.setBetas(betas);
            f_irw.setPartmu();
            flag = 0;
          }else{
            betas = trybetas;
            ll1 = tryll;
          }
          tryiter++;
        }while((delta.norm()>1)&&(tryiter<maxiter2) && (flag==1));
        iter++;
        flag = (int)(((ll0-ll1)>eps)&&(iter<maxiter));        
      }while(flag==1);
      if(ll0<ll1){
        betas=betas0;
        f_irw.setBetas(betas);
        f_irw.setPartmu();      
      }        
      f_irw.setDir(4);
      ll1 = f_irw(nlin_b1);
//    }    
    //std::cout << nlin_b0 << " " << opt_b0  << " " << f_irw.getB0() << " " << nlin_b1 << " " << opt_b1  << " " << f_irw.getB1() << std::endl;
    Rinits[0]=nlin_liphi[0];
    Rinits[1]=nlin_lphi[0];    
    Rinits[2]=nlin_b0[0];    
    Rinits[3]=nlin_b1[0];    
    for(int i=nnlin;i<np;i++)Rinits[i]=betas[i-nnlin];
    Rinits[np] = ll1;

    return ll1;
}
}
 
