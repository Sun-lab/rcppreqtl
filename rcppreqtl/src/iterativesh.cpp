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

// to use this library just use the namespace "cppoptlib"
namespace cppoptlib {

template<typename T>
class TRECASEirwsh_no : public Problem<T> {  
  public:      
    // this is just the objective (NOT optional)
    T value(const Vector<T> &x) {
      double out = 0, out2 = 0;  
      double liphi,lphi,b0;

      liphi = this->liphi;
      lphi = -this->liphi;
      b0 = this->b0;

      switch ( this->dir ) {
      case 1:
      case 2:
        liphi = x[0];
        lphi = -x[0];
        break;
      default:
        b0 = 0;
        break;
      }      

      double vals[4] = {0,0,0,0};
      //vals[1] = log1p(exp(b0))-log1p(1);
      //vals[2] = vals[1];
      //vals[3] = b0;
  
      for(int i=0;i<this->nInd;i++){
        this->tmn[i] = this->partmu[i]; + vals[this->thp[i]];
      }    
      
      out = TREC_lik(this->trc,this->trcf,this->tmn,liphi,this->nInd);
      //calc prob using given b0
      //vals[0] = 0;  //b1 in full model
      //vals[1] = b0; //b0+b1 in full model;
      //vals[2] = -b0;//-b0+b1 in full model
      //vals[3] = 0;  //b1 in full model
      //for(int i=0;i<4;i++){
      //  vals[i] = 1/(1+exp(-vals[i]));
      //}
      for(int i=0;i<this->nAse;i++){
      //  this->prob[i] = vals[this->thp[i]];
        this->prob[i] = 0.5;
      }          
      
      out2 = ASE_lik(this->asn, this->asnp, this->prob, lphi, this->nAse);
//      out2 = ASE_lik2(this->asn, this->asnp, this->prob, lphi, this->nAse);
      out += out2;
      return out;
    }
};

template<typename T>
class TRECASEirwsh_b0 : public Problem<T> {  
  public:      
    // this is just the objective (NOT optional)
    T value(const Vector<T> &x) {
      double out = 0, out2 = 0;  
      double liphi,lphi,b0;

      liphi = this->liphi;
      lphi = -this->liphi;
      b0 = this->b0;

      switch ( this->dir ) {
      case 1:
      case 2:
        liphi = x[0];
        lphi = -x[0];
        break;
      default:
        b0 = x[0];
        break;
      }      

      double vals[4] = {0,0,0,0};
      vals[1] = log1p(exp(b0))-log1p(1);
      vals[2] = vals[1];
      vals[3] = b0;
  
      for(int i=0;i<this->nInd;i++){
        this->tmn[i] = this->partmu[i] + vals[this->thp[i]];
      }    
      
      out = TREC_lik(this->trc,this->trcf,this->tmn, liphi,this->nInd);
      //calc prob using given b0
      vals[0] = 0;  //b1 in full model
      vals[1] = b0; //b0+b1 in full model;
      vals[2] = -b0;//-b0+b1 in full model
      vals[3] = 0;  //b1 in full model
      for(int i=0;i<4;i++){
        vals[i] = 1/(1+exp(-vals[i]));
      }
      for(int i=0;i<this->nAse;i++){
        this->prob[i] = vals[this->thp[i]];
      }          
      
      out2 = ASE_lik(this->asn, this->asnp, this->prob, lphi, this->nAse);
//      out2 = ASE_lik2(this->asn, this->asnp, this->prob, lphi, this->nAse);
      out += out2;
      return out;
    }
};
template<typename T>
class TRECASEirwsh_b1 : public Problem<T> {  
  public:      
    // this is just the objective (NOT optional)
    T value(const Vector<T> &x) {
      double out = 0, out2 = 0;  
      double liphi,lphi,b0, b1;

      liphi = this->liphi;
      lphi = -this->liphi;
      b0 = this->b0;
      b1 = this->b1;

      switch ( this->dir ) {
      case 1:
      case 2:
        liphi = x[0];
        lphi = -x[0];
        break;
      default:
        b1 = x[0];
        break;
      }

      double vals[4] = {0,0,0,0};
      vals[1] = log1p(exp(b0+b1))-log1p(exp(b1));
      vals[2] = log1p(exp(b0-b1))-log1p(exp(-b1));
      vals[3] = b0;
  
      for(int i=0;i<this->nInd;i++){
        this->tmn[i] = this->partmu[i] + vals[this->thp[i]];
      }    
      
      out = TREC_lik(this->trc,this->trcf,this->tmn,liphi,this->nInd);
      //calc prob using given b0,b1
      vals[0] = b1;  //b1 in full model
      vals[1] =  b0 + b1; //b0+b1 in full model;
      vals[2] = -b0 + b1;//-b0+b1 in full model
      vals[3] = b1;  //b1 in full model
      for(int i=0;i<4;i++){
        vals[i] = 1/(1+exp(-vals[i]));
      }
      for(int i=0;i<this->nAse;i++){
        this->prob[i] = vals[this->thp[i]];
      }          
      
      out2 = ASE_lik(this->asn, this->asnp, this->prob, lphi, this->nAse);
//      std::cout << liphi << " " << b0 << " " << b1 << " " << out << " " << out2 << std::endl;
//      out2 = ASE_lik2(this->asn, this->asnp, this->prob, lphi, this->nAse);
      out += out2;
      return out;
    }
};
template<typename T>
class TRECASEirwsh_b01 : public Problem<T> {  
  public:      
    // this is just the objective (NOT optional)
    T value(const Vector<T> &x) {
      double out = 0, out2 = 0;  
      double liphi, lphi, b0, b1;

      liphi = this->liphi;        
      lphi = -this->liphi;        
      b0 = this->b0;
      b1 = this->b1;
      switch ( this->dir ) {
      case 1:
      case 2:
        liphi = x[0];
        lphi = -x[0];
        break;
      case 3:
        b0 = x[0];
        break;
      case 4:
        b1 = x[0];
        break;
      default:
        b0 = x[0];
        b1 = x[1];
        break;
      }      
      double vals[4] = {0,0,0,0};
      vals[1] = log1p(exp(b0+b1))-log1p(exp(b1));
      vals[2] = log1p(exp(b0-b1))-log1p(exp(-b1));
      vals[3] = b0;
  
      for(int i=0;i<this->nInd;i++){
        this->tmn[i] = this->partmu[i] + vals[this->thp[i]];
      }    
      
      out = TREC_lik(this->trc, this->trcf, this->tmn, liphi, this->nInd);
      //calc prob using given b0,b1
/*      
      vals[0] = b1;  //b1 in full model
      vals[1] =  b0 + b1; //b0+b1 in full model;
      vals[2] = -b0 + b1;//-b0+b1 in full model
      vals[3] = b1;  //b1 in full model
*/      
      vals[0] = b1;  //b1 in full model
      vals[1] =  b0 + b1; //b0+b1 in full model;
      vals[2] = -b0 + b1;//-b0+b1 in full model
      vals[3] = b1;  //b1 in full model

      for(int i=0;i<4;i++){
        vals[i] = 1/(1+exp(-vals[i]));
      }
      for(int i=0;i<this->nAse;i++){
        this->prob[i] = vals[this->thp[i]];
      }          
      
      out2 = ASE_lik(this->asn, this->asnp, this->prob, lphi, this->nAse);
//      out2 = ASE_lik2(this->asn, this->asnp, this->prob, lphi, this->nAse);
      out += out2;
//      if((this->dir)==5){
//      std::cout << liphi << " " << lphi << " " << x.transpose() << " " << (this->betas).transpose() << " " << out << std::endl;
//      }
      return out;
    }
};
}

extern "C" {


//full TRECASE, fitting b0 - additive
//as.integer(nInd), as.integer(nGenes), as.integer(thp),as.double(trc),as.double(asn),as.double(asnp),
//  as.integer(nPar),optval=as.double(inits),as.double(Xmatr), as.integer(est.inits)

double max_trecase1irwsh_b0(int* RnInd, int* RnAse, int* Rngenes, int* Rthp, double* Rtrc, double* Rasn, double* Rasnp, 
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
    nlin_lphi[0] =  -inits[0];
    nlin_b0[0] = inits[2];
    
    cppoptlib::Vector<double> betas0(nb);betas0 = betas;
    cppoptlib::Vector<double> trybetas(nb);trybetas = betas;    
    cppoptlib::Vector<double> delta(nb);delta.setZero(); //prev iter betas
    

    cppoptlib::Matrix<double> pre;    

    double ll1=1e6,ll0=1e6,ll2=1e6,tryll=1e6,ll_old=1e6;
    int flag = 1, flag0 = 1, flag1 = 1;

    //do optim using minimizer
    cppoptlib::TRECASEirwsh_b0 <double> f_irw;
      opt_liphi = nlin_liphi;
      opt_lphi =  -nlin_liphi;
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
      f_irw.setLphi(-nlin_liphi[0]);

      if(estinit){
        k = cppoptlib::vlog1p(trc,n);
        betas = (X.transpose()*X).inverse()*(X.transpose()*k);        
        f_irw.setBetas(betas);
        f_irw.setPartmu();        
      }
      iter = 0;
      f_irw.setDir(3);      //dir = 1 for liphi, dir = 2 for lphi
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
        if(res2 > res1){
          opt_liphi=nlin_liphi;
          opt_lphi= -nlin_liphi;
        }else{
          opt_lphi= -opt_liphi;
        }

        if((f_irw(opt_liphi)-ll0)>eps){
          //if after optimization we get worse result get back to old overdispersion
          opt_liphi = nlin_liphi;
          opt_lphi =  -nlin_liphi;
          opt_b0 = nlin_b0;
        }else{
          //otherwise update overdispersion and update mean structure
          nlin_liphi = opt_liphi;
          nlin_lphi =  -opt_liphi;
          nlin_b0 = opt_b0;
        }
        f_irw.setLiphi(nlin_liphi[0]);
        f_irw.setLphi( -nlin_liphi[0]);
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
        f_irw.setLphi( -nlin_liphi[0]);
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
    Rinits[1]=-nlin_liphi[0];    
    Rinits[2]=nlin_b0[0];    
    for(int i=nnlin;i<np;i++)Rinits[i]=betas[i-nnlin];
    Rinits[np] = ll1;

    return ll1;
}


//full TRECASE, fitting b1 - poo
//as.integer(nInd), as.integer(thp),as.double(trc),as.double(asn),as.double(asnp),
//  as.integer(nPar),optval=as.double(inits),as.double(Xmatr), as.integer(est.inits), as.integer(b0)
double max_trecase1irwsh_b1(int* RnInd, int* RnAse, int* Rngenes, int* Rthp, double* Rtrc, double* Rasn, double* Rasnp, 
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
    nlin_lphi[0] =  -inits[0];
    nlin_b1[0] = inits[2];
    
    cppoptlib::Vector<double> betas0(nb);betas0 = betas;
    cppoptlib::Vector<double> trybetas(nb);trybetas = betas;    
    cppoptlib::Vector<double> delta(nb);delta.setZero(); //prev iter betas
    

    cppoptlib::Matrix<double> pre;    

    double ll1=1e6,ll0=1e6,ll2=1e6,tryll=1e6,ll_old=1e6;
    int flag = 1, flag0 = 1, flag1 = 1;

    int genej = 0;
    //do optim using minimizer
    cppoptlib::TRECASEirwsh_b1 <double> f_irw;
      opt_liphi = nlin_liphi;
      opt_lphi =  -nlin_liphi;
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
      f_irw.setLphi( -nlin_liphi[0]);
      f_irw.setB1(nlin_b1[0]);
      f_irw.setB0((double)b0(genej,0));          

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
        if(res2 > res1){
          opt_liphi=nlin_liphi;
          opt_lphi= -nlin_liphi;
        }else{
          opt_lphi= -opt_liphi;
        }

        if((f_irw(opt_liphi)-ll0)>eps){
          //if after optimization we get worse result get back to old overdispersion
          opt_liphi = nlin_liphi;
          opt_lphi =  -nlin_liphi;
          opt_b1 = nlin_b1;
        }else{
          //otherwise update overdispersion and update mean structure
          nlin_liphi = opt_liphi;
          nlin_lphi =  -opt_liphi;
          nlin_b1 = opt_b1;
        }
        f_irw.setLiphi(nlin_liphi[0]);
        f_irw.setLphi( -nlin_liphi[0]);
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
        f_irw.setLphi( -nlin_liphi[0]);
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
    Rinits[1]=-nlin_liphi[0];    
    Rinits[2]=nlin_b1[0];    
    for(int i=nnlin;i<np;i++)Rinits[i]=betas[i-nnlin];
    Rinits[np] = ll1;
    
    return ll1;
}

//full TRECASE, fitting b1 - poo
//as.integer(nInd), as.integer(thp),as.double(trc),as.double(asn),as.double(asnp),
//  as.integer(nPar),optval=as.double(inits),as.double(Xmatr), as.integer(est.inits), as.integer(b0)
double max_trecase1irwsh_b01(int* RnInd, int* RnAse, int* Rngenes, int* Rthp, double* Rtrc, double* Rasn, double* Rasnp, 
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
    nlin_lphi[0] =  -inits[0];
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
    cppoptlib::TRECASEirwsh_b01 <double> f_irw;
      opt_liphi = nlin_liphi;
      opt_lphi =  -nlin_liphi;
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
      f_irw.setLphi( -nlin_liphi[0]);
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
        if(res2 > res1){
          opt_liphi=nlin_liphi;
          opt_lphi= -nlin_liphi;
        }else{
          opt_lphi= -opt_liphi;
        }


        if((f_irw(opt_liphi)-ll0)>eps){
          //if after optimization we get worse result get back to old overdispersion
          opt_liphi = nlin_liphi;
          opt_lphi =  -nlin_liphi;
          opt_b0 = nlin_b0;
          opt_b1 = nlin_b1;
        }else{
          //otherwise update overdispersion and update mean structure
          nlin_liphi = opt_liphi;
          nlin_lphi =  -opt_liphi;
          nlin_b0 = opt_b0;
          nlin_b1 = opt_b1;
        }
        f_irw.setLiphi(nlin_liphi[0]);
        f_irw.setLphi( -nlin_liphi[0]);
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
        f_irw.setLphi( -nlin_liphi[0]);
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
    Rinits[1]=-nlin_liphi[0];    
    Rinits[2]=nlin_b0[0];    
    Rinits[3]=nlin_b1[0];    
    for(int i=nnlin;i<np;i++)Rinits[i]=betas[i-nnlin];
    Rinits[np] = ll1;

    return ll1;
}




//no b0
double max_trecase1irwsh(int* RnInd, int* RnAse, int* Rngenes, int* Rthp, double* Rtrc, double* Rasn, double* Rasnp, 
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
    inits[2] = 0;
    cppoptlib::Vector<double> betas(nb);for(int i=nnlin;i<np;i++)betas[i-nnlin] = Rinits[i];
    nlin_liphi[0] = inits[0];
    nlin_lphi[0] = inits[1];
    nlin_b0[0] = inits[2];
    //std::cout << nlin_liphi << " " << nlin_lphi << std::endl;

    cppoptlib::Vector<double> betas0(nb);betas0 = betas;
    cppoptlib::Vector<double> trybetas(nb);trybetas = betas;    
    cppoptlib::Vector<double> delta(nb);delta.setZero(); //prev iter betas
    

    cppoptlib::Matrix<double> pre;    

    double ll1=1e6,ll0=1e6,ll2=1e6,tryll=1e6,ll_old=1e6;
    int flag = 1, flag0 = 1, flag1 = 1;

    //do optim using minimizer
    //cppoptlib::TRECASEirw_b0 <double> f_irw;
    cppoptlib::TRECASEirwsh_no <double> f_irw;
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
        //f_irw.setDir(3);        
        //res1 = local_min (nlin_b0[0]+lbnd_b, nlin_b0[0]+ubnd_b, tol, f_irw, opt_b0[0] );
          res1 = ll0;
        //if(res1 < ll0){
        //  f_irw.setB0(opt_b0[0]);
        //}else{
          opt_b0=nlin_b0;
        //}


        f_irw.setDir(1);        
        lbnd_phi = LPHI + opt_liphi[0];
        if(lbnd_phi<LPHIU)lbnd_phi=LPHIU;
        ubnd_phi = UPHI + opt_liphi[0];
        if(ubnd_phi>UPHIU)ubnd_phi=UPHIU;
        res2 = local_min (lbnd_phi, ubnd_phi, tol, f_irw, opt_liphi[0] );

        if(res2 > res1){
          opt_liphi=nlin_liphi;
          opt_lphi= -nlin_liphi;
        }else{
          opt_liphi= opt_liphi;
          opt_lphi= -opt_liphi;
        }
        
        if((f_irw(opt_liphi)-ll0)>eps){
          //if after optimization we get worse result get back to old overdispersion
          opt_liphi = nlin_liphi;
          opt_lphi =  -nlin_liphi;
          opt_b0 = nlin_b0;
        }else{
          //otherwise update overdispersion and update mean structure
          nlin_liphi = opt_liphi;
          nlin_lphi =  -opt_liphi;
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
 
}
 
