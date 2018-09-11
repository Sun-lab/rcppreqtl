
//#define LBND -20
//#define UBND 20
#define LBND -10
#define UBND 10
#define LBND2 -10
#define UBND2 10
#define LPHI -4e0
#define UPHI 4e0
#define LPHIU -8e0
#define UPHIU 8e0
#define LPHI2 1e-3
#define UPHI2 1e3
#define MLIK 1e10
#define N 1e0


namespace cppoptlib {
double ASE_lik(Vector<double> asn, Vector<double> asnp, Vector<double> prob, double lphi, int n);
double ASE_lik2(Vector<double> asn, Vector<double> asnp, Vector<double> prob, double lphi, int n);
double TREC_lik(Vector<double> y, Vector<double> yf, Vector<double> lmu, double liphi, int n);
Vector<double> wt(Vector<double> lmu, int n, double iphi);
Vector<double> vexp(Vector<double> v, int n);
Vector<double> vlog1p(Vector<double> v, int n);

double max_trec1irw0(int* RnInd, int* Rngenes, int* Rthp, double* Rtrc, int* RnPar, double* Rinits, double* RXmatr, int* Restinit);
double max_trec1irw_b0(int* RnInd, int* Rngenes, int* Rthp, double* Rtrc, int* RnPar, double* Rinits, double* RXmatr, int* Restinit);
double max_trec1irw_b1(int* RnInd, int* Rngenes, int* Rthp, double* Rtrc, int* RnPar, double* Rinits, double* RXmatr, int* Restinit, double* Rb0);

double max_trecase1irw_b0(int* RnInd, int* RnAse, int* Rngenes, int* Rthp, double* Rtrc, double* Rasn, double* Rasnp, 
                          int* RnPar, double* Rinits, double* RXmatr, int* Restinit);
double max_trecase1irw_b1(int* RnInd, int* RnAse, int* Rngenes, int* Rthp, double* Rtrc, double* Rasn, double* Rasnp, 
                          int* RnPar, double* Rinits, double* RXmatr, int* Restinit,  double* Rb0);
double max_trecase1irw_b01(int* RnInd, int* RnAse, int* Rngenes, int* Rthp, double* Rtrc, double* Rasn, double* Rasnp, 
                          int* RnPar, double* Rinits, double* RXmatr, int* Restinit);

double max_trecase1irwsh_b0(int* RnInd, int* RnAse, int* Rngenes, int* Rthp, double* Rtrc, double* Rasn, double* Rasnp, 
                          int* RnPar, double* Rinits, double* RXmatr, int* Restinit);
double max_trecase1irwsh_b1(int* RnInd, int* RnAse, int* Rngenes, int* Rthp, double* Rtrc, double* Rasn, double* Rasnp, 
                          int* RnPar, double* Rinits, double* RXmatr, int* Restinit,  double* Rb0);
double max_trecase1irwsh_b01(int* RnInd, int* RnAse, int* Rngenes, int* Rthp, double* Rtrc, double* Rasn, double* Rasnp, 
                          int* RnPar, double* Rinits, double* RXmatr, int* Restinit);

template<typename T>
class TRECASEirw_b0 : public Problem<T> {  
  public:      
    // this is just the objective (NOT optional)
    T value(const Vector<T> &x) {
      double out = 0, out2 = 0;  
      double liphi,lphi,b0;

      liphi = this->liphi;
      lphi = this->lphi;
      b0 = this->b0;

      switch ( this->dir ) {
      case 1:
        liphi = x[0];
        break;
      case 2:
        lphi = x[0];
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
      
      out = TREC_lik(this->trc,this->trcf,this->tmn,liphi,this->nInd);
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
class TRECASEirw_b1 : public Problem<T> {  
  public:      
    // this is just the objective (NOT optional)
    T value(const Vector<T> &x) {
      double out = 0, out2 = 0;  
      double liphi,lphi,b1;

      liphi = this->liphi;
      lphi = this->lphi;
      b1 = this->b1;

      switch ( this->dir ) {
      case 1:
        liphi = x[0];
        break;
      case 2:
        lphi = x[0];
        break;
      default:
        b1 = x[0];
        break;
      }

      double vals[4] = {0,0,0,0};
      vals[1] = log1p(exp(this->b0+b1))-log1p(exp(b1));
      vals[2] = log1p(exp(this->b0-b1))-log1p(exp(-b1));
      vals[3] = this->b0;
  
      for(int i=0;i<this->nInd;i++){
        this->tmn[i] = this->partmu[i] + vals[this->thp[i]];
      }    
      
      out = TREC_lik(this->trc,this->trcf,this->tmn,liphi,this->nInd);
      //calc prob using given b0,b1
      vals[0] = b1;  //b1 in full model
      vals[1] =  this->b0 + b1; //b0+b1 in full model;
      vals[2] = - (this->b0) + b1;//-b0+b1 in full model
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
      return out;
    }
};
template<typename T>
class TRECASEirw_b01 : public Problem<T> {  
  public:      
    // this is just the objective (NOT optional)
    T value(const Vector<T> &x) {
      double out = 0, out2 = 0;  
      double liphi, lphi, b0, b1;

      liphi = this->liphi;        
      lphi = this->lphi;        
      b0 = this->b0;
      b1 = this->b1;
      switch ( this->dir ) {
      case 1:
        liphi = x[0];
        break;
      case 2:
        lphi = x[0];
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

template<typename T>
class TRECirw0 : public Problem<T> {
  public:
    double out = 0;  
    // this is just the objective (NOT optional)
    T value(const Vector<T> &x) {
    for(int i=0;i<this->nInd;i++){
      this->tmn[i] = this->partmu[i];
    }
    out = TREC_lik(this->trc,this->trcf,this->tmn,x[0],this->nInd);
    return out;
    }
};


template<typename T>
class TRECirw_b0 : public Problem<T> {  
  public:      
    // this is just the objective (NOT optional)
    T value(const Vector<T> &x) {
      double out = 0;  
      double liphi,b0;

      liphi = this->liphi;        
      b0 = this->b0;
      
      switch ( this->dir ) {
      case 1:
        liphi = x[0];
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
      
      out = TREC_lik(this->trc,this->trcf,this->tmn,liphi,this->nInd);
      return out;
    }
};


/*
  ends <- -log1p(exp(-b1f))
  etas[ind.lst[[4]]] <-  b0f #BxB, fem
  etas[ind.lst[[1]]] <- -b1f  + log1p(exp(b0f + b1f)) -log1p(exp(-b1f))#AxB, fem
  etas[ind.lst[[2]]] <-  log1p(exp(b0f - b1f)) -log1p(exp(-b1f))#BxA, fem
*/
template<typename T>
class TRECirw_b1 : public Problem<T> {
  public:
    
    // this is just the objective (NOT optional)
    T value(const Vector<T> &x) {
      double out = 0;  
      double liphi,b1;

      liphi = this->liphi;        
      b1 = this->b1;
      
      switch ( this->dir ) {
      case 1:
        liphi = x[0];
        break;
      default:
        b1 = x[0];
        break;
      }
//    double b1 = x[1];
    double vals[4] = {0,0,0,0};
    vals[1] = log1p(exp(this->b0+b1))-log1p(exp(b1));
    vals[2] = log1p(exp(this->b0-b1))-log1p(exp(-b1));
    vals[3] = this->b0;

    for(int i=0;i<this->nInd;i++){
      this->tmn[i] = this->partmu[i] + vals[this->thp[i]];
    }    
    out = TREC_lik(this->trc,this->trcf,this->tmn,liphi,this->nInd);
    return out;
    }
};
}
