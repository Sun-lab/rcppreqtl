# include <string.h>
# include <ctype.h>
# include <stdio.h>   
# include <stdlib.h>  
# include <math.h>  
# include <time.h>
# include <iomanip>
# include "function.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////

void nrerror(const char* error_text)
{
    fprintf(stderr,"Numerical Recipes run-time error...\n");
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    exit(1);
}
  
int *ivector(int nl, int nh)
{
    int *v;

    v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
    if (!v) nrerror("allocation failure in ivector()");
    return v-nl;
}

char *cvector(int nl, int nh) {
    char *v;

    v=(char *)malloc((unsigned) (nh-nl+1)*sizeof(char));
    if (!v) nrerror("allocation failure in cvector()");
    return v-nl;
}

double *dvector(int nl, int nh)
{
    double *v;

    v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
    if (!v) nrerror("allocation failure in vector()");
    return v-nl;
}
 
char **cmatrix(int nrl, int nrh, int ncl, int nch)
{
    int i;
	char **m;

    m=(char **)malloc((unsigned) (nrh-nrl+1)*sizeof(char*));
    if (!m) nrerror("allocation failure 1 in cmatrix()");
    m -= nrl;

    for(i=nrl;i<=nrh;i++) {
        m[i]=(char *)malloc((unsigned) (nch-ncl+1)*sizeof(char));
        if (!m[i]) nrerror("allocation failure 2 in imatrix()");
        m[i] -= ncl;
    }
    return m;
}

double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
    int i;
    double **m;

    m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
    if (!m) nrerror("allocation failure 1 in dmatrix()");
    m -= nrl;

    for(i=nrl;i<=nrh;i++) {
        m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
        if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
        m[i] -= ncl;
    }
    return m;
}

int **imatrix(int nrl, int nrh, int ncl, int nch)
{
    int i,**m;

    m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
    if (!m) nrerror("allocation failure 1 in imatrix()");
    m -= nrl;

    for(i=nrl;i<=nrh;i++) {
        m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
        if (!m[i]) nrerror("allocation failure 2 in imatrix()");
        m[i] -= ncl;
    }
    return m;
}

double ***dmatrix3(int nrl, int nrh, int ncl, int nch, int n3l, int n3h)
{
    int i, j;
    double ***m;

    m=(double ***) malloc((unsigned) (nrh-nrl+1)*sizeof(double**));
    if (!m) nrerror("error in dmatrix3");
    m -= nrl;
    for(i=nrl; i<=nrh;i++) {
        m[i]=(double **) malloc((unsigned) (nch-ncl+1)*sizeof(double*));
        if (!m[i]) nrerror("error in dmatrix3");
        m[i] -= ncl;
        for(j=ncl; j<=nch; j++){
            m[i][j]=(double *) malloc((unsigned) (n3h-n3l+1)*sizeof(double));
            if (!m[i][j]) nrerror("error in dmatrix3");
	        m[i][j] -= n3l;
		}      
	}
    return m;
}

void free_dmatrix3(double ***m, int nrl, int nrh, int ncl, int nch, int n3l, int n3h) 
{
	int i,j;
	for(i=nrh;i>=nrl;i--)
   	   for(j=nch;j>=ncl;j--)
		free((char*) (m[i][j]+n3l));
	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_ivector(int *v, int nl, int nh)
{
    free((char*) (v+nl));
}

void free_cvector(char *v, int nl, int nh) 
{ 
	free((char*) (v+nl));
}

void free_dvector(double *v, int nl, int nh)
{
    free((char*) (v+nl));
}

void free_cmatrix(char **m, int nrl, int nrh, int ncl, int nch)
{
    int i;

    for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
    free((char*) (m+nrl));
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
    int i;

    for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
    free((char*) (m+nrl));
}

void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
{
    int i;

    for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
    free((char*) (m+nrl));
}
            
/////////////////////////////////////////////////////////////////////////
//obtain p-value

#define NR_END 1
#define FREE_ARG char*

#include <cmath>
#include <limits>
#include <stddef.h>

using namespace std;

double gammq(const double a, const double x)
{
	double gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0)
		nrerror("Invalid arguments in routine gammq");
	if (x < a+1.0) {
		gser(gamser,a,x,gln);
		return 1.0-gamser;
	} else {
		gcf(gammcf,a,x,gln);
		return gammcf;
	}
}


void gser(double &gamser, const double a, const double x, double &gln)
{
	const int ITMAX=1000;
	const double EPS=numeric_limits<double>::epsilon();
	int n;
	double sum,del,ap;

	gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) 
			nrerror("x less than 0 in routine gser");
		   gamser=0.0;
		  return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=0;n<ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				gamser=sum*exp(-x+a*log(x)-gln);
				return;
			}
		}
	    nrerror("a too large, ITMAX too small in routine gser");

		return;
	}
}



void gcf(double &gammcf, const double a, const double x, double &gln)
{
	const int ITMAX=1000;
	const double EPS=numeric_limits<double>::epsilon();
	const double FPMIN=numeric_limits<double>::min()/EPS;
	int i;
	double an,b,c,d,del,h;

	gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) <= EPS) break;
	}
	if (i > ITMAX) 
		nrerror("a too large, ITMAX too small in gcf");
	gammcf=exp(-x+a*log(x)-gln)*h;
}



double gammln(const double xx)
{
	int j;
	double x,y,tmp,ser;
	static const double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
		-0.5395239384953e-5};

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<6;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}


////////////////////////////////////////////////////////////////////////
//inverse

void ludcmp(double **a,int n,int *indx) 
{
    int i,imax=0,j,k;
    double big,dum,sum,temp;
    double *vv, TINY;

        TINY=1.0e-40;
    vv=dvector(1,n);
    for (i=1;i<=n;i++) {
        big=0.0;
        for (j=1;j<=n;j++)
            if ((temp=fabs(a[i][j])) > big) big=temp;
        if (big == 0.0) nrerror("Sigmatrix in ludcmp");
        vv[i]=1.0/big;
    }
    for (j=1;j<=n;j++) {
        for (i=1;i<j;i++) {
            sum=a[i][j];
            for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
        }
        big=0.0;
        for (i=j;i<=n;i++) {
            sum=a[i][j];
            for (k=1;k<j;k++)
                sum -= a[i][k]*a[k][j];
            a[i][j]=sum;
            if ( (dum=vv[i]*fabs(sum)) >= big) {
                big=dum;
                imax=i;
            }
        }
        if (j != imax) {
            for (k=1;k<=n;k++) {
                dum=a[imax][k];
                a[imax][k]=a[j][k];
                a[j][k]=dum;
            }
            vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (a[j][j] == 0.0) a[j][j]=TINY;
        if (j != n) {
            dum=1.0/(a[j][j]);
            for (i=j+1;i<=n;i++) a[i][j] *= dum;
        }
    }
    free_dvector(vv,1,n);
}
 
void lubksb(double **a, int n, int *indx, double *b)
{
    int i,ii=0,ip,j;
    double sum;

    for (i=1;i<=n;i++) {
        ip=indx[i];
        sum=b[ip];
        b[ip]=b[i];
        if (ii)
            for (j=ii;j<=i-1;j++)
                     sum -= a[i][j]*b[j];
        else if (sum) ii=i;
        b[i]=sum;
    }
    for (i=n;i>=1;i--) {
        sum=b[i];
        for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
        b[i]=sum/a[i][i];
    }
}
void invv(double **fma, int nn)  /* any nn x nn matrix's inverse */ 
{
    int  i, j, k, *indx;
    double *COL, **ipi;

    indx=ivector(1, nn);    COL=dvector(1, nn);     ipi=dmatrix(1, nn, 1, nn);

    ludcmp(fma, nn, indx);  
    for(i=1; i<=nn; i++){
        for(j=1; j<=nn; j++) COL[j]=0.0;
        COL[i]=1.0;     
        lubksb(fma, nn, indx, COL);         
        for(j=1; j<=nn; j++) ipi[j][i]=COL[j];
    }
    for(k=1; k<=nn; k++)
        for(i=1; i<=nn; i++) fma[k][i]=ipi[k][i];

    free_ivector(indx, 1, nn);  
    free_dvector(COL, 1, nn);   
    free_dmatrix(ipi, 1, nn, 1, nn);
 
}   // end of invv.c // 

//////////////////////////////////////////////////////////////////////////////////////
//Link list


void makelist_unlinklist_ASEnull(int**hap1, int**hap2, int*nhap, int*rm, int n){
        
    for (int i=0;i<n;i++){

		rm[i] = 0;

		hap1[i][0] = 0;
		hap2[i][0] = 0;
 		nhap[i] = 1;       
   }
}

void makelist_unlinklist(int**hap1, int**hap2, int*nhap, int*rm, int*Gu, int**H, int n, int M, int*pow2){
    
    int pow1 = pow2[M-1];
    
    for (int i=0;i<n;i++){
        
		rm[i] = 0;

		if (Gu[i]==1){      //heterozygous at target SNP
			hap1[i][0] = H[i][0];		hap1[i][1] = H[i][0]+pow1;
			hap2[i][0] = H[i][1]+pow1;	hap2[i][1] = H[i][1];
			nhap[i] = 2;
		}
		else if (Gu[i]==0) {//homozygous at target SNP
			hap1[i][0] = H[i][0];
			hap2[i][0] = H[i][1];
			nhap[i] = 1;
		}
		else if (Gu[i]==2){
			hap1[i][0] = H[i][0]+pow1;
			hap2[i][0] = H[i][1]+pow1;
			nhap[i] = 1;
		}
    }
}//makelist_unlinklist


void makelist_unlinklist_Gp(int**hap1, int**hap2, int*nhap, int*rm, int n, int**Gp){
    
    
    for (int i=0;i<n;i++){
        
		rm[i] = 0;

		hap1[i][0] = Gp[i][0];
		hap2[i][0] = Gp[i][1];
		nhap[i] = 1;
    }
}//makelist_unlinklist


void reduce_list_unlinklist(int**hap1, int**hap2, int*nhap, int*rm,double thresh,double *pai, int n){
	
	for (int i = 0; i < n; i++) {

		rm[i] = 0;
		
		if (nhap[i] == 2) {
			int rm_node_0 = (pai[hap1[i][0]] < thresh || pai[hap2[i][0]] < thresh);
			int rm_node_1 = (pai[hap1[i][1]] < thresh || pai[hap2[i][1]] < thresh);

			if (rm_node_0 && !rm_node_1) {
				hap1[i][0] = hap1[i][1];
				hap2[i][0] = hap2[i][1];
				nhap[i] = 1;
			}
			else if (!rm_node_0 && rm_node_1) {
				nhap[i] = 1;
			}
			else if (rm_node_0 && rm_node_1) {
				nhap[i] = 0;
				rm[i] = 1;
			}
		}
		else if (nhap[i] == 1) {
			int rm_node_0 = (pai[hap1[i][0]] < thresh || pai[hap2[i][0]] < thresh);
			if (rm_node_0) {
				nhap[i] = 0;
				rm[i] = 1;
			}
		}
		else if (nhap[i] == 0) {
			rm[i] = 1;
		}
	}
}//reduce_list_unlinklist


//end of Link List
//////////////////////////////////////////////////////////////////////////////////////

inline bool my_isnan(double x){
   return x != x;
}

double length (double*a,int q){
    double s=0;
    int i;
    for (i=0;i<q;i++)  s+=fabs(a[i]);
    return s;
}


void reduce_hap(int *hap_left,double *pai,double threshold_F, int K){
    int i=0,j;
    for (j=0;j<K;j++){
        if (pai[j]>=threshold_F){
            hap_left[i]=j;
            i++;
		}
	}
}

int explosion_E(double**E, int p) {

	int i, j;
    double small, t, sum;

    for (i = 1; i <= p; i++) {
        small = 1.0/EPS0; 
        sum   = 0.0;
        for (j = 1; j <= p; j++) {
            sum += fabs(E[i][j]);
            if (( t = fabs(E[i][j]) ) < small) {
                small = t;
			}
		}
        if (small >= 1.0/EPS0 || my_isnan(sum) || sum < EPS0*EPS0)  return 1;
	}

    return 0;
}// explosion_E

int inproper_std(double**E, int p) {
     
	int i;
    double thresh = 1.0/EPS0;
	
	for (i = 1; i <= p; i++) {
        if (E[i][i] < 0.0 || E[i][i] > thresh) return 1;
	}

    return 0;
}// inproper_std


int A1_unlinklist(double*pai, int**hap1, int**hap2, double**prob, int*nhap, int*rm, 
	   int Y1, int Y2, int Y0, int *STATUS, int design, int K, int n) {
	//design: 1, case-control, Y0 arbitrary
	//exter: -1: no external; 1: external

	//Y1   Y2   Y0   exter
	//------------------
	//0    0    1   -1: control
	//0    1    1   -1: study
	//-1  -1   -1   1: external
	//0    1    1   1: study + external

	//design: 2, cohort, 3, cross-sectional, Y1, Y2 arbitrary
	//int*STATUS: artifical, defined outside
	//Y1   Y2   Y0   exter
	//-------------------
	//-1  -1     1   -1: study
	//-1  -1    -1   1: external
	//-1  -1     1   1: study + exernal

    int l, i, j, nST, iter;
	double loglike;
     
    double *diff = dvector(0, K-1);
    double *s    = dvector(0, K-1);

    //initialization
    for (j = 0; j < K; j++) diff[j] = 1; 

	l = 1;
	iter = 0;
	
    //EM algorithm
    while (length(diff, K) > TOL6 && iter < ITMAX0){
        iter++; 

        //E-step: evaluate pikl

		loglike = 0.0;

		nST = 0;
		if (Y0 == 1) {    
			for (i = 0; i < n; i++) {		//study
				if (rm[i]) continue;

				if ((design == 1) ? (STATUS[i] == Y1 || STATUS[i] == Y2) : 1) {
					nST++;
					double denm = 0;
					for (j = 0; j < nhap[i]; j++) {
						denm   += pai[hap1[i][j]]*pai[hap2[i][j]];
					}
					if (denm == 0) denm = 1;
					for (j = 0; j < nhap[i]; j++) {
						prob[i][j] = pai[hap1[i][j]]*pai[hap2[i][j]]/denm;
					}
					loglike += log(denm);
				}
				

			}

		}
    
        //M-step: evaluate pai_com(m+1)

        for (j = 0; j < K; j++) s[j] = 0;

		if (Y0 == 1) {      
			
			for (i = 0; i < n; i++) {										//study
				if (rm[i]) continue;

				if ((design == 1) ? (STATUS[i] == Y1 || STATUS[i] == Y2) : 1) {
					for (j = 0; j < nhap[i]; j++) {
						s[hap1[i][j]]+=prob[i][j] * 0.5;
						s[hap2[i][j]]+=prob[i][j] * 0.5;
					}
				}
			}
		}


        //difference
        for (j = 0; j < K; j++) {
            s[j]    = s[j] / nST;
            diff[j] = s[j] - pai[j];
            pai[j]  = s[j];
        }
    }//while

	if (iter == ITMAX0) l = 0;

    free_dvector(diff, 0, K-1); 
    free_dvector(s,    0, K-1);

	return l;

}//A1_unlinklist


void get_hap_freq_unlinklist(PARA*p, DATA*d) {

	int i, j;

	for (i = 0; i < p->K; i++) p->pai[i] = 1.0/p->K;
	int *STATUS=NULL, design=3;
	A1_unlinklist(p->pai, d->hap1, d->hap2, d->prob, d->nhap, d->rm, -1, -1, 1, STATUS, design, p->K, d->n); //study + external

	reduce_list_unlinklist(d->hap1, d->hap2, d->nhap, d->rm, p->THRESHOLD_L, p->pai, d->n);

	
	//if (d->ext) d->list_ex = reduce_list_ex(d->list_ex, p->THRESHOLD_L, p->pai);

	//-------------- KK, hap_freq, h -------------------//

	//1). replace p->pai by pai0 2). add comments
	p->KK = 0; 
    for (i = 0; i < p->K; i++) {
        if (p->pai[i] >= p->THRESHOLD_L) p->hap_freq[p->KK++] = i; 
	}

	for (i = 0; i < p->KK; i++) {     //order hap_freq
        int max_pos = i;
		for (j = i+1; j < p->KK; j++) {
            if (p->pai[p->hap_freq[j]] > p->pai[p->hap_freq[max_pos]]) max_pos = j;
		}
		int tmp = p->hap_freq[i];
		p->hap_freq[i] = p->hap_freq[max_pos];
        p->hap_freq[max_pos] = tmp;
	}

	for (i = 0; i < p->KK; i++) p->h[p->hap_freq[i]] = i;         //h[.]

	//recalculate p->pai
	A1_unlinklist(p->pai, d->hap1, d->hap2, d->prob, d->nhap, d->rm, -1, -1, 1, STATUS, design, p->K, d->n); //study + external

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////                     core programs                      ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void get_initial_para(int*ll_gene, double*vTREC_gene, double*vASE_gene, DATA*d, PARA*p, int TREC_fitted_gene, int ASE_fitted_gene) {

	int i;

	////////////////////////////////////////////////
	// get initial values: 
	// ll_gene, p->disp_TREC, vASE_gene, vTREC_gene
	// with a simplified link list
	////////////////////////////////////////////////
	
    makelist_unlinklist_ASEnull(d->hap1, d->hap2, d->nhap, d->rm, d->n);
	//d->list_ex = NULL;  //list_ex does not contribute to estimation of initial para

	p->M = 1;
	p->K = 2;
	get_hap_freq_unlinklist(p, d);

	double muT = 0.0, varT = 0.0;
	for (i = 0; i < d->n; i++) {
		muT += d->T[i];
		varT+= 1.0*d->T[i]*d->T[i];
	}
	muT /= d->n;
	varT = varT/d->n - muT*muT;


	double phi_GP = 1.0 - sqrt(muT/varT);
	double mu_GP = muT*(1.0-phi_GP);
	
	if (TREC_fitted_gene) {

		double *stdTREC_gene = new double[1+p->R];
		double loglike_null_TREC;

		int disp_TREC = 0; //no dispersion

		//--------------------------------------//
		//----------- fit Poission  ------------//
		//--------------------------------------//

		//initialization
		vTREC_gene[0] = log(muT);   //intercept
		vTREC_gene[1] = 0.1;        //log(cntAll)
		for (int j = 2; j < p->R; j++) vTREC_gene[j] = 0.0; 

		ll_gene[0] = fit_TReC_null(vTREC_gene, stdTREC_gene, &loglike_null_TREC, *d, p->R, disp_TREC, p->disp_TREC_model);

        if (ll_gene[0] == 4)		printf("getting initial para for Poisson model failed!: E exploded during iteration!\n");
		else if (ll_gene[0] == 3)	printf("getting initial para for Poisson model failed!: halving algorithm cannot increase log-likelihood!\n");
		else if (ll_gene[0] == 2)	printf("getting initial para for Poisson model failed!: Maximum iteration is reached!\n");

		//--------------------------------------//
		//---------- test dispersion  ----------//
		// significant when: 1). no dispersion
		//					 2). SNP effect exist
		//--------------------------------------//

		if (varT-muT > 0) { 
			const double SCORE_THRESH = 2.576; //5%: 1.96; 1%: 2.576
			double score_disp_TREC;

			if (strcmp(p->disp_TREC_model, "NB")==0) {
				get_score_disp_TREC(&score_disp_TREC, vTREC_gene, *d, p->R);
			}
			else if (strcmp(p->disp_TREC_model, "GP")==0) {
				score_disp_TREC = 3; //???
			}
		
			if (fabs(score_disp_TREC) >= SCORE_THRESH) { //dispersion exist!

				//printf("dispersion (Poisson) exists!\n");

				disp_TREC = 1;

				//--------------------------------------//
				//-------------- fix NB  ---------------//
				//--------------------------------------//

				//initialization
				for (int j = p->R-1; j >= 0; j--) {
					vTREC_gene[j+1] = vTREC_gene[j];
				}
				if (strcmp(p->disp_TREC_model, "NB")==0) {
					vTREC_gene[0] = 1.0/((varT-muT)/muT/muT);//transformed phi(TReC): 1/phi
				}
				else if (strcmp(p->disp_TREC_model, "GP")==0) {
					vTREC_gene[0] = phi_GP;
					vTREC_gene[1] = log(mu_GP);
				}


				ll_gene[0] = fit_TReC_null(vTREC_gene, stdTREC_gene, &loglike_null_TREC, *d, p->R, disp_TREC, p->disp_TREC_model);


				if (ll_gene[0] == 4)		printf("getting initial para for NB model failed!: E exploded during iteration!\n");
				else if (ll_gene[0] == 3)	printf("getting initial para for NB model failed!: halving algorithm cannot increase log-likelihood!\n");
				else if (ll_gene[0] == 2)	printf("getting initial para for NB model failed!: Maximum iteration is reached!\n");
			}
		}//varT-muT > 0

		p->disp_TREC = disp_TREC;

		delete [] stdTREC_gene;

	}

	if (ASE_fitted_gene) {

		double*stdASE_gene = new double[1]; 
		double loglike_null_ASE;

		//initialization: grid search
		vASE_gene[0] = 0.00;
		ll_gene[1] = fit_ASE_null(vASE_gene, stdASE_gene, &loglike_null_ASE, *d, *p);

		while (vASE_gene[0] < 0.5 && ll_gene[1] != 0) {
			vASE_gene[0] += 0.025;
			ll_gene[1] = fit_ASE_null(vASE_gene, stdASE_gene, &loglike_null_ASE, *d, *p);
		}

		if (ll_gene[1] != 0) {
			vASE_gene[0] = 0.00;
			while (vASE_gene[0] > -0.2 && ll_gene[1] != 0) {
				vASE_gene[0] -= 0.025;
				ll_gene[1] = fit_ASE_null(vASE_gene, stdASE_gene, &loglike_null_ASE, *d, *p);
			}
		}

		if (ll_gene[1] == 5)		printf("getting initial para for ASE model failed!: New v[AD] is below boundary!\n");
        else if (ll_gene[1] == 4)	printf("getting initial para for ASE model failed!: E exploded during iteration!\n");
		else if (ll_gene[1] == 3)	printf("getting initial para for ASE model failed!: halving algorithm cannot increase log-likelihood!\n");
		else if (ll_gene[1] == 2)	printf("getting initial para for ASE model failed!: Maximum iteration is reached!\n");

		delete [] stdASE_gene;
	}

}//get_initial_para

void get_null_para_loglike(RESULTS*r, DATA*d, PARA*p) {

	///////////////////////////////////////////////////////////////
	//  fit null model: with final link list
	//  get loglike_null[3]
	///////////////////////////////////////////////////////////////

	if (p->TREC_fitted) {
		p->TREC = 1;
		p->ASE = 0;
		get_index(p->ASE, p->TREC, p);

		r->ll_null[0] = fit_TReC_null(r->vTREC, r->std_vec, &(r->loglike_null[0]), *d, p->R, p->disp_TREC, p->disp_TREC_model);

		if		(r->ll_null[0] == 4) printf("fitting null NB model failed!: E exploded during iteration!\n");
		else if (r->ll_null[0] == 3) printf("fitting null NB model failed!: halving algorithm cannot increase log-likelihood!\n");
		else if (r->ll_null[0] == 2) printf("fitting null NB model failed!: Maximum iteration is reached!\n");
		else if (r->ll_null[0] == 1) printf("fitting null NB model failed!: EM algorithm converged but Louis estimator is not positive-definite!\n");
	}

	if (p->ASE_fitted) {
		p->TREC = 0;
		p->ASE = 1;
		get_index(p->ASE, p->TREC, p);

		r->ll_null[1] = fit_ASE_null(r->vASE, r->std_vec, &(r->loglike_null[1]), *d, *p);

		if		(r->ll_null[1] == 4) printf("fitting null ASE model failed!: E exploded during iteration!\n");
		else if (r->ll_null[1] == 3) printf("fitting null ASE model failed!: halving algorithm cannot increase log-likelihood!\n");
		else if (r->ll_null[1] == 2) printf("fitting null ASE model failed!: Maximum iteration is reached!\n");
		else if (r->ll_null[1] == 1) printf("fitting null ASE model failed!: EM algorithm converged but Louis estimator is not positive-definite!\n");
	}

	if (p->ASE_fitted & p->TREC_fitted) {

		r->ll_null[2] = r->ll_null[0] + r->ll_null[1];
		r->loglike_null[2] = r->loglike_null[0] + r->loglike_null[1];
	}

}//get_null_para_loglike



void Rvec(double*vecR, int R, double*X) {
	vecR[0] = 1.0;
	for (int i = 0; i < R-1; i++) {
		vecR[1+i] = X[i];
	}
}

void Bvec(double*vecB, int B, double*X) {
	if (B == 1) {//MAIN
		vecB[0] = 1.0;
	}
	else if (B == 2) {//INT, MAIN  //need to be revised for multivariate X
		vecB[0] = X[0];  
		vecB[1] = 1.0;  
	}
}

void get_index(int ASE, int TREC, PARA*p) {
	//input : p->B, p->R, p->disp_TREC, p->KK
	//output: p->AD, p->TD, p->REG, p->FREQ, p->nfix, p->nall

	p->ASE = ASE;
	p->TREC = TREC;

	if (p->ASE & !p->TREC) {
		p->AD = p->B;      //2					 
		p->FREQ = p->AD+1; //3					 
		p->nfix = p->FREQ; //3
		p->nall = p->nfix + p->KK;
	}
	else if (p->ASE & p->TREC) {
		p->AD = p->B;				//2	//p->AD starting point of ASE dispersion para
		p->TD = p->AD+1;			//3	//p->TD starting point of TReC dispersion para
		p->REG = p->TD+p->disp_TREC;//4	//p->REG starting point of other reg of TReC excluding SNPs
		p->FREQ = p->REG+p->R;		//13//p->FREQ starting point of hap frequency
		p->nfix = p->FREQ;			//13//p->nfix fixed para excluding hap frequency
		p->nall = p->nfix + p->KK;
	}	
	else if (!p->ASE & p->TREC) {
		p->TD = p->B;			//2
		p->REG = p->TD+p->disp_TREC;//3
		p->nfix = p->REG+p->R;	//12
		p->nall = p->nfix;
	}
}

double MDfun(double *pai, int M, int K){

	//////////////////////////////////////// 
	// new "MD" measure for this manuscript
	//////////////////////////////////////// 

	int k, l;

	int pow2 = K;

	double p1 = 0.0;
	//double EX = 0.0;
    double EX2 = 0.0;
	
	//faster version
	for (k = 0; k < K; k++) {

		double s = 2*pai[k+pow2]*pai[k];

		if (s > TOL6) {
			p1 += s;
			//EX += pai[k+pow2]*pai[k];
			EX2+= pai[k+pow2]*pai[k]*pai[k+pow2]*pai[k]/s;
		}
	}

	for (k = 0; k < K; k++) {
		for (l = k+1; l < K; l++) {

			double s = pai[k+pow2]*pai[l] + pai[k]*pai[l+pow2];

			if (s > TOL6) {
				p1 += 2*s;
				//EX += pai[k+pow2]*pai[l]+pai[l+pow2]*pai[k];
				EX2+= (pai[k+pow2]*pai[l]*pai[k+pow2]*pai[l]+pai[l+pow2]*pai[k]*pai[l+pow2]*pai[k])/s;
			}
		}
	}


	//EX  /= p1;
	EX2 /= p1;

	double MD = (EX2 - 0.25)/0.25;

	return MD;

}//MDfun


double corr2(int*Y, int*X, int n) {
	int n_true = 0;
	double XY = 0.0, X2 = 0.0, Y2 = 0.0, Xmean = 0.0, Ymean = 0.0;
	for (int i = 0; i < n; i++) {
		if (Y[i] != 9 && X[i] != 9) {
			int Xr = (X[i]/3) + (X[i]%3);
			int Yr = (Y[i]/3) + (Y[i]%3);
			XY += Xr * Yr;
			X2 += Xr * Xr;
			Y2 += Yr * Yr;
			Xmean += Xr;
			Ymean += Yr;
			n_true++;
		}
	}
	Xmean /= n_true;
	Ymean /= n_true;

	return (XY-n_true*Xmean*Ymean)*(XY-n_true*Xmean*Ymean)/(X2-n_true*Xmean*Xmean)/(Y2-n_true*Ymean*Ymean);
}//corr


int screen_SNP_LD(int*map_rf, int*pwin_r, int nsubj, int tranSNPStart, int tranSNPEnd, int**haplo, double R2_thresh){
	

	const double MISS_THRESH = 0.1;

	int pwin = tranSNPEnd - tranSNPStart + 1;

	int i, j;

	//----------------------------------------------------------
	//----- delete high-LD, high-missing, monomorphic SNPs -----
	//----------------------------------------------------------
	// pwin_r
	// SNP_out[pwin], map_rf[pwin_r]
	//----------------------------------------------------------

	int*SNP_out = new int[pwin]; 	
	for (j = 0; j < pwin; j++) SNP_out[j] = 0;

	*pwin_r = 0;
	for (j = 0; j < pwin; j++) {
		
		double miss_j = 0.0, freq_j = 0.0;
		for (i = 0; i < nsubj; i++) {
			int haplotype = haplo[tranSNPStart+j][i];
			int genotype  = (haplotype>2) ? (haplotype-2) : haplotype;
			if (genotype == 9) miss_j += 1.0;
			else freq_j += genotype;		
		}
		miss_j /= nsubj;
		freq_j /= 2*nsubj;
		if (miss_j > MISS_THRESH || freq_j < TOL6 || freq_j > 1.0-TOL6) { //delete high-missing, monomorphic
			SNP_out[j] = 1; 
			continue;
		}           

		for (int k = MAX(0,j-10); k < j; k++) {
			
			if (SNP_out[k] == 1) continue;

			double rho2 = corr2(haplo[tranSNPStart+j], haplo[tranSNPStart+k], nsubj);
			if (rho2 > R2_thresh) { //delete high-LD
				SNP_out[j] = 1; 
				break;
			}				

		}//k

		if (SNP_out[j] == 0) map_rf[(*pwin_r)++] = tranSNPStart+j;

	}//j
	

	delete [] SNP_out;
	
	return 0;

}//screen_SNP_LD


int get_tagSNP(int*tagSNP, int*M, double*MDmax, int iSNP, int nsubj, int*map_rf, int pwin_r, int**haplo, int Mmax, int*pow2){
	
	int i, j;

	int**H = imatrix(0, nsubj-1, 0, Mmax-1);
	for (i = 0; i < nsubj; i++) H[i][0] = haplo[iSNP][i];

	
	if (pwin_r+1 < Mmax) {

		*M = pwin_r+1;

		for (int iM = 0; iM < pwin_r; iM++) tagSNP[iM] = map_rf[iM];

		for (i = 0; i < nsubj; i++) {
			for (int k = 0; k < pwin_r; k++) {
				H[i][k+1] = haplo[map_rf[k]][i];
			}
		}

		get_MD(MDmax, H, nsubj, *M, pow2); 

		return 0;
	}


	*M = Mmax;

	//---------------------------------------------	
	//---- get marginal R2 for individual SNPs ----
	//---------------------------------------------
	// R2_measure[pwin_r] 
	//---------------------------------------------
 
	double*R2_measure = new double[pwin_r]; 
	
	int*tarSNP  = new int[nsubj];
	for (i = 0; i < nsubj; i++) tarSNP[i] = haplo[iSNP][i];

	int*tranSNP = new int[nsubj];
	for (j = 0; j < pwin_r; j++) {

		for (i = 0; i < nsubj; i++) tranSNP[i] = haplo[map_rf[j]][i];

		double rho2 = corr2(tarSNP, tranSNP, nsubj);
		R2_measure[j] = rho2;
	}
	delete [] tranSNP;

	//---------------------------------------------	
	//----------- sort & get top-R2 SNPs ----------
	//---------------------------------------------	
	// index_sort[pwin_r]
	//---------------------------------------------	
	
	int*index_sort = new int[pwin_r];
	for (j = 0; j < pwin_r; j++) index_sort[j] = j;

	for (j = 0; j < pwin_r; j++) {
		int max_pos = j;
		for (int k = j+1; k < pwin_r; k++) {
			if (R2_measure[k] > R2_measure[max_pos]) {
				max_pos = k;
			}
		}
		double R2_tmp = R2_measure[j];
		R2_measure[j] = R2_measure[max_pos];
		R2_measure[max_pos] = R2_tmp;
		int index_tmp = index_sort[j];
		index_sort[j] = index_sort[max_pos];
		index_sort[max_pos] = index_tmp;
	}

	//----------------------------------------------------------------------
	//----------------- exhaustive search of haplotypes --------------------
	//----------------------------------------------------------------------
	// max_tag[p.M]
	//----------------------------------------------------------------------

	const int CANDmax = 30;

	int pwin_search = (pwin_r > CANDmax) ? CANDmax : pwin_r;

	int*r = new int[*M-1];
	int*rmax = new int[*M-1];

	*MDmax = -1/TOL6;
	double MD;

	for (r[0] = 0; r[0] <= pwin_search-4; r[0]++) {
		for (r[1] = r[0]+1; r[1] <= pwin_search-3; r[1]++) {
			for (r[2] = r[1]+1; r[2] <= pwin_search-2; r[2]++) {
				for (r[3] = r[2]+1; r[3] <= pwin_search-1; r[3]++) {

					for (i = 0; i < nsubj; i++) {
						for (int k = 0; k < *M-1; k++) {
							H[i][k+1] = haplo[map_rf[index_sort[r[k]]]][i];
						}
					}

				 	get_MD(&MD, H, nsubj, *M, pow2); 

					if (MD > *MDmax + TOL6) {
						*MDmax = MD;
						for (int k = 0; k < *M-1; k++) {
							rmax[k] = r[k];
							tagSNP[k] = map_rf[index_sort[rmax[k]]];
						}
						
					}
				}
			}
		}
	}//go over


	delete [] r;
	delete [] rmax;
	free_imatrix(H, 0, nsubj-1, 0, Mmax-1);

	//--------------------- shutting down ------------------------//

	
	delete [] tarSNP;
	delete [] R2_measure; 
	delete [] index_sort;

	return 0;

}//get_tagSNP

int get_MD(double*MD, int**H, int nsubj, int M, int*pow2) {

	int j;

	int K = pow2[M-1] * 2;
	double*pai = new double[K];
	for (j = 0; j < K; j++) pai[j] = 0.0;

	double s = 0.5/nsubj;
	for (int i = 0; i < nsubj; i++) {

		int phase1 = 0, phase2 = 0; 

		for (j = 0; j < M; j++) {
			phase1 += (H[i][j]/3) * pow2[M-1-j];
			phase2 += (H[i][j]%3) * pow2[M-1-j];
		}
		pai[phase1] += s;
		pai[phase2] += s;
	}

	//double denm = 0.0;
	//for (j = 0; j < K; j++) denm   += pai[j];
	//for (j = 0; j < K; j++) pai[j] /= denm;

	*MD = MDfun(pai, M, pow2[M-1]); 

	delete [] pai;

	return 0;

}//get_MD



int fit_ASE_null(double*v0, double*std0, double*loglike, DATA d, PARA p) {

	///////////////////////////////// 
	// only estimate dispersion para
	///////////////////////////////// 
    
	int AD = 0;
	int nfix = 1;
	p.FREQ = 1; 
	int nall = nfix + p.KK;

	int p0, i, j, k, r, ll = 0, iter = 0;

	double*diff  = new double[nfix];
	double*v     = new double[nall];
	double*vnew  = new double[nfix];
    double*vold  = new double[nfix];
	double*U     = new double[nfix];
	double**E    = dmatrix(1, nfix, 1, nfix);
    
	//initialization

	v[AD] = v0[AD];		//theta (ASE)
	for (j = 0; j < p.KK; j++) v[p.FREQ+j] = p.pai[p.hap_freq[j]];	//pi: hap freq

    //ASE variables
	int max_M, max_N;
    double pn, theta, boundary;
	double *r0, *p0_r0, *l_r0, *l_r0_r, *p0_r0_r, *l_r0_r_2, *p0_r0_r_2;
    
	max_M = 0;
	max_N = 0;
	for (p0 = 0; p0 < d.n; p0++) {
		if (d.rm[p0]) continue;
		if (d.M[p0] > max_M)         max_M = d.M[p0];
		if (d.N[p0]-d.M[p0] > max_M) max_M = d.N[p0]-d.M[p0];
		if (d.N[p0] > max_N)         max_N = d.N[p0];
	}    

	r0          = new double[max_N+1];
	//beta-binomial distribution
	p0_r0       = new double[max_M+1]; //sum_r log(p0+r0)
	l_r0        = new double[max_N+1]; //sum_r log(1+r0)    
	//theta
	l_r0_r      = new double[max_N+1]; //\sum_r r/(1+r0)
	p0_r0_r     = new double[max_M+1]; //\sum_r r/(p0+r0)
	//theta2
	l_r0_r_2    = new double[max_N+1]; //\sum_r [r/(1+r0)]^2
	p0_r0_r_2   = new double[max_M+1]; //\sum_r [r/(p0+r0)]^2       

	//--------------- EM-algorithm: within M-step, one-step Newton-Raphson --------------//

	double loglike_old = 1;
	*loglike = 2;
	for (j = 0; j < nfix; j++) U[j] = 0.0; 
    U[0] = 1.0;

	while ((fabs(loglike_old-*loglike)/fabs(loglike_old) > TOL6 || length(U, nfix) > 0.001) && iter < ITMAX0) {
	//while (length(diff, nfix) > TOL6 && iter < ITMAX0) {

		for (j = 0; j < nfix; j++) vold[j] = v[j];
		loglike_old = *loglike;
        
        //save variables

		pn = 0.5;
		theta = v[AD];

		boundary = -100; //MAX(-pn/(max_M-1),-1.0/(max_N-1));
		
		for (r = 1; r <= max_N; r++) r0[r] = (r-1)*theta;

		//beta-binomial distribution
		p0_r0[0] = 0;
		l_r0[0] = 0;
        
		for (r = 1; r <= max_M; r++) {
			p0_r0[r]    = p0_r0[r-1]  + log(  pn + r0[r]);
		}
		for (r = 1; r <= max_N; r++) {
			l_r0[r]   = l_r0[r-1] + log(1 + r0[r]);
		}
        
		//----------------------- likelihood -----------------------//

		*loglike = 0.0;

		for (p0 = 0; p0 < d.n; p0++) {
			if (d.rm[p0]) continue;
             
            d.sum[p0] = 0.0;

			int index1 = d.M[p0];
			int index2 = d.N[p0] - d.M[p0];
			int index3 = d.N[p0];

			for (int pp = 0; pp < d.nhap[p0]; pp++) {
				d.sum[p0] += p.pai[d.hap1[p0][pp]]*p.pai[d.hap2[p0][pp]];
			}

            //likelihood
			*loglike -= p0_r0[index1] + p0_r0[index2] - l_r0[index3] + log(d.sum[p0]);

			if (index3 > 1) boundary = MAX(boundary, -1.0/(index3-1));
			if (index1 > 1) boundary = MAX(boundary, -pn/(index1-1));
			if (index2 > 1) boundary = MAX(boundary, -pn/(index2-1));
            
		}//person


		//if (p.ASE) {
		//	if (d.ext) link_orig_ext(&d, &p, loglike, v);
		//}

		
		//----------------------- Newton Raphson -----------------------//
        
		//theta
		l_r0_r[0] = 0;
		p0_r0_r[0] = 0;
		//theta2
		l_r0_r_2[0] = 0;
		p0_r0_r_2[0] = 0;
        
		for (r = 1; r <= max_M; r++) {
			p0_r0_r[r]      = p0_r0_r[r-1]  + (r-1)/  (pn + r0[r]);   
			p0_r0_r_2[r]    = p0_r0_r_2[r-1]  + (r-1)*(r-1)/  ((pn+r0[r])*  (pn+r0[r]));
		}      
		for (r = 1; r <= max_N; r++) {
			l_r0_r[r]     = l_r0_r[r-1] + (r-1)/(1+r0[r]);
			l_r0_r_2[r]   = l_r0_r_2[r-1] + (r-1)*(r-1)/((1+r0[r])*(1+r0[r]));
		}

        //U, E initialation
        //U = l', E = -l"
		for (i = 0; i < nfix; i++) {
			U[i] = 0;
			for (j = i; j < nfix; j++) E[i+1][j+1] = 0;
		}
		
		for (p0 = 0; p0 < d.n; p0++) {
			if (d.rm[p0]) continue;
            
			//default
			int index1 = d.M[p0];
			int index2 = d.N[p0] - d.M[p0];
			int index3 = d.N[p0];
            
			U[AD]         += (p0_r0_r[index1] + p0_r0_r[index2] - l_r0_r[index3]);
			E[AD+1][AD+1] -= (-p0_r0_r_2[index1] - p0_r0_r_2[index2] + l_r0_r_2[index3]);

		}//person


		for (k = 0; k < nfix; k++) {
            for (j = 0; j < k; j++) E[k+1][j+1] = E[j+1][k+1];
		}

        if ((ll = explosion_E(E, nfix))) {ll = 4; goto end;}
        invv(E, nfix);
        if ((ll = explosion_E(E, nfix))) {ll = 4; goto end;}

        for (k = 0; k < nfix; k++) {
			diff[k] = 0.0;
            for (j = 0; j < nfix; j++) {
				diff[k] += E[k+1][j+1] * U[j];
			}
        }
		
		double step = 1.0;
		while (v[AD]+diff[AD]*step < boundary && step > TOL4) {
			step /= 2.0;
		}

		if (v[AD]+diff[AD]*step < boundary && step <= TOL4) {
			ll = 5; goto end;
		}

		if (fabs(loglike_old-*loglike)/fabs(loglike_old) < TOL6 && length(U, nfix) < 0.001) {
			break;
		}

		//--------------------------------
		// halving algorithm
		//--------------------------------

		step *= 2.0;

		double logLnew;
        do {   //do loop

			step /= 2.0;

            for (i = 0; i < nfix; i++) {
				vnew[i] = v[i] + diff[i]*step;
			}

			pn = 0.5;
			theta = vnew[AD];
			
			for (r = 1; r <= max_N; r++) r0[r] = (r-1)*theta;

			//beta-binomial distribution
			p0_r0[0] = 0;
			l_r0[0] = 0;
        
			for (r = 1; r <= max_M; r++) {
				p0_r0[r]    = p0_r0[r-1] + log(  pn + r0[r]);
			}
			for (r = 1; r <= max_N; r++) {
				l_r0[r]   = l_r0[r-1] + log(1 + r0[r]);
			}
        
			//----------------------- likelihood -----------------------//

			logLnew = 0.0;
			for (p0 = 0; p0 < d.n; p0++) {
				if (d.rm[p0]) continue;

				d.sum[p0] = 0.0;

				int index1 = d.M[p0];
				int index2 = d.N[p0] - d.M[p0];
				int index3 = d.N[p0];

				for (int pp = 0; pp < d.nhap[p0]; pp++) {
					d.sum[p0] += p.pai[d.hap1[p0][pp]]*p.pai[d.hap2[p0][pp]];
				}

				//likelihood
				logLnew -= p0_r0[index1] + p0_r0[index2] - l_r0[index3] + log(d.sum[p0]);  
            
			}//person

		                          
			//if (p.ASE) {
			//	if (d.ext) link_orig_ext(&d, &p, &logLnew, vnew);
			//}

		} while (*loglike <= logLnew - TOL6 && step > TOL4); //inner_while

		if (*loglike <= logLnew - TOL6 && step <= TOL4) { 
			ll = 3; goto end; 
		}


		//--------------------------------
		// end of halving algorithm
		//--------------------------------


		for (k = 0; k < nfix; k++) v[k] += diff[k] * step;
			
		//diff

		for (j = 0; j < nfix; j++) diff[j] = v[j] - vold[j];
		

		(iter)++;
	
	}//while

	if ((iter == ITMAX0) || (v[AD] < boundary)) {
		ll = 2; 
		goto end;
	}

    for (i = 0; i < nfix; i++)	v0[i] = v[i];
    for (i = 0; i < nfix; i++)	std0[i] = sqrt(E[i+1][i+1]);

	*loglike = -(*loglike);  //Likelihood

end:

	delete [] diff; 
    delete [] U; 
	delete [] v; 
	delete [] vnew;
    delete [] vold;
    free_dmatrix(E, 1, nfix, 1, nfix);
    
	delete [] r0;
	//beta-binomial distribution
	delete [] p0_r0;
	delete [] l_r0;
	//theta
	delete [] l_r0_r;
	delete [] p0_r0_r;
	//theta2
	delete [] l_r0_r_2;
	delete [] p0_r0_r_2;

	return ll;

}//fit_ASE_null


int fit_TReC_null(double*v0, double*std0, double*loglike, DATA d, int R, int disp_TREC, char*disp_TREC_model) {

	//disp_TREC + other reg of TReC excluding SNP effect

	int TD = 0;
	int REG = disp_TREC; 
	int nfix = REG+R;
    
	int p0, i, j, k, ll = 0, iter = 0;

	double*diff  = new double[nfix];
	double*v     = new double[nfix];
	double*vnew  = new double[nfix];
    double*vold  = new double[nfix];
	double*U     = new double[nfix];
	double**E    = dmatrix(1, nfix, 1, nfix);
    
	//initialization

	for (j = 0; j < nfix; j++) v[j] = v0[j];


	//TREC variables
	int ifault = 0;
	double *vecR = new double[R];
	double *mu_i = new double[d.n];

	//--------------- EM-algorithm: within M-step, one-step Newton-Raphson --------------//

	double loglike_old = 1;
	*loglike = 2;
	for (j = 0; j < nfix; j++) U[j] = 0.0; 
    U[0] = 1.0;

	while ((fabs(loglike_old-*loglike)/fabs(loglike_old) > TOL6 || length(U, nfix) > 0.001) && iter < ITMAX0) {

		for (j = 0; j < nfix; j++) vold[j] = v[j];
		loglike_old = *loglike;

		//----------------------- likelihood -----------------------//
        
        //save variables
        
		double phi;
		
		if (disp_TREC) {
			phi = v[TD]; 
		}

		d.n_new = 0;
		*loglike = 0.0;
		for (p0 = 0; p0 < d.n; p0++) {
			if (d.rm[p0]) continue;     
		
            d.sum[p0] = 0.0;
            
			//TReC
			Rvec(vecR, R, d.X[p0]);
			double tmp = 0.0;
			for (j = 0; j < R; j++) {
				tmp += v[REG+j] * vecR[j]; 
			}
			double etmp = exp(tmp);

			mu_i[p0] = etmp;

			double mu = mu_i[p0];
			double T = d.T[p0];

			if (disp_TREC && strcmp(disp_TREC_model, "NB")==0) {
				d.sum[p0] += lgamma(T+phi) - d.lgamma_int[(int)T] - (phi+T)*log(phi+mu) + T*tmp;
			}
			else if (disp_TREC && strcmp(disp_TREC_model, "GP")==0) {
				d.sum[p0] += tmp + (T-1)*log(mu+T*phi) - mu - T*phi;
			}
			else if (!disp_TREC) {
				d.sum[p0] += T*tmp - mu;
			}

            //likelihood
			*loglike -= d.sum[p0];  
			d.n_new++;
            
		}//person

		//likelihood (TReC)
		if (disp_TREC && strcmp(disp_TREC_model, "NB")==0) {
			*loglike -= d.n_new * (- lgamma(phi) + phi*log(phi));
		}

		//----------------------- Newton Raphson -----------------------//
        
        //U, E initialation
        //U = l', E = -l"
		for (i = 0; i < nfix; i++) {
			U[i] = 0;
			for (j = i; j < nfix; j++) E[i+1][j+1] = 0;
		}
		
		for (p0 = 0; p0 < d.n; p0++) {
            if (d.rm[p0]) continue;

			Rvec(vecR, R, d.X[p0]);

			double T = d.T[p0];
			double mu = mu_i[p0];

			if (disp_TREC && strcmp(disp_TREC_model, "NB")==0) {
				for (j = 0; j < R; j++) {
					U[REG+j] += (T/mu - (phi+T)/(phi+mu)) * mu * vecR[j];
				}
				U[TD] += digama(T+phi, &ifault) - (phi+T)/(phi+mu) - log(phi+mu);
				if (ifault == 1) {ll = 4; goto end;}

				E[TD+1][TD+1] += -trigamma(T+phi, &ifault) + (mu-T)/((phi+mu)*(phi+mu)) + 1/(phi+mu);
				if (ifault == 1) {ll = 4; goto end;}
				for (j = 0; j < R; j++) {
					for (k = j; k < R; k++) {
						E[REG+j+1][REG+k+1] += (T/(mu*mu) - (phi+T)/((phi+mu)*(phi+mu))) *mu*mu*vecR[j]*vecR[k] - (T/mu - (phi+T)/(phi+mu)) *mu*vecR[j]*vecR[k];
					}
					E[TD+1][REG+j+1] += (1/(phi+mu) - (phi+T)/((phi+mu)*(phi+mu))) * mu*vecR[j];
				}
			}
			if (disp_TREC && strcmp(disp_TREC_model, "GP")==0) {
				double mu_T_phi = mu+T*phi;
				for (j = 0; j < R; j++) {
					U[REG+j] += (1.0/mu + (T-1)/mu_T_phi - 1) * mu * vecR[j];
				}
				U[TD] += T*(T-1)/mu_T_phi - T;
				E[TD+1][TD+1] += T*T*(T-1)/mu_T_phi/mu_T_phi;
				for (j = 0; j < R; j++) {
					for (k = j; k < R; k++) {
						E[REG+j+1][REG+k+1] += (1.0/(mu*mu) + (T-1)/mu_T_phi/mu_T_phi) *mu*mu*vecR[j]*vecR[k] - (1.0/mu + (T-1)/mu_T_phi - 1) *mu*vecR[j]*vecR[k];
					}
					E[TD+1][REG+j+1] += T*(T-1)/mu_T_phi/mu_T_phi * mu*vecR[j];
				}
			}
			else if (!disp_TREC) {
				for (j = 0; j < R; j++) {
					U[REG+j] += (T/mu - 1) * mu * vecR[j];
					for (k = j; k < R; k++) {
						E[REG+j+1][REG+k+1] += T*vecR[j]*vecR[k] - (T/mu - 1) *mu*vecR[j]*vecR[k];
					}
				}
			}
		}//person

		if (disp_TREC && strcmp(disp_TREC_model, "NB")==0) {
			//common terms
			U[TD]           += d.n_new * ( - digama(phi, &ifault) + 1 + log(phi) );
			if (ifault == 1) {ll = 4; goto end;}
			E[TD+1][TD+1] += d.n_new * (trigamma(phi, &ifault) - 1/phi );
			if (ifault == 1) {ll = 4; goto end;}
		}

		for (k = 0; k < nfix; k++) {
            for (j = 0; j < k; j++) E[k+1][j+1] = E[j+1][k+1];
		}

        if ((ll = explosion_E(E, nfix))) {ll = 4; goto end;}
        invv(E, nfix);
        if ((ll = explosion_E(E, nfix))) {ll = 4; goto end;}

        for (k = 0; k < nfix; k++) {
			diff[k] = 0.0;
            for (j = 0; j < nfix; j++) {
				diff[k] += E[k+1][j+1] * U[j];
			}
        }
		
		double step = 1.0;

		if (disp_TREC && strcmp(disp_TREC_model, "NB")==0) {
			while (v[TD]+diff[TD]*step < 0.0 && step > TOL4) {
				step /= 2.0;
			}

			if (v[TD]+diff[TD]*step < 0.0 && step <= TOL4) {
				ll = 5; goto end;
			}
		}

		if (fabs(loglike_old-*loglike)/fabs(loglike_old) < TOL6 && length(U, nfix) < 0.001) {
			break;
		}

		//--------------------------------
		// halving algorithm
		//--------------------------------

		step *= 2.0;

		double logLnew;
        do {   //do loop

			step /= 2.0;

            for (i = 0; i < nfix; i++) {
				vnew[i] = v[i] + diff[i]*step;
			}

			if (disp_TREC) {
				phi = vnew[TD]; 
			}

			logLnew = 0.0;

			for (p0 = 0; p0 < d.n; p0++) {
                if (d.rm[p0]) continue;

				d.sum[p0] = 0.0;
            
				//TReC
				Rvec(vecR, R, d.X[p0]);
				double tmp = 0.0;
				for (j = 0; j < R; j++) {
					tmp += vnew[REG+j] * vecR[j]; 
				}
				double etmp = exp(tmp);

				mu_i[p0] = etmp;

				double mu = mu_i[p0];
				double T = d.T[p0];

				if (disp_TREC && strcmp(disp_TREC_model, "NB")==0) {
					d.sum[p0] += lgamma(T+phi) - d.lgamma_int[(int)T] - (phi+T)*log(phi+mu) + T*tmp;
				}
				else if (disp_TREC && strcmp(disp_TREC_model, "GP")==0) {
					d.sum[p0] += tmp + (T-1)*log(mu+T*phi) - mu - T*phi;
				}
				else if (!disp_TREC) {
					d.sum[p0] += T*tmp - mu;
				}

				//likelihood
				logLnew -= d.sum[p0];  
            
			}//person

			//likelihood (TReC)
			if (disp_TREC && strcmp(disp_TREC_model, "NB")==0) {
				logLnew -= d.n_new * (- lgamma(phi) + phi*log(phi));
			}

		} while (*loglike <= logLnew - TOL6 && step > TOL4); //inner_while

		if (*loglike <= logLnew - TOL6 && step <= TOL4) { ll = 3; goto end; }


		//--------------------------------
		// end of halving algorithm
		//--------------------------------

		for (k = 0; k < nfix; k++) v[k] += diff[k] * step;

		//diff

		for (j = 0; j < nfix; j++) diff[j] = v[j] - vold[j];
		
		(iter)++;
	
	}//while

	if (iter == ITMAX0) {
		ll = 2; 
		goto end;
	}

    for (i = 0; i < nfix; i++)	v0[i] = v[i];
    for (i = 0; i < nfix; i++)	std0[i] = sqrt(E[i+1][i+1]);

	*loglike = -(*loglike);  //Likelihood

end:

	delete [] diff; 
    delete [] U; 
	delete [] v; 
	delete [] vnew;
    delete [] vold;
    free_dmatrix(E, 1, nfix, 1, nfix);

	delete [] vecR;
	delete [] mu_i;

	return ll;

}//fit_TReC_null 

int get_score_disp_TREC(double*score_disp_TREC, double*v, DATA d, int R) {

	//Dean's score test for over-dispersion

	int REG = 0;
    
	double *vecR = new double[R];
		
	double scoreNum = 0.0;
	double scoreDen = 0.0;

	for (int p0 = 0; p0 < d.n; p0++) {
        if (d.rm[p0]) continue;    

		Rvec(vecR, R, d.X[p0]);

		double mu = 0.0;
		for (int j = 0; j < R; j++) {
			mu += v[REG+j] * vecR[j]; 
		}
		mu = exp(mu);

		double T = d.T[p0];

		scoreNum += (T-mu) * (T-mu) - T;
		scoreDen += mu*mu;

	}//person

	*score_disp_TREC = scoreNum/sqrt(2.0*scoreDen); 

	delete [] vecR;

	return 0;

}//get_score_disp_TREC



int fit_TReCASE_unphased(double*v0, double*std0, double*loglike, DATA d, PARA p, double*vASE, double*vTREC) {

	int i, j, k, r, pp, p0;
    
	int ll = 0, iter = 0;

	double*diff  = new double[p.nall];
	double*v     = new double[p.nall];
	double*vnew  = new double[p.nall];
    double*vold  = new double[p.nall];
	double*U     = new double[p.nall];
	double**E    = dmatrix(1, p.nall, 1, p.nall);
    
	//initialization
    
    for (j = 0; j < p.B; j++) v[j] = 0.0;								//beta: SNP-related reg
	if (p.ASE) {
		v[p.AD] = vASE[0];												//theta: ASE dispersion
		for (j = 0; j < p.KK; j++) v[p.FREQ+j] = p.pai[p.hap_freq[j]];	//pi: hap freq
	}
	if (p.TREC) {
		if (p.disp_TREC) v[p.TD] = vTREC[0];							//phi: TReC dispersion   
		for (j = 0; j < p.R;  j++) v[p.REG+j] = vTREC[p.disp_TREC+j];	//eta: other reg of TReC excluding SNPs
	}

	//ASE & TRECT variables
	double *vecB = new double[p.B];
	double *exp_beta = new double[d.n];

    //ASE variables
	int pow1, max_M, max_N;
	double *Hi = NULL, *Hikl = NULL;
    double pn = 0.5, theta, boundary;
	double *p1, *r0, *p0_r0, *l_r0, *l_r0_r, *p0_r0_r, *l_r0_r_2, *p0_r0_r_2;
    
	if (p.ASE) {

		pow1 = p.pow2[p.M-1];

		max_M = 0;
		max_N = 0;
		for (int p0 = 0; p0 < d.n; p0++) {
	     if (d.rm[p0]) continue;    
			if (d.M[p0] > max_M)               max_M = d.M[p0];
			if (d.N[p0]-d.M[p0] > max_M)	   max_M = d.N[p0]-d.M[p0];
			if (d.N[p0] > max_N)               max_N = d.N[p0];
		}    

		p1	  = new double[d.n];

		r0          = new double[max_N+1];
		//beta-binomial distribution
		p0_r0       = new double[max_M+1]; //sum_r log(p0+r0)
		l_r0        = new double[max_N+1]; //sum_r log(1+r0)    
		//theta
		l_r0_r      = new double[max_N+1]; //\sum_r r/(1+r0)
		p0_r0_r     = new double[max_M+1]; //\sum_r r/(p0+r0)
		//theta2
		l_r0_r_2    = new double[max_N+1]; //\sum_r [r/(1+r0)]^2
		p0_r0_r_2   = new double[max_M+1]; //\sum_r [r/(p0+r0)]^2       
	}//ASE

	//TREC variables
	int ifault = 0;
	double *vecR, phi, *mu_i, *tmu_i;
	if (p.TREC) {
		vecR  = new double[p.R];
		mu_i  = new double[d.n];
		tmu_i = new double[d.n];
	}//TREC

	//--------------- EM-algorithm: within M-step, one-step Newton-Raphson --------------//

	double loglike_old = 1;
	*loglike = 2;
	for (j = 0; j < p.nall; j++) U[j] = 0.0; 
    U[0] = 1.0;

	while ((fabs(loglike_old-*loglike)/fabs(loglike_old) > TOL6 || length(U, p.nfix) > 0.001) && iter < ITMAX0) {
	//while (length(diff, p.nall) > TOL6 && iter < ITMAX0) {

		for (j = 0; j < p.nall; j++) vold[j] = v[j];
		loglike_old = *loglike;

		//----------------------- E-step -----------------------//
        
        //save variables that are common to individuals
		
		if (p.ASE) {

			theta = v[p.AD];
		
			for (r = 1; r <= max_N; r++) r0[r] = (r-1)*theta;

			//beta-binomial distribution
			p0_r0[0] = 0; 	for (r = 1; r <= max_M; r++) p0_r0[r]  = p0_r0[r-1]  + log(pn + r0[r]);
			l_r0[0]  = 0;	for (r = 1; r <= max_N; r++)  l_r0[r]  =  l_r0[r-1]  + log( 1 + r0[r]);
			
			boundary = -100; //MAX(-pn/(max_M-1),-1.0/(max_N-1));
		}//ASE
        
		if (p.TREC) {
			if (p.disp_TREC) phi = v[p.TD]; 
		}

		d.n_new = 0.0;
		*loglike = 0.0;
		for (p0 = 0; p0 < d.n; p0++) {
			if (d.rm[p0]) continue;    

			//exp_beta[d.n]

			exp_beta[p0] = 1.0;

			if (d.G[p0] >= 1) {
				Bvec(vecB, p.B, d.X[p0]);

				double beta_X = 0.0;
				for (j = 0; j < p.B; j++) {
					beta_X += v[j] * vecB[j];
				}
				exp_beta[p0] = exp(beta_X);
			}

            d.sum[p0] = 0.0;

			//ASE
			if (p.ASE){

				int index1 = d.M[p0];
				int index2 = d.N[p0] - d.M[p0];
				int index3 = d.N[p0];

				if (d.G[p0] == 1 && d.N[p0] > 0) { //two pairs

					//p1[d.n]
					p1[p0] = exp_beta[p0]/(1+exp_beta[p0]);

					for (int pp = 0; pp < d.nhap[p0]; pp++) {  
 
						if (d.hap1[p0][pp] / pow1 == 0) {     
							index1 = d.N[p0] - d.M[p0];
							index2 = d.M[p0];
						}
						else {
							index1 = d.M[p0];
							index2 = d.N[p0] - d.M[p0];
						}

						if (index3 > 1) boundary = MAX(boundary, -1.0/(index3-1));
						if (index1 > 1) boundary = MAX(boundary, -p1[p0]/(index1-1));
						if (index2 > 1) boundary = MAX(boundary, -(1-p1[p0])/(index2-1));

						//beta-binomial distribution
						double p_r0 = 0;	for (r = 1; r <= index1; r++)   p_r0 = p_r0   + log(  p1[p0] + r0[r]);
						double l_p_r0 = 0;  for (r = 1; r <= index2; r++) l_p_r0 = l_p_r0 + log(1-p1[p0] + r0[r]);

						d.w[p0][pp] = p_r0 + l_p_r0 - l_r0[index3] + log(v[p.FREQ+p.h[d.hap1[p0][pp]]]*v[p.FREQ+p.h[d.hap2[p0][pp]]]);
					}//pp

					if (d.nhap[p0]==2) {
						double w1 = d.w[p0][0];
						double w2 = d.w[p0][1];

						d.w[p0][0] = 1.0/(1.0+exp(w2-w1));
						d.w[p0][1] = 1.0/(1.0+exp(w1-w2));

						d.sum[p0] = w1 + log(1+exp(w2-w1));
					}
					else {
						d.sum[p0] = d.w[p0][0];
						d.w[p0][0] = 1.0;
					}

				}
				else if (d.G[p0] != 1 && d.N[p0] > 0) { //l pair
					d.w[p0][0] = 1.0;
					d.sum[p0] += p0_r0[index1] + p0_r0[index2] - l_r0[index3] + log(v[p.FREQ+p.h[d.hap1[p0][0]]]*v[p.FREQ+p.h[d.hap2[p0][0]]]);

					if (index3 > 1) boundary = MAX(boundary, -1.0/(index3-1));
					if (index1 > 1) boundary = MAX(boundary, -pn/(index1-1));
					if (index2 > 1) boundary = MAX(boundary, -pn/(index2-1));

				}
				else if (d.N[p0] == 0) {                   //can be two pairs or l pair
					for (pp = 0; pp < d.nhap[p0]; pp++) {
						d.w[p0][pp] = v[p.FREQ+p.h[d.hap1[p0][pp]]]*v[p.FREQ+p.h[d.hap2[p0][pp]]];
						d.sum[p0] += d.w[p0][pp];
					}
					for (pp = 0; pp < d.nhap[p0]; pp++) d.w[p0][pp] /= d.sum[p0];
					d.sum[p0] = log(d.sum[p0]);
				}

			}//ASE
            
			//TReC
			if (p.TREC) {

				double exp_Z;
				if (d.G[p0] == 0)		exp_Z = 1;
				else if (d.G[p0] == 1) exp_Z = 0.5*(1+exp_beta[p0]);
				else						exp_Z = exp_beta[p0];

				Rvec(vecR, p.R, d.X[p0]);
				double exp_eta = 0.0;
				for (j = 0; j < p.R; j++) {
					exp_eta += v[p.REG+j] * vecR[j]; 
				}
				exp_eta = exp(exp_eta);

				mu_i[p0] = exp_eta * exp_Z;
				tmu_i[p0] = 0.5 * d.G[p0] * exp_eta * exp_beta[p0];

				double mu = mu_i[p0];
				double T = d.T[p0];

				if (p.disp_TREC && strcmp(p.disp_TREC_model, "NB")==0) {
					d.sum[p0] += lgamma(T+phi) - d.lgamma_int[(int)T] - (phi+T)*log(phi+mu) + T*log(mu);
				}
				else if (p.disp_TREC && strcmp(p.disp_TREC_model, "GP")==0) {
					d.sum[p0] += log(mu) + (T-1)*log(mu+T*phi) - mu - T*phi;
				}
				else if (!p.disp_TREC) {
					d.sum[p0] += T*log(mu) - mu;
				}
			}//TReC

            //likelihood
			*loglike -= d.sum[p0];  
			d.n_new++;
            
		}//person

		//likelihood (TReC)
		if (p.TREC && p.disp_TREC && strcmp(p.disp_TREC_model, "NB")==0) {
			*loglike -= d.n_new * (- lgamma(phi) + phi*log(phi));
		}
		
		//----------------------- M-step -----------------------//
        
		if (p.ASE) {

			l_r0_r[0] = 0;   //theta
			p0_r0_r[0] = 0; 
			l_r0_r_2[0] = 0; //theta2
			p0_r0_r_2[0] = 0;
			for (r = 1; r <= max_M; r++) {
				p0_r0_r[r]      = p0_r0_r[r-1]  + (r-1)/  (pn + r0[r]);   
				p0_r0_r_2[r]    = p0_r0_r_2[r-1]  + (r-1)*(r-1)/  ((pn+r0[r])*  (pn+r0[r]));
			}      
			for (r = 1; r <= max_N; r++) {
				l_r0_r[r]     = l_r0_r[r-1] + (r-1)/(1+r0[r]);
				l_r0_r_2[r]   = l_r0_r_2[r-1] + (r-1)*(r-1)/((1+r0[r])*(1+r0[r]));
			}
		}//ASE

        //U, E initialation
        //U = l', E = -l"
		for (i = 0; i < p.nall; i++) {
			U[i] = 0;
			for (j = i; j < p.nall; j++) E[i+1][j+1] = 0;
		}
		
		for (p0 = 0; p0 < d.n; p0++) {
	      if (d.rm[p0]) continue;    
            
			Bvec(vecB, p.B, d.X[p0]);

			//ASE
			if (p.ASE){
				//default
				int index1 = d.M[p0];
				int index2 = d.N[p0] - d.M[p0];
				int index3 = d.N[p0];
            
				if (d.G[p0] == 1 && d.N[p0] > 0) {            //N>0, G=1
					for (int pp = 0; pp < d.nhap[p0]; pp++) {   
                
						if (d.hap1[p0][pp] / pow1 == 0) {     
							index1 = d.N[p0] - d.M[p0];
							index2 = d.M[p0];
						}
						else {
							index1 = d.M[p0];
							index2 = d.N[p0] - d.M[p0];
						}

						double p_r0_1 = 0;			//beta
						double l_p_r0_1 = 0;
						double p_r0_r = 0;			//theta
						double l_p_r0_r = 0;
						double p_r0_1_2 = 0;		//beta2
						double l_p_r0_1_2 = 0; 			
						double p_r0_r_d2 = 0;		//beta*theta
						double l_p_r0_r_d2 = 0;			
						double p_r0_r_2 = 0;		//theta2
						double l_p_r0_r_2 = 0;
					
						double p1_tmp = p1[p0];
						double l_p1_tmp = 1-p1_tmp;
						for (r = 1; r <= index1; r++) {
							p_r0_1       = p_r0_1   + 1/(p1_tmp + r0[r]);
							p_r0_r       = p_r0_r   + (r-1)/(p1_tmp + r0[r]);
							p_r0_1_2     = p_r0_1_2 + 1/((p1_tmp+r0[r])*(p1_tmp+r0[r]));
							p_r0_r_d2    = p_r0_r_d2+ (r-1)/((p1_tmp+r0[r])*(p1_tmp+r0[r]));
							p_r0_r_2     = p_r0_r_2 + (r-1)*(r-1)/((p1_tmp+r0[r])*(p1_tmp+r0[r]));
						} 
						for (r = 1; r <= index2; r++) {
							l_p_r0_1     = l_p_r0_1   + 1/(l_p1_tmp + r0[r]);
							l_p_r0_r     = l_p_r0_r   + (r-1)/(l_p1_tmp + r0[r]);
							l_p_r0_1_2   = l_p_r0_1_2 + 1/((l_p1_tmp+r0[r])*(l_p1_tmp+r0[r]));
							l_p_r0_r_d2  = l_p_r0_r_d2+ (r-1)/((l_p1_tmp+r0[r])*(l_p1_tmp+r0[r]));
							l_p_r0_r_2   = l_p_r0_r_2 + (r-1)*(r-1)/((l_p1_tmp+r0[r])*(l_p1_tmp+r0[r]));
						} 

						double exp_beta_tmp = exp_beta[p0];
						for (k = 0; k < p.B; k++) {
							U[k]           += d.w[p0][pp] * ( p_r0_1    - l_p_r0_1) * p1_tmp * l_p1_tmp * vecB[k];
							E[k+1][p.AD+1] -= d.w[p0][pp] * (-p_r0_r_d2 + l_p_r0_r_d2) * p1_tmp * l_p1_tmp * vecB[k];
							for (int l = k; l < p.B; l++) {
								E[k+1][l+1] -= d.w[p0][pp] * (-p_r0_1_2  - l_p_r0_1_2) * p1_tmp*p1_tmp*l_p1_tmp*l_p1_tmp * vecB[k]*vecB[l]
								             + d.w[p0][pp] * ( p_r0_1    - l_p_r0_1) * exp_beta_tmp * (1-exp_beta_tmp)/((1+exp_beta_tmp)*(1+exp_beta_tmp)*(1+exp_beta_tmp)) 
											         * vecB[k]*vecB[l];
							}
						}
						U[p.AD]           += d.w[p0][pp] * ( p_r0_r    + l_p_r0_r - l_r0_r[index3]);
						E[p.AD+1][p.AD+1] -= d.w[p0][pp] * (-p_r0_r_2  - l_p_r0_r_2 + l_r0_r_2[index3]);
					}
				}
				else if (d.G[p0] != 1 && d.N[p0] > 0) {      //N>0, G=0, 2
					for (int pp = 0; pp < d.nhap[p0]; pp++) {     
						U[p.AD]           += d.w[p0][pp] * (p0_r0_r[index1] + p0_r0_r[index2] - l_r0_r[index3]);
						E[p.AD+1][p.AD+1] -= d.w[p0][pp] * (-p0_r0_r_2[index1] - p0_r0_r_2[index2] + l_r0_r_2[index3]);
					}
				}
			}//ASE

			//TREC
			if (p.TREC) {
				Rvec(vecR, p.R, d.X[p0]);

				double T = d.T[p0];
				double mu = mu_i[p0];
				double tmu = tmu_i[p0];

				if (p.disp_TREC && strcmp(p.disp_TREC_model, "NB")==0) {

					for (j = 0; j < p.B; j++) {
						U[j]           += (T/mu - (phi+T)/(phi+mu)) * tmu * vecB[j];
						for (k = j; k < p.B; k++) {
							E[j+1][k+1]+= ((T/(mu*mu) - (phi+T)/((phi+mu)*(phi+mu))) *tmu*tmu - (T/mu - (phi+T)/(phi+mu)) * tmu) * vecB[j]*vecB[k];
						}
						E[j+1][p.TD+1] += (1/(phi+mu) - (phi+T)/((phi+mu)*(phi+mu))) * tmu * vecB[j];
						for (k = 0; k < p.R; k++) {
							E[j+1][p.REG+k+1] += (T/(mu*mu) - (phi+T)/((phi+mu)*(phi+mu))) * tmu*mu*vecB[j]*vecR[k] - (T/mu - (phi+T)/(phi+mu)) * tmu*vecB[j]*vecR[k];
						}
					}
					for (j = 0; j < p.R; j++) {
						U[p.REG+j] += (T/mu - (phi+T)/(phi+mu)) * mu * vecR[j];
						for (k = j; k < p.R; k++) {
							E[p.REG+j+1][p.REG+k+1] += (T/(mu*mu) - (phi+T)/((phi+mu)*(phi+mu))) *mu*mu*vecR[j]*vecR[k] - (T/mu - (phi+T)/(phi+mu)) *mu*vecR[j]*vecR[k];
						}
						E[p.TD+1][p.REG+j+1] += (1/(phi+mu) - (phi+T)/((phi+mu)*(phi+mu))) * mu*vecR[j];
					}

					U[p.TD] += digama(T+phi, &ifault) - (phi+T)/(phi+mu) - log(phi+mu);
					if (ifault == 1) {ll = 4; goto end;}
					E[p.TD+1][p.TD+1] += -trigamma(T+phi, &ifault) + (mu-T)/((phi+mu)*(phi+mu)) + 1/(phi+mu);
					if (ifault == 1) {ll = 4; goto end;}

				}//with dispersion
				if (p.disp_TREC && strcmp(p.disp_TREC_model, "GP")==0) {

					double mu_T_phi = mu+T*phi;
					for (j = 0; j < p.B; j++) {
						U[j]           += (1.0/mu + (T-1)/mu_T_phi - 1) * tmu * vecB[j];
						for (k = j; k < p.B; k++) {
							E[j+1][k+1]+= ((1.0/(mu*mu) + (T-1)/mu_T_phi/mu_T_phi) *tmu*tmu - (1.0/mu + (T-1)/mu_T_phi - 1) *tmu) * vecB[j]*vecB[k];
						}
						E[j+1][p.TD+1] += T*(T-1)/mu_T_phi/mu_T_phi * tmu * vecB[j];
						for (k = 0; k < p.R; k++) {
							E[j+1][p.REG+k+1] += ((1.0/(mu*mu) + (T-1)/mu_T_phi/mu_T_phi) *tmu*mu - (1.0/mu + (T-1)/mu_T_phi - 1) *tmu)*vecB[j]*vecR[k];
						}
					}
					for (j = 0; j < p.R; j++) {
						U[p.REG+j] += (1.0/mu + (T-1)/mu_T_phi - 1) * mu * vecR[j];
						for (k = j; k < p.R; k++) {
							E[p.REG+j+1][p.REG+k+1] += ((1.0/(mu*mu) + (T-1)/mu_T_phi/mu_T_phi) *mu*mu - (1.0/mu + (T-1)/mu_T_phi - 1) *mu)*vecR[j]*vecR[k];
						}
						E[p.TD+1][p.REG+j+1] += T*(T-1)/mu_T_phi/mu_T_phi * mu*vecR[j];
					}

					U[p.TD] += T*(T-1)/mu_T_phi - T;
					E[p.TD+1][p.TD+1] += T*T*(T-1)/mu_T_phi/mu_T_phi;

				}//with dispersion
				else if (!p.disp_TREC) {

					for (j = 0; j < p.B; j++) {
						U[j]           += (T/mu - 1) * tmu * vecB[j];
						for (k = j; k < p.B; k++) {
							E[j+1][k+1]+= (T/(mu*mu) *tmu*tmu - (T/mu - 1) * tmu) * vecB[j]*vecB[k];
						}
						for (k = 0; k < p.R; k++) {
							E[j+1][p.REG+k+1] += T/(mu*mu) * tmu*mu*vecB[j]*vecR[k] - (T/mu - 1) * tmu*vecB[j]*vecR[k];
						}
					}
					for (j = 0; j < p.R; j++) {
						U[p.REG+j] += (T/mu - 1) * mu * vecR[j];
						for (k = j; k < p.R; k++) {
							E[p.REG+j+1][p.REG+k+1] += T*vecR[j]*vecR[k] - (T/mu - 1) *mu*vecR[j]*vecR[k];
						}
					}

				}//no dispersion

			}//TREC
		}//person

		if (p.TREC && p.disp_TREC && strcmp(p.disp_TREC_model, "NB")==0) {
			//common terms
			U[p.TD]           += d.n_new * ( - digama(phi, &ifault) + 1 + log(phi) );
			if (ifault == 1) {ll = 4; goto end;}
			E[p.TD+1][p.TD+1] += d.n_new * (trigamma(phi, &ifault) - 1/phi );
			if (ifault == 1) {ll = 4; goto end;}
		}

		for (k = 0; k < p.nfix; k++) {
            for (j = 0; j < k; j++) E[k+1][j+1] = E[j+1][k+1];
		}

        if ((ll = explosion_E(E, p.nfix))) { ll = 4; goto end;}
        invv(E, p.nfix);
        if ((ll = explosion_E(E, p.nfix))) { ll = 4; goto end;}

        for (k = 0; k < p.nfix; k++) {
			diff[k] = 0.0;
            for (j = 0; j < p.nfix; j++) {
				diff[k] += E[k+1][j+1] * U[j];
			}
        }
		
		double step = 1.0;
		if (p.TREC && p.disp_TREC && strcmp(p.disp_TREC_model, "NB")==0) {
			while (v[p.TD]+diff[p.TD]*step < 0.0 && step > TOL4) {
				step /= 2.0;
			}
			if (v[p.TD]+diff[p.TD]*step < 0.0 && step <= TOL4) {
				ll = 5; goto end;
			}

		}
				
		if (p.ASE) {
			//renew \pi
			for (j = p.FREQ; j < p.FREQ+p.KK; j++) v[j] = 0; 

			for (int p0 = 0; p0 < d.n; p0++) {
				if (d.rm[p0]) continue;    
				for (int pp = 0; pp < d.nhap[p0]; pp++) {    
					v[p.FREQ+p.h[d.hap1[p0][pp]]] += d.w[p0][pp] * 0.5;
					v[p.FREQ+p.h[d.hap2[p0][pp]]] += d.w[p0][pp] * 0.5;
				}
			}

			for (j = p.FREQ; j < p.FREQ+p.KK; j++) v[j] /= d.n_new;
		}


	    if ((fabs(loglike_old-*loglike)/fabs(loglike_old) < TOL6 && length(U, p.nfix) < 0.001)) {
			break;
		}

		//--------------------------------
		// halving algorithm
		//--------------------------------

		step *= 2.0;

		double logLnew;
        do {   //do loop

			step /= 2.0;

            for (i = 0; i < p.nfix; i++) {
				vnew[i] = v[i] + diff[i]*step;
			}

			//save variables that are common to individuals
			
			if (p.ASE) {

				theta = vnew[p.AD];
			
				for (r = 1; r <= max_N; r++) r0[r] = (r-1)*theta;

				//beta-binomial distribution
				p0_r0[0] = 0; 	for (r = 1; r <= max_M; r++) p0_r0[r]  = p0_r0[r-1]  + log(pn + r0[r]);
				l_r0[0]  = 0;	for (r = 1; r <= max_N; r++)  l_r0[r]  =  l_r0[r-1]  + log( 1 + r0[r]);
				
				boundary = -100; //MAX(-pn/(max_M-1),-1.0/(max_N-1));

			}//ASE
        
			if (p.TREC) {
				if (p.disp_TREC) phi = vnew[p.TD]; 
			}

			logLnew = 0.0;
			for (int p0 = 0; p0 < d.n; p0++) {
				if (d.rm[p0]) continue;    

				//exp_beta[d.n]

				exp_beta[p0] = 1.0;

				if (d.G[p0] >= 1) {
					Bvec(vecB, p.B, d.X[p0]);

					double beta_X = 0.0;
					for (j = 0; j < p.B; j++) {
						beta_X += vnew[j] * vecB[j];
					}
					exp_beta[p0] = exp(beta_X);
				}

				d.sum[p0] = 0.0;

				//ASE
				if (p.ASE){

					int index1 = d.M[p0];
					int index2 = d.N[p0] - d.M[p0];
					int index3 = d.N[p0];

					if (d.G[p0] == 1 && d.N[p0] > 0) { //two pairs

						//p1[d.n]
						p1[p0] = exp_beta[p0]/(1+exp_beta[p0]);

						for (int pp = 0; pp < d.nhap[p0]; pp++) {  
 
							if (d.hap1[p0][pp] / pow1 == 0) {     
								index1 = d.N[p0] - d.M[p0];
								index2 = d.M[p0];
							}
							else {
								index1 = d.M[p0];
								index2 = d.N[p0] - d.M[p0];
							}

							if (index3 > 1) boundary = MAX(boundary, -1.0/(index3-1));
							if (index1 > 1) boundary = MAX(boundary, -p1[p0]/(index1-1));
							if (index2 > 1) boundary = MAX(boundary, -(1-p1[p0])/(index2-1));

							//beta-binomial distribution
							double p_r0 = 0;	for (r = 1; r <= index1; r++)   p_r0 = p_r0   + log(  p1[p0] + r0[r]);
							double l_p_r0 = 0;  for (r = 1; r <= index2; r++) l_p_r0 = l_p_r0 + log(1-p1[p0] + r0[r]);

							d.w[p0][pp] = p_r0 + l_p_r0 - l_r0[index3] + log(v[p.FREQ+p.h[d.hap1[p0][pp]]]*v[p.FREQ+p.h[d.hap2[p0][pp]]]);
						}//pp

						if (d.nhap[p0]==2) {
							double w1 = d.w[p0][0];
							double w2 = d.w[p0][1];

							d.w[p0][0] = 1.0/(1.0+exp(w2-w1));
							d.w[p0][1] = 1.0/(1.0+exp(w1-w2));

							d.sum[p0] = w1 + log(1+exp(w2-w1));
						}
						else {
							d.sum[p0] = d.w[p0][0];
							d.w[p0][0] = 1.0;
						}
					}
					else if (d.G[p0] != 1 && d.N[p0] > 0) { //l pair
						d.w[p0][0] = 1.0;
						d.sum[p0] += p0_r0[index1] + p0_r0[index2] - l_r0[index3] + log(v[p.FREQ+p.h[d.hap1[p0][0]]]*v[p.FREQ+p.h[d.hap2[p0][0]]]);

						if (index3 > 1) boundary = MAX(boundary, -1.0/(index3-1));
						if (index1 > 1) boundary = MAX(boundary, -pn/(index1-1));
						if (index2 > 1) boundary = MAX(boundary, -pn/(index2-1));
	
					}
					else if (d.N[p0] == 0) {                   //can be two pairs or l pair
						int pp;
						for (pp = 0; pp < d.nhap[p0]; pp++) {    
							d.w[p0][pp] = v[p.FREQ+p.h[d.hap1[p0][pp]]]*v[p.FREQ+p.h[d.hap2[p0][pp]]];
							d.sum[p0] += d.w[p0][pp];
						}
						for (pp = 0; pp < d.nhap[p0]; pp++)  d.w[p0][pp] /= d.sum[p0];
						d.sum[p0] = log(d.sum[p0]);
					}

				}//ASE

				//TReC
				if (p.TREC) {

					double exp_Z;
					if (d.G[p0] == 0)		exp_Z = 1;
					else if (d.G[p0] == 1) exp_Z = 0.5*(1+exp_beta[p0]);
					else						exp_Z = exp_beta[p0];

					Rvec(vecR, p.R, d.X[p0]);
					double exp_eta = 0.0;
					for (j = 0; j < p.R; j++) {
						exp_eta += vnew[p.REG+j] * vecR[j]; 
					}
					exp_eta = exp(exp_eta);

					mu_i[p0] = exp_eta * exp_Z;
					tmu_i[p0] = 0.5 * d.G[p0] * exp_eta * exp_beta[p0];

					double mu = mu_i[p0];
					double T = d.T[p0];

					if (p.disp_TREC && strcmp(p.disp_TREC_model, "NB")==0) {
						d.sum[p0] += lgamma(T+phi) - d.lgamma_int[(int)T] - (phi+T)*log(phi+mu) + T*log(mu);
					}
					else if (p.disp_TREC && strcmp(p.disp_TREC_model, "GP")==0) {
						d.sum[p0] += log(mu) + (T-1)*log(mu+T*phi) - mu - T*phi;
					}
					else if (!p.disp_TREC) {
						d.sum[p0] += T*log(mu) - mu;
					}
				}//TReC

				//likelihood
				logLnew -= d.sum[p0];  
            
			}//person

			//likelihood (TReC)
			if (p.TREC && p.disp_TREC && strcmp(p.disp_TREC_model, "NB")==0) {
				logLnew -= d.n_new * (- lgamma(phi) + phi*log(phi));
			}

		} while ((p.ASE && vnew[p.AD] < boundary || *loglike <= logLnew - TOL6) && step > TOL4); //inner_while

		if ((p.ASE && vnew[p.AD] < boundary || *loglike <= logLnew - TOL6) && step <= TOL4) { ll = 3; goto end; }


		//--------------------------------
		// end of halving algorithm
		//--------------------------------

		for (k = 0; k < p.nfix; k++) v[k] += diff[k] * step;

		//diff

		for (j = 0; j < p.nall; j++) diff[j] = v[j] - vold[j];
		

		(iter)++;
	
	}//while

	//printf("iter = %d  n_new = %d\n", iter, d.n_new);

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//premature convergence by loglike criteria, U is still not zero, resulting in E not positive definite
	//////////////////////////////////////////////////////////////////////////////////////////////////////

	if (iter == ITMAX0) {
		ll = 2; 
		goto end;
	}
	if (p.ASE) {
		if (v[p.AD] < boundary) {
			ll = 2;
			goto end;
		}
	}

	//EM algorithm converged
    for (i = 0; i < p.nfix; i++)	v0[i] = v[i];
	if (p.ASE) {
		for (i = p.FREQ; i < p.FREQ+p.KK; i++) v0[i] = v[i];
	}

	*loglike = -(*loglike);  //Likelihood


	if (p.ASE) {

		//------------------------------ std: Louis Estimator -----------------------------------//

		//E_Louis

		//recover the information matrix
		invv(E, p.FREQ);

		Hi   = new double[p.FREQ+p.KK-1];
		Hikl = new double[p.FREQ+p.KK-1];

		for (int p0 = 0; p0 < d.n; p0++) {
			if (d.rm[p0]) continue;    

			for (j = 0; j < p.FREQ+p.KK-1; j++) Hi[j] = 0;
        
			Bvec(vecB, p.B, d.X[p0]);

			//default
			int index1 = d.M[p0];
			int index2 = d.N[p0] - d.M[p0];
			int index3 = d.N[p0];
        
			for (int pp = 0; pp < d.nhap[p0]; pp++) {     
            
				for (j = 0; j < p.FREQ+p.KK-1; j++) Hikl[j] = 0.0;

				//beta & theta
				if (d.G[p0] == 1 && d.N[p0] > 0) {            //N>0, G=1
                
					if (d.hap1[p0][pp] / pow1 == 0) {     
						index1 = d.N[p0] - d.M[p0];
						index2 = d.M[p0];
					}
					else {
						index1 = d.M[p0];
						index2 = d.N[p0] - d.M[p0];
					}

					double p_r0_1 = 0;			//beta
					double l_p_r0_1 = 0;
					double p_r0_r = 0;			//theta
					double l_p_r0_r = 0;
					
					double p1_tmp = p1[p0];
					double l_p1_tmp = 1-p1_tmp;
					for (r = 1; r <= index1; r++) {
						p_r0_1       = p_r0_1   + 1/(p1_tmp + r0[r]);
						p_r0_r       = p_r0_r   + (r-1)/(p1_tmp + r0[r]);
					} 
					for (r = 1; r <= index2; r++) {
						l_p_r0_1     = l_p_r0_1   + 1/(l_p1_tmp + r0[r]);
						l_p_r0_r     = l_p_r0_r   + (r-1)/(l_p1_tmp + r0[r]);
					} 

					for (j = 0; j < p.B; j++) Hikl[j] += ( p_r0_1 - l_p_r0_1) * p1_tmp * l_p1_tmp * vecB[j];
					Hikl[p.AD] += ( p_r0_r    + l_p_r0_r - l_r0_r[index3]);
				}
				else if (d.G[p0] != 1 && d.N[p0] > 0) {      //N>0, G=0, 2
					Hikl[p.AD] += (p0_r0_r[index1] + p0_r0_r[index2] - l_r0_r[index3]);
				}

				//\pi
				if (p.h[d.hap1[p0][pp]] != 0) {
					Hikl[p.FREQ+p.h[d.hap1[p0][pp]]-1] += 1 / v[p.FREQ+p.h[d.hap1[p0][pp]]];
					E[p.FREQ+p.h[d.hap1[p0][pp]]][p.FREQ+p.h[d.hap1[p0][pp]]] += d.w[p0][pp]/(v[p.FREQ+p.h[d.hap1[p0][pp]]]*v[p.FREQ+p.h[d.hap1[p0][pp]]]);
				}
				if (p.h[d.hap2[p0][pp]] != 0) {
					Hikl[p.FREQ+p.h[d.hap2[p0][pp]]-1] += 1 / v[p.FREQ+p.h[d.hap2[p0][pp]]];
					E[p.FREQ+p.h[d.hap2[p0][pp]]][p.FREQ+p.h[d.hap2[p0][pp]]] += d.w[p0][pp]/(v[p.FREQ+p.h[d.hap2[p0][pp]]]*v[p.FREQ+p.h[d.hap2[p0][pp]]]);
				}
				int counter = (p.h[d.hap1[p0][pp]]==0) + (p.h[d.hap2[p0][pp]]==0);
				if (counter > 0) {
					for (i = 0; i < p.KK-1; i++) {
						Hikl[p.FREQ+i] -= counter / v[p.FREQ];
						for (j = i; j < p.KK-1; j++) {
							E[p.FREQ+i+1][p.FREQ+j+1] += d.w[p0][pp]*counter/(v[p.FREQ]*v[p.FREQ]);
						}
					}
				}
			
				for (k = 0; k < p.FREQ+p.KK-1; k++) {
					Hi[k] += d.w[p0][pp] * Hikl[k];
					for (int l = k; l < p.FREQ+p.KK-1; l++) 
						E[k+1][l+1] -= d.w[p0][pp] * Hikl[k] * Hikl[l]; 
				}

  			}//(k,l)

			for (k = 0; k < p.FREQ+p.KK-1; k++) {
				for (int l = k; l < p.FREQ+p.KK-1; l++) {
					E[k+1][l+1] += Hi[k] * Hi[l]; 
				}
			} 
	
		}//person

		delete [] Hi;
		delete [] Hikl;


		for (j = 1; j <= p.FREQ+p.KK-1; j++) {
			for (k = 1; k < j; k++) {
				E[j][k] = E[k][j];
			}
		}

		if ((ll = explosion_E(E, p.FREQ+p.KK-1))) {ll=1; goto end;}
		invv(E, p.FREQ+p.KK-1);
		if ((ll = explosion_E(E, p.FREQ+p.KK-1))) {ll=1; goto end;}
		if ((ll = inproper_std(E, p.FREQ))) {ll=1; goto end;}

    }//ASE

	if ((ll = inproper_std(E, p.nfix))) {ll=1; goto end;}

    for (i = 0; i < p.nfix; i++)	std0[i] = sqrt(E[i+1][i+1]);


	if (p.ASE) {
		for (i = p.FREQ+1; i < p.FREQ+p.KK; i++) std0[i] = sqrt(E[i][i]);

		std0[p.FREQ] = 0.0;
		for (i = p.FREQ+1; i < p.FREQ+p.KK; i++) {
			for (j = p.FREQ+1; j < p.FREQ+p.KK; j++) {
				std0[p.FREQ] += E[i][j];
			}
		}
		std0[p.FREQ] = sqrt(std0[p.FREQ]);
	}

end:

	delete [] diff; 
    delete [] U; 
	delete [] v; 
	delete [] vnew;
    delete [] vold;
    free_dmatrix(E, 1, p.nall, 1, p.nall);
    
	delete [] vecB;
	delete [] exp_beta;

	if (p.ASE) {
		delete [] p1;
		delete [] r0;
		delete [] p0_r0;
		delete [] l_r0;    
		delete [] l_r0_r;
		delete [] p0_r0_r;
		delete [] l_r0_r_2;
		delete [] p0_r0_r_2;  
	}

	if (p.TREC){
		delete [] vecR;
		delete [] mu_i;
		delete [] tmu_i;
	}

	return ll;

}//fit_TReCASE_unphased


void evaluate_effiScore_trans(int*ll, double*U, double*V, double**Uei, DATA d, PARA p, double*vTRECASE) {

	int saturate_fitted = 1;

	int p0, i, j, k, r, nui=0;
    
	*ll = 0;

	double loglike; 
	double*v     = new double[p.nall];
	double*Usum  = new double[p.nall];
	double**Ui   = dmatrix(0, d.n-1, 0, p.nall-1);
	double**E    = dmatrix(1, p.nall, 1, p.nall);
	double**E12_22=NULL, **E12=NULL, **E22=NULL;

	if (saturate_fitted) {

		//E22 = E, Uei = U1i
		E12_22 = dmatrix(1, p.B, 1, p.nall);
		E12  = dmatrix(1, p.B, 1, p.nall); 
	}

    
	//initialization

	if (saturate_fitted) {
		for (j = 0; j < p.nall; j++) v[j] = vTRECASE[j];
	}

	//ASE & TRECT variables
	double *vecB = new double[p.B];
	double *exp_beta = new double[d.n];

    //ASE variables
	int pow1, max_M, max_N;
	double *Hi = NULL, *Hikl = NULL;
    double pn = 0.5, theta, boundary;
	double *p1, *r0, *p0_r0, *l_r0, *l_r0_r, *p0_r0_r, *l_r0_r_2, *p0_r0_r_2;
    
	if (p.ASE) {

		pow1 = p.pow2[p.M-1];

		max_M = 0;
		max_N = 0;
		for (p0 = 0; p0 < d.n; p0++) {
			if (d.rm[p0]) continue;    
			if (d.M[p0] > max_M)			max_M = d.M[p0];
			if (d.N[p0]-d.M[p0] > max_M)	max_M = d.N[p0]-d.M[p0];
			if (d.N[p0] > max_N)            max_N = d.N[p0];
		}    

		p1	  = new double[d.n];

		r0          = new double[max_N+1];
		//beta-binomial distribution
		p0_r0       = new double[max_M+1]; //sum_r log(p0+r0)
		l_r0        = new double[max_N+1]; //sum_r log(1+r0)    
		//theta
		l_r0_r      = new double[max_N+1]; //\sum_r r/(1+r0)
		p0_r0_r     = new double[max_M+1]; //\sum_r r/(p0+r0)
		//theta2
		l_r0_r_2    = new double[max_N+1]; //\sum_r [r/(1+r0)]^2
		p0_r0_r_2   = new double[max_M+1]; //\sum_r [r/(p0+r0)]^2       
	}//ASE

	//TREC variables
	int ifault = 0;
	double *vecR=NULL, phi=0, *mu_i=NULL, *tmu_i=NULL;
	if (p.TREC) {
		vecR  = new double[p.R];
		mu_i  = new double[d.n];
		tmu_i = new double[d.n];
	}//TREC


	if (1) {
        
		////////////////////////////////////////////////////
        // save variables that are common to individuals
		////////////////////////////////////////////////////
		
		if (p.ASE) {

			theta = v[p.AD];
		
			for (r = 1; r <= max_N; r++) r0[r] = (r-1)*theta;

			//beta-binomial distribution
			p0_r0[0] = 0; 	for (r = 1; r <= max_M; r++) p0_r0[r]  = p0_r0[r-1]  + log(pn + r0[r]);
			l_r0[0]  = 0;	for (r = 1; r <= max_N; r++)  l_r0[r]  =  l_r0[r-1]  + log( 1 + r0[r]);
			
			boundary = -100; //MAX(-pn/(max_M-1),-1.0/(max_N-1));

		}//ASE
        
		if (p.TREC) {
			if (p.disp_TREC) phi = v[p.TD]; 
		}

		d.n_new = 0.0;
		loglike = 0.0;

		for (p0 = 0; p0 < d.n; p0++) {
			if (d.rm[p0]) continue;    

			//exp_beta[d.n]

			exp_beta[p0] = 1.0;

			if (d.G[p0] >= 1) {
				Bvec(vecB, p.B, d.X[p0]);

				double beta_X = 0.0;
				for (j = 0; j < p.B; j++) {
					beta_X += v[j] * vecB[j];
				}
				exp_beta[p0] = exp(beta_X);
			}

            d.sum[p0] = 0.0;

			//ASE
			if (p.ASE){

				int index1 = d.M[p0];
				int index2 = d.N[p0] - d.M[p0];
				int index3 = d.N[p0];

				if (d.G[p0] == 1 && d.N[p0] > 0) { //two pairs

					//p1[d.n]
					p1[p0] = exp_beta[p0]/(1+exp_beta[p0]);

					for (int pp = 0; pp < d.nhap[p0]; pp++) {
 
						if (d.hap1[p0][pp] / pow1 == 0) {     
							index1 = d.N[p0] - d.M[p0];
							index2 = d.M[p0];
						}
						else {
							index1 = d.M[p0];
							index2 = d.N[p0] - d.M[p0];
						}

						if (index3 > 1) boundary = MAX(boundary, -1.0/(index3-1));
						if (index1 > 1) boundary = MAX(boundary, -p1[p0]/(index1-1));
						if (index2 > 1) boundary = MAX(boundary, -(1-p1[p0])/(index2-1));

						//beta-binomial distribution
						double p_r0 = 0;	for (r = 1; r <= index1; r++)   p_r0 = p_r0   + log(  p1[p0] + r0[r]);
						double l_p_r0 = 0;  for (r = 1; r <= index2; r++) l_p_r0 = l_p_r0 + log(1-p1[p0] + r0[r]);

						d.w[p0][pp] = p_r0 + l_p_r0 - l_r0[index3] + log(v[p.FREQ+p.h[d.hap1[p0][pp]]]*v[p.FREQ+p.h[d.hap2[p0][pp]]]);
					}//pp

					if (d.nhap[p0]==2) {
						double w1 = d.w[p0][0];
						double w2 = d.w[p0][1];

						d.w[p0][0] = 1.0/(1.0+exp(w2-w1));
						d.w[p0][1] = 1.0/(1.0+exp(w1-w2));

						d.sum[p0] = w1 + log(1+exp(w2-w1));
					}
					else {
						d.sum[p0] = d.w[p0][0];
						d.w[p0][0] = 1.0;
					}
				}
				else if (d.G[p0] != 1 && d.N[p0] > 0) { //l pair
					int pp = 0;
					d.w[p0][pp] = 1.0;
					d.sum[p0] = p0_r0[index1] + p0_r0[index2] - l_r0[index3] + log(v[p.FREQ+p.h[d.hap1[p0][pp]]]*v[p.FREQ+p.h[d.hap2[p0][pp]]]);

					if (index3 > 1) boundary = MAX(boundary, -1.0/(index3-1));
					if (index1 > 1) boundary = MAX(boundary, -pn/(index1-1));
					if (index2 > 1) boundary = MAX(boundary, -pn/(index2-1));

				}
				else if (d.N[p0] == 0) {                   //can be two pairs or l pair
					int pp;
					for (pp = 0; pp < d.nhap[p0]; pp++) {
						d.w[p0][pp] = v[p.FREQ+p.h[d.hap1[p0][pp]]]*v[p.FREQ+p.h[d.hap2[p0][pp]]];
						d.sum[p0] += d.w[p0][pp];
					}
					for (pp = 0; pp < d.nhap[p0]; pp++) d.w[p0][pp] /= d.sum[p0];
					d.sum[p0] = log(d.sum[p0]);
				}

			}//ASE
            
			//TReC
			if (p.TREC) {

				double exp_Z;
				if (d.G[p0] == 0)		exp_Z = 1;
				else if (d.G[p0] == 1) exp_Z = 0.5*(1+exp_beta[p0]);
				else						exp_Z = exp_beta[p0];

				Rvec(vecR, p.R, d.X[p0]);
				double exp_eta = 0.0;
				for (j = 0; j < p.R; j++) {
					exp_eta += v[p.REG+j] * vecR[j]; 
				}
				exp_eta = exp(exp_eta);

				mu_i[p0] = exp_eta * exp_Z;
				tmu_i[p0] = 0.5 * d.G[p0] * exp_eta * exp_beta[p0];

				double mu = mu_i[p0];
				double T = d.T[p0];

				if (p.disp_TREC && strcmp(p.disp_TREC_model, "NB")==0) {
					d.sum[p0] += lgamma(T+phi) - d.lgamma_int[(int)T] - (phi+T)*log(phi+mu) + T*log(mu);
				}
				else if (p.disp_TREC && strcmp(p.disp_TREC_model, "GP")==0) {
					d.sum[p0] += log(mu) + (T-1)*log(mu+T*phi) - mu - T*phi;
				}
				else if (!p.disp_TREC) {
					d.sum[p0] += T*log(mu) - mu;
				}
			}//TReC

            //likelihood
			loglike -= d.sum[p0];  
			d.n_new++;
            
		}//person

		//likelihood (TReC)
		if (p.TREC && p.disp_TREC && strcmp(p.disp_TREC_model, "NB")==0) {
			loglike -= d.n_new * (- lgamma(phi) + phi*log(phi));
		}
		
		loglike = -loglike;  //Likelihood


		////////////////////////////////////////////////////
		// evaluate efficient scores
		////////////////////////////////////////////////////


        //U, E initialation
        //U = l', E = -l"
		for (j = 0; j < p.nall; j++) {
			for (p0 = 0; p0 < d.n; p0++) Ui[p0][j] = 0;
			for (k = j; k < p.nall; k++) E[j+1][k+1] = 0;
		}
				

		//ASE
		if (p.ASE){

			l_r0_r[0] = 0;   //theta
			p0_r0_r[0] = 0; 
			l_r0_r_2[0] = 0; //theta2
			p0_r0_r_2[0] = 0;
			for (r = 1; r <= max_M; r++) {
				p0_r0_r[r]      = p0_r0_r[r-1]  + (r-1)/  (pn + r0[r]);   
				p0_r0_r_2[r]    = p0_r0_r_2[r-1]  + (r-1)*(r-1)/  ((pn+r0[r])*  (pn+r0[r]));
			}      
			for (r = 1; r <= max_N; r++) {
				l_r0_r[r]     = l_r0_r[r-1] + (r-1)/(1+r0[r]);
				l_r0_r_2[r]   = l_r0_r_2[r-1] + (r-1)*(r-1)/((1+r0[r])*(1+r0[r]));
			}


			Hi   = new double[p.FREQ+p.KK-1];
			Hikl = new double[p.FREQ+p.KK-1];

			for (p0 = 0; p0 < d.n; p0++) {
				if (d.rm[p0]) continue;    
            
				for (j = 0; j < p.FREQ+p.KK-1; j++) Hi[j] = 0;

				Bvec(vecB, p.B, d.X[p0]);

				//default
				int index1 = d.M[p0];
				int index2 = d.N[p0] - d.M[p0];
				int index3 = d.N[p0];
          
				for (int pp = 0; pp < d.nhap[p0]; pp++) {    
	
					for (j = 0; j < p.FREQ+p.KK-1; j++) Hikl[j] = 0.0;

					if (d.G[p0] == 1 && d.N[p0] > 0) {            //N>0, G=1
                
						if (d.hap1[p0][pp] / pow1 == 0) {     
							index1 = d.N[p0] - d.M[p0];
							index2 = d.M[p0];
						}
						else {
							index1 = d.M[p0];
							index2 = d.N[p0] - d.M[p0];
						}

						double p_r0_1 = 0;			//beta
						double l_p_r0_1 = 0;
						double p_r0_r = 0;			//theta
						double l_p_r0_r = 0;
						double p_r0_1_2 = 0;		//beta2
						double l_p_r0_1_2 = 0; 			
						double p_r0_r_d2 = 0;		//beta*theta
						double l_p_r0_r_d2 = 0;			
						double p_r0_r_2 = 0;		//theta2
						double l_p_r0_r_2 = 0;
					
						double p1_tmp = p1[p0];
						double l_p1_tmp = 1-p1_tmp;
						for (r = 1; r <= index1; r++) {
							p_r0_1       = p_r0_1   + 1/(p1_tmp + r0[r]);
							p_r0_r       = p_r0_r   + (r-1)/(p1_tmp + r0[r]);
							p_r0_1_2     = p_r0_1_2 + 1/((p1_tmp+r0[r])*(p1_tmp+r0[r]));
							p_r0_r_d2    = p_r0_r_d2+ (r-1)/((p1_tmp+r0[r])*(p1_tmp+r0[r]));
							p_r0_r_2     = p_r0_r_2 + (r-1)*(r-1)/((p1_tmp+r0[r])*(p1_tmp+r0[r]));
						} 
						for (r = 1; r <= index2; r++) {
							l_p_r0_1     = l_p_r0_1   + 1/(l_p1_tmp + r0[r]);
							l_p_r0_r     = l_p_r0_r   + (r-1)/(l_p1_tmp + r0[r]);
							l_p_r0_1_2   = l_p_r0_1_2 + 1/((l_p1_tmp+r0[r])*(l_p1_tmp+r0[r]));
							l_p_r0_r_d2  = l_p_r0_r_d2+ (r-1)/((l_p1_tmp+r0[r])*(l_p1_tmp+r0[r]));
							l_p_r0_r_2   = l_p_r0_r_2 + (r-1)*(r-1)/((l_p1_tmp+r0[r])*(l_p1_tmp+r0[r]));
						} 


						double exp_beta_tmp = exp_beta[p0];
						for (k = 0; k < p.B; k++) {
							Hikl[k]			+=		   ( p_r0_1    - l_p_r0_1) * p1_tmp * l_p1_tmp * vecB[k];
							Ui[p0][k]	+= d.w[p0][pp] * ( p_r0_1    - l_p_r0_1) * p1_tmp * l_p1_tmp * vecB[k];
							E[k+1][p.AD+1]	-= d.w[p0][pp] * (-p_r0_r_d2 + l_p_r0_r_d2) * p1_tmp * l_p1_tmp * vecB[k];
							for (int l = k; l < p.B; l++) {
								E[k+1][l+1] -= d.w[p0][pp] * (-p_r0_1_2  - l_p_r0_1_2) * p1_tmp*p1_tmp*l_p1_tmp*l_p1_tmp * vecB[k]*vecB[l]
								             + d.w[p0][pp] * ( p_r0_1    - l_p_r0_1) * exp_beta_tmp * (1-exp_beta_tmp)/((1+exp_beta_tmp)*(1+exp_beta_tmp)*(1+exp_beta_tmp)) 
											         * vecB[k]*vecB[l];
							}
						}
						Hikl[p.AD]		  +=		 ( p_r0_r    + l_p_r0_r - l_r0_r[index3]);
						Ui[p0][p.AD] += d.w[p0][pp] * ( p_r0_r    + l_p_r0_r - l_r0_r[index3]);
						E[p.AD+1][p.AD+1] -= d.w[p0][pp] * (-p_r0_r_2  - l_p_r0_r_2 + l_r0_r_2[index3]);
					} //N>0, G=1

					else if (d.G[p0] != 1 && d.N[p0] > 0) {      //N>0, G=0, 2
						Hikl[p.AD]		  +=		 (p0_r0_r[index1] + p0_r0_r[index2] - l_r0_r[index3]);
						Ui[p0][p.AD] += d.w[p0][pp] * (p0_r0_r[index1] + p0_r0_r[index2] - l_r0_r[index3]);
						E[p.AD+1][p.AD+1] -= d.w[p0][pp] * (-p0_r0_r_2[index1] - p0_r0_r_2[index2] + l_r0_r_2[index3]);
					}

					//\pi
					if (p.h[d.hap1[p0][pp]] != 0) {
						Hikl[p.FREQ+p.h[d.hap1[p0][pp]]-1]					+=     1/v[p.FREQ+p.h[d.hap1[p0][pp]]];
						Ui[p0][p.FREQ+p.h[d.hap1[p0][pp]]-1]				+= d.w[p0][pp]/v[p.FREQ+p.h[d.hap1[p0][pp]]];
						E[p.FREQ+p.h[d.hap1[p0][pp]]][p.FREQ+p.h[d.hap1[p0][pp]]]	+= d.w[p0][pp]/(v[p.FREQ+p.h[d.hap1[p0][pp]]]*v[p.FREQ+p.h[d.hap1[p0][pp]]]);
					}
					if (p.h[d.hap2[p0][pp]] != 0) {
						Hikl[p.FREQ+p.h[d.hap2[p0][pp]]-1]					+=      1/v[p.FREQ+p.h[d.hap2[p0][pp]]];
						Ui[p0][p.FREQ+p.h[d.hap2[p0][pp]]-1]				+= d.w[p0][pp]/v[p.FREQ+p.h[d.hap2[p0][pp]]];
						E[p.FREQ+p.h[d.hap2[p0][pp]]][p.FREQ+p.h[d.hap2[p0][pp]]]	+= d.w[p0][pp]/(v[p.FREQ+p.h[d.hap2[p0][pp]]]*v[p.FREQ+p.h[d.hap2[p0][pp]]]);
					}
					int counter = (p.h[d.hap1[p0][pp]]==0) + (p.h[d.hap2[p0][pp]]==0);
					if (counter > 0) {
						for (i = 0; i < p.KK-1; i++) {
							Hikl[p.FREQ+i]					-=       counter/v[p.FREQ];
							Ui[p0][p.FREQ+i]			-= d.w[p0][pp]*counter/v[p.FREQ];
							for (j = i; j < p.KK-1; j++) {
								E[p.FREQ+i+1][p.FREQ+j+1]	+= d.w[p0][pp]*counter/(v[p.FREQ]*v[p.FREQ]);
							}
						}
					}

					for (k = 0; k < p.FREQ+p.KK-1; k++) {
						Hi[k] += d.w[p0][pp] * Hikl[k];
						for (int l = k; l < p.FREQ+p.KK-1; l++) {
							E[k+1][l+1] -= d.w[p0][pp] * Hikl[k] * Hikl[l]; 
						}
					}

				}//(hk,hl)

				for (k = 0; k < p.FREQ+p.KK-1; k++) {
					for (int l = k; l < p.FREQ+p.KK-1; l++) {
						E[k+1][l+1] += Hi[k] * Hi[l]; 
					}
				} 

			}//person

			delete [] Hi;
			delete [] Hikl;

			//E_Louis_ext
			//if (d.ext) Louis_orig_ext(E, U, d, p, v); //need to get Ui instead of U

			if (saturate_fitted) {
				for (j = 0; j < p.B; j++) {
					for (p0 = 0; p0 < d.n; p0++)		Uei[p0][j] = Ui[p0][j];
					for (k = 0; k < p.nall; k++)	E12[j+1][k+1] = E[j+1][k+1]; 
				}
			}

			if (v[p.AD] < boundary) {
				*ll = 2;
				goto end;
			}


		}//ASE


		if (p.TREC) {


			double c;
            
            if (p.disp_TREC && strcmp(p.disp_TREC_model, "NB")==0) {
				c = - digama(phi, &ifault) + 1 + log(phi);
			}

			for (p0 = 0; p0 < d.n; p0++) {
				if (d.rm[p0]) continue;    
            
				Bvec(vecB, p.B, d.X[p0]);

				Rvec(vecR, p.R, d.X[p0]);

				double T = d.T[p0];
				double mu = mu_i[p0];
				double tmu = tmu_i[p0];

				if (p.disp_TREC && strcmp(p.disp_TREC_model, "NB")==0) {

					for (j = 0; j < p.B; j++) {
						Ui[p0][j] += (T/mu - (phi+T)/(phi+mu)) * tmu * vecB[j];
						for (k = j; k < p.B; k++) {
							E[j+1][k+1]+= ((T/(mu*mu) - (phi+T)/((phi+mu)*(phi+mu))) *tmu*tmu - (T/mu - (phi+T)/(phi+mu)) * tmu) * vecB[j]*vecB[k];
						}
						E[j+1][p.TD+1] += (1/(phi+mu) - (phi+T)/((phi+mu)*(phi+mu))) * tmu * vecB[j];
						for (k = 0; k < p.R; k++) {
							E[j+1][p.REG+k+1] += (T/(mu*mu) - (phi+T)/((phi+mu)*(phi+mu))) * tmu*mu*vecB[j]*vecR[k] - (T/mu - (phi+T)/(phi+mu)) * tmu*vecB[j]*vecR[k];
						}
					}
					for (j = 0; j < p.R; j++) {
						Ui[p0][p.REG+j] += (T/mu - (phi+T)/(phi+mu)) * mu * vecR[j];
						for (k = j; k < p.R; k++) {
							E[p.REG+j+1][p.REG+k+1] += (T/(mu*mu) - (phi+T)/((phi+mu)*(phi+mu))) *mu*mu*vecR[j]*vecR[k] - (T/mu - (phi+T)/(phi+mu)) *mu*vecR[j]*vecR[k];
						}
						E[p.TD+1][p.REG+j+1] += (1/(phi+mu) - (phi+T)/((phi+mu)*(phi+mu))) * mu*vecR[j];
					}

					Ui[p0][p.TD] += digama(T+phi, &ifault) - (phi+T)/(phi+mu) - log(phi+mu) + c;
					if (ifault == 1) {*ll = 4; goto end;}
					E[p.TD+1][p.TD+1] += -trigamma(T+phi, &ifault) + (mu-T)/((phi+mu)*(phi+mu)) + 1/(phi+mu);
					if (ifault == 1) {*ll = 4; goto end;}

				}//with dispersion
				if (p.disp_TREC && strcmp(p.disp_TREC_model, "GP")==0) {

					double mu_T_phi = mu+T*phi;
					for (j = 0; j < p.B; j++) {
						Ui[p0][j]           += (1.0/mu + (T-1)/mu_T_phi - 1) * tmu * vecB[j];
						for (k = j; k < p.B; k++) {
							E[j+1][k+1]+= ((1.0/(mu*mu) + (T-1)/mu_T_phi/mu_T_phi) *tmu*tmu - (1.0/mu + (T-1)/mu_T_phi - 1) *tmu) * vecB[j]*vecB[k];
						}
						E[j+1][p.TD+1] += T*(T-1)/mu_T_phi/mu_T_phi * tmu * vecB[j];
						for (k = 0; k < p.R; k++) {
							E[j+1][p.REG+k+1] += ((1.0/(mu*mu) + (T-1)/mu_T_phi/mu_T_phi) *tmu*mu - (1.0/mu + (T-1)/mu_T_phi - 1) *tmu)*vecB[j]*vecR[k];
						}
					}
					for (j = 0; j < p.R; j++) {
						Ui[p0][p.REG+j] += (1.0/mu + (T-1)/mu_T_phi - 1) * mu * vecR[j];
						for (k = j; k < p.R; k++) {
							E[p.REG+j+1][p.REG+k+1] += ((1.0/(mu*mu) + (T-1)/mu_T_phi/mu_T_phi) *mu*mu - (1.0/mu + (T-1)/mu_T_phi - 1) *mu)*vecR[j]*vecR[k];
						}
						E[p.TD+1][p.REG+j+1] += T*(T-1)/mu_T_phi/mu_T_phi * mu*vecR[j];
					}

					Ui[p0][p.TD] += T*(T-1)/mu_T_phi - T;
					E[p.TD+1][p.TD+1] += T*T*(T-1)/mu_T_phi/mu_T_phi;

				}//with dispersion
				else if (!p.disp_TREC) {


					for (j = 0; j < p.B; j++) {
						Ui[p0][j] += (T/mu - 1) * tmu * vecB[j];
						for (k = j; k < p.B; k++) {
							E[j+1][k+1]+= (T/(mu*mu) *tmu*tmu - (T/mu - 1) * tmu) * vecB[j]*vecB[k];
						}
						for (k = 0; k < p.R; k++) {
							E[j+1][p.REG+k+1] += T/(mu*mu) * tmu*mu*vecB[j]*vecR[k] - (T/mu - 1) * tmu*vecB[j]*vecR[k];
						}
					}
					for (j = 0; j < p.R; j++) {
						Ui[p0][p.REG+j] += (T/mu - 1) * mu * vecR[j];
						for (k = j; k < p.R; k++) {
							E[p.REG+j+1][p.REG+k+1] += T*vecR[j]*vecR[k] - (T/mu - 1) *mu*vecR[j]*vecR[k];
						}
					}

				}//no dispersion

			}//person

			if (p.disp_TREC && strcmp(p.disp_TREC_model, "NB")==0) {
				//common terms
				E[p.TD+1][p.TD+1] += d.n_new * (trigamma(phi, &ifault) - 1/phi );
				if (ifault == 1) {*ll = 4; goto end;}
			}

		}//TREC


	}//if (1)


	for (j = 1; j <= p.nall; j++) {
		for (k = 1; k < j; k++) {
			E[j][k] = E[k][j];
		}
	}


	//------------ efficient score -------------//


	if (saturate_fitted) {

		nui = p.FREQ+p.KK-1;

		if ((*ll = explosion_E(E, nui))) {*ll = 3; goto end;}
		invv(E, nui);
        if ((*ll = explosion_E(E, nui))) {*ll = 3; goto end;}

		//E12_22 = E12 %*% E22
		for (k = 1; k <= p.B; k++) {
			for (i = 1; i <= nui; i++) {
				E12_22[k][i] = 0.0;
				for (j = 1; j <= nui; j++) E12_22[k][i] += E12[k][j] * E[j][i];
			}
		}

		//Uei
		for (p0 = 0; p0 < d.n; p0++) {
			if (d.rm[p0]) continue;
			for (k = 0; k < p.B; k++) {
				for (j = 0; j < nui; j++) {
					Uei[p0][k] -= E12_22[k+1][j+1] * Ui[p0][j];
				}
			}
		}

	}//saturate


	//V: p.B = 1;
	*U = 0.0;
	*V = 0.0;
	for (j = 0; j < nui; j++) Usum[j] = 0.0;

	for (p0 = 0; p0 < d.n; p0++) {
		if (d.rm[p0]) continue;
		*U += Uei[p0][0];
		*V += Uei[p0][0] * Uei[p0][0];

		for (j = 0; j < nui; j++) Usum[j] += Ui[p0][j];
	}

	*V = sqrt(*V - (*U)*(*U)/d.n);

end:

	delete [] v; 
	delete [] Usum;
  	free_dmatrix(Ui,  0, d.n-1, 0, p.B-1);
	free_dmatrix(E,   1, p.nall, 1, p.nall);
	if (saturate_fitted) {
		free_dmatrix(E12_22, 1, p.B, 1, p.nall);
		free_dmatrix(E12,  1, p.B, 1, p.nall);
	}

    
	delete [] vecB;
	delete [] exp_beta;

	if (p.ASE) {
		delete [] p1;
		delete [] r0;
		delete [] p0_r0;
		delete [] l_r0;    
		delete [] l_r0_r;
		delete [] p0_r0_r;
		delete [] l_r0_r_2;
		delete [] p0_r0_r_2;  
	}

	if (p.TREC){
		delete [] vecR;
		delete [] mu_i;
		delete [] tmu_i;
	}

}//evaluate_effiScore_trans


//****************************************************************************80

double digama ( double x, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    DIGAMA calculates DIGAMMA ( X ) = d ( LOG ( GAMMA ( X ) ) ) / dX
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by Jose Bernardo.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Jose Bernardo,
//    Algorithm AS 103:
//    Psi ( Digamma ) Function,
//    Applied Statistics,
//    Volume 25, Number 3, 1976, pages 315-317.
//
//  Parameters:
//
//    Input, double X, the argument of the digamma function.
//    0 < X.
//
//    Output, int *IFAULT, error flag.
//    0, no error.
//    1, X <= 0.
//
//    Output, double DIGAMA, the value of the digamma function at X.
//
{
  double c = 8.5;
  double d1 = -0.5772156649;
  double r;
  double s = 0.00001;
  double s3 = 0.08333333333;
  double s4 = 0.0083333333333;
  double s5 = 0.003968253968;
  double value;
  double y;
//
//  Check the input.
//
  if ( x <= 0.0 )
  {
    value = 0.0;
    *ifault = 1;
    return value;
  }
//
//  Initialize.
//
  *ifault = 0;
  y = x;
  value = 0.0;
//
//  Use approximation if argument <= S.
//
  if ( y <= s )
  {
    value = d1 - 1.0 / y;
    return value;
  }
//
//  Reduce to DIGAMA(X + N) where (X + N) >= C.
//
  while ( y < c )
  {
    value = value - 1.0 / y;
    y = y + 1.0;
  }
//
//  Use Stirling's (actually de Moivre's) expansion if argument > C.
//
  r = 1.0 / y;
  value = value + log ( y ) - 0.5 * r;
  r = r * r;
  value = value - r * ( s3 - r * ( s4 - r * s5 ) );

  return value;
}
//****************************************************************************80

void psi_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    PSI_VALUES returns some values of the Psi or Digamma function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      PolyGamma[x]
//
//    or
//
//      Polygamma[0,x]
//
//    PSI(X) = d ln ( Gamma ( X ) ) / d X = Gamma'(X) / Gamma(X)
//
//    PSI(1) = -Euler's constant.
//
//    PSI(X+1) = PSI(X) + 1 / X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 11

  double fx_vec[N_MAX] = { 
     -0.5772156649015329E+00,  
     -0.4237549404110768E+00,  
     -0.2890398965921883E+00,  
     -0.1691908888667997E+00,  
     -0.6138454458511615E-01,  
      0.3648997397857652E-01,  
      0.1260474527734763E+00,  
      0.2085478748734940E+00,  
      0.2849914332938615E+00,  
      0.3561841611640597E+00,  
      0.4227843350984671E+00 };

  double x_vec[N_MAX] = { 
     1.0E+00,  
     1.1E+00,  
     1.2E+00,  
     1.3E+00,  
     1.4E+00,  
     1.5E+00,  
     1.6E+00,  
     1.7E+00,  
     1.8E+00,  
     1.9E+00,  
     2.0E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80


//****************************************************************************80

double trigamma ( double x, int *ifault )

//****************************************************************************80
//
//  Purpose:
//
//    TRIGAMMA calculates trigamma(x) = d**2 log(gamma(x)) / dx**2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 January 2008
//
//  Author:
//
//    Original FORTRAN77 version by BE Schneider.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    BE Schneider,
//    Algorithm AS 121:
//    Trigamma Function,
//    Applied Statistics,
//    Volume 27, Number 1, pages 97-99, 1978.
//
//  Parameters:
//
//    Input, double X, the argument of the trigamma function.
//    0 < X.
//
//    Output, int *IFAULT, error flag.
//    0, no error.
//    1, X <= 0.
//
//    Output, double TRIGAMMA, the value of the trigamma function at X.
//
{
  double a = 0.0001;
  double b = 5.0;
  double b2 =  0.1666666667;
  double b4 = -0.03333333333;
  double b6 =  0.02380952381;
  double b8 = -0.03333333333;
  double value;
  double y;
  double z;
//
//  Check the input.
//
  if ( x <= 0.0 )
  {
    *ifault = 1;
    value = 0.0;
    return value;
  }

  *ifault = 0;
  z = x;
//
//  Use small value approximation if X <= A.
//
  if ( x <= a )
  {
    value = 1.0 / x / x;
    return value;
  }
//
//  Increase argument to ( X + I ) >= B.
//
  value = 0.0;

  while ( z < b )
  {
    value = value + 1.0 / z / z;
    z = z + 1.0;
  }
//
//  Apply asymptotic formula if argument is B or greater.
//
  y = 1.0 / z / z;

  value = value + 0.5 *
      y + ( 1.0
    + y * ( b2
    + y * ( b4
    + y * ( b6
    + y *   b8 )))) / z;

  return value;
}
//****************************************************************************80

void trigamma_values ( int *n_data, double *x, double *fx )

//****************************************************************************80
//
//  Purpose:
//
//    TRIGAMMA_VALUES returns some values of the TriGamma function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      PolyGamma[1,x]
//
//    TriGamma(X) = d^2 ln ( Gamma ( X ) ) / d X^2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *X, the argument of the function.
//
//    Output, double *FX, the value of the function.
//
{
# define N_MAX 11

  double fx_vec[N_MAX] = {
    0.1644934066848226E+01,
    0.1433299150792759E+01,
    0.1267377205423779E+01,
    0.1134253434996619E+01,
    0.1025356590529597E+01,
    0.9348022005446793E+00,
    0.8584318931245799E+00,
    0.7932328301639984E+00,
    0.7369741375017002E+00,
    0.6879720582426356E+00,
    0.6449340668482264E+00 };

  double x_vec[N_MAX] = {
     1.0E+00,
     1.1E+00,
     1.2E+00,
     1.3E+00,
     1.4E+00,
     1.5E+00,
     1.6E+00,
     1.7E+00,
     1.8E+00,
     1.9E+00,
     2.0E+00 };

  if ( *n_data < 0 )
  {
    *n_data = 0;
  }

  *n_data = *n_data + 1;

  if ( N_MAX < *n_data )
  {
    *n_data = 0;
    *x = 0.0;
    *fx = 0.0;
  }
  else
  {
    *x = x_vec[*n_data-1];
    *fx = fx_vec[*n_data-1];
  }

  return;
# undef N_MAX
}

// lgamma.cpp -- log gamma function of real argument.
//      Algorithms and coefficient values from "Computation of Special
//      Functions", Zhang and Jin, John Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
//  Returns log(gamma) of real argument.
//  NOTE: Returns 1e308 if argument is 0 or negative.
//

#define EQTL_PI 3.1415926535897932384626433

double lgamma(double x)
{
    double x0,x2,xp,gl,gl0;
    int n,k;
    static double a[] = {
        8.333333333333333e-02,
       -2.777777777777778e-03,
        7.936507936507937e-04,
       -5.952380952380952e-04,
        8.417508417508418e-04,
       -1.917526917526918e-03,
        6.410256410256410e-03,
       -2.955065359477124e-02,
        1.796443723688307e-01,
       -1.39243221690590};
    
    x0 = x;
    if (x <= 0.0) return 1e308;
    else if ((x == 1.0) || (x == 2.0)) return 0.0;
    else if (x <= 7.0) {
        n = (int)(7-x);
        x0 = x+n;
    }
    x2 = 1.0/(x0*x0);
    xp = 2.0*EQTL_PI;
    gl0 = a[9];
    for (k=8;k>=0;k--) {
        gl0 = gl0*x2 + a[k];
    }
    gl = gl0/x0+0.5*log(xp)+(x0-0.5)*log(x0)-x0;
    if (x <= 7.0) {
        for (k=1;k<=n;k++) {
            gl -= log(x0-1.0);
            x0 -= 1.0;
        }
    }
    return gl;
}
       


int cmpr_double(const void *a, const void *b) { 
	int c;
	if (*(double*)a-*(double*)b < 0.0) c = -1;
	else if (*(double*)a-*(double*)b > 0.0) c = 1;
	else c = 0;

    return c;
}


void evaluate_effiScore_TREC_TRECASE(int*ll, int*r_n_test, double*r_stat_mar, DATA d, PARA p, double***Uei, double*vASE, double*vTREC) {


	if (p.TREC)				ll[0] = 0;
	if (p.TREC && p.ASE)	ll[1] = 0;

	
	int p0, i, j, k, r, nui=0;
    
	double loglike; 
	double*v     = new double[p.nall];
	double**Ui   = dmatrix(0, d.n-1, 0, p.nall-1);
	double**E    = dmatrix(1, p.nall, 1, p.nall);

	//Uei = U1i - E12E22^{-1}U2i
	double**E12_22	= dmatrix(1, p.B, 1, p.nall-p.B);
	double**E12		= dmatrix(1, p.B, 1, p.nall-p.B);
	double**E22		= dmatrix(1, p.nall-p.B, 1, p.nall-p.B);
	double***E11	= dmatrix3(0, 1, 1, p.B, 1, p.B);

	double*Usum = new double[2];
	double**Vsum = dmatrix(0, 1, 0, 1);
	
	//initialization

	for (j = 0; j < p.B; j++) v[j] = 0.0;								//beta: SNP-related reg
	if (p.ASE) {
		v[p.AD] = vASE[0];												//theta: ASE dispersion
		for (j = 0; j < p.KK; j++) v[p.FREQ+j] = p.pai[p.hap_freq[j]];	//pi: hap freq
	}
	if (p.TREC) {
		if (p.disp_TREC) v[p.TD] = vTREC[0];							//phi: TReC dispersion   
		for (j = 0; j < p.R;  j++) v[p.REG+j] = vTREC[p.disp_TREC+j];	//eta: other reg of TReC excluding SNPs
	}


	//ASE & TRECT variables
	double *vecB = new double[p.B];
	double *exp_beta = new double[d.n];

    //ASE variables
	int pow1, max_M, max_N;
	double *Hi = NULL, *Hikl = NULL;
    double pn = 0.5, theta, boundary;
	double *p1, *r0, *p0_r0, *l_r0, *l_r0_r, *p0_r0_r, *l_r0_r_2, *p0_r0_r_2;
    
	if (p.ASE) {

		pow1 = p.pow2[p.M-1];

		max_M = 0;
		max_N = 0;
		for (p0 = 0; p0 < d.n; p0++) {
			if (d.rm[p0]) continue;    
			if (d.M[p0] > max_M)			max_M = d.M[p0];
			if (d.N[p0]-d.M[p0] > max_M)	max_M = d.N[p0]-d.M[p0];
			if (d.N[p0] > max_N)            max_N = d.N[p0];
		}    

		p1	  = new double[d.n];

		r0          = new double[max_N+1];
		//beta-binomial distribution
		p0_r0       = new double[max_M+1]; //sum_r log(p0+r0)
		l_r0        = new double[max_N+1]; //sum_r log(1+r0)    
		//theta
		l_r0_r      = new double[max_N+1]; //\sum_r r/(1+r0)
		p0_r0_r     = new double[max_M+1]; //\sum_r r/(p0+r0)
		//theta2
		l_r0_r_2    = new double[max_N+1]; //\sum_r [r/(1+r0)]^2
		p0_r0_r_2   = new double[max_M+1]; //\sum_r [r/(p0+r0)]^2       
	}//ASE

	//TREC variables
	int ifault = 0;
	double *vecR=NULL, phi=0, *mu_i=NULL, *tmu_i=NULL;
	if (p.TREC) {
		vecR  = new double[p.R];
		mu_i  = new double[d.n];
		tmu_i = new double[d.n];
	}//TREC

	
	////////////////////////////////////////////////////
    // save variables that are common to individuals
	////////////////////////////////////////////////////
		
	if (p.ASE) {

		theta = v[p.AD];
		
		for (r = 1; r <= max_N; r++) r0[r] = (r-1)*theta;

		//beta-binomial distribution
		p0_r0[0] = 0; 	for (r = 1; r <= max_M; r++) p0_r0[r]  = p0_r0[r-1]  + log(pn + r0[r]);
		l_r0[0]  = 0;	for (r = 1; r <= max_N; r++)  l_r0[r]  =  l_r0[r-1]  + log( 1 + r0[r]);
			
		boundary = -100; //MAX(-pn/(max_M-1),-1.0/(max_N-1));

	}//ASE
        
	if (p.TREC) {
		if (p.disp_TREC) phi = v[p.TD]; 
	}

	d.n_new = 0.0;
	loglike = 0.0;

	for (p0 = 0; p0 < d.n; p0++) {
		if (d.rm[p0]) continue;    

		//exp_beta[d.n]

		exp_beta[p0] = 1.0;

		if (d.G[p0] >= 1) {
			Bvec(vecB, p.B, d.X[p0]);

			double beta_X = 0.0;
			for (j = 0; j < p.B; j++) {
				beta_X += v[j] * vecB[j];
			}
			exp_beta[p0] = exp(beta_X);
		}

        d.sum[p0] = 0.0;

		//ASE
		if (p.ASE){

			int index1 = d.M[p0];
			int index2 = d.N[p0] - d.M[p0];
			int index3 = d.N[p0];

			if (d.G[p0] == 1 && d.N[p0] > 0) { //two pairs

				//p1[d.n]
				p1[p0] = exp_beta[p0]/(1+exp_beta[p0]);

				for (int pp = 0; pp < d.nhap[p0]; pp++) {
 
					if (d.hap1[p0][pp] / pow1 == 0) {     
						index1 = d.N[p0] - d.M[p0];
						index2 = d.M[p0];
					}
					else {
						index1 = d.M[p0];
						index2 = d.N[p0] - d.M[p0];
					}

					if (index3 > 1) boundary = MAX(boundary, -1.0/(index3-1));
					if (index1 > 1) boundary = MAX(boundary, -p1[p0]/(index1-1));
					if (index2 > 1) boundary = MAX(boundary, -(1-p1[p0])/(index2-1));

					//beta-binomial distribution
					double p_r0 = 0;	for (r = 1; r <= index1; r++)   p_r0 = p_r0   + log(  p1[p0] + r0[r]);
					double l_p_r0 = 0;  for (r = 1; r <= index2; r++) l_p_r0 = l_p_r0 + log(1-p1[p0] + r0[r]);

					d.w[p0][pp] = p_r0 + l_p_r0 - l_r0[index3] + log(v[p.FREQ+p.h[d.hap1[p0][pp]]]*v[p.FREQ+p.h[d.hap2[p0][pp]]]);
				}//pp

				if (d.nhap[p0]==2) {
					double w1 = d.w[p0][0];
					double w2 = d.w[p0][1];

					d.w[p0][0] = 1.0/(1.0+exp(w2-w1));
					d.w[p0][1] = 1.0/(1.0+exp(w1-w2));

					d.sum[p0] = w1 + log(1+exp(w2-w1));
				}
				else {
					d.sum[p0] = d.w[p0][0];
					d.w[p0][0] = 1.0;
				}
			}
			else if (d.G[p0] != 1 && d.N[p0] > 0) { //l pair
				int pp = 0;
				d.w[p0][pp] = 1.0;
				d.sum[p0] = p0_r0[index1] + p0_r0[index2] - l_r0[index3] + log(v[p.FREQ+p.h[d.hap1[p0][pp]]]*v[p.FREQ+p.h[d.hap2[p0][pp]]]);

				if (index3 > 1) boundary = MAX(boundary, -1.0/(index3-1));
				if (index1 > 1) boundary = MAX(boundary, -pn/(index1-1));
				if (index2 > 1) boundary = MAX(boundary, -pn/(index2-1));

			}
			else if (d.N[p0] == 0) {                   //can be two pairs or l pair
				int pp;
				for (pp = 0; pp < d.nhap[p0]; pp++) {
					d.w[p0][pp] = v[p.FREQ+p.h[d.hap1[p0][pp]]]*v[p.FREQ+p.h[d.hap2[p0][pp]]];
					d.sum[p0] += d.w[p0][pp];
				}
				for (pp = 0; pp < d.nhap[p0]; pp++) d.w[p0][pp] /= d.sum[p0];
				d.sum[p0] = log(d.sum[p0]);
			}

		}//ASE
            
		//TReC
		if (p.TREC) {

			double exp_Z;
			if (d.G[p0] == 0)			exp_Z = 1;
			else if (d.G[p0] == 1)		exp_Z = 0.5*(1+exp_beta[p0]);
			else						exp_Z = exp_beta[p0];

			Rvec(vecR, p.R, d.X[p0]);
			double exp_eta = 0.0;
			for (j = 0; j < p.R; j++) {
				exp_eta += v[p.REG+j] * vecR[j]; 
			}
			exp_eta = exp(exp_eta);

			mu_i[p0] = exp_eta * exp_Z;
			tmu_i[p0] = 0.5 * d.G[p0] * exp_eta * exp_beta[p0];

			double mu = mu_i[p0];
			double T = d.T[p0];

			if (p.disp_TREC && strcmp(p.disp_TREC_model, "NB")==0) {
				d.sum[p0] += lgamma(T+phi) - d.lgamma_int[(int)T] - (phi+T)*log(phi+mu) + T*log(mu);
			}
			else if (p.disp_TREC && strcmp(p.disp_TREC_model, "GP")==0) {
				d.sum[p0] += log(mu) + (T-1)*log(mu+T*phi) - mu - T*phi;
			}
			else if (!p.disp_TREC) {
				d.sum[p0] += T*log(mu) - mu;
			}
		}//TReC

           //likelihood
		loglike -= d.sum[p0];  
		d.n_new++;
            
	}//person

	//likelihood (TReC)
	if (p.TREC && p.disp_TREC && strcmp(p.disp_TREC_model, "NB")==0) {
		loglike -= d.n_new * (- lgamma(phi) + phi*log(phi));
	}
		
	loglike = -loglike;  //Likelihood


	if (d.n_new < d.n) printf("d.n_new (%d) < d.n (%d) !!!\n", d.n_new, d.n);

	////////////////////////////////////////////////////
	// evaluate efficient scores
	////////////////////////////////////////////////////


    //U, E initialation
    //U = l', E = -l"
	for (j = 0; j < p.nall; j++) {
		for (p0 = 0; p0 < d.n; p0++) Ui[p0][j] = 0;
		for (k = j; k < p.nall; k++) E[j+1][k+1] = 0;
	}


	//=================================================================
	//	 TReC
	//=================================================================


	if (p.TREC) {

		double c;
            
        if (p.disp_TREC && strcmp(p.disp_TREC_model, "NB")==0) {
			c = - digama(phi, &ifault) + 1 + log(phi);
		}

		for (p0 = 0; p0 < d.n; p0++) {
			if (d.rm[p0]) continue;    
            
			Bvec(vecB, p.B, d.X[p0]);

			Rvec(vecR, p.R, d.X[p0]);

			double T = d.T[p0];
			double mu = mu_i[p0];
			double tmu = tmu_i[p0];

			if (p.disp_TREC && strcmp(p.disp_TREC_model, "NB")==0) {

				for (j = 0; j < p.B; j++) {
					Ui[p0][j] += (T/mu - (phi+T)/(phi+mu)) * tmu * vecB[j];
					for (k = j; k < p.B; k++) {
						E[j+1][k+1]+= ((T/(mu*mu) - (phi+T)/((phi+mu)*(phi+mu))) *tmu*tmu - (T/mu - (phi+T)/(phi+mu)) * tmu) * vecB[j]*vecB[k];
					}
					E[j+1][p.TD+1] += (1/(phi+mu) - (phi+T)/((phi+mu)*(phi+mu))) * tmu * vecB[j];
					for (k = 0; k < p.R; k++) {
						E[j+1][p.REG+k+1] += (T/(mu*mu) - (phi+T)/((phi+mu)*(phi+mu))) * tmu*mu*vecB[j]*vecR[k] - (T/mu - (phi+T)/(phi+mu)) * tmu*vecB[j]*vecR[k];
					}
				}
				for (j = 0; j < p.R; j++) {
					Ui[p0][p.REG+j] += (T/mu - (phi+T)/(phi+mu)) * mu * vecR[j];
					for (k = j; k < p.R; k++) {
						E[p.REG+j+1][p.REG+k+1] += (T/(mu*mu) - (phi+T)/((phi+mu)*(phi+mu))) *mu*mu*vecR[j]*vecR[k] - (T/mu - (phi+T)/(phi+mu)) *mu*vecR[j]*vecR[k];
					}
					E[p.TD+1][p.REG+j+1] += (1/(phi+mu) - (phi+T)/((phi+mu)*(phi+mu))) * mu*vecR[j];
				}

				Ui[p0][p.TD] += digama(T+phi, &ifault) - (phi+T)/(phi+mu) - log(phi+mu) + c;
				if (ifault == 1) {ll[0] = 4; goto end;}
				E[p.TD+1][p.TD+1] += -trigamma(T+phi, &ifault) + (mu-T)/((phi+mu)*(phi+mu)) + 1/(phi+mu);
				if (ifault == 1) {ll[0] = 4; goto end;}

			}//with dispersion
			if (p.disp_TREC && strcmp(p.disp_TREC_model, "GP")==0) {

				double mu_T_phi = mu+T*phi;
				for (j = 0; j < p.B; j++) {
					Ui[p0][j]           += (1.0/mu + (T-1)/mu_T_phi - 1) * tmu * vecB[j];
					for (k = j; k < p.B; k++) {
						E[j+1][k+1]+= ((1.0/(mu*mu) + (T-1)/mu_T_phi/mu_T_phi) *tmu*tmu - (1.0/mu + (T-1)/mu_T_phi - 1) *tmu) * vecB[j]*vecB[k];
					}
					E[j+1][p.TD+1] += T*(T-1)/mu_T_phi/mu_T_phi * tmu * vecB[j];
					for (k = 0; k < p.R; k++) {
						E[j+1][p.REG+k+1] += ((1.0/(mu*mu) + (T-1)/mu_T_phi/mu_T_phi) *tmu*mu - (1.0/mu + (T-1)/mu_T_phi - 1) *tmu)*vecB[j]*vecR[k];
					}
				}
				for (j = 0; j < p.R; j++) {
					Ui[p0][p.REG+j] += (1.0/mu + (T-1)/mu_T_phi - 1) * mu * vecR[j];
					for (k = j; k < p.R; k++) {
						E[p.REG+j+1][p.REG+k+1] += ((1.0/(mu*mu) + (T-1)/mu_T_phi/mu_T_phi) *mu*mu - (1.0/mu + (T-1)/mu_T_phi - 1) *mu)*vecR[j]*vecR[k];
					}
					E[p.TD+1][p.REG+j+1] += T*(T-1)/mu_T_phi/mu_T_phi * mu*vecR[j];
				}

				Ui[p0][p.TD] += T*(T-1)/mu_T_phi - T;
				E[p.TD+1][p.TD+1] += T*T*(T-1)/mu_T_phi/mu_T_phi;

			}//with dispersion
			else if (!p.disp_TREC) {


				for (j = 0; j < p.B; j++) {
					Ui[p0][j] += (T/mu - 1) * tmu * vecB[j];
					for (k = j; k < p.B; k++) {
						E[j+1][k+1]+= (T/(mu*mu) *tmu*tmu - (T/mu - 1) * tmu) * vecB[j]*vecB[k];
					}
					for (k = 0; k < p.R; k++) {
						E[j+1][p.REG+k+1] += T/(mu*mu) * tmu*mu*vecB[j]*vecR[k] - (T/mu - 1) * tmu*vecB[j]*vecR[k];
					}
				}
				for (j = 0; j < p.R; j++) {
					Ui[p0][p.REG+j] += (T/mu - 1) * mu * vecR[j];
					for (k = j; k < p.R; k++) {
						E[p.REG+j+1][p.REG+k+1] += T*vecR[j]*vecR[k] - (T/mu - 1) *mu*vecR[j]*vecR[k];
					}
				}

			}//no dispersion

		}//person

		if (p.disp_TREC && strcmp(p.disp_TREC_model, "NB")==0) {
			//common terms
			E[p.TD+1][p.TD+1] += d.n_new * (trigamma(phi, &ifault) - 1/phi );
			if (ifault == 1) {ll[0] = 4; goto end;}
		}


		//------------ TReC: efficient score -------------//


		for (j = 1; j <= p.nfix; j++) {
			for (k = 1; k < j; k++) {
				E[j][k] = E[k][j];
			}
		}

		nui = p.disp_TREC + p.R;

		//E12, E22
		for (j = 1; j <= nui; j++) {
			for (k = 1; k <= nui; k++) E22[j][k] = E[p.B+p.ASE+j][p.B+p.ASE+k];
			for (k = 1; k <= p.B; k++) E12[k][j] = E[k][p.B+p.ASE+j];
		}

		if ((ll[0] = explosion_E(E22, nui))) {ll[0] = 3; goto end;}
		invv(E22, nui);
		if ((ll[0] = explosion_E(E22, nui))) {ll[0] = 3; goto end;}


		//E12_22 = E12 %*% E22
		for (k = 1; k <= p.B; k++) {
			for (i = 1; i <= nui; i++) {
				E12_22[k][i] = 0.0;
				for (j = 1; j <= nui; j++) E12_22[k][i] += E12[k][j] * E22[j][i];
			}
		}

		//Uei
		for (p0 = 0; p0 < d.n; p0++) {
			if (d.rm[p0]) {
				for (k = 0; k < p.B; k++) Uei[0][p0][k] = 0.0;
				continue;
			}
			for (k = 0; k < p.B; k++) {
				Uei[0][p0][k] = Ui[p0][k];
				for (j = 0; j < nui; j++) Uei[0][p0][k] -= E12_22[k+1][j+1] * Ui[p0][p.B+p.ASE+j];
			}
		}

		//E11 Fisher
		for (k = 1; k <= p.B; k++) {
			for (i = 1; i <= p.B; i++) {
				E11[0][k][i] = E[k][i];
				for (j = 1; j <= nui; j++) E11[0][k][i] -= E12_22[k][j] * E12[i][j];
			}
		}


	}//TREC

		
	//=================================================================
	//   ASE
	//=================================================================

	if (p.ASE){

		l_r0_r[0] = 0;   //theta
		p0_r0_r[0] = 0; 
		l_r0_r_2[0] = 0; //theta2
		p0_r0_r_2[0] = 0;
		for (r = 1; r <= max_M; r++) {
			p0_r0_r[r]      = p0_r0_r[r-1]  + (r-1)/  (pn + r0[r]);   
			p0_r0_r_2[r]    = p0_r0_r_2[r-1]  + (r-1)*(r-1)/  ((pn+r0[r])*  (pn+r0[r]));
		}      
		for (r = 1; r <= max_N; r++) {
			l_r0_r[r]     = l_r0_r[r-1] + (r-1)/(1+r0[r]);
			l_r0_r_2[r]   = l_r0_r_2[r-1] + (r-1)*(r-1)/((1+r0[r])*(1+r0[r]));
		}


		Hi   = new double[p.FREQ+p.KK-1];
		Hikl = new double[p.FREQ+p.KK-1];

		for (p0 = 0; p0 < d.n; p0++) {
			if (d.rm[p0]) continue;    
            
			for (j = 0; j < p.FREQ+p.KK-1; j++) Hi[j] = 0;

			Bvec(vecB, p.B, d.X[p0]);

			//default
			int index1 = d.M[p0];
			int index2 = d.N[p0] - d.M[p0];
			int index3 = d.N[p0];
          
			for (int pp = 0; pp < d.nhap[p0]; pp++) {    
	
				for (j = 0; j < p.FREQ+p.KK-1; j++) Hikl[j] = 0.0;

				if (d.G[p0] == 1 && d.N[p0] > 0) {            //N>0, G=1
                
					if (d.hap1[p0][pp] / pow1 == 0) {     
						index1 = d.N[p0] - d.M[p0];
						index2 = d.M[p0];
					}
					else {
						index1 = d.M[p0];
						index2 = d.N[p0] - d.M[p0];
					}

					double p_r0_1 = 0;			//beta
					double l_p_r0_1 = 0;
					double p_r0_r = 0;			//theta
					double l_p_r0_r = 0;
					double p_r0_1_2 = 0;		//beta2
					double l_p_r0_1_2 = 0; 			
					double p_r0_r_d2 = 0;		//beta*theta
					double l_p_r0_r_d2 = 0;			
					double p_r0_r_2 = 0;		//theta2
					double l_p_r0_r_2 = 0;
					
					double p1_tmp = p1[p0];
					double l_p1_tmp = 1-p1_tmp;
					for (r = 1; r <= index1; r++) {
						p_r0_1       = p_r0_1   + 1/(p1_tmp + r0[r]);
						p_r0_r       = p_r0_r   + (r-1)/(p1_tmp + r0[r]);
						p_r0_1_2     = p_r0_1_2 + 1/((p1_tmp+r0[r])*(p1_tmp+r0[r]));
						p_r0_r_d2    = p_r0_r_d2+ (r-1)/((p1_tmp+r0[r])*(p1_tmp+r0[r]));
						p_r0_r_2     = p_r0_r_2 + (r-1)*(r-1)/((p1_tmp+r0[r])*(p1_tmp+r0[r]));
					} 
					for (r = 1; r <= index2; r++) {
						l_p_r0_1     = l_p_r0_1   + 1/(l_p1_tmp + r0[r]);
						l_p_r0_r     = l_p_r0_r   + (r-1)/(l_p1_tmp + r0[r]);
						l_p_r0_1_2   = l_p_r0_1_2 + 1/((l_p1_tmp+r0[r])*(l_p1_tmp+r0[r]));
						l_p_r0_r_d2  = l_p_r0_r_d2+ (r-1)/((l_p1_tmp+r0[r])*(l_p1_tmp+r0[r]));
						l_p_r0_r_2   = l_p_r0_r_2 + (r-1)*(r-1)/((l_p1_tmp+r0[r])*(l_p1_tmp+r0[r]));
					} 


					double exp_beta_tmp = exp_beta[p0];
					for (k = 0; k < p.B; k++) {
						Hikl[k]			+=		   ( p_r0_1    - l_p_r0_1) * p1_tmp * l_p1_tmp * vecB[k];
						Ui[p0][k]	+= d.w[p0][pp] * ( p_r0_1    - l_p_r0_1) * p1_tmp * l_p1_tmp * vecB[k];
						E[k+1][p.AD+1]	-= d.w[p0][pp] * (-p_r0_r_d2 + l_p_r0_r_d2) * p1_tmp * l_p1_tmp * vecB[k];
						for (int l = k; l < p.B; l++) {
							E[k+1][l+1] -= d.w[p0][pp] * (-p_r0_1_2  - l_p_r0_1_2) * p1_tmp*p1_tmp*l_p1_tmp*l_p1_tmp * vecB[k]*vecB[l]
							             + d.w[p0][pp] * ( p_r0_1    - l_p_r0_1) * exp_beta_tmp * (1-exp_beta_tmp)/((1+exp_beta_tmp)*(1+exp_beta_tmp)*(1+exp_beta_tmp)) 
								         * vecB[k]*vecB[l];
						}
					}
					Hikl[p.AD]		  +=		 ( p_r0_r    + l_p_r0_r - l_r0_r[index3]);
					Ui[p0][p.AD] += d.w[p0][pp] * ( p_r0_r    + l_p_r0_r - l_r0_r[index3]);
					E[p.AD+1][p.AD+1] -= d.w[p0][pp] * (-p_r0_r_2  - l_p_r0_r_2 + l_r0_r_2[index3]);
				} //N>0, G=1

				else if (d.G[p0] != 1 && d.N[p0] > 0) {      //N>0, G=0, 2
					Hikl[p.AD]		  +=		 (p0_r0_r[index1] + p0_r0_r[index2] - l_r0_r[index3]);
					Ui[p0][p.AD] += d.w[p0][pp] * (p0_r0_r[index1] + p0_r0_r[index2] - l_r0_r[index3]);
					E[p.AD+1][p.AD+1] -= d.w[p0][pp] * (-p0_r0_r_2[index1] - p0_r0_r_2[index2] + l_r0_r_2[index3]);
				}

				//\pi
				if (p.h[d.hap1[p0][pp]] != 0) {
					Hikl[p.FREQ+p.h[d.hap1[p0][pp]]-1]					+=     1/v[p.FREQ+p.h[d.hap1[p0][pp]]];
					Ui[p0][p.FREQ+p.h[d.hap1[p0][pp]]-1]				+= d.w[p0][pp]/v[p.FREQ+p.h[d.hap1[p0][pp]]];
					E[p.FREQ+p.h[d.hap1[p0][pp]]][p.FREQ+p.h[d.hap1[p0][pp]]]	+= d.w[p0][pp]/(v[p.FREQ+p.h[d.hap1[p0][pp]]]*v[p.FREQ+p.h[d.hap1[p0][pp]]]);
				}
				if (p.h[d.hap2[p0][pp]] != 0) {
					Hikl[p.FREQ+p.h[d.hap2[p0][pp]]-1]					+=      1/v[p.FREQ+p.h[d.hap2[p0][pp]]];
					Ui[p0][p.FREQ+p.h[d.hap2[p0][pp]]-1]				+= d.w[p0][pp]/v[p.FREQ+p.h[d.hap2[p0][pp]]];
					E[p.FREQ+p.h[d.hap2[p0][pp]]][p.FREQ+p.h[d.hap2[p0][pp]]]	+= d.w[p0][pp]/(v[p.FREQ+p.h[d.hap2[p0][pp]]]*v[p.FREQ+p.h[d.hap2[p0][pp]]]);
				}
				int counter = (p.h[d.hap1[p0][pp]]==0) + (p.h[d.hap2[p0][pp]]==0);
				if (counter > 0) {
					for (i = 0; i < p.KK-1; i++) {
						Hikl[p.FREQ+i]					-=       counter/v[p.FREQ];
						Ui[p0][p.FREQ+i]			-= d.w[p0][pp]*counter/v[p.FREQ];
						for (j = i; j < p.KK-1; j++) {
							E[p.FREQ+i+1][p.FREQ+j+1]	+= d.w[p0][pp]*counter/(v[p.FREQ]*v[p.FREQ]);
						}
					}
				}

				for (k = 0; k < p.FREQ+p.KK-1; k++) {
					Hi[k] += d.w[p0][pp] * Hikl[k];
					for (int l = k; l < p.FREQ+p.KK-1; l++) {
						E[k+1][l+1] -= d.w[p0][pp] * Hikl[k] * Hikl[l]; 
					}
				}

			}//(hk,hl)

			for (k = 0; k < p.FREQ+p.KK-1; k++) {
				for (int l = k; l < p.FREQ+p.KK-1; l++) {
					E[k+1][l+1] += Hi[k] * Hi[l]; 
				}
			} 

		}//person

		delete [] Hi;
		delete [] Hikl;

		//E_Louis_ext
		//if (d.ext) Louis_orig_ext(E, U, d, p, v); //need to get Ui instead of U

		if (v[p.AD] < boundary) {
			ll[1] = 2;
			goto end;
		}


		//------------ TReCASE: efficient score -------------//


		for (j = 1; j <= p.nall; j++) {
			for (k = 1; k < j; k++) {
				E[j][k] = E[k][j];
			}
		}

		nui = p.FREQ+p.KK-1-p.B;


		//E12, E22
		for (j = 1; j <= nui; j++) {
			for (k = 1; k <= nui; k++) E22[j][k] = E[p.B+j][p.B+k];
			for (k = 1; k <= p.B; k++) E12[k][j] = E[k][p.B+j];
		}

		if ((ll[1] = explosion_E(E22, nui))) {ll[1] = 3; goto end;}
		invv(E22, nui);
		if ((ll[1] = explosion_E(E22, nui))) {ll[1] = 3; goto end;}


		//E12_22 = E12 %*% E22
		for (k = 1; k <= p.B; k++) {
			for (i = 1; i <= nui; i++) {
				E12_22[k][i] = 0.0;
				for (j = 1; j <= nui; j++) E12_22[k][i] += E12[k][j] * E22[j][i];
			}
		}

		//Uei
		for (p0 = 0; p0 < d.n; p0++) {
			if (d.rm[p0]) {
				for (k = 0; k < p.B; k++) Uei[1][p0][k] = 0.0;
				continue;
			}
			for (k = 0; k < p.B; k++) {
				Uei[1][p0][k] = Ui[p0][k];
				for (j = 0; j < nui; j++) {
					Uei[1][p0][k] -= E12_22[k+1][j+1] * Ui[p0][p.B+j];
				}
			}
		}

		//E11 Fisher
		for (k = 1; k <= p.B; k++) {
			for (i = 1; i <= p.B; i++) {
				E11[1][k][i] = E[k][i];
				for (j = 1; j <= nui; j++) {
					E11[1][k][i] -= E12_22[k][j] * E12[i][j];
				}
			}
		}


	}//ASE


	//===================
	//get pvalue_orig
	//===================

	//p.B should be 1


	*r_n_test = 0;
	if (ll[0] == 0) (*r_n_test)++;
	if (ll[1] == 0) (*r_n_test)++;


	for (k = 0; k < (*r_n_test); k++) {

		double Umax = 0.0, Umin = 0.0;
				
		Usum[k] = 0.0;
		Vsum[k][k] = 0.0;
		for (i = 0; i < d.n; i++) {
			Usum[k] += Uei[k][i][0];
			Vsum[k][k] += Uei[k][i][0]*Uei[k][i][0];

			//if (Uei[k][i][0] < Umin)		Umin = Uei[k][i][0];
			//else if (Uei[k][i][0] > Umax)	Umax = Uei[k][i][0];
		}

		//Usum[k] += - Umin - Umax;
		//Vsum[k][k] += - Umin*Umin - Umax*Umax;

		r_stat_mar[k] = fabs(Usum[k] / sqrt(Vsum[k][k]-Usum[k]*Usum[k]/(d.n))); //debugged: Uei tends to be big for heterozygotes; Vsum[k][k] tends to be big; r_stat_mar[1] < r_stat_mar[0]

	}//k


end:

	delete [] v; 
  	free_dmatrix(Ui,  0, d.n-1, 0, p.B-1);
	free_dmatrix(E,   1, p.nall, 1, p.nall);
	free_dmatrix(E12_22, 1, p.B, 1, p.nall-p.B);
	free_dmatrix(E12,  1, p.B, 1, p.nall-p.B);
	free_dmatrix(E22,  1, p.nall-p.B, 1, p.nall-p.B);
	free_dmatrix3(E11, 0, 1, 1, p.B, 1, p.B);

	delete [] Usum;
	free_dmatrix(Vsum, 0, 1, 0, 1);
    
	delete [] vecB;
	delete [] exp_beta;

	if (p.ASE) {
		delete [] p1;
		delete [] r0;
		delete [] p0_r0;
		delete [] l_r0;    
		delete [] l_r0_r;
		delete [] p0_r0_r;
		delete [] l_r0_r_2;
		delete [] p0_r0_r_2;  
	}

	if (p.TREC){
		delete [] vecR;
		delete [] mu_i;
		delete [] tmu_i;
	}

}//evaluate_effiScore_TREC_TRECASE

int binsearch_int_2(int key, int *str, int max, int left) {
	//find the approximate position

	int position;	
	int begin = 0; 
	int end = max - 1;
	int cond = 0;

	while(begin <= end) {
		position = (begin + end) / 2;
		if((cond = (str[position] - key)) == 0)
			return position;
		else if (cond < 0)
			begin = position + 1;
		else
			end = position - 1;
	}

	if (left == 1) return (end < 0 ? 0 : end);
	else return (begin > max-1 ? max-1 : begin);
}


double kiss(unsigned long *seed)
{
    double r;
    //static long int seed[0]=seed[0];//seeds for simulation
    //static long int seed[1]=seed[1];//seeds for simulation
    //static long int seed[2]=seed[2];//seeds for simulation
    
    seed[0] = 171 * (seed[0] % 177) - 2 * (seed[0] / 177);
    seed[1] = 172 * (seed[1] % 176) - 35 * (seed[1] / 176);
    seed[2] = 170 * (seed[2] % 178) - 63 * (seed[2] / 178);
    
    if (seed[0] < 0) seed[0] = seed[0] + 30269;
    if (seed[1] < 0) seed[1] = seed[1] + 30307;
    if (seed[2] < 0) seed[2] = seed[2] + 30323;
    
    r = fmod(double(seed[0]) / 30269. + double(seed[1]) / 30307. +  double(seed[2]) / 30323., 1.0);
    
    return r;
	
}

int runiform_n_rand(int n,double rand)
{
  int i;

  i = (int)floor(n*rand);
  if (i == n)
    return(i-1);
  else 
    return(i);
}


void permute_sample_int(int *v,int len,double*rand, int*rm) 
{
	int i,j;
	int x;
  
	for (i=len;i>0;i--) {
		if (rm[i-1]) {
			continue;
		}

		j = runiform_n_rand(i,rand[i-1]);

		if (rm[j]==1) {	
			int jtmp = j;
			while (jtmp>=1) {
				jtmp--;
				if (rm[jtmp]==0) break;
			}
			if (rm[jtmp]==1) {
				jtmp = j;
				while (jtmp<=i-2) {
					jtmp++;
					if (rm[jtmp]==0) break;
				}
			}
			j = jtmp;
		}

		x = v[j];
		v[j] = v[i-1];
		v[i-1] = x;
	}
}


void permutation (double stat_perm_pmin[2][5000], int nPERM, DATA d, PARA p, double*vTREC, double*vASE) {

	int p0, i, j, k, l, r, nui=0;
    
	//----------------------------------
	// fixed over permutation
	//----------------------------------

	int ll[2] = {9, 9};
	if (p.TREC)				ll[0] = 0;
	if (p.TREC && p.ASE)	ll[1] = 0;


	double*rand1 = new double[d.n];
	double*rand2 = new double[d.n];


	int dim;
	double *Usum, **Vsum;

	int *permi	 = new int[d.n];
	double*v     = new double[p.nall];				
	double***E22 = dmatrix3(0, 1, 1, p.nall-p.B, 1, p.nall-p.B);
	double***E11 = dmatrix3(0, 1, 1, p.B, 1, p.B);

	//initialization

	for (j = 0; j < p.B; j++) v[j] = 0.0;								//beta: SNP-related reg
	if (p.ASE) {
		v[p.AD] = vASE[0];												//theta: ASE dispersion
		for (j = 0; j < p.KK; j++) v[p.FREQ+j] = p.pai[p.hap_freq[j]];	//pi: hap freq
	}
	if (p.TREC) {
		if (p.disp_TREC) v[p.TD] = vTREC[0];							//phi: TReC dispersion   
		for (j = 0; j < p.R;  j++) v[p.REG+j] = vTREC[p.disp_TREC+j];	//eta: other reg of TReC excluding SNPs
	}

	//ASE & TRECT variables
	double *vecB = new double[p.B];

    //ASE variables
	int pow1, max_M, max_N;
    double pn = 0.5, theta, boundary;
	double *r0, *p0_r0, *l_r0, *l_r0_r, *p0_r0_r, *l_r0_r_2, *p0_r0_r_2, *p0_r0_1, *p0_r0_1_2, *p0_r0_r_d2;
	double *Hi = NULL, ***Hikl = NULL;

	if (p.ASE) {

		pow1 = p.pow2[p.M-1];

		max_M = 0;
		max_N = 0;
		for (p0 = 0; p0 < d.n; p0++) {
			if (d.rm[p0]) continue;    
			if (d.M[p0] > max_M)			max_M = d.M[p0];
			if (d.N[p0]-d.M[p0] > max_M)	max_M = d.N[p0]-d.M[p0];
			if (d.N[p0] > max_N)            max_N = d.N[p0];
		}    

		r0          = new double[max_N+1];
		//beta-binomial distribution
		p0_r0       = new double[max_M+1]; //sum_r log(p0+r0)
		l_r0        = new double[max_N+1]; //sum_r log(1+r0)    
		//theta
		l_r0_r      = new double[max_N+1]; //\sum_r r/(1+r0)
		p0_r0_r     = new double[max_M+1]; //\sum_r r/(p0+r0)
		//theta2
		l_r0_r_2    = new double[max_N+1]; //\sum_r [r/(1+r0)]^2
		p0_r0_r_2   = new double[max_M+1]; //\sum_r [r/(p0+r0)]^2    
		//beta
		p0_r0_1     = new double[max_M+1];
		//beta2
		p0_r0_1_2   = new double[max_M+1];
		//beta*theta
		p0_r0_r_d2  = new double[max_M+1];


		Hi   = new double[p.FREQ+p.KK-1];
		Hikl = dmatrix3(0, d.n-1, 0, 1, 0, p.FREQ+p.KK-2);
	
	}//ASE

	//TREC variables
	int ifault = 0;
	double *vecR=NULL, phi=0, *mu_i=NULL, *exp_eta;
	if (p.TREC) {
		vecR  = new double[p.R];
		mu_i  = new double[d.n];
		exp_eta = new double[d.n];
	}//TREC

	//----------------------------------
	// update per permutation
	//----------------------------------

	int iperm;

	double***Uei= dmatrix3(0, 1, 0, d.n-1, 0, p.B-1);
	double**Ve  = dmatrix(1, p.B, 1, p.B);
	double *U   = new double[p.nall];
	double**Ui  = dmatrix(0, d.n-1, 0, p.nall-1);
	double**E   = dmatrix(1, p.nall, 1, p.nall);
	double**E12_22 = dmatrix(1, p.B, 1, p.nall-p.B);

    
	//TREC variables
	double c;


        
	////////////////////////////////////////////////////
    // save variables 
	//--------------------------------------------------
	// theta, r0, p0_r0, l_r0, boundary, d.w[p0][pp]
	// d.n_new
	// phi, exp_eta, mu_i
	////////////////////////////////////////////////////
		
	if (p.ASE) {

		theta = v[p.AD];
		
		for (r = 1; r <= max_N; r++) r0[r] = (r-1)*theta;

		//beta-binomial distribution
		p0_r0[0] = 0; 	for (r = 1; r <= max_M; r++) p0_r0[r]  = p0_r0[r-1]  + log(pn + r0[r]);
		l_r0[0]  = 0;	for (r = 1; r <= max_N; r++)  l_r0[r]  =  l_r0[r-1]  + log( 1 + r0[r]);
			
		boundary = -100; //MAX(-pn/(max_M-1),-1.0/(max_N-1));

	}//ASE
        
	if (p.TREC) {
		if (p.disp_TREC) phi = v[p.TD]; 
	}

	d.n_new = 0.0;
	for (p0 = 0; p0 < d.n; p0++) {
		if (d.rm[p0]) continue;    

		//ASE
		if (p.ASE){

			if (d.N[p0] > 1)		boundary = MAX(boundary, -1.0/(d.N[p0]-1));
			if (d.M[p0] > 1)		boundary = MAX(boundary, -pn/(d.M[p0]-1));
			if (d.N[p0]-d.M[p0]>1)	boundary = MAX(boundary, -pn/(d.N[p0] - d.M[p0]-1));

			int pp;
			double sum = 0.0;
			for (pp = 0; pp < d.nhap[p0]; pp++) {
				d.w[p0][pp] = v[p.FREQ+p.h[d.hap1[p0][pp]]]*v[p.FREQ+p.h[d.hap2[p0][pp]]];
				sum += d.w[p0][pp];
			}
			for (pp = 0; pp < d.nhap[p0]; pp++) d.w[p0][pp] /= sum;

		}//ASE
            
		//TReC
		if (p.TREC) {

			Rvec(vecR, p.R, d.X[p0]);
			exp_eta[p0] = 0.0;
			for (j = 0; j < p.R; j++) {
				exp_eta[p0] += v[p.REG+j] * vecR[j]; 
			}
			exp_eta[p0] = exp(exp_eta[p0]);

			mu_i[p0] = exp_eta[p0];

		}//TReC

        //likelihood
		d.n_new++;
            
	}//person

	if (d.n_new < d.n) printf("d.n_new (%d) < d.n (%d) !!!\n", d.n_new, d.n);


	////////////////////////////////////////////////////
    // evaluating Ui[nui], E[nui][nui]
	////////////////////////////////////////////////////


	//U, E initialation
    //U = l', E = -l"
	for (j = p.B; j < p.nall; j++) {
		for (p0 = 0; p0 < d.n; p0++) Ui[p0][j] = 0;
		for (k = j; k < p.nall; k++) E[j+1][k+1] = 0;
	}
		

	if (p.TREC) {

            
        if (p.disp_TREC && strcmp(p.disp_TREC_model, "NB")==0) {
			c = - digama(phi, &ifault) + 1 + log(phi);
		}

		for (p0 = 0; p0 < d.n; p0++) {
			if (d.rm[p0]) continue;    
            
			Rvec(vecR, p.R, d.X[p0]);

			double T = d.T[p0];
			double mu = mu_i[p0];

			if (p.disp_TREC && strcmp(p.disp_TREC_model, "NB")==0) {

				for (j = 0; j < p.R; j++) {
					Ui[p0][p.REG+j] += (T/mu - (phi+T)/(phi+mu)) * mu * vecR[j];
					for (k = j; k < p.R; k++) {
						E[p.REG+j+1][p.REG+k+1] += (T/(mu*mu) - (phi+T)/((phi+mu)*(phi+mu))) *mu*mu*vecR[j]*vecR[k] - (T/mu - (phi+T)/(phi+mu)) *mu*vecR[j]*vecR[k];
					}
					E[p.TD+1][p.REG+j+1] += (1/(phi+mu) - (phi+T)/((phi+mu)*(phi+mu))) * mu*vecR[j];
				}

				Ui[p0][p.TD] += digama(T+phi, &ifault) - (phi+T)/(phi+mu) - log(phi+mu) + c;
				if (ifault == 1) {ll[0] = 4; }
				E[p.TD+1][p.TD+1] += -trigamma(T+phi, &ifault) + (mu-T)/((phi+mu)*(phi+mu)) + 1/(phi+mu);
				if (ifault == 1) {ll[0] = 4; }

			}//with dispersion
			if (p.disp_TREC && strcmp(p.disp_TREC_model, "GP")==0) {

				double mu_T_phi = mu+T*phi;
				for (j = 0; j < p.R; j++) {
					Ui[p0][p.REG+j] += (1.0/mu + (T-1)/mu_T_phi - 1) * mu * vecR[j];
					for (k = j; k < p.R; k++) {
						E[p.REG+j+1][p.REG+k+1] += ((1.0/(mu*mu) + (T-1)/mu_T_phi/mu_T_phi) *mu*mu - (1.0/mu + (T-1)/mu_T_phi - 1) *mu)*vecR[j]*vecR[k];
					}
					E[p.TD+1][p.REG+j+1] += T*(T-1)/mu_T_phi/mu_T_phi * mu*vecR[j];
				}

				Ui[p0][p.TD]		+= T*(T-1)/mu_T_phi - T;
				E[p.TD+1][p.TD+1]	+= T*T*(T-1)/mu_T_phi/mu_T_phi;

			}//with dispersion
			else if (!p.disp_TREC) {

				for (j = 0; j < p.R; j++) {
					Ui[p0][p.REG+j] += (T/mu - 1) * mu * vecR[j];
					for (k = j; k < p.R; k++) {
						E[p.REG+j+1][p.REG+k+1] += T*vecR[j]*vecR[k] - (T/mu - 1) *mu*vecR[j]*vecR[k];
					}
				}

			}//no dispersion

		}//person

		if (p.disp_TREC && strcmp(p.disp_TREC_model, "NB")==0) {
			//common terms
			E[p.TD+1][p.TD+1] += d.n_new * (trigamma(phi, &ifault) - 1/phi );
			if (ifault == 1) {ll[0] = 4;}
		}


		//////////////////////////////////////////////
		// save fixed E22 (TREC)
		//////////////////////////////////////////////


		nui = p.disp_TREC + p.R;


		//E22
		for (j = 1; j <= nui; j++) {
			for (k = j; k <= nui; k++) {
				E22[0][j][k] = E[p.B+p.ASE+j][p.B+p.ASE+k];
				E22[0][k][j] = E22[0][j][k];
			}
		}

		if ((ll[0] = explosion_E(E22[0], nui))) {ll[0] = 3; goto end;}
		invv(E22[0], nui);
		if ((ll[0] = explosion_E(E22[0], nui))) {ll[0] = 3; goto end;}


	}//TREC



	//ASE
	if (p.ASE){

		l_r0_r[0] = 0;   //theta
		p0_r0_r[0] = 0; 
		l_r0_r_2[0] = 0; //theta2
		p0_r0_r_2[0] = 0;
		p0_r0_1[0] = 0; //beta
		p0_r0_1_2[0] = 0; //beta2
		p0_r0_r_d2[0] = 0; //beta*theta
		for (r = 1; r <= max_M; r++) {
			p0_r0_r[r]      = p0_r0_r[r-1]  + (r-1)/  (pn + r0[r]);   
			p0_r0_r_2[r]    = p0_r0_r_2[r-1]  + (r-1)*(r-1)/  ((pn+r0[r])*  (pn+r0[r]));
			p0_r0_1[r]      = p0_r0_1[r-1]   + 1/(pn + r0[r]);
			p0_r0_1_2[r]    = p0_r0_1_2[r-1] + 1/((pn+r0[r])*(pn+r0[r]));
			p0_r0_r_d2[r]   = p0_r0_r_d2[r-1]+ (r-1)/((pn+r0[r])*(pn+r0[r]));
		}      
		for (r = 1; r <= max_N; r++) {
			l_r0_r[r]     = l_r0_r[r-1] + (r-1)/(1+r0[r]);
			l_r0_r_2[r]   = l_r0_r_2[r-1] + (r-1)*(r-1)/((1+r0[r])*(1+r0[r]));
		}


		for (p0 = 0; p0 < d.n; p0++) {
			if (d.rm[p0]) continue;    

            
			for (j = p.B; j < p.FREQ+p.KK-1; j++) Hi[j] = 0;

			int index1 = d.M[p0];
			int index2 = d.N[p0] - d.M[p0];
			int index3 = d.N[p0];
          
			Ui[p0][p.AD]	  += (p0_r0_r[index1] + p0_r0_r[index2] - l_r0_r[index3]);
			E[p.AD+1][p.AD+1] -= (-p0_r0_r_2[index1] - p0_r0_r_2[index2] + l_r0_r_2[index3]);

			for (int pp = 0; pp < d.nhap[p0]; pp++) {    
	
				for (j = p.B; j < p.FREQ+p.KK-1; j++) Hikl[p0][pp][j] = 0.0;

				if (d.N[p0] > 0) {
					Hikl[p0][pp][p.AD] += (p0_r0_r[index1] + p0_r0_r[index2] - l_r0_r[index3]);
					if (d.nhap[p0]==1) Hikl[p0][1][p.AD] = Hikl[p0][0][p.AD]; //Hikl[p0][1][p.AD] will be called later
				}

				//\pi
				int h1 = p.h[d.hap1[p0][pp]];
				int h2 = p.h[d.hap2[p0][pp]];

				if (h1 != 0) {
					Hikl[p0][pp][p.FREQ+h1-1]	+=			 1/v[p.FREQ+h1];
					Ui[p0][p.FREQ+h1-1]			+= d.w[p0][pp]/v[p.FREQ+h1];
					E[p.FREQ+h1][p.FREQ+h1]		+= d.w[p0][pp]/(v[p.FREQ+h1]*v[p.FREQ+h1]);
				}
				if (h2 != 0) {
					Hikl[p0][pp][p.FREQ+h2-1]	+=			 1/v[p.FREQ+h2];
					Ui[p0][p.FREQ+h2-1]			+= d.w[p0][pp]/v[p.FREQ+h2];
					E[p.FREQ+h2][p.FREQ+h2]		+= d.w[p0][pp]/(v[p.FREQ+h2]*v[p.FREQ+h2]);
				}
				int counter = (h1==0) + (h2==0);
				if (counter > 0) {
					for (i = 0; i < p.KK-1; i++) {
						Hikl[p0][pp][p.FREQ+i]		-=             counter/v[p.FREQ];
						Ui[p0][p.FREQ+i]			-= d.w[p0][pp]*counter/v[p.FREQ];
						for (j = i; j < p.KK-1; j++) {
							E[p.FREQ+i+1][p.FREQ+j+1]	+= d.w[p0][pp]*counter/(v[p.FREQ]*v[p.FREQ]);
						}
					}
				}

				for (k = p.B; k < p.FREQ+p.KK-1; k++) {
					Hi[k] += d.w[p0][pp] * Hikl[p0][pp][k];
					for (int l = k; l < p.FREQ+p.KK-1; l++) {
						E[k+1][l+1] -= d.w[p0][pp] * Hikl[p0][pp][k] * Hikl[p0][pp][l]; 
					}
				}

			}//(hk,hl)

			for (k = p.B; k < p.FREQ+p.KK-1; k++) {
				for (int l = k; l < p.FREQ+p.KK-1; l++) {
					E[k+1][l+1] += Hi[k] * Hi[l]; 
				}
			} 

		}//person

		//E_Louis_ext
		//if (d.ext) Louis_orig_ext(E, U, d, p, v); //need to get Ui instead of U


		//////////////////////////////////////////////
		// save fixed E22 (TRECASE)
		//////////////////////////////////////////////


		nui = p.FREQ+p.KK-1-p.B;

		//E22
		for (j = 1; j <= nui; j++) {
			for (k = j; k <= nui; k++) {
				E22[1][j][k] = E[p.B+j][p.B+k];
				E22[1][k][j] = E22[1][j][k];
			}
		}

		if ((ll[1] = explosion_E(E22[1], nui))) {ll[1] = 3; goto end;}
		invv(E22[1], nui);
		if ((ll[1] = explosion_E(E22[1], nui))) {ll[1] = 3; goto end;}


	}//ASE


	//p.B should be 1

	dim = 0;
	if (ll[0] == 0) dim++;
	if (ll[1] == 0) dim++;

	Usum = new double[dim];
	Vsum = dmatrix(0, dim-1, 0, dim-1);


	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	//                                             permutation
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////


	for (i = 0; i < d.n; i++) permi[i] = i;


	for (iperm = 0; iperm < nPERM; iperm++) {

		for (i = 0; i < d.n; i++) {
			rand1[i] = kiss(seed);
			rand2[i] = kiss(seed);
		}


		permute_sample_int(permi, d.n, rand1, d.rm);
		
		////////////////////////////////////////////////////
		// evaluate U[p.B], E[p.B][p.B]
		////////////////////////////////////////////////////


        //U, E initialation
        //U = l', E = -l"
		for (j = 0; j < p.B; j++) {
			for (p0 = 0; p0 < d.n; p0++) Ui[p0][j] = 0;
			for (k = j; k < p.nall; k++) E[j+1][k+1] = 0;
		}
		

		if (p.TREC) {

			for (p0 = 0; p0 < d.n; p0++) {
				if (d.rm[p0]) continue;    
            
				Bvec(vecB, p.B, d.X[p0]);
				Rvec(vecR, p.R, d.X[p0]); //debugged: line missed


				double T = d.T[p0];
				double mu = mu_i[p0];
				double tmu = 0.5 * d.G[permi[p0]] * exp_eta[p0];


				if (p.disp_TREC && strcmp(p.disp_TREC_model, "NB")==0) 
				{
					for (j = 0; j < p.B; j++) {
						Ui[p0][j] += (T/mu - (phi+T)/(phi+mu)) * tmu * vecB[j];
						for (k = j; k < p.B; k++) {
							E[j+1][k+1]+= ((T/(mu*mu) - (phi+T)/((phi+mu)*(phi+mu))) *tmu*tmu - (T/mu - (phi+T)/(phi+mu)) * tmu) * vecB[j]*vecB[k];
						}
						E[j+1][p.TD+1] += (1/(phi+mu) - (phi+T)/((phi+mu)*(phi+mu))) * tmu * vecB[j];
						for (k = 0; k < p.R; k++) {
							E[j+1][p.REG+k+1] += (T/(mu*mu) - (phi+T)/((phi+mu)*(phi+mu))) * tmu*mu*vecB[j]*vecR[k] - (T/mu - (phi+T)/(phi+mu)) * tmu*vecB[j]*vecR[k];
						}
					}

				}//with dispersion
				if (p.disp_TREC && strcmp(p.disp_TREC_model, "GP")==0) 
				{
					double mu_T_phi = mu+T*phi;
					for (j = 0; j < p.B; j++) {
						Ui[p0][j] += (1.0/mu + (T-1)/mu_T_phi - 1) * tmu * vecB[j];
						for (k = j; k < p.B; k++) {
							E[j+1][k+1]+= ((1.0/(mu*mu) + (T-1)/mu_T_phi/mu_T_phi) *tmu*tmu - (1.0/mu + (T-1)/mu_T_phi - 1) *tmu) * vecB[j]*vecB[k];
						}
						E[j+1][p.TD+1] += T*(T-1)/mu_T_phi/mu_T_phi * tmu * vecB[j];
						for (k = 0; k < p.R; k++) {
							E[j+1][p.REG+k+1] += ((1.0/(mu*mu) + (T-1)/mu_T_phi/mu_T_phi) *tmu*mu - (1.0/mu + (T-1)/mu_T_phi - 1) *tmu)*vecB[j]*vecR[k];
						}
					}
				}//with dispersion
				else if (!p.disp_TREC) 
				{
					for (j = 0; j < p.B; j++) {
						Ui[p0][j] += (T/mu - 1) * tmu * vecB[j];
						for (k = j; k < p.B; k++) {
							E[j+1][k+1]+= (T/(mu*mu) *tmu*tmu - (T/mu - 1) * tmu) * vecB[j]*vecB[k];
						}
						for (k = 0; k < p.R; k++) {
							E[j+1][p.REG+k+1] += T/(mu*mu) * tmu*mu*vecB[j]*vecR[k] - (T/mu - 1) * tmu*vecB[j]*vecR[k];
						}
					}

				}//no dispersion

			}//person

	
			//////////////////////////////////////////////
			// evaluate efficient score (TREC)
			//////////////////////////////////////////////

			
			nui = p.disp_TREC + p.R;


			//E12_22 = E12 %*% E22
			for (k = 1; k <= p.B; k++) {
				for (i = 1; i <= nui; i++) {
					E12_22[k][i] = 0.0;
					for (j = 1; j <= nui; j++) E12_22[k][i] += E[k][p.B+p.ASE+j] * E22[0][j][i];
				}
			}

			//Uei
			for (p0 = 0; p0 < d.n; p0++) {
				if (d.rm[p0]) {
					for (k = 0; k < p.B; k++) Uei[0][p0][k] = 0.0;
					continue;
				}
				for (k = 0; k < p.B; k++) {
					Uei[0][p0][k] = Ui[p0][k];
					for (j = 0; j < nui; j++) Uei[0][p0][k] -= E12_22[k+1][j+1] * Ui[p0][p.B+p.ASE+j];
				}
			}

			//E11 Fisher
			for (k = 1; k <= p.B; k++) {
				for (i = 1; i <= p.B; i++) {
					E11[0][k][i] = E[k][i];
					for (j = 1; j <= nui; j++) E11[0][k][i] -= E12_22[k][j] * E[i][p.B+p.ASE+j];
				}
			}


		}//TREC


		//ASE
		if (p.ASE){

			for (p0 = 0; p0 < d.n; p0++) {
				if (d.rm[p0]) continue;    

				for (j = 0; j < p.FREQ+p.KK-1; j++) Hi[j] = 0;

				Bvec(vecB, p.B, d.X[p0]);

				//default
				int index1 = d.M[p0];
				int index2 = d.N[p0] - d.M[p0];
				int index3 = d.N[p0];
          
				int switch_hap = (rand2[p0] < 0.5); 
				
				for (int pp = 0; pp < d.nhap[permi[p0]]; pp++) {    
	
					for (j = 0; j < p.B; j++) Hikl[p0][pp][j] = 0.0;

					if (d.G[permi[p0]] == 1 && d.N[p0] > 0) {            //N>0, G=1
                
						int d_hap1_permip0_pp = d.hap1[permi[p0]][pp];
						if (switch_hap) d_hap1_permip0_pp = d.hap2[permi[p0]][pp];


						if (d_hap1_permip0_pp / pow1 == 0) {     
							index1 = d.N[p0] - d.M[p0];
							index2 = d.M[p0];
						}
						else {
							index1 = d.M[p0];
							index2 = d.N[p0] - d.M[p0];
						}


						double exp_beta_tmp = 1.0;
						double p1_tmp = pn;
						double l_p1_tmp = pn;

						for (k = 0; k < p.B; k++) {
							Hikl[p0][pp][k]	+=						( p0_r0_1[index1]    - p0_r0_1[index2]) * p1_tmp * l_p1_tmp * vecB[k];
							Ui[p0][k]		+= d.w[permi[p0]][pp] * ( p0_r0_1[index1]    - p0_r0_1[index2]) * p1_tmp * l_p1_tmp * vecB[k];
							E[k+1][p.AD+1]	-= d.w[permi[p0]][pp] * (-p0_r0_r_d2[index1] + p0_r0_r_d2[index2]) * p1_tmp * l_p1_tmp * vecB[k];
							for (int l = k; l < p.B; l++) {
								E[k+1][l+1] -= d.w[permi[p0]][pp] * (-p0_r0_1_2[index1]  - p0_r0_1_2[index2]) * p1_tmp*p1_tmp*l_p1_tmp*l_p1_tmp * vecB[k]*vecB[l]
								             + d.w[permi[p0]][pp] * ( p0_r0_1[index1]    - p0_r0_1[index2]) * exp_beta_tmp * (1-exp_beta_tmp)/((1+exp_beta_tmp)*(1+exp_beta_tmp)*(1+exp_beta_tmp)) * vecB[k]*vecB[l];
							}
						}
					} //N>0, G=1

					for (k = 0; k < p.B+1; k++)					Hi[k] += d.w[permi[p0]][pp] * Hikl[p0][pp][k];
					for (k = p.FREQ; k < p.FREQ+p.KK-1; k++)	Hi[k] += d.w[permi[p0]][pp] * Hikl[permi[p0]][pp][k];
					

					for (k = 0; k < p.B; k++) {
						for (l = k; l <= p.B; l++) {
							E[k+1][l+1] -= d.w[permi[p0]][pp] * Hikl[p0][pp][k] * Hikl[p0][pp][l]; 
						}
						for (l = p.B+1; l < p.FREQ+p.KK-1; l++) {
							E[k+1][l+1] -= d.w[permi[p0]][pp] * Hikl[p0][pp][k] * Hikl[permi[p0]][pp][l]; 
						}
					}

				}//(hk,hl)

				for (k = 0; k < p.B; k++) {
					for (int l = k; l < p.FREQ+p.KK-1; l++) {
						E[k+1][l+1] += Hi[k] * Hi[l]; 
					}
				} 

			}//person

			//E_Louis_ext
			//if (d.ext) Louis_orig_ext(E, U, d, p, v); //need to get Ui instead of U


			//////////////////////////////////////////////
			// evaluate efficient score (TRECASE)
			//////////////////////////////////////////////


			nui = p.FREQ+p.KK-1-p.B;

			//E12_22 = E12 %*% E22
			for (k = 1; k <= p.B; k++) {
				for (i = 1; i <= nui; i++) {
					E12_22[k][i] = 0.0;
					for (j = 1; j <= nui; j++) E12_22[k][i] += E[k][p.B+j] * E22[1][j][i];
				}
			}

			//Uei
			for (p0 = 0; p0 < d.n; p0++) {
				if (d.rm[p0]) {
					for (k = 0; k < p.B; k++) Uei[1][p0][k] = 0.0;
					continue;
				}
				for (k = 0; k < p.B; k++) {
					Uei[1][p0][k] = Ui[p0][k];
					for (j = 0; j < p.FREQ-p.B; j++)					Uei[1][p0][k] -= E12_22[k+1][j+1] * Ui[p0][p.B+j];
					for (j = p.FREQ-p.B; j < p.FREQ+p.KK-1-p.B; j++)	Uei[1][p0][k] -= E12_22[k+1][j+1] * Ui[permi[p0]][p.B+j];
				}
			}

			//E11 Fisher
			for (k = 1; k <= p.B; k++) {
				for (i = 1; i <= p.B; i++) {
					E11[1][k][i] = E[k][i];
					for (j = 1; j <= nui; j++) E11[1][k][i] -= E12_22[k][j] * E[i][p.B+j];
				}
			}


		}//ASE


		//U as a check
		for (j = 0; j < p.nall; j++) U[j] = 0.0;
		for (p0 = 0; p0 < d.n; p0++) {
			if (d.rm[p0]) continue;
			for (j = 0; j < p.nall; j++) U[j] += Ui[p0][j];
		}

		//===================
		//get pvalue[iperm]
		//===================

		//p.B should be 1

		double r_stat_mar[2];


		for (k = 0; k < dim; k++) { //dim

			double Umax = 0.0, Umin = 0.0;
				
			Usum[k] = 0.0;
			Vsum[k][k] = 0.0;
			for (i = 0; i < d.n; i++) {
				Usum[k] += Uei[k][i][0];
				Vsum[k][k] += Uei[k][i][0]*Uei[k][i][0];

				if (Uei[k][i][0] < Umin)		Umin = Uei[k][i][0];
				else if (Uei[k][i][0] > Umax)	Umax = Uei[k][i][0];
			}

			r_stat_mar[k] = fabs(Usum[k] / sqrt(Vsum[k][k]-Usum[k]*Usum[k]/(d.n))); 
			//debugged: Uei tends to be big for heterozygotes; Vsum[k][k] tends to be big; r_stat_mar[1] < r_stat_mar[0]
		
			
			if (r_stat_mar[k] > stat_perm_pmin[k][iperm]) {
				stat_perm_pmin[k][iperm] = r_stat_mar[k];
			}
		}//k
		 
	}//iperm


	delete [] Usum;
	free_dmatrix(Vsum, 0, dim-1, 0, dim-1);


end:

	delete [] rand1;
	delete [] rand2;

	delete [] permi;
	delete [] v; 
  	delete [] U;
	free_dmatrix(Ui,  0, d.n-1, 0, p.nall-1);
	free_dmatrix(Ve,  1, p.B,   1, p.B);
	free_dmatrix(E,   1, p.nall, 1, p.nall);
	free_dmatrix(E12_22, 1, p.B, 1, p.nall-p.B);
	free_dmatrix3(E22,  0, 1, 1, p.nall-p.B, 1, p.nall-p.B);
	free_dmatrix3(E11,  0, 1, 1, p.B, 1, p.B);
    free_dmatrix3(Uei,  0, 1, 0, d.n-1, 0, p.B-1);

	delete [] vecB;

	if (p.ASE) {

		delete [] Hi;
		free_dmatrix3(Hikl, 0, d.n-1, 0, 1, 0, p.FREQ+p.KK-2);

		delete [] r0;
		delete [] p0_r0;
		delete [] l_r0;    
		delete [] l_r0_r;
		delete [] p0_r0_r;
		delete [] l_r0_r_2;
		delete [] p0_r0_r_2;  
		delete [] p0_r0_1;
		delete [] p0_r0_1_2;
		delete [] p0_r0_r_d2;
	}

	if (p.TREC){
		delete [] vecR;
		delete [] mu_i;
		delete [] exp_eta;
	}


}//permutation


double corr_int(int*Y, int*X, int n) {
	int n_true = 0;
	double XY = 0.0, X2 = 0.0, Y2 = 0.0, Xmean = 0.0, Ymean = 0.0;
	for (int i = 0; i < n; i++) {
		if (Y[i] != 9 && X[i] != 9) {
			XY += X[i] * Y[i];
			X2 += X[i] * X[i];
			Y2 += Y[i] * Y[i];
			Xmean += 1.0*X[i];
			Ymean += 1.0*Y[i];
			n_true++;
		}
	}
	Xmean /= n_true;
	Ymean /= n_true;

	return (XY-n_true*Xmean*Ymean)/sqrt((X2-n_true*Xmean*Xmean)*(Y2-n_true*Ymean*Ymean));
}//corr


int eQTL_mapping_per_gene(FILE*foutpmin[2], int chr_working, int itx, char**txName, int*txStart, int*txEnd, 
						  char**rsSNP, int*posSNP, double*freqSNP, double MAF_cutoff, double corr2_pruning,
						  int targetSNPStart, int targetSNPEnd, int tranSNPStart, int tranSNPEnd, DATA d, PARA p, RESULTS r, int**haplo) {

	//-------------------------------------------------------------------
	//--------------------- get initial para values ---------------------
	//-------------------------------------------------------------------

	int ll = 0;

	r.itx = itx;

	int i, j;	

	//get "ASE_fitted_gene"

	int TREC_fitted_gene = 1;
	int ASE_fitted_gene = 1;

	const int min_nAsCount		= 5;	//minimum of total number of reads
	const int min_nSubj_asCount = 5;	//minimum # of samples with >= min_nT reads
	const int min_nSubj_hetSNP  = 5;	//minimum # of samples with heterzygous genotypes

	int nsubj_asCount = 0;
	for (i = 0; i < d.n; i++) if (d.N[i] >= min_nAsCount) nsubj_asCount++;
	if(nsubj_asCount <= min_nSubj_asCount) ASE_fitted_gene = 0;	
	if (tranSNPEnd < tranSNPStart) ASE_fitted_gene = 0;


	//"get_initial_para"

	double mu_T = 0.0, var_T = 0.0;
	for (i = 0; i < d.n; i++) {
		mu_T       += d.T[i];
		var_T += 1.0*d.T[i]*d.T[i];
	}
	mu_T /= d.n;
	var_T = var_T/d.n - mu_T*mu_T;

	double phi_T_orig = 1.0/((var_T-mu_T)/mu_T/mu_T);


	r.ll_gene[0] = 9;
	r.ll_gene[1] = 9;
	double*vTREC_gene = new double[1+p.R];
	double*vASE_gene = new double[1];

	get_initial_para(r.ll_gene, vTREC_gene, vASE_gene, &d, &p, TREC_fitted_gene, ASE_fitted_gene);

	if (r.ll_gene[1] != 0) ASE_fitted_gene = 0;
	

	if (r.ll_gene[0] == 0 && targetSNPStart <= targetSNPEnd) {

		// screen targetSNP & tranSNP by "frequency", "missingness" and "LD"


		int pwin_r_tarSNP = 0;
		int*map_rf_tarSNP = new int[targetSNPEnd - targetSNPStart + 1];
		int*map_fr_tarSNP = new int[targetSNPEnd - targetSNPStart + 1];
		screen_SNP_LD(map_rf_tarSNP, &pwin_r_tarSNP, d.n, targetSNPStart, targetSNPEnd, haplo, corr2_pruning);
		for (j = 0; j < targetSNPEnd - targetSNPStart + 1; j++) map_fr_tarSNP[j] = 0;
		for (j = 0; j < pwin_r_tarSNP; j++) map_fr_tarSNP[map_rf_tarSNP[j]-targetSNPStart] = 1;
		


		int pwin_r_tranSNP = 0, *map_rf_tranSNP = NULL;		
		if (ASE_fitted_gene) {
			map_rf_tranSNP = new int[tranSNPEnd - tranSNPStart + 1]; 
			screen_SNP_LD(map_rf_tranSNP, &pwin_r_tranSNP, d.n, tranSNPStart, tranSNPEnd, haplo, corr2_pruning);
		}

		
		// define result variable

		double MD_pmin[2] = {-1,-1}, freq_pmin[2] = {0,0}, stat_orig_pmin[2] = {-1,-1};
		int    iSNP_pmin[2] = {-1,-1};

	
		const int nPERM = 5000;

		double stat_perm_pmin[2][nPERM] = {0.0};


		unsigned long seed_tmp0 = seed[0];
		unsigned long seed_tmp1 = seed[1];
		unsigned long seed_tmp2 = seed[2];


		for (int iSNP = targetSNPStart; iSNP <= targetSNPEnd; iSNP++) { // targetSNPStart ~ targetSNPEnd

			if (freqSNP[iSNP] < MAF_cutoff || freqSNP[iSNP] > 1.0-MAF_cutoff) continue;

			if (map_fr_tarSNP[iSNP-targetSNPStart]==0) continue;

			r.iSNP = iSNP;

			int in_range = (tranSNPStart <= iSNP && iSNP <= tranSNPEnd); 
			int known_Gp = in_range;


			//-------------------------------------------------------------------
			//-------------------- copy for actual analysis -------------------
			//-------------------------------------------------------------------
			// d.T, d.M, d.G, d.Gp, d.H
			// ASE_fitted, d.M, d.K
			//-------------------------------------------------------------------

			double freq = 0.0;
			int G0 = 0, G1 = 0, G2 = 0;

			for (i = 0; i < d.n; i++) {

				int tmp = haplo[iSNP][i];

				d.G[i] = (tmp>2) ? (tmp-2) : tmp;
				d.Gp[i][0] = tmp / 3;
				d.Gp[i][1] = tmp % 3; 

				freq += d.G[i];
				if (d.G[i]==0) G0++;
				else if (d.G[i]==1) G1++;
				else if (d.G[i]==2) G2++;

			}

			freq /= 2*d.n;

			//get "ASE_fitted_gene"

			p.TREC_fitted = TREC_fitted_gene;
			p.ASE_fitted  = ASE_fitted_gene;

			int nsubj_hetSNP_asCount = 0;
			if (p.ASE_fitted == 1) {
				nsubj_hetSNP_asCount = 0;
				for (i = 0; i < d.n; i++) {
					if (d.G[i] == 1 && d.N[i] >= min_nAsCount) nsubj_hetSNP_asCount++;
				}
				if (nsubj_hetSNP_asCount <= min_nSubj_hetSNP) p.ASE_fitted = 0;	
			}


			const int Mmax = 5;

			if (p.ASE_fitted && !known_Gp) {

				int *tagSNP = new int[Mmax-1];

				get_tagSNP(tagSNP, &(p.M), &(r.MD), iSNP, d.n, map_rf_tranSNP, pwin_r_tranSNP, haplo, Mmax, p.pow2);
				
				p.K = int(pow(2,p.M)); 


				for (i = 0; i < d.n; i++) {
					d.H[i][0] = 0;
					d.H[i][1] = 0;
					for (int k = 0; k < p.M-1; k++) {
						d.H[i][0] += haplo[tagSNP[k]][i] / 3 * p.pow2[p.M-2-k]; 
						d.H[i][1] += haplo[tagSNP[k]][i] % 3 * p.pow2[p.M-2-k];
					}
				}

				delete [] tagSNP;

			}
			else if (p.ASE_fitted && known_Gp) {
				r.MD = -1;
				p.M = 1;
				p.K = 2; 

			}
			else if (!p.ASE_fitted) {
				r.MD = -100;
			}


			//-------------------------------------------------------------------
			//------------------ data structure: link list ----------------------
			//-------------------------------------------------------------------

			if (p.ASE_fitted && !known_Gp)
			{
				makelist_unlinklist(d.hap1, d.hap2, d.nhap, d.rm, d.G, d.H, d.n, p.M, p.pow2);

				get_hap_freq_unlinklist(&p, &d);

			}
			else if (p.ASE_fitted && known_Gp)
			{
				makelist_unlinklist_Gp(d.hap1, d.hap2, d.nhap, d.rm, d.n, d.Gp);

				get_hap_freq_unlinklist(&p, &d);
			}


			//-------------------------------------------------------------------
			//-------------------- original data analysis -----------------------
			//-------------------------------------------------------------------

			// initiation

			for (j = 0; j < r.nTEST; j++)	{r.ll[j] = 9; r.ll_null[j] = 9;}
			for (j = 0; j < 1+p.R; j++)		r.vTREC[j] = vTREC_gene[j];
			for (j = 0; j < 1; j++)			r.vASE[j] = vASE_gene[j];

			// core function

			get_null_para_loglike(&r, &d, &p);
			wrap_TRECASE_score(&r, &d, &p);


//			if (p.ASE_fitted) printf("%d: Rsq=%5.2lf\n", iSNP-targetSNPStart, r.MD);
//			else printf("%d: \n", iSNP-targetSNPStart);

			for (int t = 0; t < r.n_test; t++) {		
				if (r.stat_mar[t] > stat_orig_pmin[t]) 
				{
					MD_pmin[t] = r.MD;
					freq_pmin[t] = freq;
					iSNP_pmin[t] = r.iSNP;

					stat_orig_pmin[t] = r.stat_mar[t];
				}
			}

			
			//-------------------------------------------------------------------
			//--------------------       permutation      -----------------------
			//-------------------------------------------------------------------

			seed[0] = seed_tmp0;
			seed[1] = seed_tmp1;
			seed[2] = seed_tmp2;

			permutation (stat_perm_pmin, nPERM, d, p, r.vTREC, r.vASE);


		}//iSNP (target SNPs)


		for (int t = 0; t < 1+ASE_fitted_gene; t++) {
			if (iSNP_pmin[t]==-1) continue; //debugged on Z2014-02-09

			double p_orig_pmin = gammq(0.5*p.B, 0.5*stat_orig_pmin[t]*stat_orig_pmin[t]); 

			double p_perm_pmin = 0.0; 
			for (i = 0; i < nPERM; i++) {
				if (stat_perm_pmin[t][i] > stat_orig_pmin[t]) p_perm_pmin += 1.0;
			}
			p_perm_pmin /= nPERM;


			fprintf(foutpmin[t],	"%3d  %15s %9d %9d    %5d       %5d    %5d       %5d    %10s   %9d      %5.3lf", 
								chr_working, txName[itx], txStart[itx], txEnd[itx], tranSNPEnd-tranSNPStart+1, pwin_r_tranSNP, targetSNPEnd-targetSNPStart+1, pwin_r_tarSNP, rsSNP[iSNP_pmin[t]], posSNP[iSNP_pmin[t]], freq_pmin[t]);
				
			if (MD_pmin[t] >= -1)	fprintf(foutpmin[t], "       %5.3lf", fabs(MD_pmin[t]));
			else					fprintf(foutpmin[t], "         -  ");

			fprintf(foutpmin[t],	"    %8.1e %8.1e", p_orig_pmin, p_perm_pmin);
		

			//-------------------------------------------------------------------
			//-----------  cis-trans test, beta estimation  ---------------------
			//-------------------------------------------------------------------


			for (int iSNP = iSNP_pmin[t]; iSNP <= iSNP_pmin[t]; iSNP++) { // targetSNPStart ~ targetSNPEnd

				if (map_fr_tarSNP[iSNP-targetSNPStart]==0) continue;

				r.iSNP = iSNP;

				int in_range = (tranSNPStart <= iSNP && iSNP <= tranSNPEnd); 
				int known_Gp = in_range;


				//-------------------------------------------------------------------
				//-------------------- copy for actually analysis -------------------
				//-------------------------------------------------------------------
				// d.T, d.M, d.G, d.Gp, d.H
				// ASE_fitted, d.M, d.K
				//-------------------------------------------------------------------

				double freq = 0.0;
				int G0 = 0, G1 = 0, G2 = 0;

				for (i = 0; i < d.n; i++) {

					int tmp = haplo[iSNP][i];

					d.G[i] = (tmp>2) ? (tmp-2) : tmp;
					d.Gp[i][0] = tmp / 3;
					d.Gp[i][1] = tmp % 3; 

					freq += d.G[i];
					if (d.G[i]==0) G0++;
					else if (d.G[i]==1) G1++;
					else if (d.G[i]==2) G2++;

				}

				freq /= 2*d.n;

				//get "ASE_fitted_gene"

				p.TREC_fitted = TREC_fitted_gene;
				p.ASE_fitted  = ASE_fitted_gene;

				int nsubj_hetSNP_asCount = 0;
				if (p.ASE_fitted == 1) {
					nsubj_hetSNP_asCount = 0;
					for (i = 0; i < d.n; i++) {
						if (d.G[i] == 1 && d.N[i] >= min_nAsCount) nsubj_hetSNP_asCount++;
					}
					if (nsubj_hetSNP_asCount <= min_nSubj_hetSNP) p.ASE_fitted = 0;	
				}


				const int Mmax = 5;

				if (p.ASE_fitted && !known_Gp) {

					int *tagSNP = new int[Mmax-1];

					get_tagSNP(tagSNP, &(p.M), &(r.MD), iSNP, d.n, map_rf_tranSNP, pwin_r_tranSNP, haplo, Mmax, p.pow2);
					
					p.K = int(pow(2,p.M)); 


					for (i = 0; i < d.n; i++) {
						d.H[i][0] = 0;
						d.H[i][1] = 0;
						for (int k = 0; k < p.M-1; k++) {
							d.H[i][0] += haplo[tagSNP[k]][i] / 3 * p.pow2[p.M-2-k]; 
							d.H[i][1] += haplo[tagSNP[k]][i] % 3 * p.pow2[p.M-2-k];
						}
					}

					delete [] tagSNP;

				}
				else if (p.ASE_fitted && known_Gp) {
					r.MD = -1;
					p.M = 1;
					p.K = 2; 

				}
				else if (!p.ASE_fitted) {
					r.MD = -100;
				}


				//-------------------------------------------------------------------
				//------------------ data structure: link list ----------------------
				//-------------------------------------------------------------------

				if (p.ASE_fitted && !known_Gp)
				{
					makelist_unlinklist(d.hap1, d.hap2, d.nhap, d.rm, d.G, d.H, d.n, p.M, p.pow2);

					get_hap_freq_unlinklist(&p, &d);

				}
				else if (p.ASE_fitted && known_Gp)
				{
					makelist_unlinklist_Gp(d.hap1, d.hap2, d.nhap, d.rm, d.n, d.Gp);

					get_hap_freq_unlinklist(&p, &d);
				}


				//-------------------------------------------------------------------
				//-------------------- original data analysis -----------------------
				//-------------------------------------------------------------------

				// initiation

				for (j = 0; j < r.nTEST; j++)	{r.ll[j] = 9; r.ll_null[j] = 9;}
				for (j = 0; j < 1+p.R; j++)		r.vTREC[j] = vTREC_gene[j];
				for (j = 0; j < 1; j++)			r.vASE[j] = vASE_gene[j];

				// core function

				get_null_para_loglike(&r, &d, &p);
				wrap_TRECASE_score(&r, &d, &p);

				int ll_TReCASE;
				double v_TReCASE;
				if (p.ASE_fitted) {

					wrap_est_TRECASE(&r, &d, &p);
					ll_TReCASE = r.ll[2];
					v_TReCASE = r.v[2];

					if (r.ll[2]==0) {
						wrap_test_cis(&r, &d, &p);
						double stat_orig_cis = r.Ucis/r.Vcis;
						double p_orig_cis = gammq(0.5*p.B, 0.5*stat_orig_cis*stat_orig_cis); 

						fprintf(foutpmin[t], " %8.1e", p_orig_cis);
					}
					else fprintf(foutpmin[t], "     -   ");
				}
				else fprintf(foutpmin[t], "     -   ");


				if (t==0) {
					p.TREC = 1;
					p.ASE = 0;
					wrap_est_TRECASE(&r, &d, &p);
					if (r.ll[2]==0) fprintf(foutpmin[t],	"      %7.4lf", r.v[2]);
					else fprintf(foutpmin[t], "         -   ");
				}
				else if (t==1) {
					// estimate r.v[2]

					p.TREC = 1;
					p.ASE = 0;
					wrap_est_TRECASE(&r, &d, &p);
					if (r.ll[2]==0) fprintf(foutpmin[t],	"      %7.4lf", r.v[2]);
					else fprintf(foutpmin[t], "         -   ");

					if (p.ASE_fitted) {
						p.TREC = 0;
						p.ASE = 1;
						wrap_est_TRECASE(&r, &d, &p);
						if (r.ll[2]==0) fprintf(foutpmin[t],	"  %7.4lf", r.v[2]);
						else fprintf(foutpmin[t], "     -   ");

						if (ll_TReCASE==0) fprintf(foutpmin[t],	"      %7.4lf", v_TReCASE);
						else fprintf(foutpmin[t], "         -   ");
					}
				}
				
			}//iSNP (minimum-p SNP)

			fprintf(foutpmin[t], "\n");
			
		}//t

		if (ASE_fitted_gene) delete [] map_rf_tranSNP;
		delete [] map_rf_tarSNP;
		delete [] map_fr_tarSNP;

	}//r.ll_gene[0] == 0


	delete [] vTREC_gene;
	delete [] vASE_gene;


	return ll;

}//eQTL_mapping_per_gene_score




void eQTL_mapping(FILE*foutput_file[], int start_tx, int end_tx, int chr_working, int local_SNP_win, double corr2_pruning,
				  int nsubj, int ncov, double**cov,
				  int ntx, char**txName, int*txStart, int*txEnd,
				  int**toCount, int**asCountA, int**asCountB, int T_max,
				  int nSNP, char**rsSNP, int*posSNP, double*freqSNP, double MAF_cutoff, 
				  int**haplo) 
{	

	int  isubj/*subjects*/, itx/*transcripts*/, i;

	//-------------------------------------------------------------------------------
	//------------------------------ data structure ---------------------------------
	//-------------------------------------------------------------------------------


	PARA p;

	const int Mmax = 5;  //longest haplotype spanning target + transcribed SNPs
	int Kmax = int(pow(2,Mmax));
  char nb[3] = "NB";
	p.THRESHOLD_L   =  0.0008;//1.0/sim_nsubj; debugged: removing some ind causes not enough dependence between SNPs in permutaiton  
    p.pow2          = new int[Mmax];
	p.pow2[0] = 1; for (i = 1; i < Mmax; i++)  p.pow2[i] = p.pow2[i-1]*2; 
	p.h        = new int[Kmax];
	p.hap_freq = new int[Kmax];
	p.pai	   = new double[Kmax];
	p.disp_TREC = 1;			//0: Poisson data/model; 1: NB data/model and GP data/model;
	p.disp_TREC_model = nb;	//NB or GP

	p.B = 1; 
	p.R = 1+ncov; // intercept+logcntAll+3PC	 //other regression para of TReC model excluding SNPs


	DATA d;

	d.n = nsubj;
	d.T =  new int[d.n];			//total count
	d.N =  new int[d.n];			//total asCount
    d.M =  new int[d.n];			//asCount on one haplotype
	d.G =  new int[d.n];			//genotype at the target SNP
    d.Gp=  imatrix(0, d.n-1, 0, 1); //phased genotype at the target SNP
	d.H =  imatrix(0, d.n-1, 0, 1); //haplotype at the (M-1) transcribed SNPs
	d.X =  dmatrix(0, d.n-1, 0, p.R-2); //covariates (logcntAll, pc, gender, ...)
	d.hap1 =  imatrix(0, d.n-1, 0, 1); //linklist haplotype 1, maximum possibilities = 2
	d.hap2 =  imatrix(0, d.n-1, 0, 1); //linklist haplotype 2, maximum possibilities = 2
	d.prob =  dmatrix(0, d.n-1, 0, 1); //linklist prob in A1
	d.w    =  dmatrix(0, d.n-1, 0, 1); //linklist w in EM algorithm
	d.sum  =  new double[d.n];         //linklist sum
	d.nhap =  new int[d.n];			   //# of haplo possibilities
	d.rm   =  new int[d.n];            //1: remove the individual; 0: not remove

	d.T_max = T_max;
	d.lgamma_int = new double[d.T_max+1];
	d.lgamma_int[0]=0;for (i = 1; i <= d.T_max;i++) d.lgamma_int[i] = d.lgamma_int[i-1]+log(i);

	for (isubj = 0; isubj < d.n; isubj++) {
		for (int icov = 0; icov < p.R-1; icov++) {
			d.X[isubj][icov] = cov[isubj][icov];
		}
	}


	RESULTS r;
	r.nTEST = 3;
	r.vTREC = new double[1+p.R];	
	r.vASE = new double[1];
	r.U = dmatrix3(0, r.nTEST-1, 0, d.n-1, 0, p.B-1);
	r.V = dmatrix3(0, r.nTEST-1, 1, p.B, 1, p.B);


	//-------------------------------------------------------------------------------
	//--------------------------- Analysis by Transcripts   -------------------------
	//-------------------------------------------------------------------------------
		

	for (itx = start_tx; itx <= end_tx; itx++) {

		if (itx >= ntx) break;

		printf("\nGene %d\n", itx+1);

		int targetSNPStart = binsearch_int_2(txStart[itx]-local_SNP_win, posSNP, nSNP, 0); //posSNP is ordered, 1-left, 0-right
		int targetSNPEnd   = binsearch_int_2(txEnd[itx]+local_SNP_win, posSNP, nSNP, 1); 

		int tranSNPStart   = binsearch_int_2(txStart[itx], posSNP, nSNP, 0); 
		int tranSNPEnd     = binsearch_int_2(txEnd[itx], posSNP, nSNP, 1); 


		for (isubj = 0; isubj < d.n; isubj++) {
			d.T[isubj] = toCount[itx][isubj];
			d.N[isubj] = asCountA[itx][isubj] + asCountB[itx][isubj];
			d.M[isubj] = asCountA[itx][isubj];
		}


		//-------------------------------------
		//---------- core function ------------
		//-------------------------------------

		int ll = eQTL_mapping_per_gene(foutput_file, chr_working, itx, txName, txStart, txEnd, rsSNP, posSNP, freqSNP, MAF_cutoff, corr2_pruning,
			targetSNPStart, targetSNPEnd, tranSNPStart, tranSNPEnd, d, p, r, haplo);

	}//itx


	delete [] d.T;
	delete [] d.N;
	delete [] d.M;
    delete [] d.G;
    free_imatrix(d.Gp, 0, d.n-1, 0, 1);
	free_imatrix(d.H,  0, d.n-1, 0, 1);
	free_dmatrix(d.X,  0, d.n-1, 0, p.R-2);
	free_imatrix(d.hap1, 0, d.n-1, 0, 1); 
	free_imatrix(d.hap2, 0, d.n-1, 0, 1); 
	free_dmatrix(d.prob, 0, d.n-1, 0, 1); 
	free_dmatrix(d.w   , 0, d.n-1, 0, 1); 
	delete [] d.nhap;	
	delete [] d.rm;
	delete [] d.sum;
	delete [] d.lgamma_int;

	delete [] p.pow2;
	delete [] p.h;
	delete [] p.hap_freq;
	delete [] p.pai;

	delete [] r.vTREC;
	delete [] r.vASE;
	free_dmatrix3(r.U,	0, r.nTEST-1, 0, d.n-1, 0, p.B-1);
	free_dmatrix3(r.V,	0, r.nTEST-1, 1, p.B, 1, p.B);

}//eQTL_mapping


void wrap_TRECASE_score(RESULTS*r, DATA*d, PARA*p) {

	//ll_null[0,1,2]: TREC null, TRECASE null, TRECASE null
	//loglike_null[0,1,2]: TREC null, ASE null, TRECASE null
	//ll[0,1,2]: score for TREC; score for TRECASE; score for trans

	if (p->TREC_fitted && r->ll_null[0]<=1) {
		p->TREC = 1;
		if (p->ASE_fitted && r->ll_null[1]<=1) {
			p->ASE = 1;
		}
		get_index(p->ASE, p->TREC, p);

		evaluate_effiScore_TREC_TRECASE(r->ll, &r->n_test, r->stat_mar, *d, *p, r->U, r->vASE, r->vTREC);

		if		(r->ll[0] == 4)	printf("evaluate_effiScore_TRECASE failed!: ifault==1!\n");
		else if (r->ll[0] == 3) printf("inverting E22 for TREC failed!: exploded!\n");
		else if (r->ll[0] == 2) printf("inverting V for TREC failed!: exploded!\n");		

		if		(r->ll[1] == 4)	printf("evaluate_effiScore_TRECASE failed!: ifault==1!\n");
		else if (r->ll[1] == 3) printf("inverting E22 for ASE failed!: exploded!\n");
		else if (r->ll[1] == 2) printf("inverting V for ASE failed!: exploded!\n");		
	}

}//wrap_TRECASE_score



void wrap_est_TRECASE(RESULTS*r, DATA*d, PARA*p) {

	///////////////////////////////////////////////////////////////
	//  fit alternative model
	///////////////////////////////////////////////////////////////


	get_index(p->ASE, p->TREC, p);

	r->ll[2] = fit_TReCASE_unphased(r->v_vec, r->std_vec, &(r->loglike[2]), *d, *p, r->vASE, r->vTREC);

	if		(r->ll[2] == 5) printf("fitting full TReCASE model failed!: v[TD] < 0.0!\n");
	else if	(r->ll[2] == 4) printf("fitting full TReCASE model failed!: E exploded during iteration!\n");
	else if (r->ll[2] == 3) printf("fitting full TReCASE model failed!: halving algorithm cannot increase log-likelihood!\n");
	else if (r->ll[2] == 2) printf("fitting full TReCASE model failed!: Maximum iteration is reached!\n");
	else if (r->ll[2] == 1) printf("fitting full TReCASE model failed!: EM algorithm converged but Louis estimator is not positive-definite!\n");
	else if (r->ll[2] == 0) {
		r->v[2] = r->v_vec[0];
		r->std[2] = r->std_vec[0];
	}

}//wrap_est_TRECASE



void wrap_test_cis(RESULTS*r, DATA*d, PARA*p) {

	//ll_null[0,1,2]: TREC null, TRECASE null, TRECASE null
	//loglike_null[0,1,2]: TREC null, ASE null, TRECASE null
	//ll[0,1,2]: score for TREC; score for TRECASE; score for trans


	//effi score for trans 

	if (r->ll[2]==0) {

		//r->v[2], r->std[2]: actually 'score' and 'std of score'

		evaluate_effiScore_trans(&r->ll[2], &r->Ucis, &r->Vcis, r->U[2], *d, *p, r->v_vec);

		if (r->ll[2] == 4)		printf("evaluate_effiScore_trans failed!: ifault==1!\n");
		else if (r->ll[2] == 3) printf("inverting E22 for trans failed!: exploded!\n");
		else if (r->ll[2] == 2) printf("inverting V for trans failed!: exploded!\n");		

	}

}//wrap_test_cis




