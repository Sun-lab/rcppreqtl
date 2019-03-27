#include <RcppEigen.h>
static int maxArg1eqtl,maxArg2eqtl;
#define MAX(a,b) (maxArg1eqtl=(a),maxArg2eqtl=(b),(maxArg1eqtl) > (maxArg2eqtl) ?\
	(maxArg1eqtl) : (maxArg2eqtl))


#define ITMAX0 500
#define TOL4 1.0e-4
#define TOL6 1.0e-6
#define TOL8 1.0e-8
#define TOL0 1.0e-10
#define EPS0 1.0e-10
#define EPS6 1.0e-6
#define EPS40 1.0e-40


//////////////////////////////////////////////////////////////////////

void nrerror(const char* error_text);
int    *ivector(int nl, int nh);
char   *cvector(int nl, int nh);
double *dvector(int nl, int nh);
char **cmatrix(int nrl, int nrh, int ncl, int nch);
double **dmatrix(int nrl, int nrh, int ncl, int nch);
int **imatrix(int nrl, int nrh, int ncl, int nch);
double ***dmatrix3(int nrl, int nrh, int ncl, int nch, int n3l, int n3h);
void free_dmatrix3(double ***m, int nrl, int nrh, int ncl, int nch, int n3l, int n3h);
void free_ivector(int *v, int nl, int nh);
void free_cvector(char *v, int nl, int nh);
void free_dvector(double *v, int nl, int nh);
void free_cmatrix(char **m, int nrl, int nrh, int ncl, int nch);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);

//////////////////////////////////////////////////////////////////////////////////////  

double gammq(const double a, const double x);
void gser(double &gamser, const double a, const double x, double &gln);
void gcf(double &gammcf, const double a, const double x, double &gln);
double gammln(const double xx);

////////////////////////////////////////////////////////////////////////////////////// 

void    ludcmp(double **a,int n,int *indx); 
void    lubksb(double **a, int n, int *indx, double *b);
void    invv(double **fma, int nn);   
                                                   
//////////////////////////////////// study link list //////////////////////////////////

void makelist_unlinklist(int**hap1, int**hap2, int*nhap, int*rm, int*Gu, int**H, int n, int M, int*pow2);
void makelist_unlinklist_Gp(int**hap1, int**hap2, int*nhap, int*rm, int n, int**Gp);
void makelist_unlinklist_ASEnull(int**hap1, int**hap2, int*nhap, int*rm, int n);

void reduce_list_unlinklist(int**hap1, int**hap2, int*nhap, int*rm,double thresh,double *pai, int n);


//////////////////////////////////// end of link list //////////////////////////////////

typedef struct {
	int n, n_new;
	int *T;
	int *N, *M, *G, **H, **Gp; //Gp: phased G
	double **X;
	int **hap1, **hap2, *nhap, *rm;
	double **prob, **w, *sum;
	int T_max;
	double*lgamma_int;
} DATA;

typedef struct {
	int M, K;
	int TREC, ASE, TREC_fitted, ASE_fitted;
	int B, R, KK, KK_sim;
	int disp_TREC;
	char* disp_TREC_model;
	int AD, TD, REG, FREQ, nfix, nall;
    double THRESHOLD_L;
	int *h, *hap_freq, *pow2;
	double *pai;
} PARA;

//SNP-level results
typedef struct {
	int itx, iSNP;
	int n_test, nTEST;
	int ll_gene[2]; //initial parameter of 0: TREC; 1: ASE
	double *vTREC, *vASE;
	int ll[3], ll_null[3]; //in order of TREC, ASE, TRECASE
	double v[3], std[3], loglike[3], loglike_null[3];
	double v_vec[100], std_vec[100];
	double ***U, ***V; //in order of TREC, Joint, trans/cis
	double Ucis, Vcis;
	int which_stat_max;
	double stat_mar[2], MD;
} RESULTS;


//---------------------------- interface --------------------------------//


inline bool my_isnan(double x);
double length (double*a,int q);
void reduce_hap(int *hap_left,double *pai,double threshold_F, int K);
int explosion_E(double**E, int p);
int inproper_std(double**E, int p);

void get_initial_para(int*ll_gene, double*vTREC_gene, double*vASE_gene, DATA*d, PARA*p, int TREC_fitted_gene, int ASE_fitted_gene);
void get_null_para_loglike(RESULTS*r, DATA*d, PARA*p);

void Rvec(double*vecR, int R, double*X);
void Bvec(double*vecB, int B, double*X);
void get_index(int ASE, int TReC, PARA*p);
int A1_unlinklist(double*pai, int**hap1, int**hap2, double**prob, int*nhap, int*rm, 
	   int Y1, int Y2, int Y0, int *STATUS, int design, int K, int n);
void get_hap_freq_unlinklist(PARA*p, DATA*d);
int binsearch_int_2(int key, int *str, int max, int left);

double MDfun(double*pai, int M, int K);
int screen_SNP_LD(int*map_rf, int*pwin_r, int nsubj, int tranSNPStart, int tranSNPEnd, int**haplo, double R2_thresh);
int get_tagSNP(int*tagSNP, int*M, double*MDmax, int iSNP, int nsubj, int*map_rf, int pwin_r, int**haplo, int Mmax, int*pow2);
int get_MD(double*MD, int**H, int nsubj, int M, int*pow2);


///////////////////////------- random number generator -------///////////////////////

extern unsigned long seed[];

double kiss(unsigned long *);
int runiform_n_rand(int n,double rand);
void permute_sample_int(int *,int,double *,int*); 


///////////////////////------- core functions -------///////////////////////

int fit_TReCASE_unphased(double*vM, double*stdM, double*loglike, DATA d, PARA p, double*vASE, double*vTREC);
int fit_ASE_null(double*v0, double*std0, double*loglike, DATA d, PARA p) ;
int fit_TReC_null(double*v0, double*std0, double*loglike, DATA d, int R, int disp_TREC, char*disp_TREC_model);
int get_score_disp_TREC(double*score_disp_TREC, double*vTREC, DATA d, int R);
void evaluate_effiScore_TREC_TRECASE(int*ll, int*r_n_stat, double*r_stat_mar, DATA d, PARA p, double***Uei, double*vASE, double*vTREC);
void evaluate_effiScore_trans(int*ll, double*U, double*V, double**Uei, DATA d, PARA p, double*vTRECASE);
void permutation (double stat_max[2][5000], int nPERM, DATA d, PARA p, double*vTREC, double*vASE);
void wrap_TRECASE_score(RESULTS*r, DATA*d, PARA*p);
void wrap_est_TRECASE(RESULTS*r, DATA*d, PARA*p);
void wrap_test_cis(RESULTS*r, DATA*d, PARA*p);


int eQTL_mapping_per_gene(FILE*foutpmin[2], int chr_working, int itx, char**txName, int*txStart, int*txEnd, char**rsSNP, int*posSNP, double*freq_SNP, double MAF_cutoff, double corr2_pruning,
						  int targetSNPStart, int targetSNPEnd, int tranSNPStart, int tranSNPEnd, DATA d, PARA p, RESULTS r, int**haplo);
void eQTL_mapping(FILE*foutput_file[], int start_tx, int end_tx, int chr_working, int local_SNP_win, double corr2_pruning, 
				  int nsubj, int ncov, double**cov,
				  int ntx, char**txName, int*txStart, int*txEnd,
				  int**toCount, int**asCountA, int**asCountB, int T_max,
				  int nSNP, char**rsSNP, int*posSNP, double*freqSNP, double MAF_cutoff,
				  int**haplo); 

// gamma function 
double digama ( double x, int *ifault );
void psi_values ( int *n_data, double *x, double *fx );

double trigamma ( double x, int *ifault );
void trigamma_values ( int *n_data, double *x, double *fx );

double lgamma(double x);

