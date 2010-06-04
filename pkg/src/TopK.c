#include <stdlib.h>
#include <R.h>

double sump(int,double *);
int randomsamp(int, double *);

void topksamplec (double *p, int *k, int *n, int *ns, int *samp, int *seed)
{
  int i=0,j=0,m=0;
  double *p1;  //a vector of probability at one position in a list 
  double sum_p; //sum of remain probabilities at each position in a list
  int tn=*n; //total number of items 
  int tns=*ns; //number of random samples 
  int tk=*k; //length of ordered list in one sample

  p1 = (double *) malloc(sizeof(double)*tn);

  srand(*seed); //random seed
 
  while (i<tns) {
    for (j=0; j<tk; j++) {
      for (m=0;m<tn;m++) {
	p1[m] = p[m +j*tn]; //get a vector of probs from the matrix p
      }
      if (j > 0)
	for (m=0; m<j; m++) {
	  p1[samp[i*tk+m]-1] = 0; //remove probs of previously chosen items
	}
      sum_p = sump(tn,p1); 
      if (sum_p <= 1e-6) break; //remaining probability is too small, disgard current sample
      for (m=0; m<tn; m++) {
	p1[m] /= sum_p; //normalize probs
      }
      samp[i*tk+j] = randomsamp(tn,p1); //choose one item according to the prob vector
    }
    if (j == tk) i++; //increase i only after current sample is complete
  }
  free(p1);
}

double sump(int n, double *p) 
{
  int i=0;
  double sum_p=0;

  for (i=0; i<n; i++) {
    sum_p += p[i];
  }
  return(sum_p);
}

int randomsamp(int n, double *p)
{
  int i=0;
  double *ptemp;
  double rn;
  
  ptemp = (double *) malloc(sizeof(double)*(n+1));
  
  ptemp[0] = 0;
  for (i=1; i<n+1; i++) {
    ptemp[i] = ptemp[i-1]+p[i-1]; //cumulative sums
  }
  ptemp[n] = 1;

  GetRNGstate();
  rn = unif_rand(); //uniform(0,1)
  PutRNGstate();

  for (i=1; i<n+1; i++) {
    if (rn >= ptemp[i-1] && rn <= ptemp[i]) {
      free(ptemp);
      return(i);
    }
  }
  free(ptemp);
  return(-1); //this should not happen since p is normalized earlier
}

void kendallc(int *rank_a, int *rank_b, int *n, int *nb, double *p, double *dist) 
{

  int i,j,k;
  int tn=*n; 
  int tnb=*nb; 

  for (i=0;i<tnb;i++) {
    for (j=0;j<tn-1;j++) {
      for (k=j+1;k<tn;k++) {
	if (rank_a[j]<rank_a[k] && rank_b[tn*i+j]>rank_b[tn*i+k])
	  dist[i] += 1;
	else if (rank_a[j]>rank_a[k] && rank_b[tn*i+j]<rank_b[tn*i+k])
	  dist[i] += 1;
	else if (rank_a[j]==rank_a[k] || rank_b[tn*i+j]==rank_b[tn*i+k])
	  dist[i] += *p;
      }
    }
  }
}
