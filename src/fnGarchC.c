/* fnGarchC.c
   David Ardia 
   creation : 20050212 
   updated  : 20051222
   
   this function compute the variance recursively and generate random 
   sample for GARCH(p,q) and GJR(p,q) processes

   inputs :
   n -- length of the data
   order -- order of the GARCH model :
   alpha -- alpha parameters
   delta -- asymmetry parameter
   beta -- garch parameters
   h -- conditionnal variance (length n+1)
   e -- errors (length n)
   sim -- vector of simulation terms (1 or 0, length n+1)

   outputs : 
   h -- conditional variance process
   u -- errors if simulation

   processes :
   GARCH(p,q)
   h_t = alpha_0 + alpha_1 e_t-1^2 + ... + alpha_q e_t-q^2 + beta_1 h_t-1 + ... + beta_p h_t-p

   GJR(p,q)
   h_t = alpha_0 + (alpha_1 I_e>0 + delta_1 I_e<0) e_t-1^2 + ... + (alpha_q I_e>0 + delta_q Ie<0) e_t-q^2 
   + beta_1 h_t-1 + ... + beta_p h_t-p
*/

#include <math.h>
 
void fnGarchC(int *n, int *order, double *alpha, double *delta, 
	      double *beta, double *h, double *u, int *sim){
  int i, j, q, p, s;
  double ht;
  
  /* ordr of the GARCH process */
  q = order[1]-1;
  p = order[0];
  
  /*max between p and q*/
  s = q;
  if (p>q)
    s = p;
  
  /* initialisation o of the
   conditional variance */
  for(i=0;i<s;i++){
    h[i] = alpha[0];
    /* simulation */
    if(sim[i]==1)
      u[i] = sqrt(h[i]) * u[i];
  }
  
  /* arch and garch terms */
  for (i=s;i<=*n;i++){
    ht = alpha[0];
    /* arch effects */
    if(q>0){
      for (j=1;j<=q;j++){
	if (delta[j-1]==0)
	  ht = ht + alpha[j] * u[i-j] * u[i-j];
	if (delta[j-1]>0 && u[i-j]>=0)
	  ht = ht + alpha[j] * u[i-j] * u[i-j];
	if (delta[j-1]>0 && u[i-j]<0)
	  ht = ht + delta[j-1] * u[i-j] * u[i-j];
      }
    }    
    /* garch effects */
    if (p>0){
      for (j=1;j<=p;j++)
        ht = ht + beta[j-1] * h[i-j];
    }
    h[i] = ht;
    /* if simulation */
    if (sim[i]==1)
      u[i] = sqrt(ht) * u[i];
  }
}

/* symmetric filter */
void fnFilterAlphaC(int *n, double *u, double *beta, double *vstar, double *lstar){
  int i;
  /* initialisation i=0 but t=1 */
  vstar[0] = (u[0] * u[0]);
  lstar[0] = 1;
  
  for(i=1;i<*n;i++){
    vstar[i] = (u[i] * u[i]) + *beta * vstar[i-1];
    lstar[i] = 1 + *beta * lstar[i-1];
  }
}  

/* asymmetric filter */
void fnFilterAlphaAsymC(int *n, double *u, double *beta, double *vstar, 
			double *vstarstar, double *lstar){
  int i;
  /* initialisation i=0 but t=1 */
  vstar[0] = (u[0] * u[0]);
  vstarstar[0] = (u[0] * u[0]);
  lstar[0] = 0;
  
  for(i=1;i<*n;i++){
    vstar[i] = *beta * vstar[i-1];
    if(u[i]>=0)
      vstar[i] = vstar[i] + (u[i] * u[i]);
    vstarstar[i] = *beta * vstarstar[i-1];
    if(u[i]<0)
      vstarstar[i] = vstarstar[i] + (u[i] * u[i]);
    lstar[i] = 1 + *beta * lstar[i-1];
  }
}  

/* symmetric w */
void fnFilterWC(int *n, double *u, double *alpha, double *beta, double *w){
  int i;
  /* initialisation i=0 but t=1 */
  w[0] = (u[0] * u[0]) - alpha[0];
  
  for(i=1;i<*n;i++)
    w[i] = (u[i] * u[i]) - alpha[0] - (alpha[1] + *beta) * (u[i-1] * u[i-1]) + *beta * w[i-1];
}

/* asymmetric w */
void fnFilterWAsymC(int *n, double *u, double *alpha, double *beta, double *w){
  int i;
  /* initialisation i=0 but t=1 */
  w[0] = (u[0] * u[0]) - alpha[0];
  
  for(i=1;i<*n;i++){
    w[i] = (u[i] * u[i]) - alpha[0] + *beta * w[i-1];

    /* asymmetries */
    if (u[i-1]>=0) 
      w[i] = w[i] - (alpha[1] + *beta) * (u[i-1] * u[i-1]);
    if (u[i-1]<0)
      w[i] = w[i] - (alpha[2] + *beta) * (u[i-1] * u[i-1]);
  }
}

/* quasi difference between two vectors */
void fnQDiffC(int *n, double *u, double *theta, double *w, double *w0){
  int i;
  w[0] = *w0;
  
  for(i=1;i<*n;i++)
    w[i] = u[i-1] - *theta * w[i-1];
}
