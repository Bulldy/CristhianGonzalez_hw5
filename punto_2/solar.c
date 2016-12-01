#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define g 1.985286119e-29
#define N 500000
#define Np 8

double random01();
void modelo(double logM,double alph,double *log10v,double *logrr);
double likelihood(double *tobs,double *tmodel);

int main(){
  /*Initialize random number generator*/
  time_t hoy;
  srand((unsigned) time(&hoy));

  /*Initialize arrays to hold parameters*/
  double *log10m,*alpha;
  log10m=malloc(N*sizeof(double));
  alpha=malloc(N*sizeof(double));
  log10m[0]=(double)(30);
  alpha[0]=(double)(1);

  /*initialize array to hold observed data*/
  double *logv,*logr;
  logv=malloc(Np*sizeof(double));
  logr=malloc(Np*sizeof(double));
  logv[0]=log10(130.9922493297335);logr[0]=log10(0.3374869911678286);
  logv[1]=log10(54.93795930012001);logr[1]=log10(0.7225561344522969);
  logv[2]=log10(39.56084959527357);logr[2]=log10(1.0009840612227163);
  logv[3]=log10(31.10107728190474);logr[3]=log10(1.3797287907222966);
  logv[4]=log10(7.938687587513353);logr[4]=log10(5.078374222349924);
  logv[5]=log10(4.269567859177788);logr[5]=log10(9.396651032822668);
  logv[6]=log10(1.8736524847011014);logr[6]=log10(20.095695206164145);
  logv[7]=log10(1.3183924711057402);logr[7]=log10(30.027169878698956);

  /*Metropolis-Hastings*/
  double *log10v2i,*log10v2new;
  double log10mnew,alphanew,u,pi,pnew,al,delta;
  log10v2i=malloc(Np*sizeof(double));
  log10v2new=malloc(Np*sizeof(double));
  delta=0.1;

  for(int i=0;i<N-1;i++){
    modelo(log10m[i],alpha[i],log10v2i,logr);
    log10mnew=log10m[i]+(2*random01()-1)*delta;
    alphanew=alpha[i]+(2*random01()-1)*delta;
    modelo(log10mnew,alphanew,log10v2new,logr);
    
    pi=likelihood(logv,log10v2i);
    pnew=likelihood(logv,log10v2new);
    al=exp(pnew-pi);
    
    if(al>=1.0){
      log10m[i+1]=log10mnew;
      alpha[i+1]=alphanew;
    }
    else{
      u=random01();
      if(u<=al){
        log10m[i+1]=log10mnew;
	alpha[i+1]=alphanew;
      }
      else{
        log10m[i+1]=log10m[i];
	alpha[i+1]=alpha[i];
      }
    }
  }

  /*Save the Markov chains in .txt files for both parameters*/
  FILE *fpx,*fpy;
  fpx=fopen("log10M.txt","w");
  fpy=fopen("alpha.txt","w");

  for(int i=0;i<N;i++){
    fprintf(fpx,"%f ",log10m[i]);
    fprintf(fpy,"%f ",alpha[i]);
  }
  fclose(fpx);
  fclose(fpy);

  return(0);
}

/*Returns a random number between 0 and 1*/
double random01(){
  return (double)rand()/(double)RAND_MAX;
}

/*Genera el modelo a partir de los parametros*/
void modelo(double logM,double alph,double *log10v,double *logrr){
  for(int i=0;i<Np;i++){
    log10v[i]=log10(g)+logM+(1-alph)*logrr[i];
  }
}

/*Returns the chi squared likelihood of the model vs the observed data*/
double likelihood(double *tobs,double *tmodel){
  double sum=0;
  for(int i=0;i<Np;i++){
    sum=sum+((tobs[i]-tmodel[i])*(tobs[i]-tmodel[i]));
  }
  return -0.5*sum;
}
