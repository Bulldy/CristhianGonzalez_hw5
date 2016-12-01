#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define Ne 6
#define N 1000000
#define v 5

double random01();
double likelihood(double *tobs,double *tmodel);
void tiempo(double x_epi,double y_epi,double *x_est,double *y_est,double *tt);


int main(){
  /*Initialize random number generator*/
  time_t hoy;
  srand((unsigned) time(&hoy));
  
  /*Initialize arrays to hold parameters*/
  double *x,*y;
  x=malloc(N*sizeof(double));
  y=malloc(N*sizeof(double));
  x[0]=(double)(5);
  y[0]=(double)(0.5);

  /*initialize array to hold observed data*/
  double *tobs;
  tobs=malloc(Ne*sizeof(double));
  tobs[0]=3.12;
  tobs[1]=2.98;
  tobs[2]=2.84;
  tobs[3]=3.26;
  tobs[4]=3.12;
  tobs[5]=2.98;

  /*Initialize arrays with the stations' coordinates*/
  double *x_e,*y_e;
  x_e=malloc(Ne*sizeof(double));
  y_e=malloc(Ne*sizeof(double));
  x_e[0]=3.0;y_e[0]=15.0;
  x_e[1]=4.0;y_e[1]=15.0;
  x_e[2]=5.0;y_e[2]=15.0;
  x_e[3]=3.0;y_e[3]=16.0;
  x_e[4]=4.0;y_e[4]=16.0;
  x_e[5]=5.0;y_e[5]=16.0;
  
  /* Metropolis-Hastings*/
  double *ti,*tnew;
  double xnew,ynew,u,pnew,pi,alpha,delta;
  ti=malloc(Ne*sizeof(double));
  tnew=malloc(Ne*sizeof(double));
  delta=0.001;

  for(int i=0;i<N-1;i++){
    tiempo(x[i],y[i],x_e,y_e,ti);
    xnew=x[i]+(2*random01()-1)*delta;
    ynew=y[i]+(2*random01()-1)*delta;
    tiempo(xnew,ynew,x_e,y_e,tnew);
    
    pi=likelihood(tobs,ti);
    pnew=likelihood(tobs,tnew);
    alpha=exp(pnew-pi);
    
    if(alpha>=1.0){
      x[i+1]=xnew;
      y[i+1]=ynew;
    }
    else{
      u=random01();
      if(u<=alpha){
	x[i+1]=xnew;
	y[i+1]=ynew;
      }
      else{
	x[i+1]=x[i];
	y[i+1]=y[i];
      }
    }
  }
  
  /*Save the Markov chains in .txt files for x and y*/
  FILE *fpx,*fpy;
  fpx=fopen("x.txt","w");
  fpy=fopen("y.txt","w");
  
  for(int i=0;i<N;i++){
    fprintf(fpx,"%f ",x[i]);
    fprintf(fpy,"%f ",y[i]);
  }
  fclose(fpx);
  fclose(fpy);

  return(0);
}


/*Returns a random number between 0 and 1*/
double random01(){
  return (double)rand()/(double)RAND_MAX;
}

/*Returns the chi squared likelihood of the model vs the observed data*/
double likelihood(double *tobs,double *tmodel){
  double sum=0;
  for(int i=0;i<Ne;i++){
    sum=sum+((tobs[i]-tmodel[i])*(tobs[i]-tmodel[i]))/0.01;
  }
  return -0.5*sum;
}

/*Saves the resulting times of the earthquake given the new model*/
void tiempo(double x_epi,double y_epi,double *x_est,double *y_est,double *tt){
  double d;
  for(int i=0;i<Ne;i++){
    d=sqrt((x_est[i]-x_epi)*(x_est[i]-x_epi)+(y_est[i]-y_epi)*(y_est[i]-y_epi));
    tt[i]=d/v;
  }
}
