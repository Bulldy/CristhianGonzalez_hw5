#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 1000000
#define DATA_FILE "lotka_volterra_obs.dat"
#define lines 96
#define Nt 795

void DE_solve(double *x, double *y, double a, double b, double c, double d);
double random01();
double likelihood(int *ind, double *xobs, double *yobs, double *xmodel, double *ymodel);

int main()
{
  /*Initialize random number generator*/
  time_t hoy;
  srand((unsigned) time(&hoy));

  /*Initialize arrays to hold parameters on a random number*/
  double *a,*b,*c,*d;
  a=malloc(N*sizeof(double));
  b=malloc(N*sizeof(double));
  c=malloc(N*sizeof(double));
  d=malloc(N*sizeof(double));
  a[0]=(double)(20);
  b[0]=(double)(100);
  c[0]=(double)(5);
  d[0]=(double)(1);
  
  /*Initialize arrays to hold observed data*/
  FILE  *file;
  double *tt,*xobs,*yobs;
  int *ind;
  tt=malloc(lines*sizeof(double));
  xobs=malloc(lines*sizeof(double));
  yobs=malloc(lines*sizeof(double));
  ind=malloc(lines*sizeof(int));
  
  /*Save observed data from file*/
  file=fopen(DATA_FILE,"r");
  char buff[20];
  fseek(file, SEEK_SET, 0);
  fgets(buff,20,file);
  fgets(buff,20,file);
  fgets(buff,20,file);
  for(int i=0;i<lines;i++){
    fgets(buff,20,file);
    sscanf(buff,"%lf %lf %lf",&tt[i],&xobs[i],&yobs[i]);
    ind[i]=(int)(tt[i]*1000);
  }
  fclose(file);

  /*Metropolis-Hastings for the parameters with delta=0.05*/
  double *xi,*yi,*xnew,*ynew;
  double anew,bnew,cnew,dnew,u,pnew,pi,alpha;
  xi=malloc(Nt*sizeof(double));
  yi=malloc(Nt*sizeof(double));
  xnew=malloc(Nt*sizeof(double));
  ynew=malloc(Nt*sizeof(double));
  xi[0]=15;
  yi[0]=13;
  xnew[0]=15;
  ynew[0]=13;
  double delta=0.01;
  
  for(int i=0;i<N-1;i++){
    
    DE_solve(xi,yi,a[i],b[i],c[i],d[i]);
    anew=a[i]+(2*random01()-1.0)*delta;
    while(anew<=0){
      anew=a[i]+(2*random01()-1.0)*delta;
    }
    bnew=b[i]+(2*random01()-1.0)*delta;
    while(bnew<=0){
      bnew=b[i]+(2*random01()-1.0)*delta;
    }
    cnew=c[i]+(2*random01()-1.0)*delta;
    while(cnew<=0){
      cnew=c[i]+(2*random01()-1.0)*delta;
    }
    dnew=d[i]+(2*random01()-1.0)*delta;
    while(dnew<=0){
      dnew=d[i]+(2*random01()-1.0)*delta;
    }
    DE_solve(xnew,ynew,anew,bnew,cnew,dnew);

    pi=likelihood(ind,xobs,yobs,xi,yi);
    pnew=likelihood(ind,xobs,yobs,xnew,ynew);
    alpha=exp(pnew-pi);
    
    if(alpha>=1.0){
      a[i+1]=anew;
      b[i+1]=bnew;
      c[i+1]=cnew;
      d[i+1]=dnew;
    }
    else{
      u=random01();
      if(u<=alpha){
	a[i+1]=anew;
	b[i+1]=bnew;
	c[i+1]=cnew;
	d[i+1]=dnew;
      }
      else{
	a[i+1]=a[i];
	b[i+1]=b[i];
	c[i+1]=c[i];
	d[i+1]=d[i];
      }
    }
  }

  /*Save the Markov chains in files .txt for every parameter*/
  FILE *fpa,*fpb,*fpc,*fpd;
  fpa=fopen("a.txt","w");
  fpb=fopen("b.txt","w");
  fpc=fopen("c.txt","w");
  fpd=fopen("d.txt","w");
  for(int i=0;i<N;i++){
    fprintf(fpa,"%f ",a[i]);
    fprintf(fpb,"%f ",b[i]);
    fprintf(fpc,"%f ",c[i]);
    fprintf(fpd,"%f ",d[i]);
  }
  fclose(fpa);
  fclose(fpb);
  fclose(fpc);
  fclose(fpd);
      
  return(0);
}


double likelihood(int *ind, double *xobs, double *yobs, double *xmodel, double *ymodel){
  double sum=0;
  int pos;
  for(int i=0;i<lines;i++){
    pos=ind[i]-6;
    sum=sum+((xobs[i]-xmodel[pos])*(xobs[i]-xmodel[pos]))+(yobs[i]-ymodel[pos])*(yobs[i]-ymodel[pos]);
  }
  return -0.5*sum;
}

void DE_solve(double *x, double *y, double a, double b, double c, double d){
  double dt=0.001;
  double x1,x2,x3,x4,y1,y2,y3,y4;
  for(int i=0;i<Nt-1;i++){
    x1=dt*x[i]*(a-b*y[i]);
    y1=-1*dt*y[i]*(c-d*x[i]);
    x2=0.5*dt*x1*(a-b*y1);
    y2=-1*0.5*dt*y1*(c-d*x1);
    x3=0.5*dt*x2*(a-b*y2);
    y3=-1*0.5*dt*y2*(c-d*x2);
    x4=dt*x3*(a-b*y3);
    y4=-1*dt*y3*(c-d*x3);

    x[i+1]=x[i]+(x1+2*x2+2*x3+x4)/6;
    y[i+1]=y[i]+(y1+2*y2+2*y3+y4)/6;
  }

}

double random01()
{
  return (double)rand()/(double)RAND_MAX;
}
