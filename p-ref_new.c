#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"functions1.h"

//#define PNG

#define mu  0.005  /*coefficient of viscosity*/
#define g  18 /*gravity*/
#define box_x  1.5
#define box_x_r 3.0
#define box_y  1.5
#define box_y_r 1.0


#define mh 5
#define A  2.0
#define alp  1.5
#define delta  0.01
#define kk 20
#define dt 0.001
#define Tmax 1.0
/* #define kmax */



double r_a, r_a2, r_a7, r_s, h, wall_a, wall_b;


double Wa(double *x, double *y);
double w(double r);
double dist(double *x, double *y);
double grad_Wa(double r);
double phi(double s);
double f(double *x, double *u);
double dmin(double *x, double *e);
double fop(double *x, double *bar_x);
double e_vec(double *x, double *e);



typedef struct particles
{
  int N;
  int N0;
  double **X;
  double **V;
  double *rho;
  double *P;
  double *M;
  double tau;
  int k;
} particles_t;


void fill_particles(particles_t *p_pl);
void particles_init(particles_t *p_pl);
void update_position(particles_t *p_pl);
void update_velocity(particles_t *p_pl);
void gnu_particles(FILE *gp, particles_t *p_pl, int k1);
void gnu_particles_init(FILE *gp);



int main(void)
{
  FILE *gp;
  particles_t pl;
  int k, kmax, k1;
  double tau;

  tau = pl.tau = dt;
  kmax = floor(Tmax/tau);

  particles_init(&pl);
  
  fill_particles(&pl);

  gp = popen("gnuplot", "w");

  gnu_particles_init(gp);

  gnu_particles(gp,&pl,0);

  getchar();

  
  
  for(k = 1; k <= kmax; k++)
    {
      pl.k = k;
      update_position(&pl);


      if(k % kk == 0)
  	{
  	  k1 = k/kk;
  	  gnu_particles(gp,&pl,k1);
	  printf("k = %d\n", k1);
  	  /*getchar();*/
  	}
      
      fill_particles(&pl);
      
      update_velocity(&pl);
    }
  pclose(gp);
      
  return 0;
}




/* initial particles */
void particles_init(particles_t *p_pl)
{
  int N, N0;
  double **X, **V;
  double *M, *rho, *P;
  int k, l, I;
  int n1, n2;
  double a, c;
  double m;
  h = 0.1/mh;
  r_a = 4.0*h;
  r_a2 = r_a*r_a;
  r_a7 = r_a*r_a*r_a*r_a*r_a*r_a*r_a;
  r_s = -140/r_a7/M_PI;
  m = h*h;
  /*[a*b]*[c*d] = */
  n1 = 4*mh-2;
  n2 = 8*mh-2;
  a = 2*h;
  c = 2*h;

  p_pl->k = 0;
  p_pl->N = N = (n1 + 1)*(n2 + 1);
  p_pl->N0 = N0 = N+4;
  X = p_pl->X = dmatrix(1, N0, 1, 2);
  V = p_pl->V = dmatrix(1, N, 1, 2);
  rho = p_pl->rho = dvector(1, N);
  P = p_pl->P = dvector(1, N);
  p_pl->M = M = dvector(1, N);


  for(k = 1; k <= N; k++)
    {
      V[k][1] = 0.0;
      V[k][2] = 0.0;
    }


  for(k = 0; k <= n1; k++)
    {
      for(l = 0; l <= n2; l++)
  	{
  	  X[(n1+1)*l+k+1][1] = a + k*h;
  	  X[(n1+1)*l+k+1][2] = c + l*h;
  	}
    }

  /*wall point*/
  X[N+1][1] = 0.0;
  X[N+1][2] = box_y;

  X[N+2][1] = 0.0;
  X[N+2][2] = 0.0;
  
  X[N+3][1] = box_x;
  X[N+3][2] = 0.0;
  
  X[N0][1] = box_x_r;
  X[N0][2] = box_y_r;

  wall_a = (box_y_r-box_y)/(box_x_r-box_x);
  wall_b = box_y-a*box_x;


  for(k = 1; k <= N; k++)
    {
      M[k] = m;
    }
}

void fill_particles(particles_t *p_pl)
{
  int i, j;
  double **X;
  double *rho, *P, *M;
  int N;

  N = p_pl->N;
  X = p_pl->X;
  rho = p_pl->rho;
  P = p_pl->P;
  M = p_pl->M;


  for(i = 1; i <= N; i++)
    {
      rho[i] = 0.0;
      for(j = 1; j <= N; j++)
  	{
  	  rho[i] = rho[i] + M[j]*Wa(X[i], X[j]);
  	}
       P[i] = phi(rho[i]);
    }
 }

void update_position(particles_t *p_pl)
{
  int N;
  int i, j;
  double **X, **V;
  double tau;
  
  N = p_pl->N;
  tau = p_pl->tau;
  X = p_pl->X;
  V = p_pl->V;
    
  for(i = 1; i <= N; i++)
    {
      X[i][1] += tau*V[i][1];
      X[i][2] += tau*V[i][2];
    }
}

void update_velocity(particles_t *p_pl)
{
  int N;
  int i, j;
  double **X, **V;
  double *e;
  double *rho, *P, *M;
  double tau;
  double l1;
  double tmp, r2;
  double d;



  N = p_pl->N;
  tau = p_pl->tau;
  X = p_pl->X;
  V = p_pl->V;
  rho = p_pl->rho;
  P = p_pl->P;
  M = p_pl->M;

  e = dvector(1, 2);
  
  for(i = 1; i <= N; i++)
    {
      for(j = 1; j <= N; j++)
	{
	  r2 = (X[i][1]-X[j][1])*(X[i][1]-X[j][1])
	    + (X[i][2]-X[j][2])*(X[i][2]-X[j][2]);
	  if(r2 > 0.0000001)
	    {
      	      if(r2 < r_a2)
      		{
      		  tmp = tau*M[j]*(mu*((V[i][1]-V[j][1])*(X[i][1]-X[j][1])+
      				      (V[i][2]-V[j][2])*(X[i][2]-X[j][2]))/r2
      				  - (P[i]+P[j])/2)/rho[i]/rho[j]
      		    *grad_Wa(sqrt(r2));
      		  V[i][1] += tmp*(X[i][1] - X[j][1]);
      		  V[i][2] += tmp*(X[i][2] - X[j][2]);		   
      		}
	    }
	}
      d = dmin(X[i],e);
      if(d < delta)
	{
	  l1 = tau*A*(1.0 - (d/delta))*pow(d,-alp);
	  V[i][1] += e[1]*l1;
	  V[i][2] += e[2]*l1;
	}
      V[i][2] -= tau*g;
    }
}

void gnu_particles_init(FILE *gp)
{
  FILE *fp;
  fprintf(gp,"set xrange[-0.2:3.1]\n");
  fprintf(gp,"set yrange[-0.2:1.6]\n");
  fprintf(gp,"unset key\n");

  #ifdef PNG
  {
    fp = popen("mkdir p-ref-pct","w");
    fprintf(gp,"set terminal png\n");
    pclose(fp);
  }
  #else
  {
    fprintf(gp,"set terminal aqua\n");
  }
#endif
}

void gnu_particles(FILE *gp,particles_t *p_pl,int k1)
{
  FILE *fout, *fp;
  int i, j;
  int N;
  int N0;

  N = p_pl->N;
  N0 = p_pl->N0;


  #ifdef PNG

  fprintf(gp,"set output './p-ref-pct/p-ref%05d.png'\n", k1);

  #endif

  fprintf(gp,"set title 'time = %f'\n", p_pl->k*p_pl->tau);

 if((fout = fopen("p-ref.dat","w")) == NULL)
    {
      printf("No file\n");
      exit(1);
    }

  for(i = 1; i <= N; i++)
    {
      fprintf(fout,"%f %f\n", p_pl->X[i][1], p_pl->X[i][2]);
    }
  fprintf(fout,"\n \n");

  for(i = N+1; i <= N0; i++)
    {
      fprintf(fout,"%f %f\n", p_pl->X[i][1], p_pl->X[i][2]);
    }
  fprintf(fout,"\n \n");


  fclose(fout);
  fprintf(gp,"plot'p-ref.dat' index 0 w p pt 3 lt 3,'p-ref.dat' index 1 w l lw 4 lt -1\n");
  fflush(gp);
  /*usleep(200000);*/
  /*pclose(gp);*/
}
  
/* r_a wait function */
double Wa(double *x, double *y)
{
  double r2;

  r2 = (x[1]-y[1])*(x[1]-y[1]) + (x[2]-y[2])*(x[2]-y[2]);

  if(r2 >= r_a2) 
    {
      return 0.0;
    }
  else 
    {
      return w(sqrt(r2)/r_a)/r_a2;
    }
}

/* wait function */
double w(double r)

{
  return 7*(1.0-r)*(1.0-r)*(1.0-r)*(1.0-r)*(1.0+4*r)/M_PI;
}



/* z = grad_Wa(x-y) ; D_Wa/r */
double grad_Wa(double r)
{

  double f;
  f = (r_a-r);

  return r_s*f*f*f;
}


double phi(double s)
{
  double rho0 = 1.0;
  double B = 30;

  return B*(-1.0 + s/rho0);
}


/* foot of a perpendeicular line */
double fop(double *x, double *bar_x)
{
  return bar_x[1] = (x[1]+wall_a*x[2]-wall_a*wall_b)/(1+wall_a*wall_a);
  return bar_x[2] = (wall_a*x[1]+wall_a*wall_a*x[2]+wall_b)/(1+wall_a*wall_a);
}


/* distance wall */
double dmin(double *x, double *e)
{
  if(x[1] <= 0.0 || x[1] >= box_x || x[2] <= 0.0 || x[2] >= box_y)
    {
      printf("error dmin\n");
      return 100000;
    }
  else
    {
      if(x[1] <= Min((box_x-x[1]),x[2]))
	{
	  e[1] = 1.0;
	  e[2] = 0.0;
	  return x[1];
	}
      if(box_x-x[1] <= Min(x[1],x[2]))
	{
	  e[1] = -1.0;
	  e[2] = 0.0;
	  return box_x-x[1];
	}
      else
	{
	  e[1] = 0.0;
	  e[2] = 1.0;
	  return x[2];
	}
    }
}

/* distance */
double dist(double *x, double *y)
{
  return sqrt((x[1]-y[1])*(x[1]-y[1]) + (x[2]-y[2])*(x[2]-y[2]));
}

double e_vec(double*x, double *e)
{
  return 1;
}
  


#include"functions1.c"
