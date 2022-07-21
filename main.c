# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////// System Parameters //////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Constants
const double alpha = 1.0;                   // diffusivity
const double pi    = M_PI;                  // Pi

// Creating pointers to files for storing data
FILE *fp   = NULL;
char syspar[50];
FILE *fpp  = NULL;
char datadump[50];
FILE *fppp = NULL;
char error_l2[50];

int numfile    = 0;                         // counter to keep track of the number of datafiles created
const int dump = 75;                        // data will be wriiten to a datafile after every 'dump' time steps

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////// Structs ///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Contains information about the grid
struct Grid {
  int     xpts;
  double  xmin;
  double  xmax;
  double *xvec;
  int     ypts;
  double  ymin;
  double  ymax;
  double *yvec;
  int     tpts;
  double  tmin;
  double  tmax;
  double *tvec;
};

typedef struct Grid Grid;

struct Field {
  int      xpts;
  int      ypts;
  double  *data;
  double **f;
};

typedef struct Field Field;

struct Solver {
  Field *k1;
  Field *k2;
  Field *k3;
  Field *k4;
  Field *tmp;
};

typedef struct Solver Solver;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// Function Declarations ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int     main (int argc, char* argv[]);
Grid   *newGrid(int xpts, double xmin, double xmax, int ypts, double ymin, double ymax, int tpts, double tmin, double tmax);
Field  *newField(Grid *grid);
Solver *newSolver(Grid *grid);
Field  *init(Field *field, Grid *grid);
Field  *step(Field *field, Solver *rk4, int k, Grid *grid);
Field  *rhs(Field *dudt, double t, Field *u, Grid *grid);
Field  *addField(Field *out, double a, Field *x, double b, Field *y);
Field  *analyticalSolution(Field *field, int nterms, int k, Grid *grid);
Field  *err(Field *e, double *evec, Field *u_num, Field *u_ana, int k, Grid *grid);
double  l2(Field *field, Grid *grid);
void    delGrid(Grid *grid);
void    delField(Field *field);
void    delSolver(Solver *solver);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char* argv[])
{
  // Starting timer
  clock_t comp_time;
  comp_time = clock();

  // Grid parameters
  const int xpts    = 51; 
  const double xmin = 0.0;
  const double xmax = 5.0;
  const double dx   = (xmax-xmin)/(xpts-1);
  const int ypts    = 51;
  const double ymin = 0.0;
  const double ymax = 5.0;
  const double dy   = (ymax-ymin)/(ypts-1);
  const int tpts    = 1501;
  const double tmin = 0.0;
  const double tmax = 5.0;
  const double dt   = (tmax-tmin)/(tpts-1);

  const int nterms = 100;                            // # of terms in the analytical series solution
  
  // The file 'sysparams.txt' contains the grid parameters
  sprintf(syspar, "./sysparams.txt");
  fp = fopen(syspar,"w+");
  fprintf(fp, "xpts\txmin\txmax\typts\tymin\tymax\ttpts\ttmin\ttmax\n");
  fprintf(fp, "%d\t%f\t%f\t%d\t%f\t%f\t%d\t%f\t%f\n",xpts,xmin,xmax,ypts,ymin,ymax,tpts,tmin,tmax);
  fclose(fp);

  // Allocating memory
  Grid   *G     = newGrid(xpts, xmin, xmax, ypts, ymin, ymax, tpts, tmin, tmax);
  Field  *u_num = newField(G);
  Field  *u_ana = newField(G);
  Field  *e     = newField(G);
  double *evec  = malloc(tpts*sizeof(double));
  Solver *rk4   = newSolver(G);
  
  // The file 'error.dat' contains the absolute error between the numerical solution and the analytical solution
  sprintf(error_l2, "./results/error.dat");
  fppp = fopen(error_l2,"w+");
  fprintf(fppp, "%f\n%f\n",dx,dy);
  fprintf(fppp, "t\te(t)\n");

  // t = 0
  u_num = init(u_num,G);
  u_ana = analyticalSolution(u_ana,nterms,0,G);
  e     = err(e,evec,u_num,u_ana,0,G);

  // Propagating the solution forward in time for t > 0
  for (int k=1; k<tpts; k++) {
    u_num = step(u_num,rk4,k,G);
    u_ana = analyticalSolution(u_ana,nterms,k,G);
    e     = err(e,evec,u_num,u_ana,k,G);
  }
  
  // Ending timer
  comp_time = clock() - comp_time;
  double time_taken = ((double)comp_time)/CLOCKS_PER_SEC;
  printf("Time taken: %f s \n", time_taken);

  // De-allocating memory
  delGrid(G);
  delField(u_num);
  delField(u_ana);
  delField(e);
  delSolver(rk4);
  free(evec);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Grid *newGrid(int xpts, double xmin, double xmax, int ypts, double ymin, double ymax, int tpts, double tmin, double tmax)
/* Creates a new grid object */
{
  
  Grid *grid = malloc(sizeof(Grid));

  // x-coordinate
  grid->xpts = xpts;
  grid->xmin = xmin;
  grid->xmax = xmax;
  grid->xvec = malloc(xpts*sizeof(double));
  double dx  = (xmax-xmin)/(xpts-1);
  for (int i=0; i<xpts; i++) {
    grid->xvec[i] = xmin+i*dx;
  }
  
  // y-coordinate
  grid->ypts = ypts;
  grid->ymin = ymin;
  grid->ymax = ymax;
  grid->yvec = malloc(ypts*sizeof(double));
  double dy  = (ymax-ymin)/(ypts-1);
  for (int i=0; i<ypts; i++) {
    grid->yvec[i] = ymin+i*dy;
  }

  // time
  grid->tpts = tpts;
  grid->tmin = tmin;
  grid->tmax = tmax;
  grid->tvec = malloc(tpts*sizeof(double));
  double dt  = (tmax-tmin)/(tpts-1);
  for (int i=0; i<tpts; i++) {
    grid->tvec[i] = tmin+i*dt;
  }

  return grid;
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void delGrid(Grid *grid)
/* Deallocates memory associated with a grid object */
{
  
  free(grid->xvec);
  free(grid->yvec);
  free(grid->tvec);
  free(grid);
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Field *newField(Grid *grid)
/* Creates a new field object */
{
  
  Field *field = malloc(sizeof(Field));

  field->xpts = grid->xpts;
  field->ypts = grid->ypts;
  field->data = malloc(field->xpts*field->ypts*sizeof(double));
  field->f    = malloc(field->xpts*sizeof(double *));
  for (int i=0; i<grid->xpts; i++) {
    field->f[i] = field->data+i*field->ypts;
  }

  return field;
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void delField(Field *field)
/* Deallocates memory associated with a field object */
{
  
  free(field->data);
  free(field->f);
  free(field);
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Solver *newSolver(Grid *grid)
/* Creates a new solver object */
{
  
  Solver *solver = malloc(sizeof(Solver));

  solver->k1  = newField(grid);
  solver->k2  = newField(grid);
  solver->k3  = newField(grid);
  solver->k4  = newField(grid);
  solver->tmp = newField(grid);

  return solver;
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void delSolver(Solver *solver)
/* Deallocates memory associated with a solver object */
{
  
  delField(solver->k1);
  delField(solver->k2);
  delField(solver->k3);
  delField(solver->k4);
  delField(solver->tmp);
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Field *init(Field *field, Grid *grid)
/* Specifies the initial condition for the numerical simulation */
{
  
  double x,xmax,y,ymax;

  xmax = grid->xmax;
  ymax = grid->ymax;

  const double c = 16/(xmax*xmax*ymax*ymax);                     // normalizing constant for the initial condition so that the max amplitude is unity
  /*
  // storing the initial data in a text file
  printf("t = %f\n",0.0);
  sprintf(datadump, "./results/t%d.txt", numfile);
  fpp = fopen(datadump,"w+");
  fprintf(fpp, "%f\n",0.0);
  fprintf(fpp, "x\ty\tu(t,x,y)\n");
  numfile += 1;
  */
  for (int i=0; i<field->xpts; i++) {
    x = grid->xvec[i];
    for (int j=0; j<field->ypts; j++) {
      y = grid->yvec[j];
      field->f[i][j] = c*x*(xmax-x)*y*(ymax-y);
      //fprintf(fpp, "%f\t%f\t%f\n",x,y,field->f[i][j]);    
    }
  }
  
  return field;
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Field *rhs(Field *dudt, double t, Field *u, Grid *grid)
/* Computes the right hand side of the governing equation */
{
  
  double dx = (grid->xmax-grid->xmin)/(grid->xpts-1);
  double dy = (grid->ymax-grid->ymin)/(grid->ypts-1);
  
  // interior points
  for (int i=1; i<u->xpts-1; i++) {
    for (int j=1; j<u->ypts-1; j++) {
      dudt->f[i][j] = (alpha*alpha)*((u->f[i+1][j]-2*u->f[i][j]+u->f[i-1][j])/(dx*dx)+(u->f[i][j+1]-2*u->f[i][j]+u->f[i][j-1])/(dy*dy));
    }
  }
  
  return dudt;
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Field *addField(Field *out, double a, Field *x, double b, Field *y)
/* 
out = a*x+b*y
Here:
a is a scalar
x is a 2D array
b is a scalar
y is a 2D array
*/
{
  
  for (int i=0; i<out->xpts; i++) {
    for (int j=0; j<out->ypts; j++) {
      out->f[i][j] = a*x->f[i][j]+b*y->f[i][j];
    }
  }
  
  return out;
  
}
			  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Field *step(Field *field, Solver *rk4, int k, Grid *grid)
/* Steps the numerical solution forward by one time step */
{
  
  int    xpts,ypts,tpts;
  double x,xmin,xmax,y,ymin,ymax,t,tmin,tmax,dt;
  
  xpts = grid->xpts;
  xmin = grid->xmin;
  xmax = grid->xmax;
  ypts = grid->ypts;
  ymin = grid->ymin;
  ymax = grid->ymax;
  tpts = grid->tpts;
  tmin = grid->tmin;
  tmax = grid->tmax;
  
  t    = grid->tvec[k];
  dt   = (tmax-tmin)/(tpts-1);
  /*
  // Creating the file that will contain the data for this particular time step
  if (k%dump==0) {
    // printf("k = %d\n",k);
    printf("t = %f\n",t);
    sprintf(datadump, "./results/t%d.txt", numfile);
    fpp = fopen(datadump,"w+");
    fprintf(fpp, "%f\n",t);
    fprintf(fpp, "x\ty\tu(t,x,y)\n");
    numfile += 1;
  }
  */
  rk4->k1  = rhs(rk4->k1,t,field,grid);
  rk4->tmp = addField(rk4->tmp,1.0,field,0.5*dt,rk4->k1);
  rk4->k2  = rhs(rk4->k2,t+(0.5*dt),rk4->tmp,grid);
  rk4->tmp = addField(rk4->tmp,1.0,field,0.5*dt,rk4->k2);
  rk4->k3  = rhs(rk4->k3,t+(0.5*dt),rk4->tmp,grid);
  rk4->tmp = addField(rk4->tmp,1.0,field,dt,rk4->k3);
  rk4->k4  = rhs(rk4->k4,t+dt,rk4->tmp,grid);
  /*
  for (int i=0; i<(xpts); i++) {
    x = grid->xvec[i];
    for (int j=0; j<(ypts); j++) {
      y = grid->yvec[j];
      // interior points
      if ((i!=0 && i!=(xpts-1)) && (j!=0 && j!=(ypts-1))) {
        field->f[i][j] = field->f[i][j]+(dt/6)*(rk4->k1->f[i][j]+2*rk4->k2->f[i][j]+2*rk4->k3->f[i][j]+rk4->k4->f[i][j]);
      }
      // boundary conditions
      else {
        field->f[i][j] = 0;
      }       
      
      // Writing the data to a file
      if (k%dump==0) {
        fprintf(fpp, "%f\t%f\t%f\n",x,y,field->f[i][j]);
      }
        
    }
  }
  */
  // interior points
  for (int i=1; i<(xpts-1); i++) {
    x = grid->xvec[i];
    for (int j=1; j<(ypts-1); j++) {
      y = grid->yvec[j];
      field->f[i][j] = field->f[i][j]+(dt/6)*(rk4->k1->f[i][j]+2*rk4->k2->f[i][j]+2*rk4->k3->f[i][j]+rk4->k4->f[i][j]);
      /*
      if (k%dump==0) {
        fprintf(fpp, "%f\t%f\t%f\n",x,y,field->f[i][j]);
      }
      */
    }
  }
  // xmin boundary
  x = grid->xvec[0];
  for (int j=0; j<(ypts); j++) {
    y = grid->yvec[j];
    field->f[0][j] = 0;
    /*
    if (k%dump==0) {
      fprintf(fpp, "%f\t%f\t%f\n",x,y,field->f[0][j]);
    }
    */
  }
  // xmax boundary
  x = grid->xvec[xpts-1];
  for (int j=0; j<(ypts); j++) {
    y = grid->yvec[j];
    field->f[xpts-1][j] = 0;
    /*
    if (k%dump==0) {
      fprintf(fpp, "%f\t%f\t%f\n",x,y,field->f[xpts-1][j]);
    }
    */
  }
  // ymin boundary
  y = grid->yvec[0];
  for (int i=1; i<(xpts-1); i++) {
    x = grid->xvec[i];
    field->f[i][0] = 0;
    /*
    if (k%dump==0) {
      fprintf(fpp, "%f\t%f\t%f\n",x,y,field->f[i][0]);
    }
    */
  }
  // ymax boundary
  y = grid->yvec[ypts-1];
  for (int i=1; i<(xpts-1); i++) {
    x = grid->xvec[i];
    field->f[i][ypts-1] = 0;
    /*
    if (k%dump==0) {
      fprintf(fpp, "%f\t%f\t%f\n",x,y,field->f[i][ypts-1]);
    }
    */
  }
  
  return field;
  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Field *analyticalSolution(Field *field, int nterms, int k, Grid *grid)
/* Computes the analytical solution */
{

  int    xpts, ypts;
  double x,xmin,xmax,y,ymin,ymax,t,tmp1,tmp2;

  xpts = grid->xpts;
  xmin = grid->xmin;
  xmax = grid->xmax;
  ypts = grid->ypts;
  ymin = grid->ymin;
  ymax = grid->ymax;

  t    = grid->tvec[k];

  tmp1 = 0.0;
  tmp2 = 0.0;
  
  const double c = 16/(xmax*xmax*ymax*ymax);                     // normalizing constant for the initial condition so that the max amplitude is unity
  /*
  // Creating the file that will contain the data for this particular time step
  if (k%dump==0) {
    printf("t = %f\n",t);
    sprintf(datadump, "./results/t%d.txt", numfile);
    fpp = fopen(datadump,"w+");
    fprintf(fpp, "%f\n",t);
    fprintf(fpp, "x\ty\tu(t,x,y)\n");
    numfile += 1;
  }
  */
  for (int i=0; i<xpts; i++) {
    x = grid->xvec[i];
    for (int j=0; j<ypts; j++) {
      y = grid->yvec[j];     
      // series solution
      tmp1 = 0;
      tmp2 = 0;
      for (int m=1; m<nterms; m++) {
	tmp1 += ((pow(-1,m)-1)*exp(-alpha*alpha*t*((m*pi)/xmax)*((m*pi)/xmax))*sin((m*pi*x)/xmax))/(m*m*m);
      }
      for (int n=1; n<nterms; n++) {
	tmp2 += ((pow(-1,n)-1)*exp(-alpha*alpha*t*((n*pi)/ymax)*((n*pi)/ymax))*sin((n*pi*y)/ymax))/(n*n*n);
      }
      field->f[i][j] = ((16*xmax*xmax*ymax*ymax)/pow(pi,6))*c*tmp1*tmp2;      
      /*
      // Writing the data to a file
      if (k%dump==0) {
        fprintf(fpp, "%f\t%f\t%f\n",x,y,field->f[i][j]);
      }
      */
    }
  }
  
  return field;

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Field *err(Field *e, double *evec, Field *u_num, Field *u_ana, int k, Grid *grid)
/* Computes the error between the analytical solution and the numerical solution */
{
  int    xpts,ypts;
  double x,y,t;

  xpts = grid->xpts;
  ypts = grid->ypts;
  
  t    = grid->tvec[k];

  // Creating the file that will contain the data for this particular time step
  if (k%dump==0) {
    // printf("k = %d\n",k);
    printf("t = %f\n",t);
    sprintf(datadump, "./results/t%d.txt", numfile);
    fpp = fopen(datadump,"w+");
    fprintf(fpp, "%f\n",t);
    fprintf(fpp, "x\ty\te(t,x,y)\n");
    numfile += 1;
  }

  // Computing the absolute point-wise error
  for (int i=0; i<xpts; i++) {
    x = grid->xvec[i];
    for (int j=0; j<ypts; j++) {
      y = grid->yvec[j];
      e->f[i][j] = fabs(u_ana->f[i][j]-u_num->f[i][j]);   
      // Writing the data to a file
      if (k%dump==0) {
        fprintf(fpp, "%f\t%f\t%f\n",x,y,e->f[i][j]);
      }
      
    }
  }

  evec[k] = l2(e,grid);
  fprintf(fppp, "%f\t%f\n",grid->tvec[k],evec[k]);
  
  return e;
  
}

double l2(Field *field, Grid *grid)
/* Computes the L2-norm of an solution array over a 2D grid*/
{

  int    npts;
  double tot_err;

  npts    = grid->xpts*grid->ypts;
  tot_err = 0.0;
  
  double dx = (grid->xmax-grid->xmin)/(grid->xpts-1);
  double dy = (grid->ymax-grid->ymin)/(grid->ypts-1);
  
  for (int i=0; i<grid->xpts; i++) {
    for (int j=0; j<grid->ypts; j++) {
      tot_err += field->f[i][j]*field->f[i][j]*dx*dy;
    }
  }
  
  tot_err = sqrt(tot_err/npts);

  return tot_err;
  
}
