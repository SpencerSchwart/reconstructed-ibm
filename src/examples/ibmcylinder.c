#include "navier-stokes/centered.h"
#include "../ibm.h" 
#include "view.h"

#define L0 15.
#define D 0.5
#define LEVEL 10
#define MIN_LEVEL 5

int Re;
double U0 =  1.;             // inlet velocity
double t_end = 50;
coord ci = {L0/4., L0/2.};     // initial coordinates of cylinder

face vector muv[];

u.n[left] = dirichlet (U0);
u.t[left] = dirichlet (0);
p[left]   = neumann (0);
pf[left]  = neumann (0);

u.n[right] = neumann (0);
u.t[right] = neumann (0);
p[right]   = dirichlet (0);
pf[right]  = dirichlet (0);

u.n[top] = neumann (0);
p[top] = neumann (0);
pf[top] = neumann (0);

u.n[bottom] = neumann (0);
p[bottom] = neumann (0);
pf[bottom] = neumann (0);

int main() {
  size(L0);
  init_grid (1 << (LEVEL - 2));
  mu = muv;
  TOLERANCE = 1.e-6; 
  CFL = 0.5;

  Re = 40;
  run();
}


event init (t = 0) {
  vc.x = 0, vc.y = 0;
  solid (ibm, ibmf, - sq(x - ci.x - vc.x) - sq(y - ci.y - vc.y) + sq(D/2));
  refine (ibm[] < 1 && ibm[] > 0 && level < LEVEL);
  solid (ibm, ibmf, - sq(x - ci.x - vc.x) - sq(y - ci.y - vc.y) + sq(D/2));

  foreach()
    u.x[] = U0 * (1 - ibm[]);
}


event moving_cylinder (i++) {
  solid (ibm, ibmf, - sq(x - ci.x - vc.x) - sq(y - ci.y - vc.y) + sq(D/2));
}


event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*(U0)*(D)/(Re);
   boundary ((scalar *) {muv});
}


event logfile (i++) {
  coord F = ibm_force();
  double CD = (F.x)/(0.5*sq(U0)*(D));
  double CL = (F.y)/(0.5*sq(U0)*(D));
  
  fprintf (stderr, "%d %g %d %d %d %d %d %g %g\n",
          i, t, Re, mgp.i, mgp.nrelax, mgu.i, mgu.nrelax, CD, CL); // 11  
}
/*
event profile (t = t_end) {
  double delta = L0/(pow(2,LEVEL));
  char name[80];

  sprintf (name, "vprofx2-%d", Re); // x = center of cylinder
  FILE * fv1 = fopen(name, "w");
  for(double i = 0; i <= L0; i += delta) {
    foreach_point (ci.x, i) {
      fprintf (fv1, "%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]);
    }
  }
  fflush (fv1);
  fclose (fv1);

  sprintf (name, "vprofx3-%d", Re); // x = outlet
  FILE * fv2 = fopen(name, "w");
  for(double i = 0; i <= L0; i += delta) {
    foreach_point (L0-delta, i) {
      fprintf (fv2, "%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]); 
    }
  }
  fflush (fv2);
  fclose (fv2);

  sprintf (name, "vprofy1-%d", Re); // y = L0/2
  FILE * fv3 = fopen(name, "w");
  for(double i = 0; i <= L0; i += delta) {
    foreach_point (i, ci.y) {
      fprintf (fv3, "%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]);
    }
  }
  fflush (fv3);
  fclose (fv3);
}

event snapshot (t = t_end) {
  scalar omega[];
  vorticity (u, omega);

  char name[80];
  view (fov = 2, tx = -0.25, ty = -0.465,
        width = 3000, height = 1500); 
  sprintf (name, "%g-pressure-%d.png", t, Re);
  clear();
  draw_vof ("ibm", "ibmf", lw = 5, lc = {0,0,0});
  squares ("p", min = -0.5, max = 0.5, map = blue_white_red);
  save (name);

  sprintf (name, "%g-vort-%d.png", t, Re);
  clear();
  squares ("omega", min = -3, max = 3, map = blue_white_red);
  draw_vof ("ibm", "ibmf", lw = 5, lc = {0,0,0});
  save (name);

  sprintf (name, "%g-pressureiso-%d.png", t, Re);
  clear();
  draw_vof ("ibm", "ibmf", lw = 5, lc = {0,0,0});
  isoline ("p", n = 20, min = -0.5, max = 0.5, lc = {0,0,0});
  squares ("p", min = -0.5, max = 0.5, map = blue_white_red);
  save (name);

  sprintf (name, "%g-vortiso-%d.png", t, Re);
  clear();
  isoline ("omega", n = 15, min = -3, max = 3, lc = {0,0,0});
  squares ("omega", min = -3, max = 3, map = blue_white_red);
  draw_vof ("ibm", "ibmf", lw = 5,  lc = {0,0,0});
  save (name);
}
*/


event adapt (i++) {
  adapt_wavelet ({forceTotal, u}, (double[]){3e-3, 3e-3, 3e-3, 3e-3},
		 maxlevel = LEVEL, minlevel = MIN_LEVEL);
}

event stop (t = t_end) {
  static FILE * fp = fopen("perf", "w");
  timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
  fprintf (fp, "%d\t%g\t%d\t%g\n", Re, s.real, i, s.speed);
  fflush (fp);
  return 1;
}

