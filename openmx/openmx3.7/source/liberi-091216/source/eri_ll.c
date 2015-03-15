/* ----------------------------------------------------------------------
  eri_ll.c

  Low-level functions of LIBERI

  Coded by TOYODA Masayuki, June 19, 2009
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "eri_def.h"
#include "eri.h"
#include "eri_sf.h"
#include "eri_interpolate.h"
#include "eri_gtbl.h"
#include "sbt/linear/eri_linfsbt.h"
#include "sbt/log/eri_logfsbt.h"


#define LOOP_UNROLLING 1


#define WORKSPACE_AT_STACK 0
#define MAXNUMR 1024

/*--- Switches for speeding-up ---*/
#define SUMSKIP    1
#define SUMSKIP_THRESHOLD 1e-10

/* index table */
static int g_itbl_flag = 0;
static int g_itbl_j2l[ERI_LMAXMAX*ERI_LMAXMAX];
static int g_itbl_j2m[ERI_LMAXMAX*ERI_LMAXMAX];

static char* my_version = "GAMMA-091210";


struct ERI_Struct {
  int lmax;
  int lmax_gl;
  int rcsh;

  /* Gauss-Laguerre quadrature */
  int glq_n;
  double *glq_x;
  double *glq_w;
  double *glq_j;
  double *glq_dj;

  /* SBT */
  int sbttype;
  ERI_SBT_t  *sbt;

  /* Gaunt-coefficients */
  ERI_Gtbl_t *gtbl;

  ERI_Init_Misc info;


#if WORKSPACE_AT_STACK
#else
  double *ws_in;  /* [ERI_NGRIDMAX*2] */
  double *ws_out; /* [ERI_NGRIDMAX*2] */
  double *ws_sb;  /* [(ERI_LMAXMAX*2)*ERI_NGRIDMAX] */
  double *ws_dsb; /* [(ERI_LMAXMAX*2)*ERI_NGRIDMAX] */
#endif
};
 



/*----------------------------------------------------------------------
  INTERNAL SUBROUTINES
----------------------------------------------------------------------*/
static void phase(int l, double *a)
{
  double re = a[0], im = a[1];

  switch ((l%4+4)%4) {
  case 0: break;
  case 1: a[0] = -im; a[1] =  re; break;
  case 2: a[0] = -re; a[1] = -im; break;
  case 3: a[0] =  im; a[1] = -re; break;
  }
}



const char* ERI_Version(void)
{
  return my_version;
}




/*----------------------------------------------------------------------
  ERI_Init
----------------------------------------------------------------------*/
ERI_t* ERI_Init(
  int lmax,
  int lmax_gl,
  int ngrid,
  int ngl,
  int rcsh,
  const ERI_Init_Misc *info
)
{
  int j, l, m, jmax;
  double dt;
  ERI_t *solver;
  
  STEPTRACE("ERI_Init: in");
  
  jmax  = lmax*lmax;

  /* allocation */ 
  solver = (ERI_t*)malloc(sizeof(ERI_t));
  if (NULL==solver) { return NULL; }

  solver->glq_x  = (double*)malloc(sizeof(double)*ngl); 
  solver->glq_w  = (double*)malloc(sizeof(double)*ngl); 
  solver->glq_j  = (double*)malloc(sizeof(double)*ngl*lmax); 
  solver->glq_dj = (double*)malloc(sizeof(double)*ngl*lmax); 
#if WORKSPACE_AT_STACK
#else
  solver->ws_in  = (double*)malloc(sizeof(double)*ngrid*2);
  solver->ws_out = (double*)malloc(sizeof(double)*ngrid*2);
  solver->ws_sb  = (double*)malloc(sizeof(double)*lmax*ngrid*2);
  solver->ws_dsb = (double*)malloc(sizeof(double)*lmax*ngrid*2);
#endif
  if (NULL==solver->glq_x || NULL==solver->glq_w 
    || NULL==solver->glq_j || NULL==solver->glq_dj ) {
    ERI_Free(solver);
    return NULL;
  }
  
  /* j <-> (l,m) conversion table */
  if (0==g_itbl_flag) {
    for (j=0; j<jmax; j++) {
      l = (int)sqrt((double)j);
      m = j-l*(l+1);
      g_itbl_j2l[j] = l;
      g_itbl_j2m[j] = m;
    }
  }
  
  /* abscissas and weights for Gauss-Laguerre integration */
  ERI_GLQ(solver->glq_x, solver->glq_w, ngl);

  /* Gaunt coefficients */
  solver->gtbl = ERI_Gtbl_Init(lmax, rcsh);

  /* additional parameters */
  if (NULL==info) {
    /* default values */
    solver->info.sbttype = ERI_SBT_LINEAR;
    solver->info.rmax    = 100.0;
    solver->info.rho0    = -10.0;
    solver->info.qmin    = 1e-4;
    solver->info.qmax    = 1e-2;
    solver->info.nq      = 3; 
  } else {
    solver->info.sbttype = info->sbttype;
    solver->info.rmax    = info->rmax;
    solver->info.rho0    = info->rho0;
    solver->info.qmin    = info->qmin;
    solver->info.qmax    = info->qmax;
    solver->info.nq      = info->nq; 
  }
  dt = 2.0*PI/(log(solver->info.rmax)-(solver->info.rho0));

  switch (solver->info.sbttype) {
  case ERI_SBT_LINEAR:
    solver->sbt = ERI_LinFSBT_Init(ngrid, lmax, solver->info.rmax);
    break;
  case ERI_SBT_LOG:
    solver->sbt = ERI_LogFSBT_Init(lmax, ngrid, solver->info.rho0, dt, 
      solver->info.qmin, solver->info.qmax, solver->info.nq);
    break;
  }
  if (NULL==solver->sbt) {
    ERI_Free(solver);
    return NULL;
  }

  solver->lmax    = lmax;
  solver->lmax_gl = lmax_gl;
  solver->glq_n   = ngl;
  solver->rcsh    = rcsh;

  return solver;
}




/*----------------------------------------------------------------------
  ERI_Required_Size
----------------------------------------------------------------------*/
size_t ERI_Required_Size(
  int lmax,
  int lmax_gl,
  int ngrid,
  int ngl,
  int rcsh,
  const ERI_Init_Misc *info
)
{
  int nq, type, jmax;
  size_t sz;

  jmax = lmax*lmax;
 
  sz = sizeof(ERI_t)
       + ERI_Gtbl_Required_Size(lmax, 0)
       + sizeof(double)*ngl            /* glx */
       + sizeof(double)*ngl            /* glw */
       + sizeof(double)*ngl*lmax       /* glq_j */
       + sizeof(double)*ngl*lmax       /* glq_dj */
    ;

  /* Gaunt coefficients */
  switch (rcsh) {
  case ERI_SH_COMPLEX:
    sz += ERI_Gtbl_Required_Size(lmax, 0);
    break;
  case ERI_SH_REAL:
    sz += ERI_Gtbl_Required_Size(lmax, 1);
    break;
  default:
    fprintf(stderr, "*** undefined rcsh value(%d)\n", rcsh);
    abort();
  }

  type = info ? info->sbttype : ERI_SBT_LINEAR;
  switch (type) {
  case ERI_SBT_LOG:
    if (NULL==info) {
      nq = 3;
    } else {
      nq = info->nq; 
    }
    sz += ERI_LogFSBT_Required_Size(lmax, ngrid, nq);
    break;
  case ERI_SBT_LINEAR:
    sz += ERI_LinFSBT_Required_Size(ngrid, lmax);
    break;
  default:
    fprintf(stderr, "*** undefined type value(%d)\n", type);
    abort();
  }

  return sz;
}




/*----------------------------------------------------------------------
  ERI_Free
----------------------------------------------------------------------*/
void ERI_Free(ERI_t *solver)
{
  if (solver) {
    if (solver->glq_x)  { free(solver->glq_x); }
    if (solver->glq_w)  { free(solver->glq_w); }
    if (solver->glq_j)  { free(solver->glq_j); }
    if (solver->sbt)    { ERI_SBT_Free(solver->sbt); }
    if (solver->gtbl)   { ERI_Gtbl_Free(solver->gtbl); }
    if (solver->ws_in)  { free(solver->ws_in); }
    if (solver->ws_out) { free(solver->ws_out); }
    if (solver->ws_sb)  { free(solver->ws_sb); }
    if (solver->ws_dsb) { free(solver->ws_dsb); }
    free(solver); 
  }
}



int ERI_lmax(const ERI_t *solver)  { return solver->lmax; }
int ERI_lmax_gl(const ERI_t *solver) { return solver->lmax_gl; }
int ERI_ngl(const ERI_t *solver)   { return solver->glq_n; }
int ERI_rcsh(const ERI_t *solver) { return solver->rcsh; }

const double* ERI_Mesh_Array_glx(const ERI_t *solver) 
{ 
  return solver->glq_x; 
}

const double* ERI_Mesh_Array_glw(const ERI_t *solver) 
{ 
  return solver->glq_w; 
}

int ERI_ngrid(const ERI_t *solver) 
{ 
  return ERI_SBT_ngrid(solver->sbt);
}

double ERI_Mesh_r(const ERI_t *solver, int i) 
{ 
  return ERI_SBT_Mesh_r(solver->sbt, i);
}


double ERI_Mesh_k(const ERI_t *solver, int i)
{ 
  return ERI_SBT_Mesh_k(solver->sbt, i);
}


const double* ERI_Mesh_Array_r(const ERI_t *solver)
{ 
  return ERI_SBT_Mesh_Array_r(solver->sbt);
}


const double* ERI_Mesh_Array_k(const ERI_t *solver)
{ 
  return ERI_SBT_Mesh_Array_k(solver->sbt);
}


double ERI_Mesh_dr(const ERI_t *solver, int i) 
{ 
  return ERI_SBT_Mesh_dr(solver->sbt, i);
}


double ERI_Mesh_dk(const ERI_t *solver, int i)
{ 
  return ERI_SBT_Mesh_dk(solver->sbt, i);
}



/*----------------------------------------------------------------------
  SIZE Calculators  
----------------------------------------------------------------------*/
size_t ERI_Size_of_Orbital(const ERI_t *solver) 
{
  int ngrid = ERI_ngrid(solver);
  return sizeof(double)*ngrid;
}


size_t ERI_Size_of_Gamma(const ERI_t *solver)
{
  int ngrid = ERI_ngrid(solver);
  int lmax  = ERI_lmax(solver);
  int jmax  = lmax*lmax;

#if SUMSKIP
  return sizeof(double)*(ngrid*jmax + jmax);
#else
  return sizeof(double)*(ngrid*jmax);
#endif
}


size_t ERI_Size_of_Alpha(const ERI_t *solver)
{
  int ndalp;
  int ngrid = ERI_ngrid(solver);
  int lmax  = ERI_lmax(solver);
  int jmax  = lmax*lmax;

  switch (solver->rcsh) {
  case ERI_SH_REAL:    ndalp = ngrid*jmax; break;
  case ERI_SH_COMPLEX: ndalp = ngrid*jmax*2; break;
  }
 
#if SUMSKIP
  return sizeof(double)*(ndalp + jmax);
#else
  return sizeof(double)*(ndalp);
#endif
}


size_t ERI_Size_of_Overlap(const ERI_t *solver)
{
  int ndp;
  int ngrid = ERI_ngrid(solver);
  int lmax0 = ERI_lmax_gl(solver);
  int jmax0 = lmax0*lmax0;

  switch (solver->rcsh) {
  case ERI_SH_COMPLEX: ndp = 2*ngrid*jmax0; break;
  case ERI_SH_REAL:    ndp = ngrid*jmax0; break;
  }
  
#if SUMSKIP
  return sizeof(double)*(ndp + jmax0); 
#else
  return sizeof(double)*ndp;
#endif
}


size_t ERI_Size_of_GLF(const ERI_t *solver)
{
  int ndglf;
  int ngl   = ERI_ngl(solver);
  int lmax0 = ERI_lmax_gl(solver);
  int jmax0 = lmax0*lmax0;
  
  switch (solver->rcsh) {
  case ERI_SH_COMPLEX: ndglf = 2*ngl*jmax0; break;
  case ERI_SH_REAL:    ndglf = ngl*jmax0; break;
  }

#if SUMSKIP 
  return sizeof(double)*(ndglf + jmax0); 
#else
  return sizeof(double)*ndglf;
#endif
}




/*----------------------------------------------------------------------
  ERI_Orbital

  For a given orbital function, this performs sph. Bessel transform.
----------------------------------------------------------------------*/
void ERI_Transform_Orbital(
  ERI_t        *solver,
  double       *fk,     /* (OUT) [ngrid] transformed orbital */
  const double *fr,     /* (IN)  [ngrid] radial orbital function */
  const double *xr,     /* (IN)  [ngrid] radial grid for fr */
  int           n,      /* (IN)  number of radial grid points */
  int           l       /* (IN)  angular momentum */
)
{
  int i, ngrid;
  const double *rmesh;
  ERI_CSpline_t *cs;

#if WORKSPACE_AT_STACK
  double in[ERI_NGRIDMAX*2], out[ERI_NGRIDMAX*2];
#else
  double *in = solver->ws_in;
  double *out = solver->ws_out;
#endif
  
  ngrid = ERI_ngrid(solver);
  rmesh = ERI_Mesh_Array_r(solver);

  /* change radial grid */ 
  cs = ERI_CSpline_Init(xr, fr, n);
  for (i=0; i<ngrid; i++) {
    in[2*i+0] = ERI_CSpline_Eval(rmesh[i], cs);
    in[2*i+1] = 0.0;
  }
  ERI_CSpline_Free(cs);
  
  /* transform */
  ERI_SBT_Transform(solver->sbt, out, in, l, ERI_SBT_FORWARD);

  for (i=0; i<ngrid; i++) {
    fk[i] = out[2*i+0]; 
    /* imag part out[2*i+1] is thrown away */
  }
}



/*----------------------------------------------------------------------
  ERI_Gamma
  g  : Gamma factor (out)
  dg : Deribative of Gamma factor (out)
  fk : Transformed wave function 
  R  : Displacement
----------------------------------------------------------------------*/
void ERI_LL_Gamma(
  ERI_t        *solver,
  double       *g,
  double       *dg,
  const double *fk,
  double        rd
) 
{  
  int ngrid, lmax, jmax;
  int l1, l2, ik, ig, igbase, j;
  const double *kmesh;

#if SUMSKIP
  double *gmax  = NULL;
  double *dgmax = NULL;
#endif

#if WORKSPACE_AT_STACK
  double sb[(ERI_LMAXMAX*2)*ERI_NGRIDMAX];
  double dsb[(ERI_LMAXMAX*2)*ERI_NGRIDMAX];
  double in[ERI_NGRIDMAX*2];
  double out[ERI_NGRIDMAX*2];
#else
  double *sb  = solver->ws_sb;
  double *dsb = solver->ws_dsb;
  double *in  = solver->ws_in;
  double *out = solver->ws_out;
#endif

  ngrid = ERI_ngrid(solver);
  lmax  = ERI_lmax(solver);
  jmax  = lmax*lmax;
  kmesh = ERI_Mesh_Array_k(solver);

#if SUMSKIP
  gmax  =  &g[ngrid*jmax];
  for (j=0; j<jmax; j++) { gmax[j] = 0.0; }
  if (dg) {
    dgmax =  &dg[ngrid*jmax];
    for (j=0; j<jmax; j++) { dgmax[j] = 0.0; }
  } 
#endif

  /*--- the spherical Bessel functions ---*/  
  for (ik=0; ik<ngrid; ik++) {
    ERI_Spherical_Bessel(kmesh[ik]*rd, lmax, &sb[ik*lmax], &dsb[ik*lmax]);
  }

  for (l2=0; l2<lmax; l2++) {
    for (ik=0; ik<ngrid; ik++) {
      in[2*ik+0] = fk[ik]*sb[ik*lmax+l2];
      in[2*ik+1] = 0.0;
    }
    ERI_SBT_Transform_Input(solver->sbt, in, ERI_SBT_BACKWARD);
    for (l1=0; l1<lmax; l1++) {
      j = l1*lmax+l2;
      igbase = j*ngrid;
      ERI_SBT_Transform_Output(solver->sbt, out, l1);
      for (ik=0; ik<ngrid; ik++) {
        ig = igbase+ik;
        g[ig] = out[2*ik+0]*2.0/PI;
        /* imag part out[2*ik+1] is thworn away */
#if SUMSKIP
        if (fabs(g[ig])>gmax[j]) { gmax[j] = fabs(g[ig]); }
#endif
      } /* loop of ik */
    } /* loop of l1 */
  } /* loop of l2 */

  if (dg) {
    /*--- deribative ---*/
    for (l2=0; l2<lmax; l2++) {
      for (ik=0; ik<ngrid; ik++) {
        in[2*ik+0] = fk[ik]*dsb[ik*lmax+l2]*kmesh[ik];
        in[2*ik+1] = 0.0;
      }
      ERI_SBT_Transform_Input(solver->sbt, in, ERI_SBT_BACKWARD);
      for (l1=0; l1<lmax; l1++) {
        j = l1*lmax+l2;
        igbase = j*ngrid;
        ERI_SBT_Transform_Output(solver->sbt, out, l1);
        for (ik=0; ik<ngrid; ik++) {
          ig = igbase+ik;
          dg[ig] = out[2*ik+0]*2.0/PI;
          /* imag part out[2*ik+1] is thrown away */
#if SUMSKIP
          if (fabs(dg[ig])>dgmax[j]) { dgmax[j] = fabs(dg[ig]); }
#endif
        } /* loop of ik */
      } /* loop of l1 */
    } /* loop of l2 */
  } /* end if */
}




/*----------------------------------------------------------------------
  Alpha term 
 
  Alpha term.
  If you need the ferivatives as well, you should call ERI_Alpha_d.

  This routine uses the complex spherical harmonic functions.

  alpha[l1][m2](r) = 4 pi sum[l2][m2] i^(l-l1+l2) G(l,m,l1,m1,l2,m2) 
                                        Y[l2][m2](anlge R) g[l1][l2](r)
----------------------------------------------------------------------*/
void ERI_LL_Alpha(
  ERI_t        *solver,
  double       *alp,   /* (OUT) Alpha term*/
  const double *g,     /* (IN) Gamma term */ 
  double        theta, /*      in the spherical coordinates. */
  double        phi,      
  int           l,     /* (IN) Angular momentum of the orbital */
  int           m
)
{
  int ngrid, lmax, jmax;
  int j, j1, j2, ir, l1, m1, l2, m2;
  int ia, ia_base, ig, ig_base;
  double sh[2], dsht[2], dshp[2], fpg, fpgsh[2];
  
  const int *j2l   = g_itbl_j2l;
  const int *j2m   = g_itbl_j2m;

  int ndalp;

  int i, n, igtbl;
  const int *gtbl_n;
  const long *gtbl_j1j2;
  const double *gtbl_gc;
  
#if SUMSKIP
  double *amax = NULL;
  const double *gmax = NULL;
#endif
 
  ngrid = ERI_ngrid(solver);
  lmax  = ERI_lmax(solver);
  jmax  = lmax*lmax;

  gtbl_n    = ERI_Gtbl_n(solver->gtbl);
  gtbl_j1j2 = ERI_Gtbl_j1j2(solver->gtbl);
  gtbl_gc   = ERI_Gtbl_gc(solver->gtbl);

  STEPTRACE("ERI_Alpha: in ");

  switch (solver->rcsh) {
  case ERI_SH_REAL:    ndalp = ngrid*jmax; break;
  case ERI_SH_COMPLEX: ndalp = ngrid*jmax*2; break; 
  }

#if SUMSKIP
  amax  = &alp[ndalp];
  gmax  = &g[ngrid*jmax];
  for (j1=0; j1<jmax; j1++) { amax[j1] = 0.0; }
#endif

  j = l*(l+1)+m;

  /* zero reset */
  for (ia=0; ia<ndalp; ia++) { alp[ia] = 0.0; }
 
  igtbl = ERI_Gtbl_index(solver->gtbl, j);
  n = gtbl_n[j];

  switch (solver->rcsh) {
  case ERI_SH_COMPLEX:
    for (i=0; i<n; i++) {
      fpg = 4.0*PI*gtbl_gc[igtbl];
      j1 = ERI_GTBL_UNPACK_J1( gtbl_j1j2[igtbl] );
      j2 = ERI_GTBL_UNPACK_J2( gtbl_j1j2[igtbl] );
      igtbl++;
    
      l1 = j2l[j1];
      m1 = j2m[j1];
      l2 = j2l[j2];
      m2 = j2m[j2];
#if SUMSKIP
      if (gmax[l1*lmax+l2]>SUMSKIP_THRESHOLD) 
#endif
      {
        ERI_Spherical_Harmonics(l2, m2, theta, phi, sh, dsht, dshp);
        fpgsh[0] = fpg*sh[0];
        fpgsh[1] = fpg*sh[1];
        phase(l-l1+l2, fpgsh);
      
        ia_base = j1*ngrid;
        ig_base = (l1*lmax+l2)*ngrid;
        for (ir=0; ir<ngrid; ir++) {
          ia = 2*(ia_base+ir);
          ig = ig_base+ir;
          alp[ia+0] += fpgsh[0]*g[ig];
          alp[ia+1] += fpgsh[1]*g[ig];
        } 
      }
    } /* loop of i */
#if SUMSKIP
    for (j1=0; j1<jmax; j1++) {
      ia_base = j1*ngrid;
      for (ir=0; ir<ngrid; ir++) {
        ia = 2*(ia_base+ir);
        if (fabs(alp[ia+0])>amax[j1]) { amax[j1] = fabs(alp[ia+0]); }
        if (fabs(alp[ia+1])>amax[j1]) { amax[j1] = fabs(alp[ia+1]); }
      } 
    }
#endif
    break;

  case ERI_SH_REAL:
    for (i=0; i<n; i++) {
      fpg = 4.0*PI*gtbl_gc[igtbl];
      j1 = ERI_GTBL_UNPACK_J1( gtbl_j1j2[igtbl] );
      j2 = ERI_GTBL_UNPACK_J2( gtbl_j1j2[igtbl] );
      igtbl++;
    
      l1 = j2l[j1];
      m1 = j2m[j1];
      l2 = j2l[j2];
      m2 = j2m[j2];
#if SUMSKIP
      if (gmax[l1*lmax+l2]>SUMSKIP_THRESHOLD) 
#endif
      {
        ERI_Real_Spherical_Harmonics(l2, m2, theta, phi, sh, dsht, dshp);
        fpgsh[0] = fpg*sh[0];
        fpgsh[1] = 0.0;
        phase(l-l1+l2, fpgsh);
      
        ia_base = j1*ngrid;
        ig_base = (l1*lmax+l2)*ngrid;
        for (ir=0; ir<ngrid; ir++) {
          ia = ia_base+ir;
          ig = ig_base+ir;
          alp[ia] += fpgsh[0]*g[ig];
        } 
      }
    } /* loop of i */
#if SUMSKIP
    for (j1=0; j1<jmax; j1++) {
      ia_base = j1*ngrid;
      for (ir=0; ir<ngrid; ir++) {
        ia = ia_base+ir;
        if (fabs(alp[ia])>amax[j1]) { amax[j1] = fabs(alp[ia]); }
      } 
    }
#endif
    break;
  }

  STEPTRACE("ERI_Alpha: out");
}



/*----------------------------------------------------------------------
  ERI_Alpha_d

  This calculates the Alpha term and its derivatives.
  If you don't need the derivatives, you may call ERI_Alpha which is 
  faster than this routine.
----------------------------------------------------------------------*/
void ERI_LL_Alpha_d(
  ERI_t        *solver,
  double       *alp,     /* (OUT) Alpha term */
  double       *dalp[3], /* (OUT) Derivatives of Alpha term */
  const double *g,       /* (IN)  Gamma term */
  const double *dg,      /* (IN)  Derivatives of Gamma term */
  double        rd,      /* (IN)  Displacement from the expansion center */
  double        theta,   /*       in the spherical coordinates. */
  double        phi,
  int           l,       /* (IN)  Angular mementum of the orbital */
  int           m
)
{
  int ngrid, lmax, jmax;
  int i, j, ir, j1, j2, l1, m1, l2, m2;
  int ia, ia_base, ig, ig_base;

  double sh[2], dsht[2], dshp[2], fpg, fpgsh[2], fpgdsht[2], fpgdshp[2];
  double dgsh[2], gdsht[2], gdshp[2];
 
  double ct = cos(theta);
  double st = sin(theta);
  double cp = cos(phi);
  double sp = sin(phi);
  double x1, y1, z1, x2, y2, z2, x3, y3;
  
  const int *j2l   = g_itbl_j2l;
  const int *j2m   = g_itbl_j2m;

  int ndalp;
 
  int n, igtbl;
  const int* gtbl_n;
  const long* gtbl_j1j2;
  const double* gtbl_gc;

#if SUMSKIP
  double *amax, *damax[3];
  const double *gmax, *dgmax;
#endif 
  
  const double rd_thresh = 1e-10;
 
  ngrid = ERI_ngrid(solver);
  lmax  = ERI_lmax(solver);
  jmax  = lmax*lmax;

  gtbl_n    = ERI_Gtbl_n(solver->gtbl);
  gtbl_j1j2 = ERI_Gtbl_j1j2(solver->gtbl);
  gtbl_gc   = ERI_Gtbl_gc(solver->gtbl);

  switch (solver->rcsh) {
  case ERI_SH_REAL:    ndalp = ngrid*jmax; break;
  case ERI_SH_COMPLEX: ndalp = ngrid*jmax*2; break;
  }

#if SUMSKIP
  amax     = &alp[ndalp];
  damax[0] = &dalp[0][ndalp];
  damax[1] = &dalp[1][ndalp];
  damax[2] = &dalp[2][ndalp];
  gmax     = &g[ngrid*jmax];
  dgmax    = &dg[ngrid*jmax];
  for (j1=0; j1<jmax; j1++) {
    amax[j1] = 0.0;
    damax[0][j1] = 0.0;
    damax[1][j1] = 0.0;
    damax[2][j1] = 0.0;
  }
#endif 
   
  /* zero reset */ 
  for (ia=0; ia<ndalp; ia++) {
    alp[ia] = 0.0;
    dalp[0][ia] = 0.0;
    dalp[1][ia] = 0.0;
    dalp[2][ia] = 0.0;
  }

  if (fabs(rd)<rd_thresh) {
    /* The both orbitals are on a same site,
       The derivatives are not to be calculated 
       but just set to zero. */
    ERI_LL_Alpha(solver, alp, g, theta, phi, l, m);
  } else {
    /* The orbitals are on different sites.
       The derivatives are calculated as well. */

    j = l*(l+1)+m;
  
    x1 = st*cp;     y1 = st*sp;    z1 = ct;
    x2 = ct*cp/rd;  y2 = ct*sp/rd; z2 = -st/rd;
    x3 = -sp/rd/st; y3 = cp/rd/st;

    igtbl = ERI_Gtbl_index(solver->gtbl, j);
    n = gtbl_n[j];

    switch (solver->rcsh) {
    case ERI_SH_COMPLEX:
      for (i=0; i<n; i++) {
        fpg = 4.0*PI*gtbl_gc[igtbl];
        j1  = ERI_GTBL_UNPACK_J1( gtbl_j1j2[igtbl] );
        j2  = ERI_GTBL_UNPACK_J2( gtbl_j1j2[igtbl] );
        igtbl++;

        l1 = j2l[j1];
        m1 = j2m[j1];
        l2 = j2l[j2];
        m2 = j2m[j2];

        ERI_Spherical_Harmonics(l2, m2, theta, phi, sh, dsht, dshp);
        fpgsh[0]   = fpg*sh[0];
        fpgsh[1]   = fpg*sh[1];
        fpgdsht[0] = fpg*dsht[0];
        fpgdsht[1] = fpg*dsht[1];
        fpgdshp[0] = fpg*dshp[0];
        fpgdshp[1] = fpg*dshp[1];
        phase(l-l1+l2, fpgsh);
        phase(l-l1+l2, fpgdsht);
        phase(l-l1+l2, fpgdshp); 

        ia_base  = j1*ngrid;
        ig_base  = (l1*lmax+l2)*ngrid;
#if SUMSKIP
        if (gmax[l1*lmax+l2]>SUMSKIP_THRESHOLD) 
#endif
        {
          for (ir=0; ir<ngrid; ir++) {
            ia = 2*(ia_base+ir);
            ig = ig_base+ir;
            alp[ia+0] += fpgsh[0]*g[ig];
            alp[ia+1] += fpgsh[1]*g[ig];
            /* derivatives */
            gdsht[0] = g[ig]*fpgdsht[0];
            gdsht[1] = g[ig]*fpgdsht[1];
            gdshp[0] = g[ig]*fpgdshp[0];
            gdshp[1] = g[ig]*fpgdshp[1];

            dalp[0][ia+0] += x2*gdsht[0] + x3*gdshp[0];
            dalp[0][ia+1] += x2*gdsht[1] + x3*gdshp[1];
            dalp[1][ia+0] += y2*gdsht[0] + y3*gdshp[0];
            dalp[1][ia+1] += y2*gdsht[1] + y3*gdshp[1];
            dalp[2][ia+0] += z2*gdsht[0];
            dalp[2][ia+1] += z2*gdsht[1];
          } /* --- loop of ir --- */
        } 
#if SUMSKIP
        if (dgmax[l1*lmax+l2]>SUMSKIP_THRESHOLD) 
#endif
        {
          for (ir=0; ir<ngrid; ir++) {
            /* derivatives */
            ia = 2*(ia_base+ir);
            ig = ig_base+ir;
            dgsh[0] = dg[ig]*fpgsh[0];
            dgsh[1] = dg[ig]*fpgsh[1];
            dalp[0][ia+0] += x1*dgsh[0];
            dalp[0][ia+1] += x1*dgsh[1];
            dalp[1][ia+0] += y1*dgsh[0];
            dalp[1][ia+1] += y1*dgsh[1];
            dalp[2][ia+0] += z1*dgsh[0];
            dalp[2][ia+1] += z1*dgsh[1];
          } /* --- loop of ir --- */
        } 
      } /* --- loop of i ---*/
#if SUMSKIP
      for (j1=0; j1<jmax; j1++) {
        ia_base  = j1*ngrid;
        for (ir=0; ir<ngrid; ir++) {
          ia = 2*(ia_base+ir);
          if (fabs(alp[ia+0])>amax[j1]) { amax[j1] = fabs(alp[ia+0]); }
          if (fabs(alp[ia+1])>amax[j1]) { amax[j1] = fabs(alp[ia+1]); }
          if (fabs(dalp[0][ia+0])>damax[0][j1]) {
            damax[0][j1] = fabs(dalp[0][ia+0]);
          }
          if (fabs(dalp[0][ia+1])>damax[0][j1]) {
            damax[0][j1] = fabs(dalp[0][ia+1]);
          }
          if (fabs(dalp[1][ia+0])>damax[1][j1]) {
            damax[1][j1] = fabs(dalp[1][ia+0]);
          }
          if (fabs(dalp[1][ia+1])>damax[1][j1]) {
            damax[1][j1] = fabs(dalp[1][ia+1]);
          }
          if (fabs(dalp[2][ia+0])>damax[2][j1]) {
            damax[2][j1] = fabs(dalp[2][ia+0]);
          }
          if (fabs(dalp[2][ia+1])>damax[2][j1]) {
            damax[2][j1] = fabs(dalp[2][ia+1]);
          }
        } /* loop of ir */
      } /* loop of j1 */
#endif
      break;

    case ERI_SH_REAL:
      for (i=0; i<n; i++) {
        fpg = 4.0*PI*gtbl_gc[igtbl];
        j1  = ERI_GTBL_UNPACK_J1( gtbl_j1j2[igtbl] );
        j2  = ERI_GTBL_UNPACK_J2( gtbl_j1j2[igtbl] );
        igtbl++;

        l1 = j2l[j1];
        m1 = j2m[j1];
        l2 = j2l[j2];
        m2 = j2m[j2];

        ERI_Real_Spherical_Harmonics(l2, m2, theta, phi, sh, dsht, dshp);
        fpgsh[0]   = fpg*sh[0];
        fpgsh[1]   = 0.0;
        fpgdsht[0] = fpg*dsht[0];
        fpgdsht[1] = 0.0;
        fpgdshp[0] = fpg*dshp[0];
        fpgdshp[1] = 0.0;
        phase(l-l1+l2, fpgsh);
        phase(l-l1+l2, fpgdsht);
        phase(l-l1+l2, fpgdshp); 

        ia_base  = j1*ngrid;
        ig_base  = (l1*lmax+l2)*ngrid;
#if SUMSKIP
        if (gmax[l1*lmax+l2]>SUMSKIP_THRESHOLD) 
#endif
        {
          for (ir=0; ir<ngrid; ir++) {
            ia = ia_base+ir;
            ig = ig_base+ir;
            alp[ia] += fpgsh[0]*g[ig];
            /* derivatives */
            gdsht[0] = g[ig]*fpgdsht[0];
            gdshp[0] = g[ig]*fpgdshp[0];

            dalp[0][ia] += x2*gdsht[0] + x3*gdshp[0];
            dalp[1][ia] += y2*gdsht[0] + y3*gdshp[0];
            dalp[2][ia] += z2*gdsht[0];
          } /* --- loop of ir --- */
        } 
#if SUMSKIP
        if (dgmax[l1*lmax+l2]>SUMSKIP_THRESHOLD) 
#endif
        {
          for (ir=0; ir<ngrid; ir++) {
            /* derivatives */
            ia = ia_base+ir;
            ig = ig_base+ir;
            dgsh[0] = dg[ig]*fpgsh[0];
            dalp[0][ia] += x1*dgsh[0];
            dalp[1][ia] += y1*dgsh[0];
            dalp[2][ia] += z1*dgsh[0];
          } /* --- loop of ir --- */
        } 
      } /* --- loop of i ---*/
#if SUMSKIP
      for (j1=0; j1<jmax; j1++) {
        ia_base  = j1*ngrid;
        for (ir=0; ir<ngrid; ir++) {
          ia = ia_base+ir;
          if (fabs(alp[ia])>amax[j1]) { amax[j1] = fabs(alp[ia]); }
          if (fabs(dalp[0][ia])>damax[0][j1]) { damax[0][j1] = fabs(dalp[0][ia]); }
          if (fabs(dalp[1][ia])>damax[1][j1]) { damax[1][j1] = fabs(dalp[1][ia]); }
          if (fabs(dalp[2][ia])>damax[2][j1]) { damax[2][j1] = fabs(dalp[2][ia]); }
        } /* loop of ir */
      } /* loop of j1 */
#endif
      break;
    }
  }
}






/*----------------------------------------------------------------------
  ERI_LL_Overlap
  
  Product of two orbital
----------------------------------------------------------------------*/
void ERI_LL_Overlap(
  ERI_t        *solver,
  double       *p,   /* (OUT) P term */
  const double *a1,  /* (IN)  Alpha term of the orbtial #1 */
  const double *a2   /* (IN)  Alpha term of the orbital #2 */
)
{
  int ngrid, lmax, jmax, lmax0, jmax0;
  int j, j1, j2, ir, ia1, ia2, ip;
  int ia1_base, ia2_base, ip_base;
  double gnt, a1a2[2];
  
  const int *j2l   = g_itbl_j2l;
  const int *j2m   = g_itbl_j2m;

  int ndp, ndalp;

  int i, n, igtbl;
  const int* gtbl_n;
  const long* gtbl_j1j2;
  const double* gtbl_gc;

#if SUMSKIP
  const double *amax1, *amax2;
  double *pmax;
#endif
 
  ngrid = ERI_ngrid(solver);
  lmax  = ERI_lmax(solver);
  jmax  = lmax*lmax;
  lmax0 = ERI_lmax_gl(solver);
  jmax0 = lmax0*lmax0;

  gtbl_n    = ERI_Gtbl_n(solver->gtbl);
  gtbl_j1j2 = ERI_Gtbl_j1j2(solver->gtbl);
  gtbl_gc   = ERI_Gtbl_gc(solver->gtbl);
     
  STEPTRACE("ERI_LL_Overlap: in\n");

  switch (solver->rcsh) {
  case ERI_SH_COMPLEX: 
    ndalp = 2*ngrid*jmax; 
    ndp   = 2*ngrid*jmax0; 
    break;
  case ERI_SH_REAL:
    ndalp = ngrid*jmax; 
    ndp   = ngrid*jmax0; 
    break;
  }

#if SUMSKIP
  amax1 = &a1[ndalp];
  amax2 = &a2[ndalp];
  pmax  = &p[ndp];
#endif

  /* zero rest */
  for (j=0; j<ndp; j++) { p[j] = 0.0; }
#if SUMSKIP
  for (j=0; j<jmax0; j++) { pmax[j] = 0.0; }
#endif

  switch (solver->rcsh) {
  case ERI_SH_COMPLEX: 
    igtbl=0;
    for (j=0; j<jmax0; j++) {
      n = gtbl_n[j];

      ip_base = j*ngrid;
      for (i=0; i<n; i++) {
        gnt = gtbl_gc[igtbl];
        j1 = ERI_GTBL_UNPACK_J1( gtbl_j1j2[igtbl] );
        j2 = ERI_GTBL_UNPACK_J2( gtbl_j1j2[igtbl] );
        igtbl++;

        ia1_base = j1*ngrid;
        ia2_base = j2*ngrid;
#if SUMSKIP
        if (amax1[j1]*amax2[j2]>SUMSKIP_THRESHOLD) 
#endif
        {
          for (ir=0; ir<ngrid; ir++) {  
            ia1 = 2*(ia1_base+ir); 
            ia2 = 2*(ia2_base+ir); 
            ip  = 2*(ip_base+ir);
            a1a2[0] = a1[ia1+0]*a2[ia2+0] - a1[ia1+1]*a2[ia2+1];
            a1a2[1] = a1[ia1+0]*a2[ia2+1] + a1[ia1+1]*a2[ia2+0];
            p[ip+0] += gnt*a1a2[0];
            p[ip+1] += gnt*a1a2[1];
          } /*--- loop of ir ---*/
        }
      } /*--- loop of i ---*/
#if SUMSKIP
      for (ir=0; ir<ngrid; ir++) {  
        ip = 2*(ip_base+ir);
        if (fabs(p[ip+0])>pmax[j]) { pmax[j] = fabs(p[ip+0]); }
        if (fabs(p[ip+1])>pmax[j]) { pmax[j] = fabs(p[ip+1]); }
      }
#endif
    } /*--- loop of j ---*/
    break;

  case ERI_SH_REAL: 
    igtbl=0;
    for (j=0; j<jmax0; j++) {
      n = gtbl_n[j];

      ip_base = j*ngrid;
      for (i=0; i<n; i++) {
        gnt = gtbl_gc[igtbl];
        j1 = ERI_GTBL_UNPACK_J1( gtbl_j1j2[igtbl] );
        j2 = ERI_GTBL_UNPACK_J2( gtbl_j1j2[igtbl] );
        igtbl++;

        ia1_base = j1*ngrid;
        ia2_base = j2*ngrid;
#if SUMSKIP
        if (amax1[j1]*amax2[j2]>SUMSKIP_THRESHOLD) 
#endif
        {
          for (ir=0; ir<ngrid; ir++) {  
            ia1 = ia1_base+ir; 
            ia2 = ia2_base+ir; 
            ip  = ip_base+ir;
            a1a2[0] = a1[ia1]*a2[ia2];
            p[ip] += gnt*a1a2[0];
            /* imag part should be zero */
          } /*--- loop of ir ---*/
        }
      } /*--- loop of i ---*/
#if SUMSKIP
      for (ir=0; ir<ngrid; ir++) {  
        ip = ip_base+ir;
        if (fabs(p[ip])>pmax[j]) { pmax[j] = fabs(p[ip]); }
      }
#endif
    } /*--- loop of j ---*/
    break;
  }

  STEPTRACE("ERI_Overlap: out\n");
}



/*--------------------------------------------------------------------*/
void ERI_LL_Overlap_d(
  ERI_t        *solver,
  double       *p,      /* (OUT) overlap function */ 
  double       *dp[3],  /* (OUT) derivatives */
  const double *a1,     /* (IN) alpha for orb 1 */
  const double *da1[3], /* (IN) derivatives */
  const double *a2,     /* (IN) alpha for orb 2 */
  const double *da2[3], /* (IN) derivatives */
  double        x       /* (IN) */
)
{
  int ngrid, lmax, jmax, lmax0, jmax0;
  int i, j, j1, j2, ir,  ia1, ia2, ip, ixyz, l, m;
  int ia1_base, ia2_base, ip_base;
  double gnt, a1a2[2], da1a2[3][2], a1da2[3][2];
  
  const int *j2l = g_itbl_j2l;
  const int *j2m = g_itbl_j2m;
  
  int ndp, ndalp;

  int n, igtbl;
  const int* gtbl_n;
  const long* gtbl_j1j2;
  const double* gtbl_gc;

#if SUMSKIP
  double *pmax, *dpmax[3];
  const double *amax1, *damax1[3], *amax2, *damax2[3];
#endif 
  
  ngrid = ERI_ngrid(solver);
  lmax  = ERI_lmax(solver);
  jmax  = lmax*lmax;
  lmax0 = ERI_lmax_gl(solver);
  jmax0 = lmax0*lmax0;

  gtbl_n    = ERI_Gtbl_n(solver->gtbl);
  gtbl_j1j2 = ERI_Gtbl_j1j2(solver->gtbl);
  gtbl_gc   = ERI_Gtbl_gc(solver->gtbl);
  
  switch (solver->rcsh) {
  case ERI_SH_COMPLEX: 
    ndalp = 2*ngrid*jmax;
    ndp   = 2*ngrid*jmax0; 
    break;
  case ERI_SH_REAL:
    ndalp = ngrid*jmax;
    ndp   = ngrid*jmax0;
    break;
  }

#if SUMSKIP
  pmax      = &p[ndp];
  amax1     = &a1[ndalp];
  amax2     = &a2[ndalp];
  for (i=0; i<3; i++) {
    dpmax[i]  = &(dp[i][ndp]);
    damax1[i] = &(da1[i][ndalp]);
    damax2[i] = &(da2[i][ndalp]);
  }
#endif
  STEPTRACE("MCP_P_d: begin\n");

  /* zero reset */
  for (j=0; j<ndp; j++) {
    p[j] = 0.0;
    dp[0][j] = 0.0;
    dp[1][j] = 0.0;
    dp[2][j] = 0.0;
  }
#if SUMSKIP
  for (j=0; j<jmax0; j++) {
    pmax[j]  = 0.0;
    dpmax[0][j]  = 0.0;
    dpmax[1][j]  = 0.0;
    dpmax[2][j]  = 0.0;
  }
#endif 
  STEPTRACE("MCP_P_d: step 1\n");

  switch (solver->rcsh) {
  case ERI_SH_COMPLEX:
    igtbl=0;
    for (j=0; j<jmax0; j++) {
      l = j2l[j];
      m = j2m[j];
      ip_base = j*ngrid;
      n = gtbl_n[j];
      for (i=0; i<n; i++) {
        gnt = gtbl_gc[igtbl];
        j1  = ERI_GTBL_UNPACK_J1( gtbl_j1j2[igtbl] );
        j2  = ERI_GTBL_UNPACK_J2( gtbl_j1j2[igtbl] );
        igtbl++;

        {
          ia1_base = j1*ngrid;
          ia2_base = j2*ngrid;
#if SUMSKIP
          if (amax1[j1]*amax2[j2]>SUMSKIP_THRESHOLD) {
            for (ir=0; ir<ngrid; ir++) {  
              ia1 = 2*(ia1_base+ir); 
              ia2 = 2*(ia2_base+ir); 
              ip  = 2*(ip_base+ir);
              a1a2[0] = a1[ia1+0]*a2[ia2+0] - a1[ia1+1]*a2[ia2+1];
              a1a2[1] = a1[ia1+0]*a2[ia2+1] + a1[ia1+1]*a2[ia2+0];
              p[ip+0] += gnt*a1a2[0];
              p[ip+1] += gnt*a1a2[1];
            } /* loop of ir */
          }
          for (ixyz=0; ixyz<3; ixyz++) {
            if (amax2[j2]*damax1[ixyz][j1]>SUMSKIP_THRESHOLD) {
              for (ir=0; ir<ngrid; ir++) {  
                ia1 = 2*(ia1_base+ir); 
                ia2 = 2*(ia2_base+ir); 
                ip  = 2*(ip_base+ir);
                da1a2[ixyz][0] = da1[ixyz][ia1+0]*a2[ia2+0] 
                                 -da1[ixyz][ia1+1]*a2[ia2+1];
                da1a2[ixyz][1] = da1[ixyz][ia1+0]*a2[ia2+1] 
                                 +da1[ixyz][ia1+1]*a2[ia2+0];
                dp[ixyz][ip+0] += gnt*(1.0-x)*da1a2[ixyz][0];
                dp[ixyz][ip+1] += gnt*(1.0-x)*da1a2[ixyz][1];
              } /* loop of ir */
            }
          } /* loop of ixyz */

          for (ixyz=0; ixyz<3; ixyz++) {
            if (amax1[j1]*damax2[ixyz][j2]>SUMSKIP_THRESHOLD) {
              for (ir=0; ir<ngrid; ir++) {  
                ia1 = 2*(ia1_base+ir); 
                ia2 = 2*(ia2_base+ir); 
                ip  = 2*(ip_base+ir);
                a1da2[ixyz][0] = a1[ia1+0]*da2[ixyz][ia2+0] 
                                 - a1[ia1+1]*da2[ixyz][ia2+1];
                a1da2[ixyz][1] = a1[ia1+0]*da2[ixyz][ia2+1] 
                                 + a1[ia1+1]*da2[ixyz][ia2+0];
                dp[ixyz][ip+0] += -gnt*x*a1da2[ixyz][0];
                dp[ixyz][ip+1] += -gnt*x*a1da2[ixyz][1];
              } /* loop of ir */
            }
          } /* loop of ixyz */
#else /* SUMSKIP */
          for (ir=0; ir<ngrid; ir++) {  
            ia1 = 2*(ia1_base+ir); 
            ia2 = 2*(ia2_base+ir); 
            ip  = 2*(ip_base+ir);
              
            a1a2[0] = a1[ia1+0]*a2[ia2+0]-a1[ia1+1]*a2[ia2+1];
            a1a2[1] = a1[ia1+0]*a2[ia2+1]+a1[ia1+1]*a2[ia2+0];
            p[ip+0] += gnt*a1a2[0];
            p[ip+1] += gnt*a1a2[1];
              
            da1a2[0][0] = da1[0][ia1+0]*a2[ia2+0]-da1[0][ia1+1]*a2[ia2+1];
            da1a2[0][1] = da1[0][ia1+0]*a2[ia2+1]+da1[0][ia1+1]*a2[ia2+0];
            da1a2[1][0] = da1[1][ia1+0]*a2[ia2+0]-da1[1][ia1+1]*a2[ia2+1];
            da1a2[1][1] = da1[1][ia1+0]*a2[ia2+1]+da1[1][ia1+1]*a2[ia2+0];
            da1a2[2][0] = da1[2][ia1+0]*a2[ia2+0]-da1[2][ia1+1]*a2[ia2+1];
            da1a2[2][1] = da1[2][ia1+0]*a2[ia2+1]+da1[2][ia1+1]*a2[ia2+0];
              
            a1da2[0][0] = a1[ia1+0]*da2[0][ia2+0]-a1[ia1+1]*da2[0][ia2+1];
            a1da2[0][1] = a1[ia1+0]*da2[0][ia2+1]+a1[ia1+1]*da2[0][ia2+0];
            a1da2[1][0] = a1[ia1+0]*da2[1][ia2+0]-a1[ia1+1]*da2[1][ia2+1];
            a1da2[1][1] = a1[ia1+0]*da2[1][ia2+1]+a1[ia1+1]*da2[1][ia2+0];
            a1da2[2][0] = a1[ia1+0]*da2[2][ia2+0]-a1[ia1+1]*da2[2][ia2+1];
            a1da2[2][1] = a1[ia1+0]*da2[2][ia2+1]+a1[ia1+1]*da2[2][ia2+0];
 
            dp[0][ip+0] += gnt*((1.0-x)*da1a2[0][0] - x*a1da2[0][0]);
            dp[0][ip+1] += gnt*((1.0-x)*da1a2[0][1] - x*a1da2[0][1]);
            dp[1][ip+0] += gnt*((1.0-x)*da1a2[1][0] - x*a1da2[1][0]);
            dp[1][ip+1] += gnt*((1.0-x)*da1a2[1][1] - x*a1da2[1][1]);
            dp[2][ip+0] += gnt*((1.0-x)*da1a2[2][0] - x*a1da2[2][0]);
            dp[2][ip+1] += gnt*((1.0-x)*da1a2[2][1] - x*a1da2[2][1]);
          } /*--- loop of ir ---*/
#endif /* SUMSKIP */
        }
      } /*--- loop of i ---*/
#if SUMSKIP
      for (ir=0; ir<ngrid; ir++) {  
        ip = 2*(ip_base+ir);
        if (fabs(p[ip+0])>pmax[j]) { pmax[j] = fabs(p[ip+0]); }
        if (fabs(p[ip+1])>pmax[j]) { pmax[j] = fabs(p[ip+1]); }
        for (i=0; i<3; i++) {
          if (fabs(dp[i][ip+0])>dpmax[i][j]) {
            dpmax[i][j] = fabs(dp[i][ip+0]);
          }
          if (fabs(dp[i][ip+1])>dpmax[i][j]) {
            dpmax[i][j] = fabs(dp[i][ip+1]);
          }
        } /* loop of i */
      } /* loop of ir */
#endif
    } /*--- loop of j ---*/
    break;

  case ERI_SH_REAL:
    igtbl=0;
    for (j=0; j<jmax0; j++) {
      l = j2l[j];
      m = j2m[j];
      ip_base = j*ngrid;
      n = gtbl_n[j];
      for (i=0; i<n; i++) {
        gnt = gtbl_gc[igtbl];
        j1  = ERI_GTBL_UNPACK_J1( gtbl_j1j2[igtbl] );
        j2  = ERI_GTBL_UNPACK_J2( gtbl_j1j2[igtbl] );
        igtbl++;

        {
          ia1_base = j1*ngrid;
          ia2_base = j2*ngrid;
#if SUMSKIP
          if (amax1[j1]*amax2[j2]>SUMSKIP_THRESHOLD) {
            for (ir=0; ir<ngrid; ir++) {  
              ia1 = ia1_base+ir; 
              ia2 = ia2_base+ir; 
              ip  = ip_base+ir;
              a1a2[0] = a1[ia1]*a2[ia2];
              p[ip] += gnt*a1a2[0];
              /* imag part should be zero */
            } /* loop of ir */
          }
          for (ixyz=0; ixyz<3; ixyz++) {
            if (amax2[j2]*damax1[ixyz][j1]>SUMSKIP_THRESHOLD) {
              for (ir=0; ir<ngrid; ir++) {  
                ia1 = ia1_base+ir; 
                ia2 = ia2_base+ir; 
                ip  = ip_base+ir;
                da1a2[ixyz][0] = da1[ixyz][ia1]*a2[ia2];
                dp[ixyz][ip] += gnt*(1.0-x)*da1a2[ixyz][0];
                /* imag part should be zero */
              } /* loop of ir */
            }
          } /* loop of ixyz */

          for (ixyz=0; ixyz<3; ixyz++) {
            if (amax1[j1]*damax2[ixyz][j2]>SUMSKIP_THRESHOLD) {
              for (ir=0; ir<ngrid; ir++) {  
                ia1 = ia1_base+ir; 
                ia2 = ia2_base+ir; 
                ip  = ip_base+ir;
                a1da2[ixyz][0] = a1[ia1]*da2[ixyz][ia2];
                dp[ixyz][ip] += -gnt*x*a1da2[ixyz][0];
                /* imag part should be zero */
              } /* loop of ir */
            }
          } /* loop of ixyz */
#else /* SUMSKIP */
          for (ir=0; ir<ngrid; ir++) {  
            ia1 = ia1_base+ir; 
            ia2 = ia2_base+ir; 
            ip  = ip_base+ir;
              
            a1a2[0] = a1[ia1]*a2[ia2];
            p[ip] += gnt*a1a2[0];
            /* imag part should be zero */
              
            da1a2[0][0] = da1[0][ia1]*a2[ia2];
            da1a2[1][0] = da1[1][ia1]*a2[ia2];
            da1a2[2][0] = da1[2][ia1]*a2[ia2];
              
            a1da2[0][0] = a1[ia1]*da2[0][ia2];
            a1da2[1][0] = a1[ia1]*da2[1][ia2];
            a1da2[2][0] = a1[ia1]*da2[2][ia2];
 
            dp[0][ip] += gnt*((1.0-x)*da1a2[0][0] - x*a1da2[0][0]);
            /* imag part should be zero */
            dp[1][ip] += gnt*((1.0-x)*da1a2[1][0] - x*a1da2[1][0]);
            /* imag part should be zero */
            dp[2][ip] += gnt*((1.0-x)*da1a2[2][0] - x*a1da2[2][0]);
            /* imag part should be zero */
          } /*--- loop of ir ---*/
#endif /* SUMSKIP */
        }
      } /*--- loop of i ---*/
#if SUMSKIP
      for (ir=0; ir<ngrid; ir++) {  
        ip = ip_base+ir;
        if (fabs(p[ip])>pmax[j]) { pmax[j] = fabs(p[ip]); }
        for (i=0; i<3; i++) {
          if (fabs(dp[i][ip])>dpmax[i][j]) { dpmax[i][j] = fabs(dp[i][ip]); }
        } /* loop of i */
      } /* loop of ir */
#endif
    } /*--- loop of j ---*/
  }
}




void ERI_Transform_Overlap(
  ERI_t        *solver,
  double       *F,
  const double *P
)
{
  int ngrid, lmax0, jmax0;
  int j, ip, ip_base, l;
        
#if WORKSPACE_AT_STACK
  double in[ERI_NGRIDMAX*2], out[ERI_NGRIDMAX*2];
#else
  double *in  = solver->ws_in;
  double *out = solver->ws_out;
#endif

  int ndp;

#if SUMSKIP
  int ik;
  double *fmax;
  const double *pmax;
#endif
  
  ngrid = ERI_ngrid(solver);
  lmax0 = ERI_lmax_gl(solver);
  jmax0 = lmax0*lmax0;

  STEPTRACE("ERI_Transform_Overlap: in");
  
  switch (solver->rcsh) {
  case ERI_SH_COMPLEX: ndp = 2*ngrid*jmax0; break;
  case ERI_SH_REAL:    ndp = ngrid*jmax0; break;
  }

#if SUMSKIP
  pmax = &P[ndp];
  fmax = &F[ndp];
  for (j=0; j<jmax0; j++) { fmax[j] = 0.0; }
#endif
 
  switch (solver->rcsh) {
  case ERI_SH_COMPLEX: 
    for (j=0; j<jmax0; j++) {
#if SUMSKIP
      if (pmax[j]>SUMSKIP_THRESHOLD) 
#endif
      {
        l=(int)sqrt((double)j);
        ip_base = j*ngrid;
        ERI_SBT_Transform(solver->sbt, &F[2*ip_base], &P[2*ip_base], 
                          l, ERI_SBT_FORWARD);
#if SUMSKIP
        for (ik=0; ik<ngrid; ik++) {
          ip = 2*(ip_base + ik);
          if (fabs(F[ip+0])>fmax[j]) { fmax[j] = fabs(F[ip+0]); }
          if (fabs(F[ip+1])>fmax[j]) { fmax[j] = fabs(F[ip+1]); }
        }
#endif
      }
    } /* loop of j */
    break;

  case ERI_SH_REAL:  
    for (j=0; j<jmax0; j++) {
#if SUMSKIP
      if (pmax[j]>SUMSKIP_THRESHOLD) 
#endif
      {
        l=(int)sqrt((double)j);
        ip_base = j*ngrid;
        for (ik=0; ik<ngrid; ik++) {
          ip = ip_base+ik;
          in[2*ik+0] = P[ip];
          in[2*ik+1] = 0.0;
        }
        ERI_SBT_Transform(solver->sbt, out, in, l, ERI_SBT_FORWARD);
        for (ik=0; ik<ngrid; ik++) {
          ip = ip_base+ik;
          F[ip] = out[2*ik+0];
          /* imag part out[2*ir+1] is thrown away */
#if SUMSKIP
          if (fabs(F[ip])>fmax[j]) { fmax[j] = fabs(F[ip]); }
#endif
        }
      }
    } /* loop of j */
    break;
  }
  
  STEPTRACE("ERI_Transofrm_Overlap: out");
}




/*----------------------------------------------------------------------
  ERI_GL_Interpolate
----------------------------------------------------------------------*/
void ERI_GL_Interpolate(
  ERI_t        *solver,
  double       *glF, /* (IN) Overlap matrix */
  const double *F    /* (IN) Overlap matrix */
)
{
  int ngrid, ngl, lmax0, jmax0;
  int j, ik, ig;
  double k, y[2];
  
  const double *kmesh, *glx;

#if SUMSKIP
  int ndf, ndglf;
  double *glfmax;
  const double *fmax;
#endif

  ERI_CSpline_t *cs;
 
  ngrid = ERI_ngrid(solver);
  ngl   = ERI_ngl(solver);
  lmax0 = ERI_lmax_gl(solver);
  jmax0 = lmax0*lmax0;

  kmesh = ERI_Mesh_Array_k(solver);
  glx   = ERI_Mesh_Array_glx(solver);

#if SUMSKIP
  switch (solver->rcsh) {
  case ERI_SH_COMPLEX: 
    ndglf = ngl*jmax0*2;
    ndf   = ngrid*jmax0*2;
    break;
  case ERI_SH_REAL:
    ndglf = ngl*jmax0;
    ndf   = ngrid*jmax0;
    break;
  }
  glfmax = &glF[ndglf];
  fmax   = &F[ndf];
  for (j=0; j<jmax0; j++) { glfmax[j] = 0.0; }
#endif

  switch (solver->rcsh) {
  case ERI_SH_COMPLEX:

    for (j=0; j<jmax0; j++) {
#if SUMSKIP
      if (fmax[j]>SUMSKIP_THRESHOLD) 
#endif
      {
        cs = ERI_CSpline_Init_Complex(kmesh, 
          (const double*)&F[j*ngrid*2], ngrid);
        for (ik=0; ik<ngl; ik++) {
          k = glx[ik];
          ig = 2*(j*ngl+ik);
          ERI_CSpline_Eval_Complex(y, k, cs);
          glF[ig+0] = y[0];
          glF[ig+1] = y[1];
#if SUMSKIP
          if (fabs(glF[ig+0])>glfmax[j]) { glfmax[j]=fabs(glF[ig+0]); }
          if (fabs(glF[ig+1])>glfmax[j]) { glfmax[j]=fabs(glF[ig+1]); }
#endif
        } /* loop of ik */
        ERI_CSpline_Free(cs);
      }
    } /* loop of j */
    break;

  case ERI_SH_REAL:
    for (j=0; j<jmax0; j++) {
#if SUMSKIP
      if (fmax[j]>SUMSKIP_THRESHOLD) 
#endif
      {
        cs = ERI_CSpline_Init(kmesh, (const double*)&F[j*ngrid], ngrid);
        for (ik=0; ik<ngl; ik++) {
          k = glx[ik];
          ig = j*ngl+ik;
          glF[ig] = ERI_CSpline_Eval(k, cs);
#if SUMSKIP
          if (fabs(glF[ig])>glfmax[j]) { glfmax[j]=fabs(glF[ig]); }
#endif
        } /* loop of ik */
        ERI_CSpline_Free(cs);
      }
    } /* loop of j */
    break;
  }

}



/*----------------------------------------------------------------------
  ERI_Integral_GL

  If you need the derivatives, you should call ERI_I4_d.
----------------------------------------------------------------------*/
void ERI_Integral_GL(
  ERI_t        *solver,
  double        I4[2], /* (OUT) */
  const double *F1,    /* (IN) Overlap matrix */
  const double *F2,    /* (IN) */
  double        R,     /* (IN) Displacement of two expansion centers */
  double        theta,      
  double        phi,
  double        omega,  /* (IN) screening parameter */
  int           lmax1
)
{
  int ngl, lmax, lmax0, jmax0, jmax1;
  int i, j, j1, j2, ik, l, m, l1, l2;
  int gl_ibase1, gl_ibase2, igl1, igl2;
  
  double sh[2], dsht[2], dshp[2];
  double k, gnt, gnt_r[2], sum[2], FG[2], int1[2], sum2[2];

#if WORKSPACE_AT_STACK
  double tmp_sb[(ERI_LMAXMAX*2)*ERI_NGRIDMAX];
  double tmp_dsb[(ERI_LMAXMAX*2)*ERI_NGRIDMAX];
#else
  double *tmp_sb = solver->ws_sb;
  double *tmp_dsb = solver->ws_dsb;
#endif
  
  const double *glx, *glw;
  double *glj;
  
  const int *j2l = g_itbl_j2l;
  const int *j2m = g_itbl_j2m;
  
  double scr, scr_a;

  int n, igtbl;
  const int* gtbl_n;
  const long* gtbl_j1j2;
  const double* gtbl_gc;

  int lphase;

#if SUMSKIP
  int ndglf;
  const double *fmax1, *fmax2;
#endif

  ngl   = ERI_ngl(solver);
  lmax  = ERI_lmax(solver);
  lmax0 = ERI_lmax_gl(solver);
  jmax0 = lmax0*lmax0;
  jmax1 = lmax1*lmax1;

  glx  = solver->glq_x;
  glw  = solver->glq_w;
  glj  = solver->glq_j;

  gtbl_n    = ERI_Gtbl_n(solver->gtbl);
  gtbl_j1j2 = ERI_Gtbl_j1j2(solver->gtbl);
  gtbl_gc   = ERI_Gtbl_gc(solver->gtbl);

#if SUMSKIP
  switch (solver->rcsh) {
  case ERI_SH_COMPLEX: ndglf = 2*ngl*jmax0; break; 
  case ERI_SH_REAL:    ndglf = ngl*jmax0; break;
  } 
  fmax1 = &F1[ndglf];
  fmax2 = &F2[ndglf];
#endif

  /* spherical Bessel function */
  if (omega<0.0) {
    /* no screening */
    for (ik=0; ik<ngl; ik++) {
      k = glx[ik];
      ERI_Spherical_Bessel(k*R, lmax, tmp_sb, tmp_dsb);
      for (l=0; l<lmax; l++) {
        i = l*ngl + ik;
        glj[i] = glw[ik]*exp(k)*tmp_sb[l]*32.0*PI;
      }
    }
  } else {
    /* screening */
    scr_a = 0.25/omega/omega;
    for (ik=0; ik<ngl; ik++) {
      k = glx[ik];
      ERI_Spherical_Bessel(k*R, lmax, tmp_sb, tmp_dsb);
      scr = 1.0-exp(-k*k*scr_a);
      for (l=0; l<lmax; l++) {
        i = l*ngl + ik;
        glj[i] = glw[ik]*exp(k)*tmp_sb[l]*32.0*PI*scr;
      }
    }
  }

  I4[0] = 0.0;
  I4[1] = 0.0;

  switch (solver->rcsh) {
  case ERI_SH_COMPLEX:
    igtbl=0; 
    for (j=0; j<jmax1; j++) {
      l = j2l[j];
      m = j2m[j];
      ERI_Spherical_Harmonics(l, m, theta, phi, sh, dsht, dshp);
      int1[0] = 0.0;
      int1[1] = 0.0;
      n = gtbl_n[j];
      for (i=0; i<n; i++) {
        gnt = gtbl_gc[igtbl];
        j1 = ERI_GTBL_UNPACK_J1( gtbl_j1j2[igtbl] ); 
        j2 = ERI_GTBL_UNPACK_J2( gtbl_j1j2[igtbl] ); 
        igtbl++;

        if (j1<jmax1 && j2<jmax1) {
          sum2[0] = 0.0;
          sum2[1] = 0.0;
#if SUMSKIP
          if (fmax1[j1]*fmax2[j2]>SUMSKIP_THRESHOLD) 
#endif
          {
            l1 = j2l[j1];
            l2 = j2l[j2];
            sum[0] = 0.0;
            sum[1] = 0.0;
            gl_ibase1 = j1*ngl;
            gl_ibase2 = j2*ngl;
            for (ik=0; ik<ngl; ik++) {
              igl1 = 2*(gl_ibase1+ik);
              igl2 = 2*(gl_ibase2+ik);
              FG[0] = F1[igl1+0]*F2[igl2+0] - F1[igl1+1]*F2[igl2+1];
              FG[1] = F1[igl1+0]*F2[igl2+1] + F1[igl1+1]*F2[igl2+0];
              sum[0] += glj[l*ngl+ik]*FG[0];
              sum[1] += glj[l*ngl+ik]*FG[1];
            } 
            phase(l2+l-l1, sum);
            sum2[0] += sum[0];
            sum2[1] += sum[1];
          }
          int1[0] += sum2[0]*gnt; 
          int1[1] += sum2[1]*gnt; 
        } 
      } 
      I4[0] += int1[0]*sh[0] - int1[1]*sh[1];
      I4[1] += int1[1]*sh[0] + int1[0]*sh[1];
    } /* loop of j */
    break;
  
  case ERI_SH_REAL:
    igtbl=0; 
    for (j=0; j<jmax1; j++) {
      l = j2l[j];
      m = j2m[j];
      ERI_Real_Spherical_Harmonics(l, m, theta, phi, sh, dsht, dshp);
      sh[1] = 0.0;

      int1[0] = 0.0;
      int1[1] = 0.0;
      n = gtbl_n[j];
      for (i=0; i<n; i++) {
        gnt = gtbl_gc[igtbl];
        j1 = ERI_GTBL_UNPACK_J1( gtbl_j1j2[igtbl] ); 
        j2 = ERI_GTBL_UNPACK_J2( gtbl_j1j2[igtbl] ); 
        igtbl++;
        
        if (j1>=jmax1 || j2>=jmax1) { continue; }
#if SUMSKIP
        if (fmax1[j1]*fmax2[j2]<SUMSKIP_THRESHOLD) { continue; }
#endif
        l1 = j2l[j1];
        l2 = j2l[j2];
        lphase = l2+l-l1;
        if (1==lphase%2) { continue; }

        sum[0] = 0.0;
        /* sum[1] = 0.0; */

        gl_ibase1 = j1*ngl;
        gl_ibase2 = j2*ngl;
        for (ik=0; ik<ngl; ik++) {
          igl1 = gl_ibase1+ik;
          igl2 = gl_ibase2+ik;
          FG[0] = F1[igl1]*F2[igl2];
          sum[0] += glj[l*ngl+ik]*FG[0];
        }

        if (2==lphase%4) { sum[0] = -sum[0]; }
        int1[0] += sum[0]*gnt;
      } 
      I4[0] += int1[0]*sh[0];
      /* I4[1] += int1[1]*sh[0]; */
    } /* loop of j */
    break;
  } 
}



/*----------------------------------------------------------------------
  ERI_Integral_GL_d
----------------------------------------------------------------------*/
void ERI_Integral_GL_d(
  ERI_t        *solver,
  double        I4[2],        /* (OUT) integrals */
  double        dI4[4][3][2], /* (OUT) derivatives */
  const double *F1,           /* (IN) overlap matrix 1 */
  const double *F2,           /* (IN) overlap matrix 2 */
  const double *dF1[3],       /* (IN) derivatives of overlap 1 */
  const double *dF2[3],       /* (IN) derivatives of overlap 2 */
  double        R,            /* (IN) displacement */
  double        theta,      
  double        phi,
  double        cx12,         /* (IN) */
  double        cx34,
  double        delta,
  double        omega,  /* (IN) screening parameter */
  int           lmax1
)
{
  int ngl, lmax, lmax0, jmax0, jmax1;
  int i, j, j1, j2, ik, ixyz, l, m, l1, l2;
  int gl_ibase1, gl_ibase2, igl1, igl2;

  double sh[2], dsht[2], dshp[2];
  double k, gnt;
  double sum[2], dsum[2];
  double FG[2], dFG[3][2], FdG[3][2];
  double int1[2], int2[2], int3[3][2], int4[3][2];
  double q1[2], q2[2], q3[2], q4[2], q5[3][2], q6[2], q7[2];
  double ct, st, cp, sp;
  double x1, y1, z1, x2, y2, z2, x3, y3, Rd;
  
#if WORKSPACE_AT_STACK
  double tmp_sb[(ERI_LMAXMAX*2)*ERI_NGRIDMAX];
  double tmp_dsb[(ERI_LMAXMAX*2)*ERI_NGRIDMAX];
#else
  double *tmp_sb  = solver->ws_sb;
  double *tmp_dsb = solver->ws_dsb;
#endif
  
  const double *glx, *glw;
  double *glj, *gldj;
  
  const int *j2l  = g_itbl_j2l;
  const int *j2m  = g_itbl_j2m;
  
  double scr, scr_a;

  int n, igtbl;
  const int* gtbl_n;
  const long* gtbl_j1j2;
  const double* gtbl_gc;
#if 0
  int m1, m2;
#endif

#if SUMSKIP
  int ndglf;
  const double *fmax1, *fmax2;
  const double *dfmax1[3], *dfmax2[3];
#endif

  ngl   = ERI_ngl(solver);
  lmax  = ERI_lmax(solver);
  lmax0 = ERI_lmax_gl(solver);
  jmax0 = lmax0*lmax0;
  jmax1 = lmax1*lmax1;

  glx  = solver->glq_x;
  glw  = solver->glq_w;
  glj  = solver->glq_j;
  gldj = solver->glq_dj;

  gtbl_n    = ERI_Gtbl_n(solver->gtbl);
  gtbl_j1j2 = ERI_Gtbl_j1j2(solver->gtbl);
  gtbl_gc   = ERI_Gtbl_gc(solver->gtbl);

#if SUMSKIP
  switch (solver->rcsh) {
  case ERI_SH_COMPLEX: ndglf = 2*ngl*jmax0; break;
  case ERI_SH_REAL:    ndglf = ngl*jmax0; break;
  }

  fmax1 = &F1[ndglf];
  fmax2 = &F2[ndglf];
  dfmax1[0] = &dF1[0][ndglf];
  dfmax1[1] = &dF1[1][ndglf];
  dfmax1[2] = &dF1[2][ndglf];
  dfmax2[0] = &dF2[0][ndglf];
  dfmax2[1] = &dF2[1][ndglf];
  dfmax2[2] = &dF2[2][ndglf];
#endif

  Rd = sqrt(R*R+delta*delta*exp(-R*R/delta/delta));
  ct = cos(theta);
  st = sin(theta);
  cp = cos(phi);
  sp = sin(phi);

  x1 = st*cp;     y1 = st*sp;    z1 = ct;
  x2 = ct*cp/Rd;  y2 = ct*sp/Rd; z2 = -st/Rd;
  x3 = -sp/Rd/st; y3 = cp/Rd/st;

  if (fabs(theta)<1e-10) { x3 = 0.0; y3 = 0.0; }
  
  /* spherical Bessel function */
  if (omega<0.0) {
    /* no screening */
    for (ik=0; ik<ngl; ik++) {
      k = glx[ik];
      ERI_Spherical_Bessel(k*R, lmax, tmp_sb, tmp_dsb);
      for (l=0; l<lmax; l++) {
        i = l*ngl + ik;
        glj[i]  = glw[ik]*exp(k)*tmp_sb[l]*32.0*PI;
        gldj[i] = glw[ik]*exp(k)*tmp_dsb[l]*k*32.0*PI;
      }
    }
  } else {
    /* screening */
    scr_a = 0.25/omega/omega;
    for (ik=0; ik<ngl; ik++) {
      k = glx[ik];
      ERI_Spherical_Bessel(k*R, lmax, tmp_sb, tmp_dsb);
      scr = 1.0-exp(-k*k*scr_a);
      for (l=0; l<lmax; l++) {
        i = l*ngl + ik;
        glj[i]  = glw[ik]*exp(k)*tmp_sb[l]*32.0*PI*scr;
        gldj[i] = glw[ik]*exp(k)*tmp_dsb[l]*k*32.0*PI*scr;
      }
    }
  }

  I4[0] = 0.0;
  I4[1] = 0.0;
  for (i=0; i<4; i++) {
    for (ixyz=0; ixyz<3; ixyz++) {
      dI4[i][ixyz][0] = 0.0;
      dI4[i][ixyz][1] = 0.0;
    }
  }

  switch (solver->rcsh) {
  case ERI_SH_COMPLEX:
    igtbl=0; 
    for (j=0; j<jmax1; j++) {
      l = j2l[j];
      m = j2m[j];
      ERI_Spherical_Harmonics(l, m, theta, phi, sh, dsht, dshp);
      int1[0] = 0.0;
      int1[1] = 0.0;
      int2[0] = 0.0;
      int2[1] = 0.0;
      for (ixyz=0; ixyz<3; ixyz++) {
        int3[ixyz][0] = 0.0;
        int3[ixyz][1] = 0.0;
        int4[ixyz][0] = 0.0;
        int4[ixyz][1] = 0.0;
      }
      n = gtbl_n[j];
      for (i=0; i<n; i++) {
        gnt = gtbl_gc[igtbl];
        j1  = ERI_GTBL_UNPACK_J1( gtbl_j1j2[igtbl] );
        j2  = ERI_GTBL_UNPACK_J2( gtbl_j1j2[igtbl] );
        igtbl++;

        if (j1<jmax1 && j2<jmax1) {
          l1 = j2l[j1];
          l2 = j2l[j2];
          gl_ibase1 = j1*ngl;
          gl_ibase2 = j2*ngl;
#if SUMSKIP
          if (fmax1[j1]*fmax2[j2]>SUMSKIP_THRESHOLD) 
#endif
          {
            sum[0] = 0.0;
            sum[1] = 0.0;
            dsum[0] = 0.0;
            dsum[1] = 0.0;
            for (ik=0; ik<ngl; ik++) {
              igl1 = 2*(gl_ibase1+ik);
              igl2 = 2*(gl_ibase2+ik);
              FG[0] = F1[igl1+0]*F2[igl2+0] - F1[igl1+1]*F2[igl2+1];
              FG[1] = F1[igl1+0]*F2[igl2+1] + F1[igl1+1]*F2[igl2+0];
              sum[0] += glj[l*ngl+ik]*FG[0];
              sum[1] += glj[l*ngl+ik]*FG[1];
              dsum[0] += gldj[l*ngl+ik]*FG[0];
              dsum[1] += gldj[l*ngl+ik]*FG[1];
            }
            phase(l2+l-l1, sum);
            phase(l2+l-l1, dsum);

            int1[0] += sum[0]*gnt; 
            int1[1] += sum[1]*gnt; 

            int2[0] += dsum[0]*gnt; 
            int2[1] += dsum[1]*gnt; 
          }
          for (ixyz=0; ixyz<3; ixyz++) {
#if SUMSKIP
            if (dfmax1[ixyz][j1]*fmax2[j2]>SUMSKIP_THRESHOLD) 
#endif
            {
              dsum[0] = 0.0;
              dsum[1] = 0.0;
              for (ik=0; ik<ngl; ik++) {
                igl1 = 2*(gl_ibase1+ik);
                igl2 = 2*(gl_ibase2+ik);
                dFG[ixyz][0] = dF1[ixyz][igl1+0]*F2[igl2+0]
                               -dF1[ixyz][igl1+1]*F2[igl2+1];
                dFG[ixyz][1] = dF1[ixyz][igl1+0]*F2[igl2+1]
                               +dF1[ixyz][igl1+1]*F2[igl2+0];
                dsum[0] += glj[l*ngl+ik]*dFG[ixyz][0];
                dsum[1] += glj[l*ngl+ik]*dFG[ixyz][1];
              }
              phase(l2+l-l1, dsum);
              int3[ixyz][0] += dsum[0]*gnt; 
              int3[ixyz][1] += dsum[1]*gnt; 
            }
#if SUMSKIP
            if (fmax1[j1]*dfmax2[ixyz][j2]>SUMSKIP_THRESHOLD) 
#endif
            {
              dsum[0] = 0.0;
              dsum[1] = 0.0;
              for (ik=0; ik<ngl; ik++) {
                igl1 = 2*(gl_ibase1+ik);
                igl2 = 2*(gl_ibase2+ik);
                FdG[ixyz][0] = F1[igl1+0]*dF2[ixyz][igl2+0]
                               -F1[igl1+1]*dF2[ixyz][igl2+1];
                FdG[ixyz][1] = F1[igl1+0]*dF2[ixyz][igl2+1]
                               +F1[igl1+1]*dF2[ixyz][igl2+0];
                dsum[0] += glj[l*ngl+ik]*FdG[ixyz][0];
                dsum[1] += glj[l*ngl+ik]*FdG[ixyz][1];
              } 
              phase(l2+l-l1, dsum);
              int4[ixyz][0] += dsum[0]*gnt; 
              int4[ixyz][1] += dsum[1]*gnt; 
            }
          } /* loop of ixyz */
        } 
      } 
 
      q1[0] = int1[0]*sh[0] - int1[1]*sh[1];
      q1[1] = int1[1]*sh[0] + int1[0]*sh[1];

      I4[0] += q1[0];
      I4[1] += q1[1];

      q2[0] = int2[0]*sh[0] - int2[1]*sh[1];
      q2[1] = int2[1]*sh[0] + int2[0]*sh[1];

      q3[0] = int1[0]*dsht[0] - int1[1]*dsht[1];
      q3[1] = int1[1]*dsht[0] + int1[0]*dsht[1];
      q4[0] = int1[0]*dshp[0] - int1[1]*dshp[1];
      q4[1] = int1[1]*dshp[0] + int1[0]*dshp[1];
   
      q5[0][0] = x1*q2[0] + x2*q3[0] + x3*q4[0]; 
      q5[0][1] = x1*q2[1] + x2*q3[1] + x3*q4[1]; 
      q5[1][0] = y1*q2[0] + y2*q3[0] + y3*q4[0]; 
      q5[1][1] = y1*q2[1] + y2*q3[1] + y3*q4[1]; 
      q5[2][0] = z1*q2[0] + z2*q3[0];
      q5[2][1] = z1*q2[1] + z2*q3[1];

      for (ixyz=0; ixyz<3; ixyz++) {
        q6[0] = int3[ixyz][0]*sh[0] - int3[ixyz][1]*sh[1];
        q6[1] = int3[ixyz][1]*sh[0] + int3[ixyz][0]*sh[1];
        q7[0] = int4[ixyz][0]*sh[0] - int4[ixyz][1]*sh[1];
        q7[1] = int4[ixyz][1]*sh[0] + int4[ixyz][0]*sh[1];
     
        dI4[0][ixyz][0] +=      -cx12*q5[ixyz][0] + q6[0];
        dI4[0][ixyz][1] +=      -cx12*q5[ixyz][1] + q6[1];
        dI4[1][ixyz][0] += (cx12-1.0)*q5[ixyz][0] - q6[0];
        dI4[1][ixyz][1] += (cx12-1.0)*q5[ixyz][1] - q6[1];
        dI4[2][ixyz][0] +=       cx34*q5[ixyz][0] + q7[0];
        dI4[2][ixyz][1] +=       cx34*q5[ixyz][1] + q7[1];
        dI4[3][ixyz][0] += (1.0-cx34)*q5[ixyz][0] - q7[0];
        dI4[3][ixyz][1] += (1.0-cx34)*q5[ixyz][1] - q7[1];
      } 
    } /* loop of j */
    break;

  case ERI_SH_REAL:
    abort();
  }

}




void ERI_Integral_GL_PrejY(
  ERI_t        *solver,
  double       *R,     /* (IN) Displacement of two expansion centers */
  double       *theta,      
  double       *phi,
  int           numR,
  double        omega,  /* (IN) screening parameter */
  double       *prej,   /* [lmax*ngl*numR] */
  double       *preY,   /* [numR*jmax1] */
  int          *mul_j2, /* [jmax1*jmax1*jmax1] */
  double       *mul_gc, /* [jmax1*jmax1*jmax1] */
  int          *mul_n,  /* [jmax1*jmax1] */
  int          *minimalR,     /* [numR] */
  int          *num_minimalR
) 
{
  int ngl, lmax, jmax;
  int i, j, j1, j2, ik, l, m, l1, l2;
  
  double sh[2], dsht[2], dshp[2];
  double k, z;

#if WORKSPACE_AT_STACK
  double tmp_sb[(ERI_LMAXMAX*2)*ERI_NGRIDMAX];
  double tmp_dsb[(ERI_LMAXMAX*2)*ERI_NGRIDMAX];
#else
  double *tmp_sb = solver->ws_sb;
  double *tmp_dsb = solver->ws_dsb;
#endif
  
  const double *glx, *glw;
  
  const int *j2l = g_itbl_j2l;
  const int *j2m = g_itbl_j2m;
  
  double scr;

  int iR, imR, num_mR;
  int nmul;

  int mR2iR[MAXNUMR];
  double R0;

  int n, igtbl, igtbl0;
  const int* gtbl_n;
  const long* gtbl_j1j2;
  const double* gtbl_gc;

  int lphase;

  ngl  = ERI_ngl(solver);
  lmax = ERI_lmax_gl(solver);
  jmax = lmax*lmax;

  glx  = solver->glq_x;
  glw  = solver->glq_w;
  
  gtbl_n    = ERI_Gtbl_n(solver->gtbl);
  gtbl_j1j2 = ERI_Gtbl_j1j2(solver->gtbl);
  gtbl_gc   = ERI_Gtbl_gc(solver->gtbl);

  ngl  = ERI_ngl(solver);
  lmax = ERI_lmax_gl(solver);
  jmax = lmax*lmax;

  glx  = solver->glq_x;
  glw  = solver->glq_w;
  
  /* minimal set of R */
  num_mR = 0;
  for (i=0; i<numR; i++) {
    for (j=0; j<num_mR; j++) {
      if (fabs(R[mR2iR[j]]-R[i])<1e-10) { break; }
    }
    minimalR[i] = j;
    if (j==num_mR) {  
      if (num_mR>=MAXNUMR) {
        fprintf(stderr, "***ERROR in %s (%d)\n", __FILE__, __LINE__);
        fprintf(stderr, "   MAXNUMR is too small!\n");
        abort();
      }
      mR2iR[num_mR] = i;
      num_mR++;
    }
  }
  *num_minimalR = num_mR;

  /* spherical Bessel function */
  for (imR=0; imR<num_mR; imR++) {
    R0 = R[mR2iR[imR]];
    for (ik=0; ik<ngl; ik++) {
      k = glx[ik];
      ERI_Spherical_Bessel(k*R0, lmax, tmp_sb, tmp_dsb);
      if (omega<0.0) {
        scr = 1.0; /* no screening */
      } else {
        scr = 1.0-exp(-k*k*0.25/omega/omega); /* screening */
      }
      for (l=0; l<lmax; l++) {
        prej[(imR*lmax+l)*ngl+ik] = glw[ik]*exp(k)*tmp_sb[l]*scr;
      }
    }
  }

  /* Spherical Harominics */ 
  for (iR=0; iR<numR; iR++) {
    for (j=0; j<jmax; j++) {
      l = j2l[j];
      m = j2m[j];
      ERI_Real_Spherical_Harmonics(l, m, theta[iR], phi[iR], sh, dsht, dshp);
      preY[iR*jmax+j] = sh[0];
    }
  }

  igtbl0 = 0;
  for (j=0; j<jmax; j++) {
    l = j2l[j];
    n = gtbl_n[j];
    for (j1=0; j1<jmax; j1++) {
      l1 = j2l[j1];
      nmul = 0;  
      for (i=0; i<n; i++) {
        igtbl = igtbl0 + i;
        if (j1 != ERI_GTBL_UNPACK_J1( gtbl_j1j2[igtbl] )) { continue; }

        j2 = ERI_GTBL_UNPACK_J2( gtbl_j1j2[igtbl] ); 
        if (j2>=jmax) { continue; }

        l2 = j2l[j2];
        lphase = l2+l-l1;
        if (1==lphase%2) { continue; }
          
        mul_j2[(j*jmax+j1)*jmax+nmul] = j2;
        if (2==lphase%4) { 
          mul_gc[(j*jmax+j1)*jmax+nmul] = -gtbl_gc[igtbl];
        } else {
          mul_gc[(j*jmax+j1)*jmax+nmul] = gtbl_gc[igtbl];
        }
        nmul++;
      } /* i */
      mul_n [j*jmax+j1] = nmul;
    } /* j1 */
    igtbl0 += n;
  } /* j */
}




/*----------------------------------------------------------------------
  ERI_Integral_GL

  If you need the derivatives, you should call ERI_I4_d.
----------------------------------------------------------------------*/
void ERI_Integral_GL_Post(
  ERI_t        *solver,
  double       *I4,    /* (OUT) [numR] */
  const double *F1,    /* (IN) Overlap matrix */
  const double *F2,    /* (IN) */
  int           numR,
  double       *prej,   /* [lmax*ngl*numR] */
  double       *preY,   /* [numR*jmax1] */
  int          *mul_j2, /* [jmax1*jmax1*jmax1] */
  double       *mul_gc, /* [jmax1*jmax1*jmax1] */
  int          *mul_n,  /* [jmax1*jmax1] */
  int          *minimalR,    /* [numR] */
  int           num_minimalR
)
{
  int ngl, lmax, jmax, nmul, lphase;
  int iR, imR, i, j, j1, j2, ik, l, l1, l2;
  double gnt, sh, sum;
  const int *j2l = g_itbl_j2l;
  double int1[MAXNUMR];

#if SUMSKIP
  int ndglf;
  const double *fmax1, *fmax2;
#endif

  ngl  = ERI_ngl(solver);
  lmax = ERI_lmax_gl(solver);
  jmax = lmax*lmax;

#if SUMSKIP
  ndglf = ngl*jmax;
  fmax1 = &F1[ndglf];
  fmax2 = &F2[ndglf];
#endif

  for (iR=0; iR<numR; iR++) { I4[iR]=0.0; }

  for (j=0; j<jmax; j++) {
    l = j2l[j];
    for (imR=0; imR<num_minimalR; imR++) { int1[imR]=0.0; }
       
    for (j1=0; j1<jmax; j1++) {
      nmul =  mul_n[j*jmax+j1];
      for (i=0; i<nmul; i++) {
        gnt = mul_gc[(j*jmax+j1)*jmax+i];
        j2  = mul_j2[(j*jmax+j1)*jmax+i];
#if SUMSKIP
        if (fmax1[j1]*fmax2[j2]<SUMSKIP_THRESHOLD) { continue; }
#endif

#if LOOP_UNROLLING
        for (imR=0; imR<num_minimalR; imR++) {
          sum = 0.0;
          for (ik=0; ik<ngl-3; ik+=4) {
            sum += prej[((imR)*lmax+l)*ngl+ik+0]
                    *F1[j1*ngl+ik+0]*F2[j2*ngl+ik+0]
                 + prej[((imR)*lmax+l)*ngl+ik+1]
                    *F1[j1*ngl+ik+1]*F2[j2*ngl+ik+1]
                 + prej[((imR)*lmax+l)*ngl+ik+2]
                    *F1[j1*ngl+ik+2]*F2[j2*ngl+ik+2]
                 + prej[((imR)*lmax+l)*ngl+ik+3]
                    *F1[j1*ngl+ik+3]*F2[j2*ngl+ik+3];
          }
          for (; ik<ngl; ik++) {
            sum += prej[((imR)*lmax+l)*ngl+ik]
                    *F1[j1*ngl+ik]*F2[j2*ngl+ik];
          }
          int1[imR] += sum*gnt;
        } /* loop of imR */ 
#else 
        for (imR=0; imR<num_minimalR; imR++) {
          sum = 0.0;
          for (ik=0; ik<ngl; ik++) {
            sum += prej[(imR*lmax+l)*ngl+ik]
                    *F1[j1*ngl+ik]*F2[j2*ngl+ik];
          }
          int1[imR] += sum*gnt;
        } /* loop of imR */ 
#endif 
      } /* loop of i */
    } /* loop of j1 */
    
    for (iR=0; iR<numR; iR++) {
      imR = minimalR[iR];
      I4[iR] += 32.0*PI*int1[imR]*preY[iR*jmax+j];
    }
  } /* loop of j */ 
}



/*----------------------------------------------------------------------
  ERI_Integral_GL

  If you need the derivatives, you should call ERI_I4_d.
----------------------------------------------------------------------*/
void ERI_Integral_GL_Post2(
  ERI_t        *solver,
  double       *I4,    /* (OUT) [numR] */
  const double *F1,    /* (IN) Overlap matrix */
  const double *F2,    /* (IN) */
  int           numR,
  double       *prej,   /* [lmax*ngl*numR] */
  double       *preY,   /* [numR*jmax1] */
  int          *mul_j2, /* [jmax1*jmax1*jmax1] */
  double       *mul_gc, /* [jmax1*jmax1*jmax1] */
  int          *mul_n,  /* [jmax1*jmax1] */
  int          *minimalR,    /* [numR] */
  int           num_minimalR
)
{
  int ngl, lmax, jmax, nmul, lphase;
  int iR, imR, i, j, j1, j2, ik, l, l1, l2;
  double gnt, sh, sum;
  const int *j2l = g_itbl_j2l;
  double int1[MAXNUMR*ERI_LMAXMAX*ERI_LMAXMAX];

#if SUMSKIP
  int ndglf;
  const double *fmax1, *fmax2;
#endif

  ngl  = ERI_ngl(solver);
  lmax = ERI_lmax_gl(solver);
  jmax = lmax*lmax;

#if SUMSKIP
  ndglf = ngl*jmax;
  fmax1 = &F1[ndglf];
  fmax2 = &F2[ndglf];
#endif

  for (iR=0; iR<numR; iR++) { I4[iR]=0.0; }
  for (i=0; i<MAXNUMR*ERI_LMAXMAX*ERI_LMAXMAX; i++) { int1[i]=0.0; }

  for (j=0; j<jmax; j++) {
    l = j2l[j];
    for (j1=0; j1<jmax; j1++) {
      nmul =  mul_n[j*jmax+j1];
      for (i=0; i<nmul; i++) {
        gnt = mul_gc[(j*jmax+j1)*jmax+i];
        j2  = mul_j2[(j*jmax+j1)*jmax+i];
#if SUMSKIP
        if (fmax1[j1]*fmax2[j2]<SUMSKIP_THRESHOLD) { continue; }
#endif

        for (imR=0; imR<num_minimalR; imR++) {
          sum = 0.0;
          for (ik=0; ik<ngl; ik++) {
            sum += prej[(imR*lmax+l)*ngl+ik]
                    *F1[j1*ngl+ik]*F2[j2*ngl+ik];
          }
          int1[j*numR+imR] += sum*gnt;
        } /* loop of imR */ 
      } /* loop of i */
    } /* loop of j1 */
  } /* loop of l */
  
  for (j=0; j<jmax; j++) {
    for (iR=0; iR<numR; iR++) {
      imR = minimalR[iR];
      I4[iR] += 32.0*PI*int1[j*numR+imR]*preY[iR*jmax+j];
    }
  } /* loop of j */ 
}



/*----------------------------------------------------------------------
  ERI_Integral_GL_X

  If you need the derivatives, you should call ERI_I4_d.
----------------------------------------------------------------------*/
void ERI_Integral_GL_X(
  ERI_t        *solver,
  double       *X,      /* (OUT) [numR*lmax] */
  const double *F1,    /* (IN) Overlap matrix */
  const double *F2,    /* (IN) */
  int           numR,
  double       *prej,   /* [lmax*ngl*numR] */
  int          *mul_j2,      /* [jmax1*jmax1*jmax1] */
  double       *mul_gc,      /* [jmax1*jmax1*jmax1] */
  int          *mul_n,       /* [jmax1*jmax1] */
  int          *minimalR,    /* [numR] */
  int           num_minimalR
)
{
  int ngl, lmax, jmax, nmul, lphase;
  int iR, imR, i, j, j1, j2, ik, l, l1, l2;
  double gnt, sum;
  const int *j2l = g_itbl_j2l;

#if SUMSKIP
  int ndglf;
  const double *fmax1, *fmax2;
#endif

  ngl  = ERI_ngl(solver);
  lmax = ERI_lmax_gl(solver);
  jmax = lmax*lmax;

#if SUMSKIP
  ndglf = ngl*jmax;
  fmax1 = &F1[ndglf];
  fmax2 = &F2[ndglf];
#endif

  
  for (i=0; i<numR*lmax; i++) { X[i]=0.0; }

  for (j=0; j<jmax; j++) {
    l = j2l[j];
    for (j1=0; j1<jmax; j1++) {
      nmul =  mul_n[j*jmax+j1];
      for (i=0; i<nmul; i++) {
        gnt = mul_gc[(j*jmax+j1)*jmax+i];
        j2  = mul_j2[(j*jmax+j1)*jmax+i];
#if SUMSKIP
        if (fmax1[j1]*fmax2[j2]<SUMSKIP_THRESHOLD) { continue; }
#endif

        for (imR=0; imR<num_minimalR; imR++) {
          sum = 0.0;
          for (ik=0; ik<ngl; ik++) {
            sum += prej[(imR*lmax+l)*ngl+ik]
                    *F1[j1*ngl+ik]*F2[j2*ngl+ik];
          }
          X[l*numR+imR] += sum*gnt;
        } /* loop of imR */ 
      } /* loop of i */
    } /* loop of j1 */
  } /* loop of j */ 
}


/*----------------------------------------------------------------------
  ERI_Integral_GL

  If you need the derivatives, you should call ERI_I4_d.
----------------------------------------------------------------------*/
void ERI_Integral_GL_X_Post(
  ERI_t        *solver,
  double       *I4,       /* (OUT) [numR] */
  const double *X,        /* (IN)  [numR*lmax1] */
  int           numR,
  const double *preY,     /* (IN)  [numR*jmax1] */
  const int    *minimalR  /* (IN)  [numR] */
)
{
  int lmax, jmax;
  int iR, imR, j, l;
  double sum;
  const int *j2l   = g_itbl_j2l;

  lmax = ERI_lmax_gl(solver);
  jmax = lmax*lmax;

  for (iR=0; iR<numR; iR++) {
    imR = minimalR[iR];
    sum = 0.0;
    for (j=0; j<jmax; j++) {
      l = j2l[j];
      sum += X[l*numR+imR]*preY[iR*jmax+j];
    } /* loop of j */
    I4[iR] = 32.0*PI*sum;
  } /* loop of iR */ 
}

/* EOF */
