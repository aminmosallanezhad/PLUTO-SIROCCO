/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Roe Riemann solver for the HD equations.

  Solve the Riemann problem for the Euler equations of gasdynamics
  using the standard Roe solver with local characteristic decomposition.
  Eigenvectors are identical to the ones given in the book by Toro 
  and can also be derived from the maple script "eigenv.maple"
  in Src/HD/.
  The solver can be used for adiabatic and isothermal hydrodynamics.
  
  The macro ROE_AVERAGE specifies how the averaging process is done:
    - ROE_AVERAGE == YES use Roe average (default);
    - ROE_AVERAGE == NO  use arithmetic average.

  The macro CHECK_ROE_MATRIX can be used to verify that the 
  characteristic decomposition reproduces the Roe matrix.

  On input, it takes left and right primitive state vectors 
  \c stateL->v and \c stateR->v at zone edge \c i+1/2;
  On output, return flux and pressure vectors at the same interface 
  \c i+1/2 (note that the \c i refers to \c i+1/2).
  
  Also during this step, compute maximum wave propagation speed (cmax) 
  for  explicit time step computation.
   
  \b Reference:
   -   "Riemann Solver and Numerical Methods for Fluid Dynamics"
        by E.F. Toro (Chapter 11)
        
  \authors A. Mignone (andrea.mignone@unito.it)
  \date    June 27, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

#define ROE_AVERAGE       YES
#define CHECK_ROE_MATRIX  NO

/* ********************************************************************* */
void Roe_Solver (const Sweep *sweep, int beg, int end, 
                 double *cmax, Grid *grid)
/*!
 * Solve the Riemann problem using the Roe solver.
 *
 * \param[in,out] sweep   pointer to Sweep structure
 * \param[in]     beg     initial grid index
 * \param[out]    end     final grid index
 * \param[out]    cmax    1D array of maximum characteristic speeds
 * \param[in]     grid    pointer to array of Grid structures.
 *
 *********************************************************************** */
{
  int   nv, i, j, k, nn;

  const State   *stateL = &(sweep->stateL);
  const State   *stateR = &(sweep->stateR);

  double **fL = stateL->flux, **fR = stateR->flux;
  double *a2L = stateL->a2,   *a2R = stateR->a2;
  double  *pL = stateL->prs,   *pR = stateR->prs;

  double  scrh;
  double  um[NFLX], vel2;
  double  a2, a, h;
  double  dv[NFLX], eta[NFLX];
  double  Rc[NFLX][NFLX], alambda[NFLX], lambda[NFLX];
  double  delta, delta_inv, gmm1, gmm1_inv;
#if ROE_AVERAGE == YES
  double s, c, hl, hr;
#endif
  double *ql, *qr, *uL, *uR;

  double bmin, bmax, scrh1;
  double Us[NFLX];

  delta     = 1.e-7;
  delta_inv = 1.0/delta;
  #if EOS == IDEAL
   gmm1      = g_gamma - 1.0;
   gmm1_inv  = 1.0/gmm1;
  #endif

  for (i = NFLX; i--;  ) {
  for (j = NFLX; j--;  ) {
    Rc[i][j] = 0.0;
  }}

/* ----------------------------------------------------
     compute sound speed & fluxes at zone interfaces
   ---------------------------------------------------- */

  SoundSpeed2 (stateL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (stateR, beg, end, FACE_CENTER, grid);

  Flux (stateL, beg, end);
  Flux (stateR, beg, end);

  for (i = beg; i <= end; i++)  {

    uR = stateR->u[i];
    uL = stateL->u[i];

#if SHOCK_FLATTENING == MULTID   

    /* ---------------------------------------------
       HLL switching function as in Quirk (1994).
       Since the problem is related to multidimensional 
       pathologies, it works in more than 1-D only.	 
       Use the HLL flux function if the interface 
       lies within a strong shock.
       The effect of this switch is visible
       in the Mach reflection test.
      --------------------------------------------- */

    if ((sweep->flag[i] & FLAG_HLL) || (sweep->flag[i+1] & FLAG_HLL)){        
      HLL_Speed (stateL, stateR, &bmin - i, &bmax - i, i, i);
      a     = MAX(fabs(bmin), fabs(bmax));
      cmax[i] = a;
      bmin  = MIN(0.0, bmin);
      bmax  = MAX(0.0, bmax);
      scrh  = 1.0/(bmax - bmin);
      for (nv = NFLX; nv--; ){
        sweep->flux[i][nv]  = bmin*bmax*(uR[nv] - uL[nv])
                           +  bmax*fL[i][nv] - bmin*fR[i][nv];
        sweep->flux[i][nv] *= scrh;
      }
      sweep->press[i] = (bmax*pL[i] - bmin*pR[i])*scrh;
      continue;
    }
#endif

    ql = stateL->v[i];
    qr = stateR->v[i];

  /*  ----  Define Wave Jumps  ----  */

    for (nv = NFLX; nv--;   ) dv[nv] = qr[nv] - ql[nv];

    #if ROE_AVERAGE == YES    
    s       = sqrt(qr[RHO]/ql[RHO]);
    um[RHO]  = ql[RHO]*s;
    s       = 1.0/(1.0 + s); 
    c       = 1.0 - s;
  
    um[VX1] = s*ql[VX1] + c*qr[VX1];
    um[VX2] = s*ql[VX2] + c*qr[VX2];
    um[VX3] = s*ql[VX3] + c*qr[VX3];

    #if EOS == IDEAL
    vel2 = um[VX1]*um[VX1] + um[VX2]*um[VX2] + um[VX3]*um[VX3];

    hl  = 0.5*(ql[VX1]*ql[VX1] + ql[VX2]*ql[VX2] + ql[VX3]*ql[VX3]);    
    hl += a2L[i]*gmm1_inv;
     
    hr = 0.5*(qr[VX1]*qr[VX1] + qr[VX2]*qr[VX2] + qr[VX3]*qr[VX3]);
    hr += a2R[i]*gmm1_inv;

    h = s*hl + c*hr;

  /* -------------------------------------------------
       the following should be  equivalent to 
    
       scrh = dv[VX1]*dv[VX1] + dv[VX2]*dv[VX2] + dv[VX3]*dv[VX3];

       a2 = s*a2L + c*a2R + 0.5*gmm1*s*c*scrh;

       and therefore always positive.
       just work out the coefficients...
     -------------------------------------------------- */
     
      a2 = gmm1*(h - 0.5*vel2);
      a  = sqrt(a2);
    #endif /* EOS == IDEAL */
    #else
    for (nv = NFLX; nv--;   ) um[nv] = 0.5*(ql[nv] + qr[nv]);  
    #if EOS == IDEAL
    a2   = g_gamma*um[PRS]/um[RHO];
    a    = sqrt(a2);
     
    vel2 = um[VX1]*um[VX1] + um[VX2]*um[VX2] + um[VX3]*um[VX3];
    h    = 0.5*vel2 + a2/gmm1;
    #endif /* EOS == IDEAL */
    #endif /* ROE_AVERAGE == YES/NO */
  
    #if EOS == ISOTHERMAL
    a2 = 0.5*(a2L[i] + a2R[i]);
    a  = sqrt(a2);
    #endif

  /* ----------------------------------------------------------------
      define non-zero components of conservative eigenvectors Rc, 
      eigenvalues (lambda) and wave strenght eta = L.du     
     ----------------------------------------------------------------  */

  /*  ---- (u - c_s)  ----  */ 

    nn         = 0;
    lambda[nn] = um[VXn] - a;
    #if EOS == IDEAL
    eta[nn] = 0.5/a2*(dv[PRS] - dv[VXn]*um[RHO]*a);
    #elif EOS == ISOTHERMAL
    eta[nn] = 0.5*(dv[RHO] - um[RHO]*dv[VXn]/a);
    #endif

    Rc[RHO][nn] = 1.0;
    Rc[MXn][nn] = um[VXn] - a;
    Rc[MXt][nn] = um[VXt]; 
    Rc[MXb][nn] = um[VXb];
    #if EOS == IDEAL
    Rc[ENG][nn] = h - um[VXn]*a;
    #endif

  /*  ---- (u + c_s)  ----  */ 

    nn         = 1;
    lambda[nn] = um[VXn] + a;
    #if EOS == IDEAL
    eta[nn]    = 0.5/a2*(dv[PRS] + dv[VXn]*um[RHO]*a);
    #elif EOS == ISOTHERMAL
    eta[nn] = 0.5*(dv[RHO] + um[RHO]*dv[VXn]/a);
    #endif

    Rc[RHO][nn] = 1.0;
    Rc[MXn][nn] = um[VXn] + a;
    Rc[MXt][nn] = um[VXt];
    Rc[MXb][nn] = um[VXb];
    #if EOS == IDEAL
     Rc[ENG][nn] = h + um[VXn]*a;
    #endif

  /*  ----  (u)  ----  */ 
     
    #if EOS == IDEAL
    nn         = 2;
    lambda[nn] = um[VXn];
    eta[nn]    = dv[RHO] - dv[PRS]/a2;
    Rc[RHO][nn] = 1.0;
    Rc[MX1][nn] = um[VX1];
    Rc[MX2][nn] = um[VX2];
    Rc[MX3][nn] = um[VX3];
    Rc[ENG][nn] = 0.5*vel2;
    #endif
    
  /*  ----  (u)  ----  */ 

    nn++;
    lambda[nn] = um[VXn];
    eta[nn]    = um[RHO]*dv[VXt];
    Rc[MXt][nn] = 1.0;
    #if EOS == IDEAL
    Rc[ENG][nn] = um[VXt];  
    #endif

  /*  ----  (u)  ----  */ 

    nn++;
    lambda[nn] = um[VXn];
    eta[nn]    = um[RHO]*dv[VXb];
    Rc[MXb][nn] = 1.0;
    #if EOS == IDEAL
    Rc[ENG][nn] = um[VXb];  
    #endif

  /*  ----  get max eigenvalue  ----  */

    cmax[i] = fabs(um[VXn]) + a;
    g_maxMach = MAX(fabs(um[VXn]/a), g_maxMach);

    #if DIMENSIONS > 1
   
    /* ---------------------------------------------
         use the HLL flux function if the interface 
         lies within a strong shock.
         The effect of this switch is visible
         in the Mach reflection test.
      --------------------------------------------- */

      #if EOS == IDEAL
      scrh  = fabs(ql[PRS] - qr[PRS]);
      scrh /= MIN(ql[PRS],qr[PRS]);
      #elif EOS == ISOTHERMAL
      scrh  = fabs(ql[RHO] - qr[RHO]);
      scrh /= MIN(ql[RHO],qr[RHO]);
      scrh *= a*a;
      #endif
      if (scrh > 0.5 && (qr[VXn] < ql[VXn])){   /* -- tunable parameter -- */
        bmin = MIN(0.0, lambda[0]);
        bmax = MAX(0.0, lambda[1]);
        scrh1 = 1.0/(bmax - bmin);
        for (nv = NFLX; nv--;   ){
         sweep->flux[i][nv]  = bmin*bmax*(uR[nv] - uL[nv])
                           +   bmax*fL[i][nv] - bmin*fR[i][nv];
         sweep->flux[i][nv] *= scrh1;
        }
        sweep->press[i] = (bmax*pL[i] - bmin*pR[i])*scrh1;
        continue;
      } 
    #endif  /* DIMENSIONS > 1 */

    #if CHECK_ROE_MATRIX == YES
    for (nv = 0; nv < NFLX; nv++){
      um[nv] = 0.0;
      for (k = 0; k < NFLX; k++){
      for (j = 0; j < NFLX; j++){
        um[nv] += Rc[nv][k]*(k==j)*lambda[k]*eta[j];
      }}
    }
    for (nv = 0; nv < NFLX; nv++){
      scrh = fR[i][nv] - fL[i][nv] - um[nv];
      if (nv == MXn) scrh += pR[i] - pL[i];
      if (fabs(scrh) > 1.e-6){
        printLog ("! Matrix condition not satisfied %d, %12.6e\n", nv, scrh);
        Show(sweep->vL, i);
        Show(sweep->vR, i);
        exit(1);
      }
    }
    #endif

  /* -----------------------------------------------------------
                      compute Roe flux 
     ----------------------------------------------------------- */
      
    for (nv = NFLX; nv--;   ) alambda[nv]  = fabs(lambda[nv]);

  /*  ----  entropy fix  ----  */

    if (alambda[0] <= delta) {
      alambda[0] = 0.5*lambda[0]*lambda[0]/delta + 0.5*delta;
    }
    if (alambda[1] <= delta) {
      alambda[1] = 0.5*lambda[1]*lambda[1]/delta + 0.5*delta;
    }

    for (nv = NFLX; nv--;   ) {
      sweep->flux[i][nv] = fL[i][nv] + fR[i][nv];
      for (k  = NFLX; k-- ;   ) {
        sweep->flux[i][nv] -= alambda[k]*eta[k]*Rc[nv][k];
      }
      sweep->flux[i][nv] *= 0.5;
    }
    sweep->press[i] = 0.5*(pL[i] + pR[i]);
  }
}
#undef ROE_AVERAGE        
#undef CHECK_ROE_MATRIX
