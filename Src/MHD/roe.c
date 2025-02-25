/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Implementation of the Roe Riemann solver for the MHD equations.

  Solve the Riemann problem for the adiabatic and isothermal MHD 
  equations using the linearized Riemann solver of Roe.
  The implementation follows the approach outlined by 
  Cargo & Gallice (1997). 
  
  The isothermal version is recovered by taking the limit of 
  \f$ \bar{a}^2 \f$ for \f$ \gamma \to 1\f$ which gives 
  (page 451) \f$\bar{a}^2 \to {\rm g\_isoSoundSpeed2} + X\f$
  where \c X is defined as in the adiabatic case.
  This follows by imposing zero jump across the entropy wave 
  (first of Eq. 4.20, page 452) giving 
  \f$ \Delta p = (\bar{a}^2 - X)\Delta\rho\f$.
  Then all the terms like 
  \f$ [X\Delta\rho + \Delta p] \to \bar{a}^2\Delta\rho \f$.
 
  When the macro HLL_HYBRIDIZATION is set to YES, we revert to the HLL 
  solver whenever an unphysical state appear in the solution.
  
  The macro CHECK_ROE_MATRIX can be used to verify that the 
  characteristic decomposition reproduces the Roe matrix.
  This can be done also when BACKGROUND_FIELD is set to YES.

  On input, this function takes left and right primitive state vectors 
  \c stateL->v and \c stateR->v at zone edge i+1/2;
  On output, return flux and pressure vectors at the same interface 
  \c i+1/2 (note that the \c i refers to \c i+1/2).
  
  Also during this step, compute maximum wave propagation speed (cmax) 
  for  explicit time step computation.
   
  \b Reference:
    - "Roe Matrices for Ideal MHD and Systematic Construction 
       of Roe Matrices for Systems of Conservation Laws",
       P. Cargo, G. Gallice, JCP (1997), 136, 446.
       
  \authors A. Mignone (andrea.mignone@unito.it)
           C. Zanni   (zanni@oato.inaf.it)
  \date    Oct 24, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

#define sqrt_1_2  (0.70710678118654752440)
#define HLL_HYBRIDIZATION    NO
#ifndef CHECK_ROE_MATRIX
  #define CHECK_ROE_MATRIX     NO
#endif  

/* ********************************************************************* */
void Roe_Solver (const Sweep *sweep, int beg, int end, 
                 double *cmax, Grid *grid)
/*!
 * Solve Riemann problem for the adiabatic MHD equations using the 
 * Roe Riemann solver of Cargo & Gallice (1997).
 * 
 * \param[in,out] sweep   pointer to Sweep structure
 * \param[in]     beg     initial grid index
 * \param[out]    end     final grid index
 * \param[out]    cmax    1D array of maximum characteristic speeds
 * \param[in]     grid    pointer to array of Grid structures.
 *
 *********************************************************************** */
{
  int  nv, i, j, k;
  int  ifail;

  const State   *stateL = &(sweep->stateL);
  const State   *stateR = &(sweep->stateR);

  double rho, u, v, w, vel2, bx, by, bz, pr;
  double a2, a, ca2, cf2, cs2;
  double cs, ca, cf, b2;
  double alpha_f, alpha_s, beta_y, beta_z, beta_v, scrh, sBx;
  double dV[NFLX], dU[NFLX], *vL, *vR, *uL, *uR;
  double *SL, *SR;
  double Rc[NFLX][NFLX], eta[NFLX], lambda[NFLX];
  double alambda[NFLX], Uv[NFLX];

  double sqrt_rho;
  double delta, delta_inv;
  
  double g1, sl, sr, H, Hgas, HL, HR, Bx, By, Bz, X;
  double vdm, BdB, beta_dv, beta_dB;
  double bt2, Btmag, sqr_rho_L, sqr_rho_R;

  double **bgf;
  #if BACKGROUND_FIELD == YES
  double B0x, B0y, B0z, B1x, B1y, B1z;
  #endif      
  double Us[NFLX];
  double **fL = stateL->flux, **fR = stateR->flux;
  double *a2L = stateL->a2,   *a2R = stateR->a2;
  double  *pL = stateL->prs,   *pR = stateR->prs;
  
  delta    = 1.e-6;

/* --------------------------------------------------------
   0. Get background magnetic field
   -------------------------------------------------------- */

#if BACKGROUND_FIELD == YES
  GetBackgroundField (stateL, beg, end, FACE_CENTER, grid);
  bgf = stateL->Bbck;
  #if (DIVB_CONTROL == EIGHT_WAVES)
  print ("! Roe_Solver(): background field and 8wave not tested\n");
  QUIT_PLUTO(1);
  #endif
#endif

/* --------------------------------------------------------
   1. GLM pre-Rieman solver
   -------------------------------------------------------- */
   
#ifdef GLM_MHD
  GLM_Solve (sweep, beg, end, grid);
#endif

/* --------------------------------------------------------
   2. Compute sound speed & fluxes at zone interfaces
   -------------------------------------------------------- */

  SoundSpeed2 (stateL, beg, end, FACE_CENTER, grid);
  SoundSpeed2 (stateR, beg, end, FACE_CENTER, grid);

  Flux (stateL, beg, end);
  Flux (stateR, beg, end);

  SL = sweep->SL; SR = sweep->SR;

#if EOS == IDEAL
  g1 = g_gamma - 1.0;
#endif

/* --------------------------------------------------------
   3. Some eigenvectors components will always
      be zero. Set set Rc = 0 initially  
   -------------------------------------------------------- */
     
  for (k = NFLX; k--;  ) {
  for (j = NFLX; j--;  ) {
    Rc[k][j] = 0.0;
  }}

/* --------------------------------------------------------
   4. Begin main loop
   -------------------------------------------------------- */
   
  for (i = beg; i <= end; i++) {

    vL = stateL->v[i]; uL = stateL->u[i];
    vR = stateR->v[i]; uR = stateR->u[i];

  /* ------------------------------------------------------
     4a.  switch to HLL in proximity of strong shock 
     ------------------------------------------------------ */

    #if SHOCK_FLATTENING == MULTID
    if ((sweep->flag[i] & FLAG_HLL) || (sweep->flag[i+1] & FLAG_HLL)){
      HLL_Speed (stateL, stateR, SL, SR, i, i);

      scrh = MAX(fabs(SL[i]), fabs(SR[i]));
      cmax[i] = scrh;

      if (SL[i] > 0.0) {
        NFLX_LOOP(nv) sweep->flux[i][nv] = fL[i][nv];
        sweep->press[i] = pL[i];
      } else if (SR[i] < 0.0) {
        NFLX_LOOP(nv) sweep->flux[i][nv] = fR[i][nv];
        sweep->press[i] = pR[i];
      }else{
        scrh = 1.0/(SR[i] - SL[i]);
        NFLX_LOOP(nv) {
          sweep->flux[i][nv]  = SR[i]*SL[i]*(uR[nv] - uL[nv])
                             +  SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
          sweep->flux[i][nv] *= scrh;
        }
        sweep->press[i] = (SR[i]*pL[i] - SL[i]*pR[i])*scrh;
      }
      continue;
    }
    #endif

  /* ----------------------------------------------------------------
     4b. Compute jumps in primitive and conservative variables.
      
         Note on the velocity jump:
         Formally the jump in velocity should be written using
         conservative variables:
      
         Delta(v) = Delta(m/rho) = (Delta(m) - v*Delta(rho))/rho
      
         with rho and v defined by the second and first of [4.6]
         Still, since we reconstruct v and not m, the following MAPLE 
         script shows that the previous expression simplifies to
         Delta(v) = vR - vL.
      
         restart;
         sR   := sqrt(rho[R])/(sqrt(rho[R]) + sqrt(rho[L])):
         sL   := sqrt(rho[L])/(sqrt(rho[R]) + sqrt(rho[L])):
         rhoa := sL*rho[R] + sR*rho[L]:
         va   := sL*v[L]   + sR*v[R]:
         dV := (rho[R]*v[R] - rho[L]*v[L] - va*(rho[R]-rho[L]))/rhoa:
         simplify(dV);
     ---------------------------------------------------------------- */

    NFLX_LOOP(nv) { 
      dV[nv] = vR[nv] - vL[nv];
      dU[nv] = uR[nv] - uL[nv];
    }

  /* ---------------------------------
     4c. Compute Roe averages 
     --------------------------------- */

    sqr_rho_L = sqrt(vL[RHO]);
    sqr_rho_R = sqrt(vR[RHO]);

    sl = sqr_rho_L/(sqr_rho_L + sqr_rho_R);
    sr = sqr_rho_R/(sqr_rho_L + sqr_rho_R);
/*      sl = sr = 0.5;    */
    
    rho = sr*vL[RHO] + sl*vR[RHO];

    sqrt_rho = sqrt(rho);

    u = sl*vL[VXn] + sr*vR[VXn];
    v = sl*vL[VXt] + sr*vR[VXt];
    w = sl*vL[VXb] + sr*vR[VXb];

    Bx = sr*vL[BXn] + sl*vR[BXn];
    By = sr*vL[BXt] + sl*vR[BXt];
    Bz = sr*vL[BXb] + sl*vR[BXb];

    #if BACKGROUND_FIELD == YES
   /* -- Define field B0 and total B. B1 is the deviation -- */   
    B0x = bgf[i][BXn-BX1]; B1x = sr*vL[BXn] + sl*vR[BXn]; Bx = B0x + B1x;
    B0y = bgf[i][BXt-BX1]; B1y = sr*vL[BXt] + sl*vR[BXt]; By = B0y + B1y;
    B0z = bgf[i][BXb-BX1]; B1z = sr*vL[BXb] + sl*vR[BXb]; Bz = B0z + B1z;
    #else
    Bx = sr*vL[BXn] + sl*vR[BXn];
    By = sr*vL[BXt] + sl*vR[BXt];
    Bz = sr*vL[BXb] + sl*vR[BXb];
    #endif

    sBx = (Bx >= 0.0 ? 1.0 : -1.0);

    bx = Bx/sqrt_rho;
    by = By/sqrt_rho;
    bz = Bz/sqrt_rho;
    
    bt2   = by*by + bz*bz;
    b2    = bx*bx + bt2;
    Btmag = sqrt(bt2*rho);

    X  = dV[BXn]*dV[BXn] + dV[BXt]*dV[BXt] + dV[BXb]*dV[BXb];
    X /= (sqr_rho_L + sqr_rho_R)*(sqr_rho_L + sqr_rho_R)*2.0;   

    vdm = u*dU[MXn] + v*dU[MXt] + w*dU[MXb];
    #if BACKGROUND_FIELD == YES /* BdB = B1.dB1 (deviation only) */
    BdB = B1x*dU[BXn] + B1y*dU[BXt] + B1z*dU[BXb]; 
    #else
    BdB = Bx*dU[BXn] + By*dU[BXt] + Bz*dU[BXb];
    #endif
    
  /* ---------------------------------------
     4d. Compute enthalpy and sound speed.
     --------------------------------------- */

    #if EOS == ISOTHERMAL 
    a2 = 0.5*(a2L[i] + a2R[i]) + X;  /* in most cases a2L = a2R
                                         for isothermal MHD */
    #elif EOS == BAROTROPIC
    print ("! Roe_Solver(): not implemented for barotropic EOS\n");
    QUIT_PLUTO(1);
    #elif EOS == IDEAL
    vel2    = u*u + v*v + w*w;
    dV[PRS] = g1*((0.5*vel2 - X)*dV[RHO] - vdm + dU[ENG] - BdB); 
     
    HL   = (uL[ENG] + pL[i])/vL[RHO];
    HR   = (uR[ENG] + pR[i])/vR[RHO];
    H    = sl*HL + sr*HR;   /* total enthalpy */

    #if BACKGROUND_FIELD == YES
    scrh = B1x*Bx + B1y*By + B1z*Bz;
    Hgas = H - scrh/rho;   /* gas enthalpy */
    #else
    Hgas = H - b2;         /* gas enthalpy */
    #endif

    a2 = (2.0 - g_gamma)*X + g1*(Hgas - 0.5*vel2);
    if (a2 < 0.0) {
     printLog ("! Roe_Solver(): a2 = %12.6e < 0.0 !! \n",a2);
     a2 = sqrt(g_gamma*(vL[PRS] + vR[PRS])/(vL[RHO] + vR[RHO]));
     Show(stateL->v,i);
     Show(stateR->v,i);
     QUIT_PLUTO(1);
    }      
    #endif /* EOS == IDEAL */
    
  /* ------------------------------------------------------------
     4e. Compute fast and slow magnetosonic speeds.

      The following expression appearing in the definitions
      of the fast magnetosonic speed 
    
       (a^2 - b^2)^2 + 4*a^2*bt^2 = (a^2 + b^2)^2 - 4*a^2*bx^2

      is always positive and avoids round-off errors.
      
      Note that we always use the total field to compute the 
      characteristic speeds.
     ------------------------------------------------------------ */
        
    double dab2 = a2 - b2;
    ca2  = bx*bx;
    scrh = dab2*dab2 + 4.0*bt2*a2;    
    double sqDelta = sqrt(scrh);

    cf2 = 0.5*(a2 + b2 + sqDelta); 
    cs2 = a2*ca2/cf2;   /* -- same as 0.5*(a2 + b2 - scrh) -- */
    
    cf = sqrt(cf2);
    cs = sqrt(cs2);
    ca = sqrt(ca2);
    a  = sqrt(a2); 

    if (cf == cs) {
      alpha_f = 1.0;
      alpha_s = 0.0;
/*      print ("! Roe(): degenerate case\n ");   
      QUIT_PLUTO(1);                          */
    }else{ 
      alpha_f = 0.5*( dab2/sqDelta + 1.0);
      alpha_s = 0.5*(-dab2/sqDelta + 1.0);

      alpha_f = MAX(0.0, alpha_f);
      alpha_s = MAX(0.0, alpha_s);
      alpha_f = sqrt(alpha_f);
      alpha_s = sqrt(alpha_s);
    }

    if (bt2 > 1.e-18*b2) {
      beta_y = By/Btmag; 
      beta_z = Bz/Btmag;
    } else {
      beta_z = beta_y = sqrt_1_2;
    }

  /* -------------------------------------------------------------------
     4f. Compute non-zero entries of conservative eigenvectors (Rc), 
         wave strength L*dU (=eta) for all 8 (or 7) waves using the
         expressions given by Eq. [4.18]--[4.21]. 
         Fast and slow eigenvectors are multiplied by a^2 while
         jumps are divided by a^2.
      
         Notes:
         - the expression on the paper has a typo in the very last term 
           of the energy component: it should be + and not - !
         - with background field splitting: additional terms must be 
           added to the energy component for fast, slow and Alfven waves.
           To obtain energy element, conservative eigenvector (with 
           total field) must be multiplied by | 0 0 0 0 -B0y -B0z 1 |.
           Also, H - b2 does not give gas enthalpy. A term b0*btot must 
           be added and eta (wave strength) should contain total field 
           and deviation's delta.
     ------------------------------------------------------------------- */

  /* -----------------------
      Fast wave:  u - c_f
     ----------------------- */

    k = KFASTM;
    lambda[k] = u - cf;

    scrh    = alpha_s*cs*sBx;
    beta_dv = beta_y*dV[VXt] + beta_z*dV[VXb];
    beta_dB = beta_y*dV[BXt] + beta_z*dV[BXb];
    beta_v  = beta_y*v       + beta_z*w;

    Rc[RHO][k] = alpha_f;
    Rc[MXn][k] = alpha_f*lambda[k];
    Rc[MXt][k] = alpha_f*v + scrh*beta_y;
    Rc[MXb][k] = alpha_f*w + scrh*beta_z; 

    Rc[BXt][k] = alpha_s*a*beta_y/sqrt_rho;
    Rc[BXb][k] = alpha_s*a*beta_z/sqrt_rho;

    #if EOS == IDEAL
    Rc[ENG][k] =   alpha_f*(Hgas - u*cf) + scrh*beta_v
                 + alpha_s*a*Btmag/sqrt_rho;
    #if BACKGROUND_FIELD == YES
    Rc[ENG][k] -= B0y*Rc[BXt][k] + B0z*Rc[BXb][k];
    #endif

    eta[k] =   alpha_f*(X*dV[RHO] + dV[PRS]) + rho*scrh*beta_dv
             - rho*alpha_f*cf*dV[VXn]        + sqrt_rho*alpha_s*a*beta_dB;
    #elif EOS == ISOTHERMAL
    eta[k] =   alpha_f*(0.0*X + a2)*dV[RHO] + rho*scrh*beta_dv
             - rho*alpha_f*cf*dV[VXn] + sqrt_rho*alpha_s*a*beta_dB;
    #endif
    
    eta[k] *= 0.5/a2;

  /* -----------------------
      Fast wave:  u + c_f
     ----------------------- */

    k = KFASTP;
    lambda[k] = u + cf;

    Rc[RHO][k] = alpha_f;
    Rc[MXn][k] = alpha_f*lambda[k];
    Rc[MXt][k] = alpha_f*v - scrh*beta_y;
    Rc[MXb][k] = alpha_f*w - scrh*beta_z;

    Rc[BXt][k] = Rc[BXt][KFASTM];
    Rc[BXb][k] = Rc[BXb][KFASTM];

    #if EOS == IDEAL
    Rc[ENG][k] =   alpha_f*(Hgas + u*cf) - scrh*beta_v
                 + alpha_s*a*Btmag/sqrt_rho;

    #if BACKGROUND_FIELD == YES
    Rc[ENG][k] -= B0y*Rc[BXt][k] + B0z*Rc[BXb][k];
    #endif

    eta[k] =   alpha_f*(X*dV[RHO] + dV[PRS]) - rho*scrh*beta_dv
             + rho*alpha_f*cf*dV[VXn]        + sqrt_rho*alpha_s*a*beta_dB;
    #elif EOS == ISOTHERMAL
    eta[k] =   alpha_f*(0.0*X + a2)*dV[RHO] - rho*scrh*beta_dv
             + rho*alpha_f*cf*dV[VXn]      + sqrt_rho*alpha_s*a*beta_dB;
    #endif

    eta[k] *= 0.5/a2;

  /* -----------------------
      Entropy wave:  u
     ----------------------- */

    #if EOS == IDEAL
     k = KENTRP;
     lambda[k] = u;

     Rc[RHO][k] = 1.0;
     Rc[MXn][k] = u; 
     Rc[MXt][k] = v; 
     Rc[MXb][k] = w; 
     Rc[ENG][k] = 0.5*vel2 + (g_gamma - 2.0)/g1*X;

     eta[k] = ((a2 - X)*dV[RHO] - dV[PRS])/a2;
    #endif

  /* -----------------------------------------------------------------
      div.B wave (u): this wave exists when: 

       1) 8 wave formulation
       2) CT, since we always have 8 components, but it 
          carries zero jump.

      With GLM, KDIVB is replaced by KPSI_GLMM, KPSI_GLMP and these
      two waves should not enter in the Riemann solver (eta = 0.0) 
      since the 2x2 linear system formed by (B,psi) has already 
      been solved.
     ----------------------------------------------------------------- */

    #ifdef GLM_MHD
     lambda[KPSI_GLMP] =  glm_ch;
     lambda[KPSI_GLMM] = -glm_ch;
     eta[KPSI_GLMP] = eta[KPSI_GLMM] = 0.0;
    #else
     k = KDIVB;
     lambda[k] = u;
     #if DIVB_CONTROL == EIGHT_WAVES
      Rc[BXn][k] = 1.0;
      eta[k]    = dU[BXn];
     #else
      Rc[BXn][k] = eta[k] = 0.0;
     #endif
    #endif
    
   /* -----------------------
       Slow wave:  u - c_s
      ----------------------- */

    scrh = alpha_f*cf*sBx;
     
    k = KSLOWM;
    lambda[k] = u - cs;

    Rc[RHO][k] = alpha_s;
    Rc[MXn][k] = alpha_s*lambda[k];      
    Rc[MXt][k] = alpha_s*v - scrh*beta_y;
    Rc[MXb][k] = alpha_s*w - scrh*beta_z;

    Rc[BXt][k] = - alpha_f*a*beta_y/sqrt_rho;
    Rc[BXb][k] = - alpha_f*a*beta_z/sqrt_rho;

    #if EOS == IDEAL
    Rc[ENG][k] =   alpha_s*(Hgas - u*cs) - scrh*beta_v
                 - alpha_f*a*Btmag/sqrt_rho; 
    #if BACKGROUND_FIELD == YES
    Rc[ENG][k] -= B0y*Rc[BXt][k] + B0z*Rc[BXb][k];
    #endif

    eta[k] =   alpha_s*(X*dV[RHO] + dV[PRS]) - rho*scrh*beta_dv
             - rho*alpha_s*cs*dV[VXn]        - sqrt_rho*alpha_f*a*beta_dB;
    #elif EOS == ISOTHERMAL
    eta[k] =   alpha_s*(0.*X + a2)*dV[RHO] - rho*scrh*beta_dv
             - rho*alpha_s*cs*dV[VXn]      - sqrt_rho*alpha_f*a*beta_dB;
    #endif

    eta[k] *= 0.5/a2;

   /* -----------------------
       Slow wave:  u + c_s
      ----------------------- */

    k = KSLOWP;
    lambda[k] = u + cs; 

    Rc[RHO][k] = alpha_s;
    Rc[MXn][k] = alpha_s*lambda[k]; 
    Rc[MXt][k] = alpha_s*v + scrh*beta_y;
    Rc[MXb][k] = alpha_s*w + scrh*beta_z;

    Rc[BXt][k] = Rc[BXt][KSLOWM];
    Rc[BXb][k] = Rc[BXb][KSLOWM];

    #if EOS == IDEAL
    Rc[ENG][k] =   alpha_s*(Hgas + u*cs) + scrh*beta_v
                 - alpha_f*a*Btmag/sqrt_rho;
    #if BACKGROUND_FIELD == YES
    Rc[ENG][k] -= B0y*Rc[BXt][k] + B0z*Rc[BXb][k];
    #endif

    eta[k] =   alpha_s*(X*dV[RHO] + dV[PRS]) + rho*scrh*beta_dv
             + rho*alpha_s*cs*dV[VXn]        - sqrt_rho*alpha_f*a*beta_dB; 
    #elif EOS == ISOTHERMAL
    eta[k] =   alpha_s*(0.*X + a2)*dV[RHO] + rho*scrh*beta_dv
             + rho*alpha_s*cs*dV[VXn]      - sqrt_rho*alpha_f*a*beta_dB; 
    #endif

    eta[k] *= 0.5/a2;

   /* ------------------------
       Alfven wave:  u - c_a
      ------------------------ */

    k = KALFVM;
    lambda[k] = u - ca;

    Rc[MXt][k] = - rho*beta_z;  
    Rc[MXb][k] = + rho*beta_y;
    Rc[BXt][k] = - sBx*sqrt_rho*beta_z;   
    Rc[BXb][k] =   sBx*sqrt_rho*beta_y;
    #if EOS == IDEAL
    Rc[ENG][k] = - rho*(v*beta_z - w*beta_y);
    #if BACKGROUND_FIELD == YES
    Rc[ENG][k] -= B0y*Rc[BXt][k] + B0z*Rc[BXb][k];
    #endif
    #endif

    eta[k] = + beta_y*dV[VXb]               - beta_z*dV[VXt] 
             + sBx/sqrt_rho*(beta_y*dV[BXb] - beta_z*dV[BXt]);

    eta[k] *= 0.5;

   /* -----------------------
       Alfven wave:  u + c_a 
      ----------------------- */

    k = KALFVP;
    lambda[k] = u + ca;

    Rc[MXt][k] = - Rc[MXt][KALFVM];  
    Rc[MXb][k] = - Rc[MXb][KALFVM];
    Rc[BXt][k] =   Rc[BXt][KALFVM];   
    Rc[BXb][k] =   Rc[BXb][KALFVM];
    #if EOS == IDEAL
    Rc[ENG][k] = - Rc[ENG][KALFVM];
    #if BACKGROUND_FIELD == YES
    Rc[ENG][k] -= B0y*Rc[BXt][k] + B0z*Rc[BXb][k];
    #endif
    #endif

    eta[k] = - beta_y*dV[VXb]               + beta_z*dV[VXt] 
             + sBx/sqrt_rho*(beta_y*dV[BXb] - beta_z*dV[BXt]);

    eta[k] *= 0.5;

   /* -----------------------------------------
      4g. Compute maximum signal velocity
      ----------------------------------------- */

    cmax[i] = fabs (u) + cf;
    g_maxMach = MAX (fabs (u / a), g_maxMach);
    NFLX_LOOP(k) alambda[k] = fabs(lambda[k]);

   /* --------------------------------
      4h. Entropy Fix 
      -------------------------------- */
      
    if (alambda[KFASTM] < 0.5*delta) {
      alambda[KFASTM] = lambda[KFASTM]*lambda[KFASTM]/delta + 0.25*delta;
    }
    if (alambda[KFASTP] < 0.5*delta) {
      alambda[KFASTP] = lambda[KFASTP]*lambda[KFASTP]/delta + 0.25*delta;
    }

    if (alambda[KSLOWM] < 0.5*delta) {
      alambda[KSLOWM] = lambda[KSLOWM]*lambda[KSLOWM]/delta + 0.25*delta;
    }
    if (alambda[KSLOWP] < 0.5*delta) {
      alambda[KSLOWP] = lambda[KSLOWP]*lambda[KSLOWP]/delta + 0.25*delta; 
    }
   
  /*  ---------------------------------
      4i. Compute Roe numerical flux 
      --------------------------------- */

    for (nv = 0; nv < NFLX; nv++) {
      scrh = 0.0;
      for (k = 0; k < NFLX; k++) {
        scrh += alambda[k]*eta[k]*Rc[nv][k];
      }
      sweep->flux[i][nv] = 0.5*(fL[i][nv] + fR[i][nv] - scrh);
    }
    sweep->press[i] = 0.5*(pL[i] + pR[i]);
    
  /* --------------------------------------------------------
     4j. Check the Roe matrix condition, FR - FL = A*(UR - UL)
         where A*(UR - UL) = R*lambda*eta.
    -------------------------------------------------------- */

    #if CHECK_ROE_MATRIX == YES
     for (nv = 0; nv < NFLX; nv++){
       dV[nv] = fR[i][nv] - fL[i][nv]; 
       if (nv == MXn) dV[MXn] += pR[i] - pL[i];
       for (k = 0; k < NFLX; k++){
         dV[nv] -= Rc[nv][k]*eta[k]*lambda[k];
       }
       if (fabs(dV[nv]) > 1.e-4){
         printf (" ! Roe_Solver(): matrix condition not satisfied, var = %d\n", nv);
         printf (" ! Err = %12.6e\n",dV[nv]); 
         ShowVector(vL, NFLX);
         ShowVector(vR, NFLX);
         QUIT_PLUTO(1);
       }
     } 
    #endif

  /* -------------------------------------------------------------
     4k. Save max and min Riemann fan speeds for EMF computation.
     ------------------------------------------------------------- */

    sweep->SL[i]  = lambda[KFASTM];
    sweep->SR[i]  = lambda[KFASTP];
    sweep->SaL[i] = lambda[KALFVM];
    sweep->SaR[i] = lambda[KALFVP];
    #if HAVE_ENERGY
    sweep->Sc[i]  = lambda[KENTRP];
    #endif

  /* -----------------------------------------------------------------
     4l. Hybridize with HLL solver: replace occurences of unphysical 
         states (p < 0, rho < 0) with HLL Flux. Reference:
      
        "A Positive Conservative Method for MHD based based on HLL 
         and Roe methods", P. Janhunen, JCP (2000), 160, 649
     ----------------------------------------------------------------- */

    #if HLL_HYBRIDIZATION == YES
    if (SL[i] < 0.0 && SR[i] > 0.0){
      ifail = 0;    

    /* -----------------------
         check left state
       ----------------------- */

      #if EOS == ISOTHERMAL
      Uv[RHO] = uL[RHO] + (sweep->flux[i][RHO] - fL[i][RHO])/SL[i];        
      ifail  = (Uv[RHO] < 0.0);
      #else
      NFLX_LOOP(nv) {
        Uv[nv] = uL[nv] + (sweep->flux[i][nv] - fL[i][nv])/SL[i];        
      }
      Uv[MXn] += (sweep->press[i] - pL[i])/SL[i];    
 
      vel2  = Uv[MX1]*Uv[MX1] + Uv[MX2]*Uv[MX2] + Uv[MX3]*Uv[MX3];
      b2    = Uv[BX1]*Uv[BX1] + Uv[BX2]*Uv[BX2] + Uv[BX3]*Uv[BX3];    
      pr    = Uv[ENG] - 0.5*vel2/Uv[RHO] - 0.5*b2;
      ifail = (pr < 0.0) || (Uv[RHO] < 0.0);
      #endif

    /* -----------------------
         check right state
       ----------------------- */

      #if EOS == ISOTHERMAL
      Uv[RHO] = uR[RHO] + (sweep->flux[i][RHO] - fR[i][RHO])/SR[i];
      ifail  = (Uv[RHO] < 0.0);
      #else
      for (nv = NFLX; nv--;  ){
        Uv[nv] = uR[nv] + (sweep->flux[i][nv] - fR[i][nv])/SR[i];
      }
      Uv[MXn] += (sweep->press[i] - pR[i])/SR[i];

      vel2  = Uv[MX1]*Uv[MX1] + Uv[MX2]*Uv[MX2] + Uv[MX3]*Uv[MX3];
      b2    = Uv[BX1]*Uv[BX1] + Uv[BX2]*Uv[BX2] + Uv[BX3]*Uv[BX3];
      pr    = Uv[ENG] - 0.5*vel2/Uv[RHO] - 0.5*b2;
      ifail = (pr < 0.0) || (Uv[RHO] < 0.0);
      #endif

    /* -------------------------------------------------------
        Use the HLL flux function if the interface lies 
        within a strong shock. The effect of this switch is 
        visible in the Mach reflection test.
       ------------------------------------------------------- */

      #if DIMENSIONS > 1   
      #if EOS == ISOTHERMAL
      scrh  = fabs(vL[RHO] - vR[RHO]);
      scrh /= MIN(vL[RHO], vR[RHO]);
      #else       
      scrh  = fabs(vL[PRS] - vR[PRS]);
      scrh /= MIN(vL[PRS], vR[PRS]);
      #endif
      if (scrh > 1.0 && (vR[VXn] < vL[VXn])) ifail = 1;
      #endif
      
      if (ifail){
        scrh = 1.0/(SR[i] - SL[i]);
        for (nv = 0; nv < NFLX; nv++) {
          sweep->flux[i][nv] = SL[i]*SR[i]*(uR[nv] - uL[nv]) +
                               SR[i]*fL[i][nv] - SL[i]*fR[i][nv];
          sweep->flux[i][nv] *= scrh;
        }
        sweep->press[i] = (SR[i]*pL[i] - SL[i]*pR[i])*scrh;
      }
    }
    #endif  /* HLL_HYBRIDIZATION == YES */
  }

/* --------------------------------------------------------
   5. Define point and diffusive fluxes for CT
   -------------------------------------------------------- */
  
#if DIVB_CONTROL == CONSTRAINED_TRANSPORT 
  CT_Flux (sweep, beg, end, grid);
#endif

/* --------------------------------------------------------
   6. Compute Powell's source term
   -------------------------------------------------------- */
  
  #if DIVB_CONTROL == EIGHT_WAVES
  Roe_DivBSource (sweep, beg + 1, end, grid);
  #endif

/* --------------------------------------------------------
   7. Add CR flux contribution using simplified
      upwinding.
   -------------------------------------------------------- */

  #if (PARTICLES == PARTICLES_CR) && (PARTICLES_CR_FEEDBACK == YES) 
  Particles_CR_Flux (stateL, beg, end);
  Particles_CR_Flux (stateR, beg, end);

  for (i = beg; i <= end; i++){
    for (nv = NFLX; nv--; ) {
      sweep->flux[i][nv] += 0.5*(stateL->fluxCR[i][nv] + stateR->fluxCR[i][nv]);
//                          -2.0*(stateR->u[i][nv] - stateL->u[i][nv]);
    }
  }  
  #endif

}
#undef sqrt_1_2
#undef HLL_HYBRIDIZATION
#undef CHECK_ROE_MATRIX
