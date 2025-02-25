/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute the right hand side of the relativistic 
         hydro (RHD) equations in primitive form.
  
  Implements the right hand side of the quasi-linear form of the 
  relativistic hydro equations. 
  In 1D this may be written as
  \f[ 
      \partial_t{\mathbf{V}} = - A\cdot\partial_x\mathbf{V} + \mathbf{S}
  \f]
  where \f$ A \f$ is the matrix of the primitive form of the equations,
  \f$ S \f$ is the source term.

  \b Reference:
    - "The Piecewise Parabolic Method for  Multidimensional Relativistic 
       Fluid Dynamics", Mignone, Plewa and Bodo, ApJS (2005) 160,199.
 
  The function PrimRHS() implements the first term while PrimSource() 
  implements the source term part.

  \authors A. Mignone (andrea.mignone@unito.it) \n
	         L. DelZanna
  \date   Apr 05, 2024
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void PrimRHS (double *w, double *dw, double cs2, double h, double *Adw)
/*!
 * Compute the matrix-vector multiplication \f$ A(\mathbf{v})\cdot 
 * \Delta\mathbf{v} \f$ where \c A is the matrix of the quasi-linear form 
 * of the RHD equations.
 *
 * \param [in]  w    vector of primitive variables;
 * \param [in]  dw   limited (linear) slopes;
 * \param [in]  cs2  local sound speed;
 * \param [in]  h    local enthalpy;
 * \param [out] Adw  matrix-vector product.
 *
 * \note 
 * \return This function has no return value.
 ************************************************************************** */
{
  int nv;
  double rho, u1, u2, u3;
  double d_2, g2, scrh;
  double g, v1, v2, v3, gx2;

#if RECONSTRUCT_4VEL
  rho = w[RHO];
  u1  = w[VXn];
  u2  = w[VXt];
  u3  = w[VXb];

  scrh = u1*u1 + u2*u2 + u3*u3;
  g2   = 1.0 + scrh;
  g    = sqrt(g2);
  d_2  = g/(g2 - scrh*cs2);

/* Get 3vel  */

  v1 = u1/g;
  v2 = u2/g;
  v3 = u3/g;

  gx2 = 1.0/(1.0 - v1*v1);

  scrh = v1*dw[VXn] + v2*dw[VXt] + v3*dw[VXb];

  Adw[PRS] =  d_2*(rho*h*cs2*(dw[VXn] - v1*scrh)
                      + u1*(1.0 - cs2)*dw[PRS]);

  Adw[RHO] = v1*dw[RHO] - (v1*dw[PRS] - Adw[PRS])/(h*cs2);
 
  scrh = 1.0/(g*rho*h);
  d_2  = u1*dw[PRS] - g*Adw[PRS];

  Adw[VXn] = v1*dw[VXn] + scrh*(dw[PRS] + u1*d_2);
  Adw[VXt] = v1*dw[VXt] + scrh*u2*d_2;
  Adw[VXb] = v1*dw[VXb] + scrh*u3*d_2;
  
#else
  rho = w[RHO];
  v1  = w[VXn];
  v2  = w[VXt];
  v3  = w[VXb];

  g2  = v1*v1 + v2*v2 + v3*v3;
  d_2 = 1.0/(1.0 - g2*cs2);
  g2  = 1.0/(1.0 - g2);

  Adw[PRS] = d_2*(cs2*rho*h*dw[VXn] 
              + v1*(1.0 - cs2)*dw[PRS]);

  Adw[RHO] = v1*dw[RHO] - (v1*dw[PRS] - Adw[PRS])/(h*cs2);
  
  scrh = 1.0/(g2*rho*h);
  Adw[VXn] =  v1*dw[VXn] + scrh*(dw[PRS] - v1*Adw[PRS]);
  Adw[VXt] =  v1*dw[VXt] - scrh*v2*Adw[PRS];
  Adw[VXb] =  v1*dw[VXb] - scrh*v3*Adw[PRS];
#endif

#if NSCL > 0 
  NSCL_LOOP(nv)  Adw[nv] = v1*dw[nv];
#endif

}

/* *********************************************************************  */
void PrimSource (const State *state, double **src, int beg, int end, Grid *grid)
/*!
 * Compute source terms of the RHD equations in primitive variables.
 *
 *  - Geometrical sources;
 *  - Gravity;
 *
 *  The rationale for choosing during which sweep a particular source 
 *  term has to be incorporated should match the same criterion used 
 *  during the conservative update. 
 *  For instance, in polar or cylindrical coordinates, curvilinear source
 *  terms are included during the radial sweep only.
 * 
 * \param [in]  state   pointer to a Sweep structure
 * \param [out] src     array of source terms
 * \param [in]  beg     initial index of computation
 * \param [in]  end     final   index of computation
 * \param [in]  grid    pointer to a Grid structure
 *
 *********************************************************************** */
{
  int    nv, i;
  double r_1, scrh, alpha;
  double vel2, delta;
  double lor2, lor2cs, vg;
  double *v, *x1, *x2, *x3, g[3];
  double *a2 = state->a2;
  double *h  = state->h;

#if GEOMETRY == CYLINDRICAL
  x1 = grid->xgc[IDIR];
  x2 = grid->xgc[JDIR];
  x3 = grid->xgc[KDIR];
#else  
  x1 = grid->x[IDIR]; 
  x2 = grid->x[JDIR]; 
  x3 = grid->x[KDIR]; 
#endif

  for (i = beg; i <= end; i++){
  for (nv = NVAR; nv--;  ){
    src[i][nv] = 0.0;
  }}

#if GEOMETRY == CARTESIAN

#elif GEOMETRY == CYLINDRICAL 
  
  if (g_dir == IDIR){
    for (i = beg; i <= end; i++) {
 
      r_1  = 1.0/x1[i];
      v    = state->v[i];
      vel2 = v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3];

      #if RECONSTRUCT_4VEL
      scrh    = sqrt(1.0 + vel2);
      alpha   = v[VXn]*r_1*scrh/(1.0 + vel2*(1.0 - a2[i]));
      scrh    = a2[i]*alpha;
      printLog ("! Primitive source terms not yet implemented\n");
      printLog ("! with 4-vel. Please try 3-vel\n");
      QUIT_PLUTO(1);
      #else
      alpha = v[VXn]*r_1/(1.0 - a2[i]*vel2);
      scrh  = a2[i]*(1.0 - vel2)*alpha;
      #endif

      src[i][RHO] = -v[RHO]*alpha;
      src[i][VX1] = scrh*v[VX1];
      src[i][VX2] = scrh*v[VX2];
      src[i][VX3] = scrh*v[VX3];

      src[i][iVR]   +=  v[iVPHI]*v[iVPHI]*r_1;
      src[i][iVPHI] += -v[iVPHI]*v[iVR]*r_1;
      
      src[i][PRS] = -a2[i]*v[RHO]*h[i]*alpha;

    }
  }

#else 

  printLog ("! PrimRHS(): primitive source terms not available for this geometry\n");
  printLog ("!            Use RK integrators\n");
  QUIT_PLUTO(1);

#endif   
  
/* -----------------------------------------------------------
                   Add body force
   ----------------------------------------------------------- */
  
   #if (BODY_FORCE & POTENTIAL)
   #error Cannot use BodyForcePotential in relativistic flows
   #endif

   #if (BODY_FORCE & VECTOR)
   printLog ("! PrimRHS(): Cannot use BodyForcePotential in relativistic flows\n");
   QUIT_PLUTO(1);
   #endif
}
