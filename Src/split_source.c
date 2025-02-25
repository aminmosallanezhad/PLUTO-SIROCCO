/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Include source terms using operator splitting.

  The SplitSource() function handles source terms in a separate
  step using operator splitting.
  It is called from Integrate() between standard hydro advances.
  At present these source terms are one or more of the following:

  - optically thin radiative losses (cooling)
  - Diffusion operators:
    - resistivity
    - Thermal conduction
    - Viscosity
  - additional user-defined terms may also be included here.

  It handles also the time integration of radiation fields in the
  nonrelativistic radiation module

  \authors A. Mignone (andrea.mignone@unito.it)
  \date    Oct 26, 2016
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void SplitSource (Data *d, double dt, timeStep *Dts, Grid *grid)
/*!
 *  Take one step on operator-split source terms.
 *
 *  \param [in,out] d   pointer to PLUTO Data structure containing
 *                      the solution array updated from the most
 *                      recent call
 *  \param[in]      dt  the time step used to integrate the source
 *                      terms
 *  \param[in]     Dts  pointer to the time step structure
 *  \param[in]    grid  pointer to an array of grid structures
 *
 *********************************************************************** */
{

/*  ---------------------------------------------
             Cooling/Heating losses
    ---------------------------------------------  */

#if COOLING != NO
  #if COOLING == POWER_LAW  /* -- solve exactly -- */
  PowerLawCooling (d->Vc, dt, Dts, grid);
  #elif COOLING == KROME /* -- Interfaced krome solvers -- */
  KromeCooling (d, dt, Dts, grid);
  #elif COOLING == BLONDIN
  BlondinCooling (d->Vc, d, dt, Dts, grid);
  #else
  CoolingSource (d, dt, Dts, grid);
  #endif
#endif

/* ----------------------------------------------
    Parabolic terms using STS:

    - resistivity
    - thermal conduction
    - viscosity
   ---------------------------------------------- */

#if (PARABOLIC_FLUX & SUPER_TIME_STEPPING)
  STS (d, dt, Dts, grid);
#endif

#if (PARABOLIC_FLUX & RK_LEGENDRE)
  RKL (d, dt, Dts, grid);
#endif

#if RADIATION_NR
  RadSubstepping (d, dt, Dts, grid);
#endif

}
