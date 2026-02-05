/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief This file (`init.c`) contains user-supplied functions for
         problem configuration and initialization of the Line-Driven
         Disk Wind simulation.

  \details The provided functions are used to set up and configure
           the initial conditions, boundary conditions, and other
           problem-specific parameters.

  \author Nick H
  \date   December 3, 2018

  \revised_by A. Mosallanezhad (a.mosallanezhad@soton.ac.uk)
  \date       December 17, 2025
*/
/* ///////////////////////////////////////////////////////////////////// */


#include "pluto.h"

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *********************************************************************** */
{

	 double cent_mass, rho_0, rho_alpha;
	 double disk_mdot, mu, dfloor;
 	 double rmin, gm_cgs, gm_code;
 	 double temp, cs, rho;
 	 double T_iso, rho_d, rho_a;
	 double r_0, teff;
	 double s, c, R, r, z;
	 double v_k, R_1, tan_1;
	 double H, Bphi;


   /* constant parameters from pluto.ini */

 	 cent_mass = g_inputParam[CENT_MASS];    /* central mass */
	 disk_mdot = g_inputParam[DISK_MDOT];    /* disk accretion rate */
	 rho_0     = g_inputParam[RHO_0];        /* density at r_0 set to R_ISCO */
 	 rho_alpha = g_inputParam[RHO_ALPHA];    /* drop off exponent of density */
   T_iso     = g_inputParam[T_ISO];        /* isothermal temperature  */
 	 mu        = g_inputParam[MU];           /* mean molecular weight */
   dfloor    = g_inputParam[DFLOOR];       /* minimum density */


	 rmin    = g_domBeg[IDIR];
   r_0     = rmin * UNIT_LENGTH;           /* r_0 - in cgs unit */

	 gm_cgs  = CONST_G * cent_mass;
	 gm_code = gm_cgs / (UNIT_LENGTH * UNIT_VELOCITY * UNIT_VELOCITY);


	 r   = x1;
   s   = sin(x2);
   c   = cos(x2);

   R  = r * s;   /* Cylindrical radius */
   z  = r * c;


	 R_1   = 1.0 / R;
	 tan_1 = 1.0 / tan(x2);


	  if (fabs(tan(x2))< 1.e-12) {
	 	 tan_1 = 0.0;
		 R_1   = 0.0;
	  }


	 teff  = pow(3.0 * gm_cgs * disk_mdot / (8.0 * CONST_PI * CONST_sigma), 0.25);
	 teff *= pow(r_0, -0.75);
//
	 temp  = teff * pow(rmin / r, 0.75);

	 cs    = sqrt(CONST_Rgas * temp / mu) / UNIT_VELOCITY;
	 v_k   = sqrt(gm_code / r);

	 rho_a = dfloor  / UNIT_DENSITY;

	 /* calculate density in code units - the expression from PSD98 */
	 rho_0 = (rho_0 / UNIT_DENSITY) * pow(rmin / r, rho_alpha);
	 rho_d = rho_0  * exp(-1.0 * gm_code * tan_1 * tan_1 / (2.0 * cs * cs * r));



   us[iVR]   = 0.0;
	 us[iVTH]  = 0.0;
	 us[iVPHI] = 0.0;


	 if (rho_d > rho_a)  {
 	 	 us[RHO] = rho_d;
     us[TRC] = 1.0;
   } else {
     us[RHO] = rho_a;
     us[TRC] = 0.0;
   }


   us[iVPHI] = v_k * s;

 	 /* this converts temperature (in K) to pressure in code units */
 	 #if HAVE_ENERGY
 		 	us[PRS]  = us[RHO] * temp / (KELVIN * mu);
			g_gamma = g_inputParam[GAMMA];
 	 #endif


	#if EOS == ISOTHERMAL
		temp = g_inputParam[T_ISO];
	 	g_isoSoundSpeed = sqrt(CONST_Rgas * temp / mu) / UNIT_VELOCITY;
	#endif


   #if LINE_DRIVEN_WIND != NO
     dvds_setup_flag = 0;
   #endif


}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*!
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*!
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures
 *
 *********************************************************************** */
{

}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif


/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
/*!
 *  Assign user-defined boundary conditions.
 *
 *********************************************************************** */
{

	int   i, j, k, nv;

	double *r, *th, *phi, *x2_glob, *x2r;
	double dfloor, tfloor, pfloor;
	double rho_old, rho_new, vx2_old, vx2_new;
	double v[256];

	double x2max = g_domEnd[JDIR];

	r   = grid->xgc[IDIR];
	th  = grid->xgc[JDIR];
	phi = grid->xgc[KDIR];

	x2r = grid->xr[JDIR];

	x2_glob    = grid->x_glob[JDIR];

	dfloor     = g_inputParam[DFLOOR] / UNIT_DENSITY;
	tfloor     = 5.e2;
	pfloor     = dfloor * tfloor / (KELVIN * g_inputParam[MU]);



	if (side == 0){

		/* avoid too small density near the boundary */
		RBox dom_box;
		int convert_to_cons;

		TOT_LOOP(k,j,i){

			   convert_to_cons = 0;


				 /* this should be the last 'real' theta bin - before the ghost zones. */

					if (x2r[j] >= x2max)  {

						rho_old = d->Vc[RHO][k][j][i];
						vx2_old = d->Vc[VX2][k][j][i];

						Init (v, r[i], th[j], phi[k]);
					  for (nv = 0; nv < NVAR; nv++) d->Vc[nv][k][j][i]= v[nv];

						rho_new = d->Vc[RHO][k][j][i];
						vx2_new = rho_old * vx2_old / rho_new;
						d->Vc[VX2][k][j][i] = vx2_new;

				 }


				 // Check if density is negative or below the floor

				 if (d->Vc[RHO][k][j][i] < dfloor) {

						 double rho_new = (d->Vc[RHO][k][j][i] > 0) ? d->Vc[RHO][k][j][i] : 0.0;
						 if (rho_new < dfloor) {
								 rho_new = dfloor;  // Set density to floor if it’s negative or too low
						 }

						 d->Vc[RHO][k][j][i] = dfloor;

						 // #if HAVE_ENERGY
							// 	 if (d->Vc[PRS][k][j][i] < pfloor) {
							// 			 d->Vc[PRS][k][j][i] = pfloor;
							// 	}
						 // #endif

						 convert_to_cons = 1;
				 }

				 // #if HAVE_ENERGY
					// 	 if (d->Vc[PRS][k][j][i] < pfloor) {
					// 			 d->Vc[PRS][k][j][i] = pfloor;
					// 			convert_to_cons = 1;
					// 	}
				 // #endif



				 if (convert_to_cons) {
				 		RBoxDefine (i, i, j, j, k, k, CENTER, &dom_box);
				 		PrimToCons3D(d->Vc, d->Uc, &dom_box, grid);
				 }


		 	} /* DOM_LOOP() */
	  } /* if (side == 0) */



		if (side == X1_BEG){                      /* -- X1_BEG boundary -- */
		    if (box->vpos == CENTER){
		      BOX_LOOP(box,k,j,i){
						 for (nv = 0; nv < NVAR; nv++){
						    d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IBEG];
						  }

		        if (d->Vc[iVR][k][j][i] > 0.0) d->Vc[iVR][k][j][i] = 0.0;
		      }
		    }else if (box->vpos == X2FACE){
		      #ifdef STAGGERED_MHD
		       BOX_LOOP(box,k,j,i) d->Vs[BX2s][k][j][i] = d->Vs[BX2s][k][j][IBEG];
		      #endif
		    }else if (box->vpos == X3FACE){
		      // #ifdef STAGGERED_MHD
		      //  BOX_LOOP(box,k,j,i) d->Vs[BX3s][k][j][i] = d->Vs[BX3s][k][j][IBEG];
		      // #endif
		    }
		  }


			if (side == X1_END){  /* -- X1_END boundary -- */
				if (box->vpos == CENTER){
					BOX_LOOP(box, k, j, i){
						for (nv = 0; nv < NVAR; nv++) {
							d->Vc[nv][k][j][i]= d->Vc[nv][k][j][IEND];
						}

						if (d->Vc[iVR][k][j][i] < 0.0) d->Vc[iVR][k][j][i] = 0.0;

					}
				}else if (box->vpos == X2FACE){
					#ifdef STAGGERED_MHD
					 BOX_LOOP(box,k,j,i)  d->Vs[BX2s][k][j][i] = d->Vs[BX2s][k][j][IEND];
					#endif
				}else if (box->vpos == X3FACE){
					#ifdef STAGGERED_MHD
					 // BOX_LOOP(box,k,j,i)  d->Vs[BX3s][k][j][i] = d->Vs[BX3s][k][j][IEND];
					#endif
				}
			}


			if (side == X2_BEG) {          /* -- X2_BEG boundary: rotation axis -- */
			  if (box->vpos == CENTER) {
			    BOX_LOOP(box,k,j,i) {

				      /* 1) Mirror all primitive variables from the interior */
				      for (nv = 0; nv < NVAR; nv++) {
				        d->Vc[nv][k][j][i] = d->Vc[nv][k][2*JBEG - j - 1][i];
				      }

				      /* 2) Enforce correct parity for vector components at polar axis */
				      /* Vr (VX1) stays as mirrored (even). */
				      d->Vc[VX2][k][j][i] *= -1.0;  /* Vθ odd */
				      d->Vc[VX3][k][j][i] *= -1.0;  /* Vφ odd, always flipped now */

				      /* 3) “Outflow-like” BC for scalars: zero θ-gradient at the axis */
				      d->Vc[RHO][k][j][i] = d->Vc[RHO][k][JBEG][i];

							#if EOS != ISOTHERMAL
							      d->Vc[PRS][k][j][i] = d->Vc[PRS][k][JBEG][i];
							#endif
			    }
			 }
	  }


}


#if BODY_FORCE != NO
/* ********************************************************************* */
		void BodyForceVector(double *v, double *g, double x1, double x2, double x3)

		{

			double cent_mass, gm_cgs, gm_code;

			cent_mass = g_inputParam[CENT_MASS];
			gm_cgs    = CONST_G * cent_mass;
			gm_code   = gm_cgs / (UNIT_LENGTH * UNIT_VELOCITY * UNIT_VELOCITY);


			g[IDIR] = - 1.0 * gm_code / (x1 * x1);
			g[JDIR] = 0.0;
			g[KDIR] = 0.0;

		}

		/* ********************************************************************* */
		double BodyForcePotential(double x1, double x2, double x3)
		/*!
		 * Return the gravitational potential as function of the coordinates.
		 *
		 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
		 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
		 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
		 *
		 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
		 *
		 *********************************************************************** */
		{
		}


#endif
