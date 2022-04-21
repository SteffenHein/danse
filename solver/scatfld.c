/* [ file: scatfld.c ] */
/*******************************************************************************
*                                                                              *
*   Function scatfld(*)                                                        *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r3 ]                                                          *
*                                                                              *
*   DSC scattering algorithm based on non-orthogonal condensed node            *
*   with tensor material parameters and gyrotropic structure.                  *
*   [ cf. 'Finite Difference Time-Domain Approximation of Maxwell's            *
*   Equations with Non-Orthogonal Condensed TLM Mesh', International           *
*   Journal of Numerical Modelling, vol.7, 179-188 (1994) and 'TLM             *
*   Numerical Solution of Bloch's Equations for Magnetized Gyrotropic          *
*   Media', Applied Mathematical Modelling, vol.22, 221-229 (1997).]           *
*                                                                              *
*   (C) SHEIN; Bad Aibling, October 2007                   Steffen Hein        *
*   [ Update: April 05, 2022 ]                          <contact@sfenx.de>     *
*                                                                              *
*******************************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
/*----------------------------------------------------------------------------*/
# include "../math/maths.h"  /* 'my' computation environment headers */
/*----------------------------------------------------------------------------*/
# if defined ( OPTIMIZE )
   # pragma OPT_LEVEL 3
   # pragma OPTIMIZE ON
# endif
/*----------------------------------------------------------------------------*/
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
/* Edit and customize this general configuration header: */
# include "../CONFIG.H" 
/*----------------------------------------------------------------------------*/
/* Edit and customize solver configuration header SOLVER.CONF */
# include "SOLVER.CONF"
/*----------------------------------------------------------------------------*/
# include "solvtp.h"
/*----------------------------------------------------------------------------*/
/* configuration header */
# include "SCATFLD.CONF"
/*============================================================================*/
/* function body */
# include "scatfld-1.0r3.a"
/*============================================================================*/
/************* end of Maxwell field updating function scatfld(*) **************/
/************************* end of file scatfld.c ******************************/
