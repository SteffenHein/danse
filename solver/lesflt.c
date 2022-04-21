/* [ file: lesflt.c ] */
/*******************************************************************************
*                                                                              *
*   Function lesfilter(*)                                                      *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r1 ]                                                          *
*                                                                              *
*   A simple Large Eddy Simulation filter                                      *
*   [ based on cellular coarse graining ]                                      *
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
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
# if defined ( OPTIMIZE )
   # pragma OPT_LEVEL 3
   # pragma OPTIMIZE ON
# endif
/*----------------------------------------------------------------------------*/
/* Edit and customize this general configuration header: */
# include "../CONFIG.H" 
/*----------------------------------------------------------------------------*/
/* Edit and customize elfe configuration header SOLVER.CONF */
# include "SOLVER.CONF"
/*----------------------------------------------------------------------------*/
# include "solvtp.h"
/*----------------------------------------------------------------------------*/
/* configuration header */
# include "CONFLOW.CONF"
/*----------------------------------------------------------------------------*/
/* [Don't] enable mass correction [0] 1 */
/* [ default: 1 ] */
# define LES_MASSCRR 1
/*----------------------------------------------------------------------------*/
/* [Don't] enable density coarsening [0] 1 */
/* [ default: 1 ] */
# define LES_CRSNDNS 1
# if LES_CRSNDNS == 1 /* coarsening requires mass correction !!! */
   # undef LES_MASSCRR
   # define LES_MASSCRR 1
# endif /* LES_CRSNDNS == 1 */
/*============================================================================*/
# include "lesflt.h"
/*============================================================================*/
/************************ end of function lesfilter(*) ************************/
/*************************** end of file lesflt.c *****************************/
