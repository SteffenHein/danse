/* [ file: connct.c ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function cnctfld(*)                                                 *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r3 ]                                                          *
*                                                                              *
*   DSC Maxwell field connection map                                           *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 18, 2022 ]                             <contact@sfenx.de>  *
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
/* Edit and customize solver configuration header SOLVER.CONF */
# include "SOLVER.CONF"
/*----------------------------------------------------------------------------*/
# include "solvtp.h"
/*----------------------------------------------------------------------------*/
/* the configuration header */
# include "CNCTFLD.CONF"
/*============================================================================*/
/* the function body */
# include "cnctfld.h"
/*============================================================================*/
# undef CNN_CONJGT
/************ end of Maxwell field connection function cnctfld(*) *************/





/*******************************************************************************
*                                                                              *
*   ANSI C function cncthcr(*)                                                 *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r3 ]                                                          *
*                                                                              *
*   DSC diffusive heat currents connection map                                 *
*   [ convective heat currents and fluid flows are separately connected in     *
*     function conflow(*) ]                                                    *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 18, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# if DSC_HCRMDE != 0
/*----------------------------------------------------------------------------*/
/* the configuration header */
# include "CNCTHCR.CONF"
/*============================================================================*/
/* the function body */
# include "cncthcr.h"
/*============================================================================*/
/*************** end of thermal connection function cncthcr(*) ****************/
# endif /* DSC_HCRMDE != 0 */
/******************************************************************************/





/*******************************************************************************
*                                                                              *
*   ANSI C function conflow(*)                                                 *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r3 ]                                                          *
*                                                                              *
*   DSC fluid flow connection map                                              *
*   [ diffusive heat currents are separately connected in function             *
*     cncthcr(*) ]                                                             *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 18, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# if DSC_HCRMDE != 0
# if DSC_FLDMDE != 0
/*----------------------------------------------------------------------------*/
/* the configuration header */
# include "CONFLOW.CONF"
/*============================================================================*/
/* the function body: */
# include "conflow.h"
/*============================================================================*/
/************** end of fluid flow connection function conflow(*) **************/
/******************************************************************************/





/*******************************************************************************
*                                                                              *
*   Function lesfilter(*)                                                      *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r3 ]                                                          *
*                                                                              *
*   A simple Large Eddy Simulation filter                                      *
*   [ based on cellular coarse graining ]                                      *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 18, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/

/*============================================================================*/
/* the function body: */
# include "lesflt.h"
/*============================================================================*/
/****************** end of LES filter function lesfilter(*) *******************/
# endif /* DSC_FLDMDE != 0 */
# endif /* DSC_HCRMDE != 0 */
/******************************************************************************/
# if defined ( OPTIMIZE )
   # pragma OPT_LEVEL 3
   # pragma OPTIMIZE OFF
# endif
/************************** end of file connct.c ******************************/
