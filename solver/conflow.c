/* [ file: conflow.c ] */
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
# if DSC_HCRMDE != 0
# if DSC_FLDMDE != 0
/*----------------------------------------------------------------------------*/
# include "solvtp.h"
/*----------------------------------------------------------------------------*/
/* configuration header */
# include "CONFLOW.CONF"
/*============================================================================*/
/* function body */
# include "conflow.h"
/*============================================================================*/
/************** end of fluid flow connection function conflow(*) **************/
# endif /* DSC_FLDMDE != 0 */
# endif /* DSC_HCRMDE != 0 */
/************************** end of file conflow.c *****************************/
