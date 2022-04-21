/* [ file: cnctfld.c ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function cnctfld(*)                                                 *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r3 ]                                                          *
*                                                                              *
*   DSC Maxwell field connection map                                           *
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
/* cofiguration header */
# include "CNCTFLD.CONF"
/*============================================================================*/
/* function body */
# include "cnctfld.h"
/*============================================================================*/
/************ end of Maxwell field connection function cnctfld(*) *************/
/************************** end of file cnctfld.c *****************************/
