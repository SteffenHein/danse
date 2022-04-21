/* [ file: scathcr.c ] */
/*******************************************************************************
*                                                                              *
*   Function scathcr(*)                                                        *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r3 ]                                                          *
*                                                                              *
*   DSC nodal temperatures updating function                                   *
*   [ Stable version separately updating fluid flows ]                         *
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
# if DSC_HCRMDE != 0
/*----------------------------------------------------------------------------*/
# include "solvtp.h"
/*----------------------------------------------------------------------------*/
/* configuration header: */
# include "SCATHCR.CONF"
/*============================================================================*/
/* function body: */
# include "scathcr-1.0r3.a"
/*============================================================================*/
/***************** end of thermal updating function scathcr(*) ****************/
# endif /* DSC_HCRMDE != 0 */
/*************************** end of file scathcr.c ****************************/
