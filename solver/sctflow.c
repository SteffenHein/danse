/* [ file: sctflow.c ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function sctflow(*)                                                 *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r3 ]                                                          *
*                                                                              *
*   DSC nodal fluid flow updating function                                     *
*   [ Stable version, diffusive temperature updates in function scathcr(*) ]   *
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
# if DSC_FLDMDE != 0
/*----------------------------------------------------------------------------*/
# include "solvtp.h"
/*----------------------------------------------------------------------------*/
/* configuration header */
# include "SCTFLOW.CONF"
/*============================================================================*/
/* function body */
# include "sctflow-1.0r3.a"
/*============================================================================*/
/************** end of fluid flow updating function sctflow(*) ****************/
# endif /* DSC_FLDMDE != 0 */
# endif /* DSC_HCRMDE != 0 */
/******************************************************************************/
