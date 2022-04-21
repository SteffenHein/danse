/* [ file: cllsts.c ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function cllsts(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Given the new incident and stored reflected port quantities of a DSC cell  *
*   this function computes physical fields, charge, mean particle velocity,    *
*   and relativistic mass within the cell.                                     *
*   All quantities are transferred to and returned from the function via the   *
*   structure cst of type CLLSTS.                                              *
*                                                                              *
*   (C) SHEIN; Bad Aibling, October 2007                   Steffen Hein        *
*   [ Update: April 05, 2022 ]                          <contact@sfenx.de>     *
*                                                                              *
*******************************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <stddef.h>
# include <string.h>
/*----------------------------------------------------------------------------*/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
# include "../CONFIG.H"
# include "SOLVER.CONF"
/*----------------------------------------------------------------------------*/
# include "solvtp.h"
/*----------------------------------------------------------------------------*/
static CLLSTS cst = {null};
/*----------------------------------------------------------------------------*/

/*============================================================================*/

CLLSTS *\
cllsts( CLLSTS *cpt )
{
   static CLLSTS *csp = &cst;

   ( csp->rtn ) = null;
   return csp;
} 
/*============================================================================*/
