/* [ file: idealgas.c ] */
/*******************************************************************************
*                                                                              *
*   Function idealgas(*)                                                       *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r1 ]                                                          *
*                                                                              *
*   Returns nodal pressure as a function of temperature and density            *
*   [ Boyle-Mariotte law of ideal gas ]                                        *
*                                                                              *
*   (C) SHEIN; Bad Aibling, October 2007                   Steffen Hein        *
*   [ Update: April 05, 2022 ]                          <contact@sfenx.de>     *
*                                                                              *
*******************************************************************************/
# if (( CMPRSSBL != 0 )\
    &&( BSQAPRX == 0 ))
/*----------------------------------------------------------------------------*/
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
/*============================================================================*/
# include "idealgas.b"
/*============================================================================*/
# endif /* (( CMPRSSBL != 0 )\
          &&( BSQAPRX == 0 )) */
/************************ end of function idealgas(*) *************************/
/************************** end of file idealgas.c ****************************/
