/* [ file: maths.h ] */
/*******************************************************************************
*                                                                              *
*   User defined computing environment header maths.h                          *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# ifndef __maths_h /* this condition embraces the whole content of this file */

# define __maths_h
/*----------------------------------------------------------------------------*/
# if _ISOC99_SOURCE == 1
   # include <fenv.h> 
   # include <iso646.h>
# endif
/*----------------------------------------------------------------------------*/
/*
# if defined( _GNU_Linux )
   # include <features.h>
# elif defined( _SuSE )
   # include <features.h>
# endif
*/
/*----------------------------------------------------------------------------*/
# include <math.h>
# include <float.h>
# include <limits.h>
/*----------------------------------------------------------------------------*/
# endif /* ifndef __maths_h */
/*************************** end of file maths.h *****************************/
