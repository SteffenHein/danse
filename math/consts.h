/* [ file: consts.h ] */
/*******************************************************************************
*                                                                              *
*   Additional constants definition header consts.h                            *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# ifndef __consts_h /* this condition embraces the whole content of this file */

# define __consts_h
/*----------------------------------------------------------------------------*/
# ifndef HUGE_VALF
   # define HUGE_VALF (( double )( 1.000e+277 ))
# endif 
/*----------------------------------------------------------------------------*/
/* relative precision: */
# ifndef PRECISION
   # define PRECISION (( double)( 1.000e-15 ))
# endif
/*----------------------------------------------------------------------------*/
/*
*//* numeric constants:
*/
# undef NULL
# define NULL ( void * ) 0
# define null    0
# define ONE     1
# define TWO     2
# define THREE   3
# define FOUR    4
# define FIVE    5
# define SIX     6
# define SEVEN   7
# define EIGHT   8
# define NINE    9  
# define TEN    10
# define DEC    10
# define ELEVEN 11
# define TWELVE 12

# define HALF  ( 0.5000000000000000000000 )
# define THIRD ( 0.3333333333333333333333 )
# define QUART ( 0.2500000000000000000000 )
# define ZERO  ( 0.0000000000000000000000 )

# undef  PI
# define PI    ( 3.1415926535897932384626 )
/* 
*//* structure types:
*/
typedef struct
{
   double r, i;

   double arg, nrm;
}  COMPLEX;
/* 
*//* string sizes: 
*/
# define VSS_SIZE   10  /* very short string size [ number of characters ] */
# define SHS_SIZE   20  /* short string size      [ number of characters ] */
# define STS_SIZE   80  /* standard string size   [ number of characters ] */
# define LGS_SIZE  400  /* long string size       [ number of characters ] */
# define VLS_SIZE 1200  /* very long string size  [ number of characters ] */
# define VSB_SIZE   32  /* very small block size  [ bytes ] */
# define SMB_SIZE  128  /* small block size       [ bytes ] */
# define STB_SIZE  512  /* standard block size    [ bytes ] */
# define LGB_SIZE 2048  /* large block size       [ bytes ] */
/*
*//* character strings
*//* and standard string initializers [ 80 chars ]:
*/
# define CLEAR_LINE "\r%*s", 78, "                                          "

/* very short string initializer [ 10+1 chars ]: */
# define VSSSTR  "********** "

/* short string initializer [ 20+1 chars ]: */
# define SHSSTR  "******************** "

/* standard str. initializer [ 80+1 chars ]: */
# define STDSTR  "****************************************"\
                 "**************************************** "

/* long strg. initializer [ 400+1 chars ]: */

# define LNGSTR  "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "**************************************** "

/* very long string initializer [ 1200+1 chars ]: */

# define VLGSTR  "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "****************************************"\
                 "**************************************** "

# endif /* ifndef __consts_h */
/*************************** end of file consts.h *****************************/
