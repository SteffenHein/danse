/* [ file: dgtize.h ] */
# define DO_DGTIZE "dgtize(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function dgtize(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Transforms long integer decimal expression into g-adic expression          *
*   ( in *option = 'd'igitize ) and vice-versa ( in *option = 'r'estore ).     *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# ifndef DGT_MXORDR
   # define DGT_MXORDR 10
# endif
# ifndef DGT_MXBASE
   # define DGT_MXBASE 10 
# endif
/*
# pragma OPTIMIZE OFF
*/
/*
# pragma OPTIMIZE ON
# pragma OPT_LEVEL 2 
*/
/*
# include "consts.h"
*/
typedef struct
{
   long 
      n;

   char 
      base,
      order;

   char 
      dgt[DGT_MXORDR+ONE],
      cks[DGT_MXBASE+ONE];

}  DIGITS;

/*============================================================================*/

DIGITS *
dgtize( DIGITS *inp, const char *option )
{
   static DIGITS 
      rtn = {null},
     *rpt = &rtn;

   static char 
      *ptr,
     **endp = null;

   static char
      ii = null,
      jj = null;

   static char
      base = null,
      order = null;

   static long
      rmdr = null,
      power = null;

/*================= end of declaration part ==================================*/

   ptr = ( char * ) calloc( STS_SIZE, ONE );

/*----------------------------------------------------------------------------*/

   rmdr = ( inp->n );
   base = ( inp->base );
   order = ( inp->order );

   ( rpt->n ) = rmdr;
   ( rpt->base ) = base;
   ( rpt->order ) = order;

   if( *option == 'd' )/* digitize number rmdr = (inp->n) */
   {
      power = 1;
      ii = null;
      while( ii < order ) /* power = base^order */
      {
         power *= base;
         ii++ ;
      };

      if( power < rmdr )
      {
         while(( power < rmdr )&&( ii < order )) /* power = base^order */
         {
            power *= base;
            ii++ ;
         };

         fprintf( stderr, "\n\n Error message from function %s:", DO_DGTIZE );
         fprintf( stderr, "\n Transferred number %ld exceeds digitization "
            "domain %ld^%ld !!! ", rmdr, base, order );
         fprintf( stderr,
            "\n [ Transfer higher order: %ld <= order.]\n ", ii );

         exit( EXIT_FAILURE );
      };

/*............................................................................*/
/* here starts the job: */

      ii = null;
      while( ii < base ) /* reset checksum to null */
      {
         rtn.cks[ii] = null;
         ii++;
      };

      ii = order - ONE;
      while( null <= ii )
      {
         power /= base;
         jj = ( char ) ( rmdr / power );
         ( rpt->dgt[ii] ) = jj;
         ( rpt->cks[jj] )++ ;

         rmdr -= ( long ) jj*power;
         ii-- ;
      };

      return rpt;
   }
   else if ( *option == 'r' ) /* restore number */
   {
      ii = null;
      while( ii < base ) /* reset checksum to null */
      {
         ( rpt->cks[ii] ) = null;
         ii++;
      };

      rmdr = null;
      power = 1;
      ii = null;
      while( ii < order )
      {
         jj = ( inp->dgt[ii] );
         rmdr += jj*power;
         ( rpt->cks[jj] )++ ; /* compute checksum */
         power *= base;
         ii++ ;
      };
      ( rpt->n ) = rmdr;

      return rpt;
   }
   else 
   {
      fprintf( stderr, "\n\n Error message from function %s:", DO_DGTIZE );
      fprintf( stderr, "\n Unknown or unspecified option !!!" );
      fprintf( stderr, "\n [ Known options are: *option = 'd'igitize "
         "or *option = 'r'estore .] \n " );

      return null;
   };
}
/************************** end of function dgtize(*) *************************/
