/* [ file: lotos.h ] */
# define DO_LOTOS "lotos(*)"
/*******************************************************************************
*                                                                              *
*   ANSI-C function lotos(*)                                                   *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations ]          *
*   Release 1.0                                                                *
*                                                                              *
*   This function converts given integer lngint into                           *
*   an  ASCII character string  pointed to by lngstr                           *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 18, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*============================================================================*/

char *lotos ( long lngint, char length )
{
   static long
      nn = null,
      dd = null;

   static signed char
      ii    = null,
      jj    = null,
      ssize = null;

   static char 
      *lngstr;
      
   static char
      **endp;

   static short ind = null;

   lngstr = ( char *) calloc ( SHS_SIZE, ONE );

/*............................................................................*/

   nn    = ONE;
   ii    = null;
   ssize = null;
   if ( lngint < null )
   {
      lngstr[ssize] = '-';
      lngint = - lngint;
      ssize = ONE;
      ii = ONE;
   };
   do
   {
      lngstr[ssize] = 48; /* 48:  ASCII sign '0' */
      nn *= 10;
      ssize++;
   }  while (( ssize < SHS_SIZE )&&( nn <= lngint ));
   do
   {
      nn /=10;
      dd = lngint / nn ;
      lngstr[ii] = ( char ) dd + 48;        /* (+ 48) converts  digit  dd     */
      lngint -= dd*nn;                      /*        into ASCII sign 'dd'    */
      ii++;
   }  while ( ii < ssize );
   lngstr[ssize] = null;

/* write trailing ZEROs if length > ssize: */

   if ( length > ssize ) 
   {
      jj = length;

      while ( ssize >= null )
      { 
         lngstr[jj] = lngstr[ssize]; 
         jj--;
         ssize--;
      };
      while ( jj >= null )
      {
         lngstr[jj] = 48;
         jj--; 
      };
   };

   return lngstr; 
}
/*============================================================================*/
/************* end of long-to-string conversion function lotos(*) *************/
