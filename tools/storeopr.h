/* [ file: storeopr.h ] */
/*******************************************************************************
*                                                                              *
*   DSC model operation/computation mode storage function store_operts(*)      *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 18, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# ifndef LINLEN
   # define LINLEN 53
# endif
# ifndef _DEBUG
   # define _DEBUG 1
# endif
/*----------------------------------------------------------------------------*/
# define STORE_OPR( FN, MD ) \
{ \
   fprintf((FN), spformat, ( spt->tpt->name )); \
   fprintf((FN), spformat, ( spt->tpt->text )); \
 \
/* if ( null == strcmp( fleptr, tmpfle )) */ \
   if (( MD ) == 't' ) \
   { \
      fprintf((FN), "\n%-3d       %s", trf.c[null], trf.ctx[null] ); \
 \
      fprintf((FN), "\n" ); \
      ii = null; do \
      { \
         fprintf((FN), "%c", '-' ); \
      } while(( ii++ ) < ( LINLEN+14 )); \
 \
      for ( ii=ONE; ii<=trf.c[null]; ii++ ) \
         fprintf((FN), "\n%-3d <---  c%03d|%-s", \
            trf.c[ii], ii, trf.ctx[ii] ); \
   }  \
   else \
   { \
      fprintf((FN), "\n%-3d |%s", trf.c[null], trf.ctx[null] ); \
 \
      fprintf((FN), "\n" ); \
      ii = null; do \
      { \
         fprintf((FN), "%c", '-' ); \
      } while(( ii++ ) < ( LINLEN+9 )); \
 \
      for ( ii=ONE; ii<=trf.c[null]; ii++ ) \
         fprintf((FN), "\nc%03d|%s: %3d", ii, trf.ctx[ii], trf.c[ii] ); \
   }; \
   fprintf((FN), "\n" ); \
}
/*============================================================================*/

short store_operts( char *filename, char mode )
{
/* allusions: */
/*
   extern FORMSTATE *spt;
*/
/* declarations: */

   static FILE 
     *operats = NULL;

   static long 
      ll = null;

   static short
      ii = null;

   static char
     *timeptr = null,
      ptr[STS_SIZE] = {null},
      fleptr[STS_SIZE] = {null};

   static const char
     *scformat = "80%s",
     *spformat = "%s\n";

   time_t 
      nseconds = null,
     *timer = null;
/*............................................................................*/
# if _DEBUG == 1
   static char
     *spformat_debug = "\n%-17s";

   static double
      nmbr1 = ZERO,
      nmbr2 = ZERO,
      rsult = ZERO;

   static short
      counter = null;
# endif
/*----------------------------------------------------------------------------*/
/* open operations file: */

   strcpy( fleptr, filename );
   operats = fopen( fleptr, "w+" );
/*............................................................................*/
/* set stream unbuffered: */

/* setbuf( operats, NULL ); */

   ll = null;
   ii = null;
   while ( ii < IPT_MAXLBL )
   {
      fscanf( operats, scformat, ptr );

      if ( null == strncmp( ptr, "OPERATIONS", ( size_t ) SEVEN ))
      {
         ll = ftell( operats ) + ONE;
            break;
      }
      else
      {
         ii++;
         if( ii == IPT_MAXLBL )
         {
            if ( mode == 't' ) /* temporary file */
            {
               fprintf( operats, spformat, "The actually charged" );
               fprintf( operats, spformat, "---> OPERATIONS" );
            }
            else
               fprintf( operats, "%s%s\n", "OPERATIONS-", ( spt->flbl ));

            ll = ftell( operats );
         };
      };
   };

   fseek( operats, ll, SEEK_SET );

   STORE_OPR( operats, mode );
/*............................................................................*/
# if _DEBUG == 1

   fprintf( operats, "\n%-20s\n", "FL_DBL_PRECISION" );
   fprintf( operats, spformat_debug, "FLT_DIG" );
   fprintf( operats, " =% d", ( int ) FLT_DIG );
   fprintf( operats, spformat_debug, "FLT_MIN" );
   fprintf( operats, " =% e", FLT_MIN );
   fprintf( operats, spformat_debug, "FLT_MAX" );
   fprintf( operats, " =% e", FLT_MAX );
   fprintf( operats, spformat_debug, "FLT_EPSILON" );
   fprintf( operats, " =% e\n", FLT_EPSILON );
   fprintf( operats, spformat_debug, "DBL_DIG" );
   fprintf( operats, " =% d", ( int ) DBL_DIG );
   fprintf( operats, spformat_debug, "DBL_MIN" );
   fprintf( operats, " =% e", DBL_MIN );
   fprintf( operats, spformat_debug, "DBL_MAX" );
   fprintf( operats, " =% e", DBL_MAX );
   fprintf( operats, spformat_debug, "DBL_EPSILON" );
   fprintf( operats, " =% e\n", DBL_EPSILON );

   fprintf( operats, "\n%-32s\n", "DIGITS_ACCURACY [ double type ]" );

   nmbr1 = 1.0;
   nmbr2 = 1.0;
   counter = null;

   while(( nmbr1 + nmbr2 ) != nmbr1 )
   {
      ++counter;
      nmbr2 /= 10.;
   };
   fprintf( operats, spformat_debug, "in calculations:" );
   fprintf( operats, " =% d", counter );

   nmbr1 = 1.0;
   nmbr2 = 1.0;
   counter = null;

   while(1)
   {
      rsult = nmbr1 + nmbr2;

      if ( rsult == nmbr1 )
         break;

      ++counter;
      nmbr2 /= 10.;
   };
   fprintf( operats, spformat_debug, "in storage:" );
   fprintf( operats, " =% d\n", counter );

   fprintf( operats, "\n%-32s\n", "FP_ENVIRONMENT [FENV]" );

# ifdef FE_DIVBYZERO
   fprintf( operats, spformat_debug, "FE_DIVBYZERO" );
   fprintf( operats, " =% d", ( int ) FE_DIVBYZERO );
# endif
# ifdef FE_INEXACT
   fprintf( operats, spformat_debug, "FE_INEXACT" );
   fprintf( operats, " =% d", ( int ) FE_INEXACT );
# endif
# ifdef FE_INVALID
   fprintf( operats, spformat_debug, "FE_INVALID" );
   fprintf( operats, " =% d", ( int ) FE_INVALID );
# endif
# ifdef FE_OVERFLOW
   fprintf( operats, spformat_debug, "FE_OVERFLOW" );
   fprintf( operats, " =% d", ( int ) FE_OVERFLOW );
# endif
# ifdef FE_UNDERFLOW
   fprintf( operats, spformat_debug, "FE_UNDERFLOW" );
   fprintf( operats, " =% d", ( int ) FE_UNDERFLOW );
# endif
# ifdef FE_DOWNWARD
   fprintf( operats, spformat_debug, "FE_DOWNWARD" );
   fprintf( operats, " =% d", ( int ) FE_DOWNWARD );
# endif
# ifdef FE_TONEAREST
   fprintf( operats, spformat_debug, "FE_TONEAREST" );
   fprintf( operats, " =% d", ( int ) FE_TONEAREST );
# endif
# ifdef FE_TOWARDZERO
   fprintf( operats, spformat_debug, "FE_TOWARDZERO" );
   fprintf( operats, " =% d", ( int ) FE_TOWARDZERO );
# endif
# ifdef FE_UPWARD
   fprintf( operats, spformat_debug, "FE_UPWARD" );
   fprintf( operats, " =% d\n", ( int ) FE_UPWARD );
# endif
# endif /* _DEBUG == 1 */
/*............................................................................*/
   nseconds = time( timer );
   timeptr = ctime( &nseconds );

   fprintf( operats, "\nDSC model operations logfile %s", fleptr );
   fprintf( operats, " created:\n%24s", timeptr );
   fprintf( operats, "%c", EOF );

   fclose( operats );

   return null;
}
/*============================================================================*/
/********************** end of function store_operts(*) ***********************/
