/* [ file: fluidstate.b ] */
/*******************************************************************************
*                                                                              *
*   Function body fluidstate(*)                                                *
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
/* defaults: ... */
/*============================================================================*/

double
fluidstate( CLUSTER *cl1, CLUSTER *cl2 )
{
/* declarations: */
/* register type */

   unsigned char 
      aa = 't',
      bb = 'd';

   static long
      hh = null;

   static double
      prs = ZERO,
      cmprc = ZERO, /* compression coefficient */
      rmean = ZERO; /* mean density */

/* structure pointers: */
/* geometric and material parameter structure */
/*
   static struct hcrsmx
      *hsp = NULL;
*/
   CLUSTER
     *tmp = NULL,
     *dns = NULL;
/*............................................................................*/
/*
*//* presently not used [ but maybe in a future implementation ... ]
*//*
*//*   CLUSTER 
*//*      *gradnt( CLUSTER *cls, long hh ),
*//*      *gradfc( CLUSTER *cls, long hh );
*/
/*----------------------------------------------------------------------------*/
   aa = ( cl1->par );
   bb = ( cl2->par );

   switch( aa )
   {
     default:
      break;

     case 't':
     case 'T':
      switch( bb )
      {
	default:
         break;

        case 'd':
        case 'D':
/* compute pressure from temperature and density */
         tmp = cl1;
	 dns = cl2;
	 hh = ( dns->hh );
	 rmean = ( dns->hsp->rm[hh] );
	 cmprc = ( dns->hsp->cm[hh] );

         prs = (( dns->node ) - rmean )*\
            (( dns->hsp->bm[hh] )*(( tmp->node ) + 273.3 ));

         if ( 1.0e-77 < fabs( cmprc )) 
            prs /= ( cmprc * rmean );

         return prs;
         break;
      };
      break;

     case 'd':
     case 'D':
      switch( bb )
      {
        default:
         break;

        case 't':
        case 'T':
/* compute pressure from temperature and density */
	 tmp = cl2;
         dns = cl1;
	 hh = ( dns->hh );
	 rmean = ( dns->hsp->rm[hh] );
	 cmprc = ( dns->hsp->cm[hh] );

         prs = (( dns->node ) - rmean )*\
            (( dns->hsp->bm[hh] )*(( tmp->node ) + 273.3 ));

         if ( 1.0e-77 < fabs( cmprc )) 
            prs /= ( cmprc * rmean );

         return prs;
         break;
      };
      break;
   };
/*............................................................................*/
   fprintf( stderr, "\n Error in function %s :", __func__ );
   fprintf( stderr, "\n Insufficient parameter input !!!"  );
   fprintf( stderr, "\n [ lacking temperature and density "
      "information ]\n\n" );

   exit( EXIT_FAILURE );

   return ZERO;
 }
/*============================================================================*/
# endif /* (( CMPRSSBL != 0 )\
          &&( BSQAPRX == 0 )) */
/************************* end of function idelgas(*) *************************/
/*************************** end of file idealgas.b ***************************/
