/* [ file: readval.h ] */ 
# define DO_READVAL "readval(*)"
/*******************************************************************************
*                                                                              *
*   ANSI-C function readval(*)                                                 *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations ]          *
*   Release 1.0                                                                *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 18, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define READVALUES(II,FF,NN) \
{ \
   ii = (II); \
   kk = (FF); \
   while(( ii++ ) < kk ) \
   { \
      fscanf( evalfle, scformat, rptr ); \
\
      if ( TWO == ( vpt->read )) \
         fscanf( evalfle, scformat, iptr ); \
\
      mm = null; \
      while(( mm++ ) < ( state->idx[null] )) \
      { \
         if ( ii == ( state->idx[mm] )) \
         { \
            ( fpt->r[mm][(NN)] ) = strtod( rptr, endp ); \
\
            if ( TWO == ( vpt->read )) \
               ( fpt->i[mm][(NN)] ) = strtod( iptr, endp ); \
            else \
               ( fpt->i[mm][(NN)] ) = ZERO; \
         }; /* end if ii == ( state->idx[mm] ) */ \
      }; /* end while (( mm++ ) < ( state->idx[null] )); */ \
   }; /* end while (( ii++ ) <  kk ); */ \
}
/*============================================================================*/

EVALUATE * \
readval( POSTSTATE *state, char option )
{
/* allusions: */
/*
   extern POSTSTATE eval;
*/
/* declarations: */

   static FILE 
      *evalfle = NULL;

   static EVALUATE
      *vpt = NULL;

   static FFT
      *fpt = NULL;

   static long 
      hh = null,
      nn = null;

   static char
      sgn,
     *scformat = "%80s",
    **endp = null;

   static char
/*
      ptr[STS_SIZE] = {null},
*/
      rptr[STS_SIZE] = {null},
      iptr[STS_SIZE] = {null};

   static char 
      slsh = 47,  /* slash     */
      dash = 93,  /* stroke    */
     bslsh = 92;  /* back slash */

   static short
      lines = THREE,
      columns = 10,
      lbnd = 0,
      cbnd = 8;  /* <- cbnd = columns - 2 */

   static long
      ii = null,
      kk = null,
      mm = null,
      pp = null;

# if DSC_HCRMDE != 0
   static long 
      ht = null;

   static short
      cc = null;
# endif
/*----------------------------------------------------------------------------*/
   vpt = ( state->vpt );
   fpt = ( vpt->fpt );

   if ( option == 'r' )
   {
      evalfle = fopen(( state->file ), "r+" );
      fseek( evalfle, ( state->fleofs ), SEEK_SET );
      
      printf( "\n Please wait a moment !" );
      printf( "\n [ Reading data from file %s.]", ( state->file ));
      
      hh = null;
/*...........................................................................*/\
# if DSC_HCRMDE != 0
      ht = null;
# endif
/*...........................................................................*/\
      nn = null;
      while(( nn++ ) < ( vpt->n ))
      {
         if ((( vpt->ni ) <= nn )
           &&( nn <= ( vpt->nf ))
           &&( null < ( vpt->r )))
         {
/* the actual time or iteration [ Maxwell field algorithm ] */
            if ( null < ( state->fldprd ))
            {
               fscanf( evalfle, scformat, rptr );
               fscanf( evalfle, scformat, iptr ); /* string "...seconds..." */
                                                  /* or "...iterations..."  */
            };
            READVALUES( null, ( state->fldprd ), hh );
	    hh++;
         };
/*...........................................................................*/\
# if DSC_HCRMDE != 0  
         if ((( vpt->nj ) <= nn )
           &&( nn <= ( vpt->nt ))
           &&( null < ( vpt->rc )))
         {
/* the actual time [ heat and fluid algorithm ] */
            if (( state->fldprd ) < ( state->period ))
            {
               fscanf( evalfle, scformat, rptr );
               fscanf( evalfle, scformat, iptr ); /* string "...seconds..." */
            };
            READVALUES(( state->fldprd ), ( state->period ), ht );
            ht++;
         };
# endif /* DSC_HCRMDE != 0  */
/*...........................................................................*/\
      }; /* end while( nn <= ( vtp->n )) */
      
      printf( "\n\n End of file statement "
         "overtaken from file %s:", ( state->file ));

      PRBLDCLR( "" );
      while(( sgn = fgetc( evalfle )) != EOF )
         printf( "%c", sgn );
      PRNORMAL( "" );

      fclose( evalfle );
/*...........................................................................*/\
      return vpt;
   }
   else if ( option == 'd' ) /* display */
   { 
      mm = null;

      if ((( vpt->nf ) < ( vpt->ni ))
        ||(( vpt->r ) <= null ))
         goto flds_listo; /* no stored electric values to be evaluated */

      if ( null < ( vpt->nep ))
      {
         if ( *( vpt->mode_ep ) == 'a' ) 
            mm++; 

         if (( ONE < ( vpt->nep ))&&( *( vpt->mode_ep ) == 'a' ))
            printf( "\n E-field averaging over %05d ports: >-----"
               "----------------------------------->", ( vpt->nep ));
         else 
            printf( "\n E-field evaluation at %05d ports: >------"
               "----------------------------------->", ( vpt->nep ));

         lines = (( vpt->nep )/columns ) + ONE;

         for ( kk=null; kk<lines; kk++ )
         {
            printf( "\n\n ->index:" );
            pp = kk*columns;

            for ( ii=ONE; ii<=columns; ii++ )
            {
               if (( pp < ( vpt->nep ))&&( *( vpt->mode_ep ) != 'a' )) 
                  mm++;

               if (( lbnd < kk )&&( cbnd < ii ))
               {
                  if ( *( vpt->mode_ep ) != 'a' )
                  {
                     if (( pp + ONE ) < ( vpt->nep ))
                        mm += (( vpt->nep ) - pp - ONE );
                  };
                  printf( "   ***%c", bslsh );
                  printf( "%6ld%c", ( long ) mm, bslsh );

                  goto E_cells;
               };

               if ( pp < ( vpt->nep ))  
                  printf( "%6ld%c", ( long ) mm, bslsh );

               if (( vpt->nep ) <= ( ++pp ))  
                  goto E_cells;
            };

           E_cells:

            printf( "\n -> cell:" );

            pp = kk*columns;

            for ( ii=ONE; ii<=columns; ii++ )
            {
               if (( lbnd < kk )&&( cbnd < ii ))
               {
                  printf( "   ***%c", dash );
                  printf( "%6ld%c", ( long )( vpt->mep[(( vpt->nep )-ONE )] ),
                     dash );
                  goto E_ports;
               };

               if ( pp < ( vpt->nep ))   
                  printf( "%6ld%c", ( long )( vpt->mep[pp] ), dash ); 

               if (( vpt->nep ) <= ( ++pp ))  
                  goto E_ports; 
            };

           E_ports:

            printf( "\n -> port:" );
            pp = kk*columns;

            for ( ii=ONE; ii<=columns; ii++ )
            {
               if (( lbnd < kk )&&( cbnd < ii ))
               {
                  printf( "   ***%c", slsh );
                  printf( "%6ld%c\n", ( long )( vpt->pep[(( vpt->nep )-ONE )] ),
                     slsh );
                  goto E_nodes1;
               };

               if ( pp < ( vpt->nep ))   
                  printf( "%6ld%c", ( long )( vpt->pep[pp] ), slsh ); 

               if (( vpt->nep ) <= ( ++pp ))  
               {
                  printf( "\n" );
                  goto E_nodes1;
               };
            };
         };
      }; /* end if null < ( vpt->nep ) */

     E_nodes1:

      if ( null < ( vpt->nen ))
      {
         if ( *( vpt->mode_en ) == 'a' ) 
            mm++;

         if (( ONE < ( vpt->nen ))&&( *( vpt->mode_en ) == 'a' ))
            printf( "\n E-field averaging over %05d nodes: >-----"
               "----------------------------------->", ( vpt->nen ));
         else
            printf( "\n E-field evaluation at %05d nodes: >------"
               "----------------------------------->", ( vpt->nen ));

         lines = (( vpt->nen )/columns ) + ONE;

         for ( kk=null; kk<lines; kk++ )
         {
            printf( "\n\n ->index:" );
            pp = kk*columns;

            for ( ii=ONE; ii<=columns; ii++ )
            {
               if (( pp < ( vpt->nen ))&&( *( vpt->mode_en ) != 'a' )) 
                  mm++;

               if (( lbnd < kk )&&( cbnd < ii ))
               {
                  if ( *( vpt->mode_en ) != 'a' ) 
                  {
                     if (( pp + ONE ) < ( vpt->nen ))
                        mm += (( vpt->nen ) - pp - ONE );
                  };
                  printf( "   ***%c", bslsh );
                  printf( "%6ld%c", ( long ) mm, bslsh );
                  goto E_nodes2;
               };

               if ( pp < ( vpt->nen )) 
                  printf( "%6ld%c", ( long ) mm, bslsh ); 

               if (( vpt->nen ) <= ( ++pp ))  
                  goto E_nodes2;
            };

           E_nodes2:

            printf( "\n -> cell:" );
            pp = kk*columns;

            for ( ii=ONE; ii<=columns; ii++ )
            {
               if (( lbnd < kk )&&( cbnd < ii ))
               {
                  printf( "   ***%c", dash );
                  printf( "%6ld%c", ( long )( vpt->men[(( vpt->nen )-ONE )] ),
                     dash );

                  goto E_comps;
               };

               if ( pp < ( vpt->nen ))   
                  printf( "%6ld%c", ( long )( vpt->men[pp] ), dash ); 

               if (( vpt->nen ) <= ( ++pp ))  
                  goto E_comps; 
            };

           E_comps: 

            printf( "\n -> vect:" );
            pp = kk*columns;

            for ( ii=ONE; ii<=columns; ii++ )
            {
               if (( lbnd < kk )&&( cbnd < ii ))
               {
                  printf( "   ***%c", slsh );
                  printf( "%6c%c\n",
                     (( vpt->cen[(( vpt->nen )-ONE )] )+117 ), slsh );
                  goto H_field;
               };

               if ( pp < ( vpt->nen ))
                  printf( "%6c%c", 
                     (( vpt->cen[pp] )+117 ), slsh );

               if (( vpt->nen ) <= ( ++pp ))  
               {
                  printf( "\n" );
                  goto H_field;
               };
            };
         };
      }; /* end if null < ( vpt->nen ) */

     H_field:

      if ( null < ( vpt->nhp ))
      {
         if ( *( vpt->mode_hp ) == 'a' ) 
            mm++;

         if (( ONE < ( vpt->nhp ))&&( *( vpt->mode_hp ) == 'a' ))
            printf( "\n H-field averaging over %05d ports: >-----"
               "----------------------------------->", ( vpt->nhp ));
         else
            printf( "\n H-field evaluation at %05d ports: >------"
               "----------------------------------->", ( vpt->nhp ));

         lines = (( vpt->nhp )/columns ) + ONE;

         for ( kk=null; kk<lines; kk++ )
         {
            printf( "\n\n ->index:" );
            pp = kk*columns;

            for ( ii=ONE; ii<=columns; ii++ )
            {
               if (( pp < ( vpt->nhp ))&&( *( vpt->mode_hp ) != 'a' )) 
                  mm++;

               if (( lbnd < kk )&&( cbnd < ii ))
               {
                  if ( *( vpt->mode_hp ) != 'a' )
                  {
                     if (( pp + ONE ) < ( vpt->nhp ))
                        mm += (( vpt->nhp ) - pp - ONE );
                  };
                  printf( "   ***%c", bslsh );
                  printf( "%6ld%c", ( long ) mm, bslsh );

                  goto H_cells;
               };

               if ( pp < ( vpt->nhp ))   
                  printf( "%6ld%c", ( long ) mm, bslsh );

               if (( vpt->nhp ) <= ( ++pp ))  
                  goto H_cells;
            };

           H_cells:

            printf( "\n -> cell:" );
            pp = kk*columns;

            for ( ii=ONE; ii<=columns; ii++ )
            {
               if (( lbnd < kk )&&( cbnd < ii ))
               {
                  printf( "   ***%c", dash );
                  printf( "%6ld%c", ( long )( vpt->mhp[(( vpt->nhp )-ONE )] ),
                     dash );
                  goto H_ports;
               };

               if ( pp < ( vpt->nhp ))   
                  printf( "%6ld%c", ( long )( vpt->mhp[pp] ), dash ); 

               if (( vpt->nhp ) <= ( ++pp ))  
                  goto H_ports; 
            };

           H_ports:

            printf( "\n -> port:" );
            pp = kk*columns;

            for ( ii=ONE; ii<=columns; ii++ )
            {
               if (( lbnd < kk )&&( cbnd < ii ))
               {
                  printf( "   ***%c", slsh);
                  printf( "%6ld%c\n", ( long )( vpt->php[(( vpt->nhp )-ONE )] ),
                     slsh );
                  goto H_nodes1; 
               };

               if ( pp < ( vpt->nhp ))   
                  printf( "%6ld%c", ( long )( vpt->php[pp] ), slsh ); 

               if (( vpt->nhp ) <= ( ++pp ))  
               {
                  printf( "\n" );
                  goto H_nodes1; 
               };
            };
         };
      }; /* end if null < ( vpt->nhp ) */

     H_nodes1:

      if ( null < ( vpt->nhn ))
      {
         if ( *( vpt->mode_hn ) == 'a' ) 
            mm++;

         if (( ONE < ( vpt->nhn ))&&( *( vpt->mode_hn ) == 'a' ))
            printf( "\n H-field averaging over %05d nodes: >-----"
               "----------------------------------->", vpt->nhn );
         else
            printf( "\n H-field evaluation at %05d nodes: >------"
               "----------------------------------->", vpt->nhn );

         lines = (( vpt->nhn )/columns ) + ONE;

         for ( kk=null; kk<lines; kk++ )
         {
            printf( "\n\n ->index:" );
            pp = kk*columns;

            for ( ii=ONE; ii<=columns; ii++ )
            {
               if (( pp < ( vpt->nhn ))&&( *( vpt->mode_hn ) != 'a' )) 
                  mm++;

               if (( lbnd < kk )&&( cbnd < ii ))
               {
                  if ( *( vpt->mode_hn ) != 'a' )
                  {
                     if (( pp + ONE ) < ( vpt->nhn ))
                        mm += (( vpt->nhn ) - pp - ONE );
                  };
                  printf( "   ***%c", bslsh );
                  printf( "%6ld%c", ( long ) mm, bslsh );
                  goto H_nodes2;
               };

               if ( pp < ( vpt->nhn ))   
                  printf( "%6ld%c", ( long ) mm, bslsh );

               if (( vpt->nhn ) <= ( ++pp ))
                  goto H_nodes2;
            };

           H_nodes2: 

            printf( "\n -> cell:" );
            pp = kk*columns;

            for ( ii=ONE; ii<=columns; ii++ )
            {
               if (( lbnd < kk )&&( cbnd < ii ))
               {
                  printf( "   ***%c", dash );
                  printf( "%6ld%c", ( long )( vpt->mhn[(( vpt->nhn )-ONE )] ),
                     dash );

                  goto H_comps;
               };

               if ( pp < ( vpt->nhn ))   
                  printf( "%6ld%c", ( long )( vpt->mhn[pp] ), dash ); 

               if (( vpt->nhn ) <= ( ++pp ))
                  goto H_comps; 
            };

           H_comps:

            printf( "\n -> vect:" );
            pp = kk*columns;

            for ( ii=ONE; ii<=columns; ii++ )
            {
               if (( lbnd < kk )&&( cbnd < ii ))
               {
                  printf( "   ***%c", slsh );
                  printf( "%6c%c\n",
                     (( vpt->chn[(( vpt->nhn )-ONE )] )+117 ), slsh );
                  goto flds_listo;
               };

               if ( pp < ( vpt->nhn ))
                  printf( "%6c%c",
                     (( vpt->chn[pp] )+117 ), slsh ); 

               if (( vpt->nhn ) <= ( ++pp ))
               {
                  printf( "\n" );
                  goto flds_listo;
               };
            };
         };
      }; /* end if null < ( vpt->nhn ) */

     flds_listo: ;
/*............................................................................*/
# if DSC_HCRMDE != 0 /* display evaluated current nodes: */

      if ((( vpt->nt ) < ( vpt->nj ))
        ||(( vpt->rc ) <= null ))
         goto all_listo; /* no stored heat&fluid values to be evaluated */

      cc = null; do
      { 
         if ( null < ( vpt->nhc[cc] ))
         {
            printf( "\n %02d. heat current evaluation at %05d faces: >------"
               "-------------------------->", ( cc+ONE ), ( vpt->nhc[cc] ));
            lines = (( vpt->nhc[cc] )/columns ) + ONE;

            for ( kk=null; kk<lines; kk++ )
            {
               printf( "\n\n ->index:" );
               pp = kk*columns;

               for ( ii=ONE; ii<=columns; ii++ )
               {
                  if ( pp < ( vpt->nhc[cc] ))
                     mm++;

                  if (( lbnd < kk )&&( cbnd < ii ))
                  {
                     if (( pp + ONE ) < ( vpt->nhc[cc] ))
                        mm += (( vpt->nhc[cc] ) - pp - ONE );

                     printf( "   ***%c", bslsh );
                     printf( "%6ld%c", ( long ) mm, bslsh );

                     goto HC_cells;
                  };

                  if ( pp < ( vpt->nhc[cc] ))   
                     printf( "%6ld%c", ( long ) mm, bslsh );

                  if (( vpt->nhc[cc] ) <= ( ++pp ))
                     goto HC_cells;
               };

              HC_cells:

               printf( "\n -> cell:" );
               pp = kk*columns;

               for ( ii=ONE; ii<=columns; ii++ )
               {
                  if (( lbnd < kk )&&( cbnd < ii ))
                  {
                     printf( "   ***%c", dash );
                     printf( "%6ld%c",
                        ( long )( vpt->mhc[cc][(( vpt->nhc[cc] )-ONE )] ),
                           dash );

                     goto HC_faces;
                  };

                  if ( pp < ( vpt->nhc[cc] ))
                     printf( "%6ld%c", ( long )( vpt->mhc[cc][pp] ), dash );

                  if (( vpt->nhc[cc] ) <= ( ++pp ))
                     goto HC_faces;
               };

              HC_faces:

               printf( "\n -> face:" );
               pp = kk*columns;

               for ( ii=ONE; ii<=columns; ii++ )
               {
                  if (( lbnd < kk )&&( cbnd < ii ))
                  {
                     printf( "   ***%c", slsh );

                     printf( "%6ld%c\n", ( long ) \
                        ( vpt->fhc[cc][(( vpt->nhc[cc] )-ONE )] ), slsh );

                     goto TF_temps;
                  };

                  if ( pp < ( vpt->nhc[cc] ))
                     printf( "%6ld%c", ( long ) ( vpt->fhc[cc][pp] ), slsh );

                  if (( vpt->nhc[cc] ) <= ( ++pp ))
                  {
                     printf( "\n" );
                     goto TF_temps;
                  };
               }; /* end: for ( ii = ONE; ii <= columns ... ) */
            }; /* end: for ( kk = null; kk < lines ... ) */
         }; /* end if null < ( vpt->nhc[cc] ) */

        TF_temps: ;
/*............................................................................*/
# if DSC_FCTEMP == 1
         if ( null < ( vpt->ntf[cc] ))
         {
            printf( "\n %02d. face temperature evaluation at %05d faces: >--"
               "-------------------------->", ( cc+ONE ), ( vpt->ntf[cc] ));
            lines = (( vpt->ntf[cc] )/columns ) + ONE;

            for ( kk=null; kk<lines; kk++ )
            {
               printf( "\n\n ->index:" );
               pp = kk*columns;

               for ( ii=ONE; ii<=columns; ii++ )
               {
                  if ( pp < ( vpt->ntf[cc] ))
                     mm++;

                  if (( lbnd < kk )&&( cbnd < ii ))
                  {
                     if (( pp + ONE ) < ( vpt->ntf[cc] ))
                        mm += (( vpt->ntf[cc] ) - pp - ONE );

                     printf( "   ***%c", bslsh );
                     printf( "%6ld%c", ( long ) mm, bslsh );

                     goto TF_cells;
                  };

                  if ( pp < ( vpt->ntf[cc] ))   
                     printf( "%6ld%c", ( long ) mm, bslsh );

                  if (( vpt->ntf[cc] ) <= ( ++pp ))
                     goto TF_cells;
               };

              TF_cells:

               printf( "\n -> cell:" );
               pp = kk*columns;

               for ( ii=ONE; ii<=columns; ii++ )
               {
                  if (( lbnd < kk )&&( cbnd < ii ))
                  {
                     printf( "   ***%c", dash );
                     printf( "%6ld%c",
                        ( long )( vpt->mtf[cc][(( vpt->ntf[cc] )-ONE)] ),
                           dash );

                     goto TF_faces;
                  };

                  if ( pp < ( vpt->ntf[cc] ))
                     printf( "%6ld%c", ( long )( vpt->mtf[cc][pp] ), dash );

                  if (( vpt->ntf[cc] ) <= ( ++pp ))
                     goto TF_faces;
               };

              TF_faces:

               printf( "\n -> face:" );
               pp = kk*columns;

               for ( ii=ONE; ii<=columns; ii++ )
               {
                  if (( lbnd < kk )&&( cbnd < ii ))
                  {
                     printf( "   ***%c", slsh );

                     printf( "%6ld%c\n", ( long ) \
                        ( vpt->ftf[cc][(( vpt->ntf[cc] )-ONE )] ), slsh );

                     goto TN_temps;
                  };

                  if ( pp < ( vpt->ntf[cc] ))
                     printf( "%6ld%c", ( long ) ( vpt->ftf[cc][pp] ), slsh );

                  if (( vpt->ntf[cc] ) <= ( ++pp ))
                  {
                     printf( "\n" );
                     goto TN_temps;
                  };
               }; /* end: for ( ii = ONE; ii <= columns ... ) */
            }; /* end: for ( kk = null; kk < lines ... ) */
         }; /* end if null < ( vpt->ntf[cc] ) */

         TN_temps:
	 
# endif /* DSC_FCTEMP == 1 */
/*............................................................................*/
         if ( null < ( vpt->ntn[cc] ))
         {
            printf( "\n %02d. node temperature evaluation at %05d nodes: >--"
               "-------------------------->", ( cc+ONE ), ( vpt->ntn[cc] ));
            lines = (( vpt->ntn[cc] )/columns ) + ONE;

            for ( kk=null; kk<lines; kk++ )
            {
               printf( "\n\n ->index:" );
               pp = kk*columns;

               for ( ii=ONE; ii<=columns; ii++ )
               {
                  if ( pp < ( vpt->ntn[cc] ))
                     mm++;

                  if (( lbnd < kk )&&( cbnd < ii ))
                  {
                     if (( pp + ONE ) < ( vpt->ntn[cc] ))
                        mm += (( vpt->ntn[cc] ) - pp - ONE );

                     printf( "   ***%c", bslsh );
                     printf( "%6ld%c", ( long ) mm, bslsh );

                     goto TN_nodes;
                  };

                  if ( pp < ( vpt->ntn[cc] ))   
                     printf( "%6ld%c", ( long ) mm, bslsh );

                  if (( vpt->ntn[cc] ) <= ( ++pp ))
                     goto TN_nodes;
               };

              TN_nodes: 

               printf( "\n -> cell:" );
               pp = kk*columns;

               for ( ii=ONE; ii<=columns; ii++ )
               {
                  if (( lbnd < kk )&&( cbnd < ii ))
                  {
                     printf( "   ***%c", dash );
                     printf( "%6ld%c", ( long ) \
                        ( vpt->mtn[cc][(( vpt->ntn[cc] )-ONE )] ), dash );

                     goto TN_ports;
                  };

                  if ( pp < ( vpt->ntn[cc] ))
                     printf( "%6ld%c", ( long )( vpt->mtn[cc][pp] ), dash );

                  if (( vpt->ntn[cc] ) <= ( ++pp ))
                     goto TN_ports;
               };

              TN_ports:

               printf( "\n -> port:" );
               pp = kk*columns;

               for ( ii=ONE; ii<=columns; ii++ )
               {
                  if (( lbnd < kk )&&( cbnd < ii ))
                  {
                     printf( "   ***%c", slsh );
                     printf( "%6c%c\n", 45, slsh ); /* 45 = '-' */

                     goto hcrr_listo;
                  };

                  if ( pp < ( vpt->ntn[cc] ))
                     printf( "%6c%c", 45, slsh );

                  if (( vpt->ntn[cc] ) <= ( ++pp ))
                  {
                     printf( "\n" );
                     goto hcrr_listo;
                  };
               }; /* end: for ( ii = ONE; ii <= columns ... ) */
            }; /* end: for ( kk = null; kk < lines ... ) */
         }; /* end if null < ( vpt->ntn[cc] ) */

        hcrr_listo: ;
/*............................................................................*/
# if DSC_FLDMDE != 0 
/* display nodal velocities to be evaluated: */

         if ( null < ( vpt->nun[cc] ))
         {
            printf( "\n %02d. nodal velocity evaluation at %05d nodes: >----"
               "-------------------------->", ( cc+ONE ), ( vpt->nun[cc] ));
            lines = (( vpt->nun[cc] )/columns ) + ONE;

            for ( kk=null; kk<lines; kk++ )
            {
               printf( "\n\n ->index:" );
               pp = kk*columns;

               for ( ii=ONE; ii<=columns; ii++ )
               {
                  if ( pp < ( vpt->nun[cc] ))
                     mm++;

                  if (( lbnd < kk )&&( cbnd < ii ))
                  {
                     if (( pp + ONE ) < ( vpt->nun[cc] ))
                        mm += (( vpt->nun[cc] ) - pp - ONE );

                     printf( "   ***%c", bslsh );
                     printf( "%6ld%c", ( long ) mm, bslsh );

                     goto U_nodes;
                  };

                  if ( pp < ( vpt->nun[cc] ))   
                     printf( "%6ld%c", ( long ) mm, bslsh );

                  if (( vpt->nun[cc] ) <= ( ++pp ))
                     goto U_nodes;
               };

              U_nodes:

               printf( "\n -> cell:" );
               pp = kk*columns;

               for ( ii=ONE; ii<=columns; ii++ )
               {
                  if (( lbnd < kk )&&( cbnd < ii ))
                  {
                     printf( "   ***%c", dash );
                     printf( "%6ld%c",
                        ( long )( vpt->mun[cc][(( vpt->nun[cc] )-ONE )] ),
                           dash );

                     goto U_compts;
                  };

                  if ( pp < ( vpt->nun[cc] ))
                     printf( "%6ld%c", ( long )( vpt->mun[cc][pp] ), dash );

                  if (( vpt->nun[cc] ) <= ( ++pp ))
                     goto U_compts;
               };

              U_compts:

               printf( "\n -> cmpt:" );
               pp = kk*columns;

               for ( ii=ONE; ii<=columns; ii++ )
               {
                  if (( lbnd < kk )&&( cbnd < ii ))
                  {
                     printf( "   ***%c", slsh );

                     printf( "%6c%c\n",
                        (( vpt->cun[cc][( vpt->nun[cc] )-ONE] )+120 ), slsh );

                     goto flow_listo;
                  };

                  if ( pp < ( vpt->nun[cc] ))
                     printf( "%6c%c",
                        (( vpt->cun[cc][pp] )+120 ), slsh );

                  if (( vpt->nun[cc] ) <= ( ++pp ))
                  {
                     printf( "\n" );
                     goto flow_listo;
                  };
               }; /* end: for ( ii = ONE; ii <= columns ... ) */
            }; /* end: for ( kk = null; kk < lines ... ) */
         }; /* end if null < ( vpt->nun[cc] ) */

        flow_listo: ;

# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
      } while (( ++cc ) < DSC_HCRMDE );
/*............................................................................*/
      all_listo: ;
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
      
      ( vpt->rtn ) = null;
      return vpt;
   } /* end of display section */
   else /* unknown option */
   {
      ( vpt->rtn ) = ONE;
      return vpt;
   };
}
/*============================================================================*/
/*
# undef READVALUES
*/
# undef DISPEVLPTS
/************************ end of function readval(*) **************************/
