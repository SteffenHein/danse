/* [ file: values.h ] */
/*******************************************************************************
*                                                                              *
*   ANSI C subroutine values(*)                                                *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   DSC field evaluation and ( result ) file deposition function               *
*                                                                              *
*   (C) SHEIN; Bad Aibling, October 2007                   Steffen Hein        *
*   [ Update: April 05, 2022 ]                          <contact@sfenx.de>     *
*                                                                              *
*******************************************************************************/
# include "../math/maths.h"  /* 'my' computation environment headers */
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
/* A macro that gives a conveniently readable form to the string returned */
/* with [ the pointer of ] function ctime(&nseconds) */
# include "../tools/TIMEFORM.M"
/*----------------------------------------------------------------------------*/
# define VAL_DISP 1
/*----------------------------------------------------------------------------*/
/* if val.ni < VAL_INITLBL then val.ni = VAL_INITLBL */
# ifndef VAL_INITLBL
   # define VAL_INITLBL 1
# endif
/*----------------------------------------------------------------------------*/
# if DSC_HCRMDE != 0
   # ifndef DSC_FCTEMP
      # define DSC_FCTEMP 1 /* 1: store face temperatures */
   # endif               
   # define VAL_HSCSUM 2 /* 1: sum up and store skin effect heat sources      */
                 /* 2: in addition, sum up and store mesh cell heat sources   */
                 /* 3: sum up and store only mesh cell heat sources           */
   # define VAL_HSCVOL 1 /* 1: sum up and store mesh cell heat source volume  */
# endif          /* [ all operations at final time step ] */
/*----------------------------------------------------------------------------*/
static FILE *evalfle;
/*============================================================================*/

short values( const short jj )
{
/* allusions: */
/*
   extern struct solverstat solver;
   extern struct evaluation val;
*/
/* declarations: */

   static struct solverstat
      *state = &solver;

   static struct topology
      *tpt = &top;

   static short 
      ii = null,
      kk = null;

   static long
      mm = null,
      nn = null;

   static const char 
     *astrx = "***",
     *spformat = "%-27s\n",
     *dpformat = "%+.16e\n",
     *lpformat = "%ld\n",
     *scformat = "%80s";

   static signed char 
      evlfld = null;

# if DSC_HCRMDE != 0
   static signed char 
      evlhcr = null;
# endif

   static double
      mean_epr = ZERO,
      mean_epi = ZERO,
      mean_enr = ZERO,
      mean_eni = ZERO,
      mean_hpr = ZERO,
      mean_hpi = ZERO,
      mean_hnr = ZERO,
      mean_hni = ZERO;

   static time_t
      nseconds = null,
     *timer = null;

   static char 
         ptr[STS_SIZE] = {null},
         ctmptr[STS_SIZE] = {null},
         tmestr[STS_SIZE] = {null},
        *valptr = EVALUATION_FILE, /* cf. main program solver.c */
       **endp  = null;
/*............................................................................*/
# if DSC_HCRMDE != 0

# if (( VAL_HSCSUM == 2 )\
    ||( VAL_HSCSUM == 3 ))

# if VAL_HSCVOL == 1  
   static double
      vol = ZERO;
# endif

   static struct hcrsmx
     *hsp = NULL;

   static DSC_HCRRTS
     *hci = NULL;

# endif /* (( VAL_HSCSUM == 2 )||( VAL_HSCSUM == 3 )) */
/*............................................................................*/
   static short
      cc = null;
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* prototypes: */

   time_t
      time( time_t *timer );

   char
     *lotos ( long lngint, char length );
/*----------------------------------------------------------------------------*/
/* store computed values [ results ]: */

   if ( jj == - ONE )
   {
      tpt = &top;
/*............................................................................*/
/* Store electric port values: Note that all values are stored in Volt units. */
/* The pertinent electric field strengths E [Volt/m] are obtained by deviding */
/* the stored values through the length [in meters] of the respective port. */

      if (( val.ni <= ( state->nn ))
        &&(( state->nn ) <= val.nf ))
      {     
         if ( evlfld == ONE )
         {
/* the actual [ Maxwell field algorithm ] time or iteration step: */
            fprintf( evalfle, "%+.16e   ", solver.ttime );

            if ( solver.dmn == 't' )
               fprintf( evalfle, spformat, "--------seconds--------" );
            else
               fprintf( evalfle, spformat, "------iterations-------" );

            if ( null < val.nep )
            {
               for ( ii=null; ii<val.nep; ii++ )
               {
                  if ( val.read == ONE )
                  { /* aperiodic structure in time domain, e.g.*/
                     val.epr[ii] /= val.r;

                     if ( *val.mode_ep == 'a' )
                        mean_epr += val.epr[ii];
                     else
                        fprintf( evalfle, dpformat, val.epr[ii] );

                     val.epr[ii] = ZERO;
                  } /* end if val.read == ONE */
                  else /* if val.read == TWO */  
                  {   /* periodic structure */
                      /* or frequency_domain, e.g.*/
# if DSC_DOMAIN == 3
                     fprintf( evalfle, "%+.16e   ", val.epr[ii] );
                     fprintf( evalfle, "%+.16e\n", val.epi[ii] );
# else
                     val.epr[ii] /= val.r;
                     val.epi[ii] /= val.r;

                     if ( *val.mode_ep == 'a' )
                     {
                        mean_epr += val.epr[ii];
                        mean_epi += val.epi[ii];
                     }
                     else
                     {
                        fprintf( evalfle, "%+.16e   ", val.epr[ii] );
                        fprintf( evalfle, "%+.16e\n", val.epi[ii] );
                     }

                     val.epr[ii] = ZERO;
                     val.epi[ii] = ZERO;
# endif
                  }; /* end if val.read == TWO */
               };

               if ( *val.mode_ep == 'a' )
               {
                  if ( val.read == ONE )
	          { /* aperiodic structure in time domain, e.g.*/
                     mean_epr /= val.nep;                     
                     fprintf( evalfle, dpformat, mean_epr );
                     mean_epr = ZERO;
                  } 
                  else /* if val.read == TWO */
                  {
                     mean_epr /= val.nep;
                     mean_epi /= val.nep; 

                     fprintf( evalfle, "%+.16e   ", mean_epr );
                     fprintf( evalfle, "%+.16e\n", mean_epi );

                     mean_epr = ZERO;
                     mean_epi = ZERO;
                  }; /* end if val.read == TWO */
               }; /* if ( *val.mode_ep == 'a' ) */
	    }; /* end if ( null < val.nep ) */
/*............................................................................*/
/* store electric field node values: */

            if ( null < val.nen )
            {
               for ( ii=null; ii<val.nen; ii++ )
               {
                  if ( val.read == ONE ) 
                  {
                     val.enr[ii] /= val.r;

                     if ( *val.mode_en == 'a' )
                        mean_enr += val.enr[ii];
                     else
                        fprintf( evalfle, dpformat, val.enr[ii] );

                     val.enr[ii] = ZERO;
                  } 
                  else /* if val.read == TWO */ 
                  {
                     val.enr[ii] /= val.r;
                     val.eni[ii] /= val.r;

                     if ( *val.mode_en == 'a' )
                     {
                        mean_enr += val.enr[ii];
                        mean_eni += val.eni[ii];
                     }
                     else 
                     {
                        fprintf( evalfle, "%+.16e   ", val.enr[ii] );
                        fprintf( evalfle, "%+.16e\n", val.eni[ii] ); 
                     };

                     val.enr[ii] = ZERO;
                     val.eni[ii] = ZERO;
                  }; /* end if val.read == TWO */ 
               }; /* next ii */

               if ( *val.mode_en == 'a' )
               {
                  if ( val.read == ONE ) 
                  {
                     mean_enr /= val.nen;
                     fprintf( evalfle, dpformat, mean_enr );
                     mean_enr = ZERO;
                  } 
                  else /* if val.read == TWO */
                  {
                     mean_enr /= val.nen;
                     mean_eni /= val.nen;

                     fprintf( evalfle, "%+.16e   ", mean_enr );
                     fprintf( evalfle, "%+.16e\n", mean_eni );

                     mean_enr = ZERO;
                     mean_eni = ZERO;
                  }; /* end if val.read == TWO */
	       }; /* end if ( *val.mode_en == 'a' ) */
            }; /* end if null < val.nen */
/*............................................................................*/
/* Store magnetic port values. Note that all values are stored in Volt units. */
/* The pertinent magnetic field strengths H are retrieved from the following */
/* relations: */
/* H [A/m] = ( stored value ) [Volt] / ( L * sqrt( MY_VAC_/EPS_VAC )) */
/* where L is the length [in meters] of the transverse port at the same face */

            if ( null < val.nhp )
            {
               for ( ii=null; ii<val.nhp; ii++ ) 
               {
                  if ( val.read == ONE ) 
                  {
                     val.hpr[ii] /= val.r;

                     if ( *val.mode_hp == 'a' )
                        mean_hpr += val.hpr[ii];
                     else
                        fprintf( evalfle, dpformat, val.hpr[ii] );

                     val.hpr[ii] = ZERO; 
                  } /* end if val.read == ONE */
                  else /* if val.read == TWO */
                  {
/*............................................................................*/
# if DSC_DOMAIN == 3
                     fprintf( evalfle, "%+.16e   ", val.hpr[ii] );
                     fprintf( evalfle, "%+.16e\n", val.hpi[ii] );
# else
                     val.hpr[ii] /= val.r;
                     val.hpi[ii] /= val.r;

                     if ( *val.mode_hp == 'a' )
                     {
                        mean_hpr += val.hpr[ii];
                        mean_hpi += val.hpi[ii];
                     }
                     else
                     {
                        fprintf( evalfle, "%+.16e   ", val.hpr[ii] );
                        fprintf( evalfle, "%+.16e\n", val.hpi[ii] );
                     };

                     val.hpr[ii] = ZERO;
                     val.hpi[ii] = ZERO;
# endif
/*............................................................................*/
	          }; /* end if val.read == TWO */
               }; /* next ii */

               if ( *val.mode_hp == 'a' )
               {
                  if ( val.read == ONE ) 
                  {
                     mean_hpr /= val.nhp;
                     fprintf( evalfle, dpformat, mean_hpr );
                     mean_hpr = ZERO;
                  } 
                  else /* if val.read == TWO */
                  {
                     mean_hpr /= val.nhp;
                     mean_hpi /= val.nhp; 
   
                     fprintf( evalfle, "%+.16e   ", mean_hpr );
                     fprintf( evalfle, "%+.16e\n", mean_hpi );
            
                     mean_hpr = ZERO;
                     mean_hpi = ZERO;
                  }; /* end if val.read == TWO */
	       }; /* end if ( *val.mode_hp == 'a' ) */
	    }; /* end if ( null < val.nhp ) */
/*............................................................................*/
/* store magnetic field node values: */

            if ( null < val.nhn )
            {
               for ( ii=null; ii<val.nhn; ii++ )
               {
                  if ( val.read == ONE ) 
                  {
                     val.hnr[ii] /= val.r;

                     if ( *val.mode_hn == 'a' )
                        mean_hnr += val.hnr[ii];
                     else
                        fprintf( evalfle, dpformat, val.hnr[ii] );

                     val.hnr[ii] = ZERO;
                  } 
                  else /* if val.read == TWO */
                  {
                     val.hnr[ii] /= val.r;
                     val.hni[ii] /= val.r;

                     if ( *val.mode_hn == 'a' )
                     {
                        mean_hnr += val.hnr[ii];
                        mean_hni += val.hni[ii];
                     }
                     else
                     {
                        fprintf( evalfle, "%+.16e   ", val.hnr[ii] );
                        fprintf( evalfle, "%+.16e\n", val.hni[ii] );
                     };

                     val.hnr[ii] = ZERO;
                     val.hni[ii] = ZERO;
                  }; /* end if val.read == TWO */
               }; /* next ii */

               if ( *val.mode_hn == 'a' )
               {
                  if ( val.read == ONE ) 
                  {
                     mean_hnr /= val.nhn;
                     fprintf( evalfle, dpformat, mean_hnr );
                     mean_hnr = ZERO;
                  } 
                  else /* if val.read == TWO */
                  {
                     mean_hnr /= val.nhn;
                     mean_hni /= val.nhn; 

                     fprintf( evalfle, "%+.16e   ", mean_hnr );
                     fprintf( evalfle, "%+.16e\n", mean_hni );

                     mean_hnr = ZERO;
                     mean_hni = ZERO;
                  }; /* end if val.read == TWO */
	       }; /* end if ( *val.mode_hn == 'a' ) */
            }; /* end if ( null < val.nhn ) */
         }; /* end if ( null < evlfld ) */
      }; /* end if (( val.ni <= ( state->nn ))
                  &&(( state->nn ) < val.nf )) */
/*............................................................................*/
# if DSC_HCRMDE != 0 /* store temperatures and heat currents: */

      if (( val.nj <= ( state->nn ))
        &&(( state->nn ) <= val.nt ))
      {     
         if ( null < evlhcr )
         {
/* the actual [ thermal - fluid algorithm ] time: */
            fprintf( evalfle, "%+.16e   ", solver.hctme );
            fprintf( evalfle, spformat, "--------seconds--------" );

            cc = null; do
            {
/*............................................................................*/
/* store face values [ thermal - fluid ]: */

               for ( ii=null; ii<val.nhc[cc]; ii++ )
               {
                  if ( val.read == ONE )
                  {
                     val.hc[cc][ii] /= val.rc;
                     fprintf( evalfle, dpformat, val.hc[cc][ii] );

                     val.hc[cc][ii] = ZERO;
                  }
                  else /* if val.read == TWO */
                  {
                     val.hc[cc][ii] /= val.rc;

                     fprintf( evalfle, "%+.16e   ", val.hc[cc][ii] );
                     fprintf( evalfle, "%+.16e\n", ZERO );

                     val.hc[cc][ii] = ZERO;
                  }; /* end if val.read == TWO */
               }; /* end for ( ii=null; ii<val.nhc[cc]; ii++ ) */
/*............................................................................*/
# if DSC_FCTEMP == 1 /* store face temperatures: */

               for ( ii=null; ii<val.ntf[cc]; ii++ )
               {
                  if ( val.read == ONE )
                  {
                     val.tf[cc][ii] /= val.rc;

                     fprintf( evalfle, dpformat, val.tf[cc][ii] );

                     val.tf[cc][ii] = ZERO;
                  }
                  else /* if val.read == TWO */
                  {
                     val.tf[cc][ii] /= val.rc;

                     fprintf( evalfle, "%+.16e   ", val.tf[cc][ii] );
                     fprintf( evalfle, "%+.16e\n", ZERO );

                     val.tf[cc][ii] = ZERO;
                  }; /* end if val.read == TWO */
               }; /* end for ( ii=null; ii<val.ntf[cc]; ii++ ) */
# endif
/*............................................................................*/
/* store node temperatures: */

               for ( ii=null; ii<val.ntn[cc]; ii++ )
               {
                  if ( val.read == ONE )
                  {
                     val.tn[cc][ii] /= val.rc;
                     fprintf( evalfle, dpformat, val.tn[cc][ii] );

                     val.tn[cc][ii] = ZERO;
                  }
                  else /* if val.read == TWO */
                  {
                     val.tn[cc][ii] /= val.rc;

                     fprintf( evalfle, "%+.16e   ", val.tn[cc][ii] );
                     fprintf( evalfle, "%+.16e\n", ZERO );

                     val.tn[cc][ii] = ZERO;
                  }; /* end if val.read == TWO */
               }; /* next ii */
/*............................................................................*/
# if DSC_FLDMDE != 0
/* store nodal velocities: */

               for ( ii=null; ii<val.nun[cc]; ii++ )
               {
                  if ( val.read == ONE )
                  {
                     val.un[cc][ii] /= val.rc;
                     fprintf( evalfle, dpformat, val.un[cc][ii] );

                     val.un[cc][ii] = ZERO;
                  }
                  else /* if val.read == TWO */
                  {
                     val.un[cc][ii] /= val.rc;

                     fprintf( evalfle, "%+.16e   ", val.un[cc][ii] );
                     fprintf( evalfle, "%+.16e\n", ZERO );

                     val.un[cc][ii] = ZERO;
                  }; /* end if val.read == TWO */
               }; /* next ii */
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
            } while (( ++cc ) < DSC_HCRMDE );
         }; /* end if ( null < val.rc ) */
      }; /* end if (( val.nj <= ( state->nn ))
                  &&(( state->nn ) < val.nt )) */
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/

      return ONE;
   } 
/*............................................................................*/
/* open file dsc.val[joblbl[jj]], enter evaluation parameters */

   else  if ( null <= jj )
   { 
/* allusions: - [void ] */

/* declarations: */

      static long
         fleofs = null; /* offset from beginning of file dsc.val[h] */

      static const short 
         lbnd =  0,
         cbnd =  8, /* cbnd = clns - 2 */
         clns = 10; /* columns */

      static const char 
         slsh = 47,
        bslsh = 92;

      static signed char
         pp = null;

      static short
         llns = null;
/*
      static long
         mm = null,
         nn = null;
*/
/*----------------------------------------------------------------------------*/
/* memory allocations: - [void] */
/*............................................................................*/
# if DSC_LNGNAMES == 1
      strcpy( solver.val, solver.prfx );
      strcat( solver.val, valptr );
# else
      strncpy( solver.val, solver.prfx, VSS_SIZE );
      strncat( solver.val, valptr, ( SHS_SIZE - VSS_SIZE - TWO ));
# endif
/*............................................................................*/
      strcat( solver.val, lotos( joblbl[jj], null ));

      evalfle = fopen( solver.val, "r+" );

      if ( evalfle == null )
      { 
         fprintf( display, "\n Missing evaluation file %s !!!", solver.val );
         fprintf( display, "\n [ File created for error messages. ] " );

         evalfle = fopen( solver.val, "w+" );

         fprintf( evalfle, spformat, astrx );
          
         fprintf( evalfle,"\nError message from %s :\n\n", __func__ );

         for ( ii=null; ii<THREE; ii++ )
         {
            fprintf( evalfle, "%s%s\n", "unable_to_open_evaluation_file_", 
               solver.val ); 
         };
         fprintf( evalfle , "[ file created.]\n\n" ); 

         for ( ii=null; ii<THREE; ii++)
         {
            fprintf( evalfle, "%s%s\n",
               "evaluation_instructions_missing_on_file_", solver.val );
         };

         if ( null != strncmp( solver.errmsg, astrx, THREE ))
         {
            fprintf( evalfle,"\nError message from %s :\n\n", solver.fcterr );

            for ( ii=null; ii<THREE; ii++)
            {
               fprintf( evalfle, spformat, solver.errmsg );
            };
         }
         else
         {
            strncpy( solver.fcterr, __func__, SHS_SIZE );
            strncpy( solver.errmsg, solver.val, SHS_SIZE );
            strncat( solver.errmsg, "...:_error_on_opening_val.file",
                                  LGS_SIZE - SHS_SIZE );
         };

         nseconds = time( timer );
         strcpy( ctmptr, ctime( &nseconds ));
         fprintf( evalfle, "\n%s%d:\n%s\n\n\n", "abnormal end"
                           " of DSC job no.", joblbl[jj], ctmptr );
         fflush( evalfle );
         fclose ( evalfle );

         printf( "\n DSC job no.%d skipped.", joblbl[jj] );

         return null;
      };
/*............................................................................*/
      fscanf( evalfle, scformat, ptr );
      strncpy( val.name, ptr, SHS_SIZE - ONE );

      if ( null == strncmp( val.name, astrx, THREE ))
      {
         fprintf( display, "\n Missing evaluation instructions"
                           " on file %s !!! ", solver.val );

         if ( null != strncmp( solver.errmsg, astrx, THREE ))
         {
            fprintf( evalfle, "\nError message from %s :\n\n", solver.fcterr );

            for ( ii=null; ii<THREE; ii++)
            {
               fprintf( evalfle, spformat, solver.errmsg );
            };
         };

         fprintf( evalfle, "\nError message from %s :\n\n", __func__ );
         for ( ii=null; ii<THREE; ii++)
         {
            fprintf( evalfle, "%s%s\n",
               "evaluation_instructions_missing_on_file_", solver.val );
         };

         nseconds = time( timer );
         strcpy( ctmptr, ctime( &nseconds ));
         fprintf( evalfle, "\n%s%d:\n%s\n\n\n", "abnormal end"
                           " of DSC job no.", joblbl[jj], ctmptr );
         fclose ( evalfle );

         strncpy( solver.fcterr, __func__, SHS_SIZE );
         strncpy( solver.errmsg, solver.val, SHS_SIZE );
         strncat( solver.errmsg, "...:_missing_eval.instructs.",
                               LGS_SIZE - SHS_SIZE );

         fprintf( display, "\n DSC job no.%d skipped.", joblbl[jj] );

         return null;
      };
/*............................................................................*/
      if ( null == strncmp( solver.errmsg, astrx, THREE ))
      {
         if ( null != strncmp( val.name, top.name, THREE ))
         {
            fprintf( display, "\n Inconsistent system identifier "
               "in file %s !!!\n ", solver.val );
            strncpy( solver.fcterr, __func__, SHS_SIZE );
            strncpy( solver.errmsg, solver.val, SHS_SIZE ); 
            strncat( solver.errmsg, ":inconsistent_model_identifier_", 
               ( LGS_SIZE - SHS_SIZE ));  
         };
      };
/*............................................................................*/
/* operational instructions [ computation modes, etc.]: */

      fscanf( evalfle, scformat, ptr );
      strncpy( val.text, ptr, ( STS_SIZE - ONE ));

      fscanf( evalfle, scformat, ptr ); /* string "_________________________..." */
/* the number of iteration cycles: */

      fscanf( evalfle, scformat, ptr ); /* string "iteration_cycles" */
      fscanf( evalfle, scformat, ptr ); /* long integer string */                
      val.n = strtol( ptr, endp, DEC ); 

/* the first evaluated cycle [ Maxwell field ]*/

      fscanf( evalfle, scformat, ptr ); /* string "1st_evlt'd_cycle,_Mxwfld" */
      fscanf( evalfle, scformat, ptr );
      val.ni = strtol( ptr, endp, DEC );

      if ( val.ni < VAL_INITLBL )
         val.ni = VAL_INITLBL;

/* the last computed cycle [ Maxwell field ] */

      fscanf( evalfle, scformat, ptr ); /* string "Lst_cmpt'd_cycle,_Mxwfld" */
      fscanf( evalfle, scformat, ptr );
      val.nf = strtol( ptr, endp, DEC );

/* internal  repetition rate [ Maywell field ]*/

      fscanf( evalfle, scformat, ptr); /* string "iterations/cycle,_Mxwfld" */
      fscanf( evalfle, scformat, ptr );
      val.r  = strtol( ptr, endp, DEC );

      if ( val.r < null )
         val.r = null;
/*............................................................................*/
# if DSC_HCRMDE != 0
/* the first evaluated cycle [ thermal - fluid ] */

      fscanf( evalfle, scformat, ptr ); /* string "1st_evlt'd_cycle,_Th-Fld" */
      fscanf( evalfle, scformat, ptr );
      val.nj = strtol( ptr, endp, DEC ); 

/* the last computed cycle [ thermal - fluid ] */
      fscanf( evalfle, scformat, ptr ); /* string "Lst_cmpt'd_cycle,_Th-Fld" */
      fscanf( evalfle, scformat, ptr );
      val.nt = strtol( ptr, endp, DEC );

/* internal repetition rate [ thermal || fluid ]: */
      fscanf( evalfle, scformat, ptr); /* string "iterations/cycle,_Th-Fld" */
      fscanf( evalfle, scformat, ptr );
/*............................................................................*/
# if DSC_INTLCE == 1 /* interlaced inner Maxwell field and heat propagation*/
                      /* computational loop */
      val.rc = val.r;
# else
      val.rc  = strtol( ptr, endp, DEC );

      if ( val.rc < null )
         val.rc = null;
# endif
/*............................................................................*/
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
/* frequencies, times, and phases: */

      if ( solver.dmn == 't' ) /* time_domain */
         val.dt = val.r*solver.dt;       /* val.dt: cycle time step [seconds] */
      else /* if solver.dmn == 'f'requency_domain */
      {
         val.fr = solver.fr;            /* val.fr: frequency [Hertz] */
         val.ph = solver.ph;            /* val.ph: DSC phase shift [radians] */
      };
/*............................................................................*/
/* the following would need decoupled E/H field and current storage proceedrs */
/* [ suitably in different evaluation files, e.g.] - It is not used, here     */
/*
# if DSC_HCRMDE !=0
      cc = null; do
      {
         fscanf( evalfle, scformat, ptr );
         val.nc[cc] = strtol( ptr, endp, DEC );
         fscanf( evalfle, scformat, ptr );
         val.nic[cc] = strtol( ptr, endp, DEC );
         fscanf( evalfle, scformat, ptr );
         val.rc[cc] = strtol( ptr, endp, DEC );
      } while (( ++cc ) < DSC_HCRMDE );
# endif
*/
/*............................................................................*/
      fscanf( evalfle, scformat, ptr ); /* string "Maxwell_field_evaluation" */

      evlfld = null; /* initialize [set to ONE if there are evaluted ports] */

      fscanf( evalfle, scformat, ptr ); /* string "evaluated_ports" */
      fscanf( evalfle, scformat, ptr ); /* string "number..." */

      fscanf( evalfle, scformat, ptr ); /* string "E-ports" */
      fscanf( evalfle, scformat, ptr ); /* long integer string */
      val.nep = strtol( ptr, endp, DEC ); /* number of evaluated E-ports */

      if ( EVLEP < val.nep )
      {
         fprintf( display, "\n\n Message from function %s :", __func__ );
         fprintf( display, "\n Too many E ports to be evaluated "
            "in DSC mesh %s !!!", val.name );
         fprintf( display, "\n [ Maximum number is %ld = macro EVLEP "
            "in SOLVER.CONF or %s.", ( long ) EVLEP, __func__ );
         fprintf( display, "\n - Change macro only in compliance "
            "with memory resources !\n " );
         strncpy( solver.fcterr, __func__, SHS_SIZE );
         strncpy( solver.errmsg, solver.val, SHS_SIZE );
         strncat( solver.errmsg, ":too_many_E_ports_to_be_evaluated_", 
            ( LGS_SIZE - SHS_SIZE ));

         return null;
      }
      else if ( null < val.nep )
      {
         if ( null < val.r )
            evlfld = ONE;
      
         if ( ONE < val.nep )
         {
            fscanf( evalfle, scformat, ptr ); /* mean_value operation mark */

            if ( null != strstr( ptr, "average" ))
               strncpy( val.mode_ep, "average", VSS_SIZE );
            else
               strncpy( val.mode_ep, "individual", VSS_SIZE );
         };
      };

      fscanf( evalfle, scformat, ptr ); /* string "E-nodes" */
      fscanf( evalfle, scformat, ptr );
      val.nen = strtol( ptr, endp, DEC ); /* number of E-nodes evaluated */

      if ( EVLEN < val.nen )
      {
         fprintf( display, "\n\n Message from function %s :", __func__ );
         fprintf( display, "\n Too many E nodes to be evaluated "
            "in DSC mesh %s !!!", val.name );
         fprintf( display, "\n [ Maximum number is %ld = macro EVLEN "
            "in SOLVER.CONF or %s.", ( long ) EVLEN, __func__ );
         fprintf( display, "\n - Change macro only in compliance "
            "with memory resources !\n " );
         strncpy( solver.fcterr, __func__, SHS_SIZE );
         strncpy( solver.errmsg, solver.val, SHS_SIZE );
         strncat( solver.errmsg, ":too_many_E_nodes_to_be_evaluated_", 
            ( LGS_SIZE - SHS_SIZE ));

         return null;
      }
      else if ( null < val.nen )
      {
         if ( null < val.r )
            evlfld = ONE;
      
         if ( ONE < val.nen )
         {
            fscanf( evalfle, scformat, ptr ); /* mean_value operation mark */

            if ( null != strstr( ptr, "average" ))
               strncpy( val.mode_en, "average", VSS_SIZE );
            else
               strncpy( val.mode_en, "individual", VSS_SIZE );
         };
      };

      fscanf( evalfle, scformat, ptr ); /* string "H-ports" */
      fscanf( evalfle, scformat, ptr );
      val.nhp = strtol( ptr, endp, DEC ); /* number of H-ports evaluated */

      if ( EVLHP < val.nhp )
      {
         fprintf( display, "\n\n Message from function %s :", __func__ );
         fprintf( display, "\n Too many H ports to be evaluated "
            "in DSC mesh %s !!!", val.name );
         fprintf( display, "\n [ Maximum number is %ld = macro EVLHP "
            "in SOLVER.CONF or %s.", ( long ) EVLHP, __func__ );
         fprintf( display, "\n - Change macro only in compliance "
            "with memory resources !\n " );
         strncpy( solver.fcterr, __func__, SHS_SIZE );
         strncpy( solver.errmsg, solver.val, SHS_SIZE );
         strncat( solver.errmsg, ":too_many_H_ports_to_be_evaluated_", 
            ( LGS_SIZE - SHS_SIZE ));

         return null;
      }
      else if ( null < val.nhp )
      {
         if ( null < val.r )
            evlfld = ONE;
      
         if ( ONE < val.nhp )
         {
            fscanf( evalfle, scformat, ptr ); /* mean_value operation mark */

            if ( null != strstr( ptr, "average" ))
               strncpy( val.mode_hp, "average", VSS_SIZE );
            else
               strncpy( val.mode_hp, "individual", VSS_SIZE );
         };
      };

      fscanf( evalfle, scformat, ptr ); /* string "H-nodes" */
      fscanf( evalfle, scformat, ptr );
      val.nhn = strtol( ptr, endp, DEC ); /* number of H-nds to be evaluated */

      if ( EVLHN < val.nhn )
      {
         fprintf( display, "\n\n Message from function %s :", __func__ );
         fprintf( display, "\n Too many H nodes to be evaluated "
            "in DSC mesh %s !!!", val.name );
         fprintf( display, "\n [ Maximum number is %ld = macro EVLHN "
            "in SOLVER.CONF or %s.", ( long ) EVLHN, __func__ );
         fprintf( display, "\n - Change macro only in compliance "
            "with memory resources !\n " );
         strncpy( solver.fcterr, __func__, SHS_SIZE );
         strncpy( solver.errmsg, solver.val, SHS_SIZE );
         strncat( solver.errmsg, ":too_many_H_nodes_to_be_evaluated_", 
            ( LGS_SIZE - SHS_SIZE ));

         return null;
      }
      else if ( null < val.nhn )
      {
         if ( null < val.r )
            evlfld = ONE;
      
         if ( ONE < val.nhn )
         {
            fscanf( evalfle, scformat, ptr ); /* mean_value operation mark */

            if ( null != strstr( ptr, "average" ))
               strncpy( val.mode_hn, "average", VSS_SIZE );
            else
               strncpy( val.mode_hn, "individual", VSS_SIZE );
         };
      };
/*...........................................................................*/
# if DSC_HCRMDE !=0

      evlhcr = null; /* initialize [ set to ONE if there are evaluted ports] */
/*...........................................................................*/
# if DSC_FLDMDE !=0
      fscanf( evalfle, scformat, ptr ); /* string "therm-fluid_evaluation:" */
# else
      fscanf( evalfle, scformat, ptr ); /* string "thermal_evaluation:" */
# endif /* DSC_FLDMDE ... */
/*...........................................................................*/
      fscanf( evalfle, scformat, ptr ); /* string "evaluated_ports" */
      fscanf( evalfle, scformat, ptr ); /* string "number" */

      cc = null; do
      {
         fscanf( evalfle, scformat, ptr ); /* string "heat_current:" */
         fscanf( evalfle, scformat, ptr );
         val.nhc[cc] = strtol( ptr, endp, DEC );

         if ( EVLHC < val.nhc[cc] )
         {
            fprintf( display, "\n\n Message from function %s :", __func__ );
            fprintf( display, "\n Too many heat currents to be evaluated "
               "in DSC mesh %s !!!", val.name );
            fprintf( display, "\n [ Maximum number is %ld = macro EVLHC "
               "in SOLVER.CONF or %s.", ( long ) EVLHC, __func__ );
            fprintf( display, "\n - Change macro only in compliance "
               "with memory resources !\n " );
            strncpy( solver.fcterr, __func__, SHS_SIZE );
            strncpy( solver.errmsg, solver.val, SHS_SIZE );
            strncat( solver.errmsg, \
               ":too_many_heat_currents_to_be_evaluated_",
               ( LGS_SIZE - SHS_SIZE ));

            return null;
         }
         else if ( null < val.nhc[cc] )
         {
            if ( null < val.rc )
               evlhcr = ONE;
         };
/*............................................................................*/
# if DSC_FCTEMP == 1
         fscanf( evalfle, scformat, ptr ); /* "face_temperature" */
         fscanf( evalfle, scformat, ptr );
         val.ntf[cc] = strtol( ptr, endp, DEC );

         if ( EVLTF < val.ntf[cc] )
         {
            fprintf( display, "\n\n Message from function %s :", __func__ );
            fprintf( display, "\n Too many face temperatures to be evaluated "
               "in DSC mesh %s !!!", val.name );
            fprintf( display, "\n [ Maximum number is %ld = macro EVLTF "
               "in SOLVER.CONF or %s.", ( long ) EVLTF, __func__ );
            fprintf( display, "\n - Change macro only in compliance "
               "with memory resources !\n " );
            strncpy( solver.fcterr, __func__, SHS_SIZE );
            strncpy( solver.errmsg, solver.val, SHS_SIZE );
            strncat( solver.errmsg, \
               ":too_many_face_temperatures_to_be_evaluated_",
               ( LGS_SIZE - SHS_SIZE ));

            return null;
         }
         else if ( null < val.ntf[cc] )
         {
            if ( null < val.rc )
               evlhcr = ONE;
         };
# endif /* DSC_FCTEMP == 1 */
/*............................................................................*/
         fscanf( evalfle, scformat, ptr ); /* "node_temperature" */
         fscanf( evalfle, scformat, ptr );
         val.ntn[cc] = strtol( ptr, endp, DEC );

         if ( EVLTN < val.ntn[cc] )
         {
            fprintf( display, "\n\n Message from function %s :", __func__ );
            fprintf( display, "\n Too many temperature nodes to be evaluated "
               "in DSC mesh %s !!!", val.name );
            fprintf( display, "\n [ Maximum number is %ld = macro EVLTN "
               "in SOLVER.CONF or %s.", ( long ) EVLTN, __func__ );
            fprintf( display, "\n - Change macro only in compliance "
               "with memory resources !\n " );
            strncpy( solver.fcterr, __func__, SHS_SIZE );
            strncpy( solver.errmsg, solver.val, SHS_SIZE );
            strncat( solver.errmsg, ":too_many_H_nodes_to_be_evaluated_", 
               ( LGS_SIZE - SHS_SIZE ));

            return null;
         }
         else if ( null < val.ntn[cc] )
         {
            if ( null < val.rc )
               evlhcr = ONE;
         };
/*............................................................................*/
# if DSC_FLDMDE != 0
/* store nodal velocities: */

         fscanf( evalfle, scformat, ptr ); /* "nodal_velocity" */
         fscanf( evalfle, scformat, ptr );
         val.nun[cc] = strtol( ptr, endp, DEC );

         if ( EVLUN < val.nun[cc] )
         {
            fprintf( display, "\n\n Message from function %s :", __func__ );
            fprintf( display, "\n Too many nodal velocities to be evaluated "
               "in DSC mesh %s !!!", val.name );
            fprintf( display, "\n [ Maximum number is %ld = macro EVLUN "
               "in SOLVER.CONF or %s.", ( long ) EVLUN, __func__ );
            fprintf( display, "\n - Change macro only in compliance "
               "with memory resources !\n " );
            strncpy( solver.fcterr, __func__, SHS_SIZE );
            strncpy( solver.errmsg, solver.val, SHS_SIZE );
            strncat( solver.errmsg, ":too_many_H_nodes_to_be_evaluated_", 
               ( LGS_SIZE - SHS_SIZE ));

            return null;
         }
         else if ( null < val.nun[cc] )
         {
            if ( null < val.rc )
               evlhcr = ONE;
         };
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
      } while (( ++cc ) < DSC_HCRMDE );
# endif
/*............................................................................*/
      fscanf( evalfle, scformat, ptr); /* string "...LABELS_etc....", e.g. */

      for ( ii=null; ii<val.nep; ii++ ) /* enter E-ports */
      {
         fscanf( evalfle, scformat, ptr );
         val.mep[ii] = strtol( ptr, endp, DEC ); /* mesh cell index */
         fscanf( evalfle, scformat, ptr );
         val.pep[ii] = strtol( ptr, endp, DEC ); /* port index */
      };
       
      for ( ii=null; ii<val.nen; ii++ ) /* enter E-nodes */
      {
         fscanf( evalfle, scformat, ptr );
         val.men[ii] = strtol( ptr, endp, DEC ); /* mesh cell index */
         fscanf( evalfle, scformat, ptr );
         val.cen[ii] =\
            strtol( ptr, endp, DEC ); /* component index [0:u, 1:v, 2:w] */

         if ( *ptr == 'u' ) /* correct erroneous indices: */
            val.cen[ii] = null;
         if ( *ptr == 'v' ) 
            val.cen[ii] = ONE;
         if ( *ptr == 'w' ) 
            val.cen[ii] = TWO;
      };

      for ( ii=null; ii<val.nhp; ii++ ) /* enter H-ports */
      {
         fscanf( evalfle, scformat, ptr );
         val.mhp[ii] = strtol( ptr, endp, DEC ); /* mesh cell index */
         fscanf( evalfle, scformat, ptr );
         val.php[ii] = strtol( ptr, endp, DEC ); /* port index */
      };

      for ( ii=null; ii<val.nhn; ii++ ) /* enter H-nodes */
      {
         fscanf( evalfle, scformat, ptr );
         val.mhn[ii] = strtol( ptr, endp, DEC ); /* mesh cell index */
         fscanf( evalfle, scformat, ptr );
         val.chn[ii] =\
            strtol( ptr, endp, DEC ); /* component index [0:u, 1:v, 2:w] */

         if ( *ptr == 'u' ) /* correct erroneous indices: */
            val.chn[ii] = null;
         if ( *ptr == 'v' )
            val.chn[ii] = ONE;
         if ( *ptr == 'w' )
            val.chn[ii] = TWO;
      };
/*............................................................................*/
# if DSC_HCRMDE != 0
      cc = null; do
      {
         for ( ii=null; ii<val.nhc[cc]; ii++ ) /* enter heat current faces */
         {
            fscanf( evalfle, scformat, ptr );
            val.mhc[cc][ii] = strtol( ptr, endp, DEC ); /* mesh cell index */
            fscanf( evalfle, scformat, ptr );
            val.fhc[cc][ii] = strtol( ptr, endp, DEC ); /* face index */
         };
/*............................................................................*/
# if DSC_FCTEMP == 1
         for ( ii=null; ii<val.ntf[cc]; ii++ ) /* enter temperature faces */
         {
            fscanf( evalfle, scformat, ptr );
            val.mtf[cc][ii] = strtol( ptr, endp, DEC ); /* mesh cell index */
            fscanf( evalfle, scformat, ptr );
            val.ftf[cc][ii] = strtol( ptr, endp, DEC ); /* face index */
         };
# endif /* DSC_FCTEMP == 1 */
/*............................................................................*/
         for ( ii=null; ii<val.ntn[cc]; ii++ ) /* enter temperature nodes */
         {
            fscanf( evalfle, scformat, ptr );
            val.mtn[cc][ii] = strtol( ptr, endp, DEC ); /* mesh cell index */
         };
/*............................................................................*/
# if DSC_FLDMDE != 0
         for ( ii=null; ii<val.nun[cc]; ii++ ) /* enter node velocity indices */
         {                                     
            fscanf( evalfle, scformat, ptr );
            val.mun[cc][ii] = strtol( ptr, endp, DEC ); /* mesh cell index */
            fscanf( evalfle, scformat, ptr );
            val.cun[cc][ii] =\
               strtol( ptr, endp, DEC ); /* component index [0:x, 1:y, 2:z] */

            if ( *ptr == 'x' ) /* correct erroneous indices: */
               val.cun[cc][ii] = null;
            if ( *ptr == 'y' )
               val.cun[cc][ii] = ONE;
            if ( *ptr == 'z' )
               val.cun[cc][ii] = TWO;
         };
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
      } while (( ++cc ) < DSC_HCRMDE );

# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
      fscanf( evalfle, scformat, ptr );

      if ( null != strncmp( ptr, "DSC", THREE )) /* check string *clsmsg = "DSC" */
      {
         printf( "\n\n Message from function %s :", __func__ );
         printf( "\n Illegal format of file dsc.val%d !!!", joblbl[jj] );
         printf( "\n [ Cannot read closure message 'DSC ...' ] " );

         fprintf( evalfle, "\nJob skipped.\n\n" );

         fclose( evalfle );

         fprintf( display, "\n DSC job no.%d skipped.\n", joblbl[jj] );

         return null;
      };

      for ( ii=null; ii<10; ii++ ) /* file ... created ... */
      {
         fscanf( evalfle, scformat, ptr );
      };

      fleofs = ftell( evalfle );
      fseek( evalfle, fleofs, SEEK_SET );

/* stop time: */
      nseconds = time( timer );

      strncpy( ctmptr, ctime( &nseconds ), 24 );
/* same functionality as ' ctmptr = asctime( localtime( &nseconds )); '*/
/*............................................................................*/
/* eventually, print error message on evaluation file: */

      if ( null != strncmp( solver.errmsg, astrx, THREE ))
      {
         fprintf( evalfle, "\n\n%s%d,\n%s", "Abnormal end"
                           " of DSC job no.", joblbl[jj], ctmptr );

         fprintf( evalfle, "\nError message from %s :\n\n", solver.fcterr );

         ii = null; do
         {
            fprintf( evalfle, "%s", solver.errmsg );
         }  while(( ++ii ) < THREE );

         fprintf( evalfle, "\nJob skipped.\n\n" );

         fclose( evalfle );

         fprintf( display, "\n DSC job no.%d skipped.\n", joblbl[jj] );

         return null;
      };
/*............................................................................*/
/* The following macro gives a conveniently readable form to string <ctmptr> */
/* which is returned as the new string <tmestr> */

      TIMEFORM( tmestr, ctmptr );

      fprintf( display, "\n opened: evaluation file %s", solver.val );

      fprintf( evalfle,  "\n______________________________"
         "________________________________________________\n" );
      fprintf( evalfle, "%s %s%d,\n%s", "output of program SOLVER.C,",
         "job no", joblbl[jj], tmestr );
      fprintf( evalfle,  "\n______________________________"
         "________________________________________________\n" );

      if ( bnd.p <= null ) 
      {
         val.read = 1;
         fprintf( evalfle, spformat, "APERIODIC_STRUCTURE" );
      }
      else
      {
         val.read = 2;
         fprintf( evalfle, spformat, "PERIODIC_STRUCTURE_" );
      };

      fprintf( evalfle, spformat, "Maxwell_field_computation:" );

      if ( solver.dmn == 't' ) /* time domain */
      {
         fprintf( evalfle, spformat, "TIME_DOMAIN________" );
         fprintf( evalfle, "[ Timestep = %.16e seconds ]\n",
            solver.dt ); /* DSC time step */
      }
      else /* frequency domain */
      {
         val.read = 2;
         fprintf( evalfle, spformat, "FREQUENCY_DOMAIN___" );
         fprintf( evalfle, "[ Frequency = %.16e Hertz ]\n",
            val.fr ); /* frequency */
      };

      fprintf( evalfle, spformat, "Maxwell_field_excitation:" );
      fprintf( evalfle, spformat, exc.type ); /* EMfield excitation type */
/*............................................................................*/
# if DSC_HCRMDE != 0
# if DSC_FLDMDE == 0
      fprintf( evalfle, spformat, "Thermal_computation:" );
      fprintf( evalfle, spformat, "TIME_DOMAIN" );
      fprintf( evalfle, "[ Timestep = %.16e seconds ]\n",
         solver.hcdt ); /* thermal time step */
      fprintf( evalfle, spformat, "Thermal_excitation:" );
      fprintf( evalfle, spformat, exc.hctp ); /* excitation type */
# else /* if DSC_FLDMDE != 0 */
      fprintf( evalfle, spformat, "Thermal-fluid_computation:" );
      fprintf( evalfle, spformat, "TIME_DOMAIN" );
      fprintf( evalfle, "[ Timestep = %.16e seconds ]\n",
         solver.hcdt ); /* thermal time step */
      fprintf( evalfle, spformat, "Thermal-fluid_excitation:" );
      fprintf( evalfle, spformat, exc.hctp ); /* excitation type */
# endif /* DSC_FLDMDE != 0 */
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
      fprintf( evalfle, spformat, top.name ); /* DSC system identifier */
      fprintf( evalfle, spformat, top.text ); /* text, comment ...     */

      if ( solver.dmn == 't' ) /* solver.dmn == 't'ime domain */
      {
         fprintf( evalfle, spformat, "seconds" );

         if ((( val.nhp == null )&&(val.nhn == null )) 
           &&(( null < val.nep )||( null < val.nen )))
            fprintf( evalfle, spformat, "volts" );
         else if ((( val.nep == null )&&(val.nen == null ))
                &&(( null < val.nhp )||( null < val.nhn )))
            fprintf( evalfle,spformat, "amperes" );
         else if ((( null < val.nep )||( null < val.nen ))
                &&(( null < val.nhp )||( null < val.nhn )))
            fprintf( evalfle,spformat, "[volts,amperes]" );
         else
            fprintf( evalfle,spformat, "degree_C" );

         fprintf( evalfle, dpformat, ( val.ni*val.dt ));
         fprintf( evalfle, dpformat, ( val.n*val.dt ));
         fprintf( evalfle, dpformat, val.dt );
         fprintf( evalfle, lpformat, ( val.n - val.ni + ONE ));
      } /* end if solver.dmn == 't'ime_domain */
      else /* if solver.dmn == 'f'requency_domain */
      {
         fprintf( evalfle, spformat, "cycles" );

         if ((( val.nhp == null )&&(val.nhn == null )) 
           &&(( null < val.nep )||( null < val.nen )))
            fprintf( evalfle, spformat, "volts" );
         else if ((( val.nep == null )&&(val.nen == null ))
                &&(( null < val.nhp )||( null < val.nhn )))
            fprintf( evalfle,spformat, "amperes" );
         else if ((( null < val.nep )||( null < val.nen ))
                &&(( null < val.nhp )||( null < val.nhn )))
            fprintf( evalfle,spformat, "[volts,amperes]" );
         else
            fprintf( evalfle,spformat, "degree_C" );
/*............................................................................*/
# if DSC_DOMAIN == 3
         mm = val.n - val.ni;
         nn = val.ni + ( mm*val.r );

         fprintf( evalfle, lpformat, val.ni );
         fprintf( evalfle, lpformat, nn );
         fprintf( evalfle, lpformat, ( long ) val.r );
         fprintf( evalfle, lpformat, ( mm + ONE ));
# else
         mm = val.r * val.ni;
         nn = val.r * val.n;

         fprintf( evalfle, lpformat, mm );
         fprintf( evalfle, lpformat, nn );
         fprintf( evalfle, lpformat, ( long ) val.r );
         fprintf( evalfle, lpformat, ( val.n - val.ni + ONE ));
# endif /* DSC_DOMAIN == ... */
/*............................................................................*/
      };

      fflush( evalfle );
/*............................................................................*/
# if VAL_DISP ==1
      fprintf( display, "\n\n number of iteration cycles ............: %ld",
         val.n );
      if ( ONE < val.r )
         fprintf( display, "\n internal repetition rate [Maxwell fld].: %d",
            val.r );
/*............................................................................*/
# if DSC_HCRMDE != 0
# if DSC_FLDMDE == 0
      if ( ONE < val.rc )
         fprintf( display, "\n internal repetition rate [thermal].....: %d",
            val.rc );
# else
      if ( ONE < val.rc )
         fprintf( display, "\n internal repetition rate [therm-fluid].: %d",
            val.rc );
# endif /* DSC_FLDMDE != 0 */
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
      if ( null < val.ni )
         fprintf( display, "\n first evaluated cycle .................: %ld",
            val.ni );
      if ( solver.dmn == 't' ) /* solver.dmn == 't'ime_domain */ 
      {
/*............................................................................*/
# if DSC_HCRMDE == 0
         fprintf( display, "\n computed time interval.................: "
            "[ %.5e , %.5e ] seconds", ZERO, ( val.n*val.dt ));
# else
         fprintf( display, "\n computed time interval.[Maxwell fld]...: "
            "[ %.5e , %.5e ] seconds", ZERO, ( val.n*val.dt ));
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
         if ( ONE < val.r )
         {
/*............................................................................*/
# if DSC_HCRMDE == 0
            fprintf( display, "\n internal time step.....................: "
               "%.16e seconds", solver.dt );
            fprintf( display, "\n cycle time.............................: "
               "%.16e seconds", val.dt );
# else
            fprintf( display, "\n internal time step [Maxwell field].....: "
               "%.16e seconds", solver.dt );
            fprintf( display, "\n cycle time [Maxwell field].............: "
               "%.16e seconds", val.dt );
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
         }
         else if ( null < val.r )
         {
/*............................................................................*/
# if DSC_HCRMDE == 0
            fprintf( display, "\n time step..............................: "
               "%.16e seconds", val.dt ); 
# else
            fprintf( display, "\n time step [Maxwell field]..............: "
               "%.16e seconds", val.dt );
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
         };
/*............................................................................*/
# if DSC_HCRMDE != 0
         if ( ONE < val.rc )
         {
# if DSC_FLDMDE == 0
            fprintf( display, "\n internal time step [thermal]...........: "
               "%.16e seconds", solver.hcdt );
            fprintf( display, "\n cycle time [thermal]...................: "
               "%.16e seconds", ( solver.hcdt*val.rc ));
# else
            fprintf( display, "\n internal time step [thermal-fluid].....: "
               "%.16e seconds", solver.hcdt );
            fprintf( display, "\n cycle time [thermal-fluid].............: "
               "%.16e seconds", ( solver.hcdt*val.rc ));
# endif /* DSC_FLDMDE != 0 */
         }
         else if ( null < val.rc )
         {
# if DSC_FLDMDE == 0
            fprintf( display, "\n time step [thermal]....................: "
               "%.16e seconds", solver.hcdt ); 
# else /* if DSC_FLDMDE != 0 */
            fprintf( display, "\n time step [thermal-fluid]..............: "
               "%.16e seconds", solver.hcdt ); 
# endif /* DSC_FLDMDE != 0 */
         };
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
      } /* end if solver.dmn == 't'ime_domain */
      else /* if solver.dmn == 'f'requency_domain */
      {
         fprintf( display, "\n frequency..............................: "
            "%.16e Hertz", val.fr );
         fprintf( display, "\n phase shift............................: "
            "%.16e radians", val.ph );
/*............................................................................*/
# if DSC_HCRMDE != 0
# if DSC_FLDMDE == 0
         if ( ONE < val.rc )
         {
            fprintf( display, "\n internal time step [thermal]...........: "
               "%.16e seconds", solver.hcdt );
            fprintf( display, "\n cycle time [thermal]...................: "
               "%.16e seconds", ( solver.hcdt*val.rc ));
         }
         else if ( null < val.rc )
         {
            fprintf( display, "\n time step [thermal]....................: "
               "%.16e seconds", solver.hcdt ); 
         };
# else /* if DSC_FLDMDE != 0 */
         if ( ONE < val.rc )
         {
            fprintf( display, "\n internal time step [thermal-fluid].....: "
               "%.16e seconds", solver.hcdt );
            fprintf( display, "\n cycle time [thermal-fluid].............: "
               "%.16e seconds", ( solver.hcdt*val.rc ));
         }
         else if ( null < val.rc )
         {
            fprintf( display, "\n time step [thermal-fluid]..............: "
               "%.16e seconds", solver.hcdt ); 
         };
# endif /* DSC_FLDMDE != 0 */
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
      }; /* end if solver.dmn == 'f'requency_domain */
/*............................................................................*/
/* E_field evaluation mode */

      if ( null < val.nep )
      {
         if (( ONE < val.nep )&&( *val.mode_ep == 'a' ))
         {
            fprintf( display, "\n\n E-field averaging over %05d ports: "
               ">>--------------------------------------->", val.nep );
         }
         else 
         {
            fprintf( display, "\n\n E-field evaluation at %05d ports: "
               ">>---------------------------------------->", val.nep );
         };
         llns = val.nep/clns + ONE;

         for ( kk=null; kk<llns; kk++ )
         {
            fprintf( display, "\n -> cell:" );
            pp = kk*clns;

            for ( ii=ONE; ii<=clns; ii++ )
            {
               if (( lbnd < kk )&&( cbnd < ii ))
               {
                  fprintf( display, "   ***%c", bslsh );
                  fprintf( display, "%6ld%c", val.mep[( val.nep-ONE )], bslsh );
                  goto E_ports;
               };

               if ( pp < val.nep )
                  fprintf( display, "%6ld%c", val.mep[pp], bslsh );

               if ( val.nep <= ( ++pp ))  
                  goto E_ports; 
            };

           E_ports:

            fprintf( display, "\n -> port:" );
            pp = kk*clns;

            for ( ii=ONE; ii<=clns; ii++ )
            {
               if (( lbnd < kk )&&( cbnd < ii ))
               {
                  fprintf( display, "   ***%c", slsh );
                  fprintf( display, "%6ld%c", 
                    ( long ) val.pep[( val.nep-ONE )], slsh );
                  goto E_nodes;
               };

               if ( pp < val.nep )   
                  fprintf( display, "%6ld%c", ( long ) val.pep[pp], slsh );

               if ( val.nep <= ( ++pp ))  
                  goto E_nodes;
            };
         };
      };

     E_nodes:

      if ( null < val.nen )
      {
         if (( ONE < val.nen )&&( *val.mode_en == 'a' ))
         {
            fprintf( display, "\n\n E-field averaging over %05d nodes: "
               ">>--------------------------------------->", val.nen );
         }
         else
         {
            fprintf( display, "\n\n E-field evaluation at %05d nodes: "
               ">>---------------------------------------->", val.nen );
         };
         llns = val.nen/clns + ONE;

         for ( kk=null; kk<llns; kk++ )
         {
            fprintf( display, "\n -> node:" );
            pp = kk*clns;

            for ( ii=ONE; ii<=clns; ii++ )
            {
               if (( lbnd < kk )&&( cbnd < ii ))
               {
                  fprintf( display, "   ***%c", bslsh );
                  fprintf( display, "%6ld%c", val.men[( val.nen-ONE )], bslsh );
                  goto E_compts;
               };

               if ( pp < val.nen )   
                  fprintf( display, "%6ld%c", val.men[pp], bslsh ); 

               if ( val.nen <= ( ++pp ))  
                  goto E_compts; 
            };

           E_compts:  

            fprintf( display, "\n -> vect:" );
            pp = kk*clns;

            for ( ii=ONE; ii<=clns; ii++ )
            {
               if (( lbnd < kk )&&( cbnd < ii ))
               {
                  fprintf( display, "   ***%c", slsh );
                  fprintf( display, "%6c%c",
                     val.cen[( val.nen-ONE) ]+117, slsh );

                  goto H_field;
               };
               if ( pp < val.nen )   
                  fprintf( display, "%6c%c", val.cen[pp]+117, slsh );

               if ( val.nen <= ( ++pp ))  
                  goto H_field; 
            };
         };
      };
/*............................................................................*/
/* H_field evaluation mode */

     H_field:

      if ( null < val.nhp )
      {
         if (( ONE < val.nhp )&&( *val.mode_hp == 'a' ))
         {
            fprintf( display, "\n\n H-field averaging over %05d ports: "
               ">>--------------------------------------->", val.nhp );
         }
         else
         {
            fprintf( display, "\n\n H-field evaluation at %05d ports: "
               ">>---------------------------------------->", val.nhp );
         };
         llns = val.nhp/clns + ONE;

         for ( kk=null; kk<llns; kk++ )
         {
            fprintf( display, "\n -> cell:" );
            pp = kk*clns;

            for ( ii=ONE; ii<=clns; ii++ )
            {
               if (( lbnd < kk )&&( cbnd < ii ))
               {
                  fprintf( display, "   ***%c", bslsh );
                  fprintf( display, "%6ld%c", val.mhp[( val.nhp-ONE )], bslsh );
                  goto H_ports;
               };

               if ( pp < val.nhp )   
                  fprintf( display, "%6ld%c", val.mhp[pp], bslsh );

               if ( val.nhp <= ( ++pp ))  
                  goto H_ports; 
            };

           H_ports:

            fprintf( display, "\n -> port:" );
            pp = kk*clns;

            for ( ii=ONE; ii<=clns; ii++ )
            {
               if (( lbnd < kk )&&( cbnd < ii ))
               {
                  fprintf( display, "   ***%c", slsh );
                  fprintf( display, "%6ld%c",
                    ( long ) val.php[( val.nhp-ONE )], slsh );
                  goto H_nodes; 
               };

               if ( pp < val.nhp )   
                  fprintf( display, "%6ld%c", ( long ) val.php[pp], slsh );

               if ( val.nhp <= ( ++pp ))  
                  goto H_nodes;  
            };
         };
      };

     H_nodes:

      if ( null < val.nhn )
      {
         if (( ONE < val.nhn )&&( *val.mode_hn == 'a' ))
         {
            fprintf( display, "\n\n H-field averaging over %05d nodes: "
               ">>--------------------------------------->", val.nhn );
         }
         else
         {
            fprintf( display, "\n\n H-field evaluation at %05d nodes: "
               ">>---------------------------------------->", val.nhn );
         };
         llns = val.nhn/clns + ONE;

         for ( kk=null; kk<llns; kk++ )
         {
            fprintf( display, "\n -> node:" );
            pp = kk*clns;

            for ( ii=ONE; ii<=clns; ii++ )
            {
               if (( lbnd < kk )&&( cbnd < ii ))
               {
                  fprintf( display, "   ***%c", bslsh );
                  fprintf( display, "%6ld%c", val.mhn[( val.nhn-ONE )], bslsh );
                  goto H_compts;
               };

               if ( pp < val.nhn )   
                  fprintf( display, "%6ld%c", val.mhn[pp], bslsh ); 

               if ( val.nhn <= ( ++pp ))  
                  goto H_compts; 
            };

           H_compts:

            fprintf( display, "\n -> vect:" );
            pp = kk*clns;

            for ( ii=ONE; ii<=clns; ii++ )
            {
               if (( lbnd < kk )&&( cbnd < ii ))
               {
                  fprintf( display, "   ***%c", slsh );
                  fprintf( display, "%6c%c",
                     val.chn[( val.nhn-ONE )]+117, slsh );

                  goto Maxwfld_listo;
               };
               if ( pp < val.nhn )   
                  fprintf( display, "%6c%c", val.chn[pp]+117, slsh ); 

               if ( val.nhn <= ( ++pp ))  
                  goto Maxwfld_listo; 
            };
         };
      };

     Maxwfld_listo:
/*............................................................................*/
# if DSC_HCRMDE != 0

      cc = null; do
      {
/*............................................................................*/
/* display evaluated heat currents: */

         if ( null < val.nhc[cc] )
         {
            fprintf( display, "\n\n heat current evaluation at %05d faces: "
               ">>------------------------------------>", val.nhc[cc] );

            llns = ( val.nhc[cc] /clns ) + ONE;

            for ( kk=null; kk<llns; kk++ )
            {
               fprintf( display, "\n -> cell:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )&&( cbnd < ii ))
                  {
                     fprintf( display, "   ***%c", bslsh );
                     fprintf( display, "%6ld%c",
                       ( val.mhc[cc][( val.nhc[cc]-ONE )] ), bslsh );
                     goto HC_faces;
                  };

                  if ( pp < val.nhc[cc] )
                     fprintf( display, "%6ld%c", ( val.mhc[cc][pp] ), bslsh );

                  if ( val.nhc[cc] <= ( ++pp ))
                     goto HC_faces; 
               };

              HC_faces:

               fprintf( display, "\n -> face:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )&&( cbnd < ii ))
                  {
                     fprintf( display, "   ***%c", slsh );
                     fprintf( display, "%6ld%c",
                        ( long ) ( val.fhc[cc][( val.nhc[cc]-ONE )] ), slsh );

                     goto TF_cells; 
                  };

                  if ( pp < val.nhc[cc] )
                     fprintf( display, "%6ld%c",
                        ( long ) ( val.fhc[cc][pp] ), slsh );

                  if ( val.nhc[cc] <= ( ++pp ))  
                     goto TF_cells;

               }; /* end for ( ii=ONE; ii<=clns; ii++ ) */
            }; /* end for ( kk=null; kk<llns; kk++ ) */
         }; /* end if ( null < val.nhc[cc] ) */

        TF_cells:
/*............................................................................*/
# if DSC_FCTEMP == 1
/* display evaluated temperature faces: */

         if ( null < val.ntf[cc] )
         {
            fprintf( display, "\n\n temperature evaluation at %05d faces: "
               ">>------------------------------------>", val.ntf[cc] );

            llns = ( val.ntf[cc] /clns ) + ONE;

            for ( kk=null; kk<llns; kk++ )
            {
               fprintf( display, "\n -> cell:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )&&( cbnd < ii ))
                  {
                     fprintf( display, "   ***%c", bslsh );
                     fprintf( display, "%6ld%c",
                       ( long )( val.mtf[cc][( val.ntf[cc]-ONE )] ), bslsh );
                     goto TF_faces;
                  };

                  if ( pp < val.ntf[cc] )
                     fprintf( display, "%6ld%c",
                        ( long )( val.mtf[cc][pp] ), bslsh );

                  if ( val.ntf[cc] <= ( ++pp ))
                     goto TF_faces; 
               };

              TF_faces:

               fprintf( display, "\n -> face:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )&&( cbnd < ii ))
                  {
                     fprintf( display, "   ***%c", slsh );
                     fprintf( display, "%6ld%c",
                        ( long )( val.ftf[cc][( val.ntf[cc]-ONE )] ), slsh );

                     goto TN_nodes; 
                  };

                  if ( pp < val.ntf[cc] )
                     fprintf( display, "%6ld%c",
                        ( long )( val.ftf[cc][pp] ), slsh );

                  if ( val.ntf[cc] <= ( ++pp ))  
                     goto TN_nodes;

               }; /* end for ( ii=ONE; ii<=clns; ii++ ) */
            }; /* end for ( kk=null; kk<llns; kk++ ) */
         }; /* end if ( null < val.ntf[cc] ) */

        TN_nodes:

# endif /* DSC_FCETEMP == 1 */
/*............................................................................*/
/* display evaluated temperature nodes: */

         if ( null < val.ntn[cc] )
         {
            fprintf( display, "\n\n temperature evaluation at %05d nodes: "
               ">>------------------------------------>", val.ntn[cc] );

            llns = ( val.ntn[cc] /clns ) + ONE;

            for ( kk=null; kk<llns; kk++ )
            {
               fprintf( display, "\n -> cell:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )&&( cbnd < ii ))
                  {
                     fprintf( display, "   ***%c", bslsh );
                     fprintf( display, "%6ld%c",
                       ( long )( val.mtn[cc][( val.ntn[cc]-ONE )] ), bslsh );
                     goto TN_ports;
                  };

                  if ( pp < val.ntn[cc] )   
                     fprintf( display, "%6ld%c",
                        ( long )( val.mtn[cc][pp] ), bslsh );

                  if ( val.ntn[cc] <= ( ++pp ))
                     goto TN_ports; 
               };

              TN_ports:

               fprintf( display, "\n -> port:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )&&( cbnd < ii ))
                  {
                     fprintf( display, "   ***%c", slsh );
                     fprintf( display, "%6c%c", '-', slsh );

                     goto hcrr_listo; 
                  };

                  if ( pp < val.ntn[cc] )
                     fprintf( display, "%6c%c", '-', slsh );

                  if ( val.ntn[cc] <= ( ++pp ))  
                     goto hcrr_listo;

               }; /* end for ( ii=ONE; ii<=clns; ii++ ) */
            }; /* end for ( kk=null; kk<llns; kk++ ) */
         }; /* end if ( null < val.ntn[cc] ) */

        hcrr_listo: ;
/*............................................................................*/
# if DSC_FLDMDE != 0
/* display evaluated nodal velocities: */

         if ( null < val.nun[cc] )
         {
            fprintf( display, "\n\n fluid velocity evaluation at %05d nodes: "
               ">>--------------------------------->", val.nun[cc] );

            llns = ( val.nun[cc] /clns ) + ONE;

            for ( kk=null; kk<llns; kk++ )
            {
               fprintf( display, "\n -> cell:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )&&( cbnd < ii ))
                  {
                     fprintf( display, "   ***%c", bslsh );
                     fprintf( display, "%6ld%c",
                       ( long )( val.mun[cc][( val.nun[cc]-ONE )] ), bslsh );
                     goto UN_compts;
                  };

                  if ( pp < val.nun[cc] )
                     fprintf( display, "%6ld%c",
                        ( long )( val.mun[cc][pp] ), bslsh );

                  if ( val.nun[cc] <= ( ++pp ))
                     goto UN_compts; 
               };

              UN_compts:

               fprintf( display, "\n -> comp:" );
               pp = kk*clns;

               for ( ii=ONE; ii<=clns; ii++ )
               {
                  if (( lbnd < kk )&&( cbnd < ii ))
                  {
                     fprintf( display, "   ***%c", slsh );
                     fprintf( display, "%6c%c",
                        (( val.cun[cc][( val.nun[cc]-ONE )] )+120 ), slsh );

                     goto fluids_listo;
                  };

                  if ( pp < val.nun[cc] )
                     fprintf( display, "%6c%c",
                        (( val.cun[cc][pp] )+120 ), slsh );

                  if ( val.nun[cc] <= ( ++pp ))  
                     goto fluids_listo;

               }; /* end for ( ii=ONE; ii<=clns; ii++ ) */
            }; /* end for ( kk=null; kk<llns; kk++ ) */
         }; /* end if ( null < val.nun[cc] ) */

        fluids_listo: ;

# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
      } while(( ++cc ) < DSC_HCRMDE );
# endif /* DSC_HCRMDE != 0 */
/*............................................................................*/
      fprintf( display, "\n --------------------------------------"
               "----------------------------------------\n");
# endif

      return ONE;
   }
/*............................................................................*/
   else if ( jj == ( - TWO ))
   { 
/*............................................................................*/
      fprintf( evalfle, "\n'" );
/*...........................................................................*/
# if DSC_DOMAIN == 0
      if ( null < dfl.emx )
      {
         fprintf( evalfle, "Maximum gyroelectric current density: "
            "%+.12e V/m^2\n ", dfl.emx );
      };
      if ( null < dfl.mmx ) 
      {
         fprintf( evalfle, "Maximum gyromagnetic current density: "
            "%+.12e V/m^2\n ", dfl.mmx );
      };
# elif DSC_DOMAIN == 1
      if ( null < dfl.emx )
      {
         fprintf( evalfle, "Maximum gyroelectric current density: "
            "%+.12e A/m^2\n ", dfl.emx );
      };
      if ( null < dfl.mmx ) 
      {
         fprintf( evalfle, "Maximum gyromagnetic current density: "
            "%+.12e A/m^2\n ", dfl.mmx );
      };
# endif
/*............................................................................*/
# if (( DSC_HCRMDE != 0 )&&( VAL_HSCSUM != 0 ))
/* integrate heat sources over DSC mesh, then display total heat generated */
/* on final time step and pertinent mesh volume [ where heat sources are ] */

      hsp = &hcs;

      cc = null; do
      {
         if (( null != val.nhc[cc] )
           ||( null != val.ntf[cc] )
           ||( null != val.ntn[cc] ))
         {   

# if (( VAL_HSCSUM == 1 )||( VAL_HSCSUM == 2 ))
/* sum of all skin effect heat sources: */

            solver.hsk = ZERO;
            mm = null;
            while( mm < bnd.nsk )
            {
               solver.hsk += bnd.hs[mm];
               mm++;
            }; /* end while ( mm < bnd.nsk ) */
# endif /* (( VAL_HSCSUM == 1 )||( VAL_HSCSUM == 2 )) */

# if (( VAL_HSCSUM == 2 )||( VAL_HSCSUM == 3 ))
/* inegrate mesh cell heat sources and pertinent volume: */

            hci = &hcr[cc];
            ( hci->qn[null] ) = ZERO;

# if VAL_HSCVOL == 1
            vol = ZERO;
# endif
            mm = null;
            while(( mm++ ) < ( tpt->n ))
            {
               nn = ( hsp->hh[mm] );
               if (( hsp->scs[nn] ) != null )
               {
                  ( hci->qn[null] ) += (( hci->qn[mm] )*( hsp->vol[nn] ));

# if VAL_HSCVOL == 1
                  vol += ( hsp->vol[nn] );
# endif
               };
            }; /* end while ( mm < ( tpt->n )) */
# endif /* (( VAL_HSCSUM == 2 )||( VAL_HSCSUM == 3 )) */
/*............................................................................*/
/* write results into file dsc.val<*>: */

            fprintf( evalfle, "%d.heat generated at final time step:\n ",
               ( cc+ONE ));

# if (( VAL_HSCSUM == 1 )\
    ||( VAL_HSCSUM == 2 ))
            fprintf( evalfle, "skin effect:      %.12e Watts\n ",
               solver.hsk );
# endif
# if (( VAL_HSCSUM == 2 )\
    ||( VAL_HSCSUM == 3 ))
            fprintf( evalfle, "internal sources: %.12e Watts\n ",
               ( hci->qn[null] ));
# if VAL_HSCVOL == 1
            fprintf( evalfle, "[The internal sources are distributed over "
               "the volume: %.12e m^3]\n ", vol );
# endif
# endif /* (( VAL_HSCSUM == 2 )||( VAL_HSCSUM == 3 )) */
            fprintf( evalfle, "\n " );
         }; /* end if (( null != val.nhc[cc] )||...*/
      } while (( ++cc ) < DSC_HCRMDE );
# endif /* (( DSC_HCRMDE != 0 )&&( VAL_HSCSUM != 0 )) */
/*............................................................................*/
/* closure message on file dsc.val<*>: */

      nseconds = time( timer );

/* The following macro gives a conveniently readable form to string <ctmptr> */
/* which is returned as the new string <tmestr> */

      TIMEFORM( tmestr, ctmptr );

      fprintf( evalfle, "%s\n %s'", "DSC job terminated:", tmestr );

      fflush( evalfle );
      fclose( evalfle );

      if( *solver.opt == null )
      {
         fprintf( display, "\r DSC process terminated, "
            "results on file %-38s\n", solver.val );
      }
      else
      {
         fprintf( display, "\n DSC process terminated, "
            "results on file %s", solver.val );
      };

      return ONE;
   };
   return null;
}
/*============================================================================*/
# undef VAL_INITLBL
# undef VAL_DISP
/*----------------------------------------------------------------------------*/
# if DSC_HCRMDE != 0
   # undef VAL_HSCVOL
   # undef VAL_HSCSUM
# endif
/************** end of DSC field evaluation function values(*) ****************/
