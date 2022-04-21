/* [ file: inputs.h ] */
# define DO_INPUT "input(*)"
/*******************************************************************************
*                                                                              *
*   DSC model parameter input function                                         *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0r1 ]                                                          *
*                                                                              *
*   In  option  'd'ivisions  this function  enters the discrete DSC model      *
*   parameters on calling function deflt_divsns( ) , or from file div.input    *
*   [ formatted as part `DIVISIONS' of files div.log_N or par.log_N ].         *
*   These parameters may be revised in function rvise_divsns(*), on head of    *
*   the DSC model file model.h .                                               *
*                                                                              *
*   In option `p'arameters the present function enters the continuous DSC      *
*   model parameters on calling function deflt_params( ), or from a file       *
*   `par.input' [ formatted as part `PARAMETERS' of file par.log_N ].          *
*   Again, these parameters may be revised in function rvise_params( ).        *
*                                                                              *
*   (C) SHEIN; Munich, April 2007                             Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define MAX_LABEL  500
/* [ parser command: "more -e", e.g. ] */
# ifndef DSC_PAGER
   # define DSC_PAGER "more" 
# endif
/*............................................................................*/
# ifndef DSC_EDITOR
/* [ editor command: "vi", "vim", "ed", "joe", "emacs", e.g. ] */
   # define DSC_EDITOR "vim"
# endif
/*............................................................................*/
# ifndef DSC_PRINTER
/* [ line printer command: "lp", e.g. ] */
   # define DSC_PRINTER "lp"
# endif
/*............................................................................*/
# ifndef IPT_MAXLBL
   # define IPT_MAXLBL 100
# endif
/*............................................................................*/
# ifndef IPT_SCROLL
/* [ on calling text console scroll <IPT_SCROLL minus 1 > number of lines ] */
   # define IPT_SCROLL 2
# endif
/*............................................................................*/
# ifndef IPT_OPRDEF
/* [ display default operations: 1 - else: 0 ] */
   # define IPT_OPRDEF 0
# endif
/*............................................................................*/
# ifndef IPT_DIVDEF
/* [ display default divisions: 1 - else: 0 ] */
   # define IPT_DIVDEF 0
# endif
/*............................................................................*/
# ifndef IPT_PARDEF
/* [ display default parameters: 1 - else: 0 ] */
   # define IPT_PARDEF 0
# endif
/*............................................................................*/
# ifndef IPT_OPRLOG
/* operation modes logfile name: */
   # define IPT_OPRLOG "dsc.opr"
# endif
/*............................................................................*/
# ifndef IPT_DIVLOG
/* mesh divisions logfile name: */
   # define IPT_DIVLOG "dsc.div"
# endif
/*............................................................................*/
# ifndef IPT_PARLOG
/* geometric and physical parameters logfile name: */
   # define IPT_PARLOG "dsc.par"
# endif
/*............................................................................*/
# define IPT_SCRFLS 2
/* 1: use secure temporary file names [ function mkstemp(*) ] */
/* 2: use fixed file names [ function tmpnam(*) ] */
# if defined ( _BSD )
   # undef IPT_SCRFLS
   # define IPT_SCRFLS 1
# endif
/*............................................................................*/
# if defined ( _GNU_Linux )
   # undef IPT_SCRFLS
   # define IPT_SCRFLS 1 
# endif
/*............................................................................*/
# if defined ( _Linux )
   # undef IPT_SCRFLS
   # define IPT_SCRFLS 1
# endif
/*............................................................................*/
# if USE_NCURSES == 1
/* 'my_terminal' configuration: */

   # include <termcap.h>     /* terminal type header */
   static char *term;        /* terminal type string */ 

   # define CLSCREEN {\
     printf( "%s", tgetstr( "cl", null )); \
   }

   # define PRBLDCLR(a) {\
     printf( "%s%s", tgetstr( "md", null ), (a)); /* bold clear output */ \
   }

   # define PRINVERS(a) {\
     printf( "%s%s", tgetstr( "mr", null ), (a)); /* inverse */ \
   }

   # define PRNORMAL(a) {\
     printf( "%s%s", tgetstr( "me", null ), (a)); /* back to normal output */ \
   }
# else
   # define CLSCREEN { ;\
   }

   # define PRBLDCLR(a) {\
     printf( "%s", (a));\
   }

   # define PRINVERS(a) {\
     printf( "%s", (a));\
   }

   # define PRNORMAL(a) {\
     printf( "%s", (a));\
   }
# endif
/*............................................................................*/
# if IPT_SCRFLS == 1
   # include <unistd.h>
# endif
/*----------------------------------------------------------------------------*/
# include "readopr.h"
# include "readdiv.h"
# include "readpar.h"
# include "storeopr.h"
# include "storediv.h"
# include "storepar.h"
/*============================================================================*/

short input ( char *option )
{
/* allusions: */
/*
   extern FORMSTATE *spt;
   extern struct transfer trf;
   extern struct blcstruc blc;
*/
/* declarations: */

   static TXCNSL 
     *csp = NULL;

   static short
      ii = null,
      jj = null,
      ind = null,
      item = null,
      parameters = MOD_PARMTRS;

   static char 
      sshll = ONE,
      ptr[STS_SIZE] = {null},
      fleptr[STS_SIZE] = {null},
      tmpfle[STS_SIZE] = {null},
      timeptr[STS_SIZE] = {null},
     *command = NULL; 
/*
      command[2*STS_SIZE] = {null};
*/
   const char
      items0 = SIX,
      items1 = SIX,
      items2 = SIX,
     *oprlog = IPT_OPRLOG,
     *divlog = IPT_DIVLOG,
     *parlog = IPT_PARLOG,
     *timefrm = " created: %.24s  ";

   time_t
      nseconds = null,
     *timer = null;

# ifndef DSC_PAGER
   static WVGDPAR *wve = NULL;
# endif

/* prototypes: */

   time_t
      time( time_t *timer );

   char
      *ctime( const time_t *timer );
/*............................................................................*/
# if IPT_SCRFLS == 0
   char *tmpnam( char *s );
# elif IPT_SCRFLS == 1
   int mkstemp( char *s );
# endif
/*............................................................................*/

   double 
      cos( double x ), 
      sin( double x ),
     sqrt( double x );

   void
      deflt_operts( void ),
      deflt_divsns( void ),
      deflt_params( void ),
      set_z_bases( void );

   short 
      rvise_operts( void ),
      rvise_divsns( void ),
      rvise_params( void ),
      rread_operts( char *filename, char mode ),
      rread_divsns( char *filename, char mode ),
      rread_params( char *filename, char mode );

   short 
      store_operts( char *filename, char mode ),
      store_divsns( char *filename, char mode ),
      store_params( char *filename, char mode );

   TXCNSL 
      *txcnsl( TXCNSL *csp );

# ifndef DSC_PAGER

   WVGDPAR
     *wvepar( char *type, double a,  double b,
              double eps, double my, double f );
# endif
/*----------------------------------------------------------------------------*/
# if USE_NCURSES == 1
/* get the terminal info: */

   term = ( char *) getenv( "TERM" ); /* get the terminal */

   ind = tgetent( null, term );

   if( ONE != ind )
   {
      fprintf( stderr, "Error on getting the termcap info\n" );
      exit( EXIT_FAILURE );
   };
# endif
/*............................................................................*/
/* set buffer length = null: */

   ind = setvbuf( stdin, null, _IONBF, null );
   ind = setvbuf( stdout, null, _IONBF, null ); 
/*...........................................................................*/
/* memory allocations; initializations [explicit]: */
/*............................................................................*/
# if IPT_SCRFLS == 1
   strcpy( tmpfle, "/tmp/temp.XXXXXXX" );
# endif
/*............................................................................*/

   if ( sshll == ONE )
   {
      if ( null == system( NULL ))
      {
         sshll = null;

         printf("\n Warning: System calls are inactive !" );
         printf("\n [ There is no command processor ( sshll ) ]" );
         printf( "\n Please acknowledge [ Enter any character ]:" );
         scanf( "%s", ptr );
      };
   };
/*...........................................................................*/
/* option "operations": */

   if ( *option == 'o' )
   {
      if ( trf.c[null] == null )
      {
/* initialize and enter default operations */

         for ( ii=null; ii<=MOD_OPERATS; ii++ )
         {
            for ( jj=null; jj<STS_SIZE; jj++ )   
            {
               trf.ctx[ii][jj] = null;
            };
         };
/*............................................................................*/
         deflt_operts( );        /* enter default operations                  */
/*.............................*/
      };
/*............................................................................*/
      ind = rvise_operts( );      /* revise/reconfigure default divisions     */
/*..............................*/
      if ( ind == null )
      {
         if ( MOD_OPERATS < trf.c[null] )
         {
            printf( "\n\n Message from function %s :", DO_INPUT );
            printf( "\n\n Too many operations defined !!!" );
            printf( "\n [ Number %d exceeds maximum %d = macro MOD_OPERATS",
               trf.c[null], MOD_OPERATS );
            printf( "\n   in file 'model.h'.]" );
            printf( "\n   - Change macro in compliance with memory "
               "resources.\n" );
            exit( EXIT_FAILURE );
         };

         if ( trf.ref == ONE )
            strcpy( ( spt->tpt->text ), "reference_line" );
         else if ( null == strncmp(( spt->tpt->text ), "reference", FIVE ))
            strcpy(( spt->tpt->text), "unspecified_device" );
      };
/*............................................................................*/
# if defined ( IPT_OPRTMP )
/* store the actually charged operations: */
      strcpy( tmpfle, IPT_OPRTMP );
# else
/*............................................................................*/
/* assign a temporary file name: */
# if IPT_SCRFLS == 0
      tmpnam( tmpfle ); /* [ less secure ] */
# elif IPT_SCRFLS == 1
      mkstemp( tmpfle ); /* [ secure ] */
# else
      strcpy( tmpfle, "oprXXXXXXX" );
# endif                                                        
# endif                                                        
/*............................................................................*/

      strcpy( fleptr, tmpfle );

/*............................................................................*/
      store_operts( fleptr, 't' );      /* store operations on temporary file */
/*....................................*/
# if IPT_OPRDEF == 1
      command = calloc( ONE, (2*STS_SIZE ));
      strcpy( command, DSC_PAGER );
      strcat( command, " " );
      strcat( command, fleptr );
/*............................................................................*/
      system( command );              /* display presently stored operations  */
      rread_operts( fleptr, 't' );   /*                                       */
      rvise_operts( );              /*                                        */
      store_operts( fleptr, 't' ); /*    restore on temporary file            */
/*...............................*/
# endif
/*............................................................................*/
/* initialize text console and menu */

      csp = txcnsl( null );
      item = items0;

/* if IPT_SCROLL !=0: clear screen; scroll IPT_SCROLL-ONE number of lines */
      ( csp->clscr ) = IPT_SCROLL; 

     opr_menu:

      nseconds = time( timer );
      strcpy( timeptr, ctime( &nseconds ));

      ( csp->items ) = items0;
      ( csp->dfopt ) = item;
      ( csp->dflnf ) = item; /* set line feed before option menu line */

      strcpy(( csp->title ), "DSC model generation program FORMER: " );
      strncat(( csp->title ), timeptr, 24 );
      strcpy(( csp->envmt ), "INPUT" );
      strcpy(( csp->cmmnt ), "Computation modes" );
      strcpy(( csp->tasks ), "OPERATIONS" );
      strcpy(( csp->mline[1] ), "* Display the presently charged operations" );
      strcpy(( csp->mline[2] ), "* Enter new operations from file" );
      strcpy(( csp->mline[3] ), "* Edit [ and evtly. modify ] the operations" );
      strcpy(( csp->mline[4] ), "* Reload the default operations" );
      strcpy(( csp->mline[5] ), "* Print the operations" );
      strcpy(( csp->mline[6] ), "* Continue" );
      strcpy(( csp->escpe ), "End of program / escape" );
/*............................................................................*/
      csp = txcnsl( csp );   /* call the menu displaying text console         */
/*.........................*/

      item = ( csp->option );

/*............................................................................*/
      
      switch ( item )
      {
        default: /* continue program [ e.g. item = items0 ] */
         break;

        case 1:  /* display the actual configuration [ on screen ] */
/*............................................................................*/
# if defined ( DSC_PAGER )
         command = calloc( ONE, (2*STS_SIZE ));
         strcpy( command, DSC_PAGER );
         strcat( command, " " );
         strcat( command, fleptr );
/*............................................................................*/
         system( command );      /* edit [evtentually modify] operations file */
/*.............................*/
# else
         printf( "\n\n DSC model '%s' ", ( spt->tpt->name ) );
         printf( "\n comment: %s ", ( spt->tpt->text ) );

         printf( "\n\n %03d  %s\n", trf.c[null], trf.ctx[null] );

         for ( ii=ONE ; ii<=trf.c[null] ; ii++ )
            printf( "\n n%03d;%s %3d", ii, trf.ctx[ii], trf.c[ii] );
# endif
/*............................................................................*/
         item = items0;
         ( csp->clscr ) = 0; /* if N!=0: clear screen, scroll N-1 lines */
         goto opr_menu;
         break;

        case 2: /* enter another configuration from file */

         strcpy(( csp->rqfrm ), "points" );
         strcpy(( csp->rqstr ), "Please enter filename "\
            "[ Continue/Escape: null ]" );
      /* strcpy(( csp->dfstr ), "0" ); */
/*............................................................................*/
         csp = txcnsl( csp );       /*                                        */
/*................................*/
         strcpy( ptr, ( csp->instr ));

         if( *ptr == '0' )
         { 
            PRBLDCLR( "" );
            printf( "\r %*s", 78, "INPUT" );
            PRNORMAL( "" );
            printf( "\n ==================================="
               "===========================================" );
            return -ONE; 
         };
/*............................................................................*/
         rread_operts( ptr, 'o' );      /* enter the new operations           */
         rvise_operts( );              /*  revise/reconfigure ...             */
         store_operts( fleptr, 't' ); /*   restore on the temporary file      */
/*..................................*/
         item = items0;
         ( csp->clscr ) = IPT_SCROLL; 
         goto opr_menu;
         break;

        case 3:  /* edit and/or modify the actual configuration */
/*............................................................................*/
# if defined ( DSC_EDITOR )
         command = calloc( ONE, (2*STS_SIZE ));
         strcpy( command, DSC_EDITOR );
         strcat( command, " " );
         strcat( command, fleptr );
/*............................................................................*/
         system( command );              /* edit [and modify] operations file */
         rread_operts( fleptr, 't' );   /*  enter any new configuration       */
         rvise_operts( );              /*   revise/reconfigure ...            */
         store_operts( fleptr, 't' ); /*    restore on the temporary file     */
/*..................................*/
# else
         printf( "\n No editor DSC_EDITOR defined in CONFIG.H !" );
# endif
/*............................................................................*/
         item = items0;
         ( csp->clscr ) = IPT_SCROLL; 
         goto opr_menu;
         break;

        case 4:  /* reload the default configuration */

/*............................................................................*/
         deflt_operts( );               /* enter default operations           */
         rvise_operts( );              /*  revise/reconfigure ...             */
         store_operts( fleptr, 't' ); /*   restore on the temporary file      */
/*..................................*/
         item = items0;
         ( csp->clscr ) = IPT_SCROLL; 
         goto opr_menu;
         break;

        case 5:

/*............................................................................*/
# if defined ( DSC_PRINTER )
         command = calloc( ONE, (2*STS_SIZE ));
         strcpy( command, DSC_PRINTER );
         strcat( command, " " );
         strcat( command, fleptr );
/*............................................................................*/
         system( command );       /* print actually charged operations        */
/*..............................*/
# else
         printf( "\n No printer DSC_PRINTER defined in CONFIG.H !" );
# endif
/*............................................................................*/
         item = items0;
         ( csp->clscr ) = IPT_SCROLL; 
         goto opr_menu;
         break;

        case 0:  /* end of program / escape */

         remove( tmpfle );
         PRBLDCLR( "" );
         printf( "\r %*s", 78, "INPUT" );
         PRNORMAL( "" );
         return -ONE;
         break;
      }; /* end switch(*) */
/*............................................................................*/
      strcpy( fleptr, oprlog  );
      strcat( fleptr, spt->flbl );

      printf( "\n opened: Operations logfile %s", fleptr );
/*............................................................................*/
      store_operts( fleptr, 'o' );  /* store final operation istructions      */
/*................................*/

      nseconds = time( timer );
      strcpy( timeptr, ctime( &nseconds ));
      printf( "\r Operations logfile '%s'", fleptr );
      printf( timefrm, timeptr );
      printf( CLEAR_LINE );

      return null;

   }; /* end if *option == 'o'[perations] */

/*............................................................................*/
/* option "divisions": */

   if ( *option == 'd' )
   {
      if ( trf.n[null] == null )
      {
/* initialize; enter default divisions */
         for ( ii=null; ii<=MOD_DOMAINS; ii++ )
         {
            for ( jj=null; jj<STS_SIZE; jj++ )   
            {
               trf.mtx[ii][jj] = null;
            };
         };

         for ( ii=null; ii<=MOD_DIVISNS; ii++ )
         {
            for ( jj=null; jj<STS_SIZE; jj++ )   
            {
               trf.ntx[ii][jj] = null;
            };
         };
/*............................................................................*/
         deflt_divsns( );        /* enter default divisions                   */
/*.............................*/
      };
/*............................................................................*/
      ind = rvise_divsns( );      /* revise/reconfigure default divisions     */
/*..............................*/

      if ( ind == null )
      {
         if ( MOD_DOMAINS < blc.m[null] )
         {
            printf( "\n\n Message from function %s :", DO_INPUT );
            printf( "\n\n Too many z_divisions !!!" );
            printf( "\n [ Number %d exceeds maximum %d = macro MOD_DOMAINS",
               blc.m[null], MOD_DOMAINS );
            printf( "\n   in file 'model.h'.]" );
            printf( "\n   - Change macro in compliance with memory "
               "resources.\n" );
            exit( EXIT_FAILURE );
         };

         if ( MOD_DIVISNS < trf.n[null] )
         {
            printf( "\n\n Message from function %s :", DO_INPUT );
            printf( "\n\n Too many xy_divisions !!!" );
            printf( "\n [ Number %d exceeds maximum %d = macro MOD_DIVISNS",
               trf.n[null], MOD_DIVISNS );
            printf( "\n   in file 'model.h'.]" );
            printf( "\n   - Change macro in compliance with memory "
               "resources.\n" );
            exit( EXIT_FAILURE );
         };

         if ( trf.ref == ONE )
            strcpy( ( spt->tpt->text ), "reference_line" );
         else if ( null == strncmp(( spt->tpt->text ), "reference", FIVE ))
            strcpy(( spt->tpt->text), "unspecified_device" );
      };
/*............................................................................*/
# if defined ( IPT_DIVTMP )
/* store the actually charged divisions: */
      strcpy( tmpfle, IPT_DIVTMP );
# else
/*............................................................................*/
/* assign a temporary file name: */
# if IPT_SCRFLS == 0
      tmpnam( tmpfle ); /* [ less secure ] */
# elif IPT_SCRFLS == 1
      mkstemp( tmpfle ); /* [ secure ] */
# else
      strcpy( tmpfle, "divXXXXXXX" );
# endif                                                        
# endif                                                        
/*............................................................................*/

      strcpy( fleptr, tmpfle );

/*............................................................................*/
      store_divsns( fleptr, 't' );      /* store divisions on temporary file  */
/*....................................*/
# if IPT_DIVDEF == 1
      command = calloc( ONE, (2*STS_SIZE ));
      strcpy( command, DSC_PAGER );
      strcat( command, " " );
      strcat( command, fleptr );
/*............................................................................*/
      system( command );              /* edit [eventually modify] divisions   */
      rread_divsns( fleptr, 't' );   /*                                       */
      rvise_divsns( );              /*                                        */
      store_divsns( fleptr, 't' ); /*    restore on temporary file            */
/*...............................*/
# endif
/*............................................................................*/

      csp = txcnsl( null ); /* initialize menu */
      item = items1;

/* if IPT_SCROLL !=0: clear screen; scroll IPT_SCROLL-ONE number of lines */
      ( csp->clscr ) = IPT_SCROLL; 

     div_menu:

      nseconds = time( timer );
      strcpy( timeptr, ctime( &nseconds ));

      ( csp->items ) = items1;
      ( csp->dfopt ) = item;
      ( csp->dflnf ) = item; /* set line feed before option menu line */

      strcpy(( csp->title ), "DSC model generation program FORMER: " );
      strncat(( csp->title ), timeptr, 24 );
      strcpy(( csp->envmt ), "INPUT" );
      strcpy(( csp->cmmnt ), "Model topology and mesh" );
      strcpy(( csp->tasks ), "DIVISIONS" );
      strcpy(( csp->mline[1] ), "* Display the presently charged divisions" );
      strcpy(( csp->mline[2] ), "* Enter new divisions from file" );
      strcpy(( csp->mline[3] ), "* Edit [ and evtly. modify ] the divisions" );
      strcpy(( csp->mline[4] ), "* Reload the default divisions" );
      strcpy(( csp->mline[5] ), "* Print the divisions" );
      strcpy(( csp->mline[6] ), "* Continue" );
      strcpy(( csp->escpe ), "End of program / escape" );
/*............................................................................*/
      csp = txcnsl( csp );   /* call the menu displaying text console         */
/*.........................*/

      item = ( csp->option );

/*............................................................................*/
      
      switch ( item )
      {
        default: /* continue program [ e.g. item = items1 ] */
         break;

        case 1:  /* display the actual configuration [ on screen ] */
/*............................................................................*/
# if defined ( DSC_PAGER )
         command = calloc( ONE, (2*STS_SIZE ));
         strcpy( command, DSC_PAGER );
         strcat( command, " " );
         strcat( command, fleptr );
/*............................................................................*/
         system( command );            /* edit [eventually modify] divisions  */
/*...................................*/
# else
         printf( "\n\n DSC model '%s' ", ( spt->tpt->name ) );
         printf( "\n comment: %s ", ( spt->tpt->text ) );

         printf( "\n\n %03d  %s \n", blc.m[null], trf.mtx[null] );

         for ( ii=ONE ; ii<=blc.m[null] ; ii++ )
            printf( "\n m%03d;%s %3d", ii, trf.mtx[ii], blc.m[ii] );

         printf( "\n\n %03d  %s\n", trf.n[null], trf.ntx[null] );

         for ( ii=ONE ; ii<=trf.n[null] ; ii++ )
            printf( "\n n%03d;%s %3d", ii, trf.ntx[ii], trf.n[ii] );
# endif
/*............................................................................*/
         item = items1;
         ( csp->clscr ) = 0; /* if N!=0: clear screen, scroll N-1 lines */
         goto div_menu;
         break;

        case 2: /* enter another configuration from file */

         strcpy(( csp->rqfrm ), "points" );
         strcpy(( csp->rqstr ), "Please enter filename "\
            "[ Continue/Escape: null ]" );
      /* strcpy(( csp->dfstr ), "0" ); */
/*............................................................................*/
         csp = txcnsl( csp );       /*                                        */
/*................................*/
         strcpy( ptr, ( csp->instr ));

         if( *ptr == '0' )
         { 
            PRBLDCLR( "" );
            printf( "\r %*s", 78, "INPUT" );
            PRNORMAL( "" );
            printf( "\n ==================================="
               "===========================================" );
            return -ONE; 
         };
/*............................................................................*/
         rread_divsns( ptr, 'd' );      /* enter the new divisions            */
         rvise_divsns( );              /*  revise/reconfigure ...             */
         store_divsns( fleptr, 't' ); /*   restore on temporary file          */
/*..................................*/
         item = items1;
         ( csp->clscr ) = IPT_SCROLL; 
         goto div_menu;
         break;

        case 3:  /* edit and/or modify the actual configuration */

/*............................................................................*/
# if defined ( DSC_EDITOR )
         command = calloc( ONE, (2*STS_SIZE ));
         strcpy( command, DSC_EDITOR );
         strcat( command, " " );
         strcat( command, fleptr );
/*............................................................................*/
         system( command );             /* edit [eventually modify] divisions */
         rread_divsns( fleptr, 't' );  /*  reenter ...                        */
         rvise_divsns( );             /*   revise/reconfigure ...             */
         store_divsns( fleptr, 't' );/*    restore divs on temporary file     */
/*.................................*/
# else
         printf( "\n No editor DSC_EDITOR defined in CONFIG.H !" );
# endif
/*............................................................................*/
         item = items1;
         ( csp->clscr ) = IPT_SCROLL; 
         goto div_menu;
         break;

        case 4:  /* reload the default configuration */

/*............................................................................*/
         deflt_divsns( );               /* enter default divisions            */
         rvise_divsns( );              /*  revise/reconfigure ...             */
         store_divsns( fleptr, 't' ); /*   restore on temporary file          */
/*..................................*/
         item = items1;
         ( csp->clscr ) = IPT_SCROLL; 
         goto div_menu;
         break;

        case 5:

/*............................................................................*/
# if defined ( DSC_PRINTER )
         command = calloc( ONE, (2*STS_SIZE ));
         strcpy( command, DSC_PRINTER );
         strcat( command, " " );
         strcat( command, fleptr );
/*............................................................................*/
         system( command );       /* print actually charged divisions         */
/*..............................*/
# else
         printf( "\n No printer DSC_PRINTER defined in CONFIG.H !" );
# endif
/*............................................................................*/
         item = items1;
         ( csp->clscr ) = IPT_SCROLL; 
         goto div_menu;
         break;

        case 0:  /* end of program / escape */

         remove( tmpfle );
         PRBLDCLR( "" );
         printf( "\r %*s", 78, "INPUT" );
         PRNORMAL( "" );
         return -ONE;
         break;
      }; /* end switch(*) */
/*............................................................................*/
      strcpy( fleptr, divlog  );
      strcat( fleptr, spt->flbl );

      printf( "\n opened: Divisions logfile %s", fleptr );
/*............................................................................*/
      store_divsns( fleptr, 'd' );  /* store final divisions on log file      */
/*................................*/

      nseconds = time( timer );
      strcpy( timeptr, ctime( &nseconds ));
      printf( "\r Divisions logfile '%s'", fleptr );
      printf( timefrm, timeptr );
      printf( CLEAR_LINE );

      return null;

   }; /* end if *option == 'd'[ivisions] */

/*............................................................................*/
/* option "parameters": */

   if ( *option == 'p' )
   {
      if ( ( short ) trf.s[null] == null ) 
      {
/* initialize; enter default parameters: */

         trf.wgtype = ( char *) calloc ( STS_SIZE, ONE );

         for ( ii=null; ii<=MOD_PARMTRS; ii++ )
         {
            for ( jj=null; jj<STS_SIZE; jj++ )   
            {
               trf.stx[ii][jj] = null;
            };
         };

         strncpy(( spt->ppt->name ), ( spt->tpt->name ), SHS_SIZE );
         strncpy(( spt->ppt->text ), ( spt->tpt->text ), STS_SIZE );
/*............................................................................*/
         deflt_params( );          /*                                         */
/*...............................*/

         parameters = ( short ) trf.s[null];
      };

/*............................................................................*/
      ind = rvise_params( );      /*                                          */
/*..............................*/
      
      if ( ind == null )
      {
         if ( MOD_PARMTRS < parameters )
         {
            printf( "\n\n Message from function %s :", DO_INPUT );
            printf( "\n\n Too many parameters !!!" );
            printf( "\n [ Number %d exceeds maximum %d = macro MOD_PARMTRS",
               parameters, MOD_PARMTRS );
            printf( "\n   in file 'model.h'.]" );
            printf( "\n   - Change macro in compliance with memory "
               "resources.\n" );
            exit( EXIT_FAILURE );
         };

         if ( trf.ref == ONE )
            strcpy(( spt->ppt->text ), "reference_line" );
         else if ( null == strncmp(( spt->ppt->text ), "reference", FIVE ))
            strcpy(( spt->ppt->text), "unspecified_device" );
      };
/*............................................................................*/
# if defined ( IPT_PARTMP )
/* store the actually charged parameters: */
      strcpy( tmpfle, IPT_PARTMP );
# else
/*............................................................................*/
/* assign a temporary file name: */
# if IPT_SCRFLS == 0
      tmpnam( tmpfle ); /* [ less secure ] */
# elif IPT_SCRFLS == 1
      mkstemp( tmpfle ); /* [ secure ] */
# else
      strcpy( tmpfle, "parXXXXXXX" );
# endif                                                        
# endif                                                        
/*............................................................................*/

      strcpy( fleptr, tmpfle );

/*............................................................................*/
      store_params( fleptr, 't' );      /* store parameters on temporary file */
/*....................................*/
# if IPT_PARDEF == 1
      command = calloc( ONE, (2*STS_SIZE ));
      strcpy( command, DSC_PAGER );
      strcat( command, " " );
      strcat( command, fleptr );
/*............................................................................*/
      system( command );              /* edit [eventually modify] parameters  */
      rread_params( fleptr, 't' );   /*  reenter ...                          */
      rvise_params( );              /*   revise/reconfigure ...               */
      store_params( fleptr, 't' ); /*    restore on temporary file            */
/*...............................*/
# endif
/*............................................................................*/

      csp = txcnsl( null ); /* initialize menu */
      item = items2;

/* if IPT_SCROLL !=0: clear screen; scroll IPT_SCROLL-ONE number of lines */
      ( csp->clscr ) = IPT_SCROLL; 

     par_menu:

      nseconds = time( timer );
      strcpy( timeptr, ctime( &nseconds ));

      ( csp->items ) = items2;
      ( csp->dfopt ) = item;
      ( csp->dflnf ) = item; /* set line feed before option menu line */

      strcpy(( csp->title ), "DSC model generation program FORMER: " );
      strncat(( csp->title ), timeptr, 24 );
      strcpy(( csp->envmt ), "INPUT" );
      strcpy(( csp->cmmnt ), "Model geometry and media" );
      strcpy(( csp->tasks ), "PARAMETERS" );
      strcpy(( csp->mline[1] ), "* Display the presently charged parameters" );
      strcpy(( csp->mline[2] ), "* Enter new parameters from file" );
      strcpy(( csp->mline[3] ), "* Edit [ and evtly. modify ] the parameters" );
      strcpy(( csp->mline[4] ), "* Reload the default parameters" );
      strcpy(( csp->mline[5] ), "* Print the parameters" );
      strcpy(( csp->mline[6] ), "* Continue" );
      strcpy(( csp->escpe ), "End of program / escape" );
/*............................................................................*/
      csp = txcnsl( csp );   /* call the menu displaying text console         */
/*.........................*/

      item = ( csp->option );

/*............................................................................*/
      
      switch ( item )
      {
        default: /* continue program [ e.g item = items2 ] */
         break;

        case 1:  /* display the actual configuration [ on screen ] */
/*............................................................................*/
# if defined ( DSC_PAGER )
         command = calloc( ONE, (2*STS_SIZE ));
         strcpy( command, DSC_PAGER );
         strcat( command, " " );
         strcat( command, fleptr );
/*............................................................................*/
         system( command );          /* display parameters with pager         */
/*.................................*/
# else
         printf( "\n\n DSC model '%s' ", ( spt->ppt->name ));
         printf( "\n comment: %s ", ( spt->ppt->text ));
         printf( "\n %s\n", trf.stx[null] );

         for ( ii=ONE ; ii<=parameters ; ii++ )
            printf( "\n s%03d;%s % .12e ", ii, trf.stx[ii], trf.s[ii] );
      
         printf( "\n\n frequency________________"
            "_________________________[Hz]:  % .12e ", trf.fr ); 

/*  waveguide parameters: */
   
         if ( null == strncmp( trf.wgtype, "tem", THREE ))
            wve = wvepar( trf.wgtype, ZERO, ZERO, 1., 1., trf.fr );
         else if ( null == strncmp( trf.wgtype, "micro_strip", THREE ))
            wve = wvepar( trf.wgtype, ZERO, ZERO, 1., 1., trf.fr );
         else if ( null == strncmp( trf.wgtype, "strip_line", THREE ))
            wve = wvepar( trf.wgtype, ZERO, ZERO, 1., 1., trf.fr );
         else if ( null == strncmp( trf.wgtype, "coax", THREE ))
            wve = wvepar( trf.wgtype, 2.*trf.ra, 2.*trf.ri, 1., 1., trf.fr );
         else if ( null == strncmp( trf.wgtype, "elliptic", THREE ))
            wve = wvepar( trf.wgtype, 2.*trf.ra, 2.*trf.ri, 1., 1., trf.fr );
         else if ( null == strncmp( trf.wgtype, "circular", THREE ))
            wve = wvepar( trf.wgtype, trf.a, trf.a, 1., 1., trf.fr );
         else
            wve = wvepar( trf.wgtype, trf.a, trf.b, 1., 1., trf.fr );

         printf( "\n\n waveguide_type______________"
            "__________________________:   %s ", trf.wgtype );
         printf( "\n\n [ operation parameters at %.14e  GHz: ", trf.fr*1.e-9 );

         if(( wve->q ) < 1. )
         {
            printf( "\n   free wavelength .................: "
               "%.7e m ", ( wve->w0 ));
            printf( "\n   waveguide wavelength ............: "
               "%.7e m ", ( wve->wg ));
            printf( "\n   TEmn mode char.line impedance ...: "
               "%.7e Ohm ", ( wve->zh ));
            printf( "\n   phase velocity ..................: "
               "%.7e m/s ", ( wve->vp ));
            printf( "\n   group velocity ..................: "
               "%.7e m/s ] ", ( wve->vg ));
            printf( "\n\n critical_frequency________________"
               "________________[Hz]:  % .12e ", wve->fc );
         }
         else 
         {
	    printf( "\n   ---------- frequency below cutoff  !!! "
               "----------- ] " );
         };

         if( *( spt->ppt->domain ) == 'f' )
            printf( "\n\n FREQUENCY DOMAIN " );
         else
            printf( "\n\n TIME DOMAIN      " );
# endif
/*............................................................................*/
         item = items2;
         ( csp->clscr ) = 0; /* if N!=0: clear screen, scroll N-1 lines */
         goto par_menu;
         break;

        case 2: /* enter another configuration from file */

         strcpy(( csp->rqfrm ), "points" );
         strcpy(( csp->rqstr ), "Please enter filename "\
            "[ Continue/Escape: null ]" );
      /* strcpy(( csp->dfstr ), "0" ); */
/*............................................................................*/
         csp = txcnsl( csp );       /*                                        */
/*................................*/
         strcpy( ptr, ( csp->instr ));

         if( *ptr == '0' )
         { 
            PRBLDCLR( "" );
            printf( "\r %*s", 78, "INPUT" );
            PRNORMAL( "" );
            printf( "\n ==================================="
               "===========================================" );
            return -ONE; 
         };
/*............................................................................*/
         rread_params( ptr, 'p' );      /* enter the new parameters           */
         rvise_params( );              /*  revise/reconfigure ...             */
         store_params( fleptr, 't' ); /*   restore on temporary file          */
/*..................................*/
         item = items2;
         ( csp->clscr ) = IPT_SCROLL; 
         goto par_menu;
         break;

        case 3:  /* edit and/or modify the actual configuration */

/*............................................................................*/
# if defined ( DSC_EDITOR )
         command = calloc( ONE, (2*STS_SIZE ));
         strcpy( command, DSC_EDITOR );
         strcat( command, " " );
         strcat( command, fleptr );
/*............................................................................*/
         system( command );             /* edit [eventually modify] parameters*/
         rread_params( fleptr, 't' );  /*  reenter parameters                 */
         rvise_params( );             /*   revise/reconfigure ...             */
         store_params( fleptr, 't' );/*    restore on temporary file          */
/*..................................*/
# else
         printf( "\n No editor DSC_EDITOR defined in CONFIG.H !" );
# endif
/*............................................................................*/
         item = items2;
         ( csp->clscr ) = IPT_SCROLL; 
         goto par_menu;
         break;

        case 4:  /* reload the default configuration */

/*............................................................................*/
         deflt_params( );               /* enter default parameters           */
         rvise_params( );              /*  revise/reconfigure parameters      */
         store_params( fleptr, 't' ); /*   restore on temporary file          */
/*..................................*/
         item = items2;
         goto par_menu;
         break;

        case 5:

/*............................................................................*/
# if defined ( DSC_PRINTER )
         command = calloc( ONE, (2*STS_SIZE ));
         strcpy( command, DSC_PRINTER );
         strcat( command, " " );
         strcat( command, fleptr );
/*............................................................................*/
         system( command );       /* print actually charged parameters        */
/*..............................*/
# else
         printf( "\n No printer DSC_PRINTER defined in CONFIG.H !" );
# endif
/*............................................................................*/
         item = items2;
         ( csp->clscr ) = IPT_SCROLL; 
         goto par_menu;
         break;

        case 0:  /* end of program / escape */

         remove( tmpfle );
         PRBLDCLR( "" );
         printf( "\r %*s", 78, "INPUT" );
         PRNORMAL( "" );
         return -ONE;
         break;
      }; /* end switch(*) */
/*............................................................................*/
      strcpy( fleptr, parlog  );
      strcat( fleptr, spt->flbl );

      printf( "\n opened: Parameters logfile %s", fleptr );
/*............................................................................*/
      store_params( fleptr, 'p' );  /* store final parameters on log file     */
/*................................*/

      printf( CLEAR_LINE );

      nseconds = time( timer );
      strcpy( timeptr, ctime( &nseconds ));
      printf( "\r Parameters logfile '%s'", fleptr );
      printf( timefrm, timeptr );

/*............................................................................*/
/* z-base coordinates: */

/*............................................................................*/
      set_z_bases( );      /*                                                 */
/*.......................*/

      return null;

   };/* end of *option == 'p'arameters */ 

   return null;
} 
/*============================================================================*/
# undef MAX_LABEL_
/***************** end of parameter input function input(*) *******************/
