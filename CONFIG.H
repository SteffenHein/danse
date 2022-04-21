/* [ File: CONFIG.H ] */
/*----------------------------------------------------------------------------*/
/*
*//* Configuration header file for DSC program package DANSE, release 1.0r3
*//* [ Update: April 20, 2022 ]
*/
/*----------------------------------------------------------------------------*/
# define DSC_RELEASE "1.0r3"
# define DSC_UPDATE "2022.03.24"
/*----------------------------------------------------------------------------*/
# define DSC_PRFX "dsc."/* common file prefix for DSC operation files */
/*----------------------------------------------------------------------------*/
/* where the elfe daemon searches the DSC system files [ type .../dsc.xxx ]:  */
# define DSC_WRK "/users/dsc/wrk/"
/* the solver daemon's [solvd] log and pid files: */
# define DSC_LOG "/var/log/dsc.log"
# define DSC_PID "/var/run/dsc.pid"
/*----------------------------------------------------------------------------*/
/* pragmas: */
# ifdef LINES
   # pragma LINES 63
# endif
# ifdef WIDTH
   # pragma WIDTH 80
# endif
/*----------------------------------------------------------------------------*/
# define USE_NCURSES 0
# if defined ( _HPUX )
   # undef USE_NCURSES
   # define USE_NCURSES 0
# elif defined ( _BSD )
   # undef USE_NCURSES
   # define USE_NCURSES 0 /* =1: 'my_terminal' config [ for BSD, Linux, etc.] */
# elif defined ( _GNU_Linux )
   # undef USE_NCURSES
   # define USE_NCURSES 1 /* =1: 'my_terminal' config [ for GNU_Linux ] */
# elif defined ( _Debian )
   # undef USE_NCURSES
   # define USE_NCURSES 0
# elif defined ( _SuSE )
   # undef USE_NCURSES
   # define USE_NCURSES 0
# elif defined ( _Linux )
   # undef USE_NCURSES
   # define USE_NCURSES 0
# elif defined ( _Unix )
   # undef USE_NCURSES
   # define USE_NCURSES 0
# else
   # undef USE_NCURSES
   # define USE_NCURSES 0
# endif
/*----------------------------------------------------------------------------*/
# ifdef _CCBUG
   # define REDHAT_CCBUG 1
# else
   # define REDHAT_CCBUG 0
# endif
/*----------------------------------------------------------------------------*/
/* CPU time stopping mode: */
# define SYS_TIMES 1    /* 2: use function times(*)                           */
                        /* 1: use function difftime(*)    [ ANSI C standard ] */
                        /* 0: use function clock(*)       [ ANSI C standard ] */
/*----------------------------------------------------------------------------*/
/* The following macros specify your preferred pager, editor,                 */
/* and line printer: */
/*............................................................................*/
/* your pager [ command ]: "more", e.g. */
# define DSC_PAGER "more"
/* your editor [ command ]: "vi", "vim, "joe", "emacs", e.g. */
# define DSC_EDITOR "vi"
/* your printer [ command ]: "lp", "lpr", e.g. */
# define DSC_PRINTER "lpr"
/*----------------------------------------------------------------------------*/
/* Disable/enable heat propagation [0/1]: */
# define DSC_HCRMDE 1
/*............................................................................*/
# if DSC_HCRMDE != 0
/*............................................................................*/
/* Temperature unit [ 0: Kelvin, 1: DEG Celsius ] */
    # define TEMPGGE 1 
/*............................................................................*/
/* Disable/enable fluids [0/1]: */
# define DSC_FLDMDE 1
/*............................................................................*/
# if DSC_FLDMDE != 0
/*............................................................................*/
/* Disable/enable compressible flow [0/1] */
/* [ default: 1 ] */
   # define CMPRSSBL 1
/*............................................................................*/
/* Disable/enable Boussinesq approximation [0/1] */
/* [ applies only to CMPRSSBL != 0; default: 1 */
/* - not yet implemented: BSQAPRX != 1 ] */
   # define BSQAPRX 1
/*............................................................................*/
/* Disable/enable local Nusselt numbers [ at no-slip boundary faces, e.g.] */
/* [ default: 0 ] */
   # define NUSSELT 0
/*............................................................................*/
/* Disable/enable special turbulence model [0/N], 0<N */
/* [ laminar:0, Kolmogorov-Prandtl models: 1, 2 ] */
   # define TURBMOD 0
/*............................................................................*/
# endif /* DSC_FLDMDE != 0 */
/*............................................................................*/
# ifndef DSC_INTLCE      /* 1: interlaced internal Maxwell field and          */
   # define DSC_INTLCE 0 /* heat current computation loops                    */
# endif                  /* 0: separate internal loops [allows different      */
                         /* repetition rates for Maxw and hcrr computat]      */
/*............................................................................*/
# ifndef DSC_FCTEMP
   # define DSC_FCTEMP 1 /* 1: store face temperatures                        */
# endif
/*............................................................................*/
# else /* if DSC_HCRMDE == 0 */
   # define DSC_FLDMDE 0
# endif
/*----------------------------------------------------------------------------*/
# ifndef DSC_ADJCLL  /* assign neighboring cell index tpt->mn[][k] that way:  */
   # define DSC_ADJCLL 0 /* 0: to faces [ k is a face index; 0 <= k < 6 ]     */
# endif                   /* 1: to ports [ k is a port index; 0 <= k < 12 ]   */
/*----------------------------------------------------------------------------*/
# ifndef DSC_FCELBL   /* label neighboring face index tpt->fn[][] as follows: */
   # define DSC_FCELBL 0 /* 0: unsigned, started with index null              */
# endif                  /* 1: unsigned, started with index ONE               */
                         /* 2: signed, started with index ONE                 */
                         /* [ the (negative) positive sign indicating respctly*/
                         /*   (un-) alignment with neighboring cell face ]    */
/*----------------------------------------------------------------------------*/
/**************************** end of file CONFIG.H ****************************/
