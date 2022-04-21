/* [ file: linker.c ] */
/*******************************************************************************
*                                                                              *
*   ANSI C function linker(*)                                                  *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   Topological DSC mesh linking subroutine                                    *
*   supporting non-orthogonal condensed mesh                                   *
*   in time & frequency domain                                                 *
*                                                                              *
*   Given any DCS mesh system by specifying the eight ( ordered ) cell vertex  *
*   point indices                                                              *
*                    top.cm[m][k] ;         ( k=0,...,7 )                      *
*                                                                              *
*   for every mesh cell m ( top.mi <= m <= top.mf ), this subroutine           *
*   determines the neighbouring mesh cell and port indices as                  *
*                                                                              *
*                   top.mn[m][p] and top.pn[m][p],                             *
*                                                                              *
*   respectively, for all ports p ( 0 <= p <= 11 ) of cell m.                  *
*   The neigboring port index top.pn[m][p] thereby is a signed index in the    *
*   set { +-1 ,..., +-12 }. Its absolute value indicates the neighboring       *
*   port index+ONE, and the sign gives its orientation relative to p.          *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 30, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# include <stdio.h>
# include <stdlib.h>
# include <stddef.h>
# include <stdarg.h>
# include <string.h>
# include <ctype.h>
/*----------------------------------------------------------------------------*/
# include "../math/maths.h"  /* 'my' computation environment headers */
/*----------------------------------------------------------------------------*/
# ifdef OPTIMIZE
   # pragma OPTIMIZE ON
   # pragma OPT_LEVEL 2 
# endif
/*----------------------------------------------------------------------------*/
# include "../math/consts.h" /* some frequently used constants */
/*----------------------------------------------------------------------------*/
/* Edit and customize this general configuration header: */
# include "../CONFIG.H" 
/*----------------------------------------------------------------------------*/
/* Edit and customize the following header for FORMER.C configuration: */
# include "../former/FORMER.CONF"
/*----------------------------------------------------------------------------*/
/* the following [two] macros should be defined in "../CONFIG.H" */
# ifndef DSC_ADJCLL      /* assign neighboring cell index tpt->mn[][k] so:   */
   # define DSC_ADJCLL 1 /* 0: to faces [ k is a face index; 0 <= k <= 6 ]   */
# endif                   /* 1: to ports [ k is a port index; 1 <= k <= 12 ]  */
/*----------------------------------------------------------------------------*/
# ifndef DSC_FCELBL      /* label neighboring face index tpt->fn[][] so:     */
   # define DSC_FCELBL 0 /* 0: unsigned, started with index null             */
# endif                   /* 1: unsigned, started with index ONE              */
                          /* 2: signed, started with index ONE                */
                          /* [ the negative/positive sign indicating un/align-*/
                          /*   un/alignment ment with neighboring cell face ] */
/*----------------------------------------------------------------------------*/
# if USE_NCURSES
   # include <ncurses.h>
# endif 
/*----------------------------------------------------------------------------*/
/* dimensions: The following constants may have been defined in main program  */
/*            `former.c' [ or elsewhere ]. In that case, any change has to be */
/*             carried out in the original definition !                       */
# ifndef CPNTS
   # define CPNTS 40000 /* maximum number of vertex points */
# endif
# ifndef NODES
   # define NODES 30000 /* maximum number of mesh cells */
# endif
# ifndef NFCES
   # define NFCES 90000 /* maximum number of cell faces */
# endif
# ifndef MXADJ
   # define MXADJ    25 /* maximum number of cells adjacent at a point */
# endif
/*----------------------------------------------------------------------------*/
/* operation marks: */
# define LNK_TOPOLGY  0

# define LNK_DISPLAY  2 /* LNK_DISPLAY 1: Display linker(*) operation.        */
                        /* LNL_DISPLAY 2: Additional moving cursor function   */

# define LNK_VERSION  1 /* select one of equivalent versions 1,2              */
                        /* [ faster version may dependend on compiler ]       */

# define LNK_CLRNBRS  0
                        /* LNK_CLRNBRS = 1:                                   */
                        /* Initial clearing of  neighbouring cell indices.    */
                        /* I.e. reset neighbouring cell identifier of all     */
                        /* ports except boundary ports [ i.e electric or      */
                        /* magnetic walls ] to null.                          */

                        /* Warning:  Take care that  all  b o u n d a r y     */
                        /* ports are correctly set !                          */
                        /* To this end  a l l  ports  ( including  boundary   */
                        /* ports )  should be cleared at the  beginning  of   */
                        /* each run of sygrid(*) .  ( Note that  topology     */
                        /* structure 'top' is automatically initialized only  */
                        /* at the beginning of  f i r s t   elsy job .        */
                        /* Repeated calls of  sygrid(*)  in different  confi- */
                        /* gurations - without initial boundary clearing -    */
                        /* may lead to errors ! )                             */

# define LNK_CPLWLLS 1 
                        /* LNK_CPLWLLS 1: Complete walls.                     */
                        /* If the neighbouring cell index of any face or port */
                        /* is ELECTRIC_ or MAGNETIC_WALL ( -1 or -2 ,         */
                        /* respectively ) all neighbouring cell indicators    */
                        /* of faces/ports on the same face and of faces/ports */
                        /* on adjacent faces [ if such exist ] are set        */
                        /* ELECTRIC_WALL or MAGNETIC_WALL, respectively       */
                        /* [ In conflictive cases priority is given to        */
                        /*   ELECTRIC_WALL. ]                                 */
/*----------------------------------------------------------------------------*/
/* mesh topology defining structure [ cf. header 'ecrrtp.h'] */

typedef struct
{ 
   char
      rtn;

   char 
      name[STS_SIZE],
      text[STS_SIZE];

   long
      mi, mf,                        /* initial and final mesh indices        */
      ci, cf,                        /* init.& final vertex point indices     */
      
# if NFCES != 0
      fi, ff,                        /* initial and final face indices        */
      nf[NFCES+ONE][TWO],            /* nf[I][0/1]: cellis adj to face I */
# endif
# if DSC_ADJCLL == 0
      mn[NODES+ONE][FACES],          /* mn[C][F]: cell adj to cell C, face F  */
# elif DSC_ADJCLL == 1
      mn[NODES+ONE][PORTS],          /* mn[C][P]: cell adj to cell C, port P  */
# endif
      cm[NODES+ONE][CRNRS];          /* vertex point identifier               */
                   
   signed char
      fn[NODES+ONE][FACES],         /* fn[C][F]: face adj to cell C, face F  */
      pn[NODES+ONE][PORTS];         /* neighbouring port identifier          */

} TOPOLOGY;
/*----------------------------------------------------------------------------*/
# if LNK_DISPLAY == 2 
   # include "../solver/dsptyp.h"
# endif
/*----------------------------------------------------------------------------*/
static signed char 
   port1[FACES] = {null},
   port2[FACES] = {null};

void abort( void );

int 
   printf( const char *format,...),
   scanf( const char *format,...);
/*----------------------------------------------------------------------------*/
/* math. function prototypes: */

int abs( int h );
/*----------------------------------------------------------------------------*/
# define ELECTRIC_WALL  -1
# define MAGNETIC_WALL  -2
/*============================================================================*/

TOPOLOGY *
evaluate( long ii, long jj, unsigned char kk, unsigned char ll,
          unsigned char rr, TOPOLOGY *tpt )
{
/* allusions: */
/*
   extern signed char 
      port1[FACES],
      port2[FACES];
*/
/* declarations: */   

# if NFCES != 0
   static long
      ff = null;
# endif

   static signed char 
      sgn =  ONE;

/* function prototypes: - */
/*----------------------------------------------------------------------------*/

# if NFCES != 0
   ff = ( ++( tpt->ff ));

   if ( NFCES < ( tpt->ff ))
   {
      printf( "\n\n Message from function %s :", __func__ );
      printf( "\n\n Too many cell faces defined in DCS mesh %s !!!",
              ( tpt->name ));
      printf( "\n [ Maximum number is %d = macro NFCES"
              " in file %s.", NFCES, "FORMER.CONF" );
      printf( "\n - Change macro only in compliance with memory resources " );
      printf( "\n   and same macro in function scattr(*)"
              " of program SOLVER.C.] ");

      ( tpt->rtn ) = ONE;
      return tpt;
   }
   else
   {
      ( tpt->nf[ff][null] ) = ii;
      ( tpt->nf[ff][ONE] ) = jj;
   };
# endif

# if LNK_CPLWLLS == 1
# if DSC_ADJCLL == 0
   if (( tpt->mn[ii][kk] ) == ELECTRIC_WALL )
   {
      ( tpt->mn[jj][ll] ) = ELECTRIC_WALL;
      return tpt;
   };
   if (( tpt->mn[jj][ll] ) == ELECTRIC_WALL )
   {
      ( tpt->mn[ii][kk] ) = ELECTRIC_WALL;
      return tpt;
   };
   if (( tpt->mn[ii][kk] ) == MAGNETIC_WALL )
   {
      ( tpt->mn[jj][ll] ) = MAGNETIC_WALL;
      return tpt;
   };
   if (( tpt->mn[jj][ll] ) == MAGNETIC_WALL )
   {
      ( tpt->mn[ii][kk] ) = MAGNETIC_WALL;
      return tpt;
   };
# elif DSC_ADJCLL == 1
   if (( tpt->mn[ii][port1[kk]] ) == ELECTRIC_WALL )
   {
      ( tpt->mn[jj][port1[ll]] ) = ELECTRIC_WALL;
      ( tpt->mn[jj][port2[ll]] ) = ELECTRIC_WALL;
      return tpt;
   };
   if (( tpt->mn[jj][port1[ll]] ) == ELECTRIC_WALL ) 
   {
      ( tpt->mn[ii][port1[kk]] ) = ELECTRIC_WALL;
      ( tpt->mn[ii][port2[kk]] ) = ELECTRIC_WALL;
      return tpt; 
   }; 
   if (( tpt->mn[ii][port1[kk]] ) == MAGNETIC_WALL )
   {
      ( tpt->mn[jj][port1[ll]] ) = MAGNETIC_WALL;
      ( tpt->mn[jj][port2[ll]] ) = MAGNETIC_WALL;
      return tpt; 
   };
   if (( tpt->mn[jj][port1[ll]] ) == MAGNETIC_WALL ) 
   {
      ( tpt->mn[ii][port1[kk]] ) = ELECTRIC_WALL;
      ( tpt->mn[ii][port2[kk]] ) = ELECTRIC_WALL;
      return tpt;  
   };
# endif /* end if DSC_ADJCLL == 1 */
# endif /* end if LNK_CPLWLLS == 1 */
/*............................................................................*/
# if DSC_ADJCLL == 0
/* cells adjacent to faces: */

   ( tpt->mn[ii][kk] ) = jj;
   ( tpt->mn[jj][ll] ) = ii;
/*............................................................................*/
# elif DSC_ADJCLL == 1
/* cells adjacent to ports: */

   ( tpt->mn[ii][port1[kk]] ) = jj;
   ( tpt->mn[ii][port2[kk]] ) = jj;
   ( tpt->mn[jj][port1[ll]] ) = ii;
   ( tpt->mn[jj][port2[ll]] ) = ii;
# endif
/*............................................................................*/
/* canonical orientation of adjacent faces: */ 

/* sgn = +ONE: aligned orientation [ as face 0 with face 1, e.g.]             */
/* sgn = -ONE: opposed orientation [ as face 0 with face 2, e.g.]             */

   sgn = TWO*(( kk + ll ) % TWO ) - ONE;

/* faces connected to faces [ marked by an alignment sign ] */

# if DSC_FCELBL == 0
   ( tpt->fn[ii][kk] ) = ll;
   ( tpt->fn[jj][ll] ) = kk;
# elif DSC_FCELBL == 1
   ( tpt->fn[ii][kk] ) = ll + ONE;
   ( tpt->fn[jj][ll] ) = kk + ONE;
# else /* if DSC_FCELBL == 2, e.g. */
   ( tpt->fn[ii][kk] ) = sgn*( ll + ONE );
   ( tpt->fn[jj][ll] ) = sgn*( kk + ONE );
# endif /* DSC_FCELBL == ... */

/* ports connected to ports [ with the proper orientation ] */

   switch ( rr )
   {
     default: /* case null: */
      ( tpt->pn[ii][port1[kk]] ) =   sgn*( port1[ll]+ONE );
      ( tpt->pn[ii][port2[kk]] ) =       ( port2[ll]+ONE );
      ( tpt->pn[jj][port1[ll]] ) =   sgn*( port1[kk]+ONE );
      ( tpt->pn[jj][port2[ll]] ) =       ( port2[kk]+ONE );
      return tpt; 

     case ONE:
      ( tpt->pn[ii][port1[kk]] ) =       ( port2[ll]+ONE );
      ( tpt->pn[ii][port2[kk]] ) = - sgn*( port1[ll]+ONE );
      ( tpt->pn[jj][port1[ll]] ) = - sgn*( port2[kk]+ONE );
      ( tpt->pn[jj][port2[ll]] ) =       ( port1[kk]+ONE );
      return tpt; 

     case TWO:
      ( tpt->pn[ii][port1[kk]] ) = - sgn*( port1[ll]+ONE );
      ( tpt->pn[ii][port2[kk]] ) = -     ( port2[ll]+ONE );
      ( tpt->pn[jj][port1[ll]] ) = - sgn*( port1[kk]+ONE );
      ( tpt->pn[jj][port2[ll]] ) = -     ( port2[kk]+ONE );
      return tpt; 

     case THREE:
      ( tpt->pn[ii][port1[kk]] ) = -     ( port2[ll]+ONE );
      ( tpt->pn[ii][port2[kk]] ) =   sgn*( port1[ll]+ONE );
      ( tpt->pn[jj][port1[ll]] ) =   sgn*( port2[kk]+ONE );
      ( tpt->pn[jj][port2[ll]] ) = -     ( port1[kk]+ONE );

      return tpt; 
   };
}
/*============================================================================*/

TOPOLOGY *
linker( TOPOLOGY *tpt )
{
/* allusions: */
/*
   extern signed char 
      port1[FACES],
      port2[FACES];
*/
/* declarations: */

   static long 
      gg = null,
      hh = null,
      uc = null,
      lc = null,
      mi = null,
      mj = null;

# if LNK_DISPLAY == 2 
   static DSPLAY *dsp;
# endif
     
   static short 
      ii = null,
      jj = null;

   static char *ptr;

   static char **endp = null; 

   static signed char 
       kk = null,
       ll = null,
       mm = null, 
       nn = null,
       pp = null,
       qq = null,
       rr = null;

   static long int mc[CPNTS+ONE][MXADJ] = {{null}};
    
   static signed char pnt[TWO][PORTS][ANGLS] = {{{null}}},
                          face[FACES][FACES] = {{null}},
                             cntr[CPNTS+ONE] = {null};
/* function prototypes: */
    
   TOPOLOGY *
      evaluate( long ii, long jj, unsigned char kk, unsigned char ll,
                                  unsigned char rr, TOPOLOGY *tpt );
# if LNK_DISPLAY == 2 
   DSPLAY 
     *dsplay( DSPLAY *dsp );
# endif
/*----------------------------------------------------------------------------*/
/* memory allocations: */

   ptr  = ( char *) calloc( STS_SIZE , ONE );
/*............................................................................*/
/* cell vertex point labels pnt[nn][kk][rr] of face indexed kk                */
/* nn=null: canonical point labeling scheme; nn=ONE: opposed pnt.lbl.sense    */

   pnt[0][0][0] = 0;       pnt[0][1][0] = 1;       pnt[0][2][0] = 0;
   pnt[0][0][1] = 2;       pnt[0][1][1] = 3;       pnt[0][2][1] = 4;
   pnt[0][0][2] = 6;       pnt[0][1][2] = 7;       pnt[0][2][2] = 5;
   pnt[0][0][3] = 4;       pnt[0][1][3] = 5;       pnt[0][2][3] = 1;

   pnt[0][3][0] = 2;       pnt[0][4][0] = 0;       pnt[0][5][0] = 4;
   pnt[0][3][1] = 6;       pnt[0][4][1] = 1;       pnt[0][5][1] = 5;
   pnt[0][3][2] = 7;       pnt[0][4][2] = 3;       pnt[0][5][2] = 7;
   pnt[0][3][3] = 3;       pnt[0][4][3] = 2;       pnt[0][5][3] = 6; 

   for ( kk=null; kk<FACES; kk++ )
   {
      for ( mm=null; mm<ANGLS; mm++ )
      {
         qq = (( ANGLS + ONE - mm ) % ANGLS );
         pnt[ONE][kk][mm] = pnt[null][kk][qq];
      };
   };
/*............................................................................*/
/* ports of face kk:  port1[kk], port2[kk]                                    */

   port1[0] =  7;          port1[1] =  5;          port1[2] = 11;
   port1[3] =  9;          port1[4] =  3;          port1[5] =  1;
   port2[0] = 10;          port2[1] =  8;          port2[2] =  2;
   port2[3] =  0;          port2[4] =  6;          port2[5] =  4;
/*............................................................................*/
/* default adjacent face[kk][mm]                                              */
/* index kk: face on first given cell,                                        */
/* index mm: order of choice for adjacent face on second cell                 */

   face[0][0] = 1;         face[0][1] = 0;         face[0][2] = 2;
   face[0][3] = 3;         face[0][4] = 4;         face[0][5] = 5;

   face[1][0] = 0;         face[1][1] = 1;         face[1][2] = 3;
   face[1][3] = 2;         face[1][4] = 5;         face[1][5] = 4;

   face[2][0] = 3;         face[2][1] = 2;         face[2][2] = 0;
   face[2][3] = 1;         face[2][4] = 4;         face[2][5] = 5;

   face[3][0] = 2;         face[3][1] = 3;         face[3][2] = 1;
   face[3][3] = 0;         face[3][4] = 5;         face[3][5] = 4;

   face[4][0] = 5;         face[4][1] = 4;         face[4][2] = 0;
   face[4][3] = 1;         face[4][4] = 2;         face[4][5] = 3;

   face[5][0] = 4;         face[5][1] = 5;         face[5][2] = 1;
   face[5][3] = 0;         face[5][4] = 3;         face[5][5] = 2;
/*............................................................................*/

   if ( NODES < ( tpt->mf ))
   {
      printf( "\n\n Message from function %s :", __func__ );
      printf( "\n\n Too many cells used in DCS mesh %s !!!", ( tpt->name ));
      printf( "\n [ Maximum number is %d = macro NODES"
              " in file %s.", NODES, "FORMER.CONF" );
      printf( "\n - Change macro only in compliance with memory resources " );
      printf( "\n   and same macro in function scattr(*)"
              " of program SOLVER.C.] ");

      ( tpt->rtn ) = ONE;
      return tpt;
   };

   if (( tpt->mi ) != ONE )
   {
      printf( "\n\n Warning from function linker(*):" );
      printf( "\n Initial mesh index differs from one !!! "
         "[ index = %ld ]", ( tpt->mi )); 
      printf( "\n (Re-)enter correct initial mesh index. "
         "[ Escape: enter -ONE ]: ");  

      if ( *ptr == '-' )
      {
         ( tpt->rtn ) = ONE;
         return tpt;
      }
      else
      {
         scanf( "%s", ptr );
         ( tpt->mi ) = strtol( ptr, endp, DEC ); 
      };
   };

# if LNK_DISPLAY == 1

   printf( "\n DCS mesh linking function linker(*) started." );
   printf( "\n [ Please wait a moment ] " ); 

# elif LNK_DISPLAY == 2
/*............................................................................*/
/* initialize display [ running cursor function ]: */

   dsp = dsplay( null );
/*............................................................................*/
   strcpy( dsp->messge,
      "DCS mesh linking function linker(*) started. " );
   ( dsp->option ) = 's'; /* 's'tart message */
   dsplay( dsp );

   strcpy( dsp->messge,
      "|[Please wait a moment]" );
   ( dsp->option ) = 'm'; /* 'm'essage under running cursor */
   dsplay( dsp );
# endif

   uc = null;  
   lc = CPNTS;

# if NFCES != 0
   ( tpt->fi ) = ONE;
   ( tpt->ff ) = null;
# endif

# if LNK_CLRNBRS == 1
# if NFCES != 0
   for ( hh = null; hh <= NFCES; h++ )
   {
      tpt->nf[hh][null] = 0;
      tpt->nf[hh][ONE] = 0;
   };
# endif
# endif

   for ( hh=( tpt->mi ); hh<=( tpt->mf ); hh++ )
   {
      for ( kk=null; kk<CRNRS; kk++ )
      {                                       /* determine           */
         if ( uc < ( tpt->cm[hh][kk] ))
            uc = ( tpt->cm[hh][kk] );         /* highest point index */
         if (( tpt->cm[hh][kk] ) < lc )
            lc = ( tpt->cm[hh][kk] );         /* lowest point index  */
         if ( FACES <= kk )
            goto cont1;

# if LNK_CLRNBRS == 1
# if DSC_ADJCLL == 0
         if ( null < ( tpt->mn[hh][kk] ))
            ( tpt->mn[hh][kk] ) = null;
# elif DSC_ADJCLL == 1
         if ( null < ( tpt->mn[hh][port1[kk]] ))
            ( tpt->mn[hh][port1[kk]] ) = null;

         if ( null < ( tpt->mn[hh][port2[kk]] ))
            ( tpt->mn[hh][port2[kk]] ) = null;
# endif /* end if DSC_ADJCLL == 1 */
# endif /* end if LNK_CLRNBRS == 1 */

# if LNK_CPLWLLS == 1
# if DSC_ADJCLL == 1
         if (( tpt->mn[hh][port1[kk]] ) == ELECTRIC_WALL )
         {
            ( tpt->mn[hh][port2[kk]] ) = ELECTRIC_WALL;
            goto cont1;
         };
         if (( tpt->mn[hh][port2[kk]] ) == ELECTRIC_WALL )
         {
            ( tpt->mn[hh][port1[kk]] ) = ELECTRIC_WALL;
            goto cont1;
         };
         if (( tpt->mn[hh][port1[kk]] ) == MAGNETIC_WALL )
         {
            ( tpt->mn[hh][port2[kk]] ) = MAGNETIC_WALL;
            goto cont1;
         };
         if (( tpt->mn[hh][port2[kk]] ) == MAGNETIC_WALL )
         {
            ( tpt->mn[hh][port1[kk]] ) = MAGNETIC_WALL;
            goto cont1;
         };
# endif /* end if DSC_ADJCLL == 1 */
# endif /* end if LNK_CPLWLLS == 1 */

         cont1:;
      };/* next kk [ i.e. next vertex ] */
   };               

   if ( CPNTS < uc )
   {
      printf( " \n\n Message from function %s :", __func__ );
      printf( " \n\n Too many grid points used"
              " in DCS mesh %s !!!", ( tpt->name ));
      printf( "\n [ Maximum point index is %d = macro CPNTS"
              " in file %s.", CPNTS, "FORMER.CONF" );
      printf( "\n - Change macro only in compliance with memory resources.] " );

      ( tpt->rtn ) = ONE;
      return tpt;
   };

# if LNK_CLRNBRS == 1

   printf( "\n\n Message from function linker(*):");
   printf( "\n\n ( non-boundary ) ports cleared. ");

# endif

/* lc [ uc ] : minimum [maximum] cell vertex point index, reset cell counter: */
         
   for ( hh=lc; hh<=uc; hh++ ) 
      cntr[hh] = null;

/* sum up and label all cells adjacent   */
/* to point top.cm[hh]; ii=top.mi...top.mf */
   for( hh=( tpt->mi ); hh<=( tpt->mf ); hh++ )
   {
      for ( kk=null; kk<CRNRS; kk++ )
      {
         gg = ( tpt->cm[hh][kk] );
         mc[gg][cntr[gg]] = hh; 
         cntr[gg]++;

         if ( MXADJ < cntr[gg] )
         {
            printf( "\n\n Message from function %s :", __func__ );
            printf( "\n\n Too many cells adjacent at point %ld  !!!", gg);
            printf( "\n [ Maximum number is %d = macro MXADJ"
                    " in file %s.", MXADJ, "FORMER.CONF" );
            printf( "\n - Change macro only in compliance"
                    " with memory resources.] " ); 

            ( tpt->rtn ) = ONE;
            return tpt;
         };
      };
   };                                /* the total number of those meshes is   */
                                     /* cntr[top.cm[hh]]                      */
/*............................................................................*/    
# if LNK_DISPLAY == 2
   ( dsp->range ) = uc - lc + ONE;
   ( dsp->option ) = 'c';
# endif

   hh = lc; /* hh labels points */
   while( hh <= uc )
   {
      for ( ii=null; ii<cntr[hh]; ii++ ) /*label cells connected with point hh*/
      {
         mi = mc[hh][ii];                         /* mi = 1st. cell index     */
         for ( kk=null; kk<FACES; kk++ )          /* label faces on cell mi   */
         {
/*............................................................................*/
# if LNK_CPLWLLS == 1
# if DSC_ADJCLL == 0
            if (( tpt->mn[mi][kk] ) <= null )
            {
# elif DSC_ADJCLL == 1
            if (( tpt->mn[mi][port1[kk]] ) <= null )
            {
# endif /* end if DSC_ADJCLL == 1 */
# else /* if LNK_CPLWLLS != 1 */
# if DSC_ADJCLL == 0
            if (( tpt->mn[mi][kk] ) == null )
            {
# elif DSC_ADJCLL == 1
            if (( tpt->mn[mi][port1[kk]] ) == null )
            {
# endif /* end if DSC_ADJCLL == 1 */
# endif /* end if LNK_CPLWLLS != 1 */
/*............................................................................*/
               for ( jj=ii+ONE; jj<cntr[hh]; jj++ )
               {
                  mj = mc[hh][jj];               /* 2nd.[adjacent] cell index */
                  if ( mi == mj )
                  {
                     printf( "\n\n Warning message from function"
                        " %s :", __func__ );
                     printf( "\n Coincident vertices in cell %ld, "
                        "point no. %ld !!! ", mi, hh );
                  };
                  mm = null; do
                  {                               /* label priority order mm  */
                     ll = face[kk][mm];           /* of adjacent face         */
                                                  /* [ indexed ll ]           */
/*............................................................................*/
# if LNK_CPLWLLS == 1
# if DSC_ADJCLL == 0
                     if (( tpt->mn[mj][ll] ) <= null )
                     {
# elif DSC_ADJCLL == 1
                     if (( tpt->mn[mj][port1[ll]] ) <= null )
                     {
# endif /* end if DSC_ADJCLL == 1 */
# else /* if LNK_CPLWLLS != 0 */
# if DSC_ADJCLL == 0
                     if (( tpt->mn[mj][ll]] ) == null )
                     {
# elif DSC_ADJCLL == 1
                     if (( tpt->mn[mj][port1[ll]] ) == null )
                     {
# endif /* end if DSC_ADJCLL == 1 */
# endif /* end if LNK_CPLWLLS != 1 */
/*............................................................................*/
/* canonical orientation of adjacent faces: */
/* nn=0: opposed orientation */
/* nn=1: same orientation */

                        nn = (( kk + ll + ONE ) % TWO );
                         
                        rr = null; do
                        {
# if LNK_VERSION == 1 
/* version 'a' */
                           pp = null; do
                           {
                              qq = (( pp + rr ) % ANGLS );
                                              /* compare corner point indices */
                              if (( tpt->cm[mi][pnt[null][kk][pp]] ) != 
                                  ( tpt->cm[mj][pnt[nn][ll][qq]] ))
                              {
                                 if ( null < pp ) 
                                    goto next_mm;
                                 goto next_rr;
                              };
                           } while (( ++pp ) < ANGLS ); /* next angle */

                           tpt = evaluate( mi, mj,
                                 kk, ll, rr, tpt );    /* adjacent face found */
                           goto next_kk;               /* [or walls completd] */
# else
/* version 'b' */
                           pp = null;
                           qq = (( pp + rr ) % ANGLS );
                                              /* compare corner point indices */
                           while(( pp < ANGLS )
                               &&(( tpt->cm[mi][pnt[null][kk][pp]] ) == \
                                  ( tpt->cm[mj][pnt[nn][ll][qq]] )))
                           {
                              pp++; 
                              qq = (( pp + rr ) % ANGLS );
                           };

                           switch( pp )
                           {
                              case null: 
                                 break;     
                                                       
                              case FOUR: 
                                 tpt = evaluate( mi, mj, kk, ll, rr, tpt );

                                 goto next_kk; /* adjacent face found */
                                               /* [or walls completed]*/
                              default:
                                 goto next_mm;
                           };
# endif
                 next_rr:; /* next rr [ face rotation ] */
                        } while (( ++rr ) < ANGLS );
                     }; /* end if top.mn[mj]... */
           next_mm:; /* next mm [ face of cell mj ] */
                  } while (( ++mm ) < FACES );
               }; /* next jj [ second cell connected with point hh ] */
            }; /* end if top.mn[mi]... */
   next_kk:; /* next kk [ face on cell mi ] */
         }; /* next kk */
      }; /* next ii [ first cell connected with point hh ] */

# if LNK_DISPLAY == 2 
      ( dsp->state ) = hh;
      dsplay( dsp );
# endif

      hh++;
   }; /* next hh [ corner point ] */

# if LNK_DISPLAY == 1

   printf( "\r '%s' terminated.        "
           "                           ", __func__ );

# elif LNK_DISPLAY == 2
   ( dsp->messge[null] ) = null;
   ( dsp->option ) = null;
   ( dsp->state ) = null;
   ( dsp->range ) = null;

   dsplay( dsp ); /* clear display */
   printf( "\r %s terminated. ", __func__ );
# endif

   ( tpt->rtn ) = null;
   return tpt;
} 
# ifdef OPTIMIZE
   # pragma OPTIMIZE OFF
# endif

# undef LNK_DISPLAY
# undef LNK_VERSION
# undef LNK_CLRNBRS
# undef LNK_CPLWLLS
# undef ELECTRIC_WALL
# undef MAGNETIC_WALL
/***************** end of DCS mesh linker function linker(*) ******************/
