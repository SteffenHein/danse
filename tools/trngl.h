# define DO_TRNGL "trngl(*)"
/*******************************************************************************
*                                                                              *
*   ANSI C function trngl(*)                                                   *
*   [ DANSE - Discrete Approximation of the Navier Stokes Equations,           *
*     release 1.0 ]                                                            *
*                                                                              *
*   DSC mesh generation function for triangular domains                        *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: March 18, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
/*
called functions: argz(*), subarc(*), arcsct(*)
in header "../math/arcsct.h"
*/
/*----------------------------------------------------------------------------*/
# define TGL_PRINT 0
/*----------------------------------------------------------------------------*/
# define SET_TRIANGLES 0  /* SET_TRIANGLES 1: generate only triangular cells  */
                          /* [ the number of cells will be n = p*p ].         */
                          /* SET_TRIANGLES 0: generate quadrangular cells     */
                          /* within block, and triangular cells along side bc */
                          /* [ the number of cells will be n = p*(p+1)/2 ].   */
/*----------------------------------------------------------------------------*/
# define AB_POINTS PUZZ_POINTS
/*----------------------------------------------------------------------------*/
# define ELECTRIC_WALL -1
# define MAGNETIC_WALL -2
/*----------------------------------------------------------------------------*/
struct triangle 
{
   double 
    ax,ay,az,bx,by,bz,cx,cy,cz, aab,abc,aac;

   long 
    cab[AB_POINTS],
    cbc[2*AB_POINTS],
    cac[AB_POINTS],
    sab[AB_POINTS],
    sbc[2*AB_POINTS],
    sac[AB_POINTS];

   char  
    ab[AB_POINTS],
    bc[2*AB_POINTS], 
    ac[AB_POINTS];
                      
   char top, bot, trv; /* top, bot: boundary wall indicators  */
                       /* trv: trivial node indicator         */
}; 
struct triangle tgl = {ZERO};
/*----------------------------------------------------------------------------*/
# define SET_CORNER(i) \
{ \
   lbl.p++; \
 \
   if ( CPNTS < lbl.p ) \
   { \
      printf( "\n\n Error message from function %s :", __func__ ); \
      printf( "\n Corner point limit exceeded !!!" ); \
      printf( "\n [ Maximum number is %ld = macro CPNTS " \
         "in %s.", (long) CPNTS, "FORMER.CONF" ); \
      printf( "\n - Change macro only in compliance " \
         "with memory resources.]\n" ); \
      exit( EXIT_FAILURE ); \
   }; \
 \
   if ( lbl.pi <= null ) \
      lbl.pi = lbl.p; \
 \
   (i) = lbl.p; \
}
/*----------------------------------------------------------------------------*/
# define SET_NODE( ) \
{ \
   lbl.m++; \
 \
   if ( NODES < lbl.m ) \
   { \
      printf( "\n\n Error message from function %s :", __func__ ); \
      printf( "\n Mesh cell limit exceeded !!!" ); \
      printf( "\n [ Maximum number is %ld = macro NODES " \
         "in %s.", (long) NODES, "FORMER.CONF" ); \
      printf( "\n - Change macro only in compliance with memory " \
         "resources\n   and same macro in function scattr(*) " \
         "of program SOLVER.C.]\n" ); \
      exit( EXIT_FAILURE ); \
   }; \
 \
   if ( lbl.mi <= null ) \
      lbl.mi = lbl.m; \
}
/*============================================================================*/

short trngl( short p, char *option )
{
/* allusions: */
/*
   extern FORMSTATE *spt;
   extern struct labels lbl;
   extern struct triangle tgl; 
   extern struct subarc arc; 
*/
/* declarartions: */

   static double 
      px0 = ZERO,
      py0 = ZERO,
      pz0 = ZERO,
      px1 = ZERO,
      py1 = ZERO,
      pz1 = ZERO,
      qx0 = ZERO,
      qy0 = ZERO,
      qz0 = ZERO,
      qx1 = ZERO,
      qy1 = ZERO,
      qz1 = ZERO,
      rr = ZERO,
      dbc = ZERO;
     /*  ss = ZERO, */

   static long 
      rp[AB_POINTS] = {null};

   static long
      j = null,
      k = null,
      l = null,
      m = null, 
      n = null;

   static short
      h = null,
      i = null,
      q = null,
      ind = null;

   static char
      o = null;

   double 
      argz( double x, double y, short n ),
      subarc( double x, double y, double z, double w,
                                  double alfa, double r );
/*----------------------------------------------------------------------------*/

   if ( AB_POINTS - ONE < p )
   {
      printf( "\n\n Message from function %s :", __func__ );
      printf( "\n Too many divisions in triangular block "
         "generation routine !!!" );
      printf( "\n [ Maximum number on side AB is %d = macro AB_POINTS - 1 "
         "in %s.", AB_POINTS - ONE, __func__ );
      printf( "\n - Change macro in compliance with memory resources.]\n" );
      exit( EXIT_FAILURE );
   };

   lbl.ci = lbl.m + ONE;
/*............................................................................*/

   if ( *option == 'c' )         /* copy boundary point coordinates if        */
   {                             /* boundary points are defined [ indexed>0]  */
      if ( null < tgl.cab[0] )
      {
         k = tgl.cab[0];
         tgl.ax = spt->ppt->cpt->c[k][0];
         tgl.ay = spt->ppt->cpt->c[k][1];
         tgl.az = spt->ppt->cpt->c[k][2];
      };
      if ( null < tgl.cab[p] )
      {
         k = tgl.cab[p];
         tgl.bx = spt->ppt->cpt->c[k][0];
         tgl.by = spt->ppt->cpt->c[k][1];
         tgl.bz = spt->ppt->cpt->c[k][2];
      };
      if ( null < tgl.cac[0] )
      {
         k = tgl.cac[0];
         tgl.ax = spt->ppt->cpt->c[k][0];
         tgl.ay = spt->ppt->cpt->c[k][1];
         tgl.az = spt->ppt->cpt->c[k][2];
      };
      if ( null < tgl.cac[p] )
      {
         k = tgl.cac[p];
         tgl.cx = spt->ppt->cpt->c[k][0];
         tgl.cy = spt->ppt->cpt->c[k][1];
         tgl.cz = spt->ppt->cpt->c[k][2];
      };
      if ( null < tgl.cbc[0] )
      {
         k = tgl.cbc[0];
         tgl.bx = spt->ppt->cpt->c[k][0];
         tgl.by = spt->ppt->cpt->c[k][1];
         tgl.bz = spt->ppt->cpt->c[k][2];
      };
      if ( null < tgl.cbc[2*p] )
      {
         k = tgl.cbc[2*p];
         tgl.cx = spt->ppt->cpt->c[k][0];
         tgl.cy = spt->ppt->cpt->c[k][1];
         tgl.cz = spt->ppt->cpt->c[k][2];
      };
   };
/*............................................................................*/
     
   for ( h=null ; h<=p ; h++ )
   {
      rp[h] = null;
   };

   for ( h=ONE ; h<=p ; h++ )
   {
      if ( *option == 'c' )
      {
         rr = ( double ) ( h-ONE )/p;
/*............................................................................*/
         subarc( tgl.ax, tgl.ay, tgl.cx, tgl.cy, tgl.aac, rr );     /*   */
/*.....................................................................*/
         px0 = arc.rx;
         py0 = arc.ry;
         pz0 = tgl.az + rr*( tgl.cz - tgl.az );

/*............................................................................*/
         subarc( tgl.bx, tgl.by, tgl.cx, tgl.cy, tgl.abc, rr );     /*   */
/*.....................................................................*/
         qx0 = arc.rx;
         qy0 = arc.ry;
         qz0 = tgl.bz + rr*( tgl.cz - tgl.bz );

         dbc = rr*tgl.abc;

         rr = ( double ) h/p;
/*............................................................................*/
         subarc( tgl.ax, tgl.ay, tgl.cx, tgl.cy, tgl.aac, rr );     /*   */
/*.....................................................................*/
         px1 = arc.rx;
         py1 = arc.ry;
         pz1 = tgl.az + rr*( tgl.cz- tgl.az );

/*............................................................................*/
         subarc( tgl.bx, tgl.by, tgl.cx, tgl.cy, tgl.abc, rr );     /*   */
/*.....................................................................*/
         qx1 = arc.rx;
         qy1 = arc.ry;
         qz1 = tgl.bz + rr*( tgl.cz - tgl.bz );

         dbc = rr*tgl.abc - dbc; 

# if TGL_PRINT != null  
         if ( TGL_PRINT == lbl.m + ONE )
         {
            printf( "\n\n Intermediate values from function"
                    " '%s' ", __func__ );
            printf( "\n [ generating cell no. %ld ]: " , lbl.m + ONE );
            printf( "\n\n p0 = ( %.5e , %.5e , %.5e )", px0, py0, pz0);
            printf( "\n q0 = ( %.5e , %.5e , %.5e )", qx0, qy0, qz0);
            printf( "\n p1 = ( %.5e , %.5e , %.5e )", px1, py1, pz1);
            printf( "\n q1 = ( %.5e , %.5e , %.5e )", qx1, qy1, qz1);
            printf( "\n\n [ please acknowledge ] " );
            scanf( "%s", ptr);
         };
# endif
      };/* end if *option == 'c'oordinates */

      q = p-h+ONE;
       
      j = null;
      k = null;
      l = null;    
      m = null;
      n = null;
        
      for ( i=ONE; i<=q; i++ )
      {
         SET_NODE( );

     /* set_j: */

         j = k;
         if ( null < j ) 
            goto set_k;
         if ( h == ONE ) 
            j = tgl.cab[i-ONE];
         if ( null < j ) 
            goto set_k;
         if ( i == ONE ) 
            j = tgl.cac[h-ONE];
         if ( null < j ) 
            goto set_k;

         j = rp[i-ONE]; 

         if ( null < j ) 
            goto set_k;

         SET_CORNER(j);

         if ( *option == 'c' )
         {
            rr = ( double ) ( i-ONE )/q;
/*............................................................................*/
            subarc( px0, py0, qx0, qy0, tgl.aab, rr );       /*          */
/*..............................................................*/
            spt->ppt->cpt->c[j][0] = arc.rx;
            spt->ppt->cpt->c[j][1] = arc.ry;
            spt->ppt->cpt->c[j][2] = pz0 + rr*( qz0 - pz0 );
         };

        set_k:
        
         k = null;
         if ( h == ONE ) 
            k = tgl.cab[i];
         if ( null < k ) 
            goto set_l;  
         if ( i == q ) 
            k = tgl.cbc[2*(h-ONE)];
         if ( null < k ) 
            goto set_l;

         k = rp[i];

         if ( k  > null ) 
            goto set_l; 

         SET_CORNER(k);

         if ( *option == 'c' )
         {
            rr = ( double ) i /q;
/*............................................................................*/
            subarc( px0, py0, qx0, qy0, tgl.aab, rr );       /*         */
/*..............................................................*/
            spt->ppt->cpt->c[k][0] = arc.rx;
            spt->ppt->cpt->c[k][1] = arc.ry;
            spt->ppt->cpt->c[k][2] = pz0 + rr*( qz0 - pz0 );
         };
         
        set_l:

         l = m;
         if ( null < l ) 
            goto set_m;
         if ( i == ONE ) 
            l = tgl.cac[h];
         if ( null < l ) 
            goto set_m;
         if ( i == q ) 
            l = tgl.cbc[2*h];
         if ( null < l ) 
            goto set_m;

         SET_CORNER(l);

         if ( *option == 'c' )
         {
            if ( ONE < q )
            {
               rr = ( double ) ( i-ONE )/( q-ONE );
/*............................................................................*/
               subarc( px1, py1, qx1, qy1, tgl.aab, rr );       /*      */
/*.................................................................*/
               spt->ppt->cpt->c[l][0] = arc.rx;
               spt->ppt->cpt->c[l][1] = arc.ry;
               spt->ppt->cpt->c[l][2] = pz1 + rr*( qz1 - pz1 );
            }
            else
            {
               spt->ppt->cpt->c[l][0] = px1;
               spt->ppt->cpt->c[l][1] = py1;
               spt->ppt->cpt->c[l][2] = pz1;
            };
         };

        set_m: 

         m = null;

         if ( i < q )
	 {
	    if ( i == q-ONE ) 
               m = tgl.cbc[2*h];

	    if ( null < m ) 
               goto set_n;

            SET_CORNER(m);

            if( *option == 'c' )
            { 
               rr = ( double ) i /( q-ONE );
/*............................................................................*/
               subarc( px1, py1, qx1, qy1, tgl.aab, rr );       /*       */
/*.................................................................*/
               spt->ppt->cpt->c[m][0] = arc.rx;
               spt->ppt->cpt->c[m][1] = arc.ry;
               spt->ppt->cpt->c[m][2] = pz1 + rr*( qz1 - pz1 );
            };
         };

        set_n:

         if ( i == q )
         {
            n = tgl.cbc[2*h-ONE];

            if ( null < n ) 
               goto terminal;

            SET_CORNER(n);

            if ( *option == 'c' )
            {
               rr = .5;
/*............................................................................*/
               subarc( qx0, qy0 , qx1, qy1, dbc , rr );       /*         */
/*...............................................................*/
               spt->ppt->cpt->c[n][0] = arc.rx;
               spt->ppt->cpt->c[n][1] = arc.ry;
               spt->ppt->cpt->c[n][2] = qz0 + rr*( qz1 - qz0 );
            };
         };

# if SET_TRIANGLES == 1
         if ( i < q )
         {
            SET_CORNER(n);
 
            if ( *option == 'c' )
            {
               spt->ppt->cpt->c[n][0] = .5*( spt->ppt->cpt->c[k][0] +\
                                             spt->ppt->cpt->c[l][0] );
               spt->ppt->cpt->c[n][1] = .5*( spt->ppt->cpt->c[k][1] +\
                                             spt->ppt->cpt->c[l][1] );
               spt->ppt->cpt->c[n][2] = .5*( spt->ppt->cpt->c[k][2] +\
                                             spt->ppt->cpt->c[l][2] );
            };
         };
# endif
  
        terminal:

         if ( *option == 't' )
         {

# if SET_TRIANGLES == 1

            spt->tpt->cm[lbl.m][0] = j;
            spt->tpt->cm[lbl.m][1] = k;
            spt->tpt->cm[lbl.m][2] = l;
            spt->tpt->cm[lbl.m][3] = n;
# else
            if ( i < q ) 
            {
	       spt->tpt->cm[lbl.m][0] = j;
	       spt->tpt->cm[lbl.m][1] = k;
	       spt->tpt->cm[lbl.m][2] = l;
	       spt->tpt->cm[lbl.m][3] = m;
            }
            else if ( i == q )
	    {
	       spt->tpt->cm[lbl.m][0] = j;
	       spt->tpt->cm[lbl.m][1] = k;
	       spt->tpt->cm[lbl.m][2] = l;
	       spt->tpt->cm[lbl.m][3] = n;
            };
# endif
            ind = null - ( tgl.trv == 'e' ) - 2*( tgl.trv == 'm' );

# if DSC_ADJCLL == 0
            for ( o=null; o<FACES; o++ )  /* clear neighbouring face identif.*/
               spt->tpt->mn[lbl.m][(int)o] = ind;
# elif DSC_ADJCLL == 1
            for ( o=null; o<PORTS; o++ ) /* clear neighbouring port identif.*/
               spt->tpt->mn[lbl.m][(int)o] = ind;
# endif

/* electric and magnetic walls: */

            if ( ind == null ) /* non-trivial cell */
            {
               if (( i == ONE )&&( tgl.ac[h] == 'e' ))      /* face 0 */
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][0] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][7] = ELECTRIC_WALL;
                  spt->tpt->mn[lbl.m][10] = ELECTRIC_WALL;
# endif
               };

               if (( i == ONE )&&( tgl.ac[h] == 'm' ))   
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][0] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][7] = MAGNETIC_WALL;
                  spt->tpt->mn[lbl.m][10] = MAGNETIC_WALL;
# endif
               };
       
               if (( i == q )&&( tgl.bc[2*h-ONE] == 'e' ))  /* face 1 */
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][1] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][5] = ELECTRIC_WALL;
                  spt->tpt->mn[lbl.m][8] = ELECTRIC_WALL;
# endif
               };

               if (( i == q )&&( tgl.bc[2*h-ONE] == 'm' ))
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][1] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][5] = MAGNETIC_WALL;
                  spt->tpt->mn[lbl.m][8] = MAGNETIC_WALL;
# endif
               };

               if (( h == ONE )&&( tgl.ab[i] == 'e' ))      /* face 2 */
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][2] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][2] = ELECTRIC_WALL;
                  spt->tpt->mn[lbl.m][11] = ELECTRIC_WALL;
# endif
               };

               if (( h == ONE )&&( tgl.ab[i] == 'm' )) 
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][2] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][2] = MAGNETIC_WALL;
                  spt->tpt->mn[lbl.m][11] = MAGNETIC_WALL;
# endif
               };

               if (( i == q )&&( tgl.bc[2*h] == 'e' ))      /* face 3 */
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][3] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][0] = ELECTRIC_WALL;
                  spt->tpt->mn[lbl.m][9] = ELECTRIC_WALL;
# endif
               };

               if (( i == q )&&( tgl.bc[2*h] == 'm' ))  
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][3] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][0] = MAGNETIC_WALL;
                  spt->tpt->mn[lbl.m][9] = MAGNETIC_WALL;
# endif
               };

/* top, bottom1: */

               if ( tgl.bot == 'e' )                        /* face 4 */
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][4] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][3] = ELECTRIC_WALL;
                  spt->tpt->mn[lbl.m][6] = ELECTRIC_WALL;
# endif
               };

               if ( tgl.bot == 'm' )                    
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][4] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][3] = MAGNETIC_WALL;
                  spt->tpt->mn[lbl.m][6] = MAGNETIC_WALL;
# endif
               };

               if ( tgl.top == 'e' )                       /* face 5 */
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][5] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][1] = ELECTRIC_WALL;
                  spt->tpt->mn[lbl.m][4] = ELECTRIC_WALL;
# endif
               };

               if ( tgl.top == 'm' )
               {
# if DSC_ADJCLL == 0
                  spt->tpt->mn[lbl.m][5] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                  spt->tpt->mn[lbl.m][1] = MAGNETIC_WALL;
                  spt->tpt->mn[lbl.m][4] = MAGNETIC_WALL;
# endif
               };
            };/* end if (( tgl.trv != 'e' )&&( tgl.trv != 'm' )) */ 
         };/* end if ( *option == 't'opology ) */

/* set trivial node material parameters to ZERO: */

         if ( *option == 'c' ) 
         {
            if (( tgl.trv == 'e' )||( tgl.trv == 'm' ))
               spt->ppt->mpt->idx[lbl.m] = null;

# if TGL_PRINT != null 
            if ( TGL_PRINT == lbl.m )
            {
               printf( "\n\n Intermediate values from function"
                       " '%s' ", __func__ );
               printf( "\n [ generated cell no. %ld ]: " , lbl.m );
               printf( "\n\n j = %ld : ( %.5e , %.5e , %.5e ) ", j,
                  spt->ppt->cpt->c[j][0],
                  spt->ppt->cpt->c[j][1],
                  spt->ppt->cpt->c[j][2] );

               printf( "\n k = %ld : ( %.5e , %.5e , %.5e ) ", k,
                  spt->ppt->cpt->c[k][0],
                  spt->ppt->cpt->c[k][1],
                  spt->ppt->cpt->c[k][2] );

               printf( "\n l = %ld : ( %.5e , %.5e , %.5e ) ", l,
                  spt->ppt->cpt->c[l][0],
                  spt->ppt->cpt->c[l][1],
                  spt->ppt->cpt->c[l][2] );

# if TGL_TRIANGLES == 1

               printf( "\n n = %ld : ( %.5e , %.5e , %.5e ) ", n,
                  spt->ppt->cpt->c[n][0],
                  spt->ppt->cpt->c[n][1],
                  spt->ppt->cpt->c[n][2] );
# else
               if( i < q )
               {
                  printf( "\n m = %ld : ( %.5e , %.5e , %.5e ) ", m,
                     spt->ppt->cpt->c[m][0],
                     spt->ppt->cpt->c[m][1],
                     spt->ppt->cpt->c[m][2] );
               }
               else if( i == q )
               {
                  printf( "\n n = %ld : ( %.5e , %.5e , %.5e ) ", n,
                     spt->ppt->cpt->c[n][0],
                     spt->ppt->cpt->c[n][1],
                     spt->ppt->cpt->c[n][2] );
               };
# endif
            };
# endif
         };/* end if *option == 'c'oordinates */ 

# if SET_TRIANGLES == 1

         if ( i < q )
         {
            SET_NODE( );

            if ( *option == 't' )
            { 
               spt->tpt->cm[lbl.m][0] = l;
               spt->tpt->cm[lbl.m][1] = n;
               spt->tpt->cm[lbl.m][2] = m;
               spt->tpt->cm[lbl.m][3] = k;

               ind = null - ( tgl.trv == 'e' ) - 2*( tgl.trv == 'm' );

# if DSC_ADJCLL == 1
               for ( o=null; o<=FACES ; o++ )  /* clear neighb. face identf. */
                  spt->tpt->mn[lbl.m][o] = ind;
# elif DSC_ADJCLL == 1
               for ( o=null ; o<=PORTS ; o++ ) /* clear neighb. port identf. */
                  spt->tpt->mn[lbl.m][o] = ind;
# endif
               if ( ind == null )
               {
                /* spt->tpt->bottom2: */

                  if ( tgl.bot == 'e' )                    /* face 4 */
                  {
# if DSC_ADJCLL == 1
                     spt->tpt->mn[lbl.m][4] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                     spt->tpt->mn[lbl.m][3] = ELECTRIC_WALL;
                     spt->tpt->mn[lbl.m][6] = ELECTRIC_WALL;
# endif
                  };

                  if ( tgl.bot == 'm' )
                  {
# if DSC_ADJCLL == 1
                     spt->tpt->mn[lbl.m][4] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                     spt->tpt->mn[lbl.m][3] = MAGNETIC_WALL;
                     spt->tpt->mn[lbl.m][6] = MAGNETIC_WALL;
# endif
                  };

                  if ( tgl.top == 'e' )                    /* face 5 */
                  {
# if DSC_ADJCLL == 1
                     spt->tpt->mn[lbl.m][5] = ELECTRIC_WALL;
# elif DSC_ADJCLL == 1
                     spt->tpt->mn[lbl.m][1] = ELECTRIC_WALL;
                     spt->tpt->mn[lbl.m][4] = ELECTRIC_WALL;
# endif
                  };

                  if ( tgl.top == 'm' )
                  {
# if DSC_ADJCLL == 1
                     spt->tpt->mn[lbl.m][5] = MAGNETIC_WALL;
# elif DSC_ADJCLL == 1
                     spt->tpt->mn[lbl.m][1] = MAGNETIC_WALL;
                     spt->tpt->mn[lbl.m][4] = MAGNETIC_WALL;
# endif
                  };
               };/* end if (( tgl.trv != 'e' )&&( tgl.trv != 'm' )) */
            };/* end if *option == 't'opology */ 

            if ( *option == 'c' )
            {
               if (( tgl.trv == 'e' )||( tgl.trv == 'm' ))
                  spt->ppt->mpt->idx[lbl.m] = null;

# if TGL_PRINT != null 
               if ( TGL_PRINT == lbl.m )
               {
                  printf( "\n\n Intermediate values from function"
                          " '%s' ", __func__ );
                  printf( "\n [ generated cell no. %ld ]: " , lbl.m );
                  printf( "\n\n k = %ld : ( %.5e , %.5e , %.5e ) ", k,
                     spt->ppt->cpt->c[k][0],
                     spt->ppt->cpt->c[k][1],
                     spt->ppt->cpt->c[k][2] );

                  printf( "\n l = %ld : ( %.5e , %.5e , %.5e ) ", l,
                     spt->ppt->cpt->c[l][0],
                     spt->ppt->cpt->c[l][1],
                     spt->ppt->cpt->c[l][2] );

                  printf( "\n m = %ld : ( %.5e , %.5e , %.5e ) ", m,
                     spt->ppt->cpt->c[m][0],
                     spt->ppt->cpt->c[m][1],
                     spt->ppt->cpt->c[m][2] );

                  printf( "\n n = %ld : ( %.5e , %.5e , %.5e ) ", n,
                     spt->ppt->cpt->c[n][0],
                     spt->ppt->cpt->c[n][1],
                     spt->ppt->cpt->c[n][2] );

                  printf( "\n\n [ please acknowledge ] " );
                  scanf( "%s",ptr);
              };
# endif
            };/* end if ( *option == 'c'oordinates */
         };/* end if ( i < q ) */

# endif  /* # endif SET_TRIANGLES == 1 */
      
         if ( h == ONE )
         {
            tgl.sab[i-ONE] = j;
            tgl.sab[i]     = k;
         };
         if ( i == q )
         {
            tgl.sbc[2*(h-ONE)] = k;
            tgl.sbc[2*h-ONE]   = n;
            tgl.sbc[2*h]       = l;
         };
         if ( i == ONE )
         {
            tgl.sac[h-ONE] = j;
            tgl.sac[h]     = l;
         };

         rp[i-ONE] = l;/* overtake vertex point [label 2] for next column h+1 */

      };/* next i */
   };/* next h */ 

/* reset wall and point indicators, boundaries: */

   for ( h=null ; h<=p ; h++ )
   {
      tgl.ab[h] = null; /* boundary ( wall ) indicators */
      tgl.ac[h] = null;
     tgl.cab[h] = null; /* boundary point transfer index */
     tgl.cac[h] = null;
   };
   for ( h=null ; h<=2*p ; h++ )
   {
      tgl.bc[h] = null; /* boundary (- wall ) indicators */
     tgl.cbc[h] = null; /* boundary point transfer index  */
   };

   tgl.aab = ZERO; /* reset angles */
   tgl.aac = ZERO;
   tgl.abc = ZERO;

   tgl.bot = null; /* reset wall indicators, */
   tgl.top = null; /* top and bottom         */
   tgl.trv = null; /* trivial node indicator */

   if ( lbl.pf < lbl.p )
      lbl.pf = lbl.p;

   if ( lbl.mf < lbl.m )
      lbl.mf = lbl.m;
           
   lbl.cf = lbl.m;

   return ONE;
}
/*============================================================================*/
# undef ELECTRIC_WALL
# undef MAGNETIC_WALL
# undef AB_POINTS
# undef SET_TRIANGLES
# undef SET_CORNER
# undef SET_NODE
# undef TGL_PRINT
/************************ end of function trngl(*) ****************************/
