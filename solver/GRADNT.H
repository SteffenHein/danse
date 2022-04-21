/* [ file: GRADNT.H ] */
/*******************************************************************************
*                                                                              *
*   MACRO GRADND                                                               *
*                                                                              *
*   This macro returns the nodal gradient, as well as the node vector          *
*   finite differences [representing the gradient in canonical coordinates]    *
*   of any quantity within a CLUSTER structure, pointed to by CLS *CLUSTER.    *
*   HH denotes the pertinent (heat&fluid) s-parameter set.                     *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: April 13, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define GRADND( CLS, HH ) \
{ \
/*...........................................................................*/\
/* restore gradient, */\
/* in cell coordinates [ with respect to basis b[] ] */\
/* db[k] := < b[k] | grad(f) >, */\
/* and in canonical coordinates */\
/* df[k] := < e[k] | grad(f) > = adj(B^-1 )*db ] */\
/*...........................................................................*/\
   ((CLS)->hh ) = (HH); \
\
   kk = null; do \
   { \
      ((CLS)->db[kk] ) = ((CLS)->face[TWO*kk+ONE] ) - ((CLS)->face[TWO*kk] );\
   } while(( ++kk ) < THREE );\
 \
   kk = null; do \
   { \
      ((CLS)->grad[kk] ) = ZERO; \
      ll = null; do \
      { \
         ((CLS)->grad[kk] ) += \
           (( hsp->bi[(HH)][kk][ll] )*((CLS)->db[ll] )); \
      } while(( ++ll ) < THREE ); \
   } while(( ++kk ) < THREE ); \
/*..........................................................................*/\
}
/*******************************************************************************
*                                                                              *
*   MACRO GRADFC                                                               *
*                                                                              *
*   This macro returns the face gradients of any substance in structure        *
*   CLUSTER, pointed to by CLS *CLUSTER, HH being the index of the             *
*   pertinent (heat & fluid) S-parameter set.                                  *
*                                                                              *
*   (C) SHEIN; Munich, February 2007                          Steffen Hein     *
*   [ Update: April 13, 2022 ]                             <contact@sfenx.de>  *
*                                                                              *
*******************************************************************************/
# define GRADFC( CLS, HH ) \
{ \
/*...........................................................................*/\
/* restore face temperature gradients */\
/*   dtdf[k] = adj(B^-1)*dtdb */\
/*...........................................................................*/\
   ((CLS)->hh ) = (HH); \
\
   sgn = ONE; \
   ((CLS)->dvgr ) = ZERO; \
   fc = null; do \
   { \
      kk = null; do \
      { \
         ll = ( short )( fc/TWO ); \
         ((CLS)->grdf[fc][kk] ) = 2.*sgn*\
            ( hsp->bi[(HH)][kk][ll] )*(((CLS)->node ) - ((CLS)->face[fc] ));\
 \
         ++ll; ll%=THREE; /* ll := ( ll+ONE ) mod THREE */ \
         ((CLS)->grdf[fc][kk] ) +=\
            (( hsp->bi[(HH)][kk][ll] )*((CLS)->db[ll] )); \
         ++ll; ll%=THREE; \
         ((CLS)->grdf[fc][kk] ) +=\
            (( hsp->bi[(HH)][kk][ll] )*((CLS)->db[ll] )); \
/* div grad: */\
         ((CLS)->dvgr ) +=\
            (((CLS)->grdf[fc][kk] )*( hsp->f[(HH)][fc][kk] )); \
      } while(( ++kk ) < THREE ); \
      sgn *= ( -ONE ); \
   } while(( ++fc ) < FACES ); \
}
/************************ end of macro GRADNT(*) ******************************/
/*
                A student of Nico Bourbaki,                                   
                and a teacher of math in Kentucky,                           
                thought the function of sex
                must be semi-convex
                but to prove it, he wasn't quite lucky.

                                  [ from 'Mathematical Intelligencer', 1998 ]
*/
/************************** end of file GRADNT.H ******************************/
