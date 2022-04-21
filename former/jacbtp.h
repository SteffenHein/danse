/* [ file: jacbtp.h ] */
/* The structure type definition of Jacobi eigenvalue function jacobi(*)      */
/* Update April 13, 2007                                                      */
/*----------------------------------------------------------------------------*/
# ifndef JAC_MAXRNK
   # define JAC_MAXRNK 10 /* maximum order of matrix H = jac.h.r + j*jac.h.i  */
# endif
# ifndef JAC_MAXITR
   # define JAC_MAXITR 1000 /* maximum number of iterations                   */
# endif
# ifndef JAC_BOUND1
   # define JAC_BOUND1 (1.0e-277) /*                                          */
# endif
# ifndef JAC_BOUND2
   # define JAC_BOUND2 (1.3e-27) /* changed from 5.3e-28:  09-07-1999         */
# endif                          /* [ which proved to be too restrictive      */
                                 /*               for DSC model 'mod_rj99d' ] */
# ifndef JAC_RNDOFF
   # define JAC_RNDOFF (1.0e-16) /* roundoff limit for eigenvalues jac.h.r,i  */
# endif
/*----------------------------------------------------------------------------*/
struct cmatrix  
{
   double 
      r[JAC_MAXRNK][JAC_MAXRNK], i[JAC_MAXRNK][JAC_MAXRNK];
};
/*----------------------------------------------------------------------------*/
typedef struct
{
   signed char
      rtn;

   double
      max, bnd1, bnd2, mod, phi, hn, c, s,
      sn[JAC_MAXRNK];

   short
      rank, i1, i2;

   unsigned char
      skew;

   struct cmatrix
      h, e;

}  JACOBI_EV;
/************************** end of file jacbtp.h *****************************/
