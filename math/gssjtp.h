/* [ file: gssjtp.h ] */
/* The structure type definition of Gauss-Jordan algorithm functions          */
/* gssjrd(*) [line pivoting] and gssjpv(*) [full pivoting ]                   */
/* Update: May 01, 2007                                                       */
/*----------------------------------------------------------------------------*/
# ifndef GSS_DISP
   # define GSS_DISP  1     /* 1: display if coeff. mtx.gss.m[][] is singular */
# endif
# ifndef GSS_IMPLCT
   # define GSS_IMPLCT 1    /* 1: implicit pivoting  [ normalizing equations ]*/
# endif
# ifndef GSS_STRUCT
   # define GSS_STRUCT 1    /* 1: define struct gaussc gss [may be done elsw.]*/
# endif
# ifndef GSS_LONG_DBL
   # define GSS_LONG_DBL 0  /* 1: long double prexcision */
# endif
# ifndef GSS_MAXRNK
   # define GSS_MAXRNK 10   /* maximum rank of matrix M */
# endif
# ifndef GSS_SINGLR 
   # define GSS_SINGLR (1.e-277) /* bound for singularity check */
# endif
/*----------------------------------------------------------------------------*/
typedef struct
{
   signed char
      rtn;

   char 
      opt;

   short
      rank, /* rank of matrix M = ( mr[][]+j*mi[][] ) */
      neqs; /* number of linear equations to be solved [ simultaneously ] */
            /* - only used in function gssjpv(*) */

   # if GSS_LONG_DBL  == 1
      long double
   # else
      double
   # endif

         dtr,                         /* determinant  ( real  part )          */
         dti,                         /* "            ( imag. part )          */
         mr[GSS_MAXRNK][2*GSS_MAXRNK],/* r = real part of matrix              */
         mi[GSS_MAXRNK][2*GSS_MAXRNK],/* i = imag."    "  "                   */
         zr[GSS_MAXRNK][GSS_MAXRNK],
         zi[GSS_MAXRNK][GSS_MAXRNK];

} GAUSS_JRD;
/*-------------------------- end of file gssjtp.h ----------------------------*/
