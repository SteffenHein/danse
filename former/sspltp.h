/* [ file: sspltp.h ] */
/* Update: May 10, 2007 */
/*----------------------------------------------------------------------------*/
/* Macro and structure type definition of the [ DSC mesh cell system ] plot */
/* function sysplt(*) */
/*----------------------------------------------------------------------------*/
# ifdef USER_PATH
   # define SPL_PATH USER_PATH
# else
   # define SPL_PATH ""
# endif
/*............................................................................*/
# define SPL_PREFIX "gnu."
/*............................................................................*/
# define SPL_DISPLAY 1 /* 0<N: display plot sytem file creation messages */ 
                       /* 2:   display plot data file creation messages */
/*............................................................................*/
/* SPL_PLOTMEDIA - that maximum number of media */
# ifndef SPL_PLOTMEDIA
   # define SPL_PLOTMEDIA 10 /* K: plot that maximum number of media */
# endif /* not defined SPL_PLOTMEDIA */
/*............................................................................*/
/* SPL_LAYERS - that maximum number of plotted layers */
# ifndef SPL_LAYERS
   # define SPL_LAYERS 300
# endif /* not defined SPL_LAYERS */
/*............................................................................*/
/* SPL_PLOTBOTTM = [0] 1: [don't] plot walls on cell bottom */
# ifndef SPL_PLOTBOTTM 
   # define SPL_PLOTBOTTM 0
# endif /* not defined SPL_PLOTBOTTM */
/*............................................................................*/
/* SPL_PLOTTOPS = [0] 1: [don't] plot walls on cell tops */
# ifndef SPL_PLOTTOPS
   # define SPL_PLOTTOPS 0
# endif /* not defined SPL_PLOTTOPS */
/*............................................................................*/
/* SPL_PLOTTRIVL = [0] 1: [don't] plot trivial cells */
# ifndef SPL_PLOTTRIVL
   # define SPL_PLOTTRIVL 0
# endif /* not defined SPL_PLOTTRIVL */
/*............................................................................*/
/* SPL_PLOTEWLLS = 0/1/2: don't/plot/plot_dominant electric walls */
# ifndef SPL_PLOTEWLLS
   # define SPL_PLOTEWLLS 1
# endif /* not defined SPL_PLOTEWLLS */
/*............................................................................*/
/* SPL_PLOTMWLLS = 0/1/2: don't/plot/plot_dominant magnetic walls */
# ifndef SPL_PLOTMWLLS
   # define SPL_PLOTMWLLS 1
# endif /* not defined SPL_PLOTMWLLS */
/*............................................................................*/
/* SPL_RESCLE_XY_ = [0] 1: [don't] set uniform scales on xy axes */
# ifndef SPL_RESCLE_XY_
   # define SPL_RESCLE_XY_ 1 /* [0] 1: [don't] set uniform scales on xy axes */
# endif /* not defined SPL_RESCLE_XY_ */
/*............................................................................*/
/* SPL_RESCLE_XYZ = [0] 1: [don't] set uniform scales on all axes */
# ifndef SPL_RESCLE_XYZ
   # define SPL_RESCLE_XYZ 1
# endif /* not defined SPL_RESCLE_XYZ */
/*............................................................................*/
# define PLOT2_FILE "2p"
# define PLOT3_FILE "3p"
# define CELL2_FILE "2mesh"
# define CELL3_FILE "3mesh"

# define MED12_FILE "2trvl"
# define MED13_FILE "3trvl"
# define MED22_FILE "2epsr"
# define MED23_FILE "3epsr"

# define EWLL2_FILE "2ewll"
# define EWLL3_FILE "3ewll"
# define MWLL2_FILE "2mwll"
# define MWLL3_FILE "3mwll"
/*............................................................................*/
# ifndef ELECTRIC_WALL
   # define ELECTRIC_WALL -1
# endif
/*............................................................................*/
# ifndef MAGNETIC_WALL
   # define MAGNETIC_WALL -2
# endif
/*............................................................................*/
/* The structure [ type definition ] of DSC system plot function sysplt(*) */
typedef struct
{
   short
      rtn;

   short
      media,
      layer,
      toplr, /* the top layer, in general equals blc.base[blc.m[null]] */
      lay[SPL_LAYERS+ONE];

   long
      mi, mf, pm;

   char
      flbl[THREE],
      p2[SHS_SIZE],
      p3[SHS_SIZE],
      c2[SHS_SIZE],
      c3[SHS_SIZE];

# if SPL_PLOTEWLLS != 0
   char
      ew2[SHS_SIZE],
      ew3[SHS_SIZE];
# endif
# if SPL_PLOTMWLLS != 0
   char
      mw2[SHS_SIZE],
      mw3[SHS_SIZE];
# endif
# if SPL_PLOTMEDIA != 0
   char
      md2[SPL_PLOTMEDIA][SHS_SIZE],
      md3[SPL_PLOTMEDIA][SHS_SIZE];

   short
      mdidx[SPL_PLOTMEDIA];  /* media index identifier: labels */
                             /* indices of media to be plotted */
# endif

   TOPOLOGY
     *tpt;

   PARAMETERS
     *ppt;

} SYSPLT;
/*************************** end of file sspltp.h *****************************/
