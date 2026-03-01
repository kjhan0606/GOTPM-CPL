#define NSPLINE (1<<12)
typedef double FORCEARRAYTYPE;
typedef struct DIFFFORCE{
	FORCEARRAYTYPE diff[3];
	FORCEARRAYTYPE slope[3];
} DIFFFORCE;


/* These are the GPU settings.                           */
/* ##########     WARNING       ########################
 * ##########     This should be calculated manually. ##
 * ##################################################### */
#define RMAX 6.
#define RMAXSQ 36.
#define LOG10RMAX (0.7781512) /* = alog10(RMAX) */
#define LOG10RMIN (-3.0) /* log10(RMIN) */
#define RMINSQ (1.E-6) /* == RMIN*RMIN */
#define LOG10RMAXmLOG10RMIN (3.7781512) /* == log10(RMAX)-log10(RMIN) */
#define invLOG10RMAXmLOG10RMIN (0.264680) /* == 1./(log10(RMAX)-log10(RMIN)) */



#ifdef DEFINE_TREE_COMMON
#define EXTERN
#else
#define EXTERN extern
#endif
/*
EXTERN  DIFFFORCE forcecorrect[NSPLINE];
#define (forcecorrectdiff(i,j)) forcecorrect[i].diff[j]
#define (forcecorrectslope(i,j)) forcecorrect[i].slope[j]
*/
EXTERN  FORCEARRAYTYPE diff[NSPLINE][3],slope[NSPLINE][3];
#define forcecorrectdiff(a,b) diff[a][b]
#define forcecorrectslope(a,b) slope[a][b]

EXTERN  TREEPRECISION ran2nran,invran2nran;
#undef EXTERN





