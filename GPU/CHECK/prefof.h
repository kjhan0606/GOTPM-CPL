#define NBitsperByte 8
#define FoF_RESET 0x00
char FoFBITs[NBitsperByte]={0x01, 0x02,0x04 , 0x08,0x10,0x20,0x40,0x80};
#define SET_FoFflag(a,n)    (fofflag[a] |=   FoFBITs[n] )
#define UNSET_FoFflag(a,n)  (fofflag[a] &= (~FoFBITs[n]))
#define IS_FoFflag(a,n)     (fofflag[a] &    FoFBITs[n] )
#define TOGGLE_FoFflag(a,n) (fofflag[a] ^=   FoFBITs[n] )



#ifdef DEFINE_PRE_FOF
#define EXTERN
#else
#define EXTERN extern
#endif
EXTERN unsigned char *fofflag;
EXTERN int nfofflag;
#undef EXTERN
