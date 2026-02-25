#ifdef XYZDBLE
#ifdef INDEX
#define XofP(x,indx,mx,rngx) amod((x)+mod((indx),(mx))+rngx,rngx)
#define YofP(y,indx,mx,my,rngy) amod((y)+mod((indx),((mx)*(my)))/(mx)+rngy,rngy)
#define ZofP(z,indx,mx,my,rngz) amod((z)+   ((indx)/((mx)*(my)))+rngz,rngz)
#else
#error the pair definition of XYZDBL and !INDEX  is not allowed in the current version.
#endif
#endif
      integer*8 plan,iplan
