#ifdef XYZDBL
typedef double POSTYPE;
#else
typedef float POSTYPE;
#endif


typedef struct SavedHaloType{
	INDXTYPE mbp; /* mbp ID */
	IDTYPE nowhid,majorglobalid; /* nowhid is halo id at the redshift and majorglobalid is global array id to the major subhalo
									(or mbp ) */
	POSTYPE x,y,z,hx,hy,hz; /* in unit of h^-1 Mpc in comoving coordinate 
							 (x,y,z) is position of the mbp and (hx,hy,hz) is the CoM position of the hosting FoF halo */
	float mass,fofhmass; /* in unit of M_sun/h 
						  mass is the mass of MBP or halo (if one mbp is in the halo) just before entering larger halo region.
						  fofhmass is the hosting halo mass. */
	float vx,vy,vz,hvx,hvy,hvz; /* in unit of km/second 
								 (vx,vy,vz) is the peculiar velocity of the mbp &
								 (hvx,hvy,hvz) is the peculiar velocity of the hosting FoF halo. */
	short int statusflag, linkflag;
} SavedHaloType;



typedef struct SimplifiedHalo{
	SavedHaloType basichalo;
	IDTYPE nowaid,downaid,upaid,globalid,nowhid;
} SimplifiedHalo;
