#define MAXHALO 1000000
#define MAXTREEHALO 7
#define NMINMEMBERS 30


typedef long INDXTYPE;
typedef long IDTYPE;
#define MPI_IDTYPE MPI_LONG
#define MPI_INDXTYPE MPI_LONG
typedef struct HaloProp{
	INDXTYPE hid,mbp;
	float mass,vx,vy,vz;
#ifdef XYZDBL
	double x,y,z;
#else
	float x,y,z;
#endif
	float paraspin,sigv,rotation,radius,q,s,kw;
	float ang_mom[3],v[3][3];
} HaloProp;
typedef HaloProp HaloQ;


/*
 * mbp : id of most bound particle of a halo
 */
typedef struct TrimmedHalo {
	INDXTYPE mbp;
	IDTYPE aid; /* neded to sort mbp array */
	IDTYPE nexthid,nowhid,downaid,upaid;
	float mass,fofhmass;
	short int statusflag,linkflag;
}TrHalo;



typedef struct MbpRV{
#ifdef XYZDBL
	double x,y,z;
#else
	float x,y,z;
#endif
	float vx,vy,vz;
} MbpRV;



typedef struct UpLink{
	INDXTYPE mbp;
	IDTYPE aid;
	IDTYPE nexthid,nowhid,upaid;
	float mass;
	short int linkflag;
} UpLink;

typedef struct basicparticletype {
#ifdef XYZDBL
	double x,y,z;
#else
	float x,y,z;
#endif
	float vx,vy,vz;

#ifdef INDEX
	INDXTYPE indx;
#endif
} basicparticletype;




typedef struct Halo{
	    size_t np;
#ifdef XYZDBL
	    double x,y,z;
#else
	    float x,y,z;
#endif
	    float vx,vy,vz;
} Halo;




typedef struct HidMass{
	IDTYPE hid;
	float mass;
} HidMass;

#define TERMINATE -999


#define MP_RESET 0x00
#define MP_ALONE_SUB  0x01
#define MP_NEW_BORN 0x02
#define MP_ONE_OF_MULTI_SUB 0x04
#define MP_MAJOR_SATELLITE 0x08
#define MP_MINOR_SATELLITE 0x10
#define MP_TERMINATED 0x20
#define MP_DISSOLVED 0x40

#define MP_SINGLE_LINE 0x01
#define MP_MULTI_LINE 0x02
#define MP_MAJOR_LINE 0x04
#define MP_MINOR_LINE 0x08
#define MP_MISSING_LINE 0x10

#define RESET_WHOLE_FLAGS_MP(halo,ii,flag) (halo[ii].flag  = MP_RESET)



#define SET_ALONE_SUB(halo,ii,flag) (halo[ii].flag |= MP_ALONE_SUB)
#define UNSET_ALONE_SUB(halo,ii,flag) (halo[ii].flag &= (~MP_ALONE_SUB))
#define IS_ALONE_SUB(halo,ii,flag) (halo[ii].flag & MP_ALONE_SUB)
#define TOGGLE_ALONE_SUB(halo,ii,flag) (halo[ii].flag ^= MP_ALONE_SUB)

#define SET_SINGLE_LINE(halo,ii,flag) (halo[ii].flag |= MP_SINGLE_LINE)
#define UNSET_SINGLE_LINE(halo,ii,flag) (halo[ii].flag &= (~MP_SINGLE_LINE))
#define IS_SINGLE_LINE(halo,ii,flag) (halo[ii].flag & MP_SINGLE_LINE)
#define TOGGLE_SINGLE_LINE(halo,ii,flag) (halo[ii].flag ^= MP_SINGLE_LINE)

#define SET_NEW_BORN(halo,ii,flag) (halo[ii].flag |= MP_NEW_BORN)
#define UNSET_NEW_BORN(halo,ii,flag) (halo[ii].flag &= (~MP_NEW_BORN))
#define IS_NEW_BORN(halo,ii,flag) (halo[ii].flag & MP_NEW_BORN)
#define TOGGLE_NEW_BORN(halo,ii,flag) (halo[ii].flag ^= MP_NEW_BORN)

#define SET_ONE_OF_MULTI_SUB(halo,ii,flag) (halo[ii].flag |= MP_ONE_OF_MULTI_SUB)
#define UNSET_ONE_OF_MULTI_SUB(halo,ii,flag) (halo[ii].flag &= (~MP_ONE_OF_MULTI_SUB))
#define IS_ONE_OF_MULTI_SUB(halo,ii,flag) (halo[ii].flag & MP_ONE_OF_MULTI_SUB)
#define TOGGLE_ONE_OF_MULTI_SUB(halo,ii,flag) (halo[ii].flag ^= MP_ONE_OF_MULTI_SUB)

#define SET_MULTI_LINE(halo,ii,flag) (halo[ii].flag |= MP_MULTI_LINE)
#define UNSET_MULTI_LINE(halo,ii,flag) (halo[ii].flag &= (~MP_MULTI_LINE))
#define IS_MULTI_LINE(halo,ii,flag) (halo[ii].flag & MP_MULTI_LINE)
#define TOGGLE_MULTI_LINE(halo,ii,flag) (halo[ii].flag ^= MP_MULTI_LINE)

#define SET_MAJOR_LINE(halo,ii,flag) (halo[ii].flag |= MP_MAJOR_LINE)
#define UNSET_MAJOR_LINE(halo,ii,flag) (halo[ii].flag &= (~MP_MAJOR_LINE))
#define IS_MAJOR_LINE(halo,ii,flag) (halo[ii].flag & MP_MAJOR_LINE)
#define TOGGLE_MAJOR_LINE(halo,ii,flag) (halo[ii].flag ^= MP_MAJOR_LINE)

#define SET_MINOR_LINE(halo,ii,flag) (halo[ii].flag |= MP_MINOR_LINE)
#define UNSET_MINOR_LINE(halo,ii,flag) (halo[ii].flag &= (~MP_MINOR_LINE))
#define IS_MINOR_LINE(halo,ii,flag) (halo[ii].flag & MP_MINOR_LINE)
#define TOGGLE_MINOR_LINE(halo,ii,flag) (halo[ii].flag ^= MP_MINOR_LINE)


#define SET_MISSING_LINE(halo,ii,flag) (halo[ii].flag |= MP_MISSING_LINE)
#define UNSET_MISSING_LINE(halo,ii,flag) (halo[ii].flag &= (~MP_MISSING_LINE))
#define IS_MISSING_LINE(halo,ii,flag) (halo[ii].flag & MP_MISSING_LINE)
#define TOGGLE_MISSING_LINE(halo,ii,flag) (halo[ii].flag ^= MP_MISSING_LINE)


#define SET_MAJOR_SATELLITE(halo,ii,flag) (halo[ii].flag |= MP_MAJOR_SATELLITE)
#define UNSET_MAJOR_SATELLITE(halo,ii,flag) (halo[ii].flag &= (~MP_MAJOR_SATELLITE))
#define IS_MAJOR_SATELLITE(halo,ii,flag) (halo[ii].flag & MP_MAJOR_SATELLITE)
#define TOGGLE_MAJOR_SATELLITE(halo,ii,flag) (halo[ii].flag ^= MP_MAJOR_SATELLITE)

#define SET_MINOR_SATELLITE(halo,ii,flag) (halo[ii].flag |= MP_MINOR_SATELLITE)
#define UNSET_MINOR_SATELLITE(halo,ii,flag) (halo[ii].flag &= (~MP_MINOR_SATELLITE))
#define IS_MINOR_SATELLITE(halo,ii,flag) (halo[ii].flag & MP_MINOR_SATELLITE)
#define TOGGLE_MINOR_SATELLITE(halo,ii,flag) (halo[ii].flag ^= MP_MINOR_SATELLITE)

#define SET_TERMINATED(halo,ii,flag) (halo[ii].flag |= MP_TERMINATED)
#define UNSET_TERMINATED(halo,ii,flag) (halo[ii].flag &= (~MP_TERMINATED))
#define IS_TERMINATED(halo,ii,flag) (halo[ii].flag & MP_TERMINATED)
#define TOGGLE_TERMINATED(halo,ii,flag) (halo[ii].flag ^= MP_TERMINATED)

#define SET_DISSOLVED(halo,ii,flag) (halo[ii].flag |= MP_DISSOLVED)
#define UNSET_DISSOLVED(halo,ii,flag) (halo[ii].flag &= (~MP_DISSOLVED))
#define IS_DISSOLVED(halo,ii,flag) (halo[ii].flag & MP_DISSOLVED)
#define TOGGLE_DISSOLVED(halo,ii,flag) (halo[ii].flag ^= MP_DISSOLVED)


