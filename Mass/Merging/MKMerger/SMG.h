#define MSTEPS 130
#define NSTEPS 130





typedef struct MBP{
	IDTYPE gid;
	IDTYPE hid[NSTEPS];
} MBP;
/*
typedef struct HID{
	IDTYPE hid,gid;
	float sigv;
} HID;
*/


/*
typedef struct MergerMBP{
	SavedHaloType mergerhistory[NSTEPS];
}MergerMBP;
*/

typedef struct MergeHistory{
    SavedHaloType step[MSTEPS];
} MergeHistory;

/*
typedef struct SortMergeHistory{
	IDTYPE majorgid;
	MergeHistory mergehist;
} SortMergeHistory;

typedef struct SSortMergeHistory{
	IDTYPE majorgid,gid,subid, newgid, smajorgid;
	IDTYPE stepmajorgid[NSTEPS];
} SSortMergeHistory;
*/


typedef struct MajorGid{
	IDTYPE gid, newgid, mgid;
	IDTYPE majorgid[NSTEPS];
} MajorGid;

typedef struct GID{
	IDTYPE gid, newgid;
} GID;



typedef struct SigV{
	float sigv[NSTEPS];
}SigV;
typedef struct GidMergeHistory{
	IDTYPE gid;
	/*
    MergeHistory mergerhistory;
	*/
	SavedHaloType step[NSTEPS];
	SigV sigv;
} GidMergeHistory;
