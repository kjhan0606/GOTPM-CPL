#define NSTEPS 130
typedef struct MBP{
	IDTYPE gid;
	IDTYPE hid[NSTEPS];
} MBP;
typedef struct HID{
	IDTYPE hid,gid;
	float sigv;
} HID;


typedef struct MergerMBP{
	SavedHaloType mergerhistory[NSTEPS];
}MergerMBP;
