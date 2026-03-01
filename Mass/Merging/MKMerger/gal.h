#include "merger.h"
typedef struct RedInfo{
	float z_form,z_merge,z_end;
} RedInfo;

typedef struct Mhistory{
	float red,NewMstar;
}Mhistory
typedef struct TrTrHalo{
	INDXTYPE mbp;
	IDTYPE nowhid;
	float mass;
	short int statusflag,linkflag;
}TrTrHalo;

typedef struct GalModel{
	TrTrHalo merginglink;
	RedInfo redinfo;
	size_t imhistory;
}GalModel;

#define dumpTr2TrTr(a,b) do{\
	a->merginglink.mbp = b->mbp;\
	a->merginglink.nowhid = b->nowhid;\
	a->merginglink.mass = b->mass;\
	a->merginglink.statusflag = b->statusflag;\
	a->merginglink.linkflag = b->linkflag;\
} while(0)


typedef struct treeparameters{
} treeparameters;
