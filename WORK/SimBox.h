#define NDim 3
typedef struct Range{
	float min,max;
} Range;
typedef struct Box{
	Range r[NDim];
} Box;

