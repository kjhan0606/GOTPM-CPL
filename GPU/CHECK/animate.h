typedef struct ImageSize{
	int nx,ny;
} ImageSize;
#define SonyDIgitalProjection	{4096,2160}
#define HDTV720p				{1280,720}
#define HDTV1080i				{1440,1080}
#define HDTV1080p				{1920,1080}
#define Standard				{1024,768}
#define SDTV480i				{720,480}
#define SDTV576i				{720,576}
#define UXGA					{1600,1200}
#define WSXGA					{1440,900}
#define Planetarium				{3600,3600}

#define MAXFRAME 10000

typedef struct ThreeD{
	double x,y,z;
} ThreeD;
typedef struct EulerAgle{
	double alpha,beta,gamma;
} EulerAngle;
typedef struct View{
	double dTheta, dPhi;
} View;
typedef struct Viewer{
	/* frame is the frame number & nstep is the simulation stepcount */
	int frame, nstep;
	ThreeD pos,view;
	ThreeD E1,E2,E3;
	double alpha,beta,gamma;
	View wangview;
	float HOpticalDepth;
	char outfile[30];
} Viewer;
Viewer *viewer;



void ObsAnimate(pmparticletype *,int, Viewer *, int, int, int );
