/*
 * 이 코드는 찾아낸 헤일로의 peak density 를 구하는 것이다.
 * 여기에 헤일로의 물리적 성질을 첨가할것이다.
 */
/*
 * 지금 버젼은 trim_pden.c  를 수행해서 under dense peak 을 갖는 헤일로들을
 * 미리 제거하고 난 후를 가정하고 계산 한다. 22/01/2003
 * */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include "Memory.h"
double massscale,rscale;
typedef long indxtype;
typedef struct PtlPos{
    float x,y,z;
} PtlPos;
typedef struct PtlVel{
    float vx,vy,vz;
} PtlVel;
#include "defs.h"
void DenPeakCenter(vector3d *,int ,double *,double  *,double *, int);
#define MIN(A,B) ((A)<(B) ? (A): (B))
#define MAX(A,B) ((A)>(B) ? (A): (B))
#define sqr(A) ((A)*(A))
#define pow2(A) ((A)*(A))
#define pow3(A) ((A)*(A)*(A))
#define sqr(A) ((A)*(A))
double cm2km=1.E-5L;
double sim2real,real2sim;
double pntmass,r1kineticfact,r2kineticfact;
double potentfact;
double onesolarmass=1.989E33L;
double  G= 6.672E-8L;
double pc= 3.08567802E18L;
double pi= 3.1415926535L;
double rotation;
int nspace;
float r1[3],r2[3],r3[3],radius;
float a;
float amax;
float redshift;
float size;
float rng;
float omep,omepb;
float omeplam;
float hubble;
float omepk,Hsub;
int nx,ny,nz;
float *den;
void findden(int , PtlPos *, int ,int ,int );
void triaxialshape_(int *,float *,float *,float *,float *,
        float *,float *,float *);
float wc(float );
float find_r(PtlPos *, int, float ,float, float *, float *,float *);
double ang_mom[3];
float q,s;
float d[3],v[3][3];
double para_spin;
double tenergy,tmass,pvx;
double cx,cy,cz,cvx,cvy,cvz,pvr;
double sigv;
size_t Fread(void *buff,size_t size,size_t nmem,FILE *fp){
    char a,b,c,d;
    char *data;
    int i,j;
    size_t result;
    result = fread(buff,size,nmem,fp);
    data = (char *) buff;
    for(i=0;i<nmem*size;i+=4){
        a = data[i]; b = data[i+1]; c = data[i+2]; d = data[i+3];
        data[i] = d;
        data[i+1] = c;
        data[i+2] = b;
        data[i+3] = a;
    }
    return result;
}
/*
#define fread(a,b,c,d) Fread(a,b,c,d)
*/
#include "tree.h"
#include "Time.h"
#include "force_spline.h"
void halo_potent(int nend,PtlPos *r,float *penergy){
    TStruct *TREE;
    TPtlStruct *ptl;
    particle p;
    Box box;
    float theta=1.;
    int i,j,k;
    static int np;
    np = nend;
    TREE = (TStruct *) Malloc(sizeof(TStruct)*np*4,PPTR(TREE));
    ptl = (TPtlStruct *) Malloc(sizeof(TPtlStruct)*np,PPTR(ptl));
    for(i=0;i<np;i++){
        ptl[i].type = TYPE_PTL;
        ptl[i].r[0] = r[i].x; ptl[i].r[1] = r[i].y; ptl[i].r[2] = r[i].z;
        ptl[i].mass = 1;
    }

    box = findbox(ptl,np);
    Make_Tree(TREE,ptl,np,box);
    for(i=0;i<np;i++){
        p.x = ptl[i].r[0]; p.y = ptl[i].r[1]; p.z = ptl[i].r[2];
        penergy[i] = treeplumpotential(&p,theta,TREE,ptl)*potentfact;
        
    }
    Free(ptl);
    Free(TREE);
}
void halo_properties(int nend,PtlPos *r,PtlVel *vr){
    float distx,disty,distz,tmpvx,tmpvy,tmpvz;
    float distvx,distvy,distvz;
    float tmpvvx,tmpvvy,tmpvvz;
    float *kenergy,*penergy,*rv;
    float absdist,rvx,rvy,rvz;
    float *x,*y,*z;
    double mrv;
    static int np;
    int i,j,k;
    np = nend;
    cvx = cvy = cvz = 0.;
    for(i=0;i<np;i++){
        cvx += vr[i].vx; cvy += vr[i].vy; cvz += vr[i].vz;
    }
    cvx = cvx/ (float) np; cvy = cvy/(float) np; cvz = cvz/(float) np;
    kenergy = (float *) Malloc(sizeof(float)*np,PPTR(kenergy));
    penergy = (float *) Malloc(sizeof(float)*np,PPTR(penergy));
    rv = (float *) Malloc(sizeof(float)*np,PPTR(rv));
    for(i=0;i<np;i++){
        distx = (r[i].x-cx)*r1kineticfact;
        disty = (r[i].y-cy)*r1kineticfact;
        distz = (r[i].z-cz)*r1kineticfact;
        tmpvx = (vr[i].vx-cvx)*r2kineticfact;
        tmpvy = (vr[i].vy-cvy)*r2kineticfact;
        tmpvz = (vr[i].vz-cvz)*r2kineticfact;
        tmpvvx = (distx+tmpvx); tmpvvy = (disty+tmpvy); tmpvvz = (distz+tmpvz);
        kenergy[i] = 0.5*(tmpvvx*tmpvvx+tmpvvy*tmpvvy+tmpvvz*tmpvvz);
    }
    halo_potent(np,r,penergy);

    for(i=0;i<3;i++){
        ang_mom[i] = 0.L;
    } 
    mrv = 0.L;
    for(i=0;i<np;i++){
        distx = r[i].x - cx; disty = r[i].y - cy; distz = r[i].z - cz;
        distvx = vr[i].vx-cvx; distvy = vr[i].vy-cvy; distvz = vr[i].vz-cvz;
        ang_mom[0] += disty*distvz-distz*distvy;
        ang_mom[1] += distz*distvx-distx*distvz;
        ang_mom[2] += distx*distvy-disty*distvx;
        absdist = sqrt(distx*distx+disty*disty+distz*distz);
        rvx = disty*distvz-distz*distvy;
        rvy = distz*distvx-distx*distvz;
        rvz = distx*distvy-disty*distvx;
        rv[i] = sqrt(rvx*rvx+rvy*rvy+rvz*rvz)/absdist;
        mrv = MAX(mrv,rv[i]);
    }
    tenergy = 0.;
    for(i=0;i<np;i++) tenergy += penergy[i]*0.5 + kenergy[i];
    para_spin = sqrt(sqr(ang_mom[0])+sqr(ang_mom[1])+sqr(ang_mom[2])) *
        sqrt(-tenergy)/G/pntmass*sim2real*
        r2kineticfact/pow((double)np,2.5);
    /* angular momentum in the unit of M_odot kpc km/sec */
    rotation = sqrt(sqr(ang_mom[0])+sqr(ang_mom[1])+sqr(ang_mom[2])) *
        sim2real/(1.E3L*pc)*r2kineticfact*pntmass/onesolarmass;
    tmass = pntmass*np;
    {
        x = (float *) Malloc(sizeof(float)*np,PPTR(x));
        y = (float *) Malloc(sizeof(float)*np,PPTR(y));
        z = (float *) Malloc(sizeof(float)*np,PPTR(z));
        for(i=0;i<np;i++){
            x[i] = r[i].x; y[i] = r[i].y; z[i] = r[i].z;
        }
        triaxialshape_(&np,x,y,z,&q,&s,&d[0],&v[0][0]);
        Free(x);Free(y);Free(z);
    }
    Free(rv);Free(penergy);Free(kenergy);
}
void getparameter(FILE *wp){
    float bias,rsmooth;
    float npow0;
    float z_init,delta_a;
    float SPHERE_RADIUS,particle_radius,rtheta;
    fscanf(wp,"%f %f\n",&size,&hubble);
    fscanf(wp,"%f %f %f %f %f %f\n",&npow0,&omep,&omepb,&omeplam,&bias,&rsmooth);
    fscanf(wp,"%d %d %d %d\n",&nx,&ny,&nz, &nspace);
    fscanf(wp,"%f %f %f\n",&SPHERE_RADIUS,&particle_radius,&rtheta);
    fscanf(wp,"%f %f %f\n",&z_init,&delta_a,&a);
    rng = nx;
    amax = 1.+z_init;
    redshift = amax/a - 1.;
}

/*
int main(int argc, char **argv){
*/
void MAIN_(int argc, char **argv){
    FILE *fp,*wp,*paraid,*mergerwp;
    PtlPos *r;
    PtlVel *vr;
    float peak;
    float *x,*y,*z,*vx,*vy,*vz;
    int tnend,nend,i,j,k;
    int ncount;
    indxtype *indx,id;
    int mpeak,np;
    float xmin,ymin,xmax,ymax,width;
    int anend[10000];
    int haloindex;
    int index,peakindex;
    char infile[100],parafile[100];
    char wfile[100],mergerfile[100];
    int START,inid;
    long nsize;
    int halomergerindex;
    int *nowindex,*beforeindex;

    float tt1,tt2,tt3,tt4;
    int nindex,nmpeak,npeakindex;
    int nhalomergerindex,mtmp;
    int nnowhaloindex,nbeforehaloindex;
    int beforehalomergerindex;

    nowindex = beforeindex = NULL;

    if(argc != 2) {
        fprintf(stderr,"Error input : mergertree.exe stepnum\n");
        exit(99);
    }
    (void) Make_Total_Memory();
    i_potent_spline();
    den = (float *) Malloc(256*256*256*sizeof(int),PPTR(den));
    inid = atoi(argv[1]);
    nnowhaloindex = 0;
        sprintf(wfile,"n%.5dhalomass.dat",inid);
        sprintf(mergerfile,"n%.5dmerger.dat",inid);
        sprintf(infile,"HFIND.DATA.%.5d",inid);
        sprintf(parafile,"params.%.5d",inid);
        if((fp = fopen(infile,"rb")) == NULL||(paraid=fopen(parafile,"r"))==NULL){
            fprintf(stderr,"error opening file %s && %s\n",infile,parafile);
            exit(-999);
        }
        wp = fopen(wfile,"w");
        /*
        mergerwp = fopen(mergerfile,"w");
        */
        getparameter(paraid);
        fclose(paraid);
        nsize = rng*rng*rng/nspace/nspace/nspace;



        omepk = 1.-omep-omeplam;
        Hsub = sqrt(omep*pow3(amax/a)+omeplam+omepk*pow2(amax/a));
        r1kineticfact = size/hubble/rng/amax*a*100.E5*hubble*Hsub;
        r2kineticfact = size/hubble/amax*a*a*100.E5*hubble*Hsub;
        sim2real = size/hubble/rng/amax*a*1.E6*pc;
        real2sim = 1./sim2real;
        pntmass = 3./8./pi/G*sqr(100.E5*hubble)*1.E6*pc;
        pntmass = pntmass*pow3(size/hubble)/pow3(rng/nspace)*omep;
        potentfact = G*pntmass/sim2real;
        massscale = pntmass/onesolarmass*hubble;
        rscale = size/nx;
        printf("masscale=%g rscale = %g\n",massscale,rscale);
        /* since we scale position to the simulation coordiate */

        index = 0;
        haloindex = 0;
           while(fread(&mpeak,sizeof(int),1,fp) == 1){
            double mass;
            ncount ++;
            peakindex = 0;
            for(i=0;i<mpeak;i++){
                fread(&nend,sizeof(int),1,fp);
                anend[i] = nend;
                r = (PtlPos *) Malloc(nend*sizeof(PtlPos),PPTR(r));
                vr = (PtlVel *) Malloc(nend*sizeof(PtlVel),PPTR(vr));
                indx = (indxtype *) Malloc(nend*sizeof(indxtype),PPTR(indx));
                if(fread(r,sizeof(PtlPos),nend,fp)!= nend){
                    fprintf(stderr,"error reading r nend %6d in mpeak %5d\n",
                            nend,mpeak);
                    exit(-999);
                }
                if(fread(vr,sizeof(PtlVel),nend,fp)!= nend){
                    fprintf(stderr,"error reading vr nend %6d in mpeak %5d\n",
                            nend,mpeak);
                    exit(-999);
                }
                if(fread(indx,sizeof(indxtype),nend,fp) != nend){
                    fprintf(stderr,"error reading indx nend %6d in mpeak %5d\n",
                            nend,mpeak);
                    exit(-999);
                }
                if(nend <= 5) {
                    Free(indx);Free(vr);Free(r);
                    continue;
                }
                for(j=0;j<nend;j++){
                    /* WARNING */
                    /* rescale to the original simulation coordiantes */
                    r[j].x = r[j].x * rng / size;
                    r[j].y = r[j].y * rng / size;
                    r[j].z = r[j].z * rng / size;
                    vr[j].vx = vr[j].vx/r2kineticfact*1.E5;
                    vr[j].vy = vr[j].vy/r2kineticfact*1.E5;
                    vr[j].vz = vr[j].vz/r2kineticfact*1.E5;
				}
				{
					int Numnear;
					Numnear = 27;
					DenPeakCenter((vector3d*)r,nend,&cx,&cy,&cz,Numnear);
				}
                {
                    int nxwidth,nywidth,nzwidth;
                    int nfx,nfy,nfz;
                    int nbuff;
                    int iiii;
                    float h = 1.;
                    float fmax ;
                    if(1){
                        printf("halo %d has %d particles :  %d\n",
                                i,nend,mpeak);
                        {
                            /* here insert a routine for determining the
                             * characteristics of halos */
                            halo_properties(nend,r,vr);
                        }
                        {
                            float *px,*py,*pz;
                            px = (float *) Malloc(sizeof(float)*10000,PPTR(px));
                            py = (float *) Malloc(sizeof(float)*10000,PPTR(py));
                            pz = (float *) Malloc(sizeof(float)*10000,PPTR(pz));
                            r1[0] = v[0][0]; r1[1] = v[0][1]; r1[2] = v[0][2];
                            r2[0] = v[1][0]; r2[1] = v[1][1]; r2[2] = v[1][2];
                            r3[0] = v[2][0]; r3[1] = v[2][1]; r3[2] = v[2][2];
                            radius = find_r(r,nend,q,s,r1,r2,r3)*size/rng;
                            Free(pz);
                            Free(py);
                            Free(px);
                        }
                        {
                            cvx = cvy = cvz = 0;
                            for(j=0;j<nend;j++){
                                cvx += vr[j].vx;
                                cvy += vr[j].vy;
                                cvz += vr[j].vz;
                            }
                            cvx = cvx/(double)nend;
                            cvy = cvy/(double)nend;
                            cvz = cvz/(double)nend;
                            sigv = 0.;
                            for(j=0;j<nend;j++){
                                sigv += pow2(vr[j].vx-cvx)+pow2(vr[j].vy-cvy)+pow2(vr[j].vz-cvz);
                            }
                            sigv = sigv/(double) nend;
                            sigv = sqrt(sigv);
                            sigv = sigv*r2kineticfact/1.E5;
                        }
                        /*
                        fprintf(wp,"halo %d has %d particles in %d peaks:  %g\n",
                                i,nend,mpeak,fmax);
                                */
                        mass = nend*pntmass/onesolarmass;
                        /* pvr : peculiar velocity in km/sec */
                        pvr = sqrt(cvx*cvx+cvy*cvy+cvz*cvz)*r2kineticfact/1.E5;
                        if(isnan(para_spin) != 0) {
                            para_spin = -9999.;
                        }
                        {
                            double ampang;
                            ampang = pow2(ang_mom[0])+pow2(ang_mom[1])+pow2(ang_mom[2]);
                            ampang = sqrt(ampang);
                            ang_mom[0] = ang_mom[0]/ampang;
                            ang_mom[1] = ang_mom[1]/ampang;
                            ang_mom[2] = ang_mom[2]/ampang;
                        }
                        fprintf(wp,"%g %f %f %f %d %d %d %g %g %g %g %g %g ",
                                mass,cx*size/rng,cy*size/rng,cz*size/rng,index,
                                mpeak,i,para_spin,sigv,rotation,radius,q,s);
                        fprintf(wp,"%g %g %g %g %g %g ",v[0][0],v[0][1],v[0][2],
                                v[1][0],v[1][1],v[1][2]);
                        fprintf(wp,"%g %g %g  ",v[2][0],v[2][1],v[2][2]);
                        fprintf(wp,"%g %g %g ",ang_mom[0],ang_mom[1],ang_mom[2]);
                        cvx = cvx*r2kineticfact/1.E5;
                        cvy = cvy*r2kineticfact/1.E5;
                        cvz = cvz*r2kineticfact/1.E5;
                        fprintf(wp,"%g %g %g\n",cvx,cvy,cvz);
                        fflush(wp);
                        fflush(stdout);
                        /* 첫 time step 인 경우 */
                        nnowhaloindex++;
                        haloindex++;
                    }
                }
                peakindex ++;
                Free(indx);
                Free(vr);
                Free(r);
            }
			{
				/* This part will perform the FoF searching */
			}
            index ++;
        }
        if(feof(fp)){
            printf("exit normally\n");
        }
        else if(ferror(fp)){
            printf("exit abnormally\n");
        }
        nnowhaloindex = 0;
        fclose(fp);
        fclose(wp);
        /*
        fclose(mergerwp);
        */
    Free(den);
}
float piinhin3;
float h;
void findden(int np, PtlPos *r, int nxwidth,int nywidth,int nzwidth){
    int i,j,k;
    int nbuff;
    float tmp2h;
    
    piinhin3 = 1./3.141592653;
    h = 1.;
    nbuff = 2*h+1;
    tmp2h = 2.*h;
        for(i=0;i<np;i++){
                float xi,yi,zi,xf,yf,zf;
                int mxi,myi,mzi;
                int mxf,myf,mzf;
                float wr,wresult,gx,gy,gz;
                int ii,ij,ik;
                int iii,iij,iik;
                float xp,yp,zp;
                xp = r[i].x;
                yp = r[i].y;
                zp = r[i].z;
                xi = xp - tmp2h; yi = yp - tmp2h; zi = zp - tmp2h;
                xf = xp + tmp2h; yf = yp + tmp2h; zf = zp + tmp2h;
                mxi = rint(xi) + 1; myi = rint(yi) + 1; mzi = rint(zi) + 1;
                mxf = rint(xf); myf = rint(yf); mzf = rint(zf);
                for(ik = mzi;ik<=mzf;ik++)
                for(ij = myi;ij<=myf;ij++)
                for(ii = mxi;ii<=mxf;ii++) {
                        float tmpx,tmpy,tmpz;
                        gx = ii - 0.5; gy = ij - 0.5; gz = ik - 0.5;
                        tmpx = gx-xp;tmpx = tmpx*tmpx;
                        tmpy = gy-yp;tmpy = tmpy*tmpy;
                        tmpz = gz-zp;tmpz = tmpz*tmpz;
                        wr = sqrt(tmpx + tmpy + tmpz);
                        /*
                        wr = w_(&wr);
                        */
                        wresult = wc(wr);
                        den[ii+ij*nxwidth+ik*nxwidth*nywidth] += wresult;
                }
        }
}
float wc(float r){
        float w,x;
        x = r/h;
        if(x < 1.)
                w = piinhin3*(1.-1.5*x*x+0.75*x*x*x);
        else if(x < 2.)
                w = piinhin3*0.25*(2.-x)*(2.-x)*(2.-x);
        else w = 0;
        return w;
}
