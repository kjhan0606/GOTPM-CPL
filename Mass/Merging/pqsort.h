void my_MPI_Sendrecv(void *, long long , MPI_Datatype , int , int , void *, long long , MPI_Datatype , int , int , MPI_Comm , MPI_Status *);
void my_MPI_Send(void *, long long , MPI_Datatype , int , int , MPI_Comm );
void my_MPI_Recv(void *, long long , MPI_Datatype , int , int , MPI_Comm , MPI_Status *);

#define min(a,b) ((a)<(b)? (a):(b))
#define max(a,b) ((a)>(b)? (a):(b))


#define ISUNSIGNED(a) (((typeof(a))-1) >= 0)

#define PQSORT(base,TYPE,subtype,mem,nmem,Comp,min,max,Comm) \
	long pqsort##TYPE##mem(TYPE *base,long nmem,\
		int (*Comp)(const void *, const void *),\
		subtype min,subtype max,MPI_Comm Comm){\
	MPI_Status status;\
	int src,dest;\
	long nrecv,nsend;\
	long commsize,i,j,k;\
	int myid,nid;\
	subtype lmin,lmax;\
	MPI_Comm_rank(Comm,&myid);\
	MPI_Comm_size(Comm,&nid);\
	dest = (myid+nid+1)%nid;\
	src = (myid+nid-1)%nid;\
	subtype localmin = myid*(max-min)/nid;\
	subtype localmax = (myid==nid-1? max: (myid+1)*(max-min)/nid);\
	subtype ldist,rdist;\
	TYPE *left,*right;\
	commsize =10;\
	while(commsize){\
		left = (TYPE*)base;\
		right = (TYPE*)base + nmem;\
		for(;left<right;){\
			if(left->mem >= localmax || left->mem < localmin){\
				rdist = fabs(localmax-left->mem);\
				rdist = min(rdist,max-rdist);\
				ldist = fabs(localmin-left->mem);\
				ldist = min(ldist,max-ldist);\
				if(rdist <= ldist){\
					right --;\
					TYPE tmp = *left;\
					*left = *right;\
					*right = tmp;\
				}\
				else left ++;\
			}\
			else left ++;\
		}\
		TYPE *bp2send = right;\
		long  nsend = nmem-(bp2send-base);\
		MPI_Sendrecv(&nsend,1,MPI_LONG,dest,0,&nrecv,1,MPI_LONG,src,0,Comm,&status);\
		MPI_Reduce(&nsend,&commsize,1,MPI_LONG,MPI_SUM,0,Comm);\
		MPI_Bcast(&commsize,1,MPI_LONG,0,Comm);\
		if(myid==0) {\
			printf("Total commsize = %ld",commsize);\
		}\
		long net = nrecv-nsend;\
		TYPE *tmp = (TYPE*)malloc(sizeof(TYPE)*nsend);\
		i = 0;\
		for(right = bp2send;right<base+nmem;right++) tmp[i++]=*right;\
		base = (TYPE*)realloc(base,sizeof(TYPE)*(nmem+net));\
		left = base+nmem+net-nrecv;\
		MPI_Sendrecv(tmp,nsend*sizeof(TYPE),MPI_BYTE,dest,0,\
				left,nrecv*sizeof(TYPE),MPI_BYTE,src,0,Comm,&status);\
		free(tmp);\
		nmem = nmem+net;\
	}\
	commsize = 100;\
	dest = (myid+nid-1)%nid;\
	src = (myid+nid+1)%nid;\
	while(commsize){\
		left=(TYPE*)base;\
		right = (TYPE*)base + nmem;\
		for(;left<right;){\
			if(left->mem>=localmax || left->mem < localmin){\
				right--;\
				TYPE tmp = *left;\
				*left = *right;\
				*right = tmp;\
			}\
			else left ++;\
		}\
		TYPE *bp2send = right;\
		long nsend = nmem - (bp2send-base);\
		MPI_Sendrecv(&nsend,1,MPI_LONG,dest,0,&nrecv,1,MPI_LONG,src,0,Comm,&status);\
		MPI_Reduce(&nsend,&commsize,1,MPI_LONG,MPI_SUM,0,Comm);\
		MPI_Bcast(&commsize,1,MPI_LONG,0,Comm);\
		if(myid==0) {\
			printf("Total commsize = %ld",commsize);\
		}\
		long net = nrecv-nsend;\
		TYPE *tmp = (TYPE*)malloc(sizeof(TYPE)*nsend);\
		i = 0;\
		for(right = bp2send;right<base+nmem;right++) tmp[i++]=*right;\
		base = (TYPE*)realloc(base,sizeof(TYPE)*(nmem+net));\
		left = base+nmem+net-nrecv;\
		MPI_Sendrecv(tmp,nsend*sizeof(TYPE),MPI_BYTE,dest,0,\
				left,nrecv*sizeof(TYPE),MPI_BYTE,src,0,Comm,&status);\
		free(tmp);\
		nmem = nmem+net;\
	}\
	qsort(base,nmem,sizeof(TYPE),Comp);\
	return nmem;\
} 

#ifdef __MAIN_PQSORT__

#define GenerateParallelComp(TYPE,subtype,mem)  \
	int Psort_##TYPE##_##subtype##_##mem(const void *a, const void *b){\
		TYPE *aa = (TYPE *)a;\
		TYPE *bb = (TYPE *)b;\
		if(aa->mem < bb->mem) return -1;\
		else if(aa->mem > bb->mem) return 1;\
		else return 0;\
	}\
	int Psort_##TYPE##_##subtype##_##mem##_L(const void *a, const void *lb, const void *lc){\
		if(((TYPE*)lb)->mem == ((TYPE*)lc)->mem) return 1;\
		else if(Psort_##TYPE##_##subtype##_##mem(a,lc) >=0 || Psort_##TYPE##_##subtype##_##mem(a,lb) ==-1) return 1;\
		else return 0;\
	}\
	int Psort_##TYPE##_##subtype##_##mem##_R(const void *a, const void *lb, const void *lc, int nid){\
		if(((TYPE*)lb)->mem == ((TYPE*)lc)->mem) return 1;\
		else if(Psort_##TYPE##_##subtype##_##mem(a,lc) >=0 || Psort_##TYPE##_##subtype##_##mem(a,lb) ==-1) {\
			TYPE *aa = (TYPE*) a;\
			TYPE *lbb = (TYPE*) lb;\
			TYPE *lcc = (TYPE*) lc;\
			int rdist;\
			rdist = (aa->mem - lbb->mem)/(lcc->mem-lbb->mem);\
			rdist = (rdist + nid)%nid;\
			if(rdist <= nid/2) return 1;\
			else return 0;\
		}\
		else return 0;\
	} 
#endif

void *pqsort(void *, size_t *, size_t, int (*)(const void *, const void *, const void *, int ), 
		int (*)(const void *, const void *, const void *), 
		int (*)(const void *, const void *), void *, void *, MPI_Comm);

#define MyPqsort(base,nmem,n_size,TYPE,subtype,mem,Comm) do{\
	int Psort_##TYPE##_##subtype##_##mem(const void *, const void *);\
	int Psort_##TYPE##_##subtype##_##mem##_L(const void *, const void *, const void *);\
	int Psort_##TYPE##_##subtype##_##mem##_R(const void *, const void *, const void *, int );\
	int id,mid;\
	MPI_Datatype mpi_datatype;\
	MPI_Comm_rank(Comm,&id);\
	MPI_Comm_size(Comm,&mid);\
	TYPE lmin,lmax,gmin,gmax;\
	subtype almin,almax,agmin,agmax,tagmin,tagmax;\
	if(sizeof(subtype) == sizeof(char)){\
		if(ISUNSIGNED(subtype)) {\
			mpi_datatype = MPI_UNSIGNED_CHAR;\
			agmax = 0;\
			agmin = UCHAR_MAX;\
		}\
		else {\
			mpi_datatype = MPI_CHAR;\
			agmax = SCHAR_MIN;\
			agmin = SCHAR_MAX;\
		}\
	}\
	else if(sizeof(subtype) == sizeof(int)){\
		if(ISUNSIGNED(subtype)) {\
			mpi_datatype = MPI_UNSIGNED;\
			agmax = 0;\
			agmin = UINT_MAX;\
		}\
		else {\
			mpi_datatype = MPI_INT;\
			agmax = INT_MIN;\
			agmin = INT_MAX;\
		}\
	}\
	else if(sizeof(subtype) == sizeof(long)){\
		if(ISUNSIGNED(subtype)){\
			mpi_datatype = MPI_UNSIGNED_LONG;\
			agmax = 0;\
			agmin = (subtype)ULONG_MAX;\
		}\
		else {\
			mpi_datatype = MPI_LONG;\
			agmax = (subtype)LONG_MIN;\
			agmin = (subtype)LONG_MAX;\
		}\
	}\
	else if(sizeof(subtype) == sizeof(long long)){\
			mpi_datatype = MPI_LONG;\
			agmax = (subtype)LONG_MIN;\
			agmin = (subtype)LONG_MAX;\
	}\
	long long i;\
	for(i=0;i<nmem;i++){\
		if( base[i].mem < agmin ) agmin = base[i].mem;\
		if( base[i].mem > agmax ) agmax = base[i].mem;\
	}\
	if(strcmp(#subtype,"float")==0 || strcmp(#subtype,"double")==0) agmax = agmax+0.0000001*(agmax-agmin);\
	else agmax ++;\
	MPI_Reduce(&agmin,&tagmin,1,mpi_datatype,MPI_MIN,0,Comm);\
	MPI_Reduce(&agmax,&tagmax,1,mpi_datatype,MPI_MAX,0,Comm);\
	MPI_Bcast(&tagmin,1,mpi_datatype,0,Comm);\
	MPI_Bcast(&tagmax,1,mpi_datatype,0,Comm);\
	lmin.mem = id*tagmax/mid;\
	lmax.mem = (id+1)*tagmax/mid;\
	printf("-P%d has %s local min/max = %ld %ld for %ld members\n",id,#base"."#mem,lmin.mem,lmax.mem,nmem);fflush(stdout);\
	if(id == mid-1) lmax.mem = tagmax;\
	size_t mmem = nmem;\
	base = (TYPE*)pqsort(base,&mmem,n_size,Psort_##TYPE##_##subtype##_##mem##_R, Psort_##TYPE##_##subtype##_##mem##_L, \
			Psort_##TYPE##_##subtype##_##mem,&lmin,&lmax,Comm);\
	nmem = mmem;\
	printf("+P%d has %s local min/max = %ld %ld for %ld members\n",id,#base"."#mem,lmin.mem,lmax.mem,nmem);fflush(stdout);\
} while(0)

#define MyConstRangePqsort(base,nmem,n_size,TYPE,subtype,mem,Gmin,Gmax,Comm) do{\
	int Psort_##TYPE##_##subtype##_##mem(const void *, const void *);\
	int Psort_##TYPE##_##subtype##_##mem##_L(const void *, const void *, const void *);\
	int Psort_##TYPE##_##subtype##_##mem##_R(const void *, const void *, const void *, int);\
	int _id,_mid;\
	MPI_Datatype mpi_datatype;\
	MPI_Comm_rank(Comm,&_id);\
	MPI_Comm_size(Comm,&_mid);\
	TYPE lmin,lmax,gmin,gmax;\
	lmin.mem = _id*( (Gmax-Gmin-1)/_mid+1) + Gmin;\
	lmax.mem = (_id+1)*( (Gmax-Gmin-1)/_mid+1) + Gmin;\
	if(lmin.mem > Gmax) lmin.mem = Gmax;\
	if(lmax.mem > Gmax) lmax.mem = Gmax;\
	printf("-P%d has %s local min/max = %ld %ld for %ld members\n",_id,#base"."#mem,lmin.mem,lmax.mem,nmem);fflush(stdout);\
	size_t mmem = nmem;\
	base = (TYPE*)pqsort(base,&mmem,n_size,Psort_##TYPE##_##subtype##_##mem##_R, Psort_##TYPE##_##subtype##_##mem##_L, \
			Psort_##TYPE##_##subtype##_##mem,&lmin,&lmax,Comm);\
	nmem = mmem;\
	printf("+P%d has %s local min/max = %ld %ld for %ld members\n",_id,#base"."#mem,lmin.mem,lmax.mem,nmem);fflush(stdout);\
} while(0)
