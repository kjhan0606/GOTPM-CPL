     !Simple test program to print out sigma_8 as a function of the CDM density
      module MyModule2
         use CAMB
         type (CAMBdata) OutData,InData
      end module MyModule2

      program InitPower
        use CAMB
        use MyModule2
        use LambdaGeneral
        implicit none
        integer i,j,k
        real PI2,rng,maxkhfactor,wlam0
        parameter(PI2=3.1415926535d0*2.d0,maxkhfactor=5)
    
        type(CAMBparams) P !defined in ModelParams in modules.f90
        !type (CAMBdata) OutData
        type(MatterTransferData) MTrans
        real omep,omepb,omeplam,hubble,npower,maxkh,redshift,boxsize
        integer itf,in,narg
        integer error
        real minkh,dlnkh,pamp,bias8,khmax,opamp
        real amax,anow,da
        integer npoints,nmatterpk,ireadflag
        real matterpk(8192), Ts(8192), rk
        character(LEN=80) fmt,inputfile,outfile,infile
        logical OK
        double precision cubicbox
        real w0,wa,cs2

        w0 = -1
        wa = 0
        cs2 = 1



        narg = iargc()
		if(narg .ne. 1) then
		   print *,'error occurred. Please input the parameter file'
		   stop
		endif
        call getarg(1,inputfile)
!       open(1,file=inputfile,form='formatted')
!       read(1,*) boxsize,hubble
!       read(1,*) npower,omep,omepb,omeplam,bias8
!       read(1,*) rng
!       read(1,*)
!       read(1,*) amax,da,anow
!       read(1,*) 
!       read(1,*) 
!       read(1,*) 
!       read(1,*) 
!       read(1,*) i,outfile
!       close(1)
        infile = trim(inputfile)//char(0)
        call getsimparm(infile, boxsize,hubble,npower,omep,omepb, &
                 omeplam,wlam0, bias8,rng,amax,da,anow,ireadflag,outfile)
        w_lam = wlam0
        print *, 'W(DE) = ',w_lam
!        amax =amax + 1
        redshift = amax/anow - 1
       
        cubicbox= boxsize
        cubicbox = cubicbox**3

        call CAMB_SetDefParams(P)

        P%WantTransfer= .true.
        P%WantCls = .false.

        maxkh = max(PI2/boxsize*rng*maxkhfactor,500.)
        maxkh = min(maxkh,50000.)


        P%omegab  = omepb
        P%omegac  = omep - omepb
        P%omegav  = omeplam
        P%omegan  = 0.
        P%H0      = hubble*100.
        P%OutputNormalization = 1
       
        P%InitPower%ScalarPowerAmp = 2.076e-9
        P%InitPower%nn     = 1 !number of initial power spectra
        P%InitPower%an(1)  = npower !scalar spectral index
        P%InitPower%ant(1) = 0 !Not used here
        P%InitPower%rat(1) = 1 !ditto

        !these settings seem good enough for sigma8 to a percent or so
        P%Transfer%high_precision=.true.
        P%Transfer%kmax=maxkh
        P%Transfer%k_per_logint=10
        P%Transfer%num_redshifts=1
        P%Transfer%redshifts(1)=redshift


        OK= CAMB_ValidateParams(P)
        if(OK .eq. .false.) then
           print *,'###################################'
           print *,'###################################'
           print *,'Wrong set of parameters in the CAMB'
           print *,'###################################'
           print *,'###################################'
           print *,'Ok = ',OK
           stop
        endif


        !call CAMB_GetResults(P) 
        call CAMB_GetTransfers(P,OutData)
             
!       This cubic box multiplication is to normalize the amplitude in the unit of 2pi*u not k.
        pamp = (MT%sigma_8(1,1))**2*cubicbox

        print *,'outfile = ', outfile, ireadflag

        if(ireadflag .eq. 2) then
        open(1,file=outfile)
        write(1,'("Om0=", g14.7," Ob0=", g14.7, &
                 " Ode0=", g10.3," H0=", g12.5," red=", g10.3 )') omep,omepb,omeplam, hubble,redshift
        write(1,'("w0=", g14.7," wa=", g14.7," c_s^2(DE)=", g10.3)') w0,wa,cs2
        write(1,'(g14.7," =sigma8 at R=8Mpc/h"  )') MT%sigma_8(1,1)
        write(1,*) '#   k/h        Ptotm(k)   Pcmb(k)   Pbaryon(k)'
        nmatterpk = 8192
         call GetPower2(minkh,dlnkh,matterpk,nmatterpk)
         do i = 1, nmatterpk
            rk = minkh*exp((i-1)*dlnkh)
            write(1,*) rk, matterpk(i)
         enddo

        close(1)
        else
        open(1,file=outfile,form='unformatted')
        write(1) OutData%Params
        write(1) OutData%MTrans%num_q_trans
        write(1) (OutData%MTrans%q_trans(i),i=1,OutData%MTrans%num_q_trans)
        write(1) (((OutData%MTrans%TransferData(i,j,k),&
                i=1,Transfer_max), &
                j=1,OutData%MTrans%num_q_trans), &
                k=1,CP%Transfer%num_redshifts)
        write(1) ((OutData%MTrans%sigma_8(i,j), &
                i=1,CP%Transfer%num_redshifts), &
                j=1,CP%InitPower%nn)
        close(1)
        endif
       
      end program InitPower

      subroutine GetPower2(minkh,dlnkh,matterpk,nmatterpk)
        use CAMB
        use MyModule2
        use LambdaGeneral
        implicit none
        integer i
    
        type(CAMBparams) P !defined in ModelParams in modules.f90
        !type (CAMBdata) OutData
        type(MatterTransferData) MTrans
        integer itf,in
        integer error
        real minkh,dlnkh
        integer npoints,nmatterpk
        real, dimension(:,:), allocatable :: outpower
        real matterpk(nmatterpk)
        character(LEN=80) fmt
        logical OK


        minkh = 5e-5
        dlnkh = 0.02
        Mtrans = OutData%MTrans

        itf = 1
         npoints = alog(MTrans%TransferData(Transfer_kh,Mtrans%num_q_trans,itf)/minkh)/dlnkh+1
         if(npoints .gt. nmatterpk) then
            print *,'#################################################'
            print *,'#################################################'
            print *,'Small number of power spectrum array ',nmatterpk
            print *,'The npoints in CAMB is ', npoints
            print *,'#################################################'
            print *,'#################################################'
            stop
         else
            nmatterpk = npoints
         endif
         in = 1
         call Transfer_GetMatterPower(Mtrans,matterpk(1),itf,in,minkh,dlnkh,npoints)

        end subroutine GetPower2

