     !Simple test program to print out sigma_8 as a function of the CDM density
      module MyModule
         use CAMB
         type (CAMBdata) OutData !type (CAMBdata) OutData
         type(CAMBparams) P !defined in ModelParams in modules.f90
      end module MyModule
      subroutine InitPower(omep,omepb,omeplam,wlam0,wlam1,Cs2,hubble,npower,redshift, pamp,&
                  boxsize,rng)
        use CAMB
        use MyModule
        use InitialPower
        use LambdaGeneral
        implicit none
        integer i,j,k
        real PI2,rng,maxkhfactor,wlam0, wlam1,Cs2
        parameter(PI2=3.1415926535d0*2.d0,maxkhfactor=5)
         type (CAMBdata) InData 
    
        type(MatterTransferData) MTrans
        real omep,omepb,omeplam,hubble,npower,maxkh,redshift,boxsize
        integer itf,in
        integer error
        real minkh,dlnkh,pamp,khmax
        integer npoints,nmatterpk
        character(LEN=80) fmt
        logical OK
        double precision cubicbox
           cubicbox= boxsize
           cubicbox = cubicbox**3

           call CAMB_SetDefParams(P)

           w_lam = wlam0

           P%WantTransfer= .true.
           P%WantCls = .false.

           maxkh = max(PI2/boxsize*rng*maxkhfactor,1.5)
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
             
!          This cubic box multiplication is to normalize the amplitude in the unit of 2pi*u not k.
           pamp = (MT%sigma_8(1,1))**2*cubicbox
      end subroutine InitPower

! In this old CAMB version, we didn't implement power spectrum amplitude
!   (As) scaling in this routine.
      subroutine rInitPower(omep,omepb,omeplam,wlam0,wlam1,Cs2, &
                hubble,npower,redshift, As,pamp,&
                  boxsize,rng,readflag,infile,myid)
        use CAMB
        use MyModule
        use InitialPower
        use LambdaGeneral
        implicit none
        integer i,j,k
        integer, intent(in) :: readflag,myid
        character(LEN=80), intent(in) :: infile
        real PI2,rng,maxkhfactor
        parameter(PI2=3.1415926535d0*2.d0,maxkhfactor=1.5)
         type (CAMBdata) InData 
    
        type(MatterTransferData) MTrans
        real omep,omepb,omeplam,hubble,npower,maxkh,redshift,boxsize
        integer itf,in
        integer error
        real minkh,dlnkh,pamp,khmax,As
        integer npoints,nmatterpk
        character(LEN=80) fmt
        logical OK
        double precision cubicbox
        real amax,anow,ai,bini,bnow,growth
        real wlam0,wlam1, Cs2
        external growth
        cubicbox = boxsize**3


!       print *, omep,omepb,omeplam,wlam0,wlam1,Cs2
!       print *, hubble,npower,redshift, As,pamp
!       print *, boxsize, rng,readflag,myid
!       print *, 'input file', infile
        w_lam = wlam0
                    
	        open(1,file=infile,form='unformatted')
!           read(1) pamp
           read(1) OutData%Params
           ! Do not delete or modify this line because camb subroutines need it!!!
           CP = OutData%Params 
           read(1) OutData%MTrans%num_q_trans
           call Transfer_Allocate(OutData%MTrans)
           read(1) (OutData%MTrans%q_trans(i),i=1,OutData%MTrans%num_q_trans)
           read(1) (((OutData%MTrans%TransferData(i,j,k),i=1,Transfer_max),j=1,OutData%MTrans%num_q_trans), &
              k=1,CP%Transfer%num_redshifts)
           read(1) ((OutData%MTrans%sigma_8(i,j),i=1,CP%Transfer%num_redshifts),j=1,CP%InitPower%nn)

           close(1)
           MT = OutData%MTrans
           pamp = (MT%sigma_8(1,1))**2*cubicbox
           ! These two lines are for initialization of the power spectrum call
           P = CP ! Do not delete or modify this line because camb subroutines need it!!!
           call InitializePowers(CP%InitPower,CP%curv)

           if(myid.eq.0) then
           write(*,*)
           write(*,*)
           print *, '#####################################################'
           print *, '####      POWER SPECTRUM PARAMETER READED  ##########'
           print *, '#####################################################'
           print *, 'omep                  = ',CP%omegac+CP%omegab
           print *, 'omepb                 = ',CP%omegab
           print *, 'omeplam               = ',CP%omegav
           print *, 'npower                = ',CP%InitPower%an(1)
           print *, 'Hubble (km/sec.)      = ',CP%H0
           if(CP%OutputNormalization) then
           print *, 'Normalization applied   '
           else
           print *, 'Normalization unapplied '
           endif
           print *, 'simulation amplitude  = ',pamp
           print *, 'power Amplitude       = ',CP%InitPower%ScalarPowerAmp(1)
           print *, 'number of redshift    = ',CP%InitPower%nn
           print *, 'Redshift              = ',CP%Transfer%redshifts(1)
           print *, 'Highprecision         = ',CP%Transfer%high_precision
           print *, 'kmax                  = ',CP%Transfer%kmax
           print *, 'k per logint          = ',CP%Transfer%k_per_logint
           print *, '#####################################################'
           print *, '####      POWER SPECTRUM PARAMETER READED  ##########'
           print *, '#####################################################'
           if(abs(redshift -CP%Transfer%redshifts(1)).ge. 1.e-3) then
           print *, '#####################################################'
           print *, '####      POWER SPECTRUM READSHIFT CHANGED ##########'
           print *, '# From ##############################################'
           print *, 'Redshift              = ',CP%Transfer%redshifts(1)
           print *, '# To   ##############################################'
           print *, 'Redshift              = ',redshift
           print *, '#####################################################'
           endif
           write(*,*)
           write(*,*)
           endif

      end subroutine rInitPower

      subroutine GetPower(minkh,dlnkh,matterpk,nmatterpk,Ts)
        use CAMB
        use MyModule
        use transfer
!        use InitialPower
        implicit none
        integer i
    
        type(MatterTransferData) MTrans
        integer itf,in
        integer error
        real minkh,dlnkh,logminkh
        integer npoints,nmatterpk
!       real, dimension(:,:), allocatable :: outpower
        real matterpk(nmatterpk)
        real Ts(nmatterpk)
        character(LEN=80) fmt
        logical OK
        real(dl) logmink,k,h


        minkh = 5e-5
        dlnkh = 0.02
        Mtrans = OutData%MTrans


        itf = 1

!        print *, Mtrans%num_q_trans,MTrans%TransferData(Transfer_kh,Mtrans%num_q_trans,itf)

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
         h = CP%H0/100
         logmink = log(minkh)
         do i = 1, npoints
            k = exp(logmink+dlnkh*(i-1))*h
            Ts(i) = matterpk(i)/ScalarPower(k,in)/(k*pi*twopi*h**3)*ScalarPower(k,1)
         enddo

        end subroutine GetPower

