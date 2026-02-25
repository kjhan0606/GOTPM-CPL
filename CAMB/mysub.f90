     !Simple test program to print out sigma_8 as a function of the CDM density
      subroutine InitPower(omep,omepb,omeplam,hubble,npower,redshift,matterpk,nmatterpk)
        use CAMB
        implicit none
        integer i
    
        type(CAMBparams) P !defined in ModelParams in modules.f90
        type (CAMBdata) OutData
        type(MatterTransferData) MTrans
        real omep,omepb,omeplam,hubble,npower,maxkh,redshift
        integer itf,in
        integer error
        real minkh,dlnkh
        integer npoints,nmatterpk
        real, dimension(:,:), allocatable :: outpower
        real, dimension(:) :: matterpk
        character(LEN=80) fmt
        logical OK

        call CAMB_SetDefParams(P)

        P%WantTransfer= .true.
        P%WantCls = .false.

!        omep = 0.26
!        omepb = 0.0441
!        omeplam = 0.74
!        hubble = 0.72
!        npower = 0.96
!        redshift = 0.
        maxkh = 1000.

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

        write (fmt,*) CP%InitPower%nn+1
        fmt = '('//trim(adjustl(fmt))//'E15.5)'


             
        write (*,*) 'Omc = ',real(P%Omegac),'OmLam=',real(P%Omegav) &
           , 'sigma_8 = ', real(MT%sigma_8(1,1))

        minkh = 1e-4
        dlnkh = 0.02
        Mtrans = OutData%MTrans

!         npoints = MTrans%num_q_trans
!         do in=1, npoints
!            write (*, '(7E14.6)') MTrans%TransferData(Transfer_kh:Transfer_max,in,1)
!         enddo

        itf = 1
         npoints = log(MTrans%TransferData(Transfer_kh,Mtrans%num_q_trans,itf)/minkh)/dlnkh+1
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
         allocate(outpower(npoints,CP%InitPower%nn))
         do in = 1, CP%InitPower%nn
            call Transfer_GetMatterPower(Mtrans,outpower(1,in),itf,in,minkh,dlnkh,npoints)
         enddo
         do i =1, npoints
!            write (*, fmt) minkh*exp((i-1)*dlnkh),outpower(i,1:CP%InitPower%nn)
             matterpk(i) = outpower(i,1)
         enddo
         deallocate(outpower)

        end subroutine InitPower

