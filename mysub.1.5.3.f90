     !Simple test program to print out sigma_8 as a function of the CDM density
      module MyModule
         use CAMB
         type (CAMBdata) OutData !type (CAMBdata) OutData
         type(CAMBparams) P !defined in ModelParams in modules.f90
         integer n_source_points
      end module MyModule
!--------------
      module MyPkModule
        real, parameter:: m_minkh = 1e-4 ! minimum k/h for P(k)
        real, parameter:: m_dlnkh = 0.02 ! log step size of k/h for P(k)
        real, parameter:: maxnp = 8192
        integer nkh
        real rk(maxnp), pk(maxnp),h100
      end module MyPkModule
!--------------
      subroutine InitPower(omep,omepb,omeplam,DarkEnergyModel,w0,wa,cs2,&
                hubble,npower,redshift, As,pamp,boxsize,rng)
        use CAMB
        use DarkEnergyInterface
        use DarkEnergyFluid
        use DarkEnergyPPF
        use Quintessence
        use classes
        use MyModule
        use MyPkModule
        implicit none
        integer i,j,k
        real PI2,rng,maxkhfactor,wlam
        real w0,wa,cs2
        parameter(PI2=3.1415926535d0*2.d0,maxkhfactor=5)
    
        type(MatterTransferData) MT
        real omep,omepb,omeplam,hubble,npower,maxkh,redshift,boxsize
        integer itf,in
        integer error
        real minkh,dlnkh,pamp,As,khmax
        integer npoints,nmatterpk
        character(LEN=80) fmt
        character(LEN=30) DarkEnergyModel
        logical OK
        double precision cubicbox
        cubicbox= boxsize
        cubicbox = cubicbox**3

        call CAMB_SetDefParams(P)
        h100 = hubble

        if(allocated(P%DarkEnergy)) deallocate(P%DarkEnergy)
        if(DarkEnergyModel(1:5) == 'FLUID') then
           allocate(TDarkEnergyFluid::P%DarkEnergy)
           if(w0==-1 .and. wa==0) then
              P%DarkEnergy%is_cosmological_constant = .true.
              P%DarkEnergy%num_perturb_equations = 0
           else 
              P%DarkEnergy%is_cosmological_constant = .false.
              P%DarkEnergy%num_perturb_equations = 2
           endif
        else if(DarkEnergyModel(1:3) == 'PPF') then
           allocate(TDarkEnergyPPF::P%DarkEnergy)
           if(w0==-1 .and. wa==0) then
              P%DarkEnergy%is_cosmological_constant = .true.
              P%DarkEnergy%num_perturb_equations = 0
           else 
              P%DarkEnergy%is_cosmological_constant = .false.
              P%DarkEnergy%num_perturb_equations = 1
           endif
        else if(DarkEnergyModel(1:19) == 'AXIONEFFECTIVEFLUID') then
           allocate(TAxionEffectiveFluid::P%DarkEnergy)
           P%DarkEnergy%num_perturb_equations = 2
        else if (DarkEnergyModel(1:17) == 'EARLYQUINTESSENCE') then
           allocate(TQuintessence::P%DarkEnergy)
           P%DarkEnergy%num_perturb_equations = 2
        endif


        select type(TDE=>P%DarkEnergy)
        class is (TDarkEnergyPPF)
            TDE%w_lam = w0
            TDE%wa  = wa
            TDE%cs2_lam  = cs2
        class is (TDarkEnergyFLUID)
            TDE%w_lam = w0
            TDE%wa  = wa
            TDE%cs2_lam  = cs2
        end select



        call P%SetNeutrinoHierarchy(0.00064_dl, 0._dl, 3.046_dl, neutrino_hierarchy_normal)


        P%WantTransfer= .true.
        P%WantCls = .false.
        P%WantScalars = .false.
        P%WantTensors = .false.
        P%WantVectors = .false.

        maxkh = max(PI2/boxsize*rng*maxkhfactor,100.)
        maxkh = min(maxkh,50000.)


        P%ombh2 = omepb * hubble**2
        P%omch2 = (omep-omepb) * hubble**2
        P%omk = 0._dl
        P%TCMB = COBE_CMBTemp
        P%H0      = hubble*100._dl

        P%OutputNormalization = 1

    
        select type(InitPower=>P%InitPower)
        class is (TInitialPowerLaw)
           InitPower%As = As
           InitPower%ns = npower
           InitPower%r = 1
        end select

        


        !these settings seem good enough for sigma8 to a percent or so
        P%Transfer%high_precision=.true.
        P%Transfer%kmax=maxkh
        P%Transfer%k_per_logint=10
        P%Transfer%PK_num_redshifts=1
        P%Transfer%PK_redshifts(1)=redshift


        OK= CAMBparams_Validate(P)
        if(OK .eq. .false.) then
           print *,'###################################'
           print *,'###################################'
           print *,'Wrong set of parameters in the CAMB'
           print *,'###################################'
           print *,'###################################'
           print *,'Ok = ',OK
           stop
        endif


        call CAMB_GetResults(OutData,P) 
!       call CAMB_GetTransfers(P,OutData)
!       call Transfer_SaveMatterPower(OutData%MT,OutData,outfile)
             
!        This cubic box multiplication is to normalize the amplitude in the unit of 2pi*u not k.
        pamp = (OutData%MT%sigma_8(P%Transfer%PK_num_redshifts))**2*cubicbox

      end subroutine InitPower

!--------
      subroutine rInitPower(omep,omepb,omeplam,w0,wa,cs2,hubble,npower,&
                    redshift, As, pamp,&
                  boxsize,rng,readflag,infile,myid)
        use MyPkModule
        implicit none
        integer i,j,k
        integer, intent(in) :: readflag,myid
        character(LEN=80), intent(in) :: infile
        double precision cubicbox
        real amax,anow,ai,bini,bnow,growth,ns,omep,omepb,omeplam,wlam
        real hubble,npower, redshift,  As,boxsize, rng
        real sigma8
        real, intent(out) :: pamp
        real w0,wa,cs2
        external growth

        cubicbox = boxsize**3

        h100 = hubble


!       if(myid.eq.0) then
        open(1,file=infile, form='formatted')
        read(1,*)
        read(1,*)
        read(1,*) sigma8
        read(1,*)
        i = 1
10      read(1, *, end=20) rk(i), pk(i)
        i = i + 1
        goto 10
20      close(1)
        nkh = i -1
!       endif

        pamp = sigma8*cubicbox ! from k space to u space normalization

        if(myid.eq.0) then
        write(*,*)
        write(*,*)
        print *, '#####################################################'
        print *, '####      POWER SPECTRUM PARAMETER READED  ##########'
        print *, '#####################################################'
        print *, 'omep                  = ',omep
        print *, 'omepb                 = ',omepb
        print *, 'omeplam               = ',omeplam
        print *, 'w0(eos of Vac)        = ',w0
        print *, 'wa(eos of Vac)        = ',wa
        print *, 'Cs2(Sound Speed of Vac)= ',cs2
        print *, 'npower                = ',npower
        print *, 'Hubble (km/sec.)      = ',hubble*100
        print *, '#####################################################'
        write(*,*)
        write(*,*)
        endif
! pamp should be determined outside of this routine
      end subroutine rInitPower



!--------
      subroutine GetPower(minkh,dlnkh,matterpk,nmatterpk,Ts)
        use MyPkModule
        real*8 twopi, pi
        parameter(twopi=3.1415926535d0*2.d0,pi=3.141592653d0)
        real minkh,dlnkh,logminkh
        integer nmatterpk
        real matterpk(*)
        double precision logmink,h,k
        real Ts(*)
        minkh = m_minkh
        dlnkh = m_dlnkh
        nmatterpk = nkh
        logmink = log(minkh)
        h = h100
        do i = 1, nkh
           k = exp(logmink+dlnkh*(i-1))*h
           matterpk(i) = pk(i)
           Ts(i) = matterpk(i)/(k*pi*twopi*h**3)
        enddo

      end subroutine GetPower

!--------
      subroutine kjhan_Transfer_SaveMatterPower(MTrans, TState,FileNames, R, minkh,dlnkh,outpower,points, &
            MTrans0, TState0, outpower0)
        use DarkEnergyInterface
        use DarkEnergyFluid
        use DarkEnergyPPF
        use Quintessence
        use MyModule
        use constants
      integer points
      !Export files of total  matter power spectra in h^{-1} Mpc units, against k/h.
      Type(MatterTransferData), intent(in) :: MTrans
      Type(CAMBdata) :: TState
      real :: outpower(points,6)
      Type(MatterTransferData), intent(in) :: MTrans0
      Type(CAMBdata) :: TState0
      real :: outpower0(points,6)

      character(LEN=80), intent(IN) :: FileNames(*)
      character(LEN=name_tag_len) :: columns(3)
      integer itf, i, unit
      real minkh,dlnkh
      real(dl), intent(in):: R
      Type(MatterPowerData) :: PK_data
      integer ncol
      logical all21
      integer :: var1, var2
      !JD 08/13 Changes in here to PK arrays and variables
      integer itf_PK
      double precision sigma8,sigmaR2,ww,rk,drk,rkR
      real omep,omepb,omeplam,w0,wa,hubble, redshift,cs2

      hubble = TState%CP%H0/100._dl
      omepb = TState%CP%ombh2/hubble**2
      omep = TState%CP%omch2/hubble**2 + omepb
      omeplam = 1-omep
      select type(TDE=>P%DarkEnergy)
      class is (TDarkEnergyPPF)
            w0 = TDE%w_lam 
            wa = TDE%wa 
            cs2 = TDE%cs2_lam  
      class is (TDarkEnergyFLUID)
            w0 = TDE%w_lam 
            wa = TDE%wa 
            cs2 = TDE%cs2_lam 
      end select
      redshift = TState%CP%Transfer%PK_redshifts(1)

      all21 = .false.


      do itf=1, TState%CP%Transfer%PK_num_redshifts
        if (FileNames(itf) /= '') then
            if (.not. transfer_interp_matterpower ) then
                itf_PK = TState%PK_redshifts_index(itf)

                points = MTrans%num_q_trans
!               allocate(outpower(points,ncol))

                !Sources
                if (all21) then
                    call Transfer_Get21cmPowerData(MTrans, TState, &
                    PK_data, itf_PK)
                else
                    call Transfer_GetMatterPowerData(TState, MTrans, &
                    PK_data, itf_PK)
                    !JD 08/13 for nonlinear lensing of CMB + LSS compatibility
                    !Changed (CP%NonLinear/=NonLinear_None) to CP%NonLinear/=NonLinear_none .and. CP%NonLinear/=NonLinear_Lens)
                    if(TState%CP%NonLinear/=NonLinear_none .and. &
                              TState%CP%NonLinear/=NonLinear_Lens) then
                        call TState%CP%NonLinearModel%GetNonLinRatios(&
                            TState, PK_data)
                        PK_data%matpower = PK_data%matpower +  &
                                2*log(PK_data%nonlin_ratio)
                        call MatterPowerdata_getsplines(PK_data)
                    end if
                end if

                outpower(:,1) = exp(PK_data%matpower(:,1))
                !Sources
                if (all21) then
                    outpower(:,3) = exp(PK_data%vvpower(:,1))
                    outpower(:,2) = exp(PK_data%vdpower(:,1))

                    outpower(:,1) = outpower(:,1)/1d10*const_pi* &
                     const_twopi/MTrans%TransferData(Transfer_kh,:,1)**3
                    outpower(:,2) = outpower(:,2)/1d10*const_pi* &
                     const_twopi/MTrans%TransferData(Transfer_kh,:,1)**3
                    outpower(:,3) = outpower(:,3)/1d10*const_pi* &
                     const_twopi/MTrans%TransferData(Transfer_kh,:,1)**3
                end if

                call MatterPowerdata_Free(PK_Data)
                columns = ['P   ', 'P_vd','P_vv']
                unit = open_file_header(FileNames(itf), 'k/h', &
                    columns(:ncol), 15)
                do i=1,points
                    write (unit, '(*(E15.6))') MTrans%TransferData( &
                        Transfer_kh,i,1),outpower(i,:)
                end do
                close(unit)
!               deallocate(outpower)
            else
                if (all21) stop 'Transfer_SaveMatterPower: &
                    if output all assume not interpolated'
!               allocate(outpower(points,5))
!               this is for the total power spectrum
#if defined (_CDM_PK_)
#warning "cdm power spectrum is selected to save"
!               this is for the cdm power spectrum
                call Transfer_GetMatterPowerS(TState, MTrans, &
                    outpower(:,1), itf, minkh,dlnkh, points,Transfer_cdm,Transfer_cdm)
#else
!#warning "total matter power spectrum is selected to save"
!  ifundef then, total matter
                call Transfer_GetMatterPowerS(TState, MTrans, &
                    outpower(:,1), itf, minkh,dlnkh, points)
                call Transfer_GetMatterPowerS(TState, MTrans, &
                    outpower(:,2), itf, minkh,dlnkh, points, Transfer_cdm, Transfer_cdm)
                call Transfer_GetMatterPowerS(TState, MTrans, &
                    outpower(:,3), itf, minkh,dlnkh, points, Transfer_b, Transfer_b)
                ! to get growth rate for cdm
                call Transfer_Getgrowthrate(TState, MTrans, &
                    outpower(:,4), itf, minkh,dlnkh, points, Transfer_clxcdot)
                ! to get growth rate for baryon
                call Transfer_Getgrowthrate(TState, MTrans, &
                    outpower(:,5), itf, minkh,dlnkh, points, Transfer_clxbdot)
! this is for the z=0 cdm/baryon powers
                call Transfer_GetMatterPowerS(TState0, MTrans0, &
                    outpower0(:,2), itf, minkh,dlnkh, points, Transfer_cdm, Transfer_cdm)
                call Transfer_GetMatterPowerS(TState0, MTrans0, &
                    outpower0(:,3), itf, minkh,dlnkh, points, Transfer_b, Transfer_b)
#endif

                columns(1) = 'P'
                write(6,'("Now dumping P(k) to ", a20,"  ", i0)') FileNames(itf), itf
                unit = 1
                open(unit,file=FileNames(itf), form='formatted') 
                sigma8 = (MTrans%sigma_8(TState%CP%Transfer%PK_num_redshifts))**2

                sigmaR2 = 0
                do i = 2, points-1
                   rk = minkh*exp((i-1)*dlnkh)
                   drk = minkh*(exp((i-0.5)*dlnkh) - exp((i-1.5)*dlnkh))
                   rkR = rk*R
                   ww = (3.d0 /rkR**3 *(sin(rkR)-rkR*cos(rkR)))**2
                   sigmaR2 = sigmaR2 + outpower(i,1)*ww*rk**2/2.d0/const_pi**2 * drk
                enddo
                

                write(unit,'("Om0=", g14.7," Ob0=", g14.7, &
                 " Ode0=", g10.3," H0=", g12.5," red=", g10.3 )') omep,omepb,omeplam, hubble,redshift
                write(unit,'("w0=", g14.7," wa=", g14.7, &
                 " c_s^2(DE)=", g10.3)') w0,wa,cs2

                write(unit,'(g14.7," =sigma8 at R=8Mpc/h & sigma(R)=", g14.7, &
                 " at Rsmooth= ", g10.3," Mpc/h" )') sigma8, sigmaR2,  R

#if defined (_CDM_PK_)
                write(unit, *)'#   k/h        Pcdm(k)'
                do i=1,points
                    write (unit, '(*(E15.6))') minkh*exp((i-1)*dlnkh),outpower(i,1)
                end do
#else
                write(unit, *)'#   k/h      Ptotm(k)   Pcdm(k)   Pbaryon(k)   f_cdm   f_baryon  D_tot  D_c   D_b'
                do i=1,points
                    write (unit, '(*(G15.6))') minkh*exp((i-1)*dlnkh),outpower(i,1),&
                        outpower(i,2), outpower(i,3), outpower(i,4),outpower(i,5), &
                        sqrt(outpower(i,1)/outpower0(i,1)), &
                        sqrt(outpower(i,2)/outpower0(i,2)), &
                        sqrt(outpower(i,3)/outpower0(i,3))
                end do
#endif
                close(unit)
!               deallocate(outpower)
            end if

        end if
       end do

       end subroutine kjhan_Transfer_SaveMatterPower



!--------
      subroutine kjhan_SavePoissonCorrFactor(MTrans, TState, kh, poissonfact, nkh)
        use DarkEnergyInterface
        use DarkEnergyFluid
        use DarkEnergyPPF
        use Quintessence
      use MyModule
      use constants
      !Export files of total  matter power spectra in h^{-1} Mpc units, against k/h.
      Type(MatterTransferData), intent(in) :: MTrans
      Type(CAMBdata) :: TState
      integer, intent(out):: nkh
      real, dimension(:) :: poissonfact, kh
      character(LEN=name_tag_len) :: columns(3)
      integer itf, i, unit
      integer points
      real minkh,dlnkh
      Type(MatterPowerData) :: PK_data
      integer ncol
      logical all21
      integer :: var1, var2
      !JD 08/13 Changes in here to PK arrays and variables
      integer itf_PK
      double precision sigma8,sigmaR2,ww,rk,drk,rkR
      real omep,omepb,omeplam,w0,wa,hubble, redshift,cs2

      all21 = .false.

      do itf=1, TState%CP%Transfer%PK_num_redshifts
          nkh = MTrans%num_q_trans
          print *,'nkh= ', nkh
          i = TState%PK_redshifts_index(itf)
          print *,'i= ', i
          do ik = 1, nkh
             print *, 'values ', Transfer_kh,ik,i
             print *, 'contents kh', MTrans%TransferData(Transfer_kh,ik,i)
             print *, 'contents pf', MTrans%TransferData(Transfer_poisson_factor,ik,i)
!            kh(ik) = MTrans%TransferData(Transfer_kh,ik,i)
!            poissonfact(ik) = MTrans%TransferData(Transfer_poisson_factor,ik,i)
          enddo
      end do

      end subroutine kjhan_SavePoissonCorrFactor
   


!----------------------------------------------------------------------------------

!     https://cosmocoffee.info/viewtopic.php?f=11&t=212&p=548&hilit=poisson#p548
MODULE PoissonCorrection
    use, intrinsic :: ISO_C_BINDING
    use CAMB
    use DarkEnergyInterface
    use DarkEnergyFluid
    use DarkEnergyPPF
    use Quintessence
    use MyModule
    implicit none
    contains

    subroutine getPoissonCfactor(nx,boxsize, hubble, npower, omep,omepb, omeplam, &
        C_DarkEnergyModel, w0,wa, cs2_lam,  As,rng, amax, anow, nkh, kh, poissonfact)  &
        bind (C, name='Fgetpcorr')
!   use lensing
!   use constants
!   use Bispectrum
!   use CAMBmain
!   use Transfer

    logical :: call_again,background_only
    real(c_float) w0,wa,omep,omepb,omeplam,hubble, npower, maxkh, boxsize,As,rng
    real(c_float) cs2_lam,amax,anow,maxkhfactor
    integer (c_int) :: nx
    character (len=1, kind=c_char), intent (in) ::C_DarkEnergyModel(20)
    character (len=20) ::DarkEnergyModel
    integer (c_int), intent(inout) :: nkh
    real (c_float),  intent(out) :: poissonfact(nkh), kh(nkh)


    real(dl) PI2, redshift
    parameter(PI2=3.1415926535d0*2.d0,maxkhfactor=5)
    logical OK
    integer error,i,ik, i_PK, itf
    Type(MatterTransferData), pointer :: MTrans

    global_error_flag = 0
    redshift = amax/anow-1

    DarkEnergyModel = transfer(C_DarkEnergyModel, DarkEnergyModel)



    call_again = .false.
    call CAMB_SetDefParams(P)
    call SetActiveState(OutData) !  making pointers as State => OutData &&  CP => OutData%CP
    if(allocated(P%DarkEnergy)) deallocate(P%DarkEnergy)
    if(DarkEnergyModel(1:5) == 'FLUID') then
        allocate(TDarkEnergyFluid::P%DarkEnergy)
        if(w0==-1 .and. wa==0) then
           P%DarkEnergy%is_cosmological_constant = .true.
           P%DarkEnergy%num_perturb_equations = 0
        else
           P%DarkEnergy%is_cosmological_constant = .false.
           P%DarkEnergy%num_perturb_equations = 2
        endif
     else if(DarkEnergyModel(1:3) == 'PPF') then
        allocate(TDarkEnergyPPF::P%DarkEnergy)
        if(w0==-1 .and. wa==0) then
           P%DarkEnergy%is_cosmological_constant = .true.
           P%DarkEnergy%num_perturb_equations = 0
        else
           P%DarkEnergy%is_cosmological_constant = .false.
           P%DarkEnergy%num_perturb_equations = 1
        endif
     else if(DarkEnergyModel(1:19) == 'AXIONEFFECTIVEFLUID') then
        allocate(TAxionEffectiveFluid::P%DarkEnergy)
        P%DarkEnergy%num_perturb_equations = 2
     else if (DarkEnergyModel(1:17) == 'EARLYQUINTESSENCE') then
        allocate(TQuintessence::P%DarkEnergy)
        P%DarkEnergy%num_perturb_equations = 2
     endif


     select type(TDE=>P%DarkEnergy)
     class is (TDarkEnergyPPF)
         TDE%w_lam = w0
         TDE%wa  = wa
         TDE%cs2_lam  = cs2_lam
     class is (TDarkEnergyFluid)
         TDE%w_lam = w0
         TDE%wa  = wa
         TDE%cs2_lam  = cs2_lam
     end select
     call P%SetNeutrinoHierarchy(0.00064_dl, 0._dl, 3.046_dl, neutrino_hierarchy_normal)



    P%WantTransfer = .true.
    P%WantCls = .false.
    P%WantScalars = .false.
    P%WantTensors = .false.
    P%WantVectors = .false.

    maxkh = max(PI2/boxsize*rng*maxkhfactor,500.)
    maxkh = min(maxkh,50000.)


    P%ombh2 = omepb * hubble**2
    P%omch2 = (omep-omepb) * hubble**2
    P%omk = 0._dl
    P%TCMB = COBE_CMBTemp
    P%H0      = hubble*100._dl

    P%OutputNormalization = 1

    OutData%flat = .true.
    OutData%OnlyTransfer = .true.



    select type(InitPower=>P%InitPower)
    class is (TInitialPowerLaw)
       InitPower%As = As
       InitPower%ns = npower
       InitPower%r = 1
    end select

    !these settings seem good enough for sigma8 to a percent or so
    P%Transfer%high_precision=.true.
    P%Transfer%kmax=maxkh
    P%Transfer%k_per_logint=10
    P%Transfer%PK_num_redshifts=1
    P%Transfer%PK_redshifts(1)=redshift




    OutData%flat = .true.
    OK= CAMBparams_Validate(P)
    if(.not. OK) then
        print *,'Problem in validating the parameter file '
    endif

    ! copy parameter P contents to OutData%CP for calculation
    call CAMBdata_SetParams(OutData, P, error, .false., .false., .true.)


    call CAMB_GetResults(OutData,P) ! it calls cmbmain() internally.

    do itf=1, OutData%CP%Transfer%PK_num_redshifts
          nkh = OutData%MT%num_q_trans
          i = OutData%PK_redshifts_index(itf)
          do ik = 1, nkh
             kh(ik) = OutData%MT%TransferData(Transfer_kh,ik,i)
             poissonfact(ik) = OutData%MT%TransferData(Transfer_poisson_factor,ik,i)
          enddo
    end do



    end subroutine getPoissonCfactor
end MODULE PoissonCorrection



