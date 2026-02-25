     !Simple test program to print out sigma_8 as a function of the CDM density
!     module MyModule
!        use CAMB
!        type (CAMBdata) OutData,InData
!     end module MyModule

      program makecambpower
        use CAMB
        use DarkEnergyInterface
        use DarkEnergyFluid
        use DarkEnergyPPF
        use Quintessence
        use MyModule
        use handles
        use, intrinsic :: ISO_C_BINDING
        implicit none
        integer i,j,k
        real PI2,rng,maxkhfactor,w0,wa
        real cs2_lam
        real sigma8
        parameter(PI2=3.1415926535d0*2.d0,maxkhfactor=5)
        type(MatterTransferData) MTrans
        real omep,omepb,omeplam,hubble,npower,maxkh,redshift,boxsize
        integer itf,in,narg
        integer error
        real minkh,dlnkh,pamp,bias,As,khmax,opamp
        real amax,anow,da
        real (dl):: Rsmooth
        integer npoints,nmatterpk
        character(LEN=80) fmt,inputfile, infile
        character(LEN=80) :: outfile, transferfile
        character(LEN=30) DarkEnergyModel
        character(LEN=1, KIND=c_char) C_DarkEnergyModel(30)
        character(LEN=1, KIND=c_char) C_outfile(80)
        logical OK
        double precision cubicbox
        integer ntimes,nq,noutputs,ncustomsources
        TYPE(C_FUNPTR) :: c_source_func
        parameter(ntimes=2, noutputs=Transfer_Max+9, ncustomsources=0)
        real (dl):: times(ntimes)
        real (dl), allocatable:: q(:), outputs(:, :, :)
        real, allocatable :: outpower(:,:)
        type(CAMBdata) OutData0
        real, allocatable :: outpower0(:,:)



        narg = iargc()
        if(narg .ne. 1) then
           print *,'error occurred. Please input the parameter file'
           stop
        endif
        call getarg(1,inputfile)
        infile = trim(inputfile)//char(0)
        DarkEnergyModel = ' '
        outfile = ' '
        call getwsimparm(infile, boxsize,hubble,npower,omep,omepb, &
              omeplam, C_DarkEnergyModel, w0,wa,cs2_lam,bias,As,rng, &
              amax,da,&
              anow,i,c_outfile)
        DarkEnergyModel = transfer(C_DarkEnergyModel,DarkEnergyModel)
        outfile = transfer(C_outfile,outfile)

        if(DarkEnergyModel .eq. 'PPF') then
           print *, 'DarkEnergyModel= PPF'
        else if(DarkEnergyModel .eq. 'FLUID') then
           print *, 'DarkEnergyModel= FLUID'
        endif
        print *, 'w0(DE)= ',w0,'wa(DE)= ', wa,'cs2(DE)= ',cs2_lam
        print *, 'As= ', As
        redshift = amax/anow - 1
        print *, 'redshift= ', redshift
       
        cubicbox = boxsize**3

        call CAMB_SetDefParams(P)
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
           print *, 'PPF model is selected'
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
        class is (TDarkEnergyPPF )
            TDE%w_lam = w0
            TDE%wa  = wa
            TDE%cs2_lam  = cs2_lam
        class is (TDarkEnergyFluid)
            TDE%w_lam = w0
            TDE%wa  = wa
            TDE%cs2_lam  = cs2_lam
        end select


        call P%SetNeutrinoHierarchy(0.00064_dl, 0._dl, 3.046_dl, neutrino_hierarchy_normal)


        P%WantTransfer= .true.
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
        else
           print *,'Congratulations! Parameter check has been passed.'
!          print *,'output file name is ', trim(outfile)
        endif
! This is for z=0 power
        P%Transfer%PK_redshifts(1)=0
        call CAMBdata_SetParams(OutData0, P, error, .false., .false., .true.)
        call CAMB_GetResults(OutData0,P)
        minkh = 1.e-4 ! fixed in CAMB
        dlnkh = 0.02 ! fixed in CAMB
        npoints = log(OutData0%MT%TransferData(Transfer_kh,&
            OutData0%MT%num_q_trans,1)/minkh)/dlnkh+1
        allocate(outpower0(npoints,6))


! This is for target redshift
        ! copy parameter P contents to OutData%CP for calculation
        P%Transfer%PK_redshifts(1)=redshift
        call CAMBdata_SetParams(OutData, P, error, .false., .false., .true.)
        call CAMB_GetResults(OutData,P)
        print *, 'the global error : ', global_error_flag
        sigma8 = (OutData%MT%sigma_8(P%Transfer%PK_num_redshifts))**2
        print *, 'sigma8 = ', sigma8

        Rsmooth = 130._dl
!--     set the bin of kh and output power
        minkh = 1.e-4 ! fixed in CAMB
        dlnkh = 0.02 ! fixed in CAMB
        npoints = log(OutData%MT%TransferData(Transfer_kh,&
            OutData%MT%num_q_trans,1)/minkh)/dlnkh+1
        allocate(outpower(npoints,6))
        call kjhan_Transfer_SaveMatterPower(OutData%MT,OutData, &
            outfile, Rsmooth, minkh,dlnkh, outpower, npoints, &
            OutData0%MT, OutData0, outpower0)

        transferfile = 'transfer.'//trim(outfile)
        call Transfer_SaveToFiles(OutData%MT,OutData,transferfile)


       
       deallocate(outpower0)
       deallocate(outpower)
      end program makecambpower

