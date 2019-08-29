module magnetic
  
  use constants
  use AMLutils
  use ModelParams
  use CAMB
  use richsub
  use Transfer
  use InitialPower  !ALEX.

  ! Module for calculating the amplitudes of the magnetic perturbations
  ! to the energy momentum tensor from the magnetic field properties.
  
  implicit none

  !private
  
    character(len= 100) :: mi_filename = 'magnetic_integral_new.dat'
    ! column and row , contains bothe helical and non-helcial assuming infinity k,integral
    integer, parameter :: mi_nr = 120
    integer, parameter :: mi_nc = 13
!logical :: maximal_helical = .false.  !real(dl) :: AB, AH
  ! Array of the spectral indices
  real(dl) :: specind(mi_nr)
  ! Array of the calculated angular integrals
  ! In order: delta, delta-pi, pi, pi-vector, pi-tensor, h-delta,h-delta-pi,h-pi,h-pi-vec,h-pi-tens,h-pi-vec-odd,h-pi-tens-odd
  real(dl) :: intval(mi_nr,mi_nc-1)
  ! Array of second derivatives for splining
  real(dl) :: intval2(mi_nr,mi_nc-1)

  ! Array of constants for each PS
  real(dl), parameter :: psconst(mi_nc-1) = (/0.25, 0.25, 0.25, 2., 2./3., -0.5, 0.25,-0.25, -2., 4./3., -2., -2./3./)

  ! Constant for calculating the amplitudes = 1/(2*rho_gamma_0)^2
  real(dl), parameter :: amp_constant = 1.432e-12_dl

  !Define constants for specifying type of perturbation
  !integer, parameter :: mag_scal_passive = 1, mag_scal_comp = 2, &
  !     mag_vec_comp = 3, mag_tens_passive = 4, mag_tens_comp = 5
  integer, parameter :: mag_compensated = 1, mag_passive = 2!, mag_and_prim = 3

  logical :: notloaded = .true.
  
contains

  ! Initialise the module. Read files etc. contains both helical and non-helical
  subroutine magamp_init

    implicit none

    integer :: i

    open(217,file=mi_filename)

    do i=1,mi_nr
       read(217,*) specind(i), intval(i,:)
    end do

    do i=1,mi_nc-1
       call spline_double(specind,intval(:,i),mi_nr,intval2(:,i))
    end do

    close(217)

    notloaded = .false.

  end subroutine magamp_init
  
! Get integral values by splining loaded files. contains both helical and non-helical
  real(dl) function mag_spline_int(mag_ind, corr_num)
    real(dl), intent(in) :: mag_ind
    integer, intent(in) :: corr_num

    if(notloaded) call magamp_init

    mag_spline_int = spline_val(mag_ind, specind, &
         intval(:,corr_num), intval2(:,corr_num), mi_nr)

  end function mag_spline_int

  !================ Non_Helical Parts ===================
  ! Function to calculate the common ampltiude
  ! B_lambda is in units of nano Gauss.
  real(dl) function mag_amplitude(spec_ind, b_lambda)
    
    !use AMLutils
    real(dl), intent(in) :: spec_ind
    real(dl), intent(in) :: b_lambda
    
    
    ! Scale by the spectral index dependent bits, and the field amplitude
    mag_amplitude = amp_constant * b_lambda**4 * (2*pi)**(2*spec_ind + 4) &
         / GAMMA((spec_ind + 3)/2._dl)**2

  end function mag_amplitude
    
    
  ! Function to calculate the amplitude of a particular spectrum
  real(dl) function mag_psamp(spec_ind, b_lambda, corr_num)

    real(dl), intent(in) :: spec_ind
    real(dl), intent(in) :: b_lambda
    integer, intent(in) :: corr_num
    write(*,*) "non-helical table Int and table int value = ", mag_amplitude(spec_ind, b_lambda) *&
    mag_spline_int(spec_ind, corr_num)*psconst(corr_num) !useless,check
    mag_psamp = mag_amplitude(spec_ind, b_lambda) * &
         mag_spline_int(spec_ind, corr_num) * psconst(corr_num)
  end function mag_psamp

!============== Helical part ==============================
 !Amplitude
real(dl) function mag_amplitude_hel(spec_ind, b_lambda)
    !use AMLutils
    real(dl), intent(in) :: spec_ind
    real(dl), intent(in) :: b_lambda

    ! Scale by the spectral index dependent bits, and the field amplitude
    mag_amplitude_hel = amp_constant * b_lambda**4 * (2*pi)**(2*spec_ind + 4) &
                        / GAMMA((spec_ind + 4)/2._dl)**2

end function mag_amplitude_hel

 ! Function to calculate the amplitude of a particular spectrum
real(dl) function mag_psamp_hel(spec_ind, b_lambda, corr_num)
    real(dl), intent(in) :: spec_ind
    real(dl), intent(in) :: b_lambda
    integer, intent(in) :: corr_num
    write(*,*) "helical table Int and table int value = ", mag_amplitude_hel(spec_ind, b_lambda) * &
    mag_spline_int(spec_ind, corr_num)*psconst(corr_num) !useless,check
    mag_psamp_hel = mag_amplitude_hel(spec_ind, b_lambda) * &
    mag_spline_int(spec_ind, corr_num) * psconst(corr_num)
end function mag_psamp_hel

!======== This might be useless...

  ! Set the parameters needed for adding in the Alfven velocity to get
  ! a magnetic Jeans mass. Must be done separately so we can use for
  ! the non magnetic modes.
  subroutine mag_set_alfven(mag_amp, mag_ind)
    real(dl), intent(in) :: mag_amp, mag_ind

    klambda = 2._dl * pi
    blambda = mag_amp
    nb_ind = mag_ind

  end subroutine mag_set_alfven

  ! Reset back to the default parameters.
  subroutine mag_reset_alfven
    call mag_set_alfven(0._dl, -3._dl)
  end subroutine mag_reset_alfven



 subroutine cls_from_params(CP, Cls, usedata)

   ! Barely needed, only for call to ClsFromTheoryData to get
   ! normalisation

   ! Contains all the parameters for CAMB
   type(CAMBParams), intent(in) :: CP
   
   ! Array to return Cls into
   real, intent(out) :: Cls(lmin:,:)

   logical, intent(in), optional :: usedata

   logical :: ud
   
   type(CAMBdata) :: data
   type(CAMBParams) :: tcp   

   integer :: error, l

   integer, parameter :: ScalClOrder(3) = (/1,3,2/), TensClOrder(4) = (/1,4,2,3/)
   
   ud = .false.
   if(present(usedata)) ud = usedata

   tcp = CP

   if(ud) then 
      ! Initialise data type for holding transfers etc.
      call CAMB_InitCAMBdata(data)
      data%Params = tcp
   
      ! Do hard work in calculating
      call CAMB_GetTransfers(tcp, data, error)
       !Print *, "ERRORS :", error!
      ! Calculate power spectra
      call CAMB_TransfersToPowers(data)

   else
      call CAMB_GetResults(tcp, error)
   end if
   !check
   if(error /= 0) then
      write (*,*) "Eeek error."
      return
   end if


   ! Copy out Cls into the return array.
   if(CP%WantScalars) then

      if(CP%DoLensing) then
         Cls(lmin:CP%Max_l,C_Temp:C_last) =  Cl_lensed(lmin:CP%Max_l,1, C_Temp:C_last)
      else 
         Cls(lmin:CP%Max_l,C_Temp:C_last) =  Cl_scalar(lmin:CP%Max_l,1, C_Temp:C_last)
      end if

   else if(CP%WantTensors) then
      do l = 2, CP%Max_l_tensor
         Cls(l,1:4) = Cl_tensor(l,1, TensClOrder(1:4))
      end do
   else if(CP%WantVectors) then
      do l = 2, CP%Max_l
         Cls(l,1:4) =  Cl_vector(l,1, TensClOrder(1:4))
      end do
   end if
   
   ! Free up everything created.
   if(ud) call CAMB_FreeCAMBdata(data)
   
 end subroutine cls_from_params



 ! Calculate the magnetic Cls 
 ! Largely this sets up all the relevant parameters and then calls
 ! cls_from_params, returning the results.
 ! Adding the Helical contributions
  subroutine mag_cls(CP, magmode, magamp, magind, maglrat, helical_amp, helical_ind, usedata)

    type(CAMBParams) ,intent(inout) :: CP

    integer, intent(in) :: magmode
    double precision, intent(in) :: magamp, magind, maglrat
    !logical, intent(in) :: do_helical
    double precision, intent(in) :: helical_amp, helical_ind
    logical, intent(in), optional :: usedata

!other things, check 
    real(dl) :: t, DD, DP, PP, DDH, DPH, PPH
    integer :: i_t
    
    double precision :: amp, camb_ind

    double precision :: Rg, Rv

    double precision :: lrat
    real :: Cltemp(lmin:CP%Max_l,1:4)
    real :: Cltemp2(lmin:CP%Max_l,1:4)
    integer :: i
    logical :: ud

    Rg = 0.0_dl; Rv = 0._dl

    ud = .false.
    if(present(usedata)) ud = usedata

    ! Set photon and neutrino fractions
    Rg = 1._dl / (1._dl+7._dl/8. * (4._dl/11.)**(4._dl/3.) * (CP%Num_Nu_massless + CP%Num_Nu_massive))
    Rv = 1._dl - Rg

    if(magmode == 0 .or. .not. CP%WantCls) then
       write (*,*) "Can't calculate Mag Cls if mode = 0, or WantCls is false."
       return
    end if

    !Alex: COMMON STUFF for FITTING FUNCTIONS
    CP%InitPower%kDissip = (5.5d0 * 10.d0**4)**(1.d0/(magind + 5.d0)) * (magamp)**(-2.d0/(magind+5.d0)) * &
    (2*pi)**((magind+3)/(magind+5))*(CP%H0/100)**(1/(magind+5))&
    *((CP%omegab*(CP%H0/100)**2)/0.022)**(1/(magind+5))

    write(*,*) "kD = ", CP%InitPower%kDissip


    ! Set common power spectrum stuff
    CP%InitPower%nn = 1
    CP%InitPower%rat(1) = 1._dl 
    CP%InitPower%n_run(1) = 0._dl
    CP%InitPower%nt_run(1) = 0._dl !tensor 
    CP%InitPower%k_0_scalar = 2*pi
    CP%InitPower%k_0_tensor = 2*pi

    ! Setup common options
    CP%WantCls = .true.
    CP%WantTransfer = .false.
    !CP%WantScalars = .false.
    !CP%WantVectors = .false.
    !CP%WantTensors = .false.
    CP%OnlyTransfers = .false.

    ! Convert logarithm into natural log
    lrat = maglrat * log(10._dl)

    if(CP%WantTensors) then
! Ensure tensor neutrinos are on (as mag field modes don't make
! any sense otherwise)
     DoTensorNeutrinos = .true.
!AZ: Set tensor_parameterization = 3. check
     CP%InitPower%tensor_parameterization = 3
        
     if(magmode == mag_passive) then
        CP%WantTensors = .true.
        delb = 0._dl ! Set up perturbations
        pib = 0._dl
       
            if(magind < -2.d0) then
               !if(Feedbacklevel>1) 
               write(*,*) "USING INTERPOLATION TABLE for non-helical P_tpp"
               CP%InitPower%ant(1) = 2._dl*(3 + magind)
               CP%InitPower%TensorPowerAmp(1) = mag_psamp(magind, magamp, 5) *&
                                                (3._dl*Rg*(lrat + (5._dl/(8._dl*Rv) - 1)))**2!yun
               CP%InitPower%CorrType = 0
            else       
!Use first fitting formula
               !if(Feedbacklevel>1) 
               write(*,*) "USING FITTING FORMULAS for non-helical P_tpp"
               CP%InitPower%ant(1) = magind
               CP%InitPower%TensorPowerAmp(1) = mag_amplitude(magind, magamp)*psconst(5)*&
                                                (3._dl*Rg*(lrat+(5._dl/(8._dl*Rv)-1)))**2
               CP%InitPower%CorrType = 6
            end if
!tensor helical passive
            if(do_helical) then
                if(helical_ind < -2.d0) then
                    CP%InitPower%CorrType_hel = 0
                    CP%InitPower%TensorPowerAmp_hel(1) = mag_psamp_hel(helical_ind, helical_amp,10)&
                                                       *(3._dl*Rg*(lrat + (5._dl/(8._dl*Rv) - 1)))**2 
                    CP%InitPower%TensorPowerAmp_odd(1) = mag_psamp_hel(helical_ind, helical_amp,12)&
                                                       *(3._dl*Rg*(lrat + (5._dl/(8._dl*Rv) - 1)))**2 
                               
                    CP%InitPower%ant_hel1(1) = 2.d0*(helical_ind +3.d0)
                    write(*,*) "Using interpolation TABLE for helical P_tpp"
                else 
                    write(*,*) "Using fitting functions for helical P_tpp" 
                    CP%InitPower%ant_hel1(1) = helical_ind
                    CP%InitPower%TensorPowerAmp_hel(1) = mag_amplitude_hel(helical_ind, helical_amp)*psconst(10)&
                                                       *(3._dl*Rg*(lrat + (5._dl/(8._dl*Rv) - 1)))**2
                    CP%InitPower%TensorPowerAmp_odd(1) = mag_amplitude_hel(helical_ind, helical_amp)*psconst(12)&
                                                       *(3._dl*Rg*(lrat + (5._dl/(8._dl*Rv) - 1)))**2
                    CP%InitPower%CorrType_hel = 6 
                end if
             end if
                  
     else if(magmode == mag_compensated) then
          CP%WantTensors = .true.
          delb = 0._dl ! Set up perturbations
          pib = 1._dl
          tens0 = 0._dl
  
            if(magind < -2.d0) then
               !if(Feedbacklevel>1) 
               write(*,*) "USING INTERPOLATION TABLE for non-helical P_tpp"
               CP%InitPower%ant(1) = 2._dl*(3 + magind)
               CP%InitPower%TensorPowerAmp(1) = mag_psamp(magind, magamp, 5)
               CP%InitPower%CorrType = 0
               !write(*,*) "non-helical mag_amplitude = ", mag_amplitude(magind, magamp,5) !useless,check
               ! write(*,*) "CP%InitPower%an(1) = ", CP%InitPower%ant(1) !useless,check
                !write(*,*) "CP%InitPower%TensorPowerAmp(1) = ",  CP%InitPower%TensorPowerAmp(1) !useless,check
            else       
            !Use first fitting formula
               !if(Feedbacklevel>1) 
               write(*,*) "USING FITTING FORMULAS for non-helical P_tpp"
               CP%InitPower%ant(1) = magind
               CP%InitPower%TensorPowerAmp(1)=mag_amplitude(magind, magamp)*psconst(5)
               CP%InitPower%CorrType = 7
            end if
            
            if(do_helical) then
                if(helical_ind < -2.d0) then
                    CP%InitPower%CorrType_hel = 0
                    CP%InitPower%TensorPowerAmp_hel(1)=mag_psamp_hel(helical_ind, helical_amp,10) 
                    CP%InitPower%TensorPowerAmp_odd(1)=mag_psamp_hel(helical_ind, helical_amp,12) 
                    CP%InitPower%ant_hel1(1) = 2.d0*(helical_ind +3.d0)
                    write(*,*) "Using interpolation TABLE for helical passive P_tpp"
                else 
                    write(*,*) "Using fitting functions for helical passive P_tpp" 
                    CP%InitPower%ant_hel1(1) = helical_ind
                    CP%InitPower%TensorPowerAmp_hel(1)=mag_amplitude_hel(helical_ind, helical_amp)*psconst(10)
                    CP%InitPower%TensorPowerAmp_odd(1)=mag_amplitude_hel(helical_ind, helical_amp)*psconst(12)
                    CP%InitPower%CorrType_hel = 7
                end if
             end if
     end if   
!sclar    
    else if(CP%WantScalars) then
       ! Setup the parameters for calculating the Alfven velocity
       !This might be useless..
       call mag_set_alfven(magamp, magind)

       if(magmode == mag_compensated) then
          ! This section is quite trickey as we need to combine the two
          ! correlated perturbation types. It essentially fetches each
          ! set of Cls with the correct amplitudes (taking into account
          ! the cross correlation), and then sums them up.         
          ! Set up common parameters

          CP%Scalar_initial_condition = 6
          CP%WantScalars = .true.
! Delta-Delta - Delta-Pi, The density perturbed mode (delb = 1, pib = 0)
          write(*,*) "Delta-Delta - Delta-Pi"
                  
          if(magind < -2.d0) then
                write(*,*) "Using interpolation TABLE for non-helical P_SDD-P_SDPi part"
                CP%InitPower%an(1) = 1._dl + 2._dl*(3 + magind)
                CP%InitPower%ScalarPowerAmp(1)= mag_psamp(magind, magamp, 1) - mag_psamp(magind, magamp, 2)
                CP%InitPower%CorrType = 0
          else
!Use first fitting formula
                write(*,*) "USING FITTING FORMULAS for non-helical P_SDD-P_SDPi part"
                CP%InitPower%an(1) =  magind
                CP%InitPower%ScalarPowerAmp(1) = mag_amplitude(magind, magamp)
                CP%InitPower%CorrType = 1               
          end if

          if(do_helical) then
             if(helical_ind < -2.d0) then
                write(*,*) "Using interpolation TABLE for helical P_ADD-P_ADPi part"
                CP%InitPower%CorrType_hel = 0
                CP%InitPower%ScalarPowerAmp_hel(1)=mag_psamp_hel(helical_ind, helical_amp,6) -&
                                                       mag_psamp_hel(helical_ind, helical_amp,7)
                CP%InitPower%an_hel1(1) =  1.d0 + 2.d0*(helical_ind +3.d0)
             else
!Use first fitting formula
                write(*,*) "USING FITTING FORMULAS for helical P_ADD-P_ADPi part"
                CP%InitPower%an_hel1(1) = helical_ind
                CP%InitPower%ScalarPowerAmp_hel(1) = mag_amplitude_hel(helical_ind, helical_amp)
                CP%InitPower%CorrType_hel = 1
             end if                
          end if
          
          delb = 1._dl ! Set up perturbations
          pib = 0._dl
          call cls_from_params(CP,Cltemp2,ud)
          Cltemp(:,:) = Cltemp(:,:) + Cltemp2(:,:)
!------------------------------------------------The stress perturbed mode (delb = 0, pib = 1)
          write(*,*) "Pi-Pi - Delta-Pi"
          if(magind < -2.d0) then 
                write(*,*) "Using interpolation TABLE for non-helical P_SPiPi-P_SDPi "
                CP%InitPower%ScalarPowerAmp(1)= mag_psamp(magind, magamp, 3) - mag_psamp(magind, magamp, 2)
          else             
                write(*,*) "Using FITTING Functions for non-helical P_SPiPi-P_SDPi "
                CP%InitPower%CorrType = 2   
          end if 

          if(do_helical) then
             if(helical_ind < -2.d0) then
               write(*,*) "Using interpolation TABLE for helical P_APiPi-P_ADPi part"
               CP%InitPower%ScalarPowerAmp_hel(1) = mag_psamp_hel(helical_ind, helical_amp,8) -&
                                                       mag_psamp_hel(helical_ind, helical_amp,7)
             else
!Use first fitting formula
                write(*,*) "USING FITTING FORMULAS for helical P_APiPi-P_ADPi part"
                CP%InitPower%CorrType_hel = 2
              end if               
          end if
   
          delb = 0._dl ! Set up perturbations
          pib = 1._dl
          call cls_from_params(CP,Cltemp2,ud)
          Cltemp(:,:) = Cltemp(:,:) + Cltemp2(:,:)

!-------------------------------------------------
!! The combined perturbation (delb = 1, pib = 1)
          write(*,*) "Delta-Pi"
          if(magind < -2.d0) then 
                write(*,*) "Using interpolation TABLE for non-helical PSDP"
                CP%InitPower%ScalarPowerAmp(1)= mag_psamp(magind, magamp, 2)
          else
                write(*,*) "Using FITTING Functions for non-helical PSDP"
                CP%InitPower%CorrType = 3  
          end if
          
          if(do_helical) then
             if(helical_ind < -2.d0) then
               write(*,*) "Using interpolation TABLE for helical P_ADP part"
               CP%InitPower%ScalarPowerAmp_hel(1) = mag_psamp_hel(helical_ind, helical_amp,7)
             else
!Use first fitting formula
                write(*,*) "USING FITTING FORMULAS for helical P_ADP part"
                CP%InitPower%CorrType_hel = 3
              end if               
          end if
 
          delb = 1._dl ! Set up perturbations
          pib = 1._dl
          call cls_from_params(CP,Cltemp2,ud)
          Cltemp(:,:) = Cltemp(:,:) + Cltemp2(:,:)
          write (*,*) "End"
          ! Copy back into Cl_scalar array
          Cl_scalar(lmin:CP%Max_l, 1, C_Temp:C_last) = Cltemp(lmin:CP%Max_l, C_Temp:C_last)
! Reset the Alfven velocity related stuff
      
          call mag_reset_alfven

          return

       else if(magmode == mag_passive) then

          CP%Scalar_initial_condition = 1
          ! Set up perturbations
          delb = 0._dl
          pib = 0._dl          
          CP%WantScalars = .true.
          
          if(magind < -2.d0) then
                if (Feedbacklevel>1) write(*,*) "Using interpolation TABLE for non-helical passive P_Spp part"
                CP%InitPower%an(1) = 1._dl + 2._dl*(3 + magind)
                CP%InitPower%ScalarPowerAmp(1)= mag_psamp(magind, magamp, 3)*&
                                                (Rg*(lrat + (5._dl/(8._dl*Rv) - 1)))**2
                CP%InitPower%CorrType = 0
          else
!Use first fitting formula
                if (Feedbacklevel>1) write(*,*) "USING FITTING FORMULAS for non-helical passive P_Spp part"
                CP%InitPower%an(1) = magind
                CP%InitPower%ScalarPowerAmp(1)=mag_amplitude(magind, magamp)*psconst(3)*&
                                            (Rg*(lrat + (5._dl/(8._dl*Rv) - 1)))**2
                CP%InitPower%CorrType = 4               
          end if

          if(do_helical) then
             if(helical_ind < -2.d0) then
                write(*,*) "Using interpolation TABLE for helical passive P_App part"
                CP%InitPower%CorrType_hel = 0   
                CP%InitPower%an_hel1(1) = 1.d0 + 2.d0*(3.d0+helical_ind)
                CP%InitPower%ScalarPowerAmp_hel(1)=mag_psamp_hel(helical_ind, helical_amp, 8) &
                                                    * (Rg*(lrat + (5._dl/(8._dl*Rv) - 1)))**2
             else
!Use first fitting formula
                write(*,*) "USING FITTING FORMULAS for helical passive P_App part"
                CP%InitPower%an_hel1(1) = helical_ind
                CP%InitPower%ScalarPowerAmp_hel(1)=mag_amplitude_hel(helical_ind, helical_amp)*psconst(8)*&
                                                    (Rg*(lrat + (5._dl/(8._dl*Rv) - 1)))**2
                CP%InitPower%CorrType_hel = 4
             end if                
          end if
          
       end if 

    else if(CP%WantVectors) then

       if(magmode == mag_compensated) then
         vec_sig0 = 0._dl
         delb = 0._dl ! Set up perturbations
         pib = 1._dl
         CP%WantVectors = .true.
         
           if(magind < -2.d0) then !check
                !if (Feedbacklevel>1) 
                write(*,*) "Using Table for non-helical P_Svpp integrals"
                CP%InitPower%CorrType = 0
                CP%InitPower%an(1) = 1._dl + 2._dl*(3 + magind)
                CP%InitPower%ScalarPowerAmp(1) = mag_psamp(magind, magamp, 4)
                !write(*,*) "non-helical mag_amplitude_hel = ", mag_amplitude(magind, magamp) !useless,check
                !write(*,*) "CP%InitPower%an(1) = ", CP%InitPower%an(1) !useless,check
                !write(*,*) " CP%InitPower%ScalarPowerAmp(1)= ",  CP%InitPower%ScalarPowerAmp(1) !useless,check
                write(*,*) "CP%InitPower%CorrType = ", CP%InitPower%CorrType !useless,check
           else
                CP%InitPower%CorrType = 5
                if (Feedbacklevel>1) write(*,*) "USING FITTING FUNCTION for non-helical P_Svpp integrals"
                CP%InitPower%ScalarPowerAmp(1) = mag_amplitude(magind, magamp)* psconst(4)
                CP%InitPower%an(1) = magind
                write(*,*) "CP%InitPower%CorrType = ", CP%InitPower%CorrType !useless,check
           end if
           
           if(do_helical) then
             if(helical_ind < -2.d0) then !check
                write(*,*) "Using interpolation TABLE for helical P_Avpp part"
                CP%InitPower%CorrType_hel = 0
                CP%InitPower%an_hel1(1) = 1.d0 + 2.d0*(3.d0+helical_ind)
                CP%InitPower%ScalarPowerAmp_hel(1) = mag_psamp_hel(helical_ind, helical_amp, 9)
                CP%InitPower%VectorPowerAmp_odd(1) = mag_psamp_hel(helical_ind, helical_amp, 11)
                !write(*,*) "helical mag_amplitude_hel = ", mag_amplitude_hel(helical_ind, helical_amp) !useless,check
                !write(*,*) "CP%InitPower%an_hel1(1) = ", CP%InitPower%an_hel1(1) !useless,check
                !write(*,*) " CP%InitPower%ScalarPowerAmp_hel(1)= ",  CP%InitPower%ScalarPowerAmp_hel(1) !useless,check
             else
!Use first fitting formula
                write(*,*) "USING FITTING FORMULAS for helical P_Avpp part"
                CP%InitPower%CorrType_hel = 5
                CP%InitPower%an_hel1(1) = helical_ind
                CP%InitPower%ScalarPowerAmp_hel(1) = mag_amplitude_hel(helical_ind, helical_amp)*psconst(9)
                CP%InitPower%VectorPowerAmp_odd(1) = mag_amplitude_hel(helical_ind, helical_amp)*psconst(11)
                !write(*,*) "helical mag_amplitude_hel = ", mag_amplitude_hel(helical_ind, helical_amp) !useless,check
                write(*,*) "helical CP%InitPower%CorrType_hel = ", CP%InitPower%CorrType_hel !useless,check
                
             end if                
          end if 
          
       else
          write (*,*) "There are no passive vector modes."
          return
       end if
      
    end if


    ! A little debug output
    if(DebugMsgs .and. Feedbacklevel > 0) then
       write(*,*) "\t=== Mag debug ==="
       write(*,*) "\tMag amp: ", magamp, "Pert amp: ", CP%InitPower%ScalarPowerAmp(1)
       if(CP%WantTensors) then
          write(*,*) "\tMag ind: ", magind, "CAMB ind: ", CP%InitPower%ant(1)
       else
          write(*,*) "\tMag ind: ", magind, "CAMB ind: ", CP%InitPower%an(1)
       end if
       write(*,*) "\tMag prd: ", maglrat
    end if

    ! Fetch results into temporary array. Generally we'll pull the
    ! results out of the CAMB Cl_scalar etc. variables
    call cls_from_params(CP,Cltemp,ud)

    ! Set mag perturbations back to zero.
    delb = 0._dl
    pib = 0._dl

    ! Reset the Alfven velocity related stuff
    call mag_reset_alfven

  end subroutine mag_cls

end module magnetic

