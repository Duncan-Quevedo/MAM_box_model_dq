!> \file
!> The mam4_camp_interface module.

!> An interface between MAM4 and the CAMP
module mam4_camp_interface

  use mam4_state
  use shr_kind_mod, only: r8 => shr_kind_r8
#ifdef MAM4_USE_CAMP
  use camp_camp_core
  use camp_camp_state
  use camp_rxn_data
  use camp_solver_stats
  use camp_util, only: split_string
#endif

  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef MAM4_USE_CAMP
  !> Run the CAMP module for the current MAM4 state
  subroutine mam4_camp_interface_solve(camp_core, camp_state, &
          env_state, aero_data, aero_state, gas_state, del_t)

    !> CAMP core
    type(camp_core_t), intent(in) :: camp_core
    !> CAMP state
    type(camp_state_t), intent(inout) :: camp_state
    !> Working CAMP state
    !type(camp_state_t), intent(inout) :: camp_pre_aero_state
    !> Working CAMP state
    !type(camp_state_t), intent(inout) :: camp_post_aero_state
    !> Environment.
    type(env_state_t), intent(inout) :: env_state
    !> Aerosol data
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state
    type(aero_state_t), intent(inout) :: aero_state
    !> Gas data
    !type(gas_data_t), intent(in) :: gas_data
    !> Gas state
    type(gas_state_t), intent(inout) :: gas_state
    !> Photolysis calculator
!   type(photolysis_t), intent(inout) :: photolysis
    !> Time step (s)
    real(kind=r8), intent(in) :: del_t

    real(kind=r8) :: num_conc
    integer :: camp_state_size
    type(solver_stats_t), target :: solver_stats

    ! Set the CAMP environmental state.
    call env_state_set_camp_env_state(env_state, camp_state)

    ! Set the CAMP gas-phase species concentrations
    call gas_state_set_camp_conc(gas_state, env_state, camp_state)

    ! Recalculate the photolysis rate constants
    ! call photolysis%update_rate_constants()

    ! Update the number concentrations and composition for all particles
    call mam4_camp_interface_set_camp_aerosol(aero_data, aero_state, &
                                             camp_core, camp_state, &
                                             gas_state, env_state)

    ! We're modifying particle diameters, so the bin sort is now invalid
    ! aero_state%valid_sort = .false.

    ! Solve the multi-phase chemical system
    call camp_core%solve(camp_state, del_t, solver_stats = solver_stats)

    !call assert_msg(592148911, solver_stats%status_code == 0, &
    !     "Solver failed for aerosol-phase with code "// &
    !     integer_to_string(solver_stats%solver_flag))

    ! Update the PartMC gas-phase state
    call gas_state_get_camp_conc(gas_state, camp_state)

  end subroutine mam4_camp_interface_solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set the CAMP aerosol-phase species and number concentrations
  subroutine mam4_camp_interface_set_camp_aerosol(aero_data, aero_state, &
      camp_core, camp_state, gas_state, env_state)

    !> Aerosol data.
    class(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> CAMP core.
    type(camp_core_t), intent(in) :: camp_core
    !> CAMP state.
    type(camp_state_t), intent(inout) :: camp_state
    !> Gas state.
    type(gas_state_t), intent(in) :: gas_state
    !> Ambient state.
    type(env_state_t), intent(in) :: env_state

    real(kind=r8) :: numc, numc1, numc2, numc3, numc4
    real(kind=r8), parameter :: pi = 3.1415927

    numc1 = aero_data%numc1
    numc2 = aero_data%numc2
    numc3 = aero_data%numc3
    numc4 = aero_data%numc4
    numc =  numc1 + numc2 + numc3 + numc4

    !call camp_core%update_data(numc)

    camp_state%state_var((size(gas_state%vmr)+1):(size(gas_state%vmr)+6)) &
         = (pi/6) * (/ &
         aero_state%mfso41*numc1*aero_data%dgn(1)**3*exp(4.5*log(aero_data%sg(1))**2) &
         /env_state%adv_mass(6) + &
         aero_state%mfso42*numc2*aero_data%dgn(2)**3*exp(4.5*log(aero_data%sg(2))**2) &
         /env_state%adv_mass(14) + &
         aero_state%mfso43*numc3*aero_data%dgn(3)**3*exp(4.5*log(aero_data%sg(3))**2) &
         /env_state%adv_mass(21), &
         aero_state%mfpom1*numc1*aero_data%dgn(1)**3*exp(4.5*log(aero_data%sg(1))**2) &
         /env_state%adv_mass(7) + &
         aero_state%mfpom3*numc3*aero_data%dgn(3)**3*exp(4.5*log(aero_data%sg(3))**2) &
         /env_state%adv_mass(23) + & 
         aero_state%mfpom4*numc4*aero_data%dgn(4)**3*exp(4.5*log(aero_data%sg(4))**2) &
         /env_state%adv_mass(27), &
         aero_state%mfsoa1*numc1*aero_data%dgn(1)**3*exp(4.5*log(aero_data%sg(1))**2) &
         /env_state%adv_mass(8) + &
         aero_state%mfsoa2*numc2*aero_data%dgn(2)**3*exp(4.5*log(aero_data%sg(2))**2) &
         /env_state%adv_mass(15) + &
         aero_state%mfsoa3*numc3*aero_data%dgn(3)**3*exp(4.5*log(aero_data%sg(3))**2) &
         /env_state%adv_mass(24), &
         aero_state%mfbc1*numc1*aero_data%dgn(1)**3*exp(4.5*log(aero_data%sg(1))**2) &
         /env_state%adv_mass(9) + & 
         aero_state%mfbc4*numc4*aero_data%dgn(4)**3*exp(4.5*log(aero_data%sg(4))**2) &
         /env_state%adv_mass(28), &
         aero_state%mfdst1*numc1*aero_data%dgn(1)**3*exp(4.5*log(aero_data%sg(1))**2) &
         /env_state%adv_mass(10) + & 
         aero_state%mfdst3*numc3*aero_data%dgn(3)**3*exp(4.5*log(aero_data%sg(3))**2) &
         /env_state%adv_mass(19), &
         aero_state%mfncl1*numc1*aero_data%dgn(1)**3*exp(4.5*log(aero_data%sg(1))**2) &
         /env_state%adv_mass(11) + &
         aero_state%mfncl2*numc2*aero_data%dgn(2)**3*exp(4.5*log(aero_data%sg(2))**2) &
         /env_state%adv_mass(16) + &
         aero_state%mfncl3*numc3*aero_data%dgn(3)**3*exp(4.5*log(aero_data%sg(3))**2) &
         /env_state%adv_mass(20) /)

  end subroutine mam4_camp_interface_set_camp_aerosol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
end module mam4_camp_interface
