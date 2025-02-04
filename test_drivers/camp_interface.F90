!> \file
!> The mam4_camp_interface module.

!> An interface between MAM4 and the CAMP
module mam4_camp_interface

  use shr_kind_mod, only: r8 => shr_kind_r8
#ifdef MAM4_USE_CAMP
  use mam4_state
  use camp_camp_core
  use camp_camp_state
  use camp_chem_spec_data
  use camp_aero_rep_data
  use camp_aero_rep_modal_binned_mass
  use camp_constants
  use camp_util
  use json_module
  use camp_solver_stats
#endif

  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef MAM4_USE_CAMP
  !> Run the CAMP module for the current MAM4 state
  subroutine mam4_camp_interface_solve( &
          env_state, aero_state, gas_state, del_t)

    !> CAMP core
    type(camp_core_t), pointer :: camp_core
    !> CAMP state
    type(camp_state_t), pointer :: camp_state
    !> Environment.
    type(env_state_t), intent(inout) :: env_state
    !> Aerosol state
    type(aero_state_t), intent(inout) :: aero_state
    !> Gas state
    type(gas_state_t), intent(inout) :: gas_state
    !> Time step (s)
    real(kind=r8), intent(in) :: del_t
    integer :: n
    character(len=4), parameter :: aero_rep_key = "MAM4"
    character(len=16), parameter :: mode_names(4) = (/ "accumulation    ", &
                                                       "aitken          ", &
                                                       "coarse          ", &
                                                       "primary carbon  " /)
    integer(kind=i_kind) :: mode(4)
    class(aero_rep_data_t), pointer :: aero_rep_ptr
    type(aero_rep_update_data_modal_binned_mass_GMD_t) :: update_data_GMD
    type(aero_rep_update_data_modal_binned_mass_GSD_t) :: update_data_GSD

    type(solver_stats_t), target :: solver_stats

    camp_core => camp_core_t("/home/dquevedo/AMBRS/ambrs_with_mam4-camp_support/tests/config.json")
    call camp_core%initialize()

    !call camp_core%solver_initialize()
    !camp_state => camp_core%new_state()

    call assert(209301925, camp_core%get_aero_rep(aero_rep_key, aero_rep_ptr))
    select type (aero_rep_ptr)
          type is (aero_rep_modal_binned_mass_t)
            call camp_core%initialize_update_object(aero_rep_ptr, &
                                                     update_data_GMD)
            call camp_core%initialize_update_object(aero_rep_ptr, &
                                                     update_data_GSD)
            call camp_core%solver_initialize()
            camp_state => camp_core%new_state()
            do n = 1, 4
                call assert_msg(431201141, &
                      aero_rep_ptr%get_section_id(mode_names(n), mode(n)), &
                      "Could not get mode ID")
                !print *, mode_names(n), 'GMD =', aero_state%GMD(n), 'GSD =', aero_state%GSD(n)
                call update_data_GMD%set_GMD(mode(n), aero_state%GMD(n))
                call update_data_GSD%set_GSD(mode(n), aero_state%GSD(n))
                call camp_core%update_data(update_data_GMD)
                call camp_core%update_data(update_data_GSD)
            end do
        class default
            write(*,*) "wrong aerosol rep"
            stop 3
    end select

    !call camp_core%solver_initialize()
    !camp_state => camp_core%new_state()

    ! Set the CAMP environmental state.
    call env_state_set_camp_env_state(env_state, camp_state)

    ! Set the CAMP gas-phase species concentrations
    call gas_state_set_camp_conc(camp_core, gas_state, env_state, camp_state)

    
    ! Update the mass concentrations and composition for all particles
    call mam4_camp_interface_set_camp_aerosol(aero_state, &
                                              camp_core, camp_state)

    ! Solve the multi-phase chemical system
    call camp_core%solve(camp_state, del_t, solver_stats = solver_stats)
    !call solver_stats%print()

    call mam4_camp_interface_get_camp_aerosol(aero_state, &
                                          camp_core, camp_state)

    ! Update the MAM4 gas-phase state
    call gas_state_get_camp_conc(gas_state, camp_state, camp_core)

    deallocate( camp_core, camp_state )

  end subroutine mam4_camp_interface_solve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Set the CAMP aerosol-phase species and mass concentrations
  subroutine mam4_camp_interface_set_camp_aerosol(aero_state, &
      camp_core, camp_state)

    !> Aerosol state.
    type(aero_state_t), intent(in) :: aero_state
    !> CAMP core.
    type(camp_core_t), pointer, intent(inout) :: camp_core
    !> CAMP state.
    type(camp_state_t), pointer, intent(out) :: camp_state

    integer(kind=i_kind) :: idx_H2O_1, &
                            idx_H2O_2, &
                            idx_H2O_3, &
                            idx_SO4_1, &
                            idx_SO4_2, &
                            idx_SO4_3, &
                            idx_POM_1, &
                            idx_POM_2, &
                            idx_POM_3, &
                            idx_POM_4, &
                            idx_SOA_1, &
                            idx_SOA_2, &
                            idx_SOA_3, &
                            idx_SOA_4, &
                            idx_BC_1, &
                            idx_BC_3, &
                            idx_BC_4, &
                            idx_DST_1, &
                            idx_DST_3, &
                            idx_NCL_1, &
                            idx_NCL_2, &
                            idx_NCL_3
    type(chem_spec_data_t), pointer :: chem_spec_data
    class(aero_rep_data_t), pointer :: aero_rep_data
    
    if( .not.camp_core%get_aero_rep( "MAM4", aero_rep_data ) ) then
        write(*,*) "Something's gone wrong!"
        stop 3
    end if
    idx_SO4_1 = aero_rep_data%spec_state_id("accumulation.aqueous sulfate.SO4")
    idx_SO4_2 = aero_rep_data%spec_state_id("aitken.aqueous sulfate.SO4")
    idx_SO4_3 = aero_rep_data%spec_state_id("coarse.aqueous sulfate.SO4")
    idx_POM_1 = aero_rep_data%spec_state_id("accumulation.organic matter.POM")
    idx_POM_2 = aero_rep_data%spec_state_id("aitken.organic matter.POM")
    idx_POM_3 = aero_rep_data%spec_state_id("coarse.organic matter.POM")
    idx_POM_4 = aero_rep_data%spec_state_id("primary carbon.organic matter.POM")
    idx_SOA_1 = aero_rep_data%spec_state_id("accumulation.organic matter.SOA")
    idx_SOA_2 = aero_rep_data%spec_state_id("aitken.organic matter.SOA")
    idx_SOA_3 = aero_rep_data%spec_state_id("coarse.organic matter.SOA")
    idx_SOA_4 = aero_rep_data%spec_state_id("primary carbon.organic matter.SOA")
    idx_BC_1 = aero_rep_data%spec_state_id("accumulation.soot.BC")
    idx_BC_3 = aero_rep_data%spec_state_id("coarse.soot.BC")
    idx_BC_4 = aero_rep_data%spec_state_id("primary carbon.soot.BC")
    idx_DST_1 = aero_rep_data%spec_state_id("accumulation.dust.DST")
    idx_DST_3 = aero_rep_data%spec_state_id("coarse.dust.DST")
    idx_NCL_1 = aero_rep_data%spec_state_id("accumulation.other PM.NCL")
    idx_NCL_2 = aero_rep_data%spec_state_id("aitken.other PM.NCL")
    idx_NCL_3 = aero_rep_data%spec_state_id("coarse.other PM.NCL")
    
    !idx_H2O_1 = aero_rep_data%spec_state_id("accumulation.aqueous sulfate.H2O_aq")
    !idx_H2O_2 = aero_rep_data%spec_state_id("aitken.aqueous sulfate.H2O_aq")
    !idx_H2O_3 = aero_rep_data%spec_state_id("coarse.aqueous sulfate.H2O_aq")

    camp_state%state_var(idx_SO4_1) = aero_state%qso4(1)
    camp_state%state_var(idx_SO4_2) = aero_state%qso4(2)
    camp_state%state_var(idx_SO4_3) = aero_state%qso4(3)
    camp_state%state_var(idx_POM_1) = aero_state%qpom(1)
    camp_state%state_var(idx_POM_2) = aero_state%qpom(2)
    camp_state%state_var(idx_POM_3) = aero_state%qpom(3)
    camp_state%state_var(idx_POM_4) = aero_state%qpom(4)
    camp_state%state_var(idx_SOA_1) = aero_state%qsoa(1)
    camp_state%state_var(idx_SOA_2) = aero_state%qsoa(2)
    camp_state%state_var(idx_SOA_3) = aero_state%qsoa(3)
    camp_state%state_var(idx_SOA_4) = aero_state%qsoa(4)
    camp_state%state_var(idx_BC_1) = aero_state%qbc(1)
    camp_state%state_var(idx_BC_3) = aero_state%qbc(3)
    camp_state%state_var(idx_BC_4) = aero_state%qbc(4)
    camp_state%state_var(idx_DST_1) = aero_state%qdst(1)
    camp_state%state_var(idx_DST_3) = aero_state%qdst(3)
    camp_state%state_var(idx_NCL_1) = aero_state%qncl(1)
    camp_state%state_var(idx_NCL_2) = aero_state%qncl(2)
    camp_state%state_var(idx_NCL_3) = aero_state%qncl(3)

    !camp_state%state_var(idx_H2O_1) = 1.0e-3_r8
    !camp_state%state_var(idx_H2O_2) = 1.0e-3_r8
    !camp_state%state_var(idx_H2O_3) = 1.0e0-3_r8

    print *, 'Pre-solve SO4', camp_state%state_var(idx_SO4_1) + &
                              camp_state%state_var(idx_SO4_2) + &
                              camp_state%state_var(idx_SO4_3)

  end subroutine mam4_camp_interface_set_camp_aerosol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get the CAMP aerosol-phase species and mass concentrations
  subroutine mam4_camp_interface_get_camp_aerosol(aero_state, &
      camp_core, camp_state)
      
    !> Aerosol state.
    type(aero_state_t), intent(out) :: aero_state
    !> CAMP core.
    type(camp_core_t), pointer, intent(inout) :: camp_core
    !> CAMP state.
    type(camp_state_t), pointer, intent(in) :: camp_state

    integer(kind=i_kind) :: idx_H2O_1, &
                            idx_H2O_2, &
                            idx_H2O_3, &
                            idx_SO4_1, &
                            idx_SO4_2, &
                            idx_SO4_3, &
                            idx_POM_1, &
                            idx_POM_2, &
                            idx_POM_3, &
                            idx_POM_4, &
                            idx_SOA_1, &
                            idx_SOA_2, &
                            idx_SOA_3, &
                            idx_SOA_4, &
                            idx_BC_1, &
                            idx_BC_3, &
                            idx_BC_4, &
                            idx_DST_1, &
                            idx_DST_3, &
                            idx_NCL_1, &
                            idx_NCL_2, &
                            idx_NCL_3
    type(chem_spec_data_t), pointer :: chem_spec_data
    class(aero_rep_data_t), pointer :: aero_rep_data

    if( .not.camp_core%get_aero_rep( "MAM4", aero_rep_data ) ) then
        write(*,*) "Something's gone wrong!"
        stop 3
    end if
    
    idx_SO4_1 = aero_rep_data%spec_state_id("accumulation.aqueous sulfate.SO4")
    idx_SO4_2 = aero_rep_data%spec_state_id("aitken.aqueous sulfate.SO4")
    idx_SO4_3 = aero_rep_data%spec_state_id("coarse.aqueous sulfate.SO4")
    idx_POM_1 = aero_rep_data%spec_state_id("accumulation.organic matter.POM")
    idx_POM_2 = aero_rep_data%spec_state_id("aitken.organic matter.POM")
    idx_POM_3 = aero_rep_data%spec_state_id("coarse.organic matter.POM")
    idx_POM_4 = aero_rep_data%spec_state_id("primary carbon.organic matter.POM")
    idx_SOA_1 = aero_rep_data%spec_state_id("accumulation.organic matter.SOA")
    idx_SOA_2 = aero_rep_data%spec_state_id("aitken.organic matter.SOA")
    idx_SOA_3 = aero_rep_data%spec_state_id("coarse.organic matter.SOA")
    idx_SOA_4 = aero_rep_data%spec_state_id("primary carbon.organic matter.SOA")
    idx_BC_1 = aero_rep_data%spec_state_id("accumulation.soot.BC")
    idx_BC_3 = aero_rep_data%spec_state_id("coarse.soot.BC")
    idx_BC_4 = aero_rep_data%spec_state_id("primary carbon.soot.BC")
    idx_DST_1 = aero_rep_data%spec_state_id("accumulation.dust.DST")
    idx_DST_3 = aero_rep_data%spec_state_id("coarse.dust.DST")
    idx_NCL_1 = aero_rep_data%spec_state_id("accumulation.other PM.NCL")
    idx_NCL_2 = aero_rep_data%spec_state_id("aitken.other PM.NCL")
    idx_NCL_3 = aero_rep_data%spec_state_id("coarse.other PM.NCL")

    aero_state%qso4(1) = camp_state%state_var(idx_SO4_1)
    aero_state%qso4(2) = camp_state%state_var(idx_SO4_2)
    aero_state%qso4(3) = camp_state%state_var(idx_SO4_3)
    aero_state%qpom(1) = camp_state%state_var(idx_POM_1)
    aero_state%qpom(2) = camp_state%state_var(idx_POM_2)
    aero_state%qpom(3) = camp_state%state_var(idx_POM_3)
    aero_state%qpom(4) = camp_state%state_var(idx_POM_4)
    aero_state%qsoa(1) = camp_state%state_var(idx_SOA_1)
    aero_state%qsoa(2) = camp_state%state_var(idx_SOA_2)
    aero_state%qsoa(3) = camp_state%state_var(idx_SOA_3)
    aero_state%qsoa(4) = camp_state%state_var(idx_SOA_4)
    aero_state%qbc(1) = camp_state%state_var(idx_BC_1)
    aero_state%qbc(3) = camp_state%state_var(idx_BC_3)
    aero_state%qbc(4) = camp_state%state_var(idx_BC_4)
    aero_state%qdst(1) = camp_state%state_var(idx_DST_1)
    aero_state%qdst(3) = camp_state%state_var(idx_DST_3)
    aero_state%qncl(1) = camp_state%state_var(idx_NCL_1)
    aero_state%qncl(2) = camp_state%state_var(idx_NCL_2)
    aero_state%qncl(3) = camp_state%state_var(idx_NCL_3)

    print *, 'Post-solve SO4', camp_state%state_var(idx_SO4_1) + &
                               camp_state%state_var(idx_SO4_2) + &
                               camp_state%state_var(idx_SO4_3)
    print *, 'Post-solve aero_state SO4', aero_state%qso4(1) + &
                                          aero_state%qso4(2) + &
                                          aero_state%qso4(3)

  end subroutine mam4_camp_interface_get_camp_aerosol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
end module mam4_camp_interface
