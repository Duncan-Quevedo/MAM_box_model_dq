!> \file
!> The mam4_state module.

!> The state_t structures and associated subroutines.

module mam4_state
#ifdef MAM4_USE_CAMP
      use camp_camp_core
      use camp_camp_state
      use camp_constants
      use camp_chem_spec_data
      use shr_kind_mod, only: r8 => shr_kind_r8

      implicit none

      type env_state_t
        !> Temperature (K).
        real(kind=r8) :: temp
        !> Ambient pressure (Pa).
        real(kind=r8) :: press
        !> Relative humidity (0-1).
        real(kind=r8) :: RH_CLEA
      end type env_state_t
      type aero_state_t
        !> Mass fractions (0-1).
        real(kind=r8) :: qso4(4)
        real(kind=r8) :: qpom(4)
        real(kind=r8) :: qsoa(4)
        real(kind=r8) :: qbc(4)
        real(kind=r8) :: qdst(4)
        real(kind=r8) :: qncl(4)
        real(kind=r8) :: GMD(4), GSD(4)
      end type aero_state_t
      type gas_state_t
        !> Gas mixing ratios.
        real(kind=r8), allocatable :: vmr(:)
      end type gas_state_t
      integer :: idx_so2g, idx_h2so4g
      real(kind=r8) :: to_kgperm3

      contains

              subroutine env_state_set_camp_env_state(env_state, &
                                                      camp_state)

                type(env_state_t), intent(in) :: env_state
                type(camp_state_t), intent(inout) :: camp_state

                camp_state%env_states(1)%val%rel_humid = &
                        env_state%RH_CLEA
                call camp_state%env_states(1)%set_temperature_K(env_state%temp)
                call camp_state%env_states(1)%set_pressure_Pa(env_state%press)

              end subroutine env_state_set_camp_env_state

              subroutine gas_state_set_camp_conc(camp_core, &
                                                 gas_state, &
                                                 env_state, &
                                                 camp_state)
                
                type(gas_state_t), intent(in) :: gas_state
                type(env_state_t), intent(in) :: env_state
                type(camp_state_t), intent(inout) :: camp_state
                type(camp_core_t), pointer :: camp_core
                integer, parameter :: i_camp_water = 0

                integer(kind=i_kind) :: idx_SO2, &
                                        idx_H2SO4

                type(chem_spec_data_t), pointer :: chem_spec_data
            
                if( .not.camp_core%get_chem_spec_data( chem_spec_data ) ) then
                    write(*,*) "Something's gone wrong!"
                    stop 3
                end if

                idx_SO2 = chem_spec_data%gas_state_id("SO2")
                idx_H2SO4 = chem_spec_data%gas_state_id("H2SO4")

                camp_state%state_var(i_camp_water) = &
                        env_state_rh_to_y(env_state)

                camp_state%state_var(idx_SO2) =  gas_state%vmr(idx_so2g)
                camp_state%state_var(idx_H2SO4) = gas_state%vmr(idx_h2so4g)

                print *, 'MAM4 H2SO4', gas_state%vmr(idx_h2so4g)
                print *, 'CAMP H2SO4', camp_state%state_var(idx_H2SO4)
              
              end subroutine gas_state_set_camp_conc

              real(kind=r8) function  env_state_rh_to_y(env_state)

                type(env_state_t), intent(in) :: env_state

                real(kind=r8), parameter :: t_steam = 373.15d0
                real(kind=r8) :: a, water_vp

                a = 1.0d0 - t_steam / env_state%temp
                a = (((-0.1299 * a - 0.6445) * a - 1.976) * a + 13.3185) * a
                water_vp = 101325.0 * exp(a)  ! (Pa)
                env_state_rh_to_y = env_state%RH_CLEA * water_vp * 1.0e6 &
                        / env_state%press ! (ppm)

              end function env_state_rh_to_y

              subroutine gas_state_get_camp_conc(gas_state, &
                                                 camp_state, camp_core)
                
                type(gas_state_t), intent(inout) :: gas_state
                type(camp_state_t), intent(inout) :: camp_state
                type(camp_core_t), pointer :: camp_core

                integer(kind=i_kind) :: idx_SO2, &
                                        idx_H2SO4

                type(chem_spec_data_t), pointer :: chem_spec_data
            
                if( .not.camp_core%get_chem_spec_data( chem_spec_data ) ) then
                    write(*,*) "Something's gone wrong!"
                    stop 3
                end if

                idx_SO2 = chem_spec_data%gas_state_id("SO2")
                idx_H2SO4 = chem_spec_data%gas_state_id("H2SO4")

                gas_state%vmr(idx_so2g) = camp_state%state_var(idx_SO2)
                gas_state%vmr(idx_h2so4g) = camp_state%state_var(idx_H2SO4)

                print *, 'MAM4 H2SO4', gas_state%vmr(idx_h2so4g)
                print *, 'CAMP H2SO4', camp_state%state_var(idx_H2SO4)
              
              end subroutine gas_state_get_camp_conc
#endif
end module mam4_state

