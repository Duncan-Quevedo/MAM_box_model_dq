!> \file
!> The mam4_state module.

!> The state_t structures and associated subroutines.
module mam4_state

      use camp_camp_state
      use shr_kind_mod, only: r8 => shr_kind_r8

      type env_state_t
        !> Temperature (K).
        real(kind=r8) :: temp
        !> Ambient pressure (Pa).
        real(kind=r8) :: press
        !> Relative humidity (0-1).
        real(kind=r8) :: RH_CLEA
        !> Mass of air parcel
        real(kind=r8), allocatable :: adv_mass(:)
      end type env_state_t
      type aero_state_t
        !> Mass fractions (0-1) in Mode 1.
        real(kind=r8) :: mfso41
        real(kind=r8) :: mfpom1
        real(kind=r8) :: mfsoa1
        real(kind=r8) :: mfbc1
        real(kind=r8) :: mfdst1
        real(kind=r8) :: mfncl1
        !> Mass fractions (0-1) in Mode 2.
        real(kind=r8) :: mfso42
        real(kind=r8) :: mfsoa2
        real(kind=r8) :: mfncl2
        !> Mass fractions (0-1) in Mode 3.
        real(kind=r8) :: mfdst3
        real(kind=r8) :: mfncl3
        real(kind=r8) :: mfso43
        real(kind=r8) :: mfbc3
        real(kind=r8) :: mfpom3
        real(kind=r8) :: mfsoa3
        !> Mass fractions (0-1) in Mode 4.
        real(kind=r8) :: mfpom4
        real(kind=r8) :: mfbc4
      end type aero_state_t
      type gas_state_t
        !> Gas mixing ratios.
        real(kind=r8), allocatable :: vmr(:)
      end type gas_state_t
!     type gas_data_t
!       character(len=100), allocatable :: name(:)
!       integer, allocatable :: index
!       integer :: i_camp_water = 0
!     end type gas_data_t
      type aero_data_t
        real(kind=r8) :: numc1
        real(kind=r8) :: numc2
        real(kind=r8) :: numc3
        real(kind=r8) :: numc4
        real(kind=r8) :: dgn(4)
        real(kind=r8) :: sn(4)
      end type aero_data_t

      contains

              subroutine env_state_set_camp_env_state(env_state, &
                                                      camp_state)

                type(env_state_t), intent(in) :: env_state
                type(camp_state_t), intent(inout) :: camp_state

                camp_state%env_states(1)%val%rel_humid = &
                        env_state%RH_CLEA
                camp_state%env_states(1)%val%temp = env_state%temp

              end subroutine env_state_set_camp_env_state

              subroutine gas_state_set_camp_conc(gas_state, &
                                                 env_state, &
                                                 camp_state, &
                                                 !gas_data)
                
                type(gas_state_t), intent(in) :: gas_state
                type(env_state_t), intent(in) :: env_state
                type(camp_state_t), intent(inout) :: camp_state
                !class(gas_data_t), intent(in) :: gas_data
                integer, parameter :: i_camp_water = 0

                camp_state%state_var(i_camp_water) = &
                        env_state_rh_to_y(env_state)

                camp_state%state_var(1:size(gas_state%vmr)) = &
                        1.0d-3 * gas_state%vmr ! (ppm)
              
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
                                                 camp_state)
                
                type(gas_state_t), intent(in) :: gas_state
                type(camp_state_t), intent(inout) :: camp_state

                gas_state%vmr = 1.0d3 * & ! (ppb)
                            camp_state%state_var(1:size(gas_state%vmr))
              
              end subroutine gas_state_get_camp_conc

end module mam4_state
