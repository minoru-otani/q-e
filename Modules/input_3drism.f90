!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE iosys_3drism(laue, linit)
  !-----------------------------------------------------------------------------
  !
  ! ...  Copy data read from input file (in subroutine "read_input_file") and
  ! ...  stored in modules input_parameters into internal modules
  ! ...  Note: this subroutine requires nsp(ions_base), ecutrho(gvect), dual(gvecs), alat(cell_base).
  !
  USE cell_base,        ONLY : alat
  USE gvecs,            ONLY : dual
  USE gvect,            ONLY : ecutrho
  USE ions_base,        ONLY : nsp
  USE kinds,            ONLY : DP
  USE rism,             ONLY : CLOSURE_HNC, CLOSURE_KH
  USE rism3d_facade,    ONLY : starting_corr, niter, epsv, starting_epsv, mdiis_size, mdiis_step, &
                             & ecutsolv_ => ecutsolv, rism3t, rism3d_initialize, &
                             & conv_always, planar_average, laue_nfit_ => laue_nfit, &
                             & expand_r, expand_l, starting_r, starting_l, buffer_r, buffer_l, &
                             & both_hands, ireference, IREFERENCE_NULL, IREFERENCE_AVERAGE, &
                             & IREFERENCE_RIGHT, IREFERENCE_LEFT
  USE solute,           ONLY : rmax_lj_ => rmax_lj, allocate_solU, set_solU_LJ_param
  USE solvmol,          ONLY : get_nuniq_in_solVs
  !
  ! ... CONTROL namelist
  !
  USE input_parameters, ONLY : restart_mode
  !
  ! ... RISM namelist
  !
  USE input_parameters, ONLY : closure, starting3d, tempv, ecutsolv, smear3d, &
                               solute_lj, solute_epsilon, solute_sigma, rmax_lj, &
                               rism3d_maxstep, rism3d_conv_thr, mdiis3d_size, mdiis3d_step, &
                               rism3d_conv_always, rism3d_planar_average, &
                               laue_nfit, laue_expand_right, laue_expand_left, &
                               laue_starting_right, laue_starting_left, &
                               laue_buffer_right, laue_buffer_left, &
                               laue_both_hands, laue_reference
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: laue
  LOGICAL, INTENT(IN) :: linit
  !
  INTEGER :: is
  INTEGER :: nsite
  !
  INTEGER,  PARAMETER :: MDIIS_SWITCH    = 4
  REAL(DP), PARAMETER :: MDIIS_STEP_DEF1 = 0.8_DP
  REAL(DP), PARAMETER :: MDIIS_STEP_DEF2 = 0.3_DP
  !
  ! ... check starting condition.
  IF (TRIM(restart_mode) == 'restart') THEN
    IF (TRIM(starting3d) /= 'file') THEN
      CALL infomsg('input','WARNING: "starting3d" set to '//TRIM(starting3d)//' may spoil restart')
      starting3d = 'file'
    END IF
  END IF
  !
  ! ... modify ecutsolv
  IF (ecutsolv <= 0.0_DP) THEN
    ecutsolv = ecutrho * 4.0_DP / dual
  ELSE
    ecutsolv = MAX(ecutsolv, ecutrho * 4.0_DP / dual)
    ecutsolv = MIN(ecutsolv, ecutrho)
  ENDIF
  !
  ! ... modify rism3d_planar_average
  IF (laue) THEN
    rism3d_planar_average = .TRUE.
  END IF
  !
  ! ... modify laue_reference
  IF (laue .AND. TRIM(laue_reference) == 'none') THEN
    IF (laue_expand_right > 0.0_DP .AND. laue_expand_left > 0.0_DP) THEN
      laue_reference = 'average'
    ELSE IF (laue_expand_right > 0.0_DP) THEN
      laue_reference = 'right'
    ELSE IF (laue_expand_left > 0.0_DP) THEN
      laue_reference = 'left'
    END IF
  END IF
  !
  ! ... modify mdiis3d_step
  IF (mdiis3d_step < 0.0_DP) THEN
    nsite = get_nuniq_in_solVs()
    IF (nsite <= MDIIS_SWITCH) THEN
      mdiis3d_step = MDIIS_STEP_DEF1
    ELSE
      mdiis3d_step = MDIIS_STEP_DEF2
    END IF
  END IF
  !
  ! ... set from namelist. these data are already checked.
  starting_corr  = starting3d
  niter          = rism3d_maxstep
  epsv           = rism3d_conv_thr
  starting_epsv  = epsv
  mdiis_size     = mdiis3d_size
  mdiis_step     = mdiis3d_step
  ecutsolv_      = ecutsolv
  conv_always    = rism3d_conv_always
  planar_average = rism3d_planar_average
  laue_nfit_     = laue_nfit
  expand_r       = laue_expand_right
  expand_l       = laue_expand_left
  starting_r     = laue_starting_right
  starting_l     = laue_starting_left
  buffer_r       = laue_buffer_right
  buffer_l       = laue_buffer_left
  both_hands     = laue_both_hands
  rmax_lj_       = rmax_lj
  !
  IF (TRIM(laue_reference) == 'none') THEN
    ireference = IREFERENCE_NULL
  ELSE IF (TRIM(laue_reference) == 'average') THEN
    ireference = IREFERENCE_AVERAGE
  ELSE IF (TRIM(laue_reference) == 'right') THEN
    ireference = IREFERENCE_RIGHT
  ELSE IF (TRIM(laue_reference) == 'left') THEN
    ireference = IREFERENCE_LEFT
  END IF
  !
  ! ... convert units
  IF (expand_r > 0.0_DP) THEN
    expand_r = expand_r / alat
  END IF
  !
  IF (expand_l > 0.0_DP) THEN
    expand_l = expand_l / alat
  END IF
  !
  IF (starting_r /= 0.0_DP) THEN
    starting_r = starting_r / alat
  END IF
  !
  IF (starting_l /= 0.0_DP) THEN
    starting_l = starting_l / alat
  END IF
  !
  IF (buffer_r /= 0.0_DP) THEN
    buffer_r = buffer_r / alat
  END IF
  !
  IF (buffer_l /= 0.0_DP) THEN
    buffer_l = buffer_l / alat
  END IF
  !
  ! ... initialize solute
  CALL allocate_solU()
  DO is = 1, nsp
    CALL set_solU_LJ_param(is, solute_lj(is), solute_epsilon(is), solute_sigma(is))
  END DO
  !
  ! ... initialize rism3d_facade
  IF (TRIM(closure) == 'hnc') THEN
    rism3t%closure = CLOSURE_HNC
  ELSE IF (TRIM(closure) == 'kh') THEN
    rism3t%closure = CLOSURE_KH
  END IF
  rism3t%temp = tempv
  rism3t%tau  = smear3d
  !
  IF (linit) THEN
    CALL rism3d_initialize(laue)
  END IF
  !
END SUBROUTINE iosys_3drism
