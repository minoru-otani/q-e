!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE iosys_1drism(laue)
  !-----------------------------------------------------------------------------
  !
  ! ...  Copy data read from input file (in subroutine "read_input_file") and
  ! ...  stored in modules input_parameters into internal modules
  ! ...  Note: this subroutine requires pseudo_dir(io_files), omega(cell_base),
  ! ...        ecutrho(gvect), dual(gvecs).
  !
  USE cell_base,        ONLY : alat, at, omega
  USE constants,        ONLY : tpi
  USE gvecs,            ONLY : dual
  USE gvect,            ONLY : ecutrho
  USE io_files,         ONLY : molfile
  USE kinds,            ONLY : DP
  USE molecule_const,   ONLY : BOHRm3_TO_MOLCMm3, BOHRm3_TO_MOLLm1
  USE read_solv_module, ONLY : read_solvents
  USE rism,             ONLY : CLOSURE_HNC, CLOSURE_KH
  USE rism1d_facade,    ONLY : nproc_sub, nproc_switch, starting_corr, niter, epsv, bond_width, &
                             & mdiis_size, mdiis_step, rism1t, rism1d_initialize, &
                             & rism1d_activate_right, rism1d_activate_left
  USE solvmol,          ONLY : nsolV_ => nsolV, solVs, get_nsite_in_solVs
  !
  ! ... CONTROL namelist
  !
  USE input_parameters, ONLY : restart_mode
  !
  ! ... RISM namelist
  !
  USE input_parameters, ONLY : nsolv, closure, starting1d, tempv, permittivity, rmax1d, &
                               smear1d, rism1d_maxstep, rism1d_conv_thr, rism1d_bond_width, &
                               rism1d_nproc, rism1d_nproc_switch, mdiis1d_size, mdiis1d_step, &
                               laue_expand_right, laue_expand_left, laue_both_hands
  !
  ! ... SOLVENTS card
  !
  USE input_parameters, ONLY : solv_label, solv_mfile, solv_dens1, solv_dens2, solvents_unit
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: laue
  !
  INTEGER                        :: ngrid
  INTEGER                        :: nsite
  INTEGER                        :: isolV
  REAL(DP)                       :: z0
  REAL(DP)                       :: zright
  REAL(DP)                       :: zleft
  CHARACTER(LEN=10), ALLOCATABLE :: slabel(:)
  REAL(DP),          ALLOCATABLE :: sdens1(:)
  REAL(DP),          ALLOCATABLE :: sdens2(:)
  REAL(DP)                       :: pertot, per1
  REAL(DP)                       :: dentot, den1
  INTEGER                        :: ihand
  INTEGER                        :: nhand
  !
  REAL(DP),          PARAMETER   :: RMAX1D_SCALE    = 1.5_DP
  REAL(DP),          PARAMETER   :: BOND_SCALE      = 4.0_DP
  INTEGER,           PARAMETER   :: MDIIS_SWITCH    = 6
  REAL(DP),          PARAMETER   :: MDIIS_STEP_DEF1 = 0.5_DP
  REAL(DP),          PARAMETER   :: MDIIS_STEP_DEF2 = 0.1_DP
  !
  ! ... allocate memory
  ALLOCATE(slabel(nsolv))
  ALLOCATE(sdens1(nsolv))
  ALLOCATE(sdens2(nsolv))
  !
  ! ... check starting condition.
  IF (TRIM(restart_mode) == 'restart') THEN
    IF (TRIM(starting1d) /= 'file' .AND. TRIM(starting1d) /= 'fix') THEN
      CALL infomsg('input','WARNING: "starting1d" set to '//TRIM(starting1d)//' may spoil restart')
      starting1d = 'file'
    END IF
  END IF
  !
  ! ... modify rmax1d
  IF (laue) THEN
    z0 = 0.5_DP * alat * at(3, 3)
    zright =  z0 + MAX(0.0_DP, laue_expand_right)
    zleft  = -z0 - MAX(0.0_DP, laue_expand_left )
    rmax1d = MAX(rmax1d, RMAX1D_SCALE * (zright - zleft))
  ENDIF
  !
  ! ... modify rism1d_bond_width
  IF (laue .AND. rism1d_bond_width <= 0.0_DP) THEN
    rism1d_bond_width = BOND_SCALE / SQRT(ecutrho * 4.0_DP / dual)
  ENDIF
  !
  ! ... evaluate #grid
  ngrid = number_of_grids(rmax1d)
  !
  ! ... set from namelist. these data are already checked.
  nsolV_        = nsolv
  starting_corr = starting1d
  nproc_sub     = rism1d_nproc
  nproc_switch  = rism1d_nproc_switch
  niter         = rism1d_maxstep
  epsv          = rism1d_conv_thr
  bond_width    = rism1d_bond_width
  mdiis_size    = mdiis1d_size
  mdiis_step    = mdiis1d_step
  !
  ! ... set from card
  DO isolV = 1, nsolV_
    slabel( isolV) = solv_label(isolV)
    molfile(isolV) = solv_mfile(isolV)
    sdens1( isolV) = solv_dens1(isolV)
    sdens2( isolV) = solv_dens2(isolV)
  END DO
  !
  ! ... read solvents from molecule files
  CALL read_solvents()
  !
  ! ... set variables for solVs
  DO isolV = 1, nsolV_
    !
    ! ... name
    IF (LEN_TRIM(slabel(isolV)) > 0) THEN
      solVs(isolV)%name = slabel(isolv)
    ELSE
      CALL infomsg('iosys_1drism', &
      & 'default molecular name(formula) from MOL file('//TRIM(molfile(isolV))//') is used')
    END IF
    IF (LEN_TRIM(solVs(isolV)%name) <= 0) THEN
      CALL errore('iosys_1drism', 'invalid name', isolV)
    END IF
    !
    ! ... density
    IF (sdens1(isolV) > 0.0_DP) THEN
      CALL convert_dens(TRIM(ADJUSTL(solvents_unit)), isolV, sdens1(isolV))
      solVs(isolV)%density = sdens1(isolv)
    ELSE
      CALL infomsg('iosys_1drism', &
      & 'default density from MOL file('//TRIM(molfile(isolV))//') is used')
    END IF
    IF (solVs(isolV)%density <= 0.0_DP) THEN
      CALL errore('iosys_1drism', 'invalid density', isolV)
    END IF
    !
    ! ... subdensity
    IF (.NOT. laue_both_hands) THEN
      solVs(isolV)%subdensity = solVs(isolV)%density
      !
    ELSE
      IF (sdens2(isolV) > 0.0_DP) THEN
        CALL convert_dens(TRIM(ADJUSTL(solvents_unit)), isolV, sdens2(isolV))
        solVs(isolV)%subdensity = sdens2(isolv)
      ELSE
        CALL infomsg('iosys_1drism', &
        & 'default density from MOL file('//TRIM(molfile(isolV))//') is used')
      END IF
      IF (solVs(isolV)%subdensity <= 0.0_DP) THEN
        CALL errore('iosys_1drism', 'invalid density', isolV)
      END IF
    END IF
    !
  END DO
  !
  ! ... modify mdiis1d_step (this operation must be after read_solvents)
  IF (mdiis1d_step < 0.0_DP) THEN
    nsite = get_nsite_in_solVs()
    IF (nsite <= MDIIS_SWITCH) THEN
      mdiis1d_step = MDIIS_STEP_DEF1
    ELSE
      mdiis1d_step = MDIIS_STEP_DEF2
    END IF
    mdiis_step = mdiis1d_step
  END IF
  !
  ! ... modify permittivity (this operation must be after read_solvents)
  IF (permittivity <= 0.0_DP) THEN
    pertot = 0.0_DP
    dentot = 0.0_DP
    DO isolV = 1, nsolV_
      per1 = solVs(isolV)%permittivity
      den1 = 0.5_DP * (solVs(isolV)%density + solVs(isolV)%subdensity)
      IF (per1 > 0.0_DP)  THEN
        pertot = pertot + per1 * den1
        dentot = dentot + den1
      END IF
    END DO
    !
    IF (dentot > 0.0_DP) THEN
      permittivity = pertot / dentot
    ELSE
      permittivity = 0.0_DP
    END IF
  END IF
  !
  ! ... initialize rism1d_facade
  IF (laue .AND. laue_both_hands) THEN
    nhand = 2
  ELSE
    nhand = 1
  END IF
  !
  DO ihand = 1, nhand
    IF (ihand == 1) THEN
      CALL rism1d_activate_right()
    ELSE
      CALL rism1d_activate_left()
    END IF
    !
    IF (TRIM(closure) == 'hnc') THEN
      rism1t%closure = CLOSURE_HNC
    ELSE IF (TRIM(closure) == 'kh') THEN
      rism1t%closure = CLOSURE_KH
    END IF
    rism1t%temp = tempv
    rism1t%perm = permittivity
    rism1t%tau  = smear1d
  END DO
  !
  CALL rism1d_initialize(ngrid, rmax1d, nhand > 1)
  !
  ! ... deallocate memory
  DEALLOCATE(slabel)
  DEALLOCATE(sdens1)
  DEALLOCATE(sdens2)
  !
CONTAINS
  !
  SUBROUTINE convert_dens(dens_format, isolV, dens)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)    :: dens_format
    INTEGER,          INTENT(IN)    :: isolV
    REAL (DP),        INTENT(INOUT) :: dens
    !
    SELECT CASE (dens_format)
    CASE ('1/cell')
      dens = dens / omega
      !
    CASE ('mol/L')
      dens = dens / BOHRm3_TO_MOLLm1
      !
    CASE ('g/cm^3')
      dens = (dens / solVs(isolV)%mass) / BOHRm3_TO_MOLCMm3
      !
    CASE DEFAULT
      CALL errore('iosys_1drism', 'dens_format='// &
                & TRIM(dens_format)//' not implemented', isolV)
      !
    END SELECT
  END SUBROUTINE convert_dens
  !
  FUNCTION number_of_grids(rmax) RESULT(nr)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: rmax
    INTEGER              :: nr
    !
    REAL(DP) :: tpibr
    REAL(DP) :: tpibr2
    REAL(DP) :: ecutrism
    REAL(DP) :: gcutrism
    !
    REAL(DP), PARAMETER :: ECUT_SCALE = 1.1_DP
    !
    tpibr  = tpi / (2.0_DP * rmax)
    tpibr2 = tpibr * tpibr
    ecutrism = ECUT_SCALE * MAX(ecutrho, ecutrho * 4.0_DP / dual)
    !
    gcutrism = ecutrism / tpibr2
    nr = INT(SQRT(gcutrism)) + 1
  END FUNCTION number_of_grids
  !
END SUBROUTINE iosys_1drism
