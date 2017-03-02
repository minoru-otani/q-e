!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE solute
  !--------------------------------------------------------------------------
  !
  ! ... this module keeps data of solute.
  !
  USE constants,      ONLY : BOHR_RADIUS_ANGS
  USE err_rism,       ONLY : stop_by_err_rism, IERR_RISM_NULL, &
                           & IERR_RISM_LJ_UNSUPPORTED, IERR_RISM_LJ_OUT_OF_RANGE
  USE ions_base,      ONLY : nat, atm, ityp
  USE kinds,          ONLY : DP
  USE molecule_const, ONLY : RY_TO_KCALMOLm1
  USE rism,           ONLY : rism_type
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... define constants
  INTEGER, PARAMETER :: LEN_NAME = 10
  !
  ! ... maximum range of Lennard-Jones potentials
  REAL(DP)                             :: rmax_lj = 0.0_DP  ! (in sigma)
  !
  ! ... number of solute's atoms in super cell
  INTEGER                              :: solU_nat
  !
  ! ... coordinate of solute's atoms in super cell
  REAL(DP),                ALLOCATABLE :: solU_tau(:,:)  ! (in alat)
  !
  ! ... Lennard-Jones parameters of solute's atoms in unit cell
  REAL(DP),                ALLOCATABLE :: solU_ljeps(:)  ! (in Ry)
  REAL(DP),                ALLOCATABLE :: solU_ljsig(:)  ! (in bohr)
  CHARACTER(LEN=LEN_NAME), ALLOCATABLE :: solU_ljname(:)
  !
  ! ... index of atom, super cell -> unit cell
  INTEGER,                 ALLOCATABLE :: isup_to_iuni(:)
  !
  ! ... public components
  PUBLIC :: rmax_lj
  PUBLIC :: solU_nat
  PUBLIC :: solU_tau
  PUBLIC :: solU_ljeps
  PUBLIC :: solU_ljsig
  PUBLIC :: solU_ljname
  PUBLIC :: isup_to_iuni
  PUBLIC :: allocate_solU
  PUBLIC :: deallocate_solU
  PUBLIC :: update_solU
  PUBLIC :: get_solU_LJ_force
  PUBLIC :: get_solU_LJ_stress
  PUBLIC :: set_solU_LJ_param
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE allocate_solU()
    !--------------------------------------------------------------------------
    !
    ! ... initialize this module
    !
    IMPLICIT NONE
    !
    solU_nat = 0
    !
    ALLOCATE(solU_ljeps(nat))
    ALLOCATE(solU_ljsig(nat))
    ALLOCATE(solU_ljname(nat))
    !
  END SUBROUTINE allocate_solU
  !
  !--------------------------------------------------------------------------
  SUBROUTINE deallocate_solU()
    !--------------------------------------------------------------------------
    !
    ! ... finalize this module
    !
    IMPLICIT NONE
    !
    solU_nat = 0
    !
    IF (ALLOCATED(solU_tau))     DEALLOCATE(solU_tau)
    IF (ALLOCATED(solU_ljeps))   DEALLOCATE(solU_ljeps)
    IF (ALLOCATED(solU_ljsig))   DEALLOCATE(solU_ljsig)
    IF (ALLOCATED(solU_ljname))  DEALLOCATE(solU_ljname)
    IF (ALLOCATED(isup_to_iuni)) DEALLOCATE(isup_to_iuni)
    !
  END SUBROUTINE deallocate_solU
  !
  !--------------------------------------------------------------------------
  SUBROUTINE update_solU(rismt, ierr)
    !--------------------------------------------------------------------------
    !
    ! ... update solute's structure
    !
    IMPLICIT NONE
    !
    TYPE(rism_type), INTENT(INOUT) :: rismt
    INTEGER,         INTENT(OUT)   :: ierr
    !
    ! ... update solU_tau and isup_to_iuni
    IF (ALLOCATED(solU_tau))     DEALLOCATE(solU_tau)
    IF (ALLOCATED(isup_to_iuni)) DEALLOCATE(isup_to_iuni)
    !
    CALL lj_setup_solU_tau(rismt, rmax_lj, .TRUE.,  ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      RETURN
    END IF
    !
    ALLOCATE(solU_tau(3,  solU_nat))
    ALLOCATE(isup_to_iuni(solU_nat))
    !
    CALL lj_setup_solU_tau(rismt, rmax_lj, .FALSE., ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      RETURN
    END IF
    !
    ! ... update Lennard-Jones potentials
    CALL lj_setup_solU_vlj(rismt, rmax_lj, ierr)
    !
  END SUBROUTINE update_solU
  !
  !--------------------------------------------------------------------------
  SUBROUTINE get_solU_LJ_force(rismt, force, ierr)
    !--------------------------------------------------------------------------
    !
    ! ... add solute's force (Lennard-Jones)
    !
    IMPLICIT NONE
    !
    TYPE(rism_type), INTENT(IN)  :: rismt
    REAL(DP),        INTENT(OUT) :: force(3, nat)
    INTEGER,         INTENT(OUT) :: ierr
    !
    ! ... calculate Lennard-Jones force
    CALL lj_get_force(rismt, force, rmax_lj, ierr)
    !
  END SUBROUTINE get_solU_LJ_force
  !
  !--------------------------------------------------------------------------
  SUBROUTINE get_solU_LJ_stress(rismt, sigma, ierr)
    !--------------------------------------------------------------------------
    !
    ! ... add solute's stress (Lennard-Jones)
    !
    IMPLICIT NONE
    !
    TYPE(rism_type), INTENT(IN)  :: rismt
    REAL(DP),        INTENT(OUT) :: sigma(3, 3)
    INTEGER,         INTENT(OUT) :: ierr
    !
    ! ... calculate Lennard-Jones stress
    CALL lj_get_stress(rismt, sigma, rmax_lj, ierr)
    !
  END SUBROUTINE get_solU_LJ_stress
  !
  !--------------------------------------------------------------------------
  SUBROUTINE set_solU_LJ_param(is, ffname, eps, sig)
    !--------------------------------------------------------------------------
    !
    ! ... set Lennard-Jones parameters
    !
    IMPLICIT NONE
    !
    INTEGER,          INTENT(IN) :: is
    CHARACTER(LEN=*), INTENT(IN) :: ffname
    REAL(DP),         INTENT(IN) :: eps
    REAL(DP),         INTENT(IN) :: sig
    !
    INTEGER                 :: i
    INTEGER                 :: ia
    INTEGER                 :: atomn
    CHARACTER(LEN=32)       :: ffname_
    REAL(DP)                :: eps_
    REAL(DP)                :: sig_
    CHARACTER(LEN=LEN_NAME) :: name_
    INTEGER                 :: ierr
    !
    CHARACTER(LEN=1), EXTERNAL :: capital
    INTEGER,          EXTERNAL :: atomic_number
    !
    ! ... to upper case
    ffname_ = ADJUSTL(ffname)
    DO i = 1, LEN_TRIM(ffname_)
      ffname_(i:i) = capital(ffname_(i:i))
    END DO
    !
    DO ia = 1, nat
      IF (ityp(ia) /= is) THEN
        CYCLE
      END IF
      !
      ! ... ffname -> eps_ (kcal/mol)
      ! ...        -> sig_ (angstrom)
      SELECT CASE (TRIM(ffname_))
      CASE ('NONE')
        eps_  = 0.0_DP
        sig_  = 0.0_DP
        name_ = '???'
        ierr  = IERR_RISM_NULL
        !
      CASE ('UFF')
        atomn = atomic_number(trim(atm(is)))
        CALL lj_uff(atomn, eps_, sig_, ierr)
        name_ = 'UFF'
        !
      CASE ('CLAYFF')
        atomn = atomic_number(trim(atm(is)))
        CALL lj_clayff(atomn, eps_, sig_, ierr)
        name_ = 'ClayFF'
        !
      CASE ('OPLS-AA')
        atomn = atomic_number(trim(atm(is)))
        CALL lj_oplsaa(atomn, eps_, sig_, ierr)
        name_ = 'OPLS-AA'
        !
      CASE DEFAULT
        eps_  = 0.0_DP
        sig_  = 0.0_DP
        name_ = '???'
        ierr  = IERR_RISM_LJ_UNSUPPORTED
        CALL infomsg('set_solU_LJ_param', 'incorrect force field name: ' // TRIM(ADJUSTL(ffname)))
        !
      END SELECT
      !
      ! ... eps -> eps_ (kcal/mol)
      ! ... sig -> sig_ (angstrom)
      IF (ierr == IERR_RISM_NULL) THEN
        IF (eps > 0.0_DP) THEN
          eps_  = eps
          name_ = 'given'
        END IF
        !
        IF (sig > 0.0_DP) THEN
          sig_  = sig
          name_ = 'given'
        END IF
        !
        IF (eps_ <= 0.0_DP .OR. sig_ <= 0.0_DP) THEN
          ierr = IERR_RISM_LJ_OUT_OF_RANGE
        END IF
      END IF
      !
      ! ... check status
      IF (ierr /= IERR_RISM_NULL) THEN
        CALL stop_by_err_rism('set_solU_LJ_param', ierr)
      END IF
      !
      ! ... eps_ -> solU_ljeps (Ry)
      ! ... sig_ -> solU_ljsig (bohr)
      solU_ljeps( ia) = eps_ / RY_TO_KCALMOLm1
      solU_ljsig( ia) = sig_ / BOHR_RADIUS_ANGS
      solU_ljname(ia) = name_
      !
    END DO
    !
  END SUBROUTINE set_solU_LJ_param
  !
END MODULE solute
