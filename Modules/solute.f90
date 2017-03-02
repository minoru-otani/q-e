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
  USE cell_base,      ONLY : at, bg, alat
  USE constants,      ONLY : BOHR_RADIUS_ANGS
  USE err_rism,       ONLY : stop_by_err_rism, IERR_RISM_NULL, &
                           & IERR_RISM_LJ_UNSUPPORTED, IERR_RISM_LJ_OUT_OF_RANGE
  USE ions_base,      ONLY : nat, atm, nsp, ityp, tau
  USE kinds,          ONLY : DP
  USE molecule_const, ONLY : RY_TO_KCALMOLm1
  USE rism,           ONLY : rism_type
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... define constants
  INTEGER, PARAMETER :: LEN_NAME = 12
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
    INTEGER                 :: iis
    INTEGER                 :: coord
    CHARACTER(LEN=32)       :: ffname_
    REAL(DP)                :: eps_
    REAL(DP)                :: sig_
    CHARACTER(LEN=LEN_NAME) :: name
    CHARACTER(LEN=5)        :: label
    INTEGER                 :: ierr
    LOGICAL,  ALLOCATABLE   :: loxy(:)
    REAL(DP), ALLOCATABLE   :: tau_uni(:,:)
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
    IF (TRIM(ffname_) == 'CLAYFF') THEN
      ! ... index of oxygen
      ALLOCATE(loxy(nsp))
      DO iis = 1, nsp
        atomn = atomic_number(TRIM(atm(iis)))
        loxy(iis) = (atomn == 8)
      END DO
      !
      ! ... pack atoms in unit cell
      ALLOCATE(tau_uni(3, nat))
      tau_uni = tau
      CALL cryst_to_cart(nat, tau_uni, bg, -1)
      tau_uni = tau_uni - FLOOR(tau_uni)
      CALL cryst_to_cart(nat, tau_uni, at, +1)
      !
    END IF
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
        eps_ = 0.0_DP
        sig_ = 0.0_DP
        name = '???'
        ierr = IERR_RISM_NULL
        !
      CASE ('UFF')
        atomn = atomic_number(TRIM(atm(is)))
        CALL lj_uff(atomn, eps_, sig_, ierr)
        name = 'UFF'
        !
      CASE ('CLAYFF')
        atomn = atomic_number(TRIM(atm(is)))
        coord = coordination_number(ia, atomn, loxy, tau_uni)
        CALL lj_clayff(atomn, coord, eps_, sig_, label, ierr)
        name = 'ClayFF' // TRIM(label)
        !
      CASE ('OPLS-AA')
        atomn = atomic_number(TRIM(atm(is)))
        CALL lj_oplsaa(atomn, eps_, sig_, ierr)
        name = 'OPLS-AA'
        !
      CASE DEFAULT
        eps_ = 0.0_DP
        sig_ = 0.0_DP
        name = '???'
        ierr = IERR_RISM_LJ_UNSUPPORTED
        CALL infomsg('set_solU_LJ_param', 'incorrect force field name: ' // TRIM(ADJUSTL(ffname)))
        !
      END SELECT
      !
      ! ... eps -> eps_ (kcal/mol)
      ! ... sig -> sig_ (angstrom)
      IF (ierr == IERR_RISM_NULL) THEN
        IF (eps > 0.0_DP) THEN
          eps_ = eps
          name = 'given'
        END IF
        !
        IF (sig > 0.0_DP) THEN
          sig_ = sig
          name = 'given'
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
      solU_ljname(ia) = name
      !
    END DO
    !
    ! ... deallocate memory
    IF (ALLOCATED(loxy)) THEN
      DEALLOCATE(loxy)
    END IF
    IF (ALLOCATED(tau_uni)) THEN
      DEALLOCATE(tau_uni)
    END IF
    !
  END SUBROUTINE set_solU_LJ_param
  !
  !--------------------------------------------------------------------------
  FUNCTION coordination_number(ia, atomn, loxy, tuni) RESULT(cn)
    !--------------------------------------------------------------------------
    !
    ! ... calculate coordination number of oxygen atoms
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(IN) :: ia
    INTEGER,  INTENT(IN) :: atomn
    LOGICAL,  INTENT(IN) :: loxy(nsp)
    REAL(DP), INTENT(IN) :: tuni(3, nat)
    INTEGER              :: cn
    !
    INTEGER  :: ja
    INTEGER  :: js
    INTEGER  :: i1, i2, i3
    REAL(DP) :: rad_m
    REAL(DP) :: rad_mo
    REAL(DP) :: rrad_mo
    REAL(DP) :: rr
    REAL(DP) :: t1, t2, t3
    REAL(DP) :: x1, y1, z1
    REAL(DP) :: x2, y2, z2
    REAL(DP) :: x3, y3, z3
    REAL(DP) :: dx, dy, dz
    !
    INTEGER,  PARAMETER :: ICELL  = 1
    REAL(DP), PARAMETER :: RSCALE = 1.2_DP
    !
    ! ... crystal ionic radii (in angstrom)
    REAL(DP), PARAMETER :: RAD_O  = 1.26_DP  !(-2)
    REAL(DP), PARAMETER :: RAD_SI = 0.54_DP  !(+4)
    REAL(DP), PARAMETER :: RAD_AL = 0.675_DP !(+3)
    REAL(DP), PARAMETER :: RAD_MG = 0.86_DP  !(+2)
    REAL(DP), PARAMETER :: RAD_CA = 1.14_DP  !(+2)
    REAL(DP), PARAMETER :: RAD_FE = 0.92_DP  !(+2)
    REAL(DP), PARAMETER :: RAD_LI = 0.90_DP  !(+1)
    !
    cn = 0
    !
    SELECT CASE(atomn)
    CASE(14) ! Si
      rad_m = RAD_SI
    CASE(13) ! Al
      rad_m = RAD_AL
    CASE(12) ! Mg
      rad_m = RAD_MG
    CASE(20) ! Ca
      rad_m = RAD_CA
    CASE(26) ! Fe
      rad_m = RAD_FE
    CASE( 3) ! Li
      rad_m = RAD_LI
    CASE DEFAULT
      rad_m = -1.0_DP
    END SELECT
    !
    IF (rad_m <= 0.0_DP) THEN
      RETURN
    END IF
    !
    rad_mo  = RSCALE * (rad_m + RAD_O)
    rad_mo  = rad_mo / BOHR_RADIUS_ANGS
    rad_mo  = rad_mo / alat
    rrad_mo = rad_mo * rad_mo
    !
    x1 = tuni(1, ia)
    y1 = tuni(2, ia)
    z1 = tuni(3, ia)
    !
    DO ja = 1, nat
      js = ityp(ja)
      IF (js < 1 .OR. nsp < js .OR. (.NOT. loxy(js))) THEN
        CYCLE
      END IF
      !
      x2 = tuni(1, ja)
      y2 = tuni(2, ja)
      z2 = tuni(3, ja)
      !
      DO i1 = -ICELL, ICELL
        t1 = DBLE(i1)
        DO i2 = -ICELL, ICELL
          t2 = DBLE(i2)
          DO i3 = -ICELL, ICELL
            t3 = DBLE(i3)
            !
            x3 = x2 + t1 * at(1, 1) + t2 * at(1, 2) + t3 * at(1, 3)
            y3 = y2 + t1 * at(2, 1) + t2 * at(2, 2) + t3 * at(2, 3)
            z3 = z2 + t1 * at(3, 1) + t2 * at(3, 2) + t3 * at(3, 3)
            !
            dx = x1 - x3
            dy = y1 - y3
            dz = z1 - z3
            rr = dx * dx + dy * dy + dz * dz
            IF (rr < rrad_mo) THEN
              cn = cn + 1
            END IF
            !
          END DO
        END DO
      END DO
      !
    END DO
    !
  END FUNCTION coordination_number
  !
END MODULE solute
