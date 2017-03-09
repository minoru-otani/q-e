!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE err_rism
  !--------------------------------------------------------------------------
  !
  USE mp, ONLY : mp_size, mp_rank, mp_gather, mp_bcast
  !
  ! ... this module defines error code used in RISM.
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! not error
  INTEGER, PARAMETER :: IERR_RISM_NULL                = 0
  ! incorrect data type
  INTEGER, PARAMETER :: IERR_RISM_INCORRECT_DATA_TYPE = 1
  ! 1D-RISM data is not available
  INTEGER, PARAMETER :: IERR_RISM_1DRISM_IS_NOT_AVAIL = 2
  ! 1D-, 3D-RISM iteration is not converged
  INTEGER, PARAMETER :: IERR_RISM_NOT_CONVERGED       = 3
  ! Lennard-Jones parameters are unsupported
  INTEGER, PARAMETER :: IERR_RISM_LJ_UNSUPPORTED      = 4
  ! Lennard-Jones parameters are out of range
  INTEGER, PARAMETER :: IERR_RISM_LJ_OUT_OF_RANGE     = 5
  ! cannot solve dgetrf (in LAPACK)
  INTEGER, PARAMETER :: IERR_RISM_CANNOT_DGETRF       = 6
  ! cannot solve dgetrs (in LAPACK)
  INTEGER, PARAMETER :: IERR_RISM_CANNOT_DGETRS       = 7
  ! charge of solvent is not zero (in Laue-RISM)
  INTEGER, PARAMETER :: IERR_RISM_NONZERO_CHARGE      = 8
  ! fail to make potential smooth (in Laue-RISM)
  INTEGER, PARAMETER :: IERR_RISM_FAIL_SMOOTH         = 9
  ! too large box of Laue-FFT (in Laue-RISM)
  INTEGER, PARAMETER :: IERR_RISM_LARGE_LAUE_BOX      = 10
  !
  ! ... public components
  PUBLIC :: stop_by_err_rism
  PUBLIC :: merge_ierr_rism
  PUBLIC :: IERR_RISM_NULL
  PUBLIC :: IERR_RISM_INCORRECT_DATA_TYPE
  PUBLIC :: IERR_RISM_1DRISM_IS_NOT_AVAIL
  PUBLIC :: IERR_RISM_NOT_CONVERGED
  PUBLIC :: IERR_RISM_LJ_UNSUPPORTED
  PUBLIC :: IERR_RISM_LJ_OUT_OF_RANGE
  PUBLIC :: IERR_RISM_CANNOT_DGETRF
  PUBLIC :: IERR_RISM_CANNOT_DGETRS
  PUBLIC :: IERR_RISM_NONZERO_CHARGE
  PUBLIC :: IERR_RISM_FAIL_SMOOTH
  PUBLIC :: IERR_RISM_LARGE_LAUE_BOX
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE stop_by_err_rism(parent, ierr)
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: parent
    INTEGER,          INTENT(IN) :: ierr
    !
    SELECT CASE (ierr)
    CASE (IERR_RISM_INCORRECT_DATA_TYPE)
      CALL errore(' ' // TRIM(ADJUSTL(parent)) // ' ', &
         & ' in RISM, incorrect data type ', MAX(ABS(ierr), 1))
      !
    CASE (IERR_RISM_1DRISM_IS_NOT_AVAIL)
      CALL errore(' ' // TRIM(ADJUSTL(parent)) // ' ', &
         & ' in RISM, data of 1D is not available ', MAX(ABS(ierr), 1))
      !
    CASE (IERR_RISM_NOT_CONVERGED)
      CALL errore(' ' // TRIM(ADJUSTL(parent)) // ' ', &
         & ' in RISM, iteration has not been converged ', MAX(ABS(ierr), 1))
      !
    CASE (IERR_RISM_LJ_UNSUPPORTED)
      CALL errore(' ' // TRIM(ADJUSTL(parent)) // ' ', &
         & ' in RISM, specified L.J.-parameters are not supported ', MAX(ABS(ierr), 1))
      !
    CASE (IERR_RISM_LJ_OUT_OF_RANGE)
      CALL errore(' ' // TRIM(ADJUSTL(parent)) // ' ', &
         & ' in RISM, specified L.J.-parameters are out of range ', MAX(ABS(ierr), 1))
      !
    CASE (IERR_RISM_CANNOT_DGETRF)
      CALL errore(' ' // TRIM(ADJUSTL(parent)) // ' ', &
         & ' in RISM, error at lapack::dgetrf ', MAX(ABS(ierr), 1))
      !
    CASE (IERR_RISM_CANNOT_DGETRS)
      CALL errore(' ' // TRIM(ADJUSTL(parent)) // ' ', &
         & ' in RISM, error at lapack::dgetrs ', MAX(ABS(ierr), 1))
      !
    CASE (IERR_RISM_NONZERO_CHARGE)
      CALL errore(' ' // TRIM(ADJUSTL(parent)) // ' ', &
         & ' in RISM, charge of solvent is not zero ', MAX(ABS(ierr), 1))
      !
    CASE (IERR_RISM_FAIL_SMOOTH)
      CALL errore(' ' // TRIM(ADJUSTL(parent)) // ' ', &
         & ' in RISM, fail to make potential smooth ', MAX(ABS(ierr), 1))
      !
    CASE (IERR_RISM_LARGE_LAUE_BOX)
      CALL errore(' ' // TRIM(ADJUSTL(parent)) // ' ', &
         & ' in RISM, too large box of Laue-FFT than 1D-FFT ', MAX(ABS(ierr), 1))
      !
    END SELECT
    !
  END SUBROUTINE stop_by_err_rism
  !
  !--------------------------------------------------------------------------
  SUBROUTINE merge_ierr_rism(ierr, comm)
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(INOUT) :: ierr
    INTEGER, INTENT(IN)    :: comm
    !
    INTEGER              :: nproc
    INTEGER              :: iproc
    INTEGER              :: irank
    INTEGER, ALLOCATABLE :: iallerr(:)
    !
    nproc = mp_size(comm)
    irank = mp_rank(comm)
    ALLOCATE(iallerr(nproc))
    !
    CALL mp_gather(ierr, iallerr, 0, comm)
    !
    IF (irank == 0) THEN
      ierr = IERR_RISM_NULL
      DO iproc = 1, nproc
        IF (iallerr(iproc) /= IERR_RISM_NULL) THEN
          ierr = iallerr(iproc)
          EXIT
        END IF
      END DO
    END IF
    !
    CALL mp_bcast(ierr, 0, comm)
    !
    DEALLOCATE(iallerr)
    !
  END SUBROUTINE merge_ierr_rism
  !
END MODULE err_rism

