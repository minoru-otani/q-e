!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE init_1drism(rismt, alpha, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... initialize 1D-RISM
  !
  USE err_rism, ONLY : merge_ierr_rism, IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,    ONLY : DP
  USE rism,     ONLY : rism_type, ITYPE_1DRISM
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  REAL(DP),        INTENT(IN)    :: alpha
  INTEGER,         INTENT(OUT)   :: ierr
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_1DRISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... create potentials
  IF (rismt%is_intra) THEN
    CALL potential_vv(rismt, ierr)
  ELSE
    ierr = IERR_RISM_NULL
  END IF
  !
  CALL merge_ierr_rism(ierr, rismt%super_comm)
  IF (ierr /= IERR_RISM_NULL) THEN
    RETURN
  END IF
  !
  ! ... create intra-molecular correlation
  IF (rismt%is_intra) THEN
    CALL molecorr_vv(rismt, alpha, ierr)
  ELSE
    ierr = IERR_RISM_NULL
  END IF
  !
  CALL merge_ierr_rism(ierr, rismt%super_comm)
  IF (ierr /= IERR_RISM_NULL) THEN
    RETURN
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE init_1drism
