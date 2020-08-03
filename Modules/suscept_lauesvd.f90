!
! Copyright (C) 2018 National Institute of Advanced Industrial Science and Technology (AIST)
! [ This code is written by Satomichi Nishihara. ]
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE suscept_lauesvd(rismt, lhand, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... SVD of inter-site susceptibility for Laue-RISM
  !
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum
  USE rism,      ONLY : rism_type, ITYPE_LAUERISM
  USE solvmol,   ONLY : get_nuniq_in_solVs
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  LOGICAL,         INTENT(IN)    :: lhand  ! if true, right-hand. if false, left-hand.
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER               :: nq
  INTEGER               :: iq1
  INTEGER               :: iq2
  INTEGER               :: iiq2
  INTEGER               :: irz1
  INTEGER               :: irz2
  INTEGER               :: jrz
  INTEGER               :: info
  INTEGER               :: lwork
  REAL(DP)              :: x
  REAL(DP), ALLOCATABLE :: xmajor(:,:)
  REAL(DP), ALLOCATABLE :: x0(:)
  REAL(DP), ALLOCATABLE :: xu(:)
  REAL(DP), ALLOCATABLE :: xvt(:,:)
  REAL(DP), ALLOCATABLE :: work(:)
  !
  EXTERNAL :: dgesvd
  !
  ! ... number of sites in solvents
  nq = get_nuniq_in_solVs()
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%mp_site%nsite < nq) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nrzl < rismt%lfft%nrz) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... calculation only for Gxy = 0
  IF (rismt%lfft%gxystart > 1) THEN
    !NOP
  ELSE
    GOTO 1
  END IF
  !
  ! ... allocate memory
  ALLOCATE(xmajor(nq * rismt%lfft%nrz * 3, nq * rismt%lfft%nrz))
  ALLOCATE(x0    (nq * rismt%lfft%nrz))
  ALLOCATE(xu    (nq * rismt%lfft%nrz * 3))
  ALLOCATE(xvt   (nq * rismt%lfft%nrz, nq * rismt%lfft%nrz))
  xmajor = 0.0_DP
  x0     = 0.0_DP
  xu     = 0.0_DP
  xvt    = 0.0_DP
  !
  lwork = 8 * nq * rismt%lfft%nrz
  ALLOCATE(work(lwork))
  work = 0.0_DP
  !
  ! ... construct major matrix of susceptibility
  xmajor = 0.0_DP
  !
  DO iq1 = 1, nq
    !
    DO iq2 = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      iiq2 = iq2 - rismt%mp_site%isite_start + 1
      !
      DO irz1 = 1, rismt%lfft%nrz
        !
        DO irz2 = 1, rismt%lfft%nrz * 3
          !
          jrz = ABS(irz1 - irz2) + 1
          IF (jrz <= rismt%lfft%nrz) THEN
            IF (lhand) THEN
              x = rismt%xgs(jrz, iiq2, iq1)
            ELSE
              x = rismt%ygs(jrz, iiq2, iq1)
            END IF
          ELSE
            x = 0.0_DP
          END IF
          !
          xmajor((iq2 - 1) * rismt%lfft%nrz * 3 + irz2, &
               & (iq1 - 1) * rismt%lfft%nrz     + irz1) = x
          !
        END DO
        !
      END DO
      !
    END DO
    !
  END DO
  !
  IF (nq * rismt%lfft%nrz > 0) THEN
    CALL mp_sum(xmajor, rismt%mp_site%intra_sitg_comm)
  END IF
  !
  ! ... SVD of susceptibility
  CALL dgesvd('N', 'A', nq * rismt%lfft%nrz * 3, nq * rismt%lfft%nrz, &
            & xmajor, nq * rismt%lfft%nrz * 3, x0, xu, 1, &
            & xvt, nq * rismt%lfft%nrz, work, lwork, info)
  !
  ! ... deallocate memory
  DEALLOCATE(xmajor)
  DEALLOCATE(x0)
  DEALLOCATE(xu)
  DEALLOCATE(xvt)
  DEALLOCATE(work)
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE suscept_lauesvd
