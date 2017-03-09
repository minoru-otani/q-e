!
! Copyright (C) 2016 National Institute of Advanced Industrial Science and Technology (AIST)
! [ This code is written by Satomichi Nishihara. ]
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE solvation_pbc(rismt, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... correct non-periodic solvation potential to periodic one,
  ! ... which acts for electronic wave functions.
  !
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE, IERR_RISM_FAIL_SMOOTH
  USE kinds,     ONLY : DP
  USE lauefft,   ONLY : fw_lauefft_1z
  USE rism,      ONLY : rism_type, ITYPE_LAUERISM
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER                  :: ig
  INTEGER                  :: iz
  INTEGER                  :: jz
  INTEGER                  :: igxy
  INTEGER                  :: jgxy
  COMPLEX(DP), ALLOCATABLE :: aux(:)
  COMPLEX(DP), ALLOCATABLE :: rhog_pbcl(:,:)
  COMPLEX(DP), ALLOCATABLE :: vpot_pbcl(:,:)
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nrzs < rismt%cfft%dfftt%nr3) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nrzl < rismt%lfft%nrz) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%ngxy < rismt%lfft%ngxy) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%pbc_nfit < 0) THEN
    ierr = IERR_RISM_FAIL_SMOOTH
    RETURN
  END IF
  !
  IF (rismt%cfft%dfftt%nr3 <= (2 * rismt%pbc_nfit)) THEN
    ierr = IERR_RISM_FAIL_SMOOTH
    RETURN
  END IF
  !
  ! ... allocate memory
  IF (rismt%cfft%dfftt%nnr > 0) THEN
    ALLOCATE(aux(rismt%cfft%dfftt%nnr))
  END IF
  IF (rismt%cfft%dfftt%nr3 * rismt%lfft%ngxy > 0) THEN
    ALLOCATE(rhog_pbcl(rismt%cfft%dfftt%nr3, rismt%lfft%ngxy))
    ALLOCATE(vpot_pbcl(rismt%cfft%dfftt%nr3, rismt%lfft%ngxy))
  END IF
  !
  ! ...
  ! ... fitting Rho
  ! ...
  ! ... copy rhog -> rhog_pbcl
  DO igxy = 1, rismt%lfft%ngxy
    jgxy = rismt%nrzl * (igxy - 1)
    DO iz = 1, rismt%cfft%dfftt%nr3
      jz = rismt%lfft%izcell_start + iz - 1
      rhog_pbcl(iz, igxy) = rismt%rhog(jz + jgxy)
    END DO
  END DO
  !
  ! ... correct around z=+-z0
  DO igxy = 1, rismt%lfft%ngxy
    IF (rismt%cfft%dfftt%nr3 > 0) THEN
      CALL make_smooth(rhog_pbcl(:, igxy), rismt%cfft%dfftt%nr3, rismt%pbc_nfit)
    END IF
  END DO
  !
  ! ... 1D-FFT: Laue-rep. -> G-space
  IF (rismt%cfft%dfftt%nnr > 0) THEN
    CALL fw_lauefft_1z(rismt%lfft, rhog_pbcl, rismt%cfft%dfftt%nr3, 1, aux)
  END IF
  !
  ! ... copy aux -> rhog_pbc
  IF (rismt%ng > 0) THEN
    rismt%rhog_pbc = CMPLX(0.0_DP, 0.0_DP, kind=DP)
  END IF
  DO ig = 1, rismt%cfft%ngmt
    rismt%rhog_pbc(ig) = aux(rismt%cfft%nlt(ig))
  END DO
  !
  ! ...
  ! ... fitting Vpot
  ! ...
  ! ... copy vpot -> vpot_pbcl
  DO igxy = 1, rismt%lfft%ngxy
    jgxy = rismt%nrzl * (igxy - 1)
    DO iz = 1, rismt%cfft%dfftt%nr3
      jz = rismt%lfft%izcell_start + iz - 1
      vpot_pbcl(iz, igxy) = rismt%vpot(jz + jgxy)
    END DO
  END DO
  !
  ! ... correct around z=+-z0
  DO igxy = 1, rismt%lfft%ngxy
    IF (rismt%cfft%dfftt%nr3 > 0) THEN
      CALL make_smooth(vpot_pbcl(:, igxy), rismt%cfft%dfftt%nr3, rismt%pbc_nfit)
    END IF
  END DO
  !
  ! ... 1D-FFT: Laue-rep. -> G-space
  IF (rismt%cfft%dfftt%nnr > 0) THEN
    CALL fw_lauefft_1z(rismt%lfft, vpot_pbcl, rismt%cfft%dfftt%nr3, 1, aux)
  END IF
  !
  ! ... copy aux -> vpot_pbc
  IF (rismt%ng > 0) THEN
    rismt%vpot_pbc = CMPLX(0.0_DP, 0.0_DP, kind=DP)
  END IF
  DO ig = 1, rismt%cfft%ngmt
    rismt%vpot_pbc(ig) = aux(rismt%cfft%nlt(ig))
  END DO
  !
  ! ... deallocate memory
  IF (rismt%cfft%dfftt%nnr > 0) THEN
    DEALLOCATE(aux)
  END IF
  IF (rismt%cfft%dfftt%nr3 * rismt%lfft%ngxy > 0) THEN
    DEALLOCATE(rhog_pbcl)
    DEALLOCATE(vpot_pbcl)
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
CONTAINS
  !
  SUBROUTINE make_smooth(vpot, npot, nfit)
    !
    IMPLICIT NONE
    !
    COMPLEX(kind=DP), INTENT(INOUT) :: vpot(:)
    INTEGER,          INTENT(IN)    :: npot
    INTEGER,          INTENT(IN)    :: nfit
    !
    INTEGER                       :: i
    INTEGER                       :: ndim
    REAL(kind=DP)                 :: rdim
    COMPLEX(kind=DP)              :: a, b, c, d
    COMPLEX(kind=DP)              :: value1
    COMPLEX(kind=DP)              :: value2
    COMPLEX(kind=DP)              :: slope1
    COMPLEX(kind=DP)              :: slope2
    COMPLEX(kind=DP), ALLOCATABLE :: vtmp(:)
    !
    IF (nfit < 2) THEN
      RETURN
    END IF
    !
    IF (npot < (2 * nfit + 1)) THEN
      RETURN
    END IF
    !
    ndim = 2 * nfit
    ALLOCATE(vtmp(0:(ndim + 1)))
    !
    ! ... vpot -> vtmp
    vtmp(       0) = vpot(npot - nfit)
    vtmp(ndim + 1) = vpot(nfit + 1)
    DO i = 1, nfit
      vtmp(       i) = vpot(npot - nfit + i)
      vtmp(nfit + i) = vpot(              i)
    END DO
    !
    ! ... fitting by cubic function.
    value1 = vtmp(1)
    value2 = vtmp(ndim)
    slope1 = vtmp(1) - vtmp(0)
    slope2 = vtmp(ndim + 1) - vtmp(ndim)
    !
    rdim = DBLE(ndim - 1)
    a =  (2.0_DP * (value1 - value2) + rdim * (slope1 + slope2)) / rdim / rdim / rdim
    b = -(3.0_DP * (value1 - value2) + rdim * (2.0_DP * slope1 + slope2)) / rdim / rdim
    c = slope1
    d = value1
    !
    DO i = 1, ndim
      rdim = DBLE(i - 1)
      vtmp(i) = a * (rdim * rdim * rdim) + b * (rdim * rdim) + c * rdim + d
    END DO
    !
    ! ... vtmp -> vpot
    DO i = 1, nfit
      vpot(npot - nfit + i) = vtmp(       i)
      vpot(              i) = vtmp(nfit + i)
    END DO
    !
    DEALLOCATE(vtmp)
    !
  END SUBROUTINE make_smooth
  !
END SUBROUTINE solvation_pbc
