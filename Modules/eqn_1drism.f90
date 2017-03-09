!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE eqn_1drism(rismt, gmax, lhand, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... solve 1D-RISM equation, which is defined as
  ! ...   [1 - w * c * rho] * h = w * c * w  (in G-space)
  ! ...
  ! ... (F.Hirata et al., Chem. Phys. Lett. 1981, 83, 329-334)
  !
  USE constants, ONLY : K_BOLTZMANN_RY
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE, &
                      & IERR_RISM_CANNOT_DGETRF, IERR_RISM_CANNOT_DGETRS, &
                      & merge_ierr_rism
  USE kinds,     ONLY : DP
  USE rism,      ONLY : rism_type, ITYPE_1DRISM
  USE solvmol,   ONLY : solVs, get_nsite_in_solVs, isite_to_isolV
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  REAL(DP),        INTENT(IN)    :: gmax   ! maximum of g-vectors
  LOGICAL,         INTENT(IN)    :: lhand  ! if true, right-hand. if false, left-hand
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER               :: ig, jg
  INTEGER               :: iig
  INTEGER               :: nv
  INTEGER               :: iv1, iv2
  INTEGER               :: ivv
  INTEGER               :: isolV
  INTEGER               :: ilapack
  INTEGER,  ALLOCATABLE :: ipiv(:)
  REAL(DP)              :: beta
  REAL(DP), ALLOCATABLE :: rho(:)
  REAL(DP), ALLOCATABLE :: hvv(:,:)
  REAL(DP), ALLOCATABLE :: cvv(:,:)
  REAL(DP), ALLOCATABLE :: wvv(:,:)
  REAL(DP), ALLOCATABLE :: avv(:,:)
  REAL(DP), ALLOCATABLE :: bvv(:,:)
  !
  EXTERNAL :: dgemm
  EXTERNAL :: dgetrf
  EXTERNAL :: dgetrs
  !
  ! ... number of sites in solvents
  nv = get_nsite_in_solVs()
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_1DRISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nr /= rismt%ng) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nsite < (nv * (nv + 1) / 2)) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... beta = 1 / (kB * T)
  beta = 1.0_DP / K_BOLTZMANN_RY / rismt%temp
  !
  ! ... allocate working memory
  ALLOCATE(ipiv(nv))
  ALLOCATE(rho(nv))
  ALLOCATE(hvv(nv, nv))
  ALLOCATE(cvv(nv, nv))
  ALLOCATE(wvv(nv, nv))
  ALLOCATE(avv(nv, nv))
  ALLOCATE(bvv(nv, nv))
  !
  ! ... make rho
  DO iv1 = 1, nv
    isolV = isite_to_isolV(iv1)
    IF (lhand) THEN
      rho(iv1) = solVs(isolV)%density
    ELSE
      rho(iv1) = solVs(isolV)%subdensity
    END IF
  END DO
  !
  ! ... in case G = 0
  IF (rismt%mp_task%ivec_start == 1) THEN
    jg = 2
    rismt%hg(1, :) = 0.0_DP
  ELSE
    jg = 1
  END IF
  !
  ! ... in case G /= 0
  ! ... 1D-RISM equation for each ig
  DO ig = jg, rismt%ng
    !
    ! ... check gmax
    IF (gmax > 0.0_DP) THEN
      iig = rismt%mp_task%ivec_start + ig - 1
      IF (rismt%rfft%ggrid(iig) > gmax) THEN
        hvv = 0.0_DP
        GOTO 100
      END IF
    END IF
    !
    ! ... extract data at ig
    DO iv1 = 1, nv
      DO iv2 = 1, iv1
        ivv = iv1 * (iv1 - 1) / 2 + iv2
        cvv(iv2, iv1) = rismt%csg(ig, ivv) - beta * rismt%ulg(ig, ivv)
        cvv(iv1, iv2) = cvv(iv2, iv1)
        wvv(iv2, iv1) = rismt%wg(ig, ivv)
        wvv(iv1, iv2) = wvv(iv2, iv1)
      END DO
    END DO
    !
    ! ... make tvv, avv and bvv
    ! ... b -> w * c
    CALL dgemm('T', 'N', nv, nv, nv, 1.0_DP, wvv, nv, cvv, nv, 0.0_DP, bvv, nv)
    ! ... a -> 1 - b * rho
    DO iv1 = 1, nv
      avv(:, iv1) = -bvv(:, iv1) * rho(iv1)
    END DO
    DO iv1 = 1, nv
      avv(iv1, iv1) = avv(iv1, iv1) + 1.0_DP
    END DO
    ! ... h -> b * w
    CALL dgemm('N', 'N', nv, nv, nv, 1.0_DP, bvv, nv, wvv, nv, 0.0_DP, hvv, nv)
    !
    ! ... solve linear equation
    CALL dgetrf(nv, nv, avv, nv, ipiv, ilapack)
    IF (ilapack /= 0) THEN
      ierr = IERR_RISM_CANNOT_DGETRF
      GOTO 200
    END IF
    CALL dgetrs('N', nv, nv, avv, nv, ipiv, hvv, nv, ilapack)
    IF (ilapack /= 0) THEN
      ierr = IERR_RISM_CANNOT_DGETRS
      GOTO 200
    END IF
    !
    ! ... set hvv to rismt
100 CONTINUE
    DO iv1 = 1, nv
      DO iv2 = 1, iv1
        ivv = iv1 * (iv1 - 1) / 2 + iv2
        rismt%hg(ig, ivv) = hvv(iv2, iv1)
      END DO
    END DO
    !
  END DO
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
  ! ... deallocate working memory
200 CONTINUE
  DEALLOCATE(ipiv)
  DEALLOCATE(rho)
  DEALLOCATE(hvv)
  DEALLOCATE(cvv)
  DEALLOCATE(wvv)
  DEALLOCATE(avv)
  DEALLOCATE(bvv)
  !
  ! ... merge error code through all processies
  CALL merge_ierr_rism(ierr, rismt%mp_site%inter_sitg_comm)
  !
END SUBROUTINE eqn_1drism
