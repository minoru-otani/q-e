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
SUBROUTINE solvation_lauerism(rismt, charge, ireference, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... setup charge density, solvation potential and solvation energy for DFT,
  ! ... which are derived from Laue-RISM.
  ! ...
  ! ... outside of the unit-cell, charge density is approximated as
  !
  !                   ----
  !        rho(gxy,z) >    q(v) * rho(v) * h(v; gxy,z)
  !                   ----
  !                     v
  !
  ! ... Variables:
  ! ...   charge:  total charge of solvent system
  !
  USE cell_base,      ONLY : at, alat
  USE constants,      ONLY : eps6
  USE err_rism,       ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE io_global,      ONLY : stdout
  USE kinds,          ONLY : DP
  USE lauefft,        ONLY : fw_lauefft_2xy
  USE mp,             ONLY : mp_sum
  USE rism,           ONLY : rism_type, ITYPE_LAUERISM
  USE solvmol,        ONLY : get_nuniq_in_solVs, solVs, iuniq_to_nsite, &
                           & iuniq_to_isite, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  REAL(DP),        INTENT(IN)    :: charge
  INTEGER,         INTENT(IN)    :: ireference
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER                  :: nq
  INTEGER                  :: iq
  INTEGER                  :: iiq
  INTEGER                  :: iv
  INTEGER                  :: nv
  INTEGER                  :: isolV
  INTEGER                  :: iatom
  INTEGER                  :: irz
  INTEGER                  :: iirz
  INTEGER                  :: nright
  INTEGER                  :: nleft
  INTEGER                  :: igxy
  INTEGER                  :: jgxy
  INTEGER                  :: kgxy
  INTEGER                  :: izright_tail
  INTEGER                  :: izleft_tail
  REAL(DP)                 :: rhov
  REAL(DP)                 :: rhov1
  REAL(DP)                 :: rhov2
  REAL(DP)                 :: qv
  REAL(DP)                 :: qtmp
  REAL(DP)                 :: dz
  REAL(DP)                 :: area_xy
  REAL(DP)                 :: fac
  REAL(DP)                 :: fac1
  REAL(DP)                 :: fac2
  REAL(DP)                 :: charge0
  REAL(DP)                 :: vsol0
  COMPLEX(DP), ALLOCATABLE :: ggz(:,:)
  !
  INTEGER,     PARAMETER   :: RHOG_NEDGE     = 3
  REAL(DP),    PARAMETER   :: RHOG_THRESHOLD = 1.0E-5_DP
  COMPLEX(DP), PARAMETER   :: C_ZERO         = CMPLX(0.0_DP, 0.0_DP, kind=DP)
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
  ! ... allocate memory
  IF (rismt%nrzs * rismt%ngxy * rismt%nsite > 0) THEN
    ALLOCATE(ggz(rismt%nrzs * rismt%ngxy, rismt%nsite))
  END IF
  !
  ! ... set variables
  dz      = rismt%lfft%zstep * alat
  area_xy = ABS(at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1)) * alat * alat
  !
  ! ... gr -> ggz
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq = iq - rismt%mp_site%isite_start + 1
    IF (rismt%nrzs * rismt%ngxy > 0) THEN
      ggz(:, iiq) = C_ZERO
    END IF
    !
    IF (rismt%cfft%dfftt%nnr > 0 .AND. (rismt%nrzs * rismt%ngxy) > 0) THEN
      CALL fw_lauefft_2xy(rismt%lfft, rismt%gr(:, iiq), ggz(:, iiq), rismt%nrzs, 1)
    END IF
  END DO
  !
  ! ... make qsol
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    nv    = iuniq_to_nsite(iq)
    isolV = isite_to_isolV(iv)
    iatom = isite_to_iatom(iv)
    rhov1 = DBLE(nv) * solVs(isolV)%density
    rhov2 = DBLE(nv) * solVs(isolV)%subdensity
    qv    = solVs(isolV)%charge(iatom)
    !
    rismt%qsol(iiq) = 0.0_DP
    !
    IF (rismt%lfft%gxystart > 1) THEN
      fac1 = qv * rhov1 * area_xy * dz
      fac2 = qv * rhov2 * area_xy * dz
      !
      qtmp = 0.0_DP
!$omp parallel do default(shared) private(irz) reduction(+:qtmp)
      DO irz = 1, (rismt%lfft%izleft_start - 1)
        qtmp = qtmp + fac2 * (DBLE(rismt%hsgz(irz, iiq) + rismt%hlgz(irz, iiq)))
      END DO
!$omp end parallel do
      rismt%qsol(iiq) = rismt%qsol(iiq) + qtmp
      !
      qtmp = 0.0_DP
!$omp parallel do default(shared) private(irz, iirz) reduction(+:qtmp)
      DO irz = rismt%lfft%izleft_start, rismt%lfft%izleft_gedge
        iirz = irz - rismt%lfft%izcell_start + 1
        qtmp = qtmp + fac2 * DBLE(ggz(iirz, iiq))
      END DO
!$omp end parallel do
      rismt%qsol(iiq) = rismt%qsol(iiq) + qtmp
      !
      qtmp = 0.0_DP
!$omp parallel do default(shared) private(irz, iirz) reduction(+:qtmp)
      DO irz = rismt%lfft%izright_gedge, rismt%lfft%izright_end
        iirz = irz - rismt%lfft%izcell_start + 1
        qtmp = qtmp + fac1 * DBLE(ggz(iirz, iiq))
      END DO
!$omp end parallel do
      rismt%qsol(iiq) = rismt%qsol(iiq) + qtmp
      !
      qtmp = 0.0_DP
!$omp parallel do default(shared) private(irz) reduction(+:qtmp)
      DO irz = (rismt%lfft%izright_end + 1), rismt%lfft%nrz
        qtmp = qtmp + fac1 * (DBLE(rismt%hsgz(irz, iiq) + rismt%hlgz(irz, iiq)))
      END DO
!$omp end parallel do
      rismt%qsol(iiq) = rismt%qsol(iiq) + qtmp
      !
    END IF
  END DO
  !
  IF (rismt%nsite > 0) THEN
    CALL mp_sum(rismt%qsol, rismt%mp_site%intra_sitg_comm)
  END IF
  !
  ! ... make qtot
  rismt%qtot = 0.0_DP
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq = iq - rismt%mp_site%isite_start + 1
    rismt%qtot = rismt%qtot + rismt%qsol(iiq)
  END DO
  !
  CALL mp_sum(rismt%qtot, rismt%mp_site%inter_sitg_comm)
  !
  ! ... make rhog (which is Laue-rep.)
  IF (rismt%nrzl * rismt%ngxy > 0) THEN
    rismt%rhog(:) = C_ZERO
  END IF
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    nv    = iuniq_to_nsite(iq)
    isolV = isite_to_isolV(iv)
    iatom = isite_to_iatom(iv)
    rhov1 = DBLE(nv) * solVs(isolV)%density
    rhov2 = DBLE(nv) * solVs(isolV)%subdensity
    qv    = solVs(isolV)%charge(iatom)
    !
    DO igxy = 1, rismt%ngxy
      jgxy = (igxy - 1) * rismt%nrzl
      kgxy = (igxy - 1) * rismt%nrzs
      !
!$omp parallel do default(shared) private(irz)
      DO irz = 1, (rismt%lfft%izleft_start - 1)
        rismt%rhog(irz + jgxy) = rismt%rhog(irz + jgxy) &
        & + qv * rhov2 * (rismt%hsgz(irz + jgxy, iiq) + rismt%hlgz(irz + jgxy, iiq))
      END DO
!$omp end parallel do
      !
!$omp parallel do default(shared) private(irz, iirz)
      DO irz = rismt%lfft%izleft_start, rismt%lfft%izleft_gedge
        iirz = irz - rismt%lfft%izcell_start + 1
        rismt%rhog(irz + jgxy) = rismt%rhog(irz + jgxy) &
        & + qv * rhov2 * ggz(iirz + kgxy, iiq)
      END DO
!$omp end parallel do
      !
!$omp parallel do default(shared) private(irz, iirz)
      DO irz = rismt%lfft%izright_gedge, rismt%lfft%izright_end
        iirz = irz - rismt%lfft%izcell_start + 1
        rismt%rhog(irz + jgxy) = rismt%rhog(irz + jgxy) &
        & + qv * rhov1 * ggz(iirz + kgxy, iiq)
      END DO
!$omp end parallel do
      !
!$omp parallel do default(shared) private(irz)
      DO irz = (rismt%lfft%izright_end + 1), rismt%lfft%nrz
        rismt%rhog(irz + jgxy) = rismt%rhog(irz + jgxy) &
        & + qv * rhov1 * (rismt%hsgz(irz + jgxy, iiq) + rismt%hlgz(irz + jgxy, iiq))
      END DO
!$omp end parallel do
      !
    END DO
  END DO
  !
  IF (rismt%nrzl * rismt%ngxy > 0) THEN
    CALL mp_sum(rismt%rhog, rismt%mp_site%inter_sitg_comm)
  END IF
  !
  ! ... truncate rhog
  charge0 = 0.0_DP
  izright_tail = 0
  izleft_tail  = 0
  !
  IF (rismt%lfft%gxystart > 1) THEN
    izleft_tail = 1
    DO irz = 1, rismt%lfft%izleft_gedge
      IF (ABS(DBLE(rismt%rhog(irz))) < RHOG_THRESHOLD .OR. &
      &   ABS(irz - rismt%lfft%izcell_start) <= RHOG_NEDGE) THEN
        rismt%rhog(irz) = C_ZERO
      ELSE
        izleft_tail = irz
        EXIT
      END IF
    END DO
    !
    izright_tail = rismt%lfft%nrz
    DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
      iirz = rismt%lfft%nrz + rismt%lfft%izright_gedge - irz
      IF (ABS(DBLE(rismt%rhog(iirz))) < RHOG_THRESHOLD .OR. &
      &   ABS(iirz - rismt%lfft%izcell_end) <= RHOG_NEDGE) THEN
        rismt%rhog(iirz) = C_ZERO
      ELSE
        izright_tail = iirz
        EXIT
      END IF
    END DO
    !
    fac = area_xy * dz
    DO irz = izleft_tail, rismt%lfft%izleft_gedge
      charge0 = charge0 + fac * rismt%rhog(irz)
    END DO
    DO irz = rismt%lfft%izright_gedge, izright_tail
      charge0 = charge0 + fac * rismt%rhog(irz)
    END DO
  END IF
  !
  CALL mp_sum(charge0,      rismt%mp_site%intra_sitg_comm)
  CALL mp_sum(izright_tail, rismt%mp_site%intra_sitg_comm)
  CALL mp_sum(izleft_tail,  rismt%mp_site%intra_sitg_comm)
  !
  ! ... renormalize rhog
  IF (ABS(charge0 - charge) > eps6) THEN
    IF (ABS(charge0) > eps6 .AND. ABS(charge) > eps6) THEN
      IF (rismt%lfft%gxystart > 1) THEN
        !
!$omp parallel do default(shared) private(irz)
        DO irz = izleft_tail, rismt%lfft%izleft_gedge
          rismt%rhog(irz) = rismt%rhog(irz) / charge0 * charge
        END DO
!$omp end parallel do
        !
!$omp parallel do default(shared) private(irz)
        DO irz = rismt%lfft%izright_gedge, izright_tail
          rismt%rhog(irz) = rismt%rhog(irz) / charge0 * charge
        END DO
!$omp end parallel do
        !
      END IF
      !
    ELSE
      IF (rismt%lfft%gxystart > 1) THEN
        nright = MAX(0, izright_tail - rismt%lfft%izright_gedge + 1)
        nleft  = MAX(0, rismt%lfft%izleft_gedge - izleft_tail + 1)
        fac    = area_xy * dz * DBLE(nright + nleft)
        !
!$omp parallel do default(shared) private(irz)
        DO irz = izleft_tail, rismt%lfft%izleft_gedge
          rismt%rhog(irz) = rismt%rhog(irz) - (charge0 - charge) / fac
        END DO
!$omp end parallel do
        !
!$omp parallel do default(shared) private(irz)
        DO irz = rismt%lfft%izright_gedge, izright_tail
          rismt%rhog(irz) = rismt%rhog(irz) - (charge0 - charge) / fac
        END DO
!$omp end parallel do
        !
      END IF
    END IF
    !
    WRITE(stdout, '(/,5X,"solvent charge ",F10.5, &
                    & ", renormalised to ",F10.5)') charge0, charge
    !
  ELSE
    WRITE(stdout, '(/,5X,"solvent charge ",F10.5)') charge0
    !
  END IF
  !
  ! ... make vpot
  CALL solvation_esm_potential(rismt, ireference, vsol0, ierr)
  IF (ierr /= IERR_RISM_NULL) THEN
    RETURN
  END IF
  !
  ! ... make rhog_pbc and vpot_pbc
  CALL solvation_pbc(rismt, ierr)
  IF (ierr /= IERR_RISM_NULL) THEN
    RETURN
  END IF
  !
  ! ... make esol
  rismt%esol = 0.0_DP
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    nv    = iuniq_to_nsite(iq)
    isolV = isite_to_isolV(iv)
    rhov  = DBLE(nv) * solVs(isolV)%density
    ! this version only supports G.F.
    !rismt%esol = rismt%esol + rhov * rismt%usol(iiq)
    rismt%esol = rismt%esol + rhov * rismt%usol_GF(iiq)
  END DO
  !
  CALL mp_sum(rismt%esol, rismt%mp_site%inter_sitg_comm)
  !
  ! ... make vsol (contribution of reference level shifting)
  rismt%vsol = 0.5_DP * vsol0 * charge
  !
  ! ... deallocate memory
  IF (rismt%nrzs * rismt%ngxy * rismt%nsite > 0) THEN
    DEALLOCATE(ggz)
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE solvation_lauerism
