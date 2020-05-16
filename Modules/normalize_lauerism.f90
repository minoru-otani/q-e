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
SUBROUTINE normalize_lauerism(rismt, charge, expand, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... normalize total correlations to remove noise of solvent charge density,
  ! ... which is defined as
  !                   ----
  !        rho(gxy,z) >    q(v) * rho(v) * h(v; gxy,z)
  !                   ----
  !                     v
  !
  ! ... NOTE: h(gxy,z) is used, not g(gxy,z),
  ! ...       to obtain electrostatic consistent chemical potential of solvation.
  !
  ! ... Variables:
  ! ...   charge: total charge of solvent system
  ! ...   expand: use expand-cell(.TRUE.) or unit-cell(.FALSE.)
  !
  USE cell_base,      ONLY : at, alat
  USE constants,      ONLY : eps8
  USE err_rism,       ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,          ONLY : DP
  USE mp,             ONLY : mp_sum
  USE rism,           ONLY : rism_type, ITYPE_LAUERISM
  USE solvmol,        ONLY : get_nuniq_in_solVs, nsolV, solVs, iuniq_to_nsite, &
                           & iuniq_to_isite, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  REAL(DP),        INTENT(IN)    :: charge
  LOGICAL,         INTENT(IN)    :: expand
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER               :: nq
  INTEGER               :: iq
  INTEGER               :: iiq
  INTEGER               :: iv
  INTEGER               :: nv
  INTEGER               :: isolV
  INTEGER               :: iatom
  INTEGER               :: irz
  INTEGER               :: iirz
  INTEGER               :: igxy
  INTEGER               :: jgxy
  INTEGER,  ALLOCATABLE :: izright_tail(:)
  INTEGER,  ALLOCATABLE :: izleft_tail(:)
  REAL(DP)              :: qv
  REAL(DP)              :: rhov1
  REAL(DP)              :: rhov2
  REAL(DP)              :: qrho1
  REAL(DP)              :: qrho2
  REAL(DP)              :: vrho
  REAL(DP)              :: vqrho
  REAL(DP)              :: dz
  REAL(DP)              :: area_xy
  REAL(DP)              :: dvol
  REAL(DP)              :: vtmp
  REAL(DP)              :: ntmp
  REAL(DP)              :: charge0
  REAL(DP)              :: hz
  REAL(DP)              :: hr0
  REAL(DP)              :: hr1
  REAL(DP)              :: hr2
  REAL(DP), ALLOCATABLE :: vol1(:)
  REAL(DP), ALLOCATABLE :: vol2(:)
  REAL(DP), ALLOCATABLE :: hwei(:,:)
  REAL(DP), ALLOCATABLE :: nsol(:)
  REAL(DP), ALLOCATABLE :: msol(:)
  REAL(DP), ALLOCATABLE :: qsol(:)
  !
  INTEGER,     PARAMETER :: HZ_EDGE  = 3  ! to avoid noise at edges of unit-cell
  REAL(DP),    PARAMETER :: HZ_THR   = 1.0E-3_DP
  REAL(DP),    PARAMETER :: HZ_SMEAR = 2.0_DP  ! in bohr
  !
  REAL(DP),    EXTERNAL  :: qe_erfc
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
  IF (rismt%ngxy < rismt%lfft%ngxy) THEN
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
  IF (rismt%nr < rismt%cfft%dfftt%nnr) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... allocate memory
  IF (rismt%nsite > 0) THEN
    ALLOCATE(izright_tail(rismt%nsite))
    ALLOCATE(izleft_tail( rismt%nsite))
    ALLOCATE(vol1(rismt%nsite))
    ALLOCATE(vol2(rismt%nsite))
  END IF
  IF (rismt%lfft%nrz * rismt%nsite > 0) THEN
    ALLOCATE(hwei(rismt%lfft%nrz, rismt%nsite))
  END IF
  IF (rismt%nsite > 0) THEN
    ALLOCATE(nsol(rismt%nsite))
  END IF
  IF (nsolV > 0) THEN
    ALLOCATE(msol(nsolV))
    ALLOCATE(qsol(nsolV))
  END IF
  !
  ! ... set variables
  dz      = rismt%lfft%zstep * alat
  area_xy = ABS(at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1)) * alat * alat
  dvol    = area_xy * dz
  !
  ! ...
  ! ... Define domain of total correlations
  ! ......................................................................................
  !
  ! ... detect truncating positions for each solvent site
  IF (rismt%nsite > 0) THEN
    izright_tail = 0
    izleft_tail  = 0
  END IF
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq = iq - rismt%mp_site%isite_start + 1
    !
    IF (rismt%lfft%gxystart > 1) THEN
      !
      izleft_tail(iiq) = 1
      DO irz = 1, rismt%lfft%izleft_gedge
        hz   = DBLE(rismt%hsgz(irz, iiq) + rismt%hlgz(irz, iiq))
        IF (ABS(hz) < HZ_THR .OR. ABS(irz - rismt%lfft%izcell_start) <= HZ_EDGE) THEN
          ! NOP
        ELSE
          izleft_tail(iiq) = irz
          EXIT
        END IF
      END DO
      !
      izright_tail(iiq) = rismt%lfft%nrz
      DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
        iirz = rismt%lfft%nrz + rismt%lfft%izright_gedge - irz
        hz   = DBLE(rismt%hsgz(iirz, iiq) + rismt%hlgz(iirz, iiq))
        IF (ABS(hz) < HZ_THR .OR. ABS(iirz - rismt%lfft%izcell_end) <= HZ_EDGE) THEN
          ! NOP
        ELSE
          izright_tail(iiq) = iirz
          EXIT
        END IF
      END DO
      !
    END IF
    !
  END DO
  !
  IF (rismt%nsite > 0) THEN
    CALL mp_sum(izright_tail, rismt%mp_site%intra_sitg_comm)
    CALL mp_sum(izleft_tail,  rismt%mp_site%intra_sitg_comm)
  END IF
  !
  ! ... calculate weight for total correlations
  IF (rismt%lfft%nrz * rismt%nsite > 0) THEN
    hwei = 0.0_DP
  END IF
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq = iq - rismt%mp_site%isite_start + 1
    !
!$omp parallel do default(shared) private(irz)
    DO irz = 1, rismt%lfft%izleft_gedge
      hwei(irz, iiq) = 0.5_DP * qe_erfc(DBLE(izleft_tail(iiq) - irz ) * dz / HZ_SMEAR)
    END DO
!$omp end parallel do
    !
!$omp parallel do default(shared) private(irz)
    DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
      hwei(irz, iiq) = 0.5_DP * qe_erfc(DBLE(irz - izright_tail(iiq)) * dz / HZ_SMEAR)
    END DO
!$omp end parallel do
    !
  END DO
  !
  ! ... truncate hgz, hsgz, hlgz
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq = iq - rismt%mp_site%isite_start + 1
    !
    IF (expand) THEN
      ! ... expand-cell
      DO igxy = 1, rismt%lfft%ngxy
        jgxy = (igxy - 1) * rismt%nrzl
        !
!$omp parallel do default(shared) private(irz)
        DO irz = 1, rismt%lfft%izleft_gedge
          rismt%hsgz(irz + jgxy, iiq) = hwei(irz, iiq) * rismt%hsgz(irz + jgxy, iiq)
          rismt%hlgz(irz + jgxy, iiq) = hwei(irz, iiq) * rismt%hlgz(irz + jgxy, iiq)
        END DO
!$omp end parallel do
        !
!$omp parallel do default(shared) private(irz)
        DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
          rismt%hsgz(irz + jgxy, iiq) = hwei(irz, iiq) * rismt%hsgz(irz + jgxy, iiq)
          rismt%hlgz(irz + jgxy, iiq) = hwei(irz, iiq) * rismt%hlgz(irz + jgxy, iiq)
        END DO
!$omp end parallel do
        !
      END DO
      !
    ELSE
      ! ... unit-cell
      DO igxy = 1, rismt%lfft%ngxy
        jgxy = (igxy - 1) * rismt%nrzs
        !
        DO irz = rismt%lfft%izcell_start, rismt%lfft%izleft_gedge
          iirz = irz - rismt%lfft%izcell_start + 1
          rismt%hgz(iirz + jgxy, iiq) = hwei(irz, iiq) * rismt%hgz(iirz + jgxy, iiq)
        END DO
        !
        DO irz = rismt%lfft%izright_gedge, rismt%lfft%izcell_end
          iirz = irz - rismt%lfft%izcell_start + 1
          rismt%hgz(iirz + jgxy, iiq) = hwei(irz, iiq) * rismt%hgz(iirz + jgxy, iiq)
        END DO
        !
      END DO
      !
    END IF
    !
  END DO
  !
  ! ... volume of solvent domain
  ! ... NOTE: vol1 and vol2 are defined only if gxystart > 1
  IF (rismt%nsite > 0) THEN
    vol1 = 0.0_DP
    vol2 = 0.0_DP
  END IF
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq = iq - rismt%mp_site%isite_start + 1
    !
    IF (rismt%lfft%gxystart > 1) THEN
      !
      vtmp = 0.0_DP
!$omp parallel do default(shared) private(irz) reduction(+:vtmp)
      DO irz = 1, rismt%lfft%izleft_gedge
        vtmp = vtmp + dvol * hwei(irz, iiq)
      END DO
!$omp end parallel do
      vol2(iiq) = vtmp
      !
      vtmp = 0.0_DP
!$omp parallel do default(shared) private(irz) reduction(+:vtmp)
      DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
        vtmp = vtmp + dvol * hwei(irz, iiq)
      END DO
!$omp end parallel do
      vol1(iiq) = vtmp
      !
    END IF
    !
  END DO
  !
  ! ...
  ! ... Correct stoichiometry of solvent molecules
  ! ......................................................................................
  !
  ! ... calculate total numbers of solvent atoms
  ! ... NOTE: nsol is defined only if gxystart > 1
  IF (rismt%nsite > 0) THEN
    nsol(:) = 0.0_DP
  END IF
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    isolV = isite_to_isolV(iv)
    rhov1 = solVs(isolV)%density
    rhov2 = solVs(isolV)%subdensity
    !
    IF (rismt%lfft%gxystart > 1) THEN
      !
      ntmp = 0.0_DP
!$omp parallel do default(shared) private(irz) reduction(+:ntmp)
      DO irz = 1, rismt%lfft%izleft_gedge
        ntmp = ntmp + (DBLE(rismt%hsgz(irz, iiq) + rismt%hlgz(irz, iiq)) + 1.0_DP)
      END DO
!$omp end parallel do
      nsol(iiq) = nsol(iiq) + rhov2 * dvol * ntmp
      !
      ntmp = 0.0_DP
!$omp parallel do default(shared) private(irz) reduction(+:ntmp)
      DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
        ntmp = ntmp + (DBLE(rismt%hsgz(irz, iiq) + rismt%hlgz(irz, iiq)) + 1.0_DP)
      END DO
!$omp end parallel do
      nsol(iiq) = nsol(iiq) + rhov1 * dvol * ntmp
      !
    END IF
    !
  END DO
  !
  ! ... sum numbers and charges of solvent atoms in a molecule
  ! ... NOTE: msol and qsol are defined only if gxystart > 1
  IF (nsolV > 0) THEN
    msol(:) = 0.0_DP
    qsol(:) = 0.0_DP
  END IF
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    nv    = iuniq_to_nsite(iq)
    isolV = isite_to_isolV(iv)
    iatom = isite_to_iatom(iv)
    qv    = solVs(isolV)%charge(iatom)
    !
    IF (rismt%lfft%gxystart > 1) THEN
      msol(isolV) = msol(isolV) + DBLE(nv) * nsol(iiq)
      qsol(isolV) = qsol(isolV) + DBLE(nv) * qv
    END IF
  END DO
  !
  IF (nsolV > 0) THEN
    CALL mp_sum(msol, rismt%mp_site%inter_sitg_comm)
    CALL mp_sum(qsol, rismt%mp_site%inter_sitg_comm)
  END IF
  !
  DO isolV = 1, nsolV
    IF (solVs(isolV)%natom > 0) THEN
      msol(isolV) = msol(isolV) / DBLE(solVs(isolV)%natom)
      qsol(isolV) = qsol(isolV) / DBLE(solVs(isolV)%natom)
    ELSE
      msol(isolV) = 0.0_DP
      qsol(isolV) = 0.0_DP
    END IF
  END DO
  !
  ! ... renormalize hgz, hsgz, hlgz
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    isolV = isite_to_isolV(iv)
    !
    IF (rismt%lfft%gxystart > 1) THEN
      !
      IF (ABS(msol(isolV) - nsol(iiq)) > eps8) THEN
        rhov1 = solVs(isolV)%density
        rhov2 = solVs(isolV)%subdensity
        vrho  = vol1(iiq) * rhov1 + vol2(iiq) * rhov2
        !
        IF (ABS(vrho) <= eps8) THEN  ! will not be occurred
          CALL errore('normalize_lauerism', 'vrho is zero', 1)
        END IF
        hr0 = (msol(isolV) - nsol(iiq)) / vrho
        !
        IF (expand) THEN
          ! ... expand-cell
!$omp parallel do default(shared) private(irz)
          DO irz = 1, rismt%lfft%izleft_gedge
            ! correct only short-range (for convenience)
            rismt%hsgz(irz, iiq) = rismt%hsgz(irz, iiq) &
                               & + CMPLX(hwei(irz, iiq) * hr0, 0.0_DP, kind=DP)
          END DO
!$omp end parallel do
          !
!$omp parallel do default(shared) private(irz)
          DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
            ! correct only short-range (for convenience)
            rismt%hsgz(irz, iiq) = rismt%hsgz(irz, iiq) &
                               & + CMPLX(hwei(irz, iiq) * hr0, 0.0_DP, kind=DP)
          END DO
!$omp end parallel do
          !
        ELSE
          ! ... unit-cell
          DO irz = rismt%lfft%izcell_start, rismt%lfft%izleft_gedge
            iirz = irz - rismt%lfft%izcell_start + 1
            rismt%hgz(iirz, iiq) = rismt%hgz(iirz, iiq) &
                               & + CMPLX(hwei(irz, iiq) * hr0, 0.0_DP, kind=DP)
          END DO
          !
          DO irz = rismt%lfft%izright_gedge, rismt%lfft%izcell_end
            iirz = irz - rismt%lfft%izcell_start + 1
            rismt%hgz(iirz, iiq) = rismt%hgz(iirz, iiq) &
                               & + CMPLX(hwei(irz, iiq) * hr0, 0.0_DP, kind=DP)
          END DO
          !
        END IF
        !
      END IF
      !
    END IF
    !
  END DO
  !
  ! ...
  ! ... Correct total charge of solvent system
  ! ......................................................................................
  !
  ! ... calculate total charge
  ! ... NOTE: charge0 is defined only if gxystart > 1
  charge0 = 0.0_DP
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    nv    = iuniq_to_nsite(iq)
    isolV = isite_to_isolV(iv)
    !
    IF (rismt%lfft%gxystart > 1) THEN
      charge0 = charge0 + DBLE(nv) * qsol(isolV) * msol(isolV)
    END IF
  END DO
  !
  CALL mp_sum(charge0, rismt%mp_site%inter_sitg_comm)
  !
  ! ... vqrho = sum(vol * qv^2 * rhov^2)
  ! ... NOTE: vqrho is defined only if gxystart > 1
  vqrho = 0.0_DP
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    nv    = iuniq_to_nsite(iq)
    isolV = isite_to_isolV(iv)
    rhov1 = solVs(isolV)%density
    rhov2 = solVs(isolV)%subdensity
    !
    IF (rismt%lfft%gxystart > 1) THEN
      qv    = qsol(isolV)
      qrho1 = qv * qv * rhov1 * rhov1
      qrho2 = qv * qv * rhov2 * rhov2
      vqrho = vqrho + DBLE(nv) * (vol1(iiq) * qrho1 + vol2(iiq) * qrho2)
    END IF
  END DO
  !
  CALL mp_sum(vqrho, rismt%mp_site%inter_sitg_comm)
  !
  ! ... renormalize hgz, hsgz, hlgz
  IF (rismt%lfft%gxystart > 1) THEN
    !
    IF (ABS(charge - charge0) > eps8) THEN
      !
      IF (ABS(vqrho) <= eps8) THEN  ! will not be occurred
        CALL errore('normalize_lauerism', 'vqrho is zero', 1)
      END IF
      hr0 = (charge - charge0) / vqrho
      !
      DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
        iiq   = iq - rismt%mp_site%isite_start + 1
        iv    = iuniq_to_isite(1, iq)
        isolV = isite_to_isolV(iv)
        rhov1 = solVs(isolV)%density
        rhov2 = solVs(isolV)%subdensity
        qv    = qsol(isolV)
        hr1   = hr0 * qv * rhov1
        hr2   = hr0 * qv * rhov2
        !
        IF (expand) THEN
          ! ... expand-cell
!$omp parallel do default(shared) private(irz)
          DO irz = 1, rismt%lfft%izleft_gedge
            ! correct only short-range (for convenience)
            rismt%hsgz(irz, iiq) = rismt%hsgz(irz, iiq) &
                               & + CMPLX(hwei(irz, iiq) * hr2, 0.0_DP, kind=DP)
          END DO
!$omp end parallel do
          !
!$omp parallel do default(shared) private(irz)
          DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
            ! correct only short-range (for convenience)
            rismt%hsgz(irz, iiq) = rismt%hsgz(irz, iiq) &
                               & + CMPLX(hwei(irz, iiq) * hr1, 0.0_DP, kind=DP)
          END DO
!$omp end parallel do
          !
        ELSE
          ! ... unit-cell
          DO irz = rismt%lfft%izcell_start, rismt%lfft%izleft_gedge
            iirz = irz - rismt%lfft%izcell_start + 1
            rismt%hgz(iirz, iiq) = rismt%hgz(iirz, iiq) &
                               & + CMPLX(hwei(irz, iiq) * hr2, 0.0_DP, kind=DP)
          END DO
          !
          DO irz = rismt%lfft%izright_gedge, rismt%lfft%izcell_end
            iirz = irz - rismt%lfft%izcell_start + 1
            rismt%hgz(iirz, iiq) = rismt%hgz(iirz, iiq) &
                               & + CMPLX(hwei(irz, iiq) * hr1, 0.0_DP, kind=DP)
          END DO
          !
        END IF
        !
      END DO
      !
    END IF
    !
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
  ! ... deallocate memory
100 CONTINUE
  IF (rismt%nsite > 0) THEN
    DEALLOCATE(izright_tail)
    DEALLOCATE(izleft_tail)
    DEALLOCATE(vol1)
    DEALLOCATE(vol2)
  END IF
  IF (rismt%lfft%nrz * rismt%nsite > 0) THEN
    DEALLOCATE(hwei)
  END IF
  IF (rismt%nsite > 0) THEN
    DEALLOCATE(nsol)
  END IF
  IF (nsolV > 0) THEN
    DEALLOCATE(msol)
    DEALLOCATE(qsol)
  END IF
  !
END SUBROUTINE normalize_lauerism
