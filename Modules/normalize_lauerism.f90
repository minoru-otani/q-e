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
  INTEGER               :: nright
  INTEGER               :: nleft
  INTEGER               :: igxy
  INTEGER               :: jgxy
  INTEGER               :: izright_tail
  INTEGER               :: izleft_tail
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
  REAL(DP)              :: vol1
  REAL(DP)              :: vol2
  REAL(DP)              :: charge0
  REAL(DP)              :: ntmp
  REAL(DP)              :: hr0
  REAL(DP)              :: hr1
  REAL(DP)              :: hr2
  REAL(DP), ALLOCATABLE :: rhoz(:)
  REAL(DP), ALLOCATABLE :: nsol(:)
  REAL(DP), ALLOCATABLE :: msol(:)
  REAL(DP), ALLOCATABLE :: qsol(:)
  COMPLEX(DP)           :: cpad
  !
  INTEGER,     PARAMETER :: RHOZ_NEDGE     = 3  ! to avoid noise at edges of unit-cell
  REAL(DP),    PARAMETER :: RHOZ_THRESHOLD = 1.0E-5_DP
  COMPLEX(DP), PARAMETER :: C_ZERO = CMPLX( 0.0_DP, 0.0_DP, kind=DP)
  COMPLEX(DP), PARAMETER :: C_MONE = CMPLX(-1.0_DP, 0.0_DP, kind=DP)
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
  IF (rismt%nrzl > 0) THEN
    ALLOCATE(rhoz(rismt%nrzl))
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
  ! ... rhoz: planar average of rho(r), in expand-cell
  ! ... NOTE: rhoz is defined only if gxystart > 1
  IF (rismt%nrzl > 0) THEN
    rhoz(:) = 0.0_DP
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
    IF (rismt%lfft%gxystart > 1) THEN
      !
!$omp parallel do default(shared) private(irz)
      DO irz = 1, rismt%lfft%izleft_gedge
        rhoz(irz) = rhoz(irz) + qv * rhov2 * DBLE(rismt%hsgz(irz, iiq) + rismt%hlgz(irz, iiq))
      END DO
!$omp end parallel do
      !
!$omp parallel do default(shared) private(irz)
      DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
        rhoz(irz) = rhoz(irz) + qv * rhov1 * DBLE(rismt%hsgz(irz, iiq) + rismt%hlgz(irz, iiq))
      END DO
!$omp end parallel do
      !
    END IF
    !
  END DO
  !
  IF (rismt%nrzl > 0) THEN
    CALL mp_sum(rhoz, rismt%mp_site%inter_sitg_comm)
  END IF
  !
  ! ... truncate rhoz
  izright_tail = 0
  izleft_tail  = 0
  !
  IF (rismt%lfft%gxystart > 1) THEN
    !
    izleft_tail = 1
    DO irz = 1, rismt%lfft%izleft_gedge
      IF (ABS(rhoz(irz)) < RHOZ_THRESHOLD .OR. &
      &   ABS(irz - rismt%lfft%izcell_start) <= RHOZ_NEDGE) THEN
        rhoz(irz) = 0.0_DP
      ELSE
        izleft_tail = irz
        EXIT
      END IF
    END DO
    !
    izright_tail = rismt%lfft%nrz
    DO irz = rismt%lfft%izright_gedge, rismt%lfft%nrz
      iirz = rismt%lfft%nrz + rismt%lfft%izright_gedge - irz
      IF (ABS(rhoz(iirz)) < RHOZ_THRESHOLD .OR. &
      &   ABS(iirz - rismt%lfft%izcell_end) <= RHOZ_NEDGE) THEN
        rhoz(iirz) = 0.0_DP
      ELSE
        izright_tail = iirz
        EXIT
      END IF
    END DO
    !
  END IF
  !
  CALL mp_sum(izright_tail, rismt%mp_site%intra_sitg_comm)
  CALL mp_sum(izleft_tail,  rismt%mp_site%intra_sitg_comm)
  !
  nright = MAX(0, izright_tail - rismt%lfft%izright_gedge + 1)
  nleft  = MAX(0, rismt%lfft%izleft_gedge - izleft_tail + 1)
  vol1   = dvol * DBLE(nright)
  vol2   = dvol * DBLE(nleft)
  !
  ! ... truncate hgz, hsgz, hlgz
  IF (rismt%nsite > 0) THEN
    !
    IF (expand) THEN
      ! ... expand-cell
      ! ... Gxy = 0
      IF (rismt%lfft%gxystart > 1) THEN
        !
        IF (rismt%lfft%xleft) THEN
          cpad = C_ZERO
        ELSE
          cpad = C_MONE
        END IF
        !
!$omp parallel do default(shared) private(irz)
        DO irz = 1, (izleft_tail - 1)
          rismt%hsgz(irz, :) = cpad
          rismt%hlgz(irz, :) = C_ZERO
        END DO
!$omp end parallel do
        !
        IF (rismt%lfft%xright) THEN
          cpad = C_ZERO
        ELSE
          cpad = C_MONE
        END IF
        !
!$omp parallel do default(shared) private(irz)
        DO irz = (izright_tail + 1), rismt%lfft%nrz
          rismt%hsgz(irz, :) = cpad
          rismt%hlgz(irz, :) = C_ZERO
        END DO
!$omp end parallel do
        !
      END IF
      !
      ! ... Gxy /= 0
      DO igxy = rismt%lfft%gxystart, rismt%lfft%ngxy
        jgxy = (igxy - 1) * rismt%nrzl
        !
!$omp parallel do default(shared) private(irz)
        DO irz = 1, (izleft_tail - 1)
          rismt%hsgz(irz + jgxy, :) = C_ZERO
          rismt%hlgz(irz + jgxy, :) = C_ZERO
        END DO
!$omp end parallel do
        !
!$omp parallel do default(shared) private(irz)
        DO irz = (izright_tail + 1), rismt%lfft%nrz
          rismt%hsgz(irz + jgxy, :) = C_ZERO
          rismt%hlgz(irz + jgxy, :) = C_ZERO
        END DO
!$omp end parallel do
        !
      END DO
      !
    ELSE
      ! ... unit-cell
      ! ... Gxy = 0
      IF (rismt%lfft%gxystart > 1) THEN
        !
        IF (rismt%lfft%xleft) THEN
          cpad = C_ZERO
        ELSE
          cpad = C_MONE
        END IF
        !
        DO irz = rismt%lfft%izcell_start, (izleft_tail - 1)
          iirz = irz - rismt%lfft%izcell_start + 1
          rismt%hgz(iirz, :) = cpad
        END DO
        !
        IF (rismt%lfft%xright) THEN
          cpad = C_ZERO
        ELSE
          cpad = C_MONE
        END IF
        !
        DO irz = (izright_tail + 1), rismt%lfft%izcell_end
          iirz = irz - rismt%lfft%izcell_start + 1
          rismt%hgz(iirz, :) = cpad
        END DO
        !
      END IF
      !
      ! ... Gxy /= 0
      DO igxy = rismt%lfft%gxystart, rismt%lfft%ngxy
        jgxy = (igxy - 1) * rismt%nrzs
        !
        DO irz = rismt%lfft%izcell_start, (izleft_tail - 1)
          iirz = irz - rismt%lfft%izcell_start + 1
          rismt%hgz(iirz + jgxy, :) = C_ZERO
        END DO
        !
        DO irz = (izright_tail + 1), rismt%lfft%izcell_end
          iirz = irz - rismt%lfft%izcell_start + 1
          rismt%hgz(iirz + jgxy, :) = C_ZERO
        END DO
        !
      END DO
      !
    END IF
    !
  END IF
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
  ! ... sum mean numbers and charges of solvent atoms in a molecule
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
        vrho  = vol1 * rhov1 + vol2 * rhov2
        IF (ABS(vrho) <= eps8) THEN  ! will not be occurred
          CALL errore('normalize_lauerism', 'vrho is zero', 1)
        END IF
        hr0   = (msol(isolV) - nsol(iiq)) / vrho
        !
        IF (expand) THEN
          ! ... expand-cell
!$omp parallel do default(shared) private(irz)
          DO irz = izleft_tail, rismt%lfft%izleft_gedge
            ! correct only short-range (for convenience)
            rismt%hsgz(irz, iiq) = rismt%hsgz(irz, iiq) + CMPLX(hr0, 0.0_DP, kind=DP)
          END DO
!$omp end parallel do
          !
!$omp parallel do default(shared) private(irz)
          DO irz = rismt%lfft%izright_gedge, izright_tail
            ! correct only short-range (for convenience)
            rismt%hsgz(irz, iiq) = rismt%hsgz(irz, iiq) + CMPLX(hr0, 0.0_DP, kind=DP)
          END DO
!$omp end parallel do
          !
        ELSE
          ! ... unit-cell
          DO irz = rismt%lfft%izcell_start, rismt%lfft%izleft_gedge
            iirz = irz - rismt%lfft%izcell_start + 1
            rismt%hgz(iirz, iiq) = rismt%hgz(iirz, iiq) + CMPLX(hr0, 0.0_DP, kind=DP)
          END DO
          !
          DO irz = rismt%lfft%izright_gedge, rismt%lfft%izcell_end
            iirz = irz - rismt%lfft%izcell_start + 1
            rismt%hgz(iirz, iiq) = rismt%hgz(iirz, iiq) + CMPLX(hr0, 0.0_DP, kind=DP)
          END DO
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
  ! ... qrho = sum(qv^2 * rhov^2)
  ! ... NOTE: qrho1 and qrho2 are defined only if gxystart > 1
  qrho1 = 0.0_DP
  qrho2 = 0.0_DP
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
      qrho1 = qrho1 + DBLE(nv) * qv * qv * rhov1 * rhov1
      qrho2 = qrho2 + DBLE(nv) * qv * qv * rhov2 * rhov2
    END IF
  END DO
  !
  CALL mp_sum(qrho1, rismt%mp_site%inter_sitg_comm)
  CALL mp_sum(qrho2, rismt%mp_site%inter_sitg_comm)
  !
  ! ... renormalize hgz, hsgz, hlgz
  IF (rismt%lfft%gxystart > 1) THEN
    !
    IF (ABS(charge - charge0) > eps8) THEN
      !
      vqrho = vol1 * qrho1 + vol2 * qrho2
      IF (ABS(vqrho) <= eps8) THEN  ! will not be occurred
        CALL errore('normalize_lauerism', 'vqrho is zero', 1)
      END IF
      hr0   = (charge - charge0) / vqrho
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
          DO irz = izleft_tail, rismt%lfft%izleft_gedge
            ! correct only short-range (for convenience)
            rismt%hsgz(irz, iiq) = rismt%hsgz(irz, iiq) + CMPLX(hr2, 0.0_DP, kind=DP)
          END DO
!$omp end parallel do
          !
!$omp parallel do default(shared) private(irz)
          DO irz = rismt%lfft%izright_gedge, izright_tail
            ! correct only short-range (for convenience)
            rismt%hsgz(irz, iiq) = rismt%hsgz(irz, iiq) + CMPLX(hr1, 0.0_DP, kind=DP)
          END DO
!$omp end parallel do
          !
        ELSE
          ! ... unit-cell
          DO irz = rismt%lfft%izcell_start, rismt%lfft%izleft_gedge
            iirz = irz - rismt%lfft%izcell_start + 1
            rismt%hgz(iirz, iiq) = rismt%hgz(iirz, iiq) + CMPLX(hr2, 0.0_DP, kind=DP)
          END DO
          !
          DO irz = rismt%lfft%izright_gedge, rismt%lfft%izcell_end
            iirz = irz - rismt%lfft%izcell_start + 1
            rismt%hgz(iirz, iiq) = rismt%hgz(iirz, iiq) + CMPLX(hr1, 0.0_DP, kind=DP)
          END DO
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
  IF (rismt%nrzl > 0) THEN
    DEALLOCATE(rhoz)
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
