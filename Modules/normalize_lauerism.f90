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
SUBROUTINE normalize_lauerism(rismt, charge, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... normalize total correlations to remove noise from solvent charge density.
  !
  ! ... Variables:
  ! ...   charge:  total charge of solvent system
  !
  USE cell_base,      ONLY : at, alat
  USE constants,      ONLY : eps8
  USE err_rism,       ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE io_global,      ONLY : stdout
  USE kinds,          ONLY : DP
  USE mp,             ONLY : mp_sum
  USE rism,           ONLY : rism_type, ITYPE_LAUERISM
  USE solvmol,        ONLY : get_nuniq_in_solVs, solVs, iuniq_to_nsite, &
                           & iuniq_to_isite, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  REAL(DP),        INTENT(IN)    :: charge
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
  REAL(DP)              :: rhov1
  REAL(DP)              :: rhov2
  REAL(DP)              :: qv
  REAL(DP)              :: qrho1
  REAL(DP)              :: qrho2
  REAL(DP)              :: vqrho
  REAL(DP)              :: dz
  REAL(DP)              :: area_xy
  REAL(DP)              :: dvol
  REAL(DP)              :: vol1
  REAL(DP)              :: vol2
  REAL(DP)              :: charge0
  REAL(DP)              :: gr0
  REAL(DP)              :: gr1
  REAL(DP)              :: gr2
  REAL(DP), ALLOCATABLE :: gz(:,:)
  REAL(DP), ALLOCATABLE :: rhoz(:)
  !
  INTEGER,     PARAMETER :: RHOZ_NEDGE     = 3  ! to avoid noise at edges of unit-cell
  REAL(DP),    PARAMETER :: RHOZ_THRESHOLD = 1.0E-5_DP
  COMPLEX(DP), PARAMETER :: C_ZERO = CMPLX(0.0_DP, 0.0_DP, kind=DP)
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
  IF (rismt%nr < rismt%cfft%dfftt%nnr) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... allocate memory
  IF (rismt%nrzs * rismt%nsite > 0) THEN
    ALLOCATE(gz(rismt%nrzs, rismt%nsite))
  END IF
  IF (rismt%nrzl > 0) THEN
    ALLOCATE(rhoz(rismt%nrzl))
  END IF
  !
  ! ... set variables
  dz      = rismt%lfft%zstep * alat
  area_xy = ABS(at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1)) * alat * alat
  dvol    = area_xy * dz
  !
  ! ... qrho = sum(qv^2 * rhov^2)
  qrho1 = 0.0_DP
  qrho2 = 0.0_DP
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    nv    = iuniq_to_nsite(iq)
    isolV = isite_to_isolV(iv)
    iatom = isite_to_iatom(iv)
    rhov1 = solVs(isolV)%density
    rhov2 = solVs(isolV)%subdensity
    qv    = solVs(isolV)%charge(iatom)
    qrho1 = qrho1 + DBLE(nv) * qv * qv * rhov1 * rhov1
    qrho2 = qrho2 + DBLE(nv) * qv * qv * rhov2 * rhov2
  END DO
  !
  CALL mp_sum(qrho1, rismt%mp_site%inter_sitg_comm)
  CALL mp_sum(qrho2, rismt%mp_site%inter_sitg_comm)
  !
  ! ... gz: planar average of g(r), in unit-cell
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq = iq - rismt%mp_site%isite_start + 1
    IF (rismt%nrzs > 0) THEN
      gz(:, iiq) = 0.0_DP
    END IF
    !
    IF (rismt%nr > 0 .AND. rismt%nrzs > 0) THEN
      CALL planar_average(rismt%gr(:, iiq), gz(:, iiq))
    END IF
  END DO
  !
  ! ... rhoz: planar average of rho(r), in expand-cell
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
      DO irz = 1, (rismt%lfft%izleft_start - 1)
        rhoz(irz) = rhoz(irz) + qv * rhov2 * DBLE(rismt%hsgz(irz, iiq) + rismt%hlgz(irz, iiq))
      END DO
!$omp end parallel do
      !
!$omp parallel do default(shared) private(irz, iirz)
      DO irz = rismt%lfft%izleft_start, rismt%lfft%izleft_gedge
        iirz = irz - rismt%lfft%izcell_start + 1
        rhoz(irz) = rhoz(irz) + qv * rhov2 * gz(iirz, iiq)
      END DO
!$omp end parallel do
      !
!$omp parallel do default(shared) private(irz, iirz)
      DO irz = rismt%lfft%izright_gedge, rismt%lfft%izright_end
        iirz = irz - rismt%lfft%izcell_start + 1
        rhoz(irz) = rhoz(irz) + qv * rhov1 * gz(iirz, iiq)
      END DO
!$omp end parallel do
      !
!$omp parallel do default(shared) private(irz)
      DO irz = (rismt%lfft%izright_end + 1), rismt%lfft%nrz
        rhoz(irz) = rhoz(irz) + qv * rhov1 * DBLE(rismt%hsgz(irz, iiq) + rismt%hlgz(irz, iiq))
      END DO
!$omp end parallel do
      !
    END IF
  END DO
  !
  IF (rismt%nrzl > 0) THEN
    CALL mp_sum(rhoz, rismt%mp_site%inter_sitg_comm)
  END IF
  !
  ! ... truncate rhoz
  charge0 = 0.0_DP
  izright_tail = 0
  izleft_tail  = 0
  !
  IF (rismt%lfft%gxystart > 1) THEN
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
    DO irz = izleft_tail, rismt%lfft%izleft_gedge
      charge0 = charge0 + dvol * rhoz(irz)
    END DO
    DO irz = rismt%lfft%izright_gedge, izright_tail
      charge0 = charge0 + dvol * rhoz(irz)
    END DO
  END IF
  !
  CALL mp_sum(charge0,      rismt%mp_site%intra_sitg_comm)
  CALL mp_sum(izright_tail, rismt%mp_site%intra_sitg_comm)
  CALL mp_sum(izleft_tail,  rismt%mp_site%intra_sitg_comm)
  !
  ! ... truncate hsgz, hlgz, gr
  IF (rismt%nsite > 0) THEN
    !
    DO igxy = 1, rismt%ngxy
      jgxy = (igxy - 1) * rismt%nrzl
      !
      DO irz = 1, (izleft_tail - 1)
        rismt%hsgz(irz + jgxy, :) = C_ZERO
        rismt%hlgz(irz + jgxy, :) = C_ZERO
      END DO
      !
      DO irz = (izright_tail + 1), rismt%lfft%nrz
        rismt%hsgz(irz + jgxy, :) = C_ZERO
        rismt%hlgz(irz + jgxy, :) = C_ZERO
      END DO
    END DO
    !
    CALL truncate_gr(rismt%nsite, rismt%gr, izleft_tail, izright_tail)
    !
  END IF
  !
  ! ... renormalize hsgz, hlgz, gr
  IF (ABS(charge0 - charge) > eps8) THEN
    !
    nright = MAX(0, izright_tail - rismt%lfft%izright_gedge + 1)
    nleft  = MAX(0, rismt%lfft%izleft_gedge - izleft_tail + 1)
    vol1   = dvol * DBLE(nright)
    vol2   = dvol * DBLE(nleft)
    vqrho  = vol1 * qrho1 + vol2 * qrho2
    IF (ABS(vqrho) <= eps8) THEN  ! will not be occurred
      CALL errore('normalize_lauerism', 'vqrho is zero', 1)
    END IF
    gr0    = (charge0 - charge) / vqrho
    !
    DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      iiq   = iq - rismt%mp_site%isite_start + 1
      iv    = iuniq_to_isite(1, iq)
      isolV = isite_to_isolV(iv)
      iatom = isite_to_iatom(iv)
      rhov1 = solVs(isolV)%density
      rhov2 = solVs(isolV)%subdensity
      qv    = solVs(isolV)%charge(iatom)
      gr1   = gr0 * qv * rhov1
      gr2   = gr0 * qv * rhov2
      !
      IF (rismt%lfft%gxystart > 1) THEN
        !
!$omp parallel do default(shared) private(irz)
        DO irz = izleft_tail, rismt%lfft%izleft_gedge
          ! correct only short-range (for convenience)
          rismt%hsgz(irz, iiq) = rismt%hsgz(irz, iiq) - gr2
        END DO
!$omp end parallel do
        !
!$omp parallel do default(shared) private(irz)
        DO irz = rismt%lfft%izright_gedge, izright_tail
          ! correct only short-range (for convenience)
          rismt%hsgz(irz, iiq) = rismt%hsgz(irz, iiq) - gr1
        END DO
!$omp end parallel do
        !
      END IF
      !
      CALL renormalize_gr(rismt%gr(:, iiq), gr1, gr2, izleft_tail, izright_tail)
      !
    END DO
    !
    WRITE(stdout, '(/,5X,"solvent charge ",F10.5, &
                    & ", renormalised to ",F10.5)') charge0, charge
    !
  ELSE
    WRITE(stdout, '(/,5X,"solvent charge ",F10.5)') charge0
    !
  END IF
  !
  ! ... deallocate memory
  IF (rismt%nrzs * rismt%nsite > 0) THEN
    DEALLOCATE(gz)
  END IF
  IF (rismt%nrzl > 0) THEN
    DEALLOCATE(rhoz)
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
CONTAINS
  !
  !---------------------------------------------------------------------------
  SUBROUTINE planar_average(gr, gz)
    !---------------------------------------------------------------------------
    !
    ! ... planer average of g(r) -> g(z).
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)  :: gr(:)
    REAL(DP), INTENT(OUT) :: gz(:)
    !
    INTEGER :: ir
    INTEGER :: idx
    INTEGER :: idx0
    INTEGER :: i3min
    INTEGER :: i3max
    INTEGER :: i1, i2, i3
    INTEGER :: n1, n2, n3
    INTEGER :: nx1, nx2, nx3
    INTEGER :: nnr
    INTEGER :: irz
    !
    ! ... FFT box
    n1  = rismt%cfft%dfftt%nr1
    n2  = rismt%cfft%dfftt%nr2
    n3  = rismt%cfft%dfftt%nr3
    nx1 = rismt%cfft%dfftt%nr1x
    nx2 = rismt%cfft%dfftt%nr2x
    nx3 = rismt%cfft%dfftt%nr3x
    nnr = rismt%cfft%dfftt%nnr
    !
    ! ... calculate planar average
    gz(1:n3) = 0.0_DP
    !
    idx0  = nx1 * nx2 * rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
    i3min = rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
    i3max = rismt%cfft%dfftt%npp(rismt%cfft%dfftt%mype + 1) + i3min
    !
    DO ir = 1, nnr
      !
      idx = idx0 + ir - 1
      i3  = idx / (nx1 * nx2)
      IF (i3 < i3min .OR. i3 >= i3max .OR. i3 >= n3) THEN
        CYCLE
      END IF
      !
      idx = idx - (nx1 * nx2) * i3
      i2  = idx / nx1
      IF (i2 >= n2) THEN
        CYCLE
      END IF
      !
      idx = idx - nx1 * i2
      i1  = idx
      IF (i1 >= n1) THEN
        CYCLE
      END IF
      !
      IF (i3 >= (n3 - (n3 / 2))) THEN
        irz = i3 - n3
      ELSE
        irz = i3
      END IF
      irz = irz + (n3 / 2)
      irz = irz + 1
      !
      gz(irz) = gz(irz) + gr(ir)
      !
    END DO
    !
    CALL mp_sum(gz(1:n3), rismt%mp_site%intra_sitg_comm)
    !
    gz(1:n3) = gz(1:n3) / DBLE(n1 * n2)
    !
  END SUBROUTINE planar_average
  !
  !---------------------------------------------------------------------------
  SUBROUTINE truncate_gr(nsite, gr, irz_sta, irz_end)
    !---------------------------------------------------------------------------
    !
    ! ... truncate edge of g(r).
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(IN)  :: nsite
    REAL(DP), INTENT(OUT) :: gr(:,:)
    INTEGER,  INTENT(IN)  :: irz_sta
    INTEGER,  INTENT(IN)  :: irz_end
    !
    INTEGER :: ir
    INTEGER :: idx
    INTEGER :: idx0
    INTEGER :: i3min
    INTEGER :: i3max
    INTEGER :: i1, i2, i3
    INTEGER :: n1, n2, n3
    INTEGER :: nx1, nx2, nx3
    INTEGER :: nnr
    INTEGER :: irz
    INTEGER :: irz0
    !
    ! ... FFT box
    n1  = rismt%cfft%dfftt%nr1
    n2  = rismt%cfft%dfftt%nr2
    n3  = rismt%cfft%dfftt%nr3
    nx1 = rismt%cfft%dfftt%nr1x
    nx2 = rismt%cfft%dfftt%nr2x
    nx3 = rismt%cfft%dfftt%nr3x
    nnr = rismt%cfft%dfftt%nnr
    !
    ! ... calculate planar average
    idx0  = nx1 * nx2 * rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
    i3min = rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
    i3max = rismt%cfft%dfftt%npp(rismt%cfft%dfftt%mype + 1) + i3min
    irz0  = rismt%lfft%izcell_start
    !
    DO ir = 1, nnr
      !
      idx = idx0 + ir - 1
      i3  = idx / (nx1 * nx2)
      IF (i3 < i3min .OR. i3 >= i3max .OR. i3 >= n3) THEN
        CYCLE
      END IF
      !
      idx = idx - (nx1 * nx2) * i3
      i2  = idx / nx1
      IF (i2 >= n2) THEN
        CYCLE
      END IF
      !
      idx = idx - nx1 * i2
      i1  = idx
      IF (i1 >= n1) THEN
        CYCLE
      END IF
      !
      IF (i3 >= (n3 - (n3 / 2))) THEN
        irz = i3 - n3
      ELSE
        irz = i3
      END IF
      irz = irz + (n3 / 2)
      irz = irz + irz0
      !
      IF (irz < irz_sta .OR. irz_end < irz) THEN
        gr(ir, 1:nsite) = 1.0_DP
      END IF
      !
    END DO
    !
  END SUBROUTINE truncate_gr
  !
  !---------------------------------------------------------------------------
  SUBROUTINE renormalize_gr(gr, gr1, gr2, irz_sta, irz_end)
    !---------------------------------------------------------------------------
    !
    ! ... truncate edge of g(r).
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(OUT) :: gr(:)
    REAL(DP), INTENT(IN)  :: gr1
    REAL(DP), INTENT(IN)  :: gr2
    INTEGER,  INTENT(IN)  :: irz_sta
    INTEGER,  INTENT(IN)  :: irz_end
    !
    INTEGER :: ir
    INTEGER :: idx
    INTEGER :: idx0
    INTEGER :: i3min
    INTEGER :: i3max
    INTEGER :: i1, i2, i3
    INTEGER :: n1, n2, n3
    INTEGER :: nx1, nx2, nx3
    INTEGER :: nnr
    INTEGER :: irz
    INTEGER :: irz0
    INTEGER :: irz_edg1
    INTEGER :: irz_edg2
    !
    ! ... FFT box
    n1  = rismt%cfft%dfftt%nr1
    n2  = rismt%cfft%dfftt%nr2
    n3  = rismt%cfft%dfftt%nr3
    nx1 = rismt%cfft%dfftt%nr1x
    nx2 = rismt%cfft%dfftt%nr2x
    nx3 = rismt%cfft%dfftt%nr3x
    nnr = rismt%cfft%dfftt%nnr
    !
    ! ... calculate planar average
    idx0     = nx1 * nx2 * rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
    i3min    = rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
    i3max    = rismt%cfft%dfftt%npp(rismt%cfft%dfftt%mype + 1) + i3min
    irz0     = rismt%lfft%izcell_start
    irz_edg1 = rismt%lfft%izright_gedge
    irz_edg2 = rismt%lfft%izleft_gedge
    !
    DO ir = 1, nnr
      !
      idx = idx0 + ir - 1
      i3  = idx / (nx1 * nx2)
      IF (i3 < i3min .OR. i3 >= i3max .OR. i3 >= n3) THEN
        CYCLE
      END IF
      !
      idx = idx - (nx1 * nx2) * i3
      i2  = idx / nx1
      IF (i2 >= n2) THEN
        CYCLE
      END IF
      !
      idx = idx - nx1 * i2
      i1  = idx
      IF (i1 >= n1) THEN
        CYCLE
      END IF
      !
      IF (i3 >= (n3 - (n3 / 2))) THEN
        irz = i3 - n3
      ELSE
        irz = i3
      END IF
      irz = irz + (n3 / 2)
      irz = irz + irz0
      !
      IF (irz_sta <= irz .AND. irz <= irz_edg2) THEN
        gr(ir) = gr(ir) - gr2
      ELSE IF (irz_edg1 <= irz .AND. irz <= irz_end) THEN
        gr(ir) = gr(ir) - gr1
      END IF
      !
    END DO
    !
  END SUBROUTINE renormalize_gr
  !
END SUBROUTINE normalize_lauerism
