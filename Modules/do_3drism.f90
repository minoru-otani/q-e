!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE do_3drism(rismt, maxiter, rmsconv, nbox, eta, title, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... perform 3D-RISM method
  ! ...
  ! ... Variables:
  ! ...   rismt:   data structure
  ! ...   maxiter: maximum number of iterations
  ! ...   rmsconv: RMS of residual vectors to check convergence
  ! ...   nbox:    box size of MDIIS
  ! ...   eta:     step radius of MDIIS
  ! ...   title:   subtitle of calculation
  ! ...   ierr:    status of calculation
  !
  USE check_stop,     ONLY : check_stop_now, stopped_by_user
  USE control_flags,  ONLY : iverbosity, gamma_only
  USE err_rism,       ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE, IERR_RISM_NOT_CONVERGED
  USE fft_interfaces, ONLY : fwfft, invfft
  USE io_global,      ONLY : stdout
  USE kinds,          ONLY : DP
  USE mdiis,          ONLY : mdiis_type, allocate_mdiis, deallocate_mdiis, update_by_mdiis, reset_mdiis
  USE rism,           ONLY : rism_type, ITYPE_3DRISM
  !
  IMPLICIT NONE
  !
  TYPE(rism_type),  INTENT(INOUT) :: rismt
  INTEGER,          INTENT(IN)    :: maxiter
  REAL(DP),         INTENT(IN)    :: rmsconv
  INTEGER,          INTENT(IN)    :: nbox
  CHARACTER(LEN=*), INTENT(IN)    :: title
  REAL(DP),         INTENT(IN)    :: eta
  INTEGER,          INTENT(OUT)   :: ierr
  !
  INTEGER                  :: iter
  INTEGER                  :: ngrid
  LOGICAL                  :: lconv
  REAL(DP)                 :: rmscurr
  REAL(DP)                 :: rmssave
  INTEGER                  :: rmswarn
  REAL(DP),    ALLOCATABLE :: dcsr(:,:)
  COMPLEX(DP), ALLOCATABLE :: aux(:)
  TYPE(mdiis_type)         :: mdiist
  ! if mdiist is an automatic variable,
  ! pointers in mdiis_type may not work well.
  SAVE                     :: mdiist
  REAL(DP)                 :: csr_ (1, 1)
  REAL(DP)                 :: dcsr_(1, 1)
  !
  INTEGER,     PARAMETER   :: NPRINT      = 10
  INTEGER,     PARAMETER   :: MDIIS_EXT   = 3
  INTEGER,     PARAMETER   :: RMSWARN_MAX = 16
  REAL(DP),    PARAMETER   :: RMS_SMALL   = 0.95_DP
  REAL(DP),    PARAMETER   :: RMS_LARGE   = 2.00_DP
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_3DRISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nr < rismt%cfft%dfftt%nnr) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%ng < rismt%cfft%ngmt) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... allocate memory
  IF (rismt%nr * rismt%nsite > 0) THEN
    ALLOCATE(dcsr(rismt%nr, rismt%nsite))
  END IF
  IF (rismt%cfft%dfftt%nnr > 0) THEN
    ALLOCATE(aux(rismt%cfft%dfftt%nnr))
  END IF
  CALL allocate_mdiis(mdiist, nbox, rismt%nr * rismt%nsite, eta, MDIIS_EXT)
  !
  ! ... reset conditions
  lconv       = .FALSE.
  rismt%avail = .FALSE.
  rmssave     = 1.0E+99_DP
  rmswarn     = 0
  !
  ! ... start 3D-RISM iteration
  WRITE(stdout, '()')
  IF (LEN_TRIM(title) > 0) THEN
    WRITE(stdout, '(5X,"3D-RISM Calculation (",A,")")') TRIM(title)
  ELSE
    WRITE(stdout, '(5X,"3D-RISM Calculation")')
  END IF
  WRITE(stdout, '()')
  WRITE(stdout, '(5X,"convergence threshold    =",1PE10.3)') rmsconv
#ifdef __DEBUG_RISM
  !
  IF (iverbosity > 0) THEN
    CALL write_rism_type(rismt)
  END IF
#endif
  !
  DO iter = 1, maxiter
    !
    ! ... stop by user
    IF (check_stop_now()) THEN
      EXIT
    END IF
    !
    ! ... FFT: Cs(r) -> Cs(g)
    CALL fft_csr_to_csg()
    !
    ! ... 3D-RISM eq.: Cs(g) -> H(g)
    CALL eqn_3drism(rismt, ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
    ! ... FFT: H(g) -> H(r)
    CALL fft_hg_to_hr()
    !
    ! ... Closure: H(r) -> G(r)
    CALL closure(rismt, ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
    ! ... Residual: G(r) -> dCs(r)
    IF (rismt%nr * rismt%nsite > 0) THEN
      dcsr = rismt%gr - rismt%hr - 1.0_DP
    END IF
    !
    ! ... clean date out of physical range
    CALL clean_out_of_range()
    !
    ! ... calculate RMS
    ngrid = rismt%cfft%dfftt%nr1 * rismt%cfft%dfftt%nr2 * rismt%cfft%dfftt%nr3
    IF (rismt%nr * rismt%nsite > 0) THEN
      CALL rms_residual(ngrid * rismt%mp_site%nsite, rismt%nr * rismt%nsite, &
                      & dcsr, rmscurr, rismt%intra_comm)
    ELSE
      CALL rms_residual(ngrid * rismt%mp_site%nsite, rismt%nr * rismt%nsite, &
                      & dcsr_, rmscurr, rismt%intra_comm)
    END IF
    IF (rmscurr < rmsconv) THEN
      lconv = .TRUE.
    END IF
    !
    ! ... write data
    IF (iverbosity > 0 .OR. MOD(iter - 1, NPRINT) == 0 .OR. lconv .OR. iter == maxiter) THEN
      WRITE(stdout, '(5X,"iter. #",I6,"  RMS(g-h-1)=",1PE10.3,"  nbox=",I3)') &
      & iter, rmscurr, mdiist%nbox
      FLUSH(stdout)
    END IF
#ifdef __DEBUG_RISM
    !
    IF (iverbosity > 0) THEN
      CALL write_rism_type(rismt)
    END IF
#endif
    !
    ! ... converged ?
    IF (lconv) THEN
      EXIT
    END IF
    !
    ! ... check RMS to reset MDIIS
    IF (rmscurr > (RMS_SMALL * rmssave)) THEN
      rmswarn = rmswarn + 1
    ELSE
      rmswarn = 0
    END IF
    IF (rmswarn >= RMSWARN_MAX) THEN
      CALL reset_mdiis(mdiist, .TRUE.)
      rmssave = rmscurr
      rmswarn = 0
    END IF
    !
    IF (rmscurr > (RMS_LARGE * rmssave)) THEN
      CALL reset_mdiis(mdiist, .TRUE.)
      rmssave = rmscurr
      rmswarn = 0
    ELSE
      rmssave = MIN(rmscurr, rmssave)
    END IF
    !
    ! ... MDIIS: dCs(r) -> Cs(r)
    IF (rismt%nr * rismt%nsite > 0) THEN
      CALL update_by_mdiis(mdiist, rismt%csr, dcsr, rismt%intra_comm)
    ELSE
      CALL update_by_mdiis(mdiist, csr_, dcsr_, rismt%intra_comm)
    END IF
    !
  ! ... end 3D-RISM iteration
  END DO
  !
  WRITE(stdout, '()')
  WRITE(stdout, '(5X,"End of 3D-RISM calculation")')
  WRITE(stdout, '()')
  !
  ! ... iteration has been converged ?
  IF (lconv) THEN
    ! ... write convergence message
    WRITE(stdout, '()')
    WRITE(stdout, '(5X,"convergence has been achieved in ",I5," iterations")') iter
    WRITE(stdout, '()')
    FLUSH(stdout)
    !
    ! ... calculate chemical potential
    CALL chempot(rismt, ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
    ! ... calculate solvation's charge density, potential and energy
    CALL solvation_3drism(rismt, ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
    ! ... print chemical potential
    CALL print_chempot_3drism(rismt, ierr)
    IF (ierr /= IERR_RISM_NULL) THEN
      GOTO 100
    END IF
    !
    ! ... set conditions
    ierr = IERR_RISM_NULL
    rismt%avail = .TRUE.
    !
  ELSE
    ! ... write convergence message
    WRITE(stdout, '()')
    IF (stopped_by_user) THEN
      WRITE(stdout, '(5X,"convergence NOT achieved: stopped by user")')
    ELSE
      WRITE(stdout, '(5X,"convergence NOT achieved")')
    END IF
    WRITE(stdout, '()')
    FLUSH(stdout)
    !
    ! ... set conditions
    ierr = IERR_RISM_NOT_CONVERGED
    rismt%avail = .FALSE.
  END IF
#ifdef __DEBUG_RISM
  !
  IF (iverbosity > 0) THEN
    CALL write_rism_type(rismt)
  END IF
#endif
  !
  ! ... deallocate memory
100 CONTINUE
  IF (rismt%nr * rismt%nsite > 0) THEN
    DEALLOCATE(dcsr)
  END IF
  IF (rismt%cfft%dfftt%nnr > 0) THEN
    DEALLOCATE(aux)
  END IF
  CALL deallocate_mdiis(mdiist)
  !
CONTAINS
  !
  SUBROUTINE fft_csr_to_csg()
    IMPLICIT NONE
    INTEGER :: isite
    INTEGER :: jsite
    INTEGER :: ir
    INTEGER :: ig
    !
    DO isite = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      jsite = isite - rismt%mp_site%isite_start + 1
      IF (rismt%cfft%ngmt > 0) THEN
        rismt%csgz(:, jsite) = CMPLX(0.0_DP, 0.0_DP, kind=DP)
      END IF
      !
      DO ir = 1, rismt%cfft%dfftt%nnr
        aux(ir) = CMPLX(rismt%csr(ir, jsite), 0.0_DP, kind=DP)
      END DO
      IF (rismt%cfft%dfftt%nnr > 0) THEN
        CALL fwfft('Custom', aux, rismt%cfft%dfftt)
      END IF
      !
      DO ig = 1, rismt%cfft%ngmt
        rismt%csgz(ig, jsite) = aux(rismt%cfft%nlt(ig))
      END DO
    END DO
    !
  END SUBROUTINE fft_csr_to_csg
  !
  SUBROUTINE fft_hg_to_hr()
    IMPLICIT NONE
    INTEGER :: isite
    INTEGER :: jsite
    INTEGER :: ir
    INTEGER :: ig
    !
    DO isite = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      jsite = isite - rismt%mp_site%isite_start + 1
      IF (rismt%cfft%dfftt%nnr > 0) THEN
        rismt%hr(:, jsite) = 0.0_DP
      END IF
      !
      IF (rismt%cfft%dfftt%nnr > 0) THEN
        aux = CMPLX(0.0_DP, 0.0_DP, kind=DP)
      END IF
      DO ig = 1, rismt%cfft%ngmt
        aux(rismt%cfft%nlt(ig)) = rismt%hgz(ig, jsite)
      END DO
      IF (gamma_only) THEN
        DO ig = 1, rismt%cfft%ngmt
          aux(rismt%cfft%nltm(ig)) = CONJG(rismt%hgz(ig, jsite))
        END DO
      END IF
      IF (rismt%cfft%dfftt%nnr > 0) THEN
        CALL invfft('Custom', aux, rismt%cfft%dfftt)
      END IF
      !
      DO ir = 1, rismt%cfft%dfftt%nnr
        rismt%hr(ir, jsite) = DBLE(aux(ir))
      END DO
    END DO
    !
  END SUBROUTINE fft_hg_to_hr
  !
  SUBROUTINE clean_out_of_range()
    IMPLICIT NONE
    INTEGER :: ir
    INTEGER :: idx
    INTEGER :: idx0
    INTEGER :: i3min
    INTEGER :: i3max
    INTEGER :: i1, i2, i3
    !
    IF (rismt%nsite < 1) THEN
      RETURN
    END IF
    !
    idx0 = rismt%cfft%dfftt%nr1x * rismt%cfft%dfftt%nr2x &
       & * rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
    !
    i3min = rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
    i3max = rismt%cfft%dfftt%npp(rismt%cfft%dfftt%mype + 1) + i3min
    !
    DO ir = 1, rismt%cfft%dfftt%nnr
      !
      idx = idx0 + ir - 1
      i3  = idx / (rismt%cfft%dfftt%nr1x * rismt%cfft%dfftt%nr2x)
      IF (i3 < i3min .OR. i3 >= i3max .OR. i3 >= rismt%cfft%dfftt%nr3) THEN
        rismt%csr(ir, :) = 0.0_DP
        rismt%hr (ir, :) = 0.0_DP
        rismt%gr (ir, :) = 0.0_DP
        dcsr     (ir, :) = 0.0_DP
        CYCLE
      END IF
      !
      idx = idx - (rismt%cfft%dfftt%nr1x * rismt%cfft%dfftt%nr2x) * i3
      i2  = idx / rismt%cfft%dfftt%nr1x
      IF (i2 >= rismt%cfft%dfftt%nr2) THEN
        rismt%csr(ir, :) = 0.0_DP
        rismt%hr (ir, :) = 0.0_DP
        rismt%gr (ir, :) = 0.0_DP
        dcsr     (ir, :) = 0.0_DP
        CYCLE
      END IF
      !
      idx = idx - rismt%cfft%dfftt%nr1x * i2
      i1  = idx
      IF (i1 >= rismt%cfft%dfftt%nr1) THEN
        rismt%csr(ir, :) = 0.0_DP
        rismt%hr (ir, :) = 0.0_DP
        rismt%gr (ir, :) = 0.0_DP
        dcsr     (ir, :) = 0.0_DP
        CYCLE
      END IF
      !
    END DO
    !
  END SUBROUTINE clean_out_of_range
  !
END SUBROUTINE do_3drism
