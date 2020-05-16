!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE guess_3drism(rismt, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... create initial guess of 3D-RISM or Laue-RISM.
  !
  USE cell_base, ONLY : at, alat
  USE constants, ONLY : eps4, K_BOLTZMANN_RY
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_max
  USE rism,      ONLY : rism_type, ITYPE_3DRISM, ITYPE_LAUERISM
  USE solvmol,   ONLY : get_nuniq_in_solVs, solVs, &
                      & iuniq_to_isite, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER  :: nq
  INTEGER  :: iq
  INTEGER  :: iiq
  INTEGER  :: iv
  INTEGER  :: isolV
  INTEGER  :: iatom
  INTEGER  :: ir
  INTEGER  :: idx
  INTEGER  :: idx0
  INTEGER  :: i3min
  INTEGER  :: i3max
  INTEGER  :: i1, i2, i3
  LOGICAL  :: laue
  REAL(DP) :: beta
  REAL(DP) :: qv
  REAL(DP) :: vlj
  REAL(DP) :: cs0
  REAL(DP) :: csmax
  REAL(DP) :: erf0
  !
  REAL(DP), EXTERNAL  :: qe_erf
  !
  REAL(DP), PARAMETER :: VLJ_MAX  = eps4 ! Ry
  REAL(DP), PARAMETER :: CS_SCALE = 0.1_DP
  !
  ! ... number of sites in solvents
  nq = get_nuniq_in_solVs()
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_3DRISM .AND. rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%mp_site%nsite < nq) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nr < rismt%cfft%dfftt%nnr) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... if no data, return as normally done
  IF (rismt%nsite < 1) THEN
    GOTO 1
  END IF
  !
  ! ... Laue-RISM or not
  laue = .FALSE.
  IF (rismt%itype == ITYPE_LAUERISM) THEN
    laue = .TRUE.
  END IF
  !
  ! ... beta = 1 / (kB * T)
  beta = 1.0_DP / K_BOLTZMANN_RY / rismt%temp
  !
  ! ... create guess for each solvent's site
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    isolV = isite_to_isolV(iv)
    iatom = isite_to_iatom(iv)
    qv    = solVs(isolV)%charge(iatom)
    !
    idx0 = rismt%cfft%dfftt%nr1x * rismt%cfft%dfftt%nr2x &
       & * rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
    !
    i3min = rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
    i3max = rismt%cfft%dfftt%npp(rismt%cfft%dfftt%mype + 1) + i3min
    !
    ! ... set csr to be zero
    csmax = 0.0_DP
    rismt%csr(:, iiq) = 0.0_DP
    !
    ! ... set csr initially
    DO ir = 1, rismt%cfft%dfftt%nnr
      !
      idx = idx0 + ir - 1
      i3  = idx / (rismt%cfft%dfftt%nr1x * rismt%cfft%dfftt%nr2x)
      IF (i3 < i3min .OR. i3 >= i3max .OR. i3 >= rismt%cfft%dfftt%nr3) THEN
        CYCLE
      END IF
      !
      idx = idx - (rismt%cfft%dfftt%nr1x * rismt%cfft%dfftt%nr2x) * i3
      i2  = idx / rismt%cfft%dfftt%nr1x
      IF (i2 >= rismt%cfft%dfftt%nr2) THEN
        CYCLE
      END IF
      !
      idx = idx - rismt%cfft%dfftt%nr1x * i2
      i1  = idx
      IF (i1 >= rismt%cfft%dfftt%nr1) THEN
        CYCLE
      END IF
      !
      vlj = rismt%uljr(ir, iiq)
      IF (laue) THEN  ! add LJ-wall
        vlj =  vlj + rismt%uwr(ir, iiq)
      END IF
      !
      IF (vlj >= VLJ_MAX) THEN
        cs0 = beta * qv * rismt%vlr(ir)
        csmax = MAX(csmax, ABS(cs0))
        rismt%csr(ir, iiq) = cs0
      END IF
      !
    END DO
    !
    CALL mp_max(csmax, rismt%mp_site%intra_sitg_comm)
    !
    ! ... correct csr to be smooth
    DO ir = 1, rismt%cfft%dfftt%nnr
      !
      idx = idx0 + ir - 1
      i3  = idx / (rismt%cfft%dfftt%nr1x * rismt%cfft%dfftt%nr2x)
      IF (i3 < i3min .OR. i3 >= i3max .OR. i3 >= rismt%cfft%dfftt%nr3) THEN
        CYCLE
      END IF
      !
      idx = idx - (rismt%cfft%dfftt%nr1x * rismt%cfft%dfftt%nr2x) * i3
      i2  = idx / rismt%cfft%dfftt%nr1x
      IF (i2 >= rismt%cfft%dfftt%nr2) THEN
        CYCLE
      END IF
      !
      idx = idx - rismt%cfft%dfftt%nr1x * i2
      i1  = idx
      IF (i1 >= rismt%cfft%dfftt%nr1) THEN
        CYCLE
      END IF
      !
      IF (csmax > 0.0_DP) THEN
        cs0  = rismt%csr(ir, iiq)
        erf0 = qe_erf(ABS(cs0) / (CS_SCALE * csmax))
        rismt%csr(ir, iiq) = cs0 * erf0 * erf0
      END IF
      !
    END DO
    !
  END DO
  !
  ! ... correction for Laue-RISM
  IF (laue) THEN
    CALL correct_edge()
  END IF
  !
  ! ... normally done
1 CONTINUE
  ierr = IERR_RISM_NULL
  !
CONTAINS
  !
  SUBROUTINE correct_edge()
    IMPLICIT NONE
    !
    INTEGER  :: iz
    REAL(DP) :: z
    REAL(DP) :: z0
    !
    REAL(DP), PARAMETER :: Z_SCALE = 5.0_DP ! bohr
    !
    IF (rismt%nsite < 1) THEN
      RETURN
    END IF
    !
    z0 = 0.5_DP * at(3, 3)
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
        CYCLE
      END IF
      !
      idx = idx - (rismt%cfft%dfftt%nr1x * rismt%cfft%dfftt%nr2x) * i3
      i2  = idx / rismt%cfft%dfftt%nr1x
      IF (i2 >= rismt%cfft%dfftt%nr2) THEN
        CYCLE
      END IF
      !
      idx = idx - rismt%cfft%dfftt%nr1x * i2
      i1  = idx
      IF (i1 >= rismt%cfft%dfftt%nr1) THEN
        CYCLE
      END IF
      !
      IF (i3 < (rismt%cfft%dfftt%nr3 - (rismt%cfft%dfftt%nr3 / 2))) THEN
        iz = i3 + (rismt%cfft%dfftt%nr3 / 2)
      ELSE
        iz = i3 - rismt%cfft%dfftt%nr3 + (rismt%cfft%dfftt%nr3 / 2)
      END IF
      iz = iz + rismt%lfft%izcell_start
      !
      z = rismt%lfft%zleft + rismt%lfft%zoffset + rismt%lfft%zstep * DBLE(iz - 1)
      !
      IF (rismt%lfft%xright) THEN
        erf0 = qe_erf(alat * (z0 - z) / Z_SCALE)
        rismt%csr(ir, :) = rismt%csr(ir, :) * (erf0 * erf0)
      END IF
      !
      IF (rismt%lfft%xleft) THEN
        erf0 = qe_erf(alat * (z + z0) / Z_SCALE)
        rismt%csr(ir, :) = rismt%csr(ir, :) * (erf0 * erf0)
      END IF
      !
    END DO
    !
  END SUBROUTINE correct_edge
  !
END SUBROUTINE guess_3drism
