!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE closure(rismt, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... solve a Closure Equation for RISM
  !
  USE constants, ONLY : K_BOLTZMANN_RY
  USE kinds,     ONLY : DP
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE rism,      ONLY : rism_type, ITYPE_1DRISM, CLOSURE_HNC, CLOSURE_KH
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER  :: iclosure
  INTEGER  :: mr
  REAL(DP) :: beta
  !
  ! ... check data type
  IF (rismt%itype == ITYPE_1DRISM .AND. rismt%nr /= rismt%ng) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... set variables
  iclosure = rismt%closure
  mr = rismt%nr * rismt%nsite
  beta = 1.0_DP / K_BOLTZMANN_RY / rismt%temp
  !
  ! ... if no data, return as normally done.
  IF (mr < 1) THEN
    GOTO 100
  END IF
  !
  ! ... solve closure equation
  IF (iclosure == CLOSURE_HNC) THEN
    !
    CALL closure_HNC_x(mr, beta, &
       & rismt%usr(1, 1), rismt%hr(1, 1), rismt%csr(1, 1), rismt%gr(1, 1))
    !
  ELSE IF (iclosure == CLOSURE_KH) THEN
    !
    CALL closure_KH_x(mr, beta, &
       & rismt%usr(1, 1), rismt%hr(1, 1), rismt%csr(1, 1), rismt%gr(1, 1))
    !
  ELSE
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... set zero, if R = 0 (for 1D-RISM)
  IF (rismt%itype == ITYPE_1DRISM) THEN
    IF (rismt%mp_task%ivec_start == 1) THEN
      rismt%gr(1, :) = 0.0_DP
    END IF
  END IF
  !
  ! ... normally done
100 CONTINUE
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE closure
!
!---------------------------------------------------------------------------
SUBROUTINE closure_HNC_x(nr, beta, ur, hr, cr, gr)
  !---------------------------------------------------------------------------
  !
  ! ... HyperNetted-Chain model
  ! ... (J.P.Hansen et al., Theory of simple liquids. Academic Press, London, 1990.)
  ! ...
  ! ...   g(r) = exp(t(r))
  ! ...   t(r) = -beta * u(r) + h(r) - c(r)
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: nr
  REAL(DP), INTENT(IN)  :: beta
  REAL(DP), INTENT(IN)  :: ur(1:*)
  REAL(DP), INTENT(IN)  :: hr(1:*)
  REAL(DP), INTENT(IN)  :: cr(1:*)
  REAL(DP), INTENT(OUT) :: gr(1:*)
  !
  INTEGER  :: ir
  REAL(DP) :: tr
  !
  REAL(DP), PARAMETER :: MAX_EXP = 100.0_DP
  !
  DO ir = 1, nr
    tr = -beta * ur(ir) + hr(ir) - cr(ir)
    gr(ir) = EXP(MIN(tr, MAX_EXP))
  END DO
  !
END SUBROUTINE closure_HNC_x
!
!---------------------------------------------------------------------------
SUBROUTINE closure_KH_x(nr, beta, ur, hr, cr, gr)
  !---------------------------------------------------------------------------
  !
  ! ... Kovalenko and Hirata's model
  ! ... (A.Kovalenko, F.Hirata, J. Chem. Phys. 1999, 110, 10095-10112)
  ! ...
  ! ...   g(r) = exp(t(r)), if t(r) < 0
  ! ...   g(r) = 1 + t(r) , if t(r) > 0
  ! ...   t(r) = -beta * u(r) + h(r) - c(r)
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: nr
  REAL(DP), INTENT(IN)  :: beta
  REAL(DP), INTENT(IN)  :: ur(1:*)
  REAL(DP), INTENT(IN)  :: hr(1:*)
  REAL(DP), INTENT(IN)  :: cr(1:*)
  REAL(DP), INTENT(OUT) :: gr(1:*)
  !
  INTEGER  :: ir
  REAL(DP) :: tr
  !
  DO ir = 1, nr
    tr = -beta * ur(ir) + hr(ir) - cr(ir)
    IF (tr < 0.0_DP) THEN
      gr(ir) = EXP(tr)
    ELSE
      gr(ir) = 1.0_DP + tr
    END IF
  END DO
  !
END SUBROUTINE closure_KH_x

