!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE allocate_fft_3drism(fc, ecutv, laue, mp_task)
  !---------------------------------------------------------------------------
  !
  ! ... initialize 3D-FFT for 3D-RISM
  !
  USE cell_base,        ONLY : tpiba2
  USE control_flags,    ONLY : gamma_only
  USE fft_custom,       ONLY : fft_cus, ggent, gshells_custom
  USE kinds,            ONLY : DP
  USE mp_rism,          ONLY : mp_rism_task
  USE gvect,            ONLY : ecutrho
  !
  IMPLICIT NONE
  !
  TYPE(fft_cus),      INTENT(INOUT) :: fc
  REAL(DP),           INTENT(IN)    :: ecutv
  LOGICAL,            INTENT(IN)    :: laue
  TYPE(mp_rism_task), INTENT(IN)    :: mp_task
  !
  IF (fc%initialized) THEN
    CALL errore(' allocate_fft_3drism ', ' fc is already initialized ', 1)
    RETURN
  END IF
  !
  fc%dual_t = ecutrho / ecutv
  fc%ecutt  = ecutv
  fc%gcutmt = ecutv / tpiba2
  !
  CALL data_structure_3drism(fc, gamma_only, mp_task)
  !
  fc%initialized = .TRUE.
  !
  CALL ggent(fc)
  !
  IF (.NOT. laue) THEN
    CALL gshells_custom(fc, .FALSE.) ! lmovecell must be .FALSE.
  END IF
  !
END SUBROUTINE allocate_fft_3drism
