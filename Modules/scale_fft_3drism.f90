!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE scale_fft_3drism(fc, at_old, laue)
  !---------------------------------------------------------------------------
  !
  ! ... scale G-vectors of 3D-RISM, for Variable Cell
  !
  USE cell_base,  ONLY : bg
  USE fft_custom, ONLY : fft_cus, gshells_custom
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  !
  TYPE(fft_cus), INTENT(INOUT) :: fc
  REAL(DP),      INTENT(IN)    :: at_old(3, 3)
  LOGICAL,       INTENT(IN)    :: laue
  !
  INTEGER  :: ig
  REAL(DP) :: gx
  REAL(DP) :: gy
  REAL(DP) :: gz
  !
  IF (.NOT. fc%initialized) THEN
    CALL errore(' rescale_fft_3drism ', ' fc is not initialized ', 1)
    RETURN
  END IF
  !
  ! ... scale G-vectors
  CALL cryst_to_cart(fc%ngmt, fc%gt, at_old, -1)
  CALL cryst_to_cart(fc%ngmt, fc%gt, bg,     +1)
  !
  ! ... update G^2
  DO ig = 1, fc%ngmt
    gx = fc%gt(1, ig)
    gy = fc%gt(2, ig)
    gz = fc%gt(3, ig)
    fc%ggt(ig) = gx * gx + gy * gy + gz * gz
  END DO
  !
  ! ... update G-shell
  IF (.NOT. laue) THEN
    CALL gshells_custom(fc, .FALSE.) ! lmovecell must be .FALSE.
  END IF
  !
END SUBROUTINE scale_fft_3drism
