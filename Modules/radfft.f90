!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE radfft
  !--------------------------------------------------------------------------
  !
  ! ... this module perform radial FFT, which is defined as
  ! ...                       / inf
  ! ...   g * A(g) = 4 * pi * | dr sin(g * r) * r * A(r) ,
  ! ...                       / 0
  ! ... and
  ! ...                  1        / inf
  ! ...   r * A(r) = ---------- * | dg sin(g * r) * g * A(g) .
  ! ...               2 * pi^2    / 0
  !
  USE constants,   ONLY : tpi
  USE fft_scalar,  ONLY : cft_1z
  USE fft_support, ONLY : good_fft_order
  USE kinds,       ONLY : DP
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... define variables
  TYPE radfft_type
    INTEGER           :: ngrid     ! number of grids
    INTEGER           :: mgrid     ! number of grids for FFT-box
    INTEGER           :: lgrid     ! modified mgrid
    REAL(DP), POINTER :: rgrid(:)  ! grids in R-space
    REAL(DP), POINTER :: ggrid(:)  ! grids in G-space
  END TYPE radfft_type
  !
  ! ... public components
  PUBLIC :: radfft_type
  PUBLIC :: allocate_radfft
  PUBLIC :: deallocate_radfft
  PUBLIC :: fw_radfft
  PUBLIC :: inv_radfft
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE allocate_radfft(radfft0, nr, rmax)
    !--------------------------------------------------------------------------
    !
    ! ... initialize radfft_type
    !
    IMPLICIT NONE
    !
    TYPE(radfft_type),  INTENT(INOUT) :: radfft0
    INTEGER,            INTENT(IN)    :: nr
    REAL(DP),           INTENT(IN)    :: rmax
    !
    INTEGER  :: igrid
    INTEGER  :: ngrid
    REAL(DP) :: hr
    REAL(DP) :: hg
    !
    ! check nr
    IF (nr < 2) THEN
      CALL errore(' allocate_radfft ', ' too small number of grids ', 1)
    END IF
    !
    ! number of grids
    ngrid         = nr
    radfft0%ngrid = ngrid
    radfft0%mgrid = 2 * ngrid - 1
    radfft0%lgrid = good_fft_order(2 * ngrid - 1)
    !
    ! R-space
    ALLOCATE(radfft0%rgrid(ngrid))
    hr = rmax / DBLE(ngrid)
    DO igrid = 1, ngrid
      radfft0%rgrid(igrid) = hr * DBLE(igrid - 1)
    END DO
    !
    ! G-space
    ALLOCATE(radfft0%ggrid(ngrid))
    hg = (tpi / rmax) * (DBLE(ngrid) / DBLE(2 * ngrid - 1))
    DO igrid = 1, ngrid
      radfft0%ggrid(igrid) = hg * DBLE(igrid - 1)
    END DO
    !
  END SUBROUTINE allocate_radfft
  !
  !--------------------------------------------------------------------------
  SUBROUTINE deallocate_radfft(radfft0)
    !--------------------------------------------------------------------------
    !
    ! ... finalize radfft_type
    !
    IMPLICIT NONE
    !
    TYPE(radfft_type), INTENT(INOUT) :: radfft0
    !
    radfft0%ngrid = 0
    radfft0%mgrid = 0
    radfft0%lgrid = 0
    IF (ASSOCIATED(radfft0%rgrid)) DEALLOCATE(radfft0%rgrid)
    IF (ASSOCIATED(radfft0%ggrid)) DEALLOCATE(radfft0%ggrid)
    !
  END SUBROUTINE deallocate_radfft
  !
  !--------------------------------------------------------------------------
  SUBROUTINE fw_radfft(radfft0, cr, cg)
    !--------------------------------------------------------------------------
    !
    ! ... FFT R -> G
    !
    IMPLICIT NONE
    !
    TYPE(radfft_type), INTENT(IN)  :: radfft0
    REAL(DP),          INTENT(IN)  :: cr(:)
    REAL(DP),          INTENT(OUT) :: cg(:)
    !
    INTEGER                  :: igrid
    REAL(DP)                 :: dr
    REAL(DP)                 :: fac
    COMPLEX(DP), ALLOCATABLE :: crr(:)
    COMPLEX(DP), ALLOCATABLE :: cgg(:)
    !
    ! allocate memory
    ALLOCATE(crr(radfft0%lgrid))
    ALLOCATE(cgg(radfft0%lgrid))
    !
    ! cr -> crr (ungerade)
    dr  = radfft0%rgrid(2) - radfft0%rgrid(1)
    fac = dr * tpi
    DO igrid = 1, radfft0%ngrid
      crr(igrid) = CMPLX(0.0_DP, fac * radfft0%rgrid(igrid) * cr(igrid), kind=DP)
    END DO
    DO igrid = (radfft0%ngrid + 1), radfft0%mgrid
      crr(igrid) = -crr(2 * radfft0%ngrid - igrid + 1)
    END DO
    !
    ! 1D-FFT (crr -> cgg)
    CALL cft_1z(crr, 1, radfft0%mgrid, radfft0%lgrid, -1, cgg)
    !
    ! cgg -> cg (only sin)
    cg(1) = 0.0_DP
    DO igrid = 2, radfft0%ngrid
      cg(igrid) = DBLE(cgg(igrid)) / radfft0%ggrid(igrid) * DBLE(radfft0%mgrid)
    END DO
    !
    ! deallocate memory
    DEALLOCATE(crr)
    DEALLOCATE(cgg)
    !
  END SUBROUTINE fw_radfft
  !
  !--------------------------------------------------------------------------
  SUBROUTINE inv_radfft(radfft0, cg, cr)
    !--------------------------------------------------------------------------
    !
    ! ... FFT G -> R
    !
    IMPLICIT NONE
    !
    TYPE(radfft_type), INTENT(IN)  :: radfft0
    REAL(DP),          INTENT(IN)  :: cg(:)
    REAL(DP),          INTENT(OUT) :: cr(:)
    !
    INTEGER                  :: igrid
    REAL(DP)                 :: dg
    REAL(DP)                 :: fac
    COMPLEX(DP), ALLOCATABLE :: cgg(:)
    COMPLEX(DP), ALLOCATABLE :: crr(:)
    !
    ! allocate memory
    ALLOCATE(cgg(radfft0%lgrid))
    ALLOCATE(crr(radfft0%lgrid))
    !
    ! cg -> cgg (only sin)
    dg  = radfft0%ggrid(2) - radfft0%ggrid(1)
    fac = -dg / tpi / tpi
    DO igrid = 1, radfft0%ngrid
      cgg(igrid) = CMPLX(0.0_DP, fac * radfft0%ggrid(igrid) * cg(igrid), kind=DP)
    END DO
    DO igrid = (radfft0%ngrid + 1), radfft0%mgrid
      cgg(igrid) = -cgg(2 * radfft0%ngrid - igrid + 1)
    END DO
    !
    ! 1D-FFT (cgg -> crr)
    CALL cft_1z(cgg, 1, radfft0%mgrid, radfft0%lgrid, +1, crr)
    !
    ! crr -> cr (ungerade)
    cr(1) = 0.0_DP
    DO igrid = 2, radfft0%ngrid
      cr(igrid) = DBLE(crr(igrid)) / radfft0%rgrid(igrid)
    END DO
    !
    ! deallocate memory
    DEALLOCATE(cgg)
    DEALLOCATE(crr)
    !
  END SUBROUTINE inv_radfft
  !
END MODULE radfft
