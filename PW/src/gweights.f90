!
! Copyright (C) 2001-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine gweights (nks, wk, nbnd, nelec, degauss, ngauss, &
     et, ef, demet, wg, is, isk, beta, eps)
  !--------------------------------------------------------------------
  !     calculates Ef and weights with the gaussian spreading technique
  ! ... Wrapper routine: computes first Ef, then the weights
  !
  USE kinds
  implicit none
  !
  integer, intent(in) :: nks, nbnd, ngauss, is, isk(nks)
  real(DP), intent(in) :: wk (nks), et (nbnd, nks), nelec, degauss
  ! wg must be (inout) and not (out) because if is/=0 only terms for
  ! spin=is are initialized; the remaining terms should be kept, not lost
  real(DP), intent(inout) :: wg (nbnd, nks)
  real(DP), intent(inout) :: ef
  real(DP), intent(out) :: demet
  real(DP), intent(in) :: beta, eps
  !
  real(DP) :: alpha
  real(DP) :: ef_by_n, ef_new
  real(DP), external :: efermig
  
  ! Calculate the Fermi energy ef

  ef_by_n = efermig (et, nbnd, nks, nelec, wk, degauss, ngauss, is, isk)

  if (beta > 0.0_DP) then
     ef_new = beta * ef + (1.0_DP - beta) * ef_by_n
     if (eps > 0.0_DP .and. abs(ef - ef_new) < eps) then
        alpha = (eps - abs(ef - ef_new)) / eps
        alpha = (alpha ** 3) * (10.0_DP - 15.0_DP * alpha + 6.0_DP * alpha ** 2)
        ef = alpha * ef + (1.0_DP - alpha) * ef_new
     else
        ef = ef_new
     end if
  else
     ef = ef_by_n
  end if

  ! Calculate weights

  CALL gweights_only (nks, wk, is, isk, nbnd, nelec, degauss, &
     ngauss, et, ef, demet, wg)

  return
end subroutine gweights
!
!--------------------------------------------------------------------
subroutine gweights_only (nks, wk, is, isk, nbnd, nelec, degauss, &
     ngauss, et, ef, demet, wg)
  !--------------------------------------------------------------------
  !     calculates weights with the gaussian spreading technique
  !     Fermi energy is provided in input
  !
  USE kinds
  implicit none
  !
  integer, intent(in) :: nks, nbnd, ngauss, is, isk(nks)
  real(DP), intent(in) :: wk (nks), et (nbnd, nks), nelec, degauss, ef
  ! wg must be (inout) and not (out) because if is/=0 only terms for
  ! spin=is are initialized; the remaining terms should be kept, not lost
  real(DP), intent(inout) :: wg (nbnd, nks)
  real(DP), intent(out) :: demet
  !
  integer :: kpoint, ibnd
  real(DP) , external :: wgauss, w1gauss

  demet = 0.d0
  do kpoint = 1, nks
     if (is /= 0) then
        if (isk(kpoint).ne.is) cycle
     end if
     do ibnd = 1, nbnd
        ! Calculate the gaussian weights
        wg (ibnd, kpoint) = wk (kpoint) * &
                            wgauss ( (ef-et(ibnd,kpoint)) / degauss, ngauss)
        !
        ! The correct (i.e. variational) form of the band energy is 
        !    Eband = \int e N(e) de   for e<Ef , where N(e) is the DOS
        ! This differs by the term "demet" from the sum of KS eigenvalues:
        !    Eks = \sum wg(n,k) et(n,k)
        ! which is non variational. When a Fermi-Dirac function is used
        ! for a given T, the variational energy is really the free energy F,
        ! and F = E - TS , with E = non variational energy, -TS = demet
        !
        demet = demet + wk (kpoint) * &
                 degauss * w1gauss ( (ef-et(ibnd,kpoint)) / degauss, ngauss)
     enddo

  enddo
  return
end subroutine gweights_only
