!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE upf_spinorb
  !
  !! Variables needed for calculations with spin-orbit
  !
  USE upf_kinds,   ONLY : DP
  USE upf_params,  ONLY : lmaxx, lqmax 
  !
  !! FIXME: rot_ylm could be dynamically allocated
  !
  IMPLICIT NONE
  SAVE

  LOGICAL :: lspinorb
  !! if .TRUE. this is a spin-orbit calculation
  COMPLEX (DP) :: rot_ylm(lqmax,lqmax)
  !! transform real spherical harmonics into complex ones
  COMPLEX (DP), ALLOCATABLE :: fcoef(:,:,:,:,:)
  !! function needed to account for spinors.
  !
  ! GPU vars
  COMPLEX(DP), ALLOCATABLE :: fcoef_d(:,:,:,:,:)
#if defined(__CUDA)
  attributes (DEVICE) :: fcoef_d
#endif

CONTAINS

  SUBROUTINE deallocate_spinorb
     IMPLICIT NONE
     IF( ALLOCATED( fcoef ) )   DEALLOCATE( fcoef )
     IF( ALLOCATED( fcoef_d ) ) DEALLOCATE( fcoef_d )
  END SUBROUTINE deallocate_spinorb

END MODULE upf_spinorb

