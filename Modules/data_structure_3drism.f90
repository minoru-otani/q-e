!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE data_structure_3drism(fc, gamma_only, mp_task)
  !---------------------------------------------------------------------------
  !
  ! ... this routine sets the data structure for the 3D-RISM's fft array
  ! ... In the parallel case, it distributes columns to processes, too
  !
  USE cell_base,  ONLY : at, bg
  USE fft_base,   ONLY : dfftp, smap
  USE fft_custom, ONLY : fft_cus, gvec_init
  USE fft_types,  ONLY : fft_type_init
  USE gvect,      ONLY : gcutm
  USE kinds,      ONLY : DP
  USE mp_bands,   ONLY : ntask_groups
  USE mp_rism,    ONLY : mp_rism_task
  !
  IMPLICIT NONE
  !
  TYPE(fft_cus),      INTENT(INOUT) :: fc
  LOGICAL,            INTENT(IN)    :: gamma_only
  TYPE(mp_rism_task), INTENT(IN)    :: mp_task  ! must be same as intra_bgrp_comm
  !
  INTEGER :: ngs_
  INTEGER :: me
  INTEGER :: nproc
  INTEGER :: intra_comm
  INTEGER :: root
  INTEGER :: nogrp
#if defined (__MPI) && !defined (__USE_3D_FFT)
  LOGICAL :: lpara = .TRUE.
#else
  LOGICAL :: lpara = .FALSE.
#endif
  !
  me         = mp_task%me_task
  nproc      = mp_task%nproc_task
  intra_comm = mp_task%itask_comm
  root       = mp_task%root_task
  nogrp      = ntask_groups
  !
  ! ... set up fft descriptors, including parallel stuff: sticks, planes, etc.
  !
  CALL fft_type_init(fc%dfftt, smap, "rho", gamma_only, lpara, intra_comm, at, bg, fc%gcutmt, 1.0_DP)
  !
  ngs_ = fc%dfftt%ngl(fc%dfftt%mype + 1)
  IF (gamma_only) THEN
    ngs_ = (ngs_ + 1) / 2
  END IF
  !
  ! ... initialize local and global number of G-vectors
  !
  CALL gvec_init(fc, ngs_, intra_comm)
  !
END SUBROUTINE data_structure_3drism
