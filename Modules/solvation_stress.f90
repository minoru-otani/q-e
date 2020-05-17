!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE solvation_stress(rismt, sigma, rhog_ele, vloc, dvloc, ierr)
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,     ONLY : DP
  USE solute,    ONLY : get_solU_LJ_stress
  USE rism,      ONLY : rism_type, ITYPE_3DRISM, ITYPE_LAUERISM
  USE gvect,     ONLY : ngl
  USE ions_base, ONLY : nsp
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  REAL(DP),        INTENT(OUT) :: sigma(3, 3)
  COMPLEX(DP),     INTENT(IN)  :: rhog_ele(*)
  REAL(DP),        INTENT(IN)  :: vloc(ngl, nsp)
  REAL(DP),        INTENT(IN)  :: dvloc(ngl, nsp)
  INTEGER,         INTENT(OUT) :: ierr
  !
  REAL(DP) :: sigmaion(3,3)
  REAL(DP) :: sigmaionlong(3,3)
  REAL(DP) :: sigmalj (3,3)
  REAL(DP) :: sigmahar(3,3)
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_3DRISM .AND. rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF

  IF( rismt%itype == ITYPE_3DRISM ) THEN
    !
    ! ... Lennad-Jones term
    sigmalj = 0.0_DP
    CALL get_solU_LJ_stress( rismt, sigmalj, ierr )
    IF (ierr /= IERR_RISM_NULL) RETURN

    !
    ! ... electron Coulomb term
    sigmahar = 0.0_DP
    call solvation_stress_har( rismt, sigmahar, rhog_ele, ierr )
    IF (ierr /= IERR_RISM_NULL) RETURN

    !
    ! ... ions term (short and long part)
    sigmaion = 0.0_DP
    call solvation_stress_ion( rismt, sigmaion, vloc, dvloc, ierr )
    IF (ierr /= IERR_RISM_NULL) RETURN
  END IF

  IF( rismt%itype == ITYPE_LAUERISM ) THEN
    !
    ! ... Lennad-Jones term
    sigmalj = 0.0_DP
    CALL get_solU_LJ_stress( rismt, sigmalj, ierr )
    IF (ierr /= IERR_RISM_NULL) RETURN

    !
    ! ... electron Coulomb term
    sigmahar = 0.0_DP
    call solvation_esm_stress_har( rismt, sigmahar, rhog_ele, ierr )
    IF (ierr /= IERR_RISM_NULL) RETURN

    !
    ! ... ions term (short part)
    sigmaion = 0.0_DP
    call solvation_esm_stress_locshort( rismt, sigmaion, vloc, dvloc, ierr )
    IF (ierr /= IERR_RISM_NULL) RETURN

    !
    ! ... ions term (long part)
    sigmaionlong = 0.0_DP
    call solvation_esm_stress_loclong( rismt, sigmaionlong, ierr )
    IF (ierr /= IERR_RISM_NULL) RETURN

    sigmaion(:,:) = sigmaion(:,:) + sigmaionlong(:,:)
  END IF
  !
  ! ... total solvation stress
  sigma(:,:) = sigmalj(:,:) + sigmahar(:,:) + sigmaion(:,:)

!!$  write(stdout,*) "stressLJ ", real(sigmalj(:,1))
!!$  write(stdout,*) "stressLJ ", real(sigmalj(:,2))
!!$  write(stdout,*) "stressLJ ", real(sigmalj(:,3))
!!$  write(stdout,*)
!!$
!!$  write(stdout,*) "stressHar", real(sigmahar(:,1))
!!$  write(stdout,*) "stressHar", real(sigmahar(:,2))
!!$  write(stdout,*) "stressHar", real(sigmahar(:,3))
!!$  write(stdout,*)
!!$
!!$  write(stdout,*) "stressIon", real(sigmaion(:,1))
!!$  write(stdout,*) "stressIon", real(sigmaion(:,2))
!!$  write(stdout,*) "stressIon", real(sigmaion(:,3))
!!$  write(stdout,*)
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE solvation_stress
!
!---------------------------------------------------------------------------
SUBROUTINE solvation_stress_ion( rismt, sigma, vloc, dvloc, ierr )
  !---------------------------------------------------------------------------
  !
  USE err_rism,      ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE control_flags, ONLY : gamma_only
  USE kinds,         ONLY : DP
  USE rism,          ONLY : rism_type, ITYPE_3DRISM, ITYPE_LAUERISM
  use gvect,         ONLY : ngl, igtongl
  USE ions_base,     ONLY : nsp, nat, ityp, tau
  USE cell_base,     ONLY : omega, tpiba2
  USE constants,     ONLY : e2, fpi, tpi
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  REAL(DP),        INTENT(OUT) :: sigma(3, 3)
  REAL(DP),        INTENT(IN)  :: vloc(ngl, nsp)
  REAL(DP),        INTENT(IN)  :: dvloc(ngl, nsp)
  INTEGER,         INTENT(OUT) :: ierr
  !
  COMPLEX(DP), POINTER :: rhog_liq(:)
  real(DP) :: fact, evloc, evlocg, arg
  integer :: nt, is, ig, la, mu, ia
  COMPLEX(DP)  :: strf(rismt%cfft%ngmt)
#if defined(_OPENMP)
  REAL(DP) :: sgomp(3, 3)
#endif
  !
  ! ... set solvent density
  IF (rismt%itype == ITYPE_3DRISM ) THEN
    rhog_liq => rismt%rhog
  ELSE IF(rismt%itype == ITYPE_LAUERISM ) THEN
    rhog_liq => rismt%rhog_pbc
  ELSE
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... calculate stress
  sigma = 0.0_DP
  evloc = 0.0_DP

  if (gamma_only) then
    fact = 2.d0
  else
    fact = 1.d0
  end if

  do nt = 1, nsp ! element
    ! calculate structure factor of this element
    strf(:) = (0.d0,0.d0)
    do ia = 1, nat ! atom of this element
      if( ityp(ia) == nsp ) then
!$omp parallel do default(shared) private(ig,arg)
        do ig = 1, rismt%cfft%ngmt
          arg = sum(rismt%cfft%gt(:,ig)*tau(:,ia)) * tpi
          strf(ig) = strf(ig) + CMPLX(cos (arg), -sin (arg),kind=DP)
        enddo
      endif
    enddo

    ! calculate local energy
    if (rismt%cfft%gstart_t==2) then
      evloc = evloc + &
        rhog_liq(1) * strf(1) * vloc (igtongl (1), nt)
    end if

    evlocg=0.0_DP
!$omp parallel do default(shared) private(ig) reduction(+:evlocg)
    do ig = rismt%cfft%gstart_t, rismt%cfft%ngmt
      evlocg = evlocg &
        + DBLE(CONJG(rhog_liq(ig))*strf(ig)) &
        * vloc(igtongl(ig),nt) * fact
    enddo
    evloc = evloc + evlocg

    ! no G=0 contribution
!$omp parallel default(shared) private(ig,la,mu,sgomp)
#if defined(_OPENMP)
    sgomp = 0.0_DP
#endif
!$omp do
    do ig = 1, rismt%cfft%ngmt
      do la=1, 3
        do mu=1, la
#if defined(__OPENMP)
          sgomp(la,mu) = sgomp(la,mu) &
#else
          sigma(la,mu) = sigma(la,mu) &
#endif
            + DBLE(CONJG(rhog_liq(ig))*strf(ig)) &
            * 2.0d0 * dvloc(igtongl(ig),nt) &
            * tpiba2 * rismt%cfft%gt(la,ig) * rismt%cfft%gt(mu,ig) * fact
        enddo
      enddo
    enddo
!$omp end do

#if defined(__OPENMP)
!$omp critical
    sigma(1:3,1:3) = sigma(1:3,1:3) + sgomp(1:3,1:3)
!$omp end critical
#endif
!$omp end parallel
  enddo ! nt
  !
  CALL mp_sum( evloc, rismt%mp_site%inter_sitg_comm)
  CALL mp_sum( evloc, rismt%mp_site%intra_sitg_comm)
  call mp_sum( sigma, rismt%mp_site%inter_sitg_comm)
  call mp_sum( sigma, rismt%mp_site%intra_sitg_comm)
  !
  do la=1, 3
    sigma(la,la) = sigma(la,la) + evloc
    do mu=1, la-1
      sigma(mu,la) = sigma(la,mu)
    enddo
  enddo
  !
  ! ... clear solvent density
  NULLIFY(rhog_liq)
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE solvation_stress_ion
!
!---------------------------------------------------------------------------
SUBROUTINE solvation_stress_har( rismt, sigma, rhog_ele, ierr )
  !---------------------------------------------------------------------------
  !
  USE err_rism,      ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE control_flags, ONLY : gamma_only
  USE kinds,         ONLY : DP
  USE rism,          ONLY : rism_type, ITYPE_3DRISM, ITYPE_LAUERISM
  USE cell_base,     ONLY : omega, tpiba2
  USE constants,     ONLY : e2, fpi
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  REAL(DP),        INTENT(OUT) :: sigma(3, 3)
  COMPLEX(DP),     INTENT(IN)  :: rhog_ele(rismt%cfft%ngmt)
  INTEGER,         INTENT(OUT) :: ierr
  !
  COMPLEX(DP), POINTER :: rhog_liq(:)
  real(DP) :: evhar, svhar, g2
  integer :: is, ig, la, mu
#if defined(_OPENMP)
  REAL(DP) :: sgomp(3, 3)
#endif
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_3DRISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... set solvent density
  rhog_liq => rismt%rhog
  !
  ! ... calculate stress
  sigma = 0.0_DP
  evhar = 0.0_DP
!$omp parallel default(shared) private(ig,g2,svhar,la,mu,sgomp) reduction(+:evhar)
#if defined(_OPENMP)
  sgomp = 0.0_DP
#endif
!$omp do
  DO ig = rismt%cfft%gstart_t, rismt%cfft%ngmt
    g2 = rismt%cfft%ggt(ig) * tpiba2

    svhar = REAL( rhog_ele(ig) * CONJG(rhog_liq(ig)) ) / g2
    evhar = evhar + svhar

    do la=1, 3
      do mu=1, la
#if defined(__OPENMP)
          sgomp(la,mu) = sgomp(la,mu) &
#else
          sigma(la,mu) = sigma(la,mu) &
#endif
          + svhar * tpiba2 * 2 &
          * rismt%cfft%gt(la,ig) * rismt%cfft%gt(mu,ig) / g2
      enddo
    enddo
  enddo
!$omp end do

#if defined(__OPENMP)
!$omp critical
  sigma(1:3,1:3) = sigma(1:3,1:3) + sgomp(1:3,1:3)
!$omp end critical
#endif
!$omp end parallel
  !
  if (gamma_only) then
    evhar = evhar * 2.0d0
    sigma = sigma * 2.0d0
  end if

  evhar = evhar * (-omega *fpi*e2)
  sigma = sigma * (-fpi*e2)
  !
  CALL mp_sum( evhar, rismt%mp_site%inter_sitg_comm)
  CALL mp_sum( evhar, rismt%mp_site%intra_sitg_comm)
  call mp_sum( sigma, rismt%mp_site%inter_sitg_comm)
  call mp_sum( sigma, rismt%mp_site%intra_sitg_comm)
  !
  do la = 1, 3
    sigma(la,la) = sigma(la,la) - evhar/omega
    do mu=1, la-1
      sigma(mu,la) = sigma(la, mu)
    enddo
  enddo
  !
  ! ... clear solvent density
  NULLIFY(rhog_liq)
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE solvation_stress_har
