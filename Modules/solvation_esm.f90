!
! Copyright (C) 2016 National Institute of Advanced Industrial Science and Technology (AIST)
! [ This code is written by Satomichi Nishihara. ]
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE solvation_esm_potential(rismt, iref, vref, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... calculate solvation potential of ESM(BC1), from Laue-RISM.
  ! ... calculation is performed around the expanded cell.
  !
  ! ... Variables:
  ! ...   iref: reference position of solvation potential
  ! ...   vref: reference value of solvation potential
  !
  USE cell_base,     ONLY : alat, tpiba, tpiba2
  USE constants,     ONLY : tpi, fpi, e2
  USE err_rism,      ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,         ONLY : DP
  USE lauefft,       ONLY : fw_lauefft_1z_exp, inv_lauefft_1z_exp
  USE mp,            ONLY : mp_sum
  USE rism,          ONLY : rism_type, ITYPE_LAUERISM
  USE rism3d_facade, ONLY : IREFERENCE_AVERAGE, IREFERENCE_RIGHT, IREFERENCE_LEFT
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  INTEGER,         INTENT(IN)    :: iref
  REAL(DP),        INTENT(OUT)   :: vref
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER                  :: iz
  INTEGER                  :: igz
  INTEGER                  :: igxy
  INTEGER                  :: jgxy
  REAL(DP)                 :: z
  REAL(DP)                 :: zright
  REAL(DP)                 :: zleft
  REAL(DP)                 :: zstart
  REAL(DP)                 :: dz
  REAL(DP)                 :: gz
  REAL(DP)                 :: gxy
  REAL(DP)                 :: ggxy
  REAL(DP)                 :: fac1  ! factors to    -> Energy/G/G
  REAL(DP)                 :: fac2  ! convert units -> Energy*R/G
  REAL(DP)                 :: fac3  !               -> Energy*R*R
  REAL(DP)                 :: phir
  REAL(DP)                 :: phil
  REAL(DP)                 :: rho0
  REAL(DP)                 :: realr
  REAL(DP)                 :: reall
  REAL(DP)                 :: imager
  REAL(DP)                 :: imagel
  REAL(DP)                 :: vsolv
  REAL(DP)                 :: vsolu
  COMPLEX(DP)              :: coeffr
  COMPLEX(DP)              :: coeffl
  COMPLEX(DP), ALLOCATABLE :: rhogt(:,:)
  COMPLEX(DP), ALLOCATABLE :: rhogz(:)
  COMPLEX(DP), ALLOCATABLE :: vpott(:,:)
  COMPLEX(DP), ALLOCATABLE :: expigzr(:)
  COMPLEX(DP), ALLOCATABLE :: expigzl(:)
  !
  COMPLEX(DP), PARAMETER   :: C_ZERO = CMPLX(0.0_DP, 0.0_DP, kind=DP)
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nrzl < rismt%lfft%nrz) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%ngxy < rismt%lfft%ngxy) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... allocate memory
  IF (rismt%lfft%ngz_x * rismt%lfft%ngxy > 0) THEN
    ALLOCATE(rhogt(  rismt%lfft%ngz_x, rismt%lfft%ngxy))
    ALLOCATE(vpott(  rismt%lfft%ngz_x, rismt%lfft%ngxy))
    rhogt = C_ZERO
    vpott = C_ZERO
  END IF
  IF (rismt%lfft%ngz_x > 0) THEN
    ALLOCATE(rhogz(  rismt%lfft%ngz_x))
    ALLOCATE(expigzr(rismt%lfft%ngz_x))
    ALLOCATE(expigzl(rismt%lfft%ngz_x))
    rhogz = C_ZERO
    expigzr = C_ZERO
    expigzl = C_ZERO
  END IF
  !
  ! ... set variables
  zright = rismt%lfft%zright
  zleft  = rismt%lfft%zleft
  zstart = rismt%lfft%zleft + rismt%lfft%zoffset
  dz     = rismt%lfft%zstep
  fac1   = e2 * fpi / tpiba2        ! -> Energy/G/G
  fac2   = e2 * fpi * alat / tpiba  ! -> Energy/R/G
  fac3   = e2 * fpi * alat * alat   ! -> Energy*R*R
  !
  ! ... initialize reference potential
  vref = 0.0_DP
  !
  ! ... calculate exp(i*gz*zright) and exp(i*gz*zleft)
!$omp parallel do default(shared) private(igz, gz, phir, phil)
  DO igz = 1, rismt%lfft%ngz_x
    gz   = rismt%lfft%gz_x(igz)
    phir = tpi * gz * zright
    phil = tpi * gz * zleft
    expigzr(igz) = CMPLX(COS(phir), SIN(phir), kind=DP)
    expigzl(igz) = CMPLX(COS(phil), SIN(phil), kind=DP)
  END DO
!$omp end parallel do
  !
  ! ... 1D-FFT of charge: Laue-rep. -> G-space
  IF (rismt%lfft%ngz_x * rismt%lfft%ngxy > 0) THEN
    rhogt = C_ZERO
    CALL fw_lauefft_1z_exp(rismt%lfft, rismt%rhog, rismt%nrzl, rhogt, rismt%lfft%ngz_x)
  END IF
  !
  ! ... Hartree potential: part of 4pi/G^2
  IF (rismt%lfft%ngz_x * rismt%lfft%ngxy > 0) THEN
    vpott = C_ZERO
  END IF
  !
  DO igxy = rismt%lfft%gxystart, rismt%lfft%ngxy
    ggxy = rismt%lfft%ggxy(igxy)
!$omp parallel do default(shared) private(igz, gz)
    DO igz = 1, rismt%lfft%ngz_x
      gz = rismt%lfft%gz_x(igz)
      vpott(igz, igxy) = rhogt(igz, igxy) * (fac1 / (gz * gz + ggxy))
    END DO
!$omp end parallel do
  END DO
  !
  IF (rismt%lfft%gxystart > 1) THEN
    igxy = 1
!$omp parallel do default(shared) private(igz, gz)
    DO igz = 1, rismt%lfft%ngz_x
      IF (igz /= rismt%lfft%gzzero_x) THEN
        gz = rismt%lfft%gz_x(igz)
        vpott(igz, igxy) = rhogt(igz, igxy) * (fac1 / (gz * gz))
      END IF
    END DO
!$omp end parallel do
  END IF
  !
  ! ... 1D-FFT of Hartree potential: G-space -> Laue-rep.
  IF (rismt%nrzl * rismt%ngxy > 0) THEN
    rismt%vpot = C_ZERO
  END IF
  IF (rismt%lfft%ngz_x * rismt%lfft%ngxy > 0) THEN
    CALL inv_lauefft_1z_exp(rismt%lfft, vpott, rismt%lfft%ngz_x, rismt%vpot, rismt%nrzl)
  END IF
  !
  ! ...
  ! ... Hartree potential, when Gxy /= 0
  ! ...
  DO igxy = rismt%lfft%gxystart, rismt%lfft%ngxy
    !
    jgxy = rismt%nrzl * (igxy - 1)
    gxy  = rismt%lfft%gnxy(igxy)
    !
    IF (rismt%lfft%ngz_x > 0) THEN
      rhogz(:) = rhogt(:, igxy)
    END IF
    !
    ! ... coeffr means
    !      ----           exp( i*gz*zright)
    !      >    rho(g) * -------------------
    !      ----              i*gz - gxy
    !       gz
    !
    ! ... coeffl means
    !      ----           exp( i*gz*zleft)
    !      >    rho(g) * -------------------
    !      ----              i*gz + gxy
    !       gz
    !
    coeffr = C_ZERO
    coeffl = C_ZERO
    !
!$omp parallel do default(shared) private(igz, gz) reduction(+:coeffr, coeffl)
    DO igz = 1, rismt%lfft%ngz_x
      gz = rismt%lfft%gz_x(igz)
      coeffr = coeffr + rhogz(igz) * expigzr(igz) / CMPLX(-gxy, gz, kind=DP)
      coeffl = coeffl + rhogz(igz) * expigzl(igz) / CMPLX( gxy, gz, kind=DP)
    END DO
!$omp end parallel do
    !
    ! ... when zleft <= z <= zright, potential is
    !
    !            exp( gxy*(z-zright))    ----           exp( i*gz*zright)
    !  (+4pi) * ---------------------- * >    rho(g) * -------------------
    !                 2 * gxy            ----              i*gz - gxy
    !                                     gz
    !
    !            exp(-gxy*(z+zleft))     ----           exp( i*gz*zleft)
    !  (-4pi) * ---------------------- * >    rho(g) * -------------------
    !                 2 * gxy            ----              i*gz + gxy
    !                                     gz
    !
!$omp parallel do default(shared) private(iz, z)
    DO iz = 1, rismt%lfft%nrz
      z = zstart + dz * DBLE(iz - 1)
      rismt%vpot(iz + jgxy) = rismt%vpot(iz + jgxy) + fac1 * ( &
      & + (0.5_DP / gxy) * EXP( tpi * gxy * (z - zright)) * coeffr &
      & - (0.5_DP / gxy) * EXP(-tpi * gxy * (z - zleft )) * coeffl )
    END DO
!$omp end parallel do
    !
  END DO
  !
  ! ...
  ! ... Hartree potential, when Gxy = 0
  ! ...
  IF (rismt%lfft%gxystart > 1) THEN
    !
    igxy = 1
    jgxy = rismt%nrzl * (igxy - 1)
    gxy  = rismt%lfft%gnxy(igxy)
    !
    IF (rismt%lfft%ngz_x > 0) THEN
      rhogz(:) = rhogt(:, igxy)
      rho0 = rhogz(rismt%lfft%gzzero_x)
    END IF
    !
    ! ... realr means
    !      ----  Re[ rho(gz) * exp( i*gz*zright) ]
    !      >    -----------------------------------
    !      ----                gz^2
    !      gz>0
    !
    ! ... reall means
    !      ----  Re[ rho(gz) * exp( i*gz*zleft) ]
    !      >    -----------------------------------
    !      ----                gz^2
    !      gz>0
    !
    ! ... imager means
    !      ----  Im[ rho(gz) * exp( i*gz*zright) ]
    !      >    -----------------------------------
    !      ----                gz
    !      gz>0
    !
    ! ... imagel means
    !      ----  Im[ rho(gz) * exp( i*gz*zleft) ]
    !      >    -----------------------------------
    !      ----                gz
    !      gz>0
    !
    realr  = 0.0_DP
    reall  = 0.0_DP
    imager = 0.0_DP
    imagel = 0.0_DP
    !
!$omp parallel do default(shared) private(igz, gz) reduction(+:realr, reall, imager, imagel)
    DO igz = (rismt%lfft%gzzero_x + 1), rismt%lfft%ngz_x
      gz = rismt%lfft%gz_x(igz)
      realr  = realr  + DBLE( rhogz(igz) * expigzr(igz)) / gz / gz
      reall  = reall  + DBLE( rhogz(igz) * expigzl(igz)) / gz / gz
      imager = imager + AIMAG(rhogz(igz) * expigzr(igz)) / gz
      imagel = imagel + AIMAG(rhogz(igz) * expigzl(igz)) / gz
    END DO
!$omp end parallel do
    !
    ! ... when -z0 <= z <= z0, potential is
    !
    !           ----  Re[ rho(gz) * exp( i*gz*zright) ]
    !  (-4pi) * >    -----------------------------------
    !           ----                gz^2
    !           gz>0
    !
    !           ----  Re[ rho(gz) * exp( i*gz*zrleft) ]
    !  (-4pi) * >    -----------------------------------
    !           ----                gz^2
    !           gz>0
    !
    !                        ----  Im[ rho(gz) * exp( i*gz*zright) ]
    !  (+4pi) * (z-zright) * >    -----------------------------------
    !                        ----                gz
    !                        gz>0
    !
    !                        ----  Im[ rho(gz) * exp( i*gz*zleft) ]
    !  (+4pi) * (z-zleft)  * >    -----------------------------------
    !                        ----                gz
    !                        gz>0
    !
    !  (-4pi) * ((z-zright)^2 + (z-zleft)^2) * rho(0) / 4
    !
!$omp parallel do default(shared) private(iz, z)
    DO iz = 1, rismt%lfft%nrz
      z = zstart + dz * DBLE(iz - 1)
      rismt%vpot(iz + jgxy) = rismt%vpot(iz + jgxy) + CMPLX( &
      & + fac1 * ( - realr - reall )                &
      & + fac2 * ( + (z - zright) * imager          &
      &            + (z - zleft ) * imagel )        &
      & + fac3 * 0.25_DP * rho0 * (                 &
      &            - (z - zright) * (z - zright)    &
      &            - (z - zleft ) * (z - zleft ) )  &
      & , 0.0_DP, kind=DP)
    END DO
!$omp end parallel do
    !
    ! ... modify reference of solvation potential
    vsolv = 0.0_DP
    vsolu = 0.0_DP
    !
    IF (iref == IREFERENCE_AVERAGE) THEN
      vsolv = 0.0_DP
      vsolu = 0.0_DP
      !
    ELSE IF (iref == IREFERENCE_RIGHT) THEN
      vsolv = fac1 * ( + realr - reall) &
          & + fac2 * ( + zright * imager - zleft * imagel ) &
          & + fac3 * 0.25_DP * rho0 * ( + zright * zright - zleft * zleft )
      vsolu = AIMAG(rismt%vright(igxy))
      !
    ELSE IF (iref == IREFERENCE_LEFT) THEN
      vsolv = fac1 * ( - realr + reall) &
          & + fac2 * ( - zright * imager + zleft * imagel ) &
          & + fac3 * 0.25_DP * rho0 * ( - zright * zright + zleft * zleft )
      vsolu = AIMAG(rismt%vleft(igxy))
    END IF
    !
    vref = vsolv + vsolu
    !
!$omp parallel do default(shared) private(iz)
    DO iz = 1, rismt%lfft%nrz
      rismt%vpot(iz + jgxy) = rismt%vpot(iz + jgxy) - CMPLX(vref, 0.0_DP, kind=DP)
    END DO
!$omp end parallel do
    !
  END IF
  !
  CALL mp_sum(vref, rismt%mp_site%intra_sitg_comm)
  !
  ! ... deallocate memory
  IF (rismt%lfft%ngz_x * rismt%lfft%ngxy > 0) THEN
    DEALLOCATE(rhogt)
    DEALLOCATE(vpott)
  END IF
  IF (rismt%lfft%ngz_x > 0) THEN
    DEALLOCATE(rhogz)
    DEALLOCATE(expigzr)
    DEALLOCATE(expigzl)
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE solvation_esm_potential
!
!---------------------------------------------------------------------------
SUBROUTINE solvation_esm_force(rismt, alpha, force, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... calculate solvation force of ESM(BC1), from Laue-RISM.
  ! ... local potential is derived from Gaussian functions:
  ! ...
  ! ...                      1               [   |r - R|^2 ]
  ! ...   rho(r) = -------------------- * exp[- -----------]
  ! ...             pi^(3/2) * alpha^3       [    alpha^2  ]  .
  !
  ! ... Variables:
  ! ...   alpha: gaussian width (in alat units)
  ! ...   force: solvation force from local potential of ESM(BC1)
  !
  USE cell_base,     ONLY : alat
  USE constants,     ONLY : pi, tpi, e2
  USE control_flags, ONLY : gamma_only
  USE err_rism,      ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE gvect,         ONLY : eigts1, eigts2
  USE ions_base,     ONLY : nat, tau, ityp, zv
  USE kinds,         ONLY : DP
  USE mp,            ONLY : mp_sum
  USE rism,          ONLY : rism_type, ITYPE_LAUERISM
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  REAL(DP),        INTENT(IN)  :: alpha
  REAL(DP),        INTENT(OUT) :: force(1:3, 1:*)
  INTEGER,         INTENT(OUT) :: ierr
  !
  INTEGER                  :: ia
  INTEGER                  :: it
  INTEGER                  :: ipol
  INTEGER                  :: iz
  INTEGER                  :: igxy
  INTEGER                  :: jgxy
  INTEGER                  :: mx, my
  REAL(DP)                 :: z
  REAL(DP)                 :: za
  REAL(DP)                 :: zstart
  REAL(DP)                 :: dz
  REAL(DP)                 :: gx, gy
  REAL(DP)                 :: gxy
  REAL(DP)                 :: qa
  REAL(DP)                 :: mult
  REAL(DP)                 :: rterm1
  REAL(DP)                 :: rterm2
  REAL(DP)                 :: rcoeff
  REAL(DP)                 :: rhogr
  REAL(DP)                 :: rhogi
  REAL(DP)                 :: dvlocr
  REAL(DP)                 :: dvloci
  REAL(DP)                 :: forctmp(3)
  REAL(DP),    ALLOCATABLE :: forcesm(:,:)
  COMPLEX(DP)              :: ccoeff
  COMPLEX(DP)              :: strf_xy
  COMPLEX(DP), ALLOCATABLE :: dvloc(:,:)
  COMPLEX(DP), ALLOCATABLE :: rhogz(:)
  !
  COMPLEX(DP), PARAMETER   :: C_ZERO = CMPLX(0.0_DP, 0.0_DP, kind=DP)
  !
  REAL(DP), EXTERNAL :: qe_erf
  REAL(DP), EXTERNAL :: qe_erfc
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nrzl < rismt%lfft%nrz) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%ngxy < rismt%lfft%ngxy) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... allocate memory
  IF (nat > 0) THEN
    ALLOCATE(forcesm(3, nat))
  END IF
  IF (rismt%lfft%nrz > 0) THEN
    ALLOCATE(dvloc(3, rismt%lfft%nrz))
    ALLOCATE(rhogz(   rismt%lfft%nrz))
  END IF
  !
  ! ... set variables
  zstart = rismt%lfft%zleft + rismt%lfft%zoffset
  dz     = rismt%lfft%zstep
  IF (gamma_only) THEN
    mult = 2.0_DP
  ELSE
    mult = 1.0_DP
  END IF
  !
  ! ... initialize force
  forcesm = 0.0_DP
  !
  ! ...
  ! ... local potential, when Gxy /= 0
  ! ...
  DO igxy = rismt%lfft%gxystart, rismt%lfft%ngxy
    !
    jgxy = rismt%nrzl * (igxy - 1)
    gx   = rismt%lfft%gxy(1, igxy)
    gy   = rismt%lfft%gxy(2, igxy)
    gxy  = rismt%lfft%gnxy(igxy)
    mx   = rismt%lfft%millxy(1, igxy)
    my   = rismt%lfft%millxy(2, igxy)
    !
    IF (rismt%lfft%nrz > 0) THEN
      rhogz(:) = rismt%rhog((1 + jgxy):(rismt%lfft%nrz + jgxy))
    END IF
    !
    DO ia = 1, nat
      !
      it = ityp(ia)
      qa = zv(it)
      za = tau(3, ia)
      !
      ! ... ccoeff means
      !
      !      pi*exp(-i(gx*xa+gy*ya))
      !
      strf_xy = eigts1(mx, ia) * eigts2(my, ia)
      rcoeff  = -qa * e2 * pi
      ccoeff  = rcoeff * strf_xy
      !
      dvloc = C_ZERO
      !
!$omp parallel do default(shared) private(iz, z, rterm1, rterm2)
      DO iz = 1, rismt%lfft%nrz
        z = zstart + dz * DBLE(iz - 1)
        !
        ! ... rterm1 means
        !                              gxy*alpha     z-za
        !     exp( gxy*(z-za)) * erfc(----------- + -------)
        !                                  2         alpha
        !
        ! ... rterm2 means
        !                              gxy*alpha     z-za
        !     exp(-gxy*(z-za)) * erfc(----------- - -------)
        !                                  2         alpha
        !
        ! ... NOTE: to avoid overflows,
        ! ...       exp(var1)*erfc(var2) = exp(var1 + log(erfc(var2))) .
        !
        rterm1 = EXP( tpi * gxy * (z - za) + LOG(qe_erfc(0.5_DP * tpi * gxy * alpha + (z - za) / alpha)))
        rterm2 = EXP(-tpi * gxy * (z - za) + LOG(qe_erfc(0.5_DP * tpi * gxy * alpha - (z - za) / alpha)))
        !
        ! ... derive by X
        !
        !      -i*gx                              [                          gxy*alpha     z-za
        !     ------- * pi*exp(-i(gx*xa+gy*ya)) * [ exp( gxy*(z-za)) * erfc(----------- + -------)
        !       gxy                               [                              2         alpha
        !                                                                    gxy*alpha     z-za    ]
        !                                         + exp(-gxy*(z-za)) * erfc(----------- - -------) ]
        !                                                                        2         alpha   ]
        !
        ! ... derive by Y
        !
        !      -i*gx                              [                          gxy*alpha     z-za
        !     ------- * pi*exp(-i(gx*xa+gy*ya)) * [ exp( gxy*(z-za)) * erfc(----------- + -------)
        !       gxy                               [                              2         alpha
        !                                                                    gxy*alpha     z-za    ]
        !                                         + exp(-gxy*(z-za)) * erfc(----------- - -------) ]
        !                                                                        2         alpha   ]
        ! ... derive by Z
        !
        !                                 [                          gxy*alpha     z-za
        !     - pi*exp(-i(gx*xa+gy*ya)) * [ exp( gxy*(z-za)) * erfc(----------- + -------)
        !                                 [                              2         alpha
        !                                                            gxy*alpha     z-za    ]
        !                                 - exp(-gxy*(z-za)) * erfc(----------- - -------) ]
        !                                                                2         alpha   ]
        !
        dvloc(1, iz) = CMPLX(0.0_DP, -gx / gxy, kind=DP) * ccoeff * (rterm1 + rterm2)
        dvloc(2, iz) = CMPLX(0.0_DP, -gy / gxy, kind=DP) * ccoeff * (rterm1 + rterm2)
        dvloc(3, iz) = -ccoeff * (rterm1 - rterm2)
      END DO
!$omp end parallel do
      !
      forctmp = 0.0_DP
!$omp parallel do default(shared) private(iz, ipol, rhogr, rhogi, dvlocr, dvloci) reduction(+:forctmp)
      DO iz = 1, rismt%lfft%izleft_gedge
        rhogr = -DBLE( rhogz(iz))
        rhogi = -AIMAG(rhogz(iz))
        DO ipol = 1, 3
          dvlocr = DBLE( dvloc(ipol, iz))
          dvloci = AIMAG(dvloc(ipol, iz))
          forctmp(ipol) = forctmp(ipol) - mult * (dvlocr * rhogr + dvloci * rhogi)
        END DO
      END DO
!$omp end parallel do
      forcesm(:, ia) = forcesm(:, ia) + forctmp(:)
      !
      forctmp = 0.0_DP
!$omp parallel do default(shared) private(iz, ipol, rhogr, rhogi, dvlocr, dvloci) reduction(+:forctmp)
      DO iz = rismt%lfft%izright_gedge, rismt%lfft%nrz
        rhogr = -DBLE( rhogz(iz))
        rhogi = -AIMAG(rhogz(iz))
        DO ipol = 1, 3
          dvlocr = DBLE( dvloc(ipol, iz))
          dvloci = AIMAG(dvloc(ipol, iz))
          forctmp(ipol) = forctmp(ipol) - mult * (dvlocr * rhogr + dvloci * rhogi)
        END DO
      END DO
!$omp end parallel do
      forcesm(:, ia) = forcesm(:, ia) + forctmp(:)
      !
    END DO
    !
  END DO
  !
  ! ...
  ! ... local potential, when Gxy = 0
  ! ...
  IF (rismt%lfft%gxystart > 1) THEN
    !
    igxy = 1
    jgxy = rismt%nrzl * (igxy - 1)
    gx   = rismt%lfft%gxy(1, igxy)
    gy   = rismt%lfft%gxy(2, igxy)
    gxy  = rismt%lfft%gnxy(igxy)
    !
    IF (rismt%lfft%nrz > 0) THEN
      rhogz(:) = rismt%rhog((1 + jgxy):(rismt%lfft%nrz + jgxy))
    END IF
    !
    DO ia = 1, nat
      !
      it = ityp(ia)
      qa = zv(it)
      za = tau(3, ia)
      !
      dvloc = C_ZERO
      !
!$omp parallel do default(shared) private(iz, z)
      DO iz = 1, rismt%lfft%nrz
        z = zstart + dz * DBLE(iz - 1)
        !
        ! ... derive by Z
        !
        !               z-za
        !    2pi * erf(-------)
        !               alpha
        !
        dvloc(1, iz) = C_ZERO
        dvloc(2, iz) = C_ZERO
        dvloc(3, iz) = CMPLX((-qa * e2 * tpi) * qe_erf((z - za) / alpha), 0.0_DP, kind=DP)
      END DO
!$omp end parallel do
      !
      forctmp = 0.0_DP
!$omp parallel do default(shared) private(iz, ipol, rhogr, dvlocr) reduction(+:forctmp)
      DO iz = 1, rismt%lfft%izleft_gedge
        rhogr = -DBLE(rhogz(iz))
        DO ipol = 1, 3
          dvlocr = DBLE(dvloc(ipol, iz))
          forctmp(ipol) = forctmp(ipol) - dvlocr * rhogr
        END DO
      END DO
!$omp end parallel do
      forcesm(:, ia) = forcesm(:, ia) + forctmp(:)
      !
      forctmp = 0.0_DP
!$omp parallel do default(shared) private(iz, ipol, rhogr, dvlocr) reduction(+:forctmp)
      DO iz = rismt%lfft%izright_gedge, rismt%lfft%nrz
        rhogr = -DBLE(rhogz(iz))
        DO ipol = 1, 3
          dvlocr = DBLE(dvloc(ipol, iz))
          forctmp(ipol) = forctmp(ipol) - dvlocr * rhogr
        END DO
      END DO
!$omp end parallel do
      forcesm(:, ia) = forcesm(:, ia) + forctmp(:)
      !
    END DO
    !
  END IF
  !
  IF (nat > 0) THEN
    CALL mp_sum(forcesm, rismt%mp_site%intra_sitg_comm)
    force(1:3, 1:nat) = forcesm(1:3, 1:nat) * dz * alat
  END IF
  !
  ! ... deallocate memory
  IF (nat > 0) THEN
    DEALLOCATE(forcesm)
  END IF
  IF (rismt%lfft%nrz > 0) THEN
    DEALLOCATE(dvloc)
    DEALLOCATE(rhogz)
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE solvation_esm_force
!
!---------------------------------------------------------------------------
SUBROUTINE solvation_esm_stress_har( rismt, sigma, rhog_ele, ierr )
  !---------------------------------------------------------------------------
  !
  USE err_rism,      ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  use control_flags, only : gamma_only
  use kinds,         only : DP
  USE rism,          ONLY : rism_type, ITYPE_3DRISM, ITYPE_LAUERISM
  use cell_base,     only : omega, alat, tpiba, at, bg
  use constants,     only : e2, tpi, fpi
  USE lauefft,       ONLY : lauefft_type, inv_lauefft_1z, fw_lauefft_1z_exp, inv_lauefft_1z_exp
  use mp,            only : mp_sum
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  REAL(DP),        INTENT(OUT) :: sigma(3, 3)
  COMPLEX(DP),     INTENT(IN)  :: rhog_ele(rismt%cfft%dfftt%nnr) ! (gz,gp) in unit cell
  INTEGER,         INTENT(OUT) :: ierr
  !
  real(DP), parameter :: delta(2,2) = reshape( (/ 1._DP, 0._DP, 0._DP, 1._DP /), (/2,2/) )
  complex(DP), parameter :: ci = dcmplx(0._DP, 1._DP)
  !
  real(DP) :: z0, zright, zleft, zstart, dz, S
  integer  :: igp, igpz, igz, iz, jz, la, mu
  real(DP) :: gp, gz, g(2), z

  complex(DP) :: rg3
  complex(DP) :: sum1pp, sum1pm, sum1mp, sum1mm
  complex(DP) :: sum2pp, sum2pm, sum2mp, sum2mm
  complex(DP) :: sum1cr, sum1ci, sum2c, sum1sr, sum2si
  real(DP)    :: dgp_deps(2,2)  !! dgp/deps
  real(DP)    :: dgp2_deps(2,2)  !! dgp^2/deps
  real(DP)    :: dinvgp_deps(2,2)  !! dgp^-1/deps

  complex(DP), allocatable :: rhol_ele(:,:) ! (iz,igp) in unit cell
  complex(DP), allocatable :: rholl_ele(:,:) ! (iz,igp) in expanded cell
  complex(DP), allocatable :: rhogl_ele(:,:)  ! (igz,igp) in expanded cell

  complex(DP), allocatable :: Vr(:), dVr_deps(:,:,:)
  complex(DP), allocatable :: workr(:), workg(:)
  TYPE(lauefft_type) :: lauefft0
#if defined(__OPENMP)
  REAL(DP) :: sgomp(2, 2)
#endif

  ! check data type
  IF (rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF

  ! one line FFT
  lauefft0 = rismt%lfft
  lauefft0%ngxy = 1

  ! cell settings
  z0     = 0.5_DP * at(3, 3) * alat
  zright = rismt%lfft%zright * alat
  zleft  = rismt%lfft%zleft * alat
  zstart = rismt%lfft%zleft * alat + rismt%lfft%zoffset * alat
  dz     = rismt%lfft%zstep * alat
  S      = omega/(2*z0)

  ! convert electron density from (gz,gp) to (z,gp) in unit cell
  allocate( rhol_ele(rismt%cfft%dfftt%nr3, rismt%lfft%ngxy) )
  rhol_ele = 0.0_DP
  call inv_lauefft_1z(rismt%lfft, rhog_ele, rhol_ele, rismt%cfft%dfftt%nr3, 1 )

  ! expand electron density from unit cell to expanded cell
  allocate( rholl_ele(rismt%nrzl,           rismt%lfft%ngxy) )
  rholl_ele = 0.0_DP
  do igp = 1, rismt%lfft%ngxy
!$omp parallel do default(shared) private(iz, jz)
    do iz=1, rismt%cfft%dfftt%nr3
      jz = rismt%lfft%izcell_start + iz - 1
      rholl_ele(jz,igp) = rhol_ele(iz,igp)
    end do
  end do

  ! convert electron density from (z,gp) to (gz,gp) in expanded cell
  allocate( rhogl_ele(rismt%lfft%ngz_x, rismt%lfft%ngxy) )
  rhogl_ele = 0.0_DP
  CALL fw_lauefft_1z_exp(rismt%lfft, rholl_ele, rismt%nrzl, rhogl_ele, rismt%lfft%ngz_x )

  deallocate( rhol_ele, rholl_ele )

  allocate( Vr(rismt%lfft%nrz),   dVr_deps(rismt%lfft%nrz,2,2) )
  allocate( workr(rismt%lfft%nrz), workg(rismt%lfft%ngz_x) )

  ! initialize
  sigma(:,:) = 0._DP

  !****For gp!=0 case ********************
  do igp = rismt%lfft%gxystart, rismt%lfft%ngxy
    igpz = rismt%nrzl*(igp-1)
    gp = rismt%lfft%gnxy(igp) * tpiba
    g(:) = rismt%lfft%gxy(:,igp) * tpiba
    if( gp==0._DP ) cycle ! skip gp=0

    ! derivatives by strain tensor
    do la=1, 2
      do mu=1, 2
        dgp_deps(la,mu)    = -g(la)*g(mu)/gp
        dgp2_deps(la,mu)   = -g(la)*g(mu)*2._DP
        dinvgp_deps(la,mu) = +g(la)*g(mu)/gp**3
      end do
    end do

    ! summations over gz
    sum1pp=( 0._DP, 0._DP )
    sum1pm=( 0._DP, 0._DP )
    sum1mp=( 0._DP, 0._DP )
    sum1mm=( 0._DP, 0._DP )
    sum2pp=( 0._DP, 0._DP )
    sum2pm=( 0._DP, 0._DP )
    sum2mp=( 0._DP, 0._DP )
    sum2mm=( 0._DP, 0._DP )
!$omp parallel do default(shared) private(igz,gz,rg3) &
!$omp reduction(+:sum1pp,sum1pm,sum1mp,sum1mm,sum2pp,sum2pm,sum2mp,sum2mm)
    do igz=1, rismt%lfft%ngz_x
      gz = rismt%lfft%gz_x(igz) * tpiba
      rg3 = rhogl_ele(igz,igp)
      sum1pp = sum1pp + rg3*qe_exp(+ci*gz*z0)/(gp+ci*gz)
      sum1pm = sum1pm + rg3*qe_exp(+ci*gz*z0)/(gp-ci*gz)
      sum1mp = sum1mp + rg3*qe_exp(-ci*gz*z0)/(gp+ci*gz)
      sum1mm = sum1mm + rg3*qe_exp(-ci*gz*z0)/(gp-ci*gz)

      sum2pp = sum2pp + rg3*qe_exp(+ci*gz*z0)/(gp+ci*gz)**2
      sum2pm = sum2pm + rg3*qe_exp(+ci*gz*z0)/(gp-ci*gz)**2
      sum2mp = sum2mp + rg3*qe_exp(-ci*gz*z0)/(gp+ci*gz)**2
      sum2mm = sum2mm + rg3*qe_exp(-ci*gz*z0)/(gp-ci*gz)**2
    end do ! igz

    ! calculate dV(z)/deps in center cell -z0<=z<=+z0
!$omp parallel do default(shared) private(iz,z)
    do iz = rismt%lfft%izcell_start, rismt%lfft%izcell_end
      z = zstart + dz*(iz-1)

      dVr_deps(iz,:,:) = &
        +( dgp_deps(:,:)*tpi/gp**2*(1._DP-gp*(z-z0)) + delta(:,:)*tpi/gp ) &
        * exp(+gp*(z-z0)) * sum1pm &
        +( dgp_deps(:,:)*tpi/gp ) &
        * exp(+gp*(z-z0)) * sum2pm &
        +( dgp_deps(:,:)*tpi/gp**2*(1._DP+gp*(z+z0)) + delta(:,:)*tpi/gp ) &
        * exp(-gp*(z+z0)) * sum1mp &
        +( dgp_deps(:,:)*tpi/gp ) &
        * exp(-gp*(z+z0)) * sum2mp
    end do

    ! calculate dV(z)/deps in left expanded cell z<-z0
!$omp parallel do default(shared) private(iz,z)
    do iz = 1, (rismt%lfft%izcell_start - 1)
      z = zstart + dz*(iz-1)

      dVr_deps(iz,:,:) = &
        -( dgp_deps(:,:)*tpi/gp**2*(1._DP-gp*(z+z0)) + delta(:,:)*tpi/gp ) &
        * exp(+gp*(z+z0)) * sum1mm &
        -( dgp_deps(:,:)*tpi/gp ) &
        * exp(+gp*(z+z0)) * sum2mm &
        +( dgp_deps(:,:)*tpi/gp**2*(1._DP-gp*(z-z0)) + delta(:,:)*tpi/gp ) &
        * exp(+gp*(z-z0)) * sum1pm &
        +( dgp_deps(:,:)*tpi/gp ) &
        * exp(+gp*(z-z0)) * sum2pm
    end do

    ! calculate dV(z)/deps in right expanded cell +z0<z
!$omp parallel do default(shared) private(iz,z)
    do iz = (rismt%lfft%izcell_end + 1), rismt%lfft%nrz
      z = zstart + dz*(iz-1)

      dVr_deps(iz,:,:) = &
        -( dgp_deps(:,:)*tpi/gp**2*(1._DP+gp*(z-z0)) + delta(:,:)*tpi/gp ) &
        * exp(-gp*(z-z0)) * sum1pp &
        -( dgp_deps(:,:)*tpi/gp ) &
        * exp(-gp*(z-z0)) * sum2pp &
        +( dgp_deps(:,:)*tpi/gp**2*(1._DP+gp*(z+z0)) + delta(:,:)*tpi/gp ) &
        * exp(-gp*(z+z0)) * sum1mp &
        +( dgp_deps(:,:)*tpi/gp ) &
        * exp(-gp*(z+z0)) * sum2mp
    end do

    do la=1, 2
      do mu=1, 2
        ! calculate bare coulomn term in work(gz)
!$omp parallel do default(shared) private(igz,gz,rg3)
        do igz=1, rismt%lfft%ngz_x
          gz = rismt%lfft%gz_x(igz) * tpiba
          rg3 = rhogl_ele(igz,igp)

          workg(igz) = &
            - delta(la,mu)* fpi*rg3/(gp**2+gz**2) &
            - dgp2_deps(la,mu)* fpi*rg3/(gp**2+gz**2)**2
        end do ! igz

        ! convert work(gz) to work(z)
        call inv_lauefft_1z_exp(lauefft0, workg, rismt%lfft%ngz_x, &
          workr, rismt%nrzl) ! one line FFT

        ! filter  work(z)
!$omp parallel do default(shared) private(iz)
        do iz=1, rismt%lfft%nrz
          if( iz<rismt%lfft%izcell_start .or. rismt%lfft%izcell_end<iz ) then
            workr(iz) = 0.0_DP ! remove bare coulomb in left or right expanded cell.
          end if
        end do

!$omp parallel do default(shared) private(iz)
        do iz = 1, rismt%lfft%nrz
           dVr_deps(iz,la,mu) = dVr_deps(iz,la,mu) + workr(iz)
        end do

      end do ! la
    end do ! mu

    ! modifications
    if( gamma_only ) then
      dVr_deps(:,:,:) = dVr_deps(:,:,:)*2._DP
    end if

    ! calculate stress tensor
!$omp parallel default(shared) private(iz,sgomp)
#if defined(__OPENMP)
    sgomp = 0.0_DP
#endif
!$omp do
    do iz = 1, rismt%lfft%nrz
#if defined(__OPENMP)
      sgomp(1:2,1:2) = sgomp(1:2,1:2) &
         + real( CONJG(rismt%rhog(iz+igpz)) * dVr_deps(iz,1:2,1:2) )
#else
      sigma(1:2,1:2) = sigma(1:2,1:2) &
         + real( CONJG(rismt%rhog(iz+igpz)) * dVr_deps(iz,1:2,1:2) )
#endif
    end do
!$omp end do

#if defined(__OPENMP)
!$omp critical
    sigma(1:2,1:2) = sigma(1:2,1:2) + sgomp(1:2,1:2)
!$omp end critical
#endif
!$omp end parallel

  end do ! igp

  !****For gp=0 case ********************
  if( rismt%lfft%gxystart > 1 ) then
    igp = 1
    igpz = rismt%nrzl*(igp-1)
    gp  = rismt%lfft%gnxy(igp) * tpiba

    ! summations over gz
    sum1cr = ( 0._DP, 0._DP )
    sum1ci = ( 0._DP, 0._DP )
    sum2c  = ( 0._DP, 0._DP )
    sum1sr = ( 0._DP, 0._DP )
    sum2si = ( 0._DP, 0._DP )
!$omp parallel do default(shared) private(igz,gz,rg3) &
!$omp reduction(+:sum1cr,sum1ci,sum2c,sum1sr,sum2si)
    do igz=1, rismt%lfft%ngz_x
      gz = rismt%lfft%gz_x(igz) * tpiba

      if( gz==0._DP ) cycle ! skip gz=0

      rg3 = rhogl_ele(igz,igp)
      sum1ci = sum1ci + rg3*(-ci)*cos(gz*z0)/gz
      sum1cr = sum1cr + rg3*cos(gz*z0)/gz
      sum2c  = sum2c  + rg3*cos(gz*z0)/gz**2
      sum1sr = sum1sr + rg3*sin(gz*z0)/gz
      sum2si = sum2si + rg3*(-ci)*sin(gz*z0)/gz**2
    end do ! igz

    ! calculate V(z) in center cell -z0<=z<=+z0
!$omp parallel do default(shared) private(iz,z,rg3)
    do iz = rismt%lfft%izcell_start, rismt%lfft%izcell_end
      z = zstart + dz*(iz-1)

      rg3 = rhogl_ele(rismt%lfft%gzzero_x,igp)
      Vr(iz) = &
        - tpi*z**2*rg3 - tpi*z0**2*rg3 &
        - fpi*sum2c + fpi*z*sum1ci &
        - fpi*z0 * sum1sr
    end do ! iz

    ! calculate V(z) in left expanded cell z<-z0
!$omp parallel do default(shared) private(iz,z,rg3)
    do iz = 1, (rismt%lfft%izcell_start - 1)
      z = zstart + dz*(iz-1)

      rg3 = rhogl_ele(rismt%lfft%gzzero_x,igp)
      Vr(iz) = &
        + fpi*z0*z*rg3 &
        - fpi*z0*sum1cr &
        - fpi*sum2si + fpi*z*sum1sr
    end do ! iz

    ! calculate V(z) in right expanded cell +z0<z
!$omp parallel do default(shared) private(iz,z,rg3)
    do iz = (rismt%lfft%izcell_end + 1), rismt%lfft%nrz
      z = zstart + dz*(iz-1)

      rg3 = rhogl_ele(rismt%lfft%gzzero_x,igp)
      Vr(iz) = &
        - fpi*z0*z*rg3 &
        + fpi*z0*sum1cr &
        + fpi*sum2si - fpi*z*sum1sr
    end do ! iz

    ! calculate bare coulomn term in work(gz)
    workg(:) = 0.0_DP
!$omp parallel do default(shared) private(igz,gz,rg3)
    do igz=1, rismt%lfft%ngz_x
      gz = rismt%lfft%gz_x(igz) * tpiba
      if( gz==0._DP ) cycle ! skip gz=0

      rg3 = rhogl_ele(igz,igp)
      workg(igz) = fpi*rg3/gz**2
    end do ! igz

    ! convert work(gz) to work(z)
    call inv_lauefft_1z_exp(lauefft0, workg, rismt%lfft%ngz_x, &
      workr, rismt%nrzl) ! one line FFT

    ! filter  work(z)
!$omp parallel do default(shared) private(iz)
    do iz=1, rismt%lfft%nrz
      if( iz<rismt%lfft%izcell_start .or. rismt%lfft%izcell_end<iz ) then
        workr(iz) = 0.0_DP ! remove bare coulomb in left or right expanded cell.
      end if
    end do

    ! add work(z) to V(z)
!$omp parallel do default(shared) private(iz)
    do iz = 1, rismt%lfft%nrz
      Vr(iz) = Vr(iz) + workr(iz)
    end do

    ! calculate dV/deps(z)
!$omp parallel do default(shared) private(iz)
    do iz = 1, rismt%lfft%nrz
      dVr_deps(iz,:,:) = -delta(:,:) * Vr(iz)
    end do ! igz

    ! calculate stress tensor
!$omp parallel default(shared) private(iz,sgomp)
#if defined(__OPENMP)
    sgomp = 0.0_DP
#endif
!$omp do
    do iz = 1, rismt%lfft%nrz
#if defined(__OPENMP)
      sgomp(1:2,1:2) = sgomp(1:2,1:2) &
         + real( CONJG(rismt%rhog(iz+igpz)) * dVr_deps(iz,1:2,1:2) )
#else
      sigma(1:2,1:2) = sigma(1:2,1:2) &
         + real( CONJG(rismt%rhog(iz+igpz)) * dVr_deps(iz,1:2,1:2) )
#endif
    end do
!$omp end do

#if defined(__OPENMP)
!$omp critical
    sigma(1:2,1:2) = sigma(1:2,1:2) + sgomp(1:2,1:2)
!$omp end critical
#endif
!$omp end parallel
  end if ! rismt%lfft%gxystart > 1

  ! e2 means hartree -> Ry.
  sigma(:,:) = sigma(:,:) * (-e2)
  sigma(:,:) = sigma(:,:) * (1.0d0/rismt%lfft%nrz)

  call mp_sum( sigma, rismt%mp_site%inter_sitg_comm)
  call mp_sum( sigma, rismt%mp_site%intra_sitg_comm)

  deallocate( workr, workg )
  deallocate( Vr, dVr_deps )
  deallocate( rhogl_ele )
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
  return

CONTAINS
  complex(DP) function qe_exp( x )
    complex(DP), intent(in) :: x
    real(DP) :: r, i, c, s

    r = dreal(x)
    i = dimag(x)
    c = cos(i)
    s = sin(i)

    qe_exp = exp(r)*cmplx(c,s,kind=DP)

  end function qe_exp
END SUBROUTINE solvation_esm_stress_har
!
!---------------------------------------------------------------------------
SUBROUTINE solvation_esm_stress_locshort( rismt, sigma, vloc, dvloc, ierr )
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
#if defined(__OPENMP)
  REAL(DP) :: sgomp(2, 2)
#endif

  ! check data type
  IF (rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF

  ! ... set solvent density
  rhog_liq => rismt%rhog_pbc
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
#if defined(__OPENMP)
    sgomp = 0.0_DP
#endif
!$omp do
    do ig = 1, rismt%cfft%ngmt
      do la=1, 2
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
    sigma(1:2,1:2) = sigma(1:2,1:2) + sgomp(1:2,1:2)
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
  do la=1, 2
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
  return
END SUBROUTINE solvation_esm_stress_locshort

!
!---------------------------------------------------------------------------
SUBROUTINE solvation_esm_stress_loclong( rismt, sigma, ierr )
  !---------------------------------------------------------------------------
  !
  USE err_rism,      ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  use control_flags, only : gamma_only
  use kinds,         only : DP
  USE rism,          ONLY : rism_type, ITYPE_3DRISM, ITYPE_LAUERISM
  use constants,     only : sqrtpm1, pi, tpi, fpi, e2
  use cell_base,     only : omega, alat, tpiba, at
  use ions_base,     only : zv, nat, tau, ityp
  USE lauefft,       ONLY : lauefft_type, fw_lauefft_1z_exp
  use mp,            only : mp_sum
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  REAL(DP),        INTENT(OUT) :: sigma(3, 3)
  INTEGER,         INTENT(OUT) :: ierr

  complex(DP), parameter :: ci = dcmplx(0._DP, 1._DP)
  real(DP), parameter :: delta(2,2) = reshape( (/ 1._DP, 0._DP, 0._DP, 1._DP /), (/2,2/) )

  real(DP) :: z0, zright, zleft, zstart, dz, S, alpha, salp
  integer  :: igp, igpz, igz, iz, la, mu, ia
  real(DP) :: gp, gz, g(2), z
  real(DP) :: Qa, ra(2), za

  complex(DP) :: rg3
  complex(DP) :: expimgpr, experfcm, experfcp, dexperfcm_dgp, dexperfcp_dgp

  real(DP) :: dgp_deps(2,2)  !! dgp/deps
  real(DP) :: dinvgp_deps(2,2)  !! dgp^-1/deps

  complex(DP), allocatable :: Vr(:), dVr_deps(:,:,:)
#if defined(__OPENMP)
  REAL(DP) :: sgomp(2, 2)
#endif

  real(DP), external   :: qe_erf, qe_erfc

  ! check data type
  IF (rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF

  ! cell settings
  z0     = 0.5_DP * at(3, 3) * alat
  zright = rismt%lfft%zright * alat
  zleft  = rismt%lfft%zleft * alat
  zstart = rismt%lfft%zleft * alat + rismt%lfft%zoffset * alat
  dz     = rismt%lfft%zstep * alat
  S      = omega/(2*z0)
  alpha  = 1._DP
  salp   = sqrt(alpha)

  allocate( Vr(rismt%lfft%nrz),   dVr_deps(rismt%lfft%nrz,2,2) )

  ! initialize
  sigma(:,:) = 0._DP

  !****For gp!=0 case ********************
  do igp = rismt%lfft%gxystart, rismt%lfft%ngxy
    igpz = rismt%nrzl*(igp-1)
    gp = rismt%lfft%gnxy(igp) * tpiba
    g(:) = rismt%lfft%gxy(:,igp) * tpiba

    if( gp==0._DP ) cycle ! skip gp=0

    ! derivatives by strain tensor
    do la=1, 2
      do mu=1, 2
        dgp_deps(la,mu)    = -g(la)*g(mu)/gp
        dinvgp_deps(la,mu) = +g(la)*g(mu)/gp**3
      end do
    end do

    dVr_deps = ( 0._DP, 0._DP )
    ! calculate dV(z)/deps in center cell -z0<=z<=+z0
!$omp parallel do default(shared) private(iz,z,ia,Qa,ra,za) &
!$omp private(expimgpr,experfcm,experfcp,dexperfcm_dgp,dexperfcp_dgp)
    do iz = rismt%lfft%izcell_start, rismt%lfft%izcell_end
      z = zstart + dz*(iz-1)

      ! summations over all atoms
      do ia=1, nat
        Qa = (-1._DP)*zv(ityp(ia))
        ra(1:2) = tau(1:2,ia)*alat
        za = tau(3,ia)*alat

        expimgpr = qe_exp( - ci*(g(1)*ra(1) + g(2)*ra(2)) )

        experfcm = exp_erfc( -gp*(z-za), gp/2.d0/salp-salp*(z-za) )
        experfcp = exp_erfc( +gp*(z-za), gp/2.d0/salp+salp*(z-za) )

        dexperfcm_dgp = -(z-za)*exp_erfc( -gp*(z-za), gp/2.d0/salp-salp*(z-za) ) &
          - sqrtpm1/salp * exp( -(gp/2.d0/salp)**2 - (salp*(z-za))**2 )
        dexperfcp_dgp = +(z-za)*exp_erfc( +gp*(z-za), gp/2.d0/salp+salp*(z-za) ) &
          - sqrtpm1/salp * exp( -(gp/2.d0/salp)**2 - (salp*(z-za))**2 )

        dVr_deps(iz,:,:) = dVr_deps(iz,:,:) &
          + dinvgp_deps(:,:) * gp * pi/gp * Qa/S * expimgpr * experfcm &
          - delta(:,:)            * pi/gp * Qa/S * expimgpr * experfcm &
          + dgp_deps(:,:)         * pi/gp * Qa/S * expimgpr * dexperfcm_dgp

        dVr_deps(iz,:,:) = dVr_deps(iz,:,:) &
          + dinvgp_deps(:,:) * gp * pi/gp * Qa/S * expimgpr * experfcp &
          - delta(:,:)            * pi/gp * Qa/S * expimgpr * experfcp &
          + dgp_deps(:,:)         * pi/gp * Qa/S * expimgpr * dexperfcp_dgp
      end do ! ia
    end do ! iz

    ! calculate dV(z)/deps in left expanded cell z<-z0
!$omp parallel do default(shared) private(iz,z,ia,Qa,ra,za) &
!$omp private(expimgpr,experfcm,experfcp,dexperfcm_dgp,dexperfcp_dgp)
    do iz = 1, (rismt%lfft%izcell_start - 1)
      z = zstart + dz*(iz-1)

      ! summations over all atoms
      do ia=1, nat
        Qa = (-1._DP)*zv(ityp(ia))
        ra(1:2) = tau(1:2,ia)*alat
        za = tau(3,ia)*alat

        expimgpr = qe_exp( - ci*(g(1)*ra(1) + g(2)*ra(2)) )
        experfcp = 2._DP*exp( +gp*(z-za) )
        dexperfcp_dgp = -(z-za)*experfcp

        dVr_deps(iz,:,:) = dVr_deps(iz,:,:) &
          + dinvgp_deps(:,:) * gp * pi/gp * Qa/S * expimgpr * experfcp &
          - delta(:,:)            * pi/gp * Qa/S * expimgpr * experfcp &
          + dgp_deps(:,:)         * pi/gp * Qa/S * expimgpr * dexperfcp_dgp
      end do ! ia
    end do

    ! calculate dV(z)/deps in right expanded cell +z0<z
!$omp parallel do default(shared) private(iz,z,ia,Qa,ra,za) &
!$omp private(expimgpr,experfcm,experfcp,dexperfcm_dgp,dexperfcp_dgp)
    do iz = (rismt%lfft%izcell_end + 1), rismt%lfft%nrz
      z = zstart + dz*(iz-1)

      ! summations over all atoms
      do ia=1, nat
        Qa = (-1._DP)*zv(ityp(ia))
        ra(1:2) = tau(1:2,ia)*alat
        za = tau(3,ia)*alat

        expimgpr = qe_exp( - ci*(g(1)*ra(1) + g(2)*ra(2)) )
        experfcm = 2._DP*exp( -gp*(z-za) )
        dexperfcm_dgp = +(z-za)*experfcm

        dVr_deps(iz,:,:) = dVr_deps(iz,:,:) &
          + dinvgp_deps(:,:) * gp * pi/gp * Qa/S * expimgpr * experfcm &
          - delta(:,:)            * pi/gp * Qa/S * expimgpr * experfcm &
          + dgp_deps(:,:)         * pi/gp * Qa/S * expimgpr * dexperfcm_dgp
      end do ! ia
    end do

    ! modifications
    if( gamma_only ) then
      dVr_deps(:,:,:) = dVr_deps(:,:,:)*2._DP
    end if

    ! calculate stress tensor
!$omp parallel default(shared) private(iz,sgomp)
#if defined(__OPENMP)
    sgomp = 0.0_DP
#endif
!$omp do
    do iz = 1, rismt%lfft%nrz
#if defined(__OPENMP)
      sgomp(1:2,1:2) = sgomp(1:2,1:2) &
         - real( CONJG(rismt%rhog(iz+igpz)) * dVr_deps(iz,1:2,1:2) )
#else
      sigma(1:2,1:2) = sigma(1:2,1:2) &
         - real( CONJG(rismt%rhog(iz+igpz)) * dVr_deps(iz,1:2,1:2) )
#endif
    end do
!$omp end do

#if defined(__OPENMP)
!$omp critical
    sigma(1:2,1:2) = sigma(1:2,1:2) + sgomp(1:2,1:2)
!$omp end critical
#endif
!$omp end parallel
  end do ! igp

  !****For gp=0 case ********************
  if( rismt%lfft%gxystart > 1 ) then
    igp = 1
    igpz = rismt%nrzl*(igp-1)
    gp  = rismt%lfft%gnxy(igp) * tpiba

    Vr(:) = 0._DP

    ! calculate V(z) in center cell -z0<=z<=+z0
!$omp parallel do default(shared) private(iz,z,ia,Qa,ra,za)
    do iz = rismt%lfft%izcell_start, rismt%lfft%izcell_end
      z = zstart + dz*(iz-1)

      do ia=1, nat
        Qa = (-1._DP)*zv(ityp(ia))
        ra(1:2) = tau(1:2,ia)*alat
        za = tau(3,ia)*alat

        Vr(iz) = Vr(iz) - tpi * Qa/S &
            * ( (z-za)*qe_erf(salp*(z-za)) &
          + exp(-alpha*(z-za)**2)*sqrtpm1/salp )
      end do ! ia
    end do ! iz

    ! calculate V(z) in left expanded cell z<-z0
!$omp parallel do default(shared) private(iz,z,ia,Qa,ra,za)
    do iz = 1, (rismt%lfft%izcell_start - 1)
      z = zstart + dz*(iz-1)

      do ia=1, nat
        Qa = (-1._DP)*zv(ityp(ia))
        ra(1:2) = tau(1:2,ia)*alat
        za = tau(3,ia)*alat

        Vr(iz) = Vr(iz) + tpi * Qa/S * (z-za)
      end do ! ia
    end do ! iz

    ! calculate V(z) in right expanded cell +z0<z
!$omp parallel do default(shared) private(iz,z,ia,Qa,ra,za)
    do iz = (rismt%lfft%izcell_end + 1), rismt%lfft%nrz
      z = zstart + dz*(iz-1)

      do ia=1, nat
        Qa = (-1._DP)*zv(ityp(ia))
        ra(1:2) = tau(1:2,ia)*alat
        za = tau(3,ia)*alat

        Vr(iz) = Vr(iz) - tpi * Qa/S * (z-za)
      end do ! ia
    end do ! iz

    ! calculate dV/deps(z)
!$omp parallel do default(shared) private(igz)
    do iz = 1, rismt%lfft%nrz
      dVr_deps(iz,:,:) = -delta(:,:) * Vr(iz)
    end do

    ! calculate stress tensor
!$omp parallel default(shared) private(iz,sgomp)
#if defined(__OPENMP)
    sgomp = 0.0_DP
#endif
!$omp do
    do iz = 1, rismt%lfft%nrz
#if defined(__OPENMP)
      sgomp(1:2,1:2) = sgomp(1:2,1:2) &
       - real( CONJG(rismt%rhog(iz+igpz)) * dVr_deps(iz,1:2,1:2) )
#else
      sigma(1:2,1:2) = sigma(1:2,1:2) &
       - real( CONJG(rismt%rhog(iz+igpz)) * dVr_deps(iz,1:2,1:2) )
#endif
    end do
!$omp end do

#if defined(__OPENMP)
!$omp critical
    sigma(1:2,1:2) = sigma(1:2,1:2) + sgomp(1:2,1:2)
!$omp end critical
#endif
!$omp end parallel
  end if ! rismt%lfft%gxystart > 1

  ! e2 means hartree -> Ry.
  sigma(:,:) = sigma(:,:)*(e2)
  sigma(:,:) = sigma(:,:)*(1.0d0/rismt%lfft%nrz)

  call mp_sum( sigma, rismt%mp_site%inter_sitg_comm )
  call mp_sum( sigma, rismt%mp_site%intra_sitg_comm )

  deallocate( Vr, dVr_deps )
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
  return

CONTAINS
  complex(DP) function qe_exp( x )
    complex(DP), intent(in) :: x
    real(DP) :: r, i, c, s

    r = dreal(x)
    i = dimag(x)
    c = cos(i)
    s = sin(i)

    qe_exp = exp(r)*cmplx(c,s,kind=DP)

  end function qe_exp

  ! exp(x) * erfc(y)
  ! This function is to avoid INFINITY * ZERO for large positive x and y.
  real(8) function exp_erfc (x, y)
    implicit none
    real(8), intent(in) :: x, y
    real(8)             :: ym, ym2, nume, deno
    real(8), parameter  :: rtpim = 0.564189583547756279d0 ! 1 / sqrt(PI)
    real(8), parameter  :: r(0:4) = (/ &
      -2.99610707703542174d-3, -4.94730910623250734d-2, &
      -2.26956593539686930d-1, -2.78661308609647788d-1, &
      -2.23192459734184686d-2 /)
    real(8), parameter  :: s(0:4) = (/ &
      1.06209230528467918d-2, 1.91308926107829841d-1, &
      1.05167510706793207d0,  1.98733201817135256d0,  &
      1.00000000000000000d0 /)

    if( x < 709.0d0 .or. y < 4.0d0 ) then
      exp_erfc = exp(x) * qe_erfc(y)
    else
      ym  = 1d0 / y
      ym2 = ym ** 2
      nume =( ( ( r(4) * ym2 + r(3) ) * ym2 + r(2) ) * ym2 + r(1) ) * ym2 + r(0)
      deno =( ( ( s(4) * ym2 + s(3) ) * ym2 + s(2) ) * ym2 + s(1) ) * ym2 + s(0)
      exp_erfc = exp( - y**2 + x) * ym * ( rtpim + ym2 * nume / deno )
    end if

    return
  end function exp_erfc

END SUBROUTINE solvation_esm_stress_loclong
