!> This module contains UGWP v0 initialization schemes
module ugwp_common_v0
  use machine,  only: kind_phys
  implicit none

  real(kind=kind_phys), parameter :: dw2min=1.0, bnv2min=1.e-6
  real(kind=kind_phys) :: grcp = 1.0E30
  real(kind=kind_phys) :: rgrav = 1.0E30
  real(kind=kind_phys) :: rdi = 1.0E30
  real(kind=kind_phys) :: gor = 1.0E30
  real(kind=kind_phys) :: gr2 = 1.0E30
  real(kind=kind_phys) :: gocp = 1.0E30
  real(kind=kind_phys) :: rcpd = 1.0E30
  real(kind=kind_phys) :: rcpd2 = 1.0E30
  real(kind=kind_phys) :: pi2 = 1.0E30
  real(kind=kind_phys) :: omega1 = 1.0E30
  real(kind=kind_phys) :: omega2 = 1.0E30
  real(kind=kind_phys) :: rad_to_deg = 1.0E30
  real(kind=kind_phys) :: deg_to_rad = 1.0E30
  real(kind=kind_phys) :: velmin = 1.0E30
  real(kind=kind_phys) :: pi = 1.0E30
  real(kind=kind_phys) :: grav = 1.0E30
  real(kind=kind_phys) :: rd = 1.0E30
  real(kind=kind_phys) :: rv = 1.0E30
  real(kind=kind_phys) :: cpd = 1.0E30
  real(kind=kind_phys) :: fv = 1.0E30
  real(kind=kind_phys) :: arad = 1.0E30

contains

  subroutine ugwp_common_v0_init(con_pi, con_g, con_rd, con_rv, &
       con_cp, con_fvirt, con_rerth)
    real(kind=kind_phys), intent(in) :: con_pi, con_g, con_rd, con_rv
    real(kind=kind_phys), intent(in) :: con_cp, con_fvirt, con_rerth

    pi = con_pi
    grav = con_g
    rd = con_rd
    rv = con_rv
    cpd = con_cp
    fv = con_fvirt
    arad = con_rerth

    grcp = grav/cpd
    rgrav = 1.0d0/grav
    rdi  = 1.0d0/rd
    gor  = grav/rd
    gr2   = grav*gor
    gocp = grav/cpd
    rcpd = 1./cpd
    rcpd2 = 0.5*rcpd
    pi2  = pi + pi
    omega1 = pi2/86400.0
    omega2 = omega1+omega1
    rad_to_deg=180.0/pi
    deg_to_rad=pi/180.0
    velmin=sqrt(dw2min)
  end subroutine ugwp_common_v0_init
end module ugwp_common_v0
