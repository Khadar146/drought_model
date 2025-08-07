!> @file lib/cdf.f90
!! Wrapper module to expose selected CDFLIB routines
module cdf
  implicit none
  public

  interface
    function gaminv(p, a, q, status, bound)
      real(8), intent(in) :: p, a
      real(8)             :: q
      integer             :: status
      real(8)             :: bound
    end function gaminv

    function ppnd(p, status)
      real(8), intent(in) :: p
      integer             :: status
      real(8)             :: ppnd
    end function ppnd
  end interface
end module cdf