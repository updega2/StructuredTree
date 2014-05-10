subroutine impedance_driver(tmstps,Period,ff1,ff2,ff3,rho,mu,r_root,r_min,y_xt,Lr,Fr2,q,g)
  use f90_tools
  use root_imp 
  implicit none

  integer, intent(in)    :: tmstps
  real(lng), intent(in)  :: Period,ff1,ff2,ff3,rho,mu,r_root,r_min,Lr,Fr2,q,g
  real(lng), intent(out) :: y_xt(tmstps)


  call impedance (tmstps,Period,ff1,ff2,ff3,rho,mu,r_root,r_min,y_xt,Lr,Fr2,q,g)
end subroutine impedance_driver
