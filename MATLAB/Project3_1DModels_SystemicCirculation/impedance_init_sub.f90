subroutine impedance_init_driver(tmstps)
  use f90_tools
  use root_imp 
  implicit none

  integer, intent (in) :: tmstps

  call impedance_init(tmstps)

end subroutine impedance_init_driver
