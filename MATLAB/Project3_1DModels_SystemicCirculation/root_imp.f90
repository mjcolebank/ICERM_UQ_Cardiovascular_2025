!  MODULE ROOT_IMP

!***************************************************************************
!*                                                                         *
!* Module: root_imp.f90                                                    *
!* Version: 1.0                                                            *
!* By: Mette Olufsen, Math-Tech                                            *
!* Last modified: 14. Jan. 1997                                            * 
!*                                                                         *
!* Contains subroutines and functions needed to compute the impedance at   *
!* the root of a structured tree (with constant asymmetry-ratio and a      *
!* geometry where the size of all branches depends of the radius of the    *
!* branch at the root of the tree.                                         *
!*                                                                         *
!***************************************************************************
module root_imp
use f90_tools     ! Contains several tools used for when computing FFT.
implicit none

private
public impedance  ! The subroutine called from the c-program arteries.cxx
                  ! (the driver routine)
public impedance_init, impedance_close


! The number of generations in the structured tree
integer, parameter   :: Maxgen = 100

integer count, alpha_max, beta_max

! The asymmetry ratio of the structured tree
!real(lng), parameter :: asym   = 0.4048

! Exponent in radius relation.
!real(lng), parameter :: expo   = 2.76_lng

! A temporary matrix for storing root impedances in parts of the
! structured tree, those which are repeated because of the constant
! asymmetry ratio.
complex(lng)         :: Computed(0:Maxgen,0:Maxgen)

complex(lng), allocatable   :: Z_om(:)

integer localmax


contains

!***************************************************************************
!*                                                                         *
!* Function: Z0func                                                        *
!* Version: 1.0                                                            *
!* By: Mette Olufsen, Math-Tech                                            *
!* Last modified: 14. Jan. 1997                                            * 
!*                                                                         *
!* Computes the impedance as a function of Omega using a recursive         *
!* formula for a asymmetric tree with constant asymmetry-ratio.            *
!*                                                                         *
!* This subroutine is called from function comp_imp by:                    *
!*                                                                         *
!* Z0func (omega_k,alpha_pow,beta_pow,ff*,rho,mu,r_root,r_min,Lr,Fr2,q,g)  *
!*                                                                         *
!* where:                                                                  *
!*                                                                         *
!* Z0func    Returns the frequency dep. impedance (Z0) from an asym. tree. *
!* omega_k   Frequency for which to compute Z0.                            *
!* alpha_pow The power of alpha in the structured tree.                    *
!* beta_pow  The power of beta in the structured tree.                     *
!* ff*       Constants to model Eh/r = ff1*exp(ff2*r)+ff3, g/cm/s^2.       *
!* ff*       Constants to model Eh/r=-ff1*tanh(ff2*(r-ff3))+ff4, g/cm/s^2. *
!* rho       The density of the blood, g/cm^3.                             *
!* mu        Bloodplasma viscosity, g/cm/s.                                *
!* r_root    The input radius at the root of the structured tree, cm.      *
!* r_min     The minimum bottom radius, gives the order of the tree, cm.   *
!* Lr        Characteristic length scale (radius), cm.                     *
!* Fr2       The squared Froude number.                                    *
!* q         Characteristic flow-rate, cm^3/s.                             *
!* g         The gravitational force, cm/s^2.                              *
!*                                                                         *
!***************************************************************************
recursive function Z0func (omega_k,alpha_pow,beta_pow,ff1,ff2,ff3,rho,mu,r_root,r_min, &
                                    trm_rst,Lr,Fr2,q,g,par1_in,par2_in,lrr) result (Z0)
implicit none

  integer, intent(in)   :: alpha_pow, beta_pow
  real(lng), intent(in) :: omega_k,ff1,ff2,ff3,rho,r_root,r_min,trm_rst,Lr,Fr2,q,g,mu
  real(lng), intent(in) :: par1_in,par2_in,lrr
  real(lng) ::  diameter, mu_D, eta_D,mu1
  
  integer      :: generations
  real(lng)    :: nu, r, r_d, l, A, A_d, D, wom
  real(lng)    :: C_exp, Hct, Hct_term
  real(lng)    :: beta, alpha
  complex(lng) :: i, g_omega, c_omega, kappa, Z0, ZL, Zl_0, Zr_0, t1, t2

  ! Physical constants.
  i     = cmplx(0.0_lng,1.0_lng,lng) ! The complex unit.

    ! ADDED BY MJC
    ! if you want to use the traditional structured tree parameters, one of the options will be greater than 1
    ! otherwise, we will assume you pass alpha and beta in
    if (par2_in>1) then
        alpha  = ((par1_in**(par2_in/2.0)+1.0)**(-1.0/par2_in))**alpha_pow
        beta = (sqrt(par1_in)*(par1_in**(par2_in/2.0)+1.0)**(-1.0/par2_in))**beta_pow
!        beta  = ((par1_in**(par2_in/2.0)+1.0)**(-1.0/par2_in))**alpha_pow
!        alpha = (sqrt(par1_in)*(par1_in**(par2_in/2.0)+1.0)**(-1.0/par2_in))**beta_pow
    else
      alpha  = par1_in**alpha_pow
      beta = par2_in**beta_pow
    endif

    if(alpha_pow>alpha_max) then
        alpha_max = alpha_pow
    endif
    if(beta_pow>beta_max) then
        beta_max = beta_pow
    endif

  generations = alpha_pow + beta_pow  ! The current generation.
  r_d  = alpha*beta*r_root     ! Radius at root, cm.
  A_d  = pi*r_d**2.0           ! Cross-sectional area, cm^2.
  r    = r_d                   ! Radius at root, dimension-less.
  A    = A_d                   ! Cross-sectional area, dimension-less.
      
  l = lrr*r



  diameter = 2.0*r*(10.0**4.0)
  eta_D = 3.2+6*exp(-0.085*diameter)-2.44*exp(-0.06*(diameter**0.645))
  !  Terms that depend on hematocrit
  C_exp = (0.8+exp(-0.075*diameter))*(1.0/(1.0+10**(-11.0) * diameter**12.0) - 1.0) + (1.0/(1.0+10**(-11.0) * diameter**12.0))
  Hct = 0.45
  Hct_term = ((1.0 - Hct)**C_exp - 1.0)/((1.0-0.45)**C_exp - 1.0)
  mu_D  = (1.0 + (eta_D - 1.0)*Hct_term*(diameter/(diameter - 1.1))**2.0)*(diameter/(diameter - 1.1))**2.0
  mu1 = mu!*mu_D/3.2
  
  nu   = mu1/rho                ! Kinematic blood viscosity, m^2/s.
  D    = 1.0/(ff1*exp(ff2*r_d)+ff3)*3.0*A_d/2.0 ! Distensibility.
  wom  = r_d*sqrt(omega_k/nu)  ! Womersleys parameter.




  ! Temporary functions of r.
  if (wom > 3.0) then 
    g_omega = sqrt(D*A/rho)* &
              sqrt(1.0_lng-2.0/i**(0.5)/wom*(1.0_lng+1.0/2.0/wom)) 
    c_omega = sqrt(A/D/rho)* &
              sqrt(1.0_lng-2.0/i**(0.5)/wom*(1.0_lng+1.0/2.0/wom))
  else
    if (wom > 2.0) then
      g_omega = sqrt(D*A/rho)*((3.0_lng-wom)* &
                sqrt(i*wom**2.0/8.0+wom**4.0/48.0) + &
                (wom-2.0_lng)*&
                sqrt(1.0_lng-2.0/i**(0.5)/wom*(1.0_lng+1.0/2.0/wom)))
      c_omega = sqrt(A/D/rho)*((3.0_lng-wom)* &
                sqrt(i*wom**2.0/8.0+wom**4.0/48.0) + &
                (wom-2.0_lng)*&
                sqrt(1.0_lng-2.0/i**(0.5)/wom*(1.0_lng+1.0/2.0/wom)))
    else
      if (wom == 0) then
        g_omega = 0.0
        c_omega = 0.0
      else
        g_omega = sqrt(D*A/rho)*sqrt(i*wom**2/8+wom**4/48)
        c_omega = sqrt(A/D/rho)*sqrt(i*wom**2/8+wom**4/48)
      end if
    end if
  end if
 
  ! Temporary function of omega_k. 
   if (omega_k /= 0) then
    kappa  = omega_k*l/c_omega
  else
    count = count+1
    kappa   = 0.0
  end if

  ! Determine the resistance of the root of the terminal branches.
  ! if (generations >= Maxgen) 
  if (r <= r_min) then
    if (generations >= localmax) then
      localmax = generations
    end if
    !if (generations >= Maxgen) then
    !  write (*,*) 'Generations level exceeded'
    !end if
    if (omega_k == 0) then
    ZL = trm_rst
    else
    ZL = trm_rst
    end if
    !ZL = trm_rst
  else
    ! Get Z0 recursively at reduced frequencies.
    if (abs(Computed(alpha_pow+1, beta_pow)) /= 0.0) then
      Zl_0 = Computed(alpha_pow+1, beta_pow)
    else  
      Zl_0 = Z0func (omega_k,alpha_pow+1,beta_pow,ff1,ff2,ff3,rho,mu,r_root,r_min,trm_rst,Lr,Fr2,q,g,par1_in,par2_in,lrr)
    end if

    if (abs(Computed(alpha_pow, beta_pow+1)) /= 0.0) then
      Zr_0 = Computed(alpha_pow, beta_pow+1)
    else
      Zr_0 = Z0func (omega_k,alpha_pow,beta_pow+1,ff1,ff2,ff3,rho,mu,r_root,r_min,trm_rst,Lr,Fr2,q,g,par1_in,par2_in,lrr)
    end if

    ! Prediction of the resulting impedance from the recursion formula.
    ZL = rho*g*Lr/(q/Zl_0 + q/Zr_0)
  end if
 
  ! Finally get Z0(omega) as theoretically derived.
  if (g_omega /= 0.0) then 
    t1   = i*sin(kappa)/g_omega + cos(kappa)*ZL
    t2   = cos(kappa) + i*g_omega*sin(kappa)*ZL
    Z0   = q/(rho*g*Lr)*(t1/t2)
  else 
    Z0   = q/(rho*g*Lr)*(8.0*mu1*l/(A*r**2.0) + ZL)
  end if

  Computed(alpha_pow,beta_pow) = Z0

end function Z0func



!***************************************************************************
!*                                                                         *
!* Function: comp_imp                                                      *
!* Version: 1.0                                                            *
!* By: Mette Olufsen, Math-Tech                                            *
!* Last modified: 14. Jan. 1997                                            * 
!*                                                                         *
!* Compute the impedance using a recursive formula for an asymmetric       *
!* tree.                                                                   *
!*                                                                         *
!* This subroutine is called from subroutine impedance by:                 *
!*                                                                         *
!* Z_om = comp_imp (N,Omega,ff*,rho,mu,r_root,Lr,Fr2,q,g)                  *
!*                                                                         *
!* comp_imp  Returns frequency dependent impedances (Z_om) at the root of  *
!*           the structured asymmetric tree.                               *
!* N         Number of frequencies (in Omega).                             *  
!* Omega     A vector containing N frequencies.                            *
!* ff*       Constants to model Eh/r=ff1*exp(ff2*r)+ff3, g/cm/s^2.         *
!* ff*       Constants to model Eh/r=-ff1*tanh(ff2*(r-ff3))+ff4, g/cm/s^2. *
!* rho       The density of the blood, g/cm^3.                             *
!* mu        Bloodplasma viscosity, g/cm/s.                                *
!* r_root    The input radius at the root of the structured tree, cm.      *
!* Lr        Characteristic length scale (radius), cm.                     *
!* Fr2       The squared Froude number.                                    *
!* q         Characteristic flow-rate, cm^3/s.                             *
!* g         The gravitational force, cm/s^2.                              *
!*                                                                         *
!* This function routine is a driver routine that uses Z0func to compute   *
!* Z_om at each generation.                                                *
!*                                                                         *
!***************************************************************************
function comp_imp (N,Omega,trm_rst,ff1,ff2,ff3,rho,mu,r_root,r_min,Lr,Fr2,q,g,par1_in,par2_in,lrr) result (Z_om)
implicit none

  integer, intent(in)       :: N
  real(lng), intent(in)     :: Omega(:),trm_rst,ff1,ff2,ff3,rho,r_root,Lr,Fr2,q,g,r_min,mu
  real(lng), intent(in)     :: par1_in,par2_in,lrr
  complex(lng)              :: temp, Z_om(N)
  !real(lng)                 :: mu1
  integer                   :: k

  ! Omega contains tmstps+1 frequencies because when computing the
  ! impedance it is easier to invert it when it is computed for all positive
  ! frequencies and since the interval goes from [-N/2/Tper:N/2/Tper-df],
  ! we include the frequency N/2/Tper in Omega and from this we get
  ! Z(-N/2/Tper) we then end up throwing pushing Z(N/2/Tper) out.
 
  Z_om = cmplx (0.0_lng,0.0_lng,lng) ! Initialize Z_om with zeros

  ! For all the positive frequencies compute the impedance.
  ! Since we know that the system is self-adjoint we don't
  ! have to compute the negative frequencies as well. 
  do k = N/2+1, N+1
    Computed  = 0.0_lng ! For each frequency make sure we don't carry
                        ! any values with us from the computations for the
                        ! previous frequency, so make sure Computed is zero.
    ! Since Z_om only has N places leave result for the frequency k at
    ! Z_om(k-1) we will later make up for this.

    Z_om(k-1) = Z0func (Omega(k),0,0,ff1,ff2,ff3,rho,mu,r_root,r_min,trm_rst,Lr,Fr2,q,g,par1_in,par2_in,lrr)
  end do 
  
  temp = Z_om (N/2)
  
  ! Apply self-adjointness 
  Z_om(1:N/2)   = conjg(flipud(Z_om(N/2+1:N)))

  ! Shift the results for the positive frequencies one to the right
  ! to make up for the above.
  Z_om(N/2+1:N) = eoshift(Z_om(N/2+1:N),-1)    

  ! Insert Z(0,0) as described in the theoretical derivation
  Z_om(N/2+1) = temp
end function comp_imp


!***************************************************************************
!*                                                                         *
!* Subroutine: impedance                                                   *
!* Version: 1.0                                                            *
!* By: Mette Olufsen, Math-Tech                                            *
!* Last modified: 14. Jan. 1997                                            * 
!*                                                                         *
!*                                                                         *
!* Compute the Impedance as a function of space and time at the root of a  *
!* structured tree.                                                        *
!*                                                                         *
!* This subroutine is called from the c++ subroutine arteries.cxx whenever *
!* a new end-tube (of the tree of larger arteries) is constructed. We then *
!* use the resulting impedances to get a time-dependent boundary-condition.*
!* It is called by:                                                        *
!*                                                                         *
!* Impedance (tmstps,Tper,ff*,rho,mu,r_root,r_min,z_xt,Lr,Fr2,q,g)         *
!*                                                                         *
!* where:                                                                  *
!*                                                                         *
!* tmstps    Number of sample points in time (Must be a power of 2)!       *
!* Tper      The period of the simulation, s.                              *
!* ff*       Constants to model Eh/r = ff1*exp(ff2*r)+ff3, g/cm/s^2.       *
!* ff*       Constants to model Eh/r=-ff1*tanh(ff2*(r-ff3))+ff4, g/cm/s^2. *
!* rho       Blood density, g/cm^3.                                        *
!* mu        Bloodplasma viscosity, g/cm/s.                                *
!* r_root    The radius of the root of the structured tree, cm.            *
!* r_min     The terminal radius of the structured tree, cm.               *
!* z_xt      Will return the impedances at the root of the structured tree *
!*           in the time domain.                                           *
!* Lr        Characteristic length scale (radius), cm.                     *
!* Fr2       The squared Froude number.                                    *
!* q         Characteristic flow-rate, cm^3/s.                             *
!* g         The gravitational force, cm/s^2.                              *
!*                                                                         *
!***************************************************************************
subroutine impedance (tmstps,Period,ff1,ff2,ff3,rho,mu,r_root,r_min,y_xt,Lr,Fr2,q,g,par1_in,par2_in,lrr,term_ID,Z_flag)
implicit none

  integer,   intent(in)      :: tmstps
  real(lng), intent(in)      :: Period,ff1,ff2,ff3,rho,Lr,Fr2,q,g,r_root,r_min,mu
  real(lng), intent(in)      :: par1_in,par2_in,lrr
  real(lng)                  :: y_xt(tmstps)
  integer,   intent(in)      :: term_ID,Z_flag

  integer                    :: j
  real(lng)                  :: df, Freq(tmstps+1), Omega(tmstps+1), trm_rst
  complex(lng)               :: Z_hat(tmstps), Y_hat(tmstps)!,Zin

  integer, parameter                  :: nbuf = 3, f1 = 10
  !character (len=30)                  :: fn
  character (len=10)                  :: term_name
!  character (len=40), dimension(nbuf) :: buffer  ! Temporary strings
!  integer k

  ! Physical parameters
  df     = 1/Period                            ! Frequency interval. 
  Freq   = (/ (j*df, j=-tmstps/2, tmstps/2) /) ! Frequency-vector (abscissae). 
  Omega  = 2*pi*Freq                           ! Freq.-vector scaled by a factor 2pi.


  !call counting (0,0,alpha,beta,r_root,r_min,nb_terms)

  localmax = 0
  trm_rst = 0.1!0.05     ! Terminal resistance could be (nb_terms*resist)

   count = 0
   alpha_max = 0
   beta_max = 0
    write(term_name,'(I3)') term_ID
    if (Z_flag==1) then
      ! Compute the impedance at the root of the structured tree.
      Z_om =comp_imp (tmstps,Omega,trm_rst,ff1,ff2,ff3,rho,mu,r_root,r_min,Lr,Fr2,q,g,par1_in,par2_in,lrr)



      ! Transform P back to the time domain.
      ! Divide by tmstps to approximate continuous inv. Fourier transform.
      ! In particular, amplitude must be independent of resolution.
       Z_hat = Z_om
      Y_hat = 1.0/Z_om
      y_xt   = real(IFFT(bitreverse(FFTshift(Y_hat)/Period)),lng)

     ! WANT TO USE term_ID TO PRINT FILES
!     write(*,*) "Ves:",term_ID,"Alpha max",alpha_max,"beta max",beta_max
!     write(*,*) "Zc:",real(Z_om(tmstps/2+1))+real(Z_om(tmstps/2+2))+real(Z_om(tmstps/2+3))+real(Z_om(tmstps/2+4)) &
!                        +real(Z_om(tmstps/2+5))
!      write(term_name,'(I3)') term_ID
!      fn = "STnetwork_" // adjustl(trim(term_name))
!        open (f1, file=fn, action='write')
!            write(f1,*) count
!            write(f1,*) alpha_max
!            write(f1,*) beta_max
!        do k = 1, tmstps
!            write(f1,*) Z_om(k)
!            write(f1,*) real(Z_om(k))
!        end do
!        close(f1)

    else ! DONT RECOMPUTE IMPEDANCE
 !       fn = "STnetwork_" // adjustl(trim(term_name))
 !       open (f1, file=fn)
 !           read(f1,*) count
 !           read(f1,*) alpha_max
 !           read(f1,*) beta_max
 !       do k = 1, tmstps
 !           read(f1,*) Zin
 !           Z_om(k) = Zin
!            write(*,*) Zin
 !       end do
 !       close(f1)
    Y_hat = 1.0/Z_om
        y_xt   = real(IFFT(bitreverse(FFTshift(Y_hat)/Period)),lng)
    endif


  return
  
end subroutine impedance


subroutine impedance_init (tmstps)

  integer,   intent(in)  :: tmstps
 
  allocate(Z_om(tmstps))
  
end subroutine impedance_init


subroutine impedance_close 

  deallocate(Z_om)
    
end subroutine impedance_close

end module root_imp



