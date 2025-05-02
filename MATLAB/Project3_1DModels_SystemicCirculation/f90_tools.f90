!  MODULE F90_TOOLS

!***************************************************************************
!*                                                                         *
!* Module: f90_tools.f90                                                   *
!* Version: 1.0                                                            *
!* By: Mette Olufsen, Math-Tech                                            *
!* Last modified: 14. Jan. 1997                                            *
!*                                                                         *
!* Contains subroutines and functions needed to compute the fast Fourier   *
!* transform and to flip a vector up down. These functions are all library *
!* subroutines and taken directly either from the corresponding Matlab     *
!* routines or in case of FFT from any standard textbook. I will therefore *
!* not make any specific comments of the separate routines.                *
!*                                                                         *
!***************************************************************************
module f90_tools
implicit none

integer,   parameter :: lng = selected_real_kind(12,99)
real(lng), parameter :: pi  = 3.14159265358979323_lng

contains

!***************************************************************************
!*                                                                         *
!* Function: flipud                                                        *
!*                                                                         *
!* Flips the contents of the vector Z up-down and thereby creating an mirror *
!* of the second half of the vector.                                       *
!* This function is equivalent to the Matlab function with the same name.  *
!*                                                                         *
!***************************************************************************
function flipud (Z)
implicit none

  complex(lng) :: Z(:), flipud(size(Z))
  integer      :: N

  N = size(Z)
  flipud(1:N:1) = Z(N:1:-1) 
end function flipud

!***************************************************************************
!*                                                                         *
!* Function: FFT                                                           *
!*                                                                         *
!* Fast Fourier transforms the input vector Z and returns Zhat.            *
!* This function is made according to a standard description in any text-  *
!* book.                                                                   *
!*                                                                         *
!***************************************************************************
function FFT (Z) result(Zhat)
implicit none
  
! Fast Fourier Transform. 
! Requires N data-points in natural order in Z, where we assume that
! N is a power of 2.
! The function returns the bit reversed transform in Zhat.

  complex(lng) :: Z(:), Zhat(size(Z)), wa(size(Z)/2), u, w, t
  integer      :: NV2, N, M, j, l, k, ns, le, i, index
    
  N    = size(Z)
  M    = nint(log(real(N))/log(2.0))
  NV2  = N/2
  Zhat = Z

  ! Check consistency of data.
  if (N /= 2**M) then
    print *,'N <> 2**M'
    return
  end if

  ! Generates a bit reversed sequence of twiddle factors, W**k.
  j = 1
  u = (1.0_lng,0.0_lng)
  w = cmplx(cos(pi/NV2),-sin(pi/NV2),lng)
  do  l = 1, NV2-1
    wa(j) = u
    u     = u*w
    k     = NV2/2
  
    do while (k < j) 
      j = j - k
      k = k/2
    end do
    j = j+k   !J-1 is a (M-1) bit reversing of l-1.
  end do
  wa(NV2) = u

  ! FFT algorithm. The input data should appear in natural order, the
  ! output will be a bit reversed transform.
  ns = NV2
  do l = 1, m
    le = 2**(l-1)
    index = 1
    do j = 1, le
      do i = 1, NS
        t = Zhat(index+ns)*wa(j)
        Zhat(index+ns) = Zhat(index) - t
        Zhat(index)    = Zhat(index) + t
        index = index + 1
      end do
      index = index + ns
    end do
    ns = ns/2
  end do
  u = (1.0_lng,0.0_lng)/N
  do i = 1, n
    Zhat(i) = Zhat(i)*u
  end do
  return
end function FFT

!***************************************************************************
!*                                                                         *
!* Function: IFFT                                                          *
!*                                                                         *
!* Inverse fast Fourier transforms the input vector Zhat and returns Z.    *
!* This function is made according to a standard description in any text-  *
!* book.                                                                   *
!*                                                                         *
!***************************************************************************
function IFFT (Zhat) result (Z)
implicit none

  ! Inverse FFT algorithm.
  ! The input data must be bit reversed, the resulting vector will appear in
  ! in natural order.

  complex(lng) :: Zhat(:), Z(size(Zhat)), u, w, t
  integer      :: N, m, l, le, le1, i, j, ip

  N = size(Zhat)
  M = nint(log(real(N))/log(2.0))
  Z = Zhat
  do  l = 1, M
    le  = 2**l
    le1 = le/2
    u   = (1.0_lng,0.0_lng)
    w   = cmplx(cos(pi/le1),sin(pi/le1),lng)
    do j = 1, le1
      do i = j, N, le
        ip    = i + le1
        t     = Z(ip)*u
        Z(ip) = Z(i) - t
        Z(i)  = Z(i) + t
      end do
      u = u*w
    end do
  end do
  return
end function IFFT

!***************************************************************************
!*                                                                         *
!* Function: FFTshift                                                      *
!*                                                                         *
!* Shifts the input vector Z by N/2 in a circular way since we shift by    *
!* N/2 the result is a swapping of the left and right half of the vector. *
!* This function is equivalent to the Matlab function with the same name.  *
!*                                                                         *
!***************************************************************************
function FFTshift (Z) result(Y)
implicit none

  complex(lng) :: Z(:), Y(size(Z))
  integer      :: N

  N = size(Z)
  Y = cshift(Z,N/2)  
end function FFTshift


!***************************************************************************
!*                                                                         *
!* Function: bitreverse                                                    *
!*                                                                         *
!* Since the FFT/IFFT respective outputs and requires as input a vector    *
!* which is bit reversed we need this routine. It takes a vector and bit   *
!* reverses it. Then having a normalized input we get a bit-reversed out-  *
!* and having a bit-reversed input we get a normalized output.             *
!* This function is made according the description of the output/input for *
!* FFT/IFFT.                                                               *
!*                                                                         *
!***************************************************************************
function bitreverse (X) result(Y)
implicit none

  complex(lng) :: X(:), Y(size(X))
  integer      :: j, k, l, N

  N = size(X)
  j = 1

  do l = 1, N-1
    Y(j) = X(l)
    k = N/2          ! Unit

    do               ! While k < j ripple carry to the right
      if (k >= j) exit  
      j = j-k        ! Reset bit
      k = k/2        ! Next carry
    end do
  
    j = j+k          ! Add final carry
  end do
  Y(N) = X(N)
end function bitreverse

end module f90_tools
