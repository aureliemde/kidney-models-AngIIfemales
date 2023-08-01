subroutine jacobi2(fcn,n,x,fvec,fjac,ldfjac,iflag,epsfcn,numpar,pars,pt,id) 

include 'values.h'
include 'defs.h'

! passed variables
type (membrane) :: pt
integer n,ldfjac,iflag,numpar,id
double precision epsfcn
double precision x(n),fvec(n),fjac(ldfjac,n),pars(numpar)

!     **********
!
!     subroutine fdjac1
!
!     this subroutine computes a forward-difference approximation
!     to the n by n jacobian matrix associated with a specified
!     problem of n functions in n variables. if the jacobian has
!     a banded form, then function evaluations are saved by only
!     approximating the nonzero terms.
!
!     the subroutine statement is
!
!       subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,epsfcn) !
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(n,x,fvec,iflag)
!         integer n,iflag
!         double precision x(n),fvec(n)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of fdjac1.
!         in this case set iflag to a negative integer.
!
!       n is a positive integer input variable set to the number
!         of functions and variables.
!
!       x is an input array of length n.
!
!       fvec is an input array of length n which must contain the
!         functions evaluated at x.
!
!       fjac is an output n by n array which contains the
!         approximation to the jacobian matrix evaluated at x.
!
!       ldfjac is a positive integer input variable not less than n
!         which specifies the leading dimension of the array fjac.
!
!       iflag is an integer variable which can be used to terminate
!         the execution of fdjac1. see description of fcn.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!     subprograms called
!
!       minpack-supplied ... dpmpar
!
!       fortran-supplied ... dabs,dmax1,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
!

! local variables
double precision wa1(n),wa2(n)
integer i,j,k,msum
double precision eps,epsmch,h,temp,zero
double precision dpmpar
data zero /0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
      eps = dsqrt(dmax1(epsfcn,epsmch))
!
!        computation of dense approximate jacobian.
!	   
         do j = 1, n
            temp = x(j)
            h = eps*dabs(temp)
            if (h == zero) h = eps
            x(j) = temp + h
            call fcn(n,x,wa1,iflag,numpar,pars,pt,id)
!print*,"j",j,fvec(1)
!           if (iflag < 0) go to 30
            x(j) = temp
            do i = 1, n
               fjac(i,j) = (wa1(i) - fvec(i))/h
          end do
        end do

      return
      end

