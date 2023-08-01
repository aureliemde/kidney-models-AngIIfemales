!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!
!  *************************** NEWTON SOLVER *************************
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     This is the Newton solver used to determine the epithelial 
!	concentrations, volumes, and EP at the inlet (jz = 0)

subroutine qnewton1icb (tube,CPimpref,CAimpref,CBimpref,id)
      
  include 'values.h'
  include 'global.h'
  include 'defs.h'
  
  external fcn1C,fcn1OMC

  ! passed variables
  type (membrane) :: tube
  double precision CPimpref, CAimpref, CBimpref
  integer id

! local variables
  double precision AVF(NUC),A(NUC,NUC)
  double precision Cprev(NS,NC),Volprev(NC),EPprev(NC),phprev(NC)
  double precision osmol(NC)
  double precision theta(NC),Slum(NC),Slat(NC),Sbas(NC)

  integer numpar
  parameter (numpar=2*NS+8)
  double precision pars(numpar)
  
!      for fjac subroutine
	 !parameter (num=74)
  parameter (num=9+4*NS2)
  integer iflag,ml,mu
  double precision epsfcn
  double precision x(num),fvec(num),fjac(num,num),wa1(num),wa2(num)
  double precision fjacold(num,num),fjacinv(num,num)

!	for matrix inversion
  integer info, lwork
  parameter (lwork=10000)
  integer ipiv(NUC)
  double precision work(lwork)

!---------------------------------------------------------------------72
!     private variables
!---------------------------------------------------------------------72

  TOL = 1.0d-4 !!!!!!!!!! MODIFIED (OLD VALUE = 1.0d-3) !!!!

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     INITIAL GUESSES FOR CONCENTRATIONS, VOLUMES, AND POTENTIALS
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72


!---------------------------------------------------------------------72
!	From restart file
!---------------------------------------------------------------------72

  if (id==5) then
     open ( unit=61, file='CNTresults')
  elseif (id==6) then
     open ( unit=61, file='CCDresults')
  elseif (id==7) then
     open ( unit=61, file='OMCresults')
  end if
  Do I = 1,4
     read(61,200),Ci1,tube%conc(I,2),tube%conc(I,3),tube%conc(I,4),tube%conc(I,5)
  end Do
  Do I = 5,NS
     read(61,210),Ci1,tube%conc(I,2),tube%conc(I,3),tube%conc(I,4),tube%conc(I,5)
  end Do
  read(61,200),phi1,tube%ph(2),tube%ph(3),tube%ph(4),tube%ph(5)
  read(61,200),Vol1,tube%vol(2),tube%vol(3),tube%vol(4),tube%vol(5)
  read(61,200),tube%ep(1),tube%ep(2),tube%ep(3),tube%ep(4),tube%ep(5)
  close ( unit=61 )

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     store initial concentrations
!	CPREV's are used in computing iteration residue
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
       
  Do K = 1, NC
     Do I = 1, NS
        Cprev(I,K) = tube%conc(I,K)
     End Do
     phprev(K)=tube%ph(K)
     Volprev(K)=tube%vol(K)
     EPprev(K)=tube%ep(K)
  End Do
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     START ITERATION LOOP TO FIND SOLUTION
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
         
  res = 1.0d0
  iter = 0
         
25  Do While ( res .gt. TOL ) 
     
     res = 0.0d0   
     iter = iter + 1
     
!---------------------------------------------------------------------72
!	Create vector x of variables to solve for
!---------------------------------------------------------------------72

     Do I = 1,NS2
        x(1+4*(I-1))=tube%conc(I,2)
        x(2+4*(I-1))=tube%conc(I,3)
        x(3+4*(I-1))=tube%conc(I,4)
        x(4+4*(I-1))=tube%conc(I,5)
     End do
     Do J=2,5
        x(J-1+4*NS2)=tube%vol(J)
        x(J+3+4*NS2)=tube%ep(J)
     End do
     x(9+4*NS2)=tube%ep(1)

!---------------------------------------------------------------------72
!	determine vector fvec of non-linear equations to be solved
!	determine Jacobian matrix fjac corresponding to fvec
!	fjac is a (NUxNU) matrix
!---------------------------------------------------------------------72

     ldfjac = num
     iflag = 1
     ml = num
     mu = num
     epsfcn = 1.d-5

     if(id==5 .or. id==6) then
        call fcn1C(num,x,fvec,iflag,numpar,pars,tube,id)
        
        call jacobi2(fcn1C,num,x,fvec,fjac,ldfjac,iflag,epsfcn,numpar,pars,tube,id)
     elseif (id==7) then
        call fcn1OMC(num,x,fvec,iflag,numpar,pars,tube,id)
        
        call jacobi2(fcn1OMC,num,x,fvec,fjac,ldfjac,iflag,epsfcn,numpar,pars,tube,id)
     end if



!---------------------------------------------------------------------72
!     invert Jacobian matrix
!---------------------------------------------------------------------72

     call dgetrf( num, num, fjac, num, ipiv, info )
     
     call dgetri( num, fjac, num, ipiv, work, lwork, info )

!---------------------------------------------------------------------72
!     Calculate AVF: product of fjac  (inverted Jacobian) and fvec
!---------------------------------------------------------------------72
     
     Do L = 1, num
        AVF(L)=0.d0
        Do M = 1, num
           AVF(L) = AVF(L) + 0.50*fjac(L,M)*fvec(M)
        End Do
     End do
     

!---------------------------------------------------------------------72
!     Calculate updated concentrations, volumes, and potentials
!	Then update "old" values
!---------------------------------------------------------------------72


     Do I = 1,NS2-1
        tube%conc(I,2) = Cprev(I,2)-AVF(1+4*(I-1))
        tube%conc(I,3) = Cprev(I,3)-AVF(2+4*(I-1))
        tube%conc(I,4) = Cprev(I,4)-AVF(3+4*(I-1))
        tube%conc(I,5) = Cprev(I,5)-AVF(4+4*(I-1))
        res = dmax1(res,dabs(tube%conc(I,2)/Cprev(I,2)-1.d0))
        res = dmax1(res,dabs(tube%conc(I,3)/Cprev(I,3)-1.d0))
        res = dmax1(res,dabs(tube%conc(I,4)/Cprev(I,4)-1.d0))
        res = dmax1(res,dabs(tube%conc(I,5)/Cprev(I,5)-1.d0))

        do K=2,5
           if (tube%conc(I,K).le.0.0d0) then
              print*,"negative concentration in qnewton1CD",K,I,tube%conc(I,K)
           end if
        end do
     End Do

     Do J = 2,5
     tube%vol(J) = Volprev(J)-AVF(J-1+4*NS2)
     res = dmax1(res,dabs(tube%vol(J)-Volprev(J)))
     tube%ep(J) = EPprev(J)-AVF(J+3+4*NS2)
     res = dmax1(res,dabs(tube%ep(J)-EPprev(J)))

  End Do

  tube%ep(1) = EPprev(1)-AVF(9+4*NS2)
  res = dmax1(res,dabs(tube%ep(1)-EPprev(1))*0.1d0)


!	update "previous" concentrations, volumes, and potentials

  Do K = 2,5
     Do I = 1, NS
        Cprev(I,K) = tube%conc(I,K)
     End Do
     tube%ph(K) = -dlog(tube%conc(12,K)/1.d3)/dlog(10.0d0)
     phprev(K) = tube%ph(K)
     Volprev(K) = tube%vol(K)
     EPprev(K) = tube%ep(K)
  End Do
  EPprev(1)=tube%ep(1)
	 
  print *, "iter=", iter, "   res=", res
  print*

  if (iter >= 10) then
     write(*,'(a,4g12.4)')," Too many iterations ",res
     print*
     if (iter < 20) then
        TOL = 1.0d-3
     else
        TOL = res*1.010
     endif
     goto 25
  end if
  
End Do

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	                           FINAL RESULTS
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72


!---------------------------------------------------------------------72
!	 print to output file
!---------------------------------------------------------------------72


es=tube%ep(6)

500 if (id==5) then
   open ( unit=61, file='CNTresults' )
elseif (id==6) then
   open ( unit=61, file='CCDresults' )
elseif (id==7) then
   open ( unit=61, file='OMCresults' )
end if

Do I = 1,4
   write(61,200),tube%conc(I,1),tube%conc(I,2),tube%conc(I,3), tube%conc(I,4),tube%conc(I,5)
end Do
Do I = 5,NS
   write(61,210),tube%conc(I,1),tube%conc(I,2),tube%conc(I,3), tube%conc(I,4),tube%conc(I,5)
end Do
write(61,200),tube%ph(1),tube%ph(2),tube%ph(3),tube%ph(4),tube%ph(5)
write(61,200),tube%vol(1),tube%vol(2),tube%vol(3),tube%vol(4), tube%vol(5)
write(61,200),tube%ep(1)-es,tube%ep(2)-es,tube%ep(3)-es, tube%ep(4)-es,tube%ep(5)-es

close ( unit=61 )
 
!---------------------------------------------------------------------72
!	 Determine fluxes
!---------------------------------------------------------------------72

Do I = 1,NS2
   x(1+4*(I-1))=tube%conc(I,2)
   x(2+4*(I-1))=tube%conc(I,3)
   x(3+4*(I-1))=tube%conc(I,4)
   x(4+4*(I-1))=tube%conc(I,5)
End do
Do J=2,5
   x(J-1+4*NS2)=tube%vol(J)
   x(J+3+4*NS2)=tube%ep(J)
End do
x(9+4*NS2)=tube%ep(1)


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

200 format (6f12.5)
210 format (6e12.5)
300 format (A21)
305 format (A2)
310 format (A6)
320 format (A32)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

return
end subroutine qnewton1icb


