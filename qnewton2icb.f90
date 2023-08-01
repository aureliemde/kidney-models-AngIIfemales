!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!
!  *************************** NEWTON SOLVER *************************
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     This is the Newton solver used to determine the luminal and epithelial 
!	concentrations, volumes, and EP below the inlet 
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

subroutine qnewton2icb (tubeprev,tube,Lz,id,Vol0,VolEinit,VolPinit, &
						VolAinit,VolBinit,ND)
      
  include 'values.h'
  include 'global.h'
  include 'defs.h'
  
  external fcn2C, fcn2OMC

  ! passed variables
  type (membrane) :: tubeprev,tube
  integer ND
  double precision Vol0(NC), VolEinit,VolPinit
  integer Lz,id
  
  double precision Cb(NS,NC),Volb(NC),EPb(NC),phb(NC)
  double precision AVF(NDC),A(NDC,NDC)
  double precision Cprev(NS,NC),Volprev(NC),EPprev(NC),phprev(NC)
  double precision osmol(NC)
  
  integer numpar
!  parameter (numpar=5*NS+16)
  parameter (numpar=1+(NC-1)*(NS+2))
  double precision pars(numpar)

  !      for fjac subroutine
  parameter (num=11+5*NS2)
  integer iflag,ml,mu
  double precision epsfcn
  double precision x(num),fvec(num),fjac(num,num),wa1(num),wa2(num)
  
  !	for matrix inversion
  integer info, lwork
  parameter (lwork=10000)
  integer ipiv(NDC)
  double precision work(lwork)
  
  TOL = 1.0d-5
!  TOL = 1.0d-3  ! when we set TOL 0.001, IMCD doesn't converge...

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     INITIAL GUESS: ASSUME SAME VALUES AS AT PREVIOUS STEP
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     store initial concentrations
!	CPREV's are used in computing iteration residue
!---------------------------------------------------------------------72

  Do K = 1, NC-1
     Do I = 1, NS
        Cb(I,K)=tubeprev%conc(I,K)
        Cprev(I,K) = Cb(I,K)
     End Do
     
     phb(K)=tubeprev%ph(K)
     phprev(K)=phb(K)
     
     Volb(K)=tubeprev%vol(K)
     Volprev(K)=Volb(K)
     
     EPb(K)=tubeprev%ep(K)
     EPprev(K)=EPb(K)     
  End Do
         
  PMb=tubeprev%pres !PMz(Lz)
  PMprev=PMb

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     START ITERATION LOOP TO DETERMINE SOLUTION
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
  
  nwrite = 0
  res = 1.0d0
  iter = 0
  
25 Do While ( res .gt. TOL ) 
     
     res = 0.0d0   
     iter = iter + 1
!---------------------------------------------------------------------72
!	Create vector x of variables to solve for
!---------------------------------------------------------------------72

     Do I = 1,NS2
        x(1+5*(I-1))=Cb(I,1)
        x(2+5*(I-1))=Cb(I,2)
        x(3+5*(I-1))=Cb(I,3)
        x(4+5*(I-1))=Cb(I,4)
        x(5+5*(I-1))=Cb(I,5)
     End do
     Do J=1,NC-1
        x(J+5*NS2)=Volb(J)
        x(J+5+5*NS2)=EPb(J)
     End do
     x(11+5*NS2)=PMb
	
!---------------------------------------------------------------------72
!	determine vector fvec of non-linear equations to be solved
!	determine Jacobian matrix fjac corresponding to fvec
!	fjac is a (NUxNU) matrix
!---------------------------------------------------------------------72

     ldfjac = num
     iflag = 1
     ml = num
     mu = num
     epsfcn = 1.d-6
     
     do K = 1, NC-1
        do I = 1, NS
           pars(I+(K-1)*NS) = tubeprev%conc(I,K)
        end do
     end do
     do K = 1, NC-1
        pars(K+(NC-1)*NS) = tubeprev%vol(K)
     end do
     do K = 1, NC-1
        pars(K+(NC-1)*(NS+1)) = tubeprev%ep(K)
     end do
     pars(1+(NC-1)*(NS+2)) = tubeprev%pres
     
     if(id==6 .or. id==7) then
        call fcn2C(num,x,fvec,iflag,numpar,pars,tube,id)
        
        call jacobi2(fcn2C,num,x,fvec,fjac,ldfjac,iflag,epsfcn,numpar,pars,tube,id)

     elseif(id==8) then
        call fcn2OMC(num,x,fvec,iflag,numpar,pars,tube,id)
        
        call jacobi2(fcn2OMC,num,x,fvec,fjac,ldfjac,iflag,epsfcn,numpar,pars,tube,id)
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
        if (id==6) then  ! 6 = CCD, 7 = OMCD
          Do M = 1, num
             AVF(L) = AVF(L) + 0.250*fjac(L,M)*fvec(M)
          End Do
        else
          Do M = 1, num
             AVF(L) = AVF(L) + 0.250*fjac(L,M)*fvec(M)
          End Do
        end if
     End do
!---------------------------------------------------------------------72
!     Calculate updated concentrations, volumes, and potentials
!	Then update "old" values
!---------------------------------------------------------------------72

     Do I = 1,NS2
        Cb(I,1) = Cprev(I,1)-AVF(1+5*(I-1))
        Cb(I,2) = Cprev(I,2)-AVF(2+5*(I-1))
        Cb(I,3) = Cprev(I,3)-AVF(3+5*(I-1))
        Cb(I,4) = Cprev(I,4)-AVF(4+5*(I-1))
        Cb(I,5) = Cprev(I,5)-AVF(5+5*(I-1))
                   
		if (I.eq.14) then
					
			res = dmax1(res,dabs(Cb(I,1)/Cprev(I,1)-1.d0)*0.01)
			res = dmax1(res,dabs(Cb(I,2)/Cprev(I,2)-1.d0)*0.01)
			res = dmax1(res,dabs(Cb(I,3)/Cprev(I,3)-1.d0)*0.01)
			res = dmax1(res,dabs(Cb(I,4)/Cprev(I,4)-1.d0)*0.01)
			res = dmax1(res,dabs(Cb(I,5)/Cprev(I,5)-1.d0)*0.01)

		else 

			res = dmax1(res,dabs(Cb(I,1)/Cprev(I,1)-1.d0))
			res = dmax1(res,dabs(Cb(I,2)/Cprev(I,2)-1.d0))
			res = dmax1(res,dabs(Cb(I,3)/Cprev(I,3)-1.d0))
			res = dmax1(res,dabs(Cb(I,4)/Cprev(I,4)-1.d0))
			res = dmax1(res,dabs(Cb(I,5)/Cprev(I,5)-1.d0))
                
		end if

        do K=1,5
           if (Cb(I,K).le.0.0d0) then
              print*,"warning in newton2icb",I,K,Cb(I,K)
           end if
        end do

     End Do

     Do J = 1,5
        phb(J) = -dlog(Cb(12,J)/1.0d3)/dlog(10.0d0)
        volb(J) = Volprev(J)-AVF(J+5*NS2)
        res = dmax1(res,dabs(Volb(J)-Volprev(J)))
        EPb(J) = EPprev(J)-AVF(J+5+5*NS2)
        res = dmax1(res,1.d-3*dabs(EPb(J)-EPprev(J)))
     End Do
     
     PMb = PMprev-AVF(11+5*NS2)
     res = dmax1(res,dabs(PMb-PMprev))

!	update "previous" concentrations, volumes, and potentials

     Do K = 1,5
        Do I = 1, NS
           Cprev(I,K) = Cb(I,K)
        End Do
        phprev(K) = phb(K)
        Volprev(K) = Volb(K)
        EPprev(K) = EPb(K)
     End Do
     PMprev=PMb

     if (iter >= 11) then

        if (iter < 21) then
            TOL = 1.0d0-4
 !           write(*,'(a,4g12.4)'),"Qnewton2icb - iter < 21",Lz+1,res
        else if (iter < 31) then
            TOL = 1.0d-3
 !           write(*,'(a,4g12.4)'),"Qnewton2icb - iter < 31",Lz+1,res
        else
            TOL = res*1.010
            write(*,'(a,4g12.4)'),"Qnewton2icb - iter > 31",Lz+1,res,TOL
        endif

        goto 25
     end if
     
  End Do

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	 Update concentrations at Lz + 1
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  Do K = 1,5
     Do I = 1, NS
        tube%conc(I,K) = Cb(I,K)
     End Do
     tube%ph(K) = phb(K)
     tube%vol(K) = Volb(K)
     tube%ep(K) = EPb(K)
  End Do
  tube%pres = PMb


!---------------------------------------------------------------------72
!	 print to output file
!---------------------------------------------------------------------72

50 if (Lz==Nz-1) then

  if (id==6) then
     open ( unit=42, file='CNToutlet' )
  elseif (id==7) then
     open ( unit=42, file='CCDoutlet' )
  elseif (id==8) then
     open ( unit=42, file='OMCoutlet' )
  end if
  
  Do I = 1,4
     write(42,200),tube%conc(I,1),tube%conc(I,6)
  end Do
  Do I = 5,NS
     write(42,210),tube%conc(I,1),tube%conc(I,6)
  end Do
  write(42,200),tube%ph(1),tube%ph(6)
  write(42,200),tube%ep(1),tube%ep(6)
  write(42,210),tube%vol(1),tube%pres
  close ( unit=42 )
 

!---------------------------------------------------------------------72
!	 print to output file
!---------------------------------------------------------------------72
  
  if (id==6) then
     open ( unit=44, file='CNToutlet_all' )
  elseif (id==7) then
     open ( unit=44, file='CCDoutlet_all' )
  elseif (id==8) then
     open ( unit=44, file='OMCoutlet_all' )
  end if
  
  Do I = 1,4
     write(44,200),tube%conc(I,2),tube%conc(I,3),tube%conc(I,4),tube%conc(I,5)
  end Do
  Do I = 5,NS
     write(44,210),tube%conc(I,2),tube%conc(I,3),tube%conc(I,4),tube%conc(I,5)
  end Do
  write(44,200),tube%ph(2),tube%ph(3),tube%ph(4),tube%ph(5)
  write(44,200),tube%ep(2),tube%ep(3),tube%ep(4),tube%ep(5)
  write(44,210),tube%vol(2),tube%vol(3),tube%vol(4),tube%vol(5)
  close ( unit=44 )
  
  endif

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
end subroutine qnewton2icb


