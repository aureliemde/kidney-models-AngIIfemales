!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
 
!  *************************** NEWTON SOLVER *************************
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     This is the Newton solver used to determine the luminal and epithelial 
!	concentrations, volumes, and EP below the inlet 
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

subroutine qnewton2b (tubeprev,tube,Vol0,Lz,id,ND)

  include 'values.h'
  include 'global.h'
  include 'defs.h'

  external fcn2b

  ! passed variables
  type (membrane) :: tubeprev, tube
  integer ND
  double precision Vol0(NC)
  integer Lz,id

  double precision Cb(NS,NC),Volb(NC),EPb(NC),phb(NC)
  double precision AVF(ND),A(ND,ND)
  double precision Cprev(NS,NC),Volprev(NC),EPprev(NC),phprev(NC)
  double precision osmol(NC)


!      for fjac subroutine
  parameter (num=7+3*NS2)
  parameter (numpar=8+5*NS)
  integer iflag,ml,mu
  double precision epsfcn
  double precision x(num),fvec(num),fjac(num,num),wa1(num),wa2(num)
  double precision pars(numpar)

!	for matrix inversion
  integer info, lwork
  parameter (lwork=10000)
  integer ipiv(ND)
  double precision work(lwork)

  TOL = 1.0d-5

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
         
	 PMb=tubeprev%pres
	 PMprev=PMb


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     START ITERATION LOOP TO DETERMINE SOLUTION
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
  
         res = 1.0d0
         iter = 0
         
25         Do While ( res > TOL ) 
            
            res = 0.0d0   
            iter = iter + 1
!---------------------------------------------------------------------72
!	Create vector x of variables to solve for
!---------------------------------------------------------------------72

	 Do I = 1,NS2
		x(1+3*(I-1))=Cb(I,1)
		x(2+3*(I-1))=Cb(I,2)
		x(3+3*(I-1))=Cb(I,5)
	 End do
	 x(1+3*NS2) = Volb(1)
	 x(2+3*NS2) = Volb(2)
	 x(3+3*NS2) = Volb(5)
	 x(4+3*NS2) = EPb(1)
	 x(5+3*NS2) = EPb(2)
	 x(6+3*NS2) = EPb(5)
	 x(7+3*NS2) = PMb
	


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

         pars(1:NS) = tubeprev%conc(:,1)	
         pars(1+NS:2*NS) = tubeprev%conc(:,2)	
         pars(1+2*NS:3*NS) = tubeprev%conc(:,5)	
         pars(1+3*NS) = tubeprev%vol(1)
         pars(2+3*NS) = tubeprev%vol(2)
         pars(3+3*NS) = tubeprev%vol(5)
         pars(4+3*NS) = tubeprev%pres
         pars(5+3*NS:4+4*NS) = tubeprev%conc(:,6)
         pars(5+4*NS:4+5*NS) = tube%conc(:,6)
         pars(5+5*NS) = tube%ep(6)
         if (id==3) then ! mTAL
           pars(6+5*NS) = dimLA
           pars(7+5*NS) = CPimprefA
           pars(8+5*NS) = DiamA
        elseif (id==4) then  ! cTAL
           pars(6+5*NS) = dimLT
           pars(7+5*NS) = CPimprefT
           pars(8+5*NS) = DiamT
         elseif (id==5) then ! DCT
           pars(6+5*NS) = dimLD
           pars(7+5*NS) = CPimprefD
           pars(8+5*NS) = DiamD
         elseif (id==9) then ! IMCD
           pars(6+5*NS) = dimLIMC
           pars(7+5*NS) = CPimprefIMC
           pars(8+5*NS) = DiamIMC
         end if
 

	  call fcn2b(num,x,fvec,iflag,numpar,pars,tube,id)

	  call jacobi2(fcn2b,num,x,fvec,fjac,ldfjac,iflag,epsfcn,numpar,pars,tube,id)

!---------------------------------------------------------------------72
!     invert Jacobian matrix
!---------------------------------------------------------------------72


	 call dgetrf( num, num, fjac, num, ipiv, info )

     call dgetri( num, fjac, num, ipiv, work, lwork, info )

!---------------------------------------------------------------------72
!     Calculate AVF: product of fjac (inverted Jacobian) and fvec
!---------------------------------------------------------------------72

	 Do L = 1, num
		AVF(L)=0.d0
		Do M = 1, num
		    if (id == 5) then
			AVF(L) = AVF(L) + 0.250*fjac(L,M)*fvec(M) ! original = 0.500
			else
			AVF(L) = AVF(L) + 0.50*fjac(L,M)*fvec(M) ! original = 0.250 or 0.500
			end if
		End Do
	 End do

!---------------------------------------------------------------------72
!     Calculate updated concentrations, volumes, and potentials
!	Then update "old" values
!---------------------------------------------------------------------72


	 Do I = 1,NS2
			Cb(I,1) = Cprev(I,1)-AVF(1+3*(I-1))
			Cb(I,2) = Cprev(I,2)-AVF(2+3*(I-1))
			Cb(I,5) = Cprev(I,5)-AVF(3+3*(I-1))

			if (I.eq.14) then
					
			res = dmax1(res,dabs(Cb(I,1)/Cprev(I,1)-1.d0)*0.01)
			res = dmax1(res,dabs(Cb(I,2)/Cprev(I,2)-1.d0)*0.01)
			res = dmax1(res,dabs(Cb(I,5)/Cprev(I,5)-1.d0)*0.01)

			else

			res = dmax1(res,dabs(Cb(I,1)/Cprev(I,1)-1.d0))
			res = dmax1(res,dabs(Cb(I,2)/Cprev(I,2)-1.d0))
			res = dmax1(res,dabs(Cb(I,5)/Cprev(I,5)-1.d0))

			endif

			if (Cb(I,2).le.0.0d0 .or. Cb(I,5).le.0.0d0) then
				print*,"warning in newton2b",I,Cb(I,2),Cb(I,5) 
			end if


	 End Do
	 
	 phb(1) = -dlog(Cb(12,1)/1.0d3)/dlog(10.0d0)
	 phb(2) = -dlog(Cb(12,2)/1.0d3)/dlog(10.0d0)
	 phb(5) = -dlog(Cb(12,5)/1.0d3)/dlog(10.0d0)

	 volb(1) = Volprev(1)-AVF(1+3*NS2)
	 res = dmax1(res,dabs(Volb(1)-Volprev(1)))

	 volb(2) = Volprev(2)-AVF(2+3*NS2)
	 res = dmax1(res,dabs(Volb(2)-Volprev(2)))

	 volb(5) = Volprev(5)-AVF(3+3*NS2)
	 res = dmax1(res,dabs(Volb(5)-Volprev(5)))
		
	 EPb(1) = EPprev(1)-AVF(4+3*NS2)
	 res = dmax1(res,1.d-3*dabs(EPb(1)-EPprev(1)))
	 EPb(2) = EPprev(2)-AVF(5+3*NS2)
	 res = dmax1(res,1.d-3*dabs(EPb(2)-EPprev(2)))
	 EPb(5) = EPprev(5)-AVF(6+3*NS2)
	 res = dmax1(res,1.d-3*dabs(EPb(5)-EPprev(5)))
	 

	 PMb = PMprev-AVF(7+3*NS2)
	 res = dmax1(res,dabs(PMb-PMprev))


!	update "previous" concentrations, volumes, and potentials

	  Do I = 1, NS
			Cprev(I,1) = Cb(I,1)
			Cprev(I,2) = Cb(I,2)
			Cprev(I,5) = Cb(I,5)
	  End Do

	  phprev(1) = phb(1)
	  Volprev(1) = Volb(1)
	  EPprev(1) = EPb(1)

	  phprev(2) = phb(2)
	  Volprev(2) = Volb(2)
	  EPprev(2) = EPb(2)

	  phprev(5) = phb(5)
	  Volprev(5) = Volb(5)
	  EPprev(5) = EPb(5)

	  PMprev=PMb

     if (iter >= 11) then

        if (iter < 21) then
            TOL = 1.0d0-4 ! 1.0d-3
!            write(*,'(a,4g12.4)'),"Qnewton2b - iter < 21",Lz+1,res
        else if (iter < 31) then
            TOL = 1.0d-3
!            write(*,'(a,4g12.4)'),"Qnewton2b - iter < 31",Lz+1,res
        else
            TOL = res*1.010
            write(*,'(a,4g12.4)'),"Qnewton2b - iter > 31",Lz+1,res,TOL
        endif

        goto 25
     end if

       End Do

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	 Update concentrations at Lz + 1
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72


	 Do I = 1, NS
			tube%conc(I,1) = Cb(I,1)
			tube%conc(I,2) = Cb(I,2)
			tube%conc(I,5) = Cb(I,5)
	 End Do

	 tube%ph(1) = phb(1)
	 tube%vol(1) = Volb(1)
	 tube%ep(1) = EPb(1)

	 tube%ph(2) = phb(2)
	 tube%vol(2) = Volb(2)
	 tube%ep(2) = EPb(2)

	 tube%ph(5) = phb(5)
	 tube%vol(5) = Volb(5)
	 tube%ep(5) = EPb(5)

	 tube%pres = PMb


!---------------------------------------------------------------------72
!	 print luminal and interstitial results to output file
!---------------------------------------------------------------------72
   if (Lz==NZ-1) then

      if (id==2) then ! SDL
        open ( unit=12, file='SDLoutlet' )
      elseif (id==3) then  ! mTAL
         open ( unit=12, file='mTALoutlet' )
       elseif (id==4) then  ! cTAL
         open ( unit=12, file='cTALoutlet' )
       elseif (id==5) then  ! DCT
         open ( unit=12, file='DCToutlet' )
       elseif (id==9) then  ! IMCD
         open ( unit=12, file='IMCoutlet' )
       end if

	 Do I = 1,4
			write(12,200),tube%conc(I,1),tube%conc(I,6)
	 end Do
	 Do I = 5,NS
			write(12,210),tube%conc(I,1),tube%conc(I,6)
	 end Do
	 write(12,200),tube%ph(1),tube%ph(6)
	 write(12,200),tube%ep(1),tube%ep(6)
	 write(12,210),tube%vol(1),tube%pres

       close ( unit=12 )

!---------------------------------------------------------------------72
!	 print cell and LIS results to output file
!---------------------------------------------------------------------72

       if (id==2) then !SDL
         open ( unit=14, file='SDLoutlet_all' )
       else if (id==3) then !  mTAL
         open ( unit=14, file='mTALoutlet_all' )
       elseif (id==4) then ! cTAL
         open ( unit=14, file='cTALoutlet_all' )
       elseif (id==5) then ! DCT
         open ( unit=14, file='DCToutlet_all' )
       elseif (id==9) then ! IMCD
         open ( unit=14, file='IMCoutlet_all' )
       end if

	 Do I = 1,4
			write(14,200),tube%conc(I,2),tube%conc(I,5)
	 end Do
	 Do I = 5,NS
			write(14,210),tube%conc(I,2),tube%conc(I,5)
	 end Do
	 write(14,200),tube%ph(2),tube%ph(5)
	 write(14,200),tube%ep(2),tube%ep(5)
	 write(14,210),tube%vol(2),tube%vol(5)

       close ( unit=14 )

  end if

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

 200   format (6f12.5)
 210   format (6e12.5)
 300   format (A21)
 305   format (A2)
 310   format (A6)
 320   format (A32)

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

1000	 return
       end


