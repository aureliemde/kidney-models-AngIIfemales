!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
 
!  *************************** NEWTON SOLVER *************************
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     This is the Newton solver used to determine the luminal and epithelial 
!	concentrations, volumes, and EP below the inlet 
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

subroutine qnewton2SDL (sdlprev,sdl,Lz,id,Vol0,ND)

  include 'values.h'
  include 'global.h'
  include 'defs.h'

  external fcn2SDL

  ! passed variables
  type (membrane) :: sdlprev, sdl
  integer ND
!  double precision PMz(0:NZ)
  double precision Vol0(NC), VolEinit,VolPinit
  integer Lz,id

  double precision Cb(NS,NC),Volb(NC),EPb(NC),phb(NC)
  double precision AVF(ND),A(ND,ND)
  double precision Cprev(NS,NC),Volprev(NC),EPprev(NC),phprev(NC)
  double precision osmol(NC)


!      for fjac subroutine
  parameter (num=2+NS2)
  parameter (numpar=8+4*NS)
  integer iflag,ml,mu
  double precision epsfcn
  double precision x(num),fvec(num),fjac(num,num),wa1(num),wa2(num)
  double precision pars(numpar)

!	for matrix inversion
  integer info, lwork
  parameter (lwork=10000)
  integer ipiv(ND)
  double precision work(lwork)

  TOL = 1.0d-5 ! 4.0d-4

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
        Cb(I,K)=sdlprev%conc(I,K)!   Cz(I,K,Lz)
        Cprev(I,K) = Cb(I,K)
     End Do
     
     phb(K)=sdlprev%ph(K)
     phprev(K)=phb(K)
     
     Volb(K)=sdlprev%vol(K)
     Volprev(K)=Volb(K)
     
     EPb(K)=sdlprev%ep(K)
     EPprev(K)=EPb(K)

  End Do
  
  PMb=sdlprev%pres!PMz(Lz)
  PMprev=PMb
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     START ITERATION LOOP TO DETERMINE SOLUTION
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
  
  res = 1.0d0
  iter = 0
         
  Do While ( res > TOL ) 
     
     res = 0.0d0   
     iter = iter + 1
!---------------------------------------------------------------------72
!	Create vector x of variables to solve for
!---------------------------------------------------------------------72

     Do I = 1,NS2
        x(I)=Cb(I,1)
     End do
     x(1+NS2) = Volb(1)
     x(2+NS2) = PMb
	
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

     pars(1:NS) = sdlprev%conc(:,1)  !Cz(:,1,Lz)	
     pars(1+NS:2*NS) = sdlprev%conc(:,2)  !Cz(:,2,Lz)	
     pars(1+2*NS:3*NS) = sdlprev%conc(:,5)  !Cz(:,5,Lz)	
     pars(1+3*NS) = sdlprev%vol(1)
     pars(2+3*NS) = sdlprev%vol(2)
     pars(3+3*NS) = sdlprev%vol(5)
     pars(4+3*NS) = sdlprev%pres!PMz(Lz)
     pars(5+3*NS:4+4*NS) = sdlprev%conc(:,6)  !Cz(:,6,Lz)
     pars(5+4*NS) = dimLSDL
     pars(6+4*NS) = CPimprefSDL
     pars(7+4*NS) = DiamSDL
     pars(8+4*NS) = Vol0(1)

     call fcn2SDL(num,x,fvec,iflag,numpar,pars,sdl,id)

     call jacobi2(fcn2SDL,num,x,fvec,fjac,ldfjac,iflag,epsfcn,numpar,pars,sdl,id)

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
           AVF(L) = AVF(L) + 1.0*fjac(L,M)*fvec(M) ! original = 0.200
        End Do
     End do

!---------------------------------------------------------------------72
!     Calculate updated concentrations, volumes, and potentials
!	Then update "old" values
!---------------------------------------------------------------------72

     Do I = 1,NS2
        Cb(I,1) = Cprev(I,1)-AVF(I)
        
        if (I.eq.14) then           
           res = dmax1(res,dabs(Cb(I,1)/Cprev(I,1)-1.d0)*0.1)
        else           
           res = dmax1(res,dabs(Cb(I,1)/Cprev(I,1)-1.d0))
        endif
        
     End Do
	 
     phb(1) = -dlog(Cb(12,1)/1.0d3)/dlog(10.0d0)
     volb(1) = Volprev(1)-AVF(1+NS2)
     res = dmax1(res,dabs(Volb(1)-Volprev(1)))
     
     PMb = PMprev-AVF(2+NS2)
     res = dmax1(res,dabs(PMb-PMprev))
     EPb(1) = EPprev(1)-AVF(4+3*NS2)
     res = dmax1(res,1.d-3*dabs(EPb(1)-EPprev(1)))

!	update "previous" concentrations, volumes, and potentials
  
     Do I = 1, NS
        Cprev(I,1) = Cb(I,1)
        Cprev(I,2) = Cb(I,2)
        Cprev(I,5) = Cb(I,5)
     End Do
     
     PMprev=PMb
     phprev(1) = phb(1)
     Volprev(1) = Volb(1)
     EPprev(1) = EPb(1)

       if (iter >= 22) then
        if (iter < 40) then
            TOL = 1.0d-3
  !          write(*,'(a,4g12.4)'),"Qnewton2SDL - iter > 21",Lz+1,res
        else
            TOL = res*1.010
            write(*,'(a,4g12.4)'),"Qnewton2SDL - iter > 40",Lz+1,res,TOL
        endif
     end if

  End Do

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!	 Update concentrations at Lz + 1
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

  Do I = 1, NS
     sdl%conc(I,1) = Cb(I,1)
     sdl%conc(I,2) = Cb(I,2)
     sdl%conc(I,5) = Cb(I,5)
  End Do
  
  sdl%vol(1) = Volb(1)
  sdl%pres = PMb
  sdl%ph(1) = -dlog(Cb(12,1)/1.0d3)/dlog(10.0d0)


  
!---------------------------------------------------------------------72
!	 print to output file
!---------------------------------------------------------------------72
  open ( unit=12, file='SDLoutlet' )
  
  Do I = 1,4
     write(12,200),sdl%conc(I,1),sdl%conc(I,6)
  end Do
  Do I = 5,NS
     write(12,210),sdl%conc(I,1),sdl%conc(I,6)
  end Do
  write(12,200),sdl%ph(1),sdl%ph(6)
  write(12,200),sdl%ep(1),sdl%ep(6)
  write(12,210),sdl%vol(1),sdl%pres  
  
  close ( unit=12 )
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

200 format (6f12.5)
210 format (6e12.5)
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72


  return
end subroutine qnewton2SDL


