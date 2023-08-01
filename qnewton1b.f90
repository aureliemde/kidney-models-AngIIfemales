!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
 
!  *************************** NEWTON SOLVER *************************
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72


!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     This is the Newton solver used to determine the epithelial 
!	concentrations, volumes, and EP at the inlet (jz = 0)

subroutine qnewton1b (tube,CPimpref,id)

  include 'values.h'
  include 'global.h'
  include 'defs.h'
  
  external fcn1b

  !  passed variables
  type(membrane) :: tube
  double precision CPimpref
  integer id

  double precision AVF(NUA),A(NUA,NUA)
  double precision Cprev(NS,NC),Volprev(NC),EPprev(NC),phprev(NC)
  double precision osmol(NC)
  double precision theta(NC),Slum(NC),Slat(NC),Sbas(NC)
  
!      for fjac subroutine
	 !parameter (num=42)
  parameter (num=5+2*NS2)
  parameter (numpar=2*NS+2)
  integer iflag,ml,mu
  double precision epsfcn, pars(numpar)
  double precision x(NDA),fvec(num),fjac(num,num),wa1(num),wa2(num)
  double precision fjacold(num,num),fjacinv(num,num)
  
!	for matrix inversion
  integer info, lwork
  parameter (lwork=10000)
  integer ipiv(NUA)
  double precision work(lwork)

!---------------------------------------------------------------------72
!     private variables
!---------------------------------------------------------------------72

  TOL = 1.0d-4 !!!!!!!!!!! MODIFIED (OLD VALUE = 5.0d-4) !!!!

!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
!     INITIAL GUESSES FOR CONCENTRATIONS, VOLUMES, AND POTENTIALS
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72
        
        if (id==1) then
           open ( unit=11, file='mTALresults')
        elseif (id==2) then
           open ( unit=11, file='TALresults')
        elseif (id==3) then
           open ( unit=11, file='DCTresults')
        elseif (id==8) then
           open (unit=11, file='IMCDresults')
        elseif (id==-1) then
           open ( unit=11, file='SDLresultsAL')
        else
           open ( unit=11, file='MDresults')
        end if
        
        Do I = 1,4
           read(11,200),Ci1,tube%conc(I,2),tube%conc(I,5)
        end Do
        Do I = 5,NS
           read(11,210),Ci1,tube%conc(I,2),tube%conc(I,5)
        end Do
        read(11,200),phi1,tube%ph(2),tube%ph(5)
        read(11,200),Vol1,tube%vol(2),tube%vol(5)
        read(11,200),tube%ep(1),tube%ep(2),tube%ep(5)
        close ( unit=11 )
        
!    Force glucose concentrations to values near peritubular concentrations     
!	if (id == 3) then
!		tube%conc(15,2) = TotGluCM
!		tube%conc(15,5) = TotGluCM
!	end if
			
!	Assign arbitrary values to compartments A (4) and B (5)

	Do I = 1,NS
		tube%conc(I,3) = tube%conc(I,2)
		tube%conc(I,4) = tube%conc(I,2)
	End Do
  
	tube%ph(3) = tube%ph(2)
	tube%ph(4) = tube%ph(2)
  
	tube%vol(3) = 0 !VolAinitmtal
	tube%vol(4) = 0 !VolBinitmtal
  
	tube%ep(3) = tube%ep(2)
	tube%ep(4) = tube%ep(2)
  
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
!        
  res = 1.0d0
  iter = 0
         
25 Do While ( res .gt. TOL ) 
     
     res = 0.0d0   
     iter = iter + 1

!---------------------------------------------------------------------72
!	Create vector x of variables to solve for
!---------------------------------------------------------------------72

     Do I = 1,NS2
        x(1+2*(I-1))=tube%conc(I,2)
        x(2+2*(I-1))=tube%conc(I,5)
     End do
     x(1+2*NS2)=tube%vol(2)
     x(2+2*NS2)=tube%vol(5)
     x(3+2*NS2)=tube%ep(1)
     x(4+2*NS2)=tube%ep(2)
     x(5+2*NS2)=tube%ep(5)
     
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
	
     pars(1:NS) = tube%conc(1:NS,1)
     pars(1+NS:2*NS) = tube%conc(1:NS,6)
     pars(1+2*NS) = tube%ep(6)
     pars(2+2*NS) = CPimpref
      
!    print*,"x",x
     call fcn1b(num,x,fvec,iflag,numpar,pars,tube,id)
!    print*,"fvec",fvec
	 
     call jacobi2(fcn1b,num,x,fvec,fjac,ldfjac,iflag,epsfcn,numpar,pars,tube,id)

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
           AVF(L) = AVF(L) + fjac(L,M)*fvec(M)
        End Do
     End do

!---------------------------------------------------------------------72
!     Calculate updated concentrations, volumes, and potentials
!	Then update "old" values
!---------------------------------------------------------------------72

     Do I = 1,NS2
        tube%conc(I,2) = Cprev(I,2)-AVF(1+2*(I-1))
        tube%conc(I,5) = Cprev(I,5)-AVF(2+2*(I-1))
        res = dmax1(res,dabs(tube%conc(I,2)/Cprev(I,2)-1.d0))
        res = dmax1(res,dabs(tube%conc(I,5)/Cprev(I,5)-1.d0))

        if (tube%conc(I,2).le.0.0d0 .or. tube%conc(I,5).le.0.0d0) then
           print*,"warning in newton1b",I,tube%conc(I,2),tube%conc(I,5) 
        end if

!		if (I == 14) then 
!			print*,"res before glucose",res
!			else if (I == 15) then
!				print*,"res with glucose",res
!				print*,"Cell concentrations",tube%conc(I,2),Cprev(I,2)
!				print*,"LIS concentrations",tube%conc(I,5),Cprev(I,5)
!				pause
!			end if
!		end if

     End Do
     
     tube%vol(2) = Volprev(2)-AVF(1+2*NS2)
     res = dmax1(res,dabs(tube%vol(2)-Volprev(2)))
     tube%vol(5) = Volprev(5)-AVF(2+2*NS2)
     res = dmax1(res,dabs(tube%vol(5)-Volprev(5)))     
     
     tube%ep(1) = EPprev(1)-AVF(3+2*NS2)
     res = dmax1(res,dabs(tube%ep(1)-EPprev(1)))
     tube%ep(2) = EPprev(2)-AVF(4+2*NS2)
     res = dmax1(res,dabs(tube%ep(2)-EPprev(2)))
     
     tube%ep(5) = EPprev(5)-AVF(5+2*NS2)
     res = dmax1(res,dabs(tube%ep(5)-EPprev(5)))

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
        write(*,'(a,4g12.4)')," Too many iterations at inlet",res
        print*
        if (iter < 20) then
           TOL = 1.0d-3 !!!!!!!!!! MODIFIED (OLD VALUE = 1.0d-2) !!!!
        else
           TOL = res*1.010
           !	        pause
        endif
        goto 25
     end if
     
  End Do


!---------------------------------------------------------------------72
!	 print to output file
!---------------------------------------------------------------------72


  
 if (id==1) then
     open ( unit=11, file='mTALresults' )
  elseif (id==2) then
     open ( unit=11, file='TALresults' )
  elseif (id==3) then
     open ( unit=11, file='DCTresults' )
  elseif (id==4) then
     open ( unit=11, file='CNTresults' )
  elseif (id==8) then
     open ( unit=11, file='IMCDresults' )
  elseif (id==-1) then
     open ( unit=11, file='SDLresults' )
  else
     open ( unit=11, file='MDresults' )
  end if
  
  Do I = 1,4
     write(11,200),tube%conc(I,1),tube%conc(I,2),tube%conc(I,5)
  end Do
  Do I = 5,NS
     write(11,210),tube%conc(I,1),tube%conc(I,2),tube%conc(I,5)
  end Do
  write(11,200),tube%ph(1),tube%ph(2),tube%ph(5)
  write(11,200),tube%vol(1),tube%vol(2),tube%vol(5)
  write(11,200),tube%ep(1)-tube%ep(6),tube%ep(2)-tube%ep(6),  tube%ep(5)-tube%ep(6)
!  end if
  
  close ( unit=11 )
  
!---------------------------------------------------------------------72
!---------------------------------------------------------------------72

200 format (6f12.5)
210 format (6e12.5)
300 format (A21)
305 format (A2)
310 format (A6)
320 format (A32)
  
1000 print*
  
  return
end subroutine qnewton1b


