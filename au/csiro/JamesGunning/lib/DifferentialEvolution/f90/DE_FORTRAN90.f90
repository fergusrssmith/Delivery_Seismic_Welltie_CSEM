module data_type
  implicit none 
  integer(kind=4), parameter :: IB=4, RP=8
end module data_type

module data_Rosen
  use data_type
  implicit none
  integer(kind=IB), parameter ::  Dim_XC=10
end module data_Rosen

module data_HDE
	use data_type
	use data_Rosen
	implicit none
	integer(kind=IB), parameter :: NP=50, itermax=20000, strategy=2, &
                                 refresh=10000, iwrite=7
    real(kind=RP), parameter :: VTR=-1.0e-4_RP, CR_XC=0.5_RP
    real(kind=RP) :: F_XC=0.8_RP   
	real(kind=RP), dimension(Dim_XC), parameter :: XCmin=-10.0_RP, XCmax=10.0_RP
	real(kind=RP), dimension(Dim_XC) :: bestmem_XC
    integer(kind=IB) :: nfeval
    real(kind=RP) :: bestval
end module data_HDE


program Rosen
  use data_type
  use data_Rosen
  use data_HDE
  implicit none
  integer(kind=IB) :: i
  integer (kind=IB), dimension(8) :: time 
  intrinsic date_and_time
  external FTN
  open(iwrite,file='Rosen.dat')
  call date_and_time(values=time)
  write(unit=iwrite, FMT=11) time(1:3), time(5:7)

  call DE_Fortran90(FTN, Dim_XC, XCmin, XCmax, VTR, NP, itermax, F_XC,&
                CR_XC, strategy, refresh, iwrite, bestmem_XC, &
				bestval, nfeval) 
              
              write(iwrite,205) NP, nfeval
			  write(iwrite,FMT=201) F_XC, CR_XC
 		      write(iwrite,FMT=200) bestval
			  do i=1,Dim_XC
		         write(iwrite,FMT=202) i,bestmem_XC(i)
              end do
200	    format(/2x, 'Bestval=', ES14.7)
201     format(2x, 'F_XC=',F6.3, 2x, 'CR_XC=', F5.3)
202     format(2x, 'best_XC(',I3,')=',es14.7)
205     format(2x, 'NP=', I4, 2x, 'No. function call=', I9)
 	call date_and_time(values=time)
	write(unit=iwrite, FMT=10)time(1:3), time(5:7)
10	format(/1x, 'End of Program. Date:', I4, '/', I2,'/', I2, ', Time: ', I2,':',I2,':',I2) 
11	format(1x, 'Beginning of Program. Date:', I4, '/', I2,'/', I2, ', Time: ', I2,':',I2,':',I2)
end program Rosen


subroutine FTN(XC, objval)
  use data_type
  use data_Rosen
  implicit none
  real(kind=RP), dimension(Dim_XC), intent(in) :: XC
  real(kind=RP), intent(out) :: objval
  real(kind=RP), dimension(Dim_XC) :: x
  integer(kind=IB) :: i
  i=Dim_XC 
  x=XC
  objval=sum(100.0*(x(1:i-1)**2-x(2:i))**2+(1.0-x(1:i-1))**2)
end subroutine FTN



subroutine DE_Fortran90(obj, Dim_XC, XCmin, XCmax, VTR, NP, itermax, F_XC, &
                      CR_XC, strategy, refresh, iwrite, bestmem_XC, bestval, nfeval)
!.......................................................................
!    
! Differential Evolution for Optimal Control Problems
!
!.......................................................................
!  This Fortran 90 program translates from the original MATLAB 
!  version of differential evolution (DE). This FORTRAN 90 code 
!  has been tested on Compaq Visual Fortran v6.1. 
!  Any users new to the DE are encouraged to read the article of Storn and Price. 
!
!  Refences:
!  Storn, R., and Price, K.V., (1996). Minimizing the real function of the 
!    ICEC'96 contest by differential evolution. IEEE conf. on Evolutionary 
!    Comutation, 842-844.
!
!  This Fortran 90 program is written by Dr. Feng-Sheng Wang 
!  Department of Chemical Engineering, National Chung Cheng University, 
!  Chia-Yi 621, Taiwan, e-mail: chmfsw@ccunix.ccu.edu.tw
!.........................................................................
!                obj : The user provided file for evlauting the objective function.
!                      subroutine obj(xc,fitness)
!                      where "xc" is the real decision parameter vector.(input)
!                            "fitness" is the fitness value.(output)
!             Dim_XC : Dimension of the real decision parameters.
!      XCmin(Dim_XC) : The lower bound of the real decision parameters.
!      XCmax(Dim_XC) : The upper bound of the real decision parameters.
!                VTR : The expected fitness value to reach.
!                 NP : Population size.
!            itermax : The maximum number of iteration.
!               F_XC : Mutation scaling factor for real decision parameters.
!              CR_XC : Crossover factor for real decision parameters.
!           strategy : The strategy of the mutation operations is used in HDE.
!            refresh : The intermediate output will be produced after "refresh"
!                      iterations. No intermediate output will be produced if
!                      "refresh < 1".
!             iwrite : The unit specfier for writing to an external data file.
! bestmen_XC(Dim_XC) : The best real decision parameters.
!              bestval : The best objective function.
!             nfeval : The number of function call.

	 use data_type, only : IB, RP
	 implicit none
	 integer(kind=IB), intent(in) :: NP, Dim_XC, itermax, strategy,   &
	                                 iwrite, refresh
     real(kind=RP), intent(in) :: VTR, CR_XC, F_XC
	 real(kind=RP), dimension(Dim_XC), intent(in) :: XCmin, XCmax
     real(kind=RP), dimension(Dim_XC), intent(inout) :: bestmem_XC
	 real(kind=RP), intent(out) :: bestval
	 integer(kind=IB), intent(out) :: nfeval     
     real(kind=RP), dimension(NP,Dim_XC) :: pop_XC, bm_XC, mui_XC, mpo_XC,   &
	                                        popold_XC, rand_XC, ui_XC
     integer(kind=IB) :: i, ibest, ppp, iter, k
     integer(kind=IB), dimension(1) :: inumber
     integer(kind=IB), dimension(NP) :: rot, a1, a2, a3, a4, a5, rt
     integer(kind=IB), dimension(4) :: ind
     real(kind=RP) :: tempval, s, diversity, xxval
	 real(kind=RP), dimension(NP) :: val
     real(kind=RP), dimension(Dim_XC) :: bestmemit_XC, temp_C, xx_XC
     real(kind=RP), dimension(Dim_XC) :: rand_C1, rand_C2
	 external  obj
	 intrinsic max, min, random_number, mod, abs, any, all, maxloc
	 interface
        function randperm(num)
		   use data_type, only : IB	
		   implicit none
	       integer(kind=IB), intent(in) :: num
           integer(kind=IB), dimension(num) :: randperm
        end function randperm
     end interface
 !!-----Initialize a population --------------------------------------------!!

        pop_XC=0.0_RP
        do i=1,NP
           call random_number(rand_C1)
           pop_XC(i,:)=XCmin+rand_C1*(XCmax-XCmin)
		end do
 
!!--------------------------------------------------------------------------!!

!!------Evaluate fitness functions and find the best member-----------------!!
     val=0.0_RP
     nfeval=0
     ibest=1
     call obj(pop_XC(1,:), val(1))
     bestval=val(1)
     nfeval=nfeval+1
     do i=2,NP
        call obj(pop_XC(i,:), val(i))
        nfeval=nfeval+1
        if (val(i) < bestval) then
           ibest=i
	       bestval=val(i)
        end if
     end do  	 
     bestmemit_XC=pop_XC(ibest,:)
     bestmem_XC=bestmemit_XC
!!--------------------------------------------------------------------------!!

     bm_XC=0.0_RP
     rot=(/(i,i=0,NP-1)/)
     ppp=0
     iter=1 
!!------Perform evolutionary computation------------------------------------!! 

     do while (iter <= itermax)
        popold_XC=pop_XC

!!------Mutation operation--------------------------------------------------!!
        ind=randperm(4)
        a1=randperm(NP)
        rt=mod(rot+ind(1),NP)
        a2=a1(rt+1)
        rt=mod(rot+ind(2),NP)
        a3=a2(rt+1)
        rt=mod(rot+ind(3),NP)
        a4=a3(rt+1)
        rt=mod(rot+ind(4),NP)
        a5=a4(rt+1)
        bm_XC=spread(bestmemit_XC, DIM=1, NCOPIES=NP)
!---- select a mutation strategy-----------------------------------------!
		select case (strategy)
        case (1)
           ui_XC=bm_XC+F_XC*(popold_XC(a1,:)-popold_XC(a2,:))

	    case default
           ui_XC=popold_XC(a3,:)+F_XC*(popold_XC(a1,:)-popold_XC(a2,:))

	    case (3)
           ui_XC=popold_XC+F_XC*(bm_XC-popold_XC+popold_XC(a1,:)-popold_XC(a2,:))

	    case (4)
           ui_XC=bm_XC+F_XC*(popold_XC(a1,:)-popold_XC(a2,:)+popold_XC(a3,:)-popold_XC(a4,:))

		case (5)
		   ui_XC=popold_XC(a5,:)+F_XC*(popold_XC(a1,:)-popold_XC(a2,:)+popold_XC(a3,:) &
		         -popold_XC(a4,:))

        end select
!!--------------------------------------------------------------------------!!
!!------Crossover operation-------------------------------------------------!!
        call random_number(rand_XC)
           mui_XC=0.0_RP
           mpo_XC=0.0_RP
        where (rand_XC < CR_XC)
           mui_XC=1.0_RP
           mpo_XC=0.0_RP
        elsewhere
           mui_XC=0.0_RP
           mpo_XC=1.0_RP
        end where
		ui_XC=popold_XC*mpo_XC+ui_XC*mui_XC
!!--------------------------------------------------------------------------!!
!!------Evaluate fitness functions and find the best member-----------------!!
        do i=1,NP
!!------Confine each of feasible individuals in the lower-upper bound-------!!
		   ui_XC(i,:)=max(min(ui_XC(i,:),XCmax),XCmin)
           call obj(ui_XC(i,:), tempval)
	       nfeval=nfeval+1
	       if (tempval < val(i)) then
	          pop_XC(i,:)=ui_XC(i,:)
	          val(i)=tempval
	          if (tempval < bestval) then
	             bestval=tempval
		         bestmem_XC=ui_XC(i,:)
              end if
           end if
        end do

	    if( (refresh > 0) .and. (mod(iter,refresh)==0)) then
   		     write(unit=iwrite,FMT=203) iter
		     write(unit=*, FMT=203) iter	
 		     do i=1,Dim_XC
		         write(unit=iwrite, FMT=202) i, bestmem_XC(i)
				 write(*,FMT=202) i,bestmem_XC(i)
             end do
			 write(unit=iwrite, FMT=201) bestval
			 write(unit=*, FMT=201) bestval	
        end if
        iter=iter+1
        if ( bestval <= VTR .and. refresh > 0) then
           write(unit=iwrite, FMT=*) ' The best fitness is smaller than VTR'
		   write(unit=*, FMT=*) 'The best fitness is smaller than VTR' 
           exit
        endif
	 end do
!!------end the evolutionary computation------------------------------!!
201 format(2x, 'bestval=', ES14.7, /)
202 format(5x, 'bestmem_XC(', I3, ')=', ES12.5)
203 format(2x, 'No. of iteration =', I8)
end subroutine DE_Fortran90


function randperm(num)
  use data_type, only : IB, RP
  implicit none
  integer(kind=IB), intent(in) :: num
  integer(kind=IB) :: number, i, j, k
  integer(kind=IB), dimension(num) :: randperm
  real(kind=RP), dimension(num) :: rand2
  intrinsic random_number
  call random_number(rand2)
  do i=1,num
     number=1
     do j=1,num
        if (rand2(i) > rand2(j)) then
	       number=number+1
        end if
     end do
     do k=1,i-1
        if (rand2(i) == rand2(k)) then
	       number=number+1
        end if
     end do
     randperm(i)=number
  end do
  return
end function randperm