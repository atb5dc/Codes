    program 2D_laplace
    !
!****************************************************************************************
!********************!
! Author:: Ayodeji BodeOke
! Date:: 10 November 2015
! This code solve a 2D laplace equation with initial and boundary condition using point
! and line methods, given user input.
! The accuracy of the schemes are caluculated using error norms, the effect of relaxation
!on each method as
! the number of iterations and rates of convergence for each method are computed
! The different schemes have been written using the finite difference approximations of
! the PDE Ut=Uxx + Uyy.
! where 0<=x<=1, 0<=y<=1, T(0,y)=T(1,y)=0 , T(x,1)=sin(nx), dT/dy(x,0)=0
! Point Methods: Point Jacobi, Point Gauss Seidel
! Line Methods: Line Jacob SoR, Line Gauss Seidel SoR
!
****************************************************************************************
********************!
    implicit none
    INTEGER :: n_row, n_col
    INTEGER :: nx,ny
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: T_old,T_new, err_arr, T_anly,T_star,err_arr_k
    double precision :: L1_norm
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: T_line,T_rhs
    DOUBLE PRECISION :: delta_x, delta_y, length_x , length_y
    INTEGER  :: i,j,n_iter
    INTEGER  :: M,N,tt
    DOUBLE PRECISION :: err , criterion, omega
    double precision, parameter :: pi = 4.0d0*ATAN(1.0d0)
    integer :: choice
    double precision :: a, b ,c, p
    double precision :: spectral_radius, q
    Character (len=1) :: input

    !**********************************************************!
    !! This program solves the 2D heat equation using point and line methods
    length_x = 1
    length_y = 1
    delta_x = 0.025
    delta_y = 0.025
    nx = 1 + length_x/delta_x
    ny = 1 + length_y/delta_y

    !spectral analysis
    q=8

    criterion = 1e-6
    omega = 1.0d0
    allocate(T_anly(nx+1,ny+1))
    allocate(T_new(nx+1,ny+1))
    allocate(T_old(nx+1,ny+1))
    allocate(T_star(nx+1,ny+1))
    allocate(T_rhs(nx+1))
    allocate(T_line(nx+1))
    allocate(err_arr(nx+1,ny+1))
    allocate(err_arr_k(nx+1,ny+1))

    !open file for writing data
    open(101,file='Error_vs_iteration_Pjacobi.dat')
    open(102,file='Error_vs_iteration_PGS.dat')
    open(103,file='Error_vs_iteration_LGS.dat')
    open(104,file='Error_vs_iteration_LJacobi.dat')
    open(111,file='abserror_vs_niter_Pjacobi.dat')
    open(112,file='abserror_vs_niter_PGS.dat')
    open(113,file='abserror_vs_niter__LGS.dat')
    open(114,file='abserror_vs_niter__LJacobi.dat')

    open(11, file='analytical_solution.dat')
    write(11,*) 'TITLE="ANALYTICAL SOLUTION"'
    write(11,*) 'VARIABLES= "X","Y","T"'
    write(11,*) 'ZONE I=',nx+1, 'j=', ny+1,'F=point'

    open(12, file='numerical_solution_Pjacobi.dat')
    write(12,*) 'TITLE="NUMERICAL SOLUTION"'
    write(12,*) 'VARIABLES= "X","Y","T"'
    write(12,*) 'ZONE I=',nx+1, 'j=', ny+1,'F=point'

    open(13, file='absoluteerror_PJacobi.dat')
    write(13,*) 'TITLE="ABSOLUTE ERROR"'
    write(13,*) 'VARIABLES= "X","Y","abs error"'
    write(13,*) 'ZONE I=',nx+1, 'j=', ny+1,'F=point'

    open(14, file='numerical_solution_PGS.dat')
    write(14,*) 'TITLE="NUMERICAL SOLUTION"'
    write(14,*) 'VARIABLES= "X","Y","T"'
    write(14,*) 'ZONE I=',nx+1, 'j=', ny+1,'F=point'

    open(15, file='absoluteerror_PGS.dat')
    write(15,*) 'TITLE="ABSOLUTE ERROR"'
    write(15,*) 'VARIABLES= "X","Y","abs error"'
    write(15,*) 'ZONE I=',nx+1, 'j=', ny+1,'F=point'

    open(16, file='numerical_solution_LGS.dat')
    write(16,*) 'TITLE="NUMERICAL SOLUTION"'
    write(16,*) 'VARIABLES= "X","Y","T"'
    write(16,*) 'ZONE I=',nx+1, 'j=', ny+1,'F=point'

    open(17, file='absoluteerror_LGS.dat')
    write(17,*) 'TITLE="ABSOLUTE ERROR"'
    write(17,*) 'VARIABLES= "X","Y","abs error"'
    write(17,*) 'ZONE I=',nx+1, 'j=', ny+1,'F=point'

    open(18, file='numerical_solution_LJacobi.dat')
    write(18,*) 'TITLE="NUMERICAL SOLUTION"'
    write(18,*) 'VARIABLES= "X","Y","T"'
    write(18,*) 'ZONE I=',nx+1, 'j=', ny+1,'F=point'

    open(19, file='absoluteerror_LJacobi.dat')
    write(19,*) 'TITLE="ABSOLUTE ERROR"'
    write(19,*) 'VARIABLES= "X","Y","abs error"'
    write(19,*) 'ZONE I=',nx+1, 'j=', ny+1,'F=point'

    !! Analytical Solution
    do j=1,ny+1
        do i=1,nx+1
            T_anly(i,j)=((sin(pi*(i-1)*delta_x))*(cosh(pi*(1-(j-1)*delta_y))))/cosh(pi)
            write(11,*) (i-1)*delta_x, (j-1)*delta_y, T_anly(i,j)
        enddo
    enddo

    do while(.true.)
        !!user input! 1 = point jacobi sor , 2 = point GS sor, 3= line GS sor, 4= line Jacobi sor
        write(*,*) "Choose which method to use : "
        write(*,*) "1.PJ-SOR , 2.PGS-SOR , 3.LGS-SOR , 4.LJ-SOR"
        read(*,*) choice
        PRINT*, choice
        if ( choice >= 1 .and. choice <= 4  )then
        !!define initial temp values
        T_new(:,:)  = 10.0d0
        T_old(:,:)  = 10.0d0

        n_iter=0
        err=1.0d0
        err_arr_k(:,:) = 0
        do while(err > criterion)
            err_arr_k = err_arr
            !!*********************************************************************************************!!!POINT METHODS
            !! Point SOR Jacobi
            !! boundary conditions
            !top BC
            do i=1,nx+1
                T_old(i,1)=sin(pi*(i-1)*delta_x)
            enddo
            !left BC
            T_old(1,:)=0.0d0
            !right BC
            T_old(nx+1,:)=0.0d0
            !bottom BC
            T_old(:,ny+1)= T_old(:,ny)
            T_new = T_old

            if (choice==1)then !!use point jacobi
                !!! inner points
                do j=2,ny
                    do i=2,nx
                        T_star(i,j) = (T_old(i-1,j) + T_old(i+1,j) + T_old(i,j-1) + T_old(i,j+1))/4.0d0
                        T_new(i,j)= (1.d0-omega)*T_old(i,j) + omega*T_star(i,j)
                    enddo
                enddo
            elseif(choice==2)then !! use point gauss seidel
                !T_star=T_old
                do j=2,ny
                    do i=2,nx
                        T_star(i,j) = (T_new(i-1,j) + T_old(i+1,j) + T_new(i,j-1) + T_old(i,j+1))/4.0d0
                        T_new(i,j)= (1.d0-omega)*T_old(i,j) + omega*T_star(i,j)
                    enddo
                enddo

                !!*********************************************************************************************!!!LINE METHODS
                !!!!!********!!!!!!!! LINE GAUSS SEIDEL
            elseif (choice==3)then
                a = -0.25d0
                c = a
                b = 1.0d0

                do j=2,ny
                    !!RHS
                    do i=2,nx
                        T_rhs(i) = (1-omega)*T_old(i,j) + omega*(T_new(i,j-1) + T_old(i,j+1))/4.0d0
                    enddo
                    CALL tdma(a,b,c,T_rhs,T_line,2,nx)
                    do i=2,nx
                        T_new(i,j) = T_line(i)
                    enddo
                enddo
                !!!!!********!!!!!!!! LINE JACOBI
            elseif (choice==4)then
                a = -0.25d0
                c = a
                b = 1.0d0

                do j=2,ny
                    !!RHS
                    do i=2,nx
                        T_rhs(i) = (1-omega)*T_old(i,j) + omega*(T_old(i,j-1) + T_old(i,j+1))/4.0d0
                    enddo
                    CALL tdma(a,b,c,T_rhs,T_line,2,nx)
                    do i=2,nx
                        T_new(i,j) = T_line(i)
                    enddo
                enddo
            endif
               

            !!for the dowhile loop
            err_arr = T_new-T_old
            err = sum(abs(err_arr))
            n_iter=n_iter+1
            T_old = T_new

            if (choice==1) write(111,*) n_iter, log(err/(nx*ny))
            if (choice==2) write(112,*) n_iter, log(err/(nx*ny))
            if (choice==3) write(113,*) n_iter, log(err/(nx*ny))
            if (choice==4) write(114,*) n_iter, log(err/(nx*ny))

            if (choice==4) write(104,*) n_iter, err
            if (choice==3) write(103,*) n_iter, err
            if (choice==2) write(102,*) n_iter, err
            if (choice==1) write(101,*) n_iter, err
        enddo
             
        L1_norm = (sum(abs(T_new-T_anly)))/((nx-1)*(ny-1))
        spectral_radius = maxval((err_arr)/(err_arr_k))
        p = q/-log(spectral_radius)
        print*, spectral_radius , n_iter, p , n_iter/p, -log(spectral_radius)

        !write numerical solution and error to file
        do j=1,ny+1
            do i=1,nx+1
                if (choice==1)then
                    write(12,*) (i-1)*delta_x, (j-1)*delta_y, T_new(i,j)
                    write(13,*) -length_x*(i-1)*delta_x, -length_y*(j-1)*delta_y,  abs((T_new(i,j)-T_anly(i,j)))
                endif
                if(choice==2)then
                    write(14,*) (i-1)*delta_x, (j-1)*delta_y, T_new(i,j)
                    write(15,*) -length_x*(i-1)*delta_x, -length_y*(j-1)*delta_y, abs((T_new(i,j)-T_anly(i,j)))
                endif
                if(choice==3)then
                    write(16,*) (i-1)*delta_x, (j-1)*delta_y, T_new(i,j)
                    write(17,*) (i-1)*delta_x, (j-1)*delta_y, abs((T_new(i,j)-T_anly(i,j)))
                endif
                if(choice==4)then
                    write(18,*) (i-1)*delta_x, (j-1)*delta_y, T_new(i,j)
                    write(19,*) (i-1)*delta_x, (j-1)*delta_y, abs((T_new(i,j)-T_anly(i,j)))
                endif
            enddo
        enddo

        write(*,*) "Do you wish to perform another calculation (Y/N)?"
        read(*,*) input
        if (input .eq. "Y") then
            cycle
        elseif (input .eq. "N") then
            stop
        endif
        else
            write(*,*) "check your input"
            endif
    enddo


    !!close all open files
    close(101)
    close(102)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)

    end