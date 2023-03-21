Module var_and_data
    
    implicit none
    
    integer, parameter :: double=selected_real_kind(15,307)
    
    complex(kind=double) :: im=(0,1) !imaginary constant
    
    real(kind=double) :: hbar=1.0
    
    real(kind=double), parameter :: pi=3.141592653589793238462643d0
    
    integer :: Npx, Npy, Npt, Nprint  !number of x-,y-positions, time steps
    
    integer :: Npix_x,Npix_y,printEn_count  !Output wavefunction number of point on x,y and counter to print energy every # step - Fabio June 8 2020
    
    integer :: Ngrid                        !Grid update counter - Fabio Dec. 2020
    
    integer :: flag_potential, flag_incond1, flag_incond2  !flag to choose the kind of potential and initial condition of species 1 and 2 - Fabio June 8 2020 

    real(kind=double) :: Xmax, Ymax, Tmax !Maximum absolute x-,y-values, time 
    
    real(kind=double), allocatable :: dx(:,:), dy(:,:)
    
    real(kind=double) :: dt
    
    real(kind=double) :: sig1, sig2, x01, y01, x02, y02  !stddevs, center coordinates for each Gaussian

    real(kind=double) :: h_leng !Healing length of the vortex
    
    real(kind=double) :: Omega !Angular speed of the trap [rad/s]   - Fabio Nov. 2020
    
    real(kind=double) :: wx1, wy1, wx2, wy2 !omega-constants for potential function

    real(kind=double) :: m1, m2, m2ratio, N1, N2, g1, g2, g12, g21 !masses of particles, numbers of particles, interaction energy constants
    
    real(kind=double), allocatable :: Xp(:,:), Yp(:,:) !x-, y-position vectors
    
    real(kind=double), allocatable :: V(:,:,:) !potential function
    
    complex(kind=double), allocatable :: Psi(:,:,:) !wave function
    
    character(len=1024) :: datafile !name of file containing energy and volumes
    
    !Other thomas method arrays
    
    real(kind=double) :: volume1, volume2 !volume variables
    
    real(kind=double) :: total_energy, energy1, energy2, i_energy1, i_energy2, imaginary_energy !energy variables 
    
    real(kind=double) :: Bnd01,Bnd02,Bnd1,Bnd2
    
    real(kind=double) :: Lz1, Lz2, im_Lz1, im_Lz2       !angular Momentum Lz variables - Fabio Nov. 2020
    
    real(kind=double) :: time       !real time variable
    
end module var_and_data
    
    
program BEC2d_2GPE_AMR_v07
    
    use var_and_data
    implicit none
    
    integer :: t, z, k, i, j                             !time and print counters
    real(kind=double) :: g1_0, g2_0, g12_0      !local reading interaction values  
    real(kind=double),allocatable :: vtest(:)
    real(kind=double) :: xi
    
    write(*,*) 'Reading Input Parameters:...'
    open(unit=342, file = "parameters.dat", status="old", position="asis")
    !reads parameters from file
    READ(342,*) Npx, Npy, Npt                                    !space-time number of points    
    READ(342,*) Npix_x, Npix_y, Nprint, printEn_count, Ngrid     !Resolution outputfile, number of time slices to print and how "often to print energy", how often to update the grid - Fabio Dec. 2020
    READ(342,*) Xmax, Ymax, Tmax
    READ(342,*) flag_potential
    READ(342,*) flag_incond1, flag_incond2          
    !READ(342,*) dx, dy, dt     !Fabio correction - May 26 2020
    READ(342,*) Omega
    READ(342,*) h_leng,sig1, sig2, x01, y01, x02, y02     
    READ(342,*) wx1, wy1, wx2, wy2   
    READ(342,*) m1, m2, N1, N2, g1_0, g2_0, g12_0      
    
    close(342)         
    
    !Check Resolution Input-Consistency
    if (Npix_x>=Npx) then
        Npix_x=Npx
    endif
    if (Npix_y>=Npy) then
        Npix_y=Npy
    endif
    
    allocate(Xp(Npx,2))
    allocate(Yp(Npy,2))
    allocate(dx(Npx-1,2))
    allocate(dy(Npy-1,2))
    allocate(Psi(Npx,Npy,2))
    allocate(V(Npx,Npy,2))
    
    !Conmpute initial Uniform grid spacings 
    do k=1,2
        do i=1,Npx-1
            dx(i,k)=(1.*(2*Xmax)/(Npx-1))    !contructs differential constants
        enddo
        do i=1,Npy-1
            dy(i,k)=(1.*(2*Ymax)/(Npy-1))
        enddo
    enddo
    
    dt=1.*Tmax/Npt
    
    !m2=m1*m2ratio      !Removed as it was confusing, better have m1 and m2, Fabio July 31 2020
    
    !Interaction constant re-definition
    g1=g1_0*N1
    g2=g2_0*N2
    g12=g12_0*N2
    g21=g12_0*N1
    
    write(*,*) 'Initializing:...'
    call Xbuilding(Xp) !initializes position vectors
    
    call Ybuilding(Yp) 
        
    !initializes potential function
    select case (flag_potential)    
        case (1) 
            call HarmonicPotential(V)       !Harmonic potential
        case (2)
            call CilindricalPotential(V)    !Cilindrical potential
        !case (3)    
        case default
            call HarmonicPotential(V) !default is Harmonic potential  
    end select
    
    !Select Initial conditions for Psi_1
    select case (flag_incond1)    
        case (0) 
            call psi_zero(Psi,1)                        !Psi zero everywhere
        case (1) 
            call psi_gaussian_building(Psi,1)           !Gaussian wave-packet
        case (2)
            call psi_vortex_building(Psi,1)             !Gaussian wave-packet + vortex at the center
        case (3) 
            call psi_ThomasVortex_building(Psi,1)       !Thomas Fermi + vortex at the center
        case (4)
            call psi_ThomasFermi_building(Psi,1)        !Thomas Fermi approx. for Harmonic Potential
        case (5)
            call psi_constant_building(Psi,1)           !Constant Psi
        case (6)
            call psi_constantVortex_building(Psi,1)     !Constant Psi + Vortex at the center
        case default
            call psi_gaussian_building(Psi,1)           !Gaussian wave-packet
    end select
    
    !Select Initial conditions for Psi_2
    select case (flag_incond2)
        case (0) 
            call psi_zero(Psi,2)                        !Psi zero everywhere
        case (1) 
            call psi_gaussian_building(Psi,2)           !Gaussian wave-packet
        case (2)
            call psi_vortex_building(Psi,2)             !Gaussian wave-packet + vortex at the center
        case (3) 
            call psi_ThomasVortex_building(Psi,2)       !Thomas Fermi + vortex at the center
        case (4)
            call psi_ThomasFermi_building(Psi,2)        !Thomas Fermi approx. for Harmonic Potential
        case (5)
            call psi_constant_building(Psi,2)           !Constant Psi
        case (6)
            call psi_constantVortex_building(Psi,2)           !Constant Psi + Vortex at the center
        case default
            call psi_gaussian_building(Psi,2)           !Gaussian wave-packet
        end select
        
    call volumecalc(Psi,volume1,volume2) !calculates initial volumes
    if (volume1==0.) then      !avoid normalization when used the test Psi_zero - Fabio, July 2020 
        volume1=1.
    endif
    if (volume2==0.) then
        volume2=1.
    endif
    
    call normalization(Psi,volume1,volume2) !normalizes potential functions
    
    call volumecalc(Psi,volume1,volume2) !calculates to confirm normalization of wave functions
    
    call firstprinting(volume1,volume2) !prints initial conditions to files
    
    !First AMR - Adaptive Mesh Refinement
    call mesh_refinement(Xp,Yp,dx,dy,Psi)
    call fullRes_printing(Psi,1,'Psi1_firstAMR.dat')
    call fullRes_printing(Psi,2,'Psi2_firstAMR.dat')
    pause
    
    write(datafile,"(A)") "Energy_and_Norms.dat" 

    open(unit=17, file = datafile, status="replace", position="append")                             !creates file to write energy and volumes
    open(unit=18, file = "Boundary_contribution_check.dat", status="replace", position="append")    !creates file to write boundary_contribution
    open(unit=19, file = "AngularMomentum_Lz.dat", status="replace", position="append")             !creates file to write Angular momentum Lz - Fabio Nov. 2020
    
!---------TIME EVOLUTION------------------------------------------------------------------------------------------------------------------
    
    time=0.d0                           !Initial time
    
    z=0 
    
    call stateprinting_new(z,Psi)       !New printing subroutine, with controlled resolution
    
    !Initial Energy Computation and printing
    energy1=0.0
    energy2=0.0    
    i_energy1=0.0
    i_energy2=0.0    
    call energycalc(Psi,V,1,energy1,i_energy1) !calculates the energy of each species
    call energycalc(Psi,V,2,energy2,i_energy2)    
    energy1=N1*energy1                                  !Fabio July 2020
    energy2=N2*energy2
    Lz1=0.d0
    Lz2=0.d0
    im_Lz1=0.d0
    im_Lz2=0.d0
    call angular_momentum(Psi, 1, Lz1, im_Lz1)
    call angular_momentum(Psi, 2, Lz2, im_Lz2)
    !Initial values printing
    call boundary_contribution_check(Psi,Bnd01,Bnd02,Bnd1,Bnd2)
    write(17,500) energy1, volume1, energy2, volume2                      !prints energy and volumes to file
    write(18,502) time, Bnd01, Bnd1, volume1, Bnd02, Bnd2, volume2        !prints boundary contribution and volumes
    write(19,500) Lz1, Lz2, im_Lz1, im_Lz2                                !write the angular momentum of the two components - Fabio Nov. 2020
    
    !Beginning of the time-evolution
    write(*,*) 'Computing:...'
    do t=1, Npt
        
        time=time + dt                !Update time
               
        call CrankNicholson(Psi,V) !solves for Psi(t+dt)
        
        !Compute and print Energy and Norm
        if (modulo(t,printEn_count)==0) then
            !Energy Computation
            energy1=0.0
            energy2=0.0
            i_energy1=0.0
            i_energy2=0.0
            call energycalc(Psi,V,1,energy1,i_energy1)          !calculates the energy of each species
            call energycalc(Psi,V,2,energy2,i_energy2)
            total_energy=0
            imaginary_energy=0
            energy1=N1*energy1                                  !Fabio July 2020
            energy2=N2*energy2
            total_energy=energy1+energy2              
            imaginary_energy=i_energy1+i_energy2
            
            !Angular Momentum Lz Computation                   - Fabio Nov. 2020
            Lz1=0.d0
            Lz2=0.d0
            im_Lz1=0.d0
            im_Lz2=0.d0
            call angular_momentum(Psi, 1, Lz1, im_Lz1)
            call angular_momentum(Psi, 2, Lz2, im_Lz2)
            
            
            !Volume Computation
            volume1=0
            volume2=0
            call volumecalc(Psi,volume1,volume2)                !calculates volume at this time step
            
            !Boundary contribution to the norm check
            call boundary_contribution_check(Psi,Bnd01,Bnd02,Bnd1,Bnd2)
            
            !Print on file
            write(17,500) energy1, volume1, energy2, volume2                      !prints energy and volumes to file
            
            write(18,502) time, Bnd01, Bnd1, volume1, Bnd02, Bnd2, volume2        !prints boundary contribution and volumes
            
            write(19,500) Lz1, Lz2, im_Lz1, im_Lz2                                !write the angular momentum of the two components - Fabio Nov. 2020
            
            !Print status on screen/Output-file
            write(*,501) t,Npt,real(100.d0*(t*1.d0)/(Npt*1.d0))
        end if
        
        !Update the grid "Adaptive Mesh Refinement"
        if (modulo(t,Ngrid)==0) then
            call mesh_refinement(Xp,Yp,dx,dy,Psi)
        endif
        
        !Print wave-functions
        if (modulo(t,(Npt/Nprint))==0) then
            z=z+1
            call stateprinting_new(z,Psi)       !New printing subroutine, with controlled resolution   --- !prints out Psi data files for each species
        endif
    enddo !END time-loop
    !Last Printing
    call lastprinting
    close(17)
    close(18)
    write(*,*) 'DONE.'
    

	500 format (2x,ES23.15E3,3x,ES23.15E3,3x,ES23.15E3,3x,ES23.15E3)  
    501 format ('Step: ',I10,'/',I10,' (',F5.1,'% Complete)')    
    502 format (2x,ES23.15E3,3x,ES23.15E3,3x,ES23.15E3,3x,ES23.15E3,3x,ES23.15E3,3x,ES23.15E3,3x,ES23.15E3)  
    
end program BEC2d_2GPE_AMR_v07


    
!------------------------------------------------------------------------------------------------------    
!--------SUBROUTINES (ordered as called)---------------------------------------------------------------      
!------------------------------------------------------------------------------------------------------
    
subroutine Xbuilding(Xgr) !builds x-position vector

    use var_and_data
    implicit none
    
    real(kind=double), intent(inout) :: Xgr(Npx,2)
    integer :: i,k
    
    do k=1,2
        Xgr(1,k) = -Xmax
        do i=2, Npx
            Xgr(i,k) = Xgr(i-1,k) + dx(i-1,k)
        enddo
    enddo 
    
end subroutine 
    
subroutine Ybuilding(Ygr) !builds y-position vector

    use var_and_data
    implicit none    
    
    real(kind=double), intent(inout) :: Ygr(Npy,2)
    integer :: i,k
    
    do k=1,2
        Ygr(1,k) = -Ymax
        do i=2, Npy
            Ygr(i,k) = Ygr(i-1,k) + dy(i-1,k) 
        enddo
    enddo
    
end subroutine     
    
subroutine psibuilding(f) !builds both Psi species

    use var_and_data
    implicit none
    
    complex(kind=double), intent(inout) :: f(Npx,Npy,2)
    integer :: i,j,k
    real(kind=double) :: sig(2)
    real(kind=double) :: x0(2)
    real(kind=double) :: y0(2)
    real(kind=double) :: N(2)
    real(kind=double) :: current_point
    
    sig(1)=sig1
    sig(2)=sig2
    x0(1)=x01
    x0(2)=x02
    y0(1)=y01
    y0(2)=y02
    N(1)=N1
    N(2)=N2
    
    do k=1,2 
        do i=1, Npx
            do j=1, Npy
                
                if(i==1 .OR. i==Npx .OR. j==1 .OR. j==Npy) then
                    f(i,j,k)=0.0
                else               
                    f(i,j,k)=exp(-((Xp(i,k)-x0(k))**2 + (Yp(j,k)-y0(k))**2)/(sig(k)**2)) 
                endif
                
                current_point=f(i,j,k)
                
                if(current_point .LE. 10**(-10)) then
                    f(i,j,k)=0.0
                endif
                
            enddo
        enddo
    enddo
    
    
end subroutine  

subroutine psi_vortex_building(f,sp) !builds both Psi species

    use var_and_data
    implicit none
    
    integer, intent(in) :: sp
    complex(kind=double), intent(inout) :: f(Npx,Npy,2)
    integer :: i,j,k
    real(kind=double) :: sig(2)
    real(kind=double) :: x0(2)
    real(kind=double) :: y0(2)
    real(kind=double) :: N(2)
    real(kind=double) :: current_point
    real(kind=double) :: f1_r,r,theta
    
    sig(1)=sig1
    sig(2)=sig2
    x0(1)=x01
    x0(2)=x02
    y0(1)=y01
    y0(2)=y02
    N(1)=N1
    N(2)=N2
    
    k=sp
    do i=1, Npx
        do j=1, Npy
                
            if(i==1 .OR. i==Npx .OR. j==1 .OR. j==Npy) then
                f(i,j,k)=0.0
            else               
                r=sqrt( (Xp(i,k)-x0(k))**2. + (Yp(j,k)-y0(k))**2. )
                if ( ((Xp(i,k)-x0(k))==0.d0).AND.((Yp(j,k)-y0(k))==0.d0) ) then
                    theta=0.d0                                  !Singularity point of the angle (center)
                else
                    theta=atan2((Yp(j,k)-y0(k)),(Xp(i,k)-x0(k)))    !Compute the right angle in the right quadrant
                    if ( theta<0 ) then
                        theta = 2*pi + theta
                    endif
                endif
                if (r<= (sqrt(2.)*h_leng) ) then    !Here, h_leng is the healing lenght of the vortex (size of the eye)
                    f1_r=0.5*r/h_leng
                else
                    f1_r=sqrt(1 - (h_leng/r)**2. )
                endif   
                f(i,j,k)= f1_r*exp(im*theta)          !Ref. S.Stringari et. al., Phys. Rev. A 97, 063615 (2018)
                
                f(i,j,k)=f(i,j,k)*exp(-((Xp(i,k)-x0(k))**2 + (Yp(j,k)-y0(k))**2)/(sig(k)**2)) 
            endif          
                
        enddo
    enddo
        
end subroutine     

subroutine psi_gaussian_building(f,sp) !builds both Psi species

    use var_and_data
    implicit none
    
    integer, intent(in) :: sp
    complex(kind=double), intent(inout) :: f(Npx,Npy,2)
    integer :: i,j,k
    real(kind=double) :: sig(2)
    real(kind=double) :: x0(2)
    real(kind=double) :: y0(2)
    real(kind=double) :: N(2)
    real(kind=double) :: current_point
    
    sig(1)=sig1
    sig(2)=sig2
    x0(1)=x01
    x0(2)=x02
    y0(1)=y01
    y0(2)=y02
    N(1)=N1
    N(2)=N2
    
    k=sp
    do i=1, Npx
        do j=1, Npy
                
            if(i==1 .OR. i==Npx .OR. j==1 .OR. j==Npy) then
                f(i,j,k)=0.0
            else               
                f(i,j,k)=exp(-((Xp(i,k)-x0(k))**2 + (Yp(j,k)-y0(k))**2)/(sig(k)**2)) 
            endif
                
            current_point=abs(f(i,j,k))**2.d0
                
            if(current_point .LE. 10**(-10)) then
                f(i,j,k)=0.0
            endif
                
        enddo
    enddo

end subroutine         

subroutine psi_ThomasVortex_building(f,sp) !builds both Psi species

    use var_and_data
    implicit none
    
    integer, intent(in) :: sp
    complex(kind=double), intent(inout) :: f(Npx,Npy,2)
    integer :: i,j,k
    real(kind=double) :: sig(2)
    real(kind=double) :: x0(2)
    real(kind=double) :: y0(2)
    real(kind=double) :: N(2)
    real(kind=double) :: wx(2)
    real(kind=double) :: wy(2)
    real(kind=double) :: m(2)   
    real(kind=double) :: current_point,BECij
    real(kind=double) :: f1_r,r,theta,g,echem,argo,Rc,Vijk,g1_0,g2_0
    
    sig(1)=sig1
    sig(2)=sig2
    x0(1)=x01
    x0(2)=x02
    y0(1)=y01
    y0(2)=y02
    N(1)=N1
    N(2)=N2
    wx(1)=wx1
    wx(2)=wx2
    wy(1)=wy1
    wy(2)=wy2
    m(1)=m1
    m(2)=m2
    
  !assumed wx==wy
    if (sp == 1) then
        !echem=( (3*g1*N1*sqrt(m1*(wx1**2. )))/(4*sqrt(2.)) )**(2.d0/3.d0)
        echem=sqrt(g1*m1*(wx1**2.0)/pi)
        g=g1
        Rc=sqrt(2.*echem/(m1*(wx1**2.0)))
    else
        !echem=( (3*g2*N2*sqrt(m2*(wx2**2. )))/(4*sqrt(2.)) )**(2.d0/3.d0)
        echem=sqrt(g2*m2*(wx2**2.0)/pi)
        g=g2
        Rc=sqrt(2.*echem/(m2*(wx2**2.0)))
    endif
    
    print*,'Thomas-Fermi Radius:',Rc
    !pause
    
    k=sp
    do i=1, Npx
        do j=1, Npy

            if(i==1 .OR. i==Npx .OR. j==1 .OR. j==Npy) then
                f(i,j,k)=0.0
            else 
                r=sqrt( (Xp(i,k)-x0(k))**2. + (Yp(j,k)-y0(k))**2. )
                if ( ((Xp(i,k)-x0(k))==0.d0).AND.((Yp(j,k)-y0(k))==0.d0) ) then
                    theta=0.d0                                  !Singularity point of the angle (center)
                else
                    theta=atan2((Yp(j,k)-y0(k)),(Xp(i,k)-x0(k)))    !Compute the right angle in the right quadrant
                    if ( theta<0 ) then
                        theta = 2*pi + theta
                    endif
                endif
                
                !Create Vortex
                if (r<= (sqrt(2.)*h_leng) ) then    !Here, h_leng is the healing lenght of the vortex (size of the eye)
                    f1_r=0.5*r/h_leng
                else
                    f1_r=sqrt(1 - (h_leng/r)**2. )
                endif   
                f(i,j,k)= f1_r*exp(im*5*theta)          !Ref. S.Stringari et. al., Phys. Rev. A 97, 063615 (2018)
                
                Vijk= 0.5*m(k)*(wx(k)**2)*((Xp(i,k) - x0(k))**2) + 0.5*m(k)*(wy(k)**2)*((Yp(j,k) - y0(k))**2)
                argo = echem - Vijk
                if (argo>0) then
                    BECij= sqrt( argo/g )       !Thomas-Fermi 
                    !print*,argo
                    !pause
                else
                    BECij=0.0
                endif
                
                f(i,j,k)=f(i,j,k)*BECij                 !Thomas-Fermi  + vortex (to work properly: h_heal << Rthomas)
                
            endif
            
            !current_point=abs(f(i,j,k))**2.d0
            !if(current_point .LE. 10**(-10)) then
            !    f(i,j,k)=0.0
            !endif    
            !print*,i,j,f(i,j,k)
                    
        enddo
    enddo
    !pause
end subroutine 
    
subroutine psi_ThomasFermi_building(f,sp) !builds both Psi species

    use var_and_data
    implicit none
    
    integer, intent(in) :: sp
    complex(kind=double), intent(inout) :: f(Npx,Npy,2)
    integer :: i,j,k
    real(kind=double) :: sig(2)
    real(kind=double) :: x0(2)
    real(kind=double) :: y0(2)
    real(kind=double) :: N(2)
    real(kind=double) :: wx(2)
    real(kind=double) :: wy(2)
    real(kind=double) :: m(2)    
    real(kind=double) :: current_point,BECij
    real(kind=double) :: f1_r,r,theta,g,echem,argo,Rc,Vijk,g1_0,g2_0
    
    sig(1)=sig1
    sig(2)=sig2
    x0(1)=x01
    x0(2)=x02
    y0(1)=y01
    y0(2)=y02
    N(1)=N1
    N(2)=N2
    wx(1)=wx1
    wx(2)=wx2
    wy(1)=wy1
    wy(2)=wy2
    m(1)=m1
    m(2)=m2
    
    !assumed wx==wy
    if (sp == 1) then
        !echem=( (3*g1*N1*sqrt(m1*(wx1**2. )))/(4*sqrt(2.)) )**(2.d0/3.d0)
        echem=sqrt(g1*m1*(wx1**2.0)/pi)
        g=g1
        Rc=sqrt(2.*echem/(m1*(wx1**2.0)))
    else
        !echem=( (3*g2*N2*sqrt(m2*(wx2**2. )))/(4*sqrt(2.)) )**(2.d0/3.d0)
        echem=sqrt(g2*m2*(wx2**2.0)/pi)
        g=g2
        Rc=sqrt(2.*echem/(m2*(wx2**2.0)))
    endif
    
    print*,'Thomas-Fermi Radius:',Rc
    !pause
    
    k=sp
    do i=1, Npx
        do j=1, Npy

            if(i==1 .OR. i==Npx .OR. j==1 .OR. j==Npy) then
                f(i,j,k)=0.0
            else
                Vijk= 0.5*m(k)*(wx(k)**2)*((Xp(i,k) - x0(k))**2) + 0.5*m(k)*(wy(k)**2)*((Yp(j,k) - y0(k))**2)
                argo = echem - Vijk
                if (argo>0) then
                    BECij= sqrt( argo/g )       !Thomas-Fermi 
                    !print*,argo
                    !pause
                else
                    BECij=0.0
                endif
                
                f(i,j,k)=BECij                  
                
            endif
            
            !current_point=abs(f(i,j,k))**2.d0
            !if(current_point .LE. 10**(-10)) then
            !    f(i,j,k)=0.0
            !endif    
            !print*,i,j,f(i,j,k)
                    
        enddo
    enddo
    !pause
end subroutine 
    
subroutine psi_constant_building(f,sp) !builds both Psi species

    use var_and_data
    implicit none
    
    integer, intent(in) :: sp
    complex(kind=double), intent(inout) :: f(Npx,Npy,2)
    integer :: i,j,k
    real(kind=double) :: sig(2)
    real(kind=double) :: x0(2)
    real(kind=double) :: y0(2)
    real(kind=double) :: N(2)
    real(kind=double) :: Rpot(2)
    real(kind=double) :: current_point,BECij
    real(kind=double) :: f1_r,r,theta,g,echem,argo,sigma,mu
    
    sig(1)=sig1
    sig(2)=sig2
    x0(1)=x01
    x0(2)=x02
    y0(1)=y01
    y0(2)=y02
    N(1)=N1
    N(2)=N2
    
    Rpot(1)=wx1
    Rpot(2)=wx2
    sigma=0.5 !Rpot(sp)/100
    mu=2*sigma
    k=sp
    do i=1, Npx
        do j=1, Npy

            if(i==1 .OR. i==Npx .OR. j==1 .OR. j==Npy) then
                f(i,j,k)=0.0
            else 
                r=sqrt( Xp(i,k)**2. + Yp(j,k)**2. )
                if (flag_potential==2) then         !if cilindrical potential constant only inside the cilinder
                    if ( r<=(Rpot(k) - mu) ) then
                        f(i,j,k)=1.d0
                    else
                        f(i,j,k)= exp( ( -((r - (Rpot(k) - mu) )/sigma )**2.d0)/2.d0 )
                    endif
                else                    !other potentials
                    f(i,j,k)=1.d0
                endif
            endif
            
            current_point=abs(f(i,j,k))**2.d0
            if(current_point .LE. 10**(-10)) then
                f(i,j,k)=0.0
            endif    
            !print*,i,j,f(i,j,k)
                    
        enddo
    enddo
    !pause
end subroutine 

subroutine psi_zero(f,sp) !builds both Psi species

    use var_and_data
    implicit none
    
    integer, intent(in) :: sp
    complex(kind=double), intent(inout) :: f(Npx,Npy,2)
    integer :: i,j,k
    real(kind=double) :: sig(2)
    real(kind=double) :: x0(2)
    real(kind=double) :: y0(2)
    real(kind=double) :: N(2)
    real(kind=double) :: Rpot(2)
    real(kind=double) :: current_point,BECij
    real(kind=double) :: f1_r,r,theta,g,echem,argo,sigma,mu
    
    k=sp
    do i=1, Npx
        do j=1, Npy
            f(i,j,k)=0.0        
        enddo
    enddo
    !pause
end subroutine     
    
subroutine psi_constantVortex_building(f,sp) !builds both Psi species

    use var_and_data
    implicit none
    
    integer, intent(in) :: sp
    complex(kind=double), intent(inout) :: f(Npx,Npy,2)
    integer :: i,j,k
    real(kind=double) :: sig(2)
    real(kind=double) :: x0(2)
    real(kind=double) :: y0(2)
    real(kind=double) :: N(2)
    real(kind=double) :: Rpot(2)
    real(kind=double) :: current_point,BECij
    real(kind=double) :: f1_r,r,theta,g,echem,argo,sigma,mu,f1_Rpot
    
    sig(1)=sig1
    sig(2)=sig2
    x0(1)=x01
    x0(2)=x02
    y0(1)=y01
    y0(2)=y02
    N(1)=N1
    N(2)=N2
    
    Rpot(1)=wx1
    Rpot(2)=wx2
    sigma=0.5 !Rpot(sp)/100
    mu=2*sigma
    k=sp
    do i=1, Npx
        do j=1, Npy

            if(i==1 .OR. i==Npx .OR. j==1 .OR. j==Npy) then
                f(i,j,k)=0.0
            else 
                r=sqrt( Xp(i,k)**2. + Yp(j,k)**2. )
                if ( ((Xp(i,k)-x0(k))==0.d0).AND.((Yp(j,k)-y0(k))==0.d0) ) then
                    theta=0.d0                                  !Singularity point of the angle (center)
                else
                    theta=atan2((Yp(j,k)-y0(k)),(Xp(i,k)-x0(k)))    !Compute the right angle in the right quadrant
                    if ( theta<0 ) then
                        theta = 2*pi + theta
                    endif
                endif
                
                if (flag_potential==2) then         !if cilindrical potential constant only inside the cilinder
                    if ( r<=(Rpot(k) - mu) ) then
                        f(i,j,k)=1.d0
                        !Create Vortex
                        if (r<= (sqrt(2.)*h_leng) ) then    !Here, h_leng is the healing lenght of the vortex (size of the eye)
                            f1_r=0.5*r/h_leng
                        else
                            f1_r=sqrt(1 - (h_leng/r)**2. )
                        endif   
                        f(i,j,k)= f1_r*exp(im*theta)          !Ref. S.Stringari et. al., Phys. Rev. A 97, 063615 (2018)
                    else
                        f1_Rpot=sqrt(1 - (h_leng/(Rpot(k) - mu))**2. )
                        f(i,j,k)= f1_Rpot*exp( ( -((r - (Rpot(k) - mu) )/sigma )**2.d0)/2.d0 )
                    endif
                else                    !other potentials
                    !Create Vortex
                    if (r<= (sqrt(2.)*h_leng) ) then    !Here, h_leng is the healing lenght of the vortex (size of the eye)
                        f1_r=0.5*r/h_leng
                    else
                        f1_r=sqrt(1 - (h_leng/r)**2. )
                    endif   
                    f(i,j,k)= f1_r*exp(im*theta)          !Ref. S.Stringari et. al., Phys. Rev. A 97, 063615 (2018)
                endif
                !Core
            endif
            
            current_point=abs(f(i,j,k))**2.d0
            if(current_point .LE. 10**(-10)) then
                f(i,j,k)=0.0
            endif    
            !print*,i,j,f(i,j,k)
                    
        enddo
    enddo
    !pause
end subroutine 

subroutine rotate_psi(f,sp)                     !Rotate the wavefunction of an angle Omega*dt (used for rotating initial conditions) - Fabio Nov. 2020

    use var_and_data
    implicit none
    
    integer, intent(in) :: sp
    complex(kind=double), intent(inout) :: f(Npx,Npy,2)
    integer :: i,j,k
    complex(kind=double) :: Lz_ij, defx(Npx,Npy), defy(Npx,Npy)
    
    !Compute derivatives along y for all x
    do i=1, Npx
        do j=2, Npy
            defy(i,j)=(f(i,j,sp) - f(i,j-1,sp) )/dy(j,sp)
        enddo
        defy(i,1)=(f(i,2,sp) - f(i,1,sp) )/dy(j,sp)
    enddo
    
    !Compute derivatives along x for all y
    do j=1, Npy
        do i=2, Npx
            defx(i,j)=(f(i,j,sp) - f(i-1,j,sp) )/dx(i,sp)
        enddo
        defx(1,j)=(f(2,j,sp) - f(1,j,sp) )/dx(i,sp)
    enddo
    
    !Perform Rotation
    do i=1, Npx
        do j=1, Npy
            Lz_ij= -(im*hbar)*( Xp(i,sp)*defy(i,j) - Yp(j,sp)*defx(i,j) )
            f(i,j,sp)=f(i,j,sp) - (im/hbar)*Omega*dt*Lz_ij
        enddo
    enddo
    
end subroutine    
    
subroutine HarmonicPotential(f) !builds the potential function

    use var_and_data
    implicit none
    
    real(kind=double), intent(inout) :: f(Npx,Npy,2)
    integer :: i,j,k
    real(kind=double) :: wx(2)
    real(kind=double) :: wy(2)
    real(kind=double) :: m(2)       !Fabio - June 3 2020
    
    wx(1)=wx1
    wx(2)=wx2
    wy(1)=wy1
    wy(2)=wy2
    m(1)=m1
    m(2)=m2
    
    do k=1,2 
        do i=1, Npx
            do j=1, Npy
                
                f(i,j,k)= 0.5*m(k)*(wx(k)**2)*(Xp(i,k)**2) + 0.5*m(k)*(wy(k)**2)*(Yp(j,k)**2)
        
            enddo
        enddo
    enddo
    
end subroutine 

subroutine CilindricalPotential(f) !builds the potential function

    use var_and_data
    implicit none
    
    real(kind=double), intent(inout) :: f(Npx,Npy,2)
    integer :: i,j,k
    real(kind=double) :: wx(2)
    real(kind=double) :: wy(2)
    real(kind=double) :: m(2)       !Fabio - June 3 2020
    real(kind=double) :: r,Rpot
    
    wx(1)=wx1
    wx(2)=wx2
    wy(1)=wy1
    wy(2)=wy2
    m(1)=m1
    m(2)=m2
    
    if (wx1.ne.wx2) then
        write(*,*) 'Radius for cilindrical potential taken as wx1 for both species'
    endif
    
    Rpot=wx1
    
    do k=1,2 
        do i=1, Npx
            do j=1, Npy
                r=sqrt( Xp(i,k)**2. + Yp(j,k)**2. )
                if (r<=Rpot) then
                    f(i,j,k)= 0.d0
                else
                    f(i,j,k)= 10.d15
                endif
            enddo
        enddo
    enddo
    
end subroutine     

subroutine volumecalc(f,v1,v2) !calculates the volume using the 2-D trapezoidal rule

    use var_and_data
    implicit none
    
    complex(kind=double), intent(in) :: f(Npx, Npy, 2)
    real(kind=double), intent(out) :: v1, v2
    integer :: i, j, k
    real(kind=double) :: vol(2)
    v1=0
    v2=0
    vol(1)=0
    vol(2)=0
    do k=1,2
        do i=2, (Npx-2)
            do j=2, (Npy-2)
                vol(k)=vol(k)+0.25*dx(i,k)*dy(j,k)*((abs(f(i,j,k))**2)+(abs(f(i+1,j,k))**2)+(abs(f(i,j+1,k))**2)+(abs(f(i+1,j+1,k))**2))
            enddo
        enddo
    enddo
    
    v1=vol(1)
    v2=vol(2)
    
end subroutine 

subroutine boundary_contribution_check(f,v01,v02,v1,v2) !calculates the contribution of boundaries to the norm

    use var_and_data
    implicit none
    
    complex(kind=double), intent(in) :: f(Npx, Npy, 2)
    real(kind=double), intent(out) :: v01, v02, v1, v2
    integer :: i, j, k
    real(kind=double) :: vol(2),vol0(2)
    v1=0
    v2=0
    vol(1)=0
    vol(2)=0
    vol0(1)=0
    vol0(2)=0
    do k=1,2
        !Internal real Boundary [2: Nxy-1]
        do i=2, (Npx-2)
            vol(k)=vol(k)+0.5*dx(i,k)*( (abs(f(i,2,k))**2) + (abs(f(i+1,2,k))**2) )
            vol(k)=vol(k)+0.5*dx(i,k)*( (abs(f(i,Npy-1,k))**2) + (abs(f(i+1,Npy-1,k))**2) )
        enddo
        do j=2, (Npy-2)
            vol(k)=vol(k)+0.5*dy(j,k)*( (abs(f(2,j,k))**2) + (abs(f(2,j+1,k))**2) )
            vol(k)=vol(k)+0.5*dy(j,k)*( (abs(f(Npx-1,j,k))**2) + (abs(f(Npx-1,j+1,k))**2) )
        enddo
        !Outer "safety"_boundary (should be always zero)
        do i=1, (Npx-1)
            vol0(k)=vol0(k)+0.5*dx(i,k)*( (abs(f(i,1,k))**2) + (abs(f(i+1,1,k))**2) )
            vol0(k)=vol0(k)+0.5*dx(i,k)*( (abs(f(i,Npy,k))**2) + (abs(f(i+1,Npy,k))**2) )
        enddo
        do j=1, (Npy-1)
            vol0(k)=vol0(k)+0.5*dy(j,k)*( (abs(f(1,j,k))**2) + (abs(f(1,j+1,k))**2) )
            vol0(k)=vol0(k)+0.5*dy(j,k)*( (abs(f(Npx,j,k))**2) + (abs(f(Npx,j+1,k))**2) )
        enddo
    enddo
    
    v01=vol0(1)
    v02=vol0(2)
    v1=vol(1)
    v2=vol(2)
    
end subroutine     
    
    
subroutine normalization(f,v1,v2) !normalizes Psi 

    use var_and_data
    implicit none
    
    complex(kind=double), intent(inout) :: f(Npx,Npy,2)
    real(kind=double), intent(in) :: v1, v2
    integer :: i,j,k
    complex(kind=double) :: vol(2)
    
    vol(1)=v1
    vol(2)=v2
    do k=1,2
        do i=2, (Npx-1)
            do j=2, (Npy-1)
                f(i,j,k)=f(i,j,k)/(sqrt(vol(k)))
            enddo
        enddo
    enddo
    
end subroutine 
    
subroutine firstprinting(v1,v2) !prints initial Psi data files and inital volumes
    
    use var_and_data
    implicit none

    real(kind=double), intent(in) :: v1, v2
    integer :: i,j,k
    
    open(unit=6, file="Psi1_initial.dat", status="replace", position="append")
    open(unit=3, file="Psi2_initial.dat", status="replace", position="append")
    
    do i=1, Npx
        do j=1, Npy
            write(6,57) Xp(i,1), Yp(j,1), real(Psi(i,j,1)), imag(Psi(i,j,1))
            write(3,57) Xp(i,2), Yp(j,2), real(Psi(i,j,2)), imag(Psi(i,j,2))
        enddo
    enddo
    
    close(6)
    close(3)
    
    open(unit=4, file="initial_volumes.dat", status="replace", position="append")
    write(4,43) v1, v2
    close(4)
    
43 format (x,ES23.15E3,3x,ES23.15E3)
57 format (x,ES23.15E3,3x,ES23.15E3,3x,ES23.15E3,3x,ES23.15E3)
end subroutine    
    
subroutine lastprinting !prints initial Psi data files and inital volumes
    
    use var_and_data
    implicit none

    integer :: i,j,k
    
    open(unit=6, file="Psi1_final.dat", status="replace", position="append")
    open(unit=3, file="Psi2_final.dat", status="replace", position="append")
    
    do i=1, Npx
        do j=1, Npy
            write(6,57) Xp(i,1), Yp(j,1), real(Psi(i,j,1)), imag(Psi(i,j,1))
            write(3,57) Xp(i,2), Yp(j,2), real(Psi(i,j,2)), imag(Psi(i,j,2))
        enddo
    enddo
    
    close(6)
    close(3)
    
    
43 format (x,ES23.15E3,3x,ES23.15E3)
57 format (x,ES23.15E3,3x,ES23.15E3,3x,ES23.15E3,3x,ES23.15E3)
end subroutine     

subroutine CrankNicholson(f,g) !Time evolution subroutine

    use var_and_data
    implicit none
    
    complex(kind=double), intent(inout) :: f(Npx,Npy,2) !psi
    real(kind=double), intent(in) :: g(Npx,Npy,2) !potential
    complex(kind=double) :: psinext(Npx,Npy,2), psi_half(Npx,Npy,2), psi_avg(Npx,Npy,2) !arrays for next step
    integer :: i,j,k
    
    real(kind=double) :: modded_dt

    
    modded_dt = dt/1.


    call thomas(f,g,psi_half,modded_dt) !solves a tridiagonal system !!only have one thomas method that solves both 
    
    do k=1,2
        do i=1, Npx
            do j=1, Npy
                psi_avg(i,j,k)=0.5*(f(i,j,k)+psi_half(i,j,k))
                if(i==1 .OR. i==Npx .OR. j==1 .OR. j==Npy) then         !Fabio July 31st 2020
                    psi_avg(i,j,k)=0.0
                endif
            enddo
        enddo
    enddo


    call second_thomas(f,g,psi_avg,psinext,modded_dt)


    do k=1,2
        do i=1, Npx
            do j=1, Npy
                if(i==1 .OR. i==Npx .OR. j==1 .OR. j==Npy) then
                    psinext(i,j,k)=0.0
                else
                endif
                f(i,j,k)=psinext(i,j,k)
            enddo
        enddo
    enddo  
    
end subroutine


subroutine energycalc(f,w,sp,energy,im_energy) !calculates energy

use var_and_data
implicit none

real(kind=double) :: mass, PE, IE, g,gg !mass dummy, Potential and Interaction energies
integer :: m,n,i,j,sp2 !sp2 is other species of psi
real(kind=double) :: ans
integer, intent(in) :: sp !species of Psi i.e 1 or 2
real(kind=double), intent(in) :: w(Npx,Npy,2) !potential dummy
complex(kind=double), intent(in) :: f(Npx,Npy,2) !Psi dummy
complex(kind=double) :: d2f(Npx,Npy), fconj(Npx,Npy), KE !Laplacian, conjugate of Psi(sp) and kinetic energy
real(kind=double), intent(out) :: energy, im_energy !energy output

if(sp==1)then !sets values for sp, sp2, mass
    sp2=2
    mass=m1
    g=g1    !Fabio - May 27 2020
    gg=g12  !Fabio - July 16 2020
else
    sp2=1
    mass=m2
    g=g2    !Fabio - May 27 2020
    gg=g21  !Fabio - July 16 2020
endif

do i=1, Npx !creates conjugate Psi(sp) array
    do j=1, Npy
        fconj(i,j)=conjg(f(i,j,sp))
    enddo
enddo

do i=2, (Npx-1) !creates Laplacian array
    do j=2, (Npy-1)
        d2f(i,j)=((f(i+1,j,sp)+f(i-1,j,sp)-2*f(i,j,sp))/(dx(i,sp)**2) + (f(i,j+1,sp)+f(i,j-1,sp)-2*f(i,j,sp))/(dy(j,sp)**2))          !Problem def. of second derivative need to be changed for non-uniform dx and dy - Fabio Dec. 2020
    enddo
enddo

energy=0 
KE=0
PE=0
IE=0 !interaction energy

do i=1, (Npx-1) !uses 2-D trapezoidal rule
    do j=1, (Npy-1)
        do m=0,1
            do n=0,1
                KE=KE + 0.25*dx(i,sp)*dy(j,sp)*((-hbar**2)/(2*mass)*fconj(i+m,j+n)*d2f(i+m,j+n)) !Kinectic energy calculation
                PE=PE + 0.25*dx(i,sp)*dy(j,sp)*(w(i+m,j+n,sp)*abs(f(i+m,j+n,sp))**2) !Potential energy calculation
                IE=IE + 0.25*dx(i,sp)*dy(j,sp)*(0.5*g*abs(f(i+m,j+n,sp))**4 + 0.5*gg*(abs(f(i+m,j+n,sp))**2)*abs(f(i+m,j+n,sp2))**2) !Interaction energy calculation
            enddo
        enddo
    enddo
enddo

energy = Real(KE) + PE + IE !total energy
im_energy= Imag(KE)

!write(*,*) "energy", energy, Real(KE), PE, IE

end subroutine energycalc

subroutine angular_momentum(f, sp, Lz, im_Lz)                     !Compute the expectation value of the angular momentum <L_z> - Fabio Nov. 2020

    use var_and_data
    implicit none
    
    integer, intent(in) :: sp
    complex(kind=double), intent(inout) :: f(Npx,Npy,2)
    real(kind=double), intent(out) :: Lz, im_Lz !Angular Momentum Output
    integer :: i,j,k,m,n
    complex(kind=double) :: Lz_ij, Lz_tmp, defx(Npx,Npy), defy(Npx,Npy), f_Lz(Npx,Npy)
    
    !Compute derivatives along y for all x
    do i=1, Npx
        do j=2, Npy
            defy(i,j)=(f(i,j,sp) - f(i,j-1,sp) )/dy(j-1,sp)
        enddo
        defy(i,1)=(f(i,2,sp) - f(i,1,sp) )/dy(1,sp)
    enddo
    
    !Compute derivatives along x for all y
    do j=1, Npy
        do i=2, Npx
            defx(i,j)=(f(i,j,sp) - f(i-1,j,sp) )/dx(i-1,sp)
        enddo
        defx(1,j)=(f(2,j,sp) - f(1,j,sp) )/dx(1,sp)
    enddo
    
    !Compute the action of the Lz operator on f
    do i=1, Npx
        do j=1, Npy
            Lz_ij= -(im*hbar)*( Xp(i,sp)*defy(i,j) - Yp(j,sp)*defx(i,j) )
            f_Lz(i,j)= conjg(f(i,j,sp))*Lz_ij
        enddo
    enddo
    
    !Compute the expectation value of Lz (integral of f_Lz)
    Lz_tmp=0
    do i=1, (Npx-1) !uses 2-D trapezoidal rule
        do j=1, (Npy-1)
            do m=0,1
                do n=0,1
                    Lz_tmp=Lz_tmp + 0.25*dx(i,sp)*dy(j,sp)*f_Lz(i+m,j+n)
                enddo
            enddo
        enddo
    enddo
    
    Lz=real(Lz_tmp)
    im_Lz=imag(Lz_tmp)
    
    !write(*,*) "Lz ", Lz, im_Lz
    
end subroutine   
    
subroutine stateprinting(a,f) !prints out Psi datafiles
    
    use var_and_data
    implicit none
    
    integer, intent(in) :: a !printing step
    integer :: i,j, unit1, unit2
    complex(kind=double), intent(in) :: f(Npx,Npy,2) !psi
    character(len=1024) :: filename1, filename2
    character(len=1024) :: step
    
    write(step, '(I4.4)') a
    
    write(filename1,"(A,I3)") "Psi1_"//trim(adjustl(step))//".dat" 
    write(filename2,"(A,I3)") "Psi2_"//trim(adjustl(step))//".dat"
    
    unit1=111+a
    unit2=222+a
    
    open(unit=unit1, file=trim(filename1), status="unknown", position="append")
    open(unit=unit2, file=trim(filename2), status="unknown", position="append")
    
    do i=1, Npx
        do j=1, Npy
            
            write(unit1,100) Xp(i,1), Yp(j,1), Real(f(i,j,1)), Imag(f(i,j,1))
            write(unit2,100) Xp(i,2), Yp(j,2), Real(f(i,j,2)), Imag(f(i,j,2))
            
        enddo
    enddo
    
    100 format (2x,ES23.15E3,3x,ES23.15E3,3x,ES23.15E3,3x,ES23.15E3)    
    end subroutine    


subroutine stateprinting_new(a,f) !prints out Psi datafiles
    
    use var_and_data
    implicit none
    
    integer, intent(in) :: a !printing step
    integer :: i,j, unit1, unit2, ii, jj, de_ii, de_jj
    complex(kind=double), intent(in) :: f(Npx,Npy,2) !psi
    character(len=1024) :: filename1, filename2
    character(len=1024) :: step
    
    write(step, '(I4.4)') a
    
    write(filename1,"(A,I3)") "Psi1_"//trim(adjustl(step))//".dat" 
    write(filename2,"(A,I3)") "Psi2_"//trim(adjustl(step))//".dat"
    
    unit1=111+a
    unit2=222+a
    
    open(unit=unit1, file=trim(filename1))
    open(unit=unit2, file=trim(filename2))
    
    ii=1
    de_ii=int(Npx/Npix_x)         !int rounds toward zero, Check this for consistency: it works with Npx = Npix_x: To check: what about even odd? Need +1?
    de_jj=int(Npy/Npix_y)
    do i=1, Npix_x
        jj=1
        do j=1, Npix_y

            write(unit1,100) Xp(ii,1), Yp(jj,1), Real(f(ii,jj,1)), Imag(f(ii,jj,1))
            write(unit2,100) Xp(ii,2), Yp(jj,2), Real(f(ii,jj,2)), Imag(f(ii,jj,2))
            
            jj=jj + de_jj
        enddo
        ii=ii + de_ii
    enddo
    
    close(unit1)
    close(unit2)
    100 format (2x,ES23.15E3,3x,ES23.15E3,3x,ES23.15E3,3x,ES23.15E3)    
end subroutine    

subroutine fullRes_printing(f,sp,fname)  !Test initial Psi with mesh redifined 
    
    use var_and_data
    implicit none
    complex(kind=double), intent(in) :: f(Npx,Npy,2) !psi
    integer, intent(in) :: sp
    character (*), intent(in) :: fname
    
    integer :: i,j,k
    
    open(unit=6, file=trim(fname), status="replace", position="append")
    
    do i=1, Npx
        do j=1, Npy
            write(6,57) Xp(i,sp), Yp(j,sp), real(f(i,j,sp)), imag(f(i,j,sp))        
        enddo
    enddo
    
    close(6)

57 format (x,ES23.15E3,3x,ES23.15E3,3x,ES23.15E3,3x,ES23.15E3)
end subroutine    

subroutine test_1dfun_print(xv,fv,Np,fname)  !Test initial Psi with mesh redifined 
    
    use var_and_data
    implicit none
    integer, intent(in) :: Np
    real(kind=double), intent(in) :: xv(Np)
    real(kind=double), intent(in) :: fv(Np) 
    character (*), intent(in) :: fname
    
    integer :: i,j,k
    
    open(unit=6, file=trim(fname), status="replace", position="append")
    
    do i=1, Np
        write(6,59) xv(i), fv(i)         
    enddo
    
    close(6)

59 format (x,ES23.15E3,3x,ES23.15E3)
end subroutine    
    
    
!----------------------------------------------------------------------------------------------------
!--------META SUBROUTINES----------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------

subroutine SecondDerivY(h,d2yh,a,b,c) !calcuates d2/d2y
    
    use var_and_data
    implicit none

    complex(kind=double), intent(in) :: h(Npx,Npy,2) !psi
    complex(kind=double), intent(out) :: d2yh 
    integer, intent(in) :: a,b,c
    
    d2yh=(h(a,b-1,c)+h(a,b+1,c)-2*h(a,b,c))/(dy(b,c)**2)         !Problem def. of second derivative need to be changed for non-uniform dx and dy - Fabio Dec. 2020
            
    end subroutine
    
subroutine SecondDerivX(h,d2xh,a,b,c) !calculates d2/dx2
    
    use var_and_data
    implicit none
    
    complex(kind=double), intent(in) :: h(Npx,Npy,2) !psi
    complex(kind=double), intent(out) :: d2xh
    integer, intent(in) :: a,b,c
    
    d2xh=(h(a-1,b,c)+h(a+1,b,c)-2.*h(a,b,c))/(dx(a,c)**2)        !Problem def. of second derivative need to be changed for non-uniform dx and dy - Fabio Dec. 2020
                    
    end subroutine
        
    
    
subroutine thomas(f,w,g,timestep) !solves the tridiagonal system for both halfsteps

    use var_and_data
    implicit none

    complex(kind=double), intent(in) :: f(Npx,Npy,2) !psi
    real(kind=double), intent(in) :: w(Npx,Npy,2) !potential
    real(kind=double), intent(in) :: timestep !the modified time step from that inputted
    complex(kind=double) :: u(Npx,Npy,2) !psi after first solution
    complex(kind=double), intent(out) :: g(Npx,Npy,2) !psi after thomas method
    integer :: i,j,k
    integer :: m !number of steps
    integer :: sp_self, sp_opp
    real(kind=double) :: g_self, g_opp
    complex(kind=double) :: alpha,a,b,c,cminus1,d,dminus1,d2xu,d2yf 
    complex(kind=double) :: cprime(Npx,2), dprime(Npx,2)
    real(kind=double) :: mass(2)
    
    m=2
    
    mass(1)=m1
    mass(2)=m2
    
    alpha=(im/hbar)*(timestep/m)*(1./2.)
    
    !-------solves d2x2----------------------------------------------------------------------------------
    
    
    do k=1,2
!    open(unit=6, file="der.dat", status="replace", position="append")
        
       if(k==1) then
            sp_self=1; sp_opp=2
            g_self=g1; g_opp=g12
       else
            sp_self=2; sp_opp=1
            g_self=g2; g_opp=g21        !Changed to g21 - Fabio July 16 2020         
       endif
            
       do j=2,(Npy-1)
          do i=2, (Npx-1)
             b=1.+alpha*((-(hbar**2))*(-2)/(2*mass(k)*dx(i,k)**2)+0.5*(w(i,j,k)+g_self*abs(f(i,j,sp_self))**2 + g_opp*abs(f(i,j,sp_opp))**2))        !Problem def. of second derivative need to be changed for non-uniform dx and dy - Fabio Dec. 2020
             
             call SecondDerivY(f,d2yf,i,j,k)
             
             d=f(i,j,k) - alpha*((-(hbar**2))/(2*mass(k))*d2yf + 0.5*(w(i,j,k) + g_self*abs(f(i,j,sp_self))**2 + g_opp*abs(f(i,j,sp_opp))**2)*f(i,j,k))        !Problem def. of second derivative need to be changed for non-uniform dx and dy - Fabio Dec. 2020
             
		     !1*psi   - alpha*(  hbar^2/2m*d2/dy2*psi       + 1/2(    V     + g_self|psi|^2                  + g_opp|psi|^2)*psi)  
             
             if(i .LT. (Npx-1)) then
                c=(im/hbar)*(timestep/m)*(1./2.)*(-(hbar**2))/(2.*mass(k)*dx(i,k)**2)        !Problem def. of second derivative need to be changed for non-uniform dx and dy - Fabio Dec. 2020
             else
                c=0.0
             endif
             
             if(i==2) then
                a=0.0
                cprime(i,k)=c/b
                dprime(i,k)=d/b
             else
                a=(im/hbar)*(timestep/m)*(1./2.)*(-(hbar**2))/(2.*mass(k)*dx(i,k)**2)        !Problem def. of second derivative need to be changed for non-uniform dx and dy - Fabio Dec. 2020
                cprime(i,k)=c/(b-(a*cminus1))
                dprime(i,k)=(d-(a*dminus1))/(b-(a*cminus1))
             endif
             
             
             cminus1=cprime(i,k)
             dminus1=dprime(i,k)
          enddo
          u((Npx-1),j,k)=dprime((Npx-1),k)
          do i=(Npx-2),2,-1
             u(i,j,k)=dprime(i,k)-cprime(i,k)*u(i+1,j,k)                                          
          enddo
       enddo
    enddo


    
    !-------solves d2y2--------------------------------------------------------------------------------------------

    do k=1,2
        
       if(k==1) then
            sp_self=1; sp_opp=2
            g_self=g1; g_opp=g12
       else
            sp_self=2; sp_opp=1
            g_self=g2; g_opp=g21        !Changed to g21 - Fabio July 16 2020         
       endif
       
       do i=2,(Npx-1)
          do j=2, (Npy-1)
             b=1.+alpha*((-(hbar**2))*(-2)/(2*mass(k)*dy(j,k)**2)+0.5*(w(i,j,k)+g_self*abs(u(i,j,sp_self))**2 + g_opp*abs(u(i,j,sp_opp))**2))        !Problem def. of second derivative need to be changed for non-uniform dx and dy - Fabio Dec. 2020
             
             call SecondDerivX(u,d2xu,i,j,k)
             
             d=u(i,j,k) - alpha*((-(hbar**2))/(2*mass(k))*d2xu + 0.5*(w(i,j,k) + g_self*abs(u(i,j,sp_self))**2 + g_opp*abs(u(i,j,sp_opp))**2)*u(i,j,k))        !Problem def. of second derivative need to be changed for non-uniform dx and dy - Fabio Dec. 2020
             
    		!1*psi   - alpha*(  hbar^2/2m*d2/dx2*psi       + 1/2(    V     + g_self|psi|^2                  + g_opp|psi|^2)*psi)             
             
             if(j .LT. (Npy-1)) then
                c=(im/hbar)*(timestep/m)*0.5*(-(hbar**2))/(2.*mass(k)*dy(j,k)**2)        !Problem def. of second derivative need to be changed for non-uniform dx and dy - Fabio Dec. 2020
             else
                c=0.0
             endif
             
             if(j==2) then
                a=0.0
                cprime(j,k)=c/b
                dprime(j,k)=d/b
             else
                a=(im/hbar)*(timestep/m)*0.5*(-(hbar**2))/(2.*mass(k)*dy(j,k)**2)        !Problem def. of second derivative need to be changed for non-uniform dx and dy - Fabio Dec. 2020
                cprime(j,k)=c/(b-(a*cminus1))
                dprime(j,k)=(d-(a*dminus1))/(b-(a*cminus1))
             endif
             cminus1=cprime(j,k)
             dminus1=dprime(j,k)
          enddo
          
          g(i,(Npy-1),k)=dprime((Npy-1),k)
          do j=(Npy-2),2,-1
             g(i,j,k)=dprime(j,k)-cprime(j,k)*g(i,j+1,k)
          enddo
       enddo
    enddo
     
     !both tridiagonal equations solved
     
    end subroutine
 
subroutine second_thomas(f,w,h,g,timestep) !solves the tridiagonal system for both halfsteps

    use var_and_data
    implicit none

    complex(kind=double), intent(in) :: f(Npx,Npy,2), h(Npx,Npy,2) !psi
    real(kind=double), intent(in) :: w(Npx,Npy,2) !potential
    real(kind=double), intent(in) :: timestep !the modified time step from that inputted
    complex(kind=double) :: u(Npx,Npy,2) !psi after first solution
    complex(kind=double), intent(out) :: g(Npx,Npy,2) !psi after thomas method
    integer :: i,j,k
    integer :: m !number of steps
    integer :: sp_self, sp_opp
    real(kind=double) :: g_self, g_opp
    complex(kind=double) :: alpha,a,b,c,cminus1,d,dminus1,d2xu,d2yf 
    complex(kind=double) :: cprime(Npx,2), dprime(Npx,2)
    real(kind=double) :: mass(2)
    
    m=2
    
    mass(1)=m1
    mass(2)=m2
    
    alpha = (im/hbar)*(timestep/m)*(1./2.)
    
    !-------solves d2x2----------------------------------------------------------------------------------
    
    
    do k=1,2
!    open(unit=6, file="der.dat", status="replace", position="append")
        
       if(k==1) then
            sp_self=1; sp_opp=2
            g_self=g1; g_opp=g12
       else
            sp_self=2; sp_opp=1
            g_self=g2; g_opp=g21        !Changed to g21 - Fabio July 16 2020           
       endif
              
       do j=2,(Npy-1)
          do i=2, (Npx-1)
             b=1.+alpha*((-(hbar**2))*(-2)/(2*mass(k)*dx(i,k)**2)+0.5*(w(i,j,k) +g_self*abs(h(i,j,sp_self))**2 + g_opp*abs(h(i,j,sp_opp))**2))        !Problem def. of second derivative need to be changed for non-uniform dx and dy - Fabio Dec. 2020
             
             call SecondDerivY(f,d2yf,i,j,k)
                          
             d=f(i,j,k) - alpha*((-(hbar**2))/(2*mass(k))*d2yf + 0.5*(w(i,j,k) + g_self*abs(h(i,j,sp_self))**2 + g_opp*abs(h(i,j,sp_opp))**2)*f(i,j,k))        !Problem def. of second derivative need to be changed for non-uniform dx and dy - Fabio Dec. 2020
             
		    !1*psi   - alpha*(  hbar^2/2m*d2/dy2*psi    +    1/2(    V    +  g_self|psi|^2           +        g_opp|psi|^2)*psi)       
	

             if(i .LT. (Npx-1)) then
                c=(im/hbar)*(timestep/m)*(1./2.)*(-(hbar**2))/(2.*mass(k)*dx(i,k)**2)        !Problem def. of second derivative need to be changed for non-uniform dx and dy - Fabio Dec. 2020
             else
                c=0.0
             endif
             
             if(i==2) then
                a=0.0
                cprime(i,k)=c/b
                dprime(i,k)=d/b
             else
                a=(im/hbar)*(timestep/m)*(1./2.)*(-(hbar**2))/(2.*mass(k)*dx(i,k)**2)        !Problem def. of second derivative need to be changed for non-uniform dx and dy - Fabio Dec. 2020
                cprime(i,k)=c/(b-(a*cminus1))
                dprime(i,k)=(d-(a*dminus1))/(b-(a*cminus1))
             endif
             
             
             cminus1=cprime(i,k)
             dminus1=dprime(i,k)
          enddo
          u((Npx-1),j,k)=dprime((Npx-1),k)
          do i=(Npx-2),2,-1
             u(i,j,k)=dprime(i,k)-cprime(i,k)*u(i+1,j,k)                                          
          enddo
       enddo
    enddo



    
    
    !-------solves d2y2--------------------------------------------------------------------------------------------

    do k=1,2
       
       if(k==1) then
            sp_self=1; sp_opp=2
            g_self=g1; g_opp=g12
       else
            sp_self=2; sp_opp=1
            g_self=g2; g_opp=g21        !Changed to g21 - Fabio July 16 2020         
       endif
       
       do i=2,(Npx-1)
          do j=2, (Npy-1)
             b=1. + alpha*((-(hbar**2))*(-2)/(2*mass(k)*dy(j,k)**2)+0.5*(w(i,j,k)+g_self*abs(h(i,j,sp_self))**2 + g_opp*abs(h(i,j,sp_opp))**2))        !Problem def. of second derivative need to be changed for non-uniform dx and dy - Fabio Dec. 2020
             
             call SecondDerivX(u,d2xu,i,j,k)
             
             d=u(i,j,k) - alpha*((-(hbar**2))/(2*mass(k))*d2xu + 0.5*(w(i,j,k) + g_self*abs(h(i,j,sp_self))**2 + g_opp*abs(h(i,j,sp_opp))**2)*u(i,j,k))        !Problem def. of second derivative need to be changed for non-uniform dx and dy - Fabio Dec. 2020
                    
             
             if(j .LT. (Npy-1)) then
                c=(im/hbar)*(timestep/m)*0.5*(-(hbar**2))/(2.*mass(k)*dy(j,k)**2)        !Problem def. of second derivative need to be changed for non-uniform dx and dy - Fabio Dec. 2020
             else
                c=0.0
             endif
             
             if(j==2) then
                a=0.0
                cprime(j,k)=c/b
                dprime(j,k)=d/b
             else
                a=(im/hbar)*(timestep/m)*0.5*(-(hbar**2))/(2.*mass(k)*dy(j,k)**2)        !Problem def. of second derivative need to be changed for non-uniform dx and dy - Fabio Dec. 2020
                cprime(j,k)=c/(b-(a*cminus1))
                dprime(j,k)=(d-(a*dminus1))/(b-(a*cminus1))
             endif
             cminus1=cprime(j,k)
             dminus1=dprime(j,k)
          enddo
          
          g(i,(Npy-1),k)=dprime((Npy-1),k)
          do j=(Npy-2),2,-1
             g(i,j,k)=dprime(j,k)-cprime(j,k)*g(i,j+1,k)
          enddo
       enddo
    enddo
     
     !both tridiagonal equations solved
     
end subroutine

!ADAPTIVE MESH REFINEMENT SUBROUTINES

subroutine evaluate_g(Xg,g,Np,sp,xr,g_value,dxr)
    
    use var_and_data
    implicit none 
    
    integer, intent(in) :: Np,sp 
    real(kind=double), intent(in) ::  Xg(Np,2)
    real(kind=double), intent(in) ::  g(Np)
    real(kind=double), intent(in) :: xr
    real(kind=double), intent(out) :: g_value
    real(kind=double), intent(out) :: dxr
    
    integer :: i,found,saved_index
    real(kind=double) :: mr
    
    i=2
    found=0
    saved_index=0
    do while( (found==0).AND.(i<=Np) )
        if ( (Xg(i-1,sp) <= xr).AND.(xr <= Xg(i,sp)) ) then
            found=1
            saved_index=i
        endif
        i=i+1
    enddo
    
    !Return the value (-1 error)
    if (found==0) then
        g_value=-1
        dxr=1
    else
        mr=(g(saved_index) - g(saved_index-1))/(Xg(saved_index,sp) - Xg(saved_index-1,sp))
        g_value= mr*( xr - Xg(saved_index-1,sp) ) + g(saved_index-1) 
        dxr=abs(Xg(saved_index,sp) - Xg(saved_index-1,sp))
    endif
    
end subroutine   

subroutine fit_grid(Xg,Yg,Ngx,Ngy,sp,xr,yr,ix,jy,errflg)    
    
    use var_and_data
    implicit none 
    
    integer, intent(in) :: Ngx,Ngy,sp 
    real(kind=double), intent(in) ::  Xg(Ngx,2)
    real(kind=double), intent(in) ::  Yg(Ngy,2)
    real(kind=double), intent(in) :: xr
    real(kind=double), intent(in) :: yr
    integer, intent(out) :: ix
    integer, intent(out) :: jy
    integer, intent(out) :: errflg
    
    integer :: i,found_x,found_y,saved_index
    
    !X search
    i=2
    found_x=0
    saved_index=0
    do while( (found_x==0).AND.(i<=Ngx) )
        if ( (Xg(i-1,sp) <= xr).AND.(xr <= Xg(i,sp)) ) then
            found_x=1
            saved_index=i
        endif
        i=i+1
    enddo
    ix=saved_index
    
    !Y search
    i=2
    found_y=0
    saved_index=0
    do while( (found_y==0).AND.(i<=Ngy) )
        if ( (Yg(i-1,sp) <= yr).AND.(yr <= Yg(i,sp)) ) then
            found_y=1
            saved_index=i
        endif
        i=i+1
    enddo
    jy=saved_index
    
    !Return the value (-1 error)
    if (found_x*found_y==0) then
        errflg=-1
    else
        errflg=0
    endif
    
end subroutine

subroutine insert_ordered(vec,Np,count,xi)

    use var_and_data
    implicit none 
    
    integer, intent(in) :: Np,count
    real(kind=double), intent(inout) ::  vec(Np)
    real(kind=double), intent(in) ::  xi
    
    integer :: i, found, saved_index
    real(kind=double) :: prev,dummy
    
    if (count+1<=Np) then
        i=2
        found=0
        saved_index=0
        if ( xi <= vec(1) ) then    !has to be first
            found=1
            saved_index=1
        else
            do while( (found==0).AND.(i<=count) )
                if ( (vec(i-1) <= xi).AND.(xi <= vec(i)) ) then
                    found=1
                    saved_index=i
                endif
                i=i+1
            enddo
        endif
    
        if (found==0) then  !add xi to the end
            vec(count+1)=xi
        else                !add xi to the saved position and rescale the entire vector
            prev=vec(saved_index)
            vec(saved_index)=xi
            do i=saved_index+1,count
                dummy=vec(i)
                vec(i)=prev
                prev=dummy
            enddo
            vec(count+1)=prev
        endif
    else
        print*,'Insert_ordered ERROR. Vector already full!'
    endif
    
end subroutine
    
subroutine mesh_refinement(Xgr,Ygr,d_x,d_y,f) !solves the tridiagonal system for both halfsteps

    use var_and_data
    implicit none 
    
    complex(kind=double), intent(inout) :: f(Npx,Npy,2) !psi
    real(kind=double), intent(inout) ::  Xgr(Npx,2),Ygr(Npy,2)
    real(kind=double), intent(inout) ::  d_x(Npx-1,2),d_y(Npy-1,2)
    !real(kind=double), intent(out)   ::  Vp(Npx,Npy,2)
    
    real(kind=double) :: gx(Npx),gy(Npy),listX(Npx),listY(Npy),dxuni,dyuni,r0,xr,yr,pacc,p_i,vol,gval,dxr,dyr,xri,yri,ft,fu
    integer :: i,j,k,count,error,int_error,ix,jy
    complex(kind=double) :: newf(Npx,Npy),y1,y2,y3,y4
    !Test
    real(kind=double) :: xaxis(Npx),yaxis(Npy),g_ij
    
    !Loop on the 2 species to compute the new axis Xgr and Ygr
    do k=1,2
        ! Compute gx Probability density along x
        do i=1, Npx
            gx(i)=0
            do j=1, Npy
                if (i<npx) then     !derivative of psi
                    gx(i)=gx(i) +  abs((abs(f(i+1,j,k))**2.d0 - abs(f(i,j,k))**2.d0)/d_x(i,k))
                else    !i=npx
                    gx(i)=gx(i) +  abs((abs(f(i,j,k))**2.d0 - abs(f(i-1,j,k))**2.d0)/d_x(i-1,k))
                endif
                gx(i)=gx(i) + abs(f(i,j,k))**2 !Psi value 
            enddo
            gx(i)=dsqrt(gx(i)/Npy)
            if (gx(i)<1d-9) then
                gx(i)=0.d0    
            else
                gx(i)=1.d0
            endif
        enddo
        vol = 0
        do i=1,Npx!-1      !Normalize the prob. function (assumed continuous)
            vol = vol + gx(i)!0.5d0*(gx(i) + gx(i+1))*d_x(i,k) 
        enddo
        if (vol.ne.0) gx=gx/vol
        
        ! Compute gy Probability density along y
        do j=1, Npy
            gy(j)=0
            do i=1, Npx
                if (j<npy) then     !derivative of psi
                    gy(j)=gy(j) + abs((abs(f(i,j+1,k))**2.d0 - abs(f(i,j,k))**2.d0)/d_y(j,k))
                else    !j=npy
                    gy(j)=gy(j) + abs((abs(f(i,j,k))**2.d0 - abs(f(i,j-1,k))**2.d0)/d_y(j-1,k))
                endif
                gy(j)=gy(j) + abs(f(i,j,k))**2  !Psi value
            enddo
            gy(j)=dsqrt(gy(j)/Npx)
            if (gy(j)<1d-9) then
                gy(j)=0.d0    
            else
                gy(j)=1.d0
            endif
        enddo
        vol = 0
        do i=1,Npy!-1      !Normalize the prob. function (assumed discrete)
            vol = vol + gy(i)!0.5d0*(gy(i) + gy(i+1))*d_y(i,k)
        enddo
        if (vol.ne.0) gy=gy/vol
        
        do i=1,Npx
            xaxis(i)=Xgr(i,k)
        enddo
        call test_1dfun_print(xaxis,gx,Npx,'plot_gx.dat')
        do i=1,Npy
            yaxis(i)=Ygr(i,k)
        enddo
        call test_1dfun_print(yaxis,gy,Npy,'plot_gy.dat')
        
        ! Sample Npx points along x with prob. gx and compute Xgr of species k
        count=0
        !Assign 10 backbone equally spaced points
        listX(1) = -Xmax
        dxuni=2*Xmax/9
        do i=2, 10
            listX(i) = listX(i-1) + dxuni
        enddo 
        count=10
        do while (count<Npx) 
            call random_number(r0)
            xr=2.d0*Xmax*r0 - Xmax
            call evaluate_g(Xgr,gx,Npx,k,xr,gval,dxr)
            pacc=gval!*dxr               !Continuos probability in the interval where xr is: gval = (g(i-1) + g(i))/2  where x(i-1)<xr<x(i)
            if (pacc<0) print*,'boom'
            call random_number(p_i)
            if ( (p_i<=pacc).AND.(pacc>=0) ) then
                call insert_ordered(listX,Npx,count,xr)
                count=count+1
                !listX(count)=xr
            endif
             
        enddo
        
        ! Sample Npy points along y with prob. gy and compute Ygr of species k
        count=0
        !Assign 10 backbone equally spaced points
        listY(1) = -Ymax
        dyuni=2*Ymax/9
        do i=2, 10
            listY(i) = listY(i-1) + dyuni
        enddo 
        count=10
        do while (count<Npy) 
            call random_number(r0)
            yr=2.d0*Ymax*r0 - Ymax
            call evaluate_g(Ygr,gy,Npy,k,yr,gval,dyr)
            pacc=gval!*dyr               !Continuos probability in the interval where xr is: gval = (g(i-1) + g(i))/2  where y(i-1)<yr<y(i)
            if (pacc<0) print*,'boom'
            call random_number(p_i)
            if ( (p_i<=pacc).AND.(pacc>=0) ) then
                call insert_ordered(listY,Npy,count,yr)
                count=count+1
                !listY(count)=yr
            endif
             
        enddo
        
        !Interpolate the Wavefunction on the new mesh
        int_error=0
        do i=1,Npx
            do j=1,Npy
                call fit_grid(Xgr,Ygr,Npx,Npy,k,listX(i),listY(j),ix,jy,error)
                if (error==0) then
                    !newf(i,j)=(f(ix,jy,k) + f(ix-1,jy,k) + f(ix,jy-1,k) + f(ix-1,jy-1,k) )/4
                    
                    xri=listX(i)
                    yri=listY(j)
                    !Bilinear interpolation
                    y1=f(ix-1,jy-1,k)
                    y2=f(ix,jy-1,k)
                    y3=f(ix,jy,k) 
                    y4=f(ix-1,jy,k)
                    ft=(xri - Xgr(ix-1,k) )/(Xgr(ix,k) - Xgr(ix-1,k))
                    fu=(yri - Ygr(jy-1,k) )/(Ygr(jy,k) - Ygr(jy-1,k))
                    newf(i,j)=(1-ft)*(1-fu)*y1 + ft*(1-fu)*y2 + ft*fu*y3 + (1-ft)*fu*y4
                else
                    print*, "AMR error interpolating the wavefunction, point not on the grid."
                    int_error=1
                endif
            enddo
        enddo
        do i=1,Npx
            do j=1,Npy
                f(i,j,k)=newf(i,j)
            enddo
        enddo
        
        
        ! Save new Xgr and Ygr for species 1 and 2
        do i=1,Npx
            Xgr(i,k)= listX(i)   
        enddo
        do i=1,Npy
            Ygr(i,k)= listY(i)   
        enddo
        
        ! Compute d_x and d_y
        do i=1,Npx-1
            d_x(i,k)= Xgr(i+1,k) - Xgr(i,k)    
        enddo
        do i=1,Npy-1
            d_y(i,k)= Ygr(i+1,k) - Ygr(i,k)    
        enddo
        
    enddo
    
    !Compute mapping between axis of species 1 and 2
    
    
end subroutine    

    

    
