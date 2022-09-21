!potential is DLVO + Hydrophobic type
program BD
use omp_lib
implicit none
integer,parameter::natom=6400 , max_bin=120 !natom = no of particles present
real,parameter:: L=18.00 , Lx = 39.5  !Box length
integer :: nstep,  nequ  , nequ_step, nsteady , iseed, tconfig
common /BLOCK2/ nstep , nequ , nequ_step , nsteady , tconfig
real :: alpha,alph,alpc
common /BLOCK3/ alpha, alph, alpc 
integer i,j,k,count,t,run_start
real,dimension(natom)::Rx,Ry,Rz,fx,fy,fz, Rx_npc, Ry_npc, Rz_npc
common / BLOCK1 / Rx, Ry, Rz, Rx_npc, Ry_npc, Rz_npc
real ,parameter:: gama = 569 , gama2 = 106, T1 = 182 , T2 = 248.67 ,  sigma = 1.0, D1 = T1/gama , D2 = T2/gama2   
! D = diffusion Coefficient
real dens,epot , ch_sys_h, ch_dil_h,ch_sys_c, ch_dil_c, pack_fr, pack_fr_cross_h,pack_fr_cross_c
real,parameter:: pi=3.1415927 , kappa_h= 3.0 , kappa_c=2.5 , surface_pot_h = 1.2, surface_pot_c=1,radius_scale=100.0 
real :: time_step = 0.0001 
!time_Step = integration time step
real,dimension(max_bin)::bin1,bin2,bine1,bine2,bin_eq,bine_eq !change
real step_eq , step  !increment of bin
real dens_hot , dens_cold , part_hot , part_cold !particle number on hotter and colder side
double precision wtime
!charge = sqrt(beta_epsilon*T2*(1+kappa/2.0)*(1+kappa/2.0))
!--------------------
nstep=2000000 !total integration step
nequ = 50000  !Equilibrium time in Equilibrium Calculation
nequ_step=100000 !Total Equilibrium Calculation
nsteady = 1900000 !steady state time 
!-------------------------------------
alpha = 50.0
alph  = 50.0
alpc  = 300.0
!----------------------------------------

run_start = 0
step =L/(3*max_bin)
step_eq = L/(2*max_bin)
dens =natom/(Lx*L*L)
ch_dil_h = radius_scale*(1+kappa_h/2.0)*surface_pot_h
ch_dil_c = radius_scale*(1+kappa_c/2.0)*surface_pot_c
pack_fr = dens*(pi/6.0)
pack_fr_cross_h = (kappa_h*kappa_h)/(3*(1+kappa_h/2.0))/4.0
pack_fr_cross_c = (kappa_c*kappa_c)/(3*(1+kappa_c/2.0))/4.0

ch_sys_h = ch_dil_h/(1+pack_fr/pack_fr_cross_h)
ch_sys_c = ch_dil_c/(1+pack_fr/pack_fr_cross_c)

dens_hot = 0.0
dens_cold = 0.0
tconfig = 500  !skip of time to calculate g(r) , save coordinates
!iseed = 16 
iseed = 106701 !seed for random number generator
call ranint(iseed)
write(*,*)dens , kappa_h, kappa_c,pack_fr 

write ( *, '(a,i8)' ) '  The number of processors available = ', omp_get_num_procs ( )
write ( *, '(a,i8)' ) '  The number of threads available    = ', omp_get_max_threads ( )

open(1,file = 'fcc-lattice6400_disorder1.dat') !fcc lattice file
open(11,file = 'forcefact_test10.dat') !force calculation file
open(12,file = 'epot1l_eq_1.dat') !potential energy in equilibrium calculation per step file 
open(13,file = 'epot1l_soret_1.dat') !potential energy in soret effect calculation per step file 
open(2,file="gr_20_4000_1.dat") !pair_correlation function file after soret effect
open(3,file="particle_number4000_1.dat") !particle number on both sides
open(4,file="flux_x1.dat") !flux in x direction
open(7,file = 'final-coordinate1.dat')
open(5,file="Information1.dat") !flux in x direction
open(8,file='Equilibrium Coordinates1.dat')
open(9,file='Steadystate-Coordinates1.dat')
write(5,*) 'number of particle :' , natom ,'Length of box :' , L , 'density of initial homogeneous system ' , dens
write(5,*) 'constant surface potential hot' , surface_pot_h , 'constant surface potential cold' , surface_pot_c
write(5,*) 'packing fraction :' , pack_fr , 'crossover packinf fraction :' , pack_fr_cross_h
write(5,*) 'radius of particle in simulation in bjerrum length scale' , radius_scale
write(5,*) 'charge in dilution-hot: ' ,ch_dil_h , 'charge in the system-hot: ', ch_sys_h
write(5,*) 'charge in dilution-cold: ' ,ch_dil_c , 'charge in the system-cold: ', ch_sys_c
write(5,*) 'inverese debye screen length-hot' ,kappa_h , 'inverese debye screen length-cold' ,kappa_c
write(5,*) 'gama cold region:' ,gama , 'gama hot region:' ,gama2 , 'integration time step : ' , time_step  
write(5,*) 'cold temp' ,T1 ,'hot temp:',T2 , 'seed random number :' , iseed
write(5,*) 'total integration step' , nstep , 'Total Equilibrium Calculation :' , nequ_step
write(5,*) 'Equilibrium time in Equilibrium Calculation: ' , nequ , 'steady state time : ' , nsteady
write(5,*) 'skip of time to calculate g(r) , save coordinates :' , tconfig
 
 
do k=1,natom
   read(1,*) Rx_npc(k) , Ry_npc(k) , Rz_npc(k) , count  !reading lattice structure coordinate from file
   Rx(k) = Rx_npc(k) - Lx*anint(Rx_npc(k)/Lx)
   Ry(k) = Ry_npc(k) - L*anint(Ry_npc(k)/L)
   Rz(k) = Rz_npc(k) - L*anint(Rz_npc(k)/L)
enddo

!write(*,*) count
do i =1,max_bin  
   bin1(i) = 0.0  !initalizing bin
   bin2(i) = 0.0
   bin_eq(i) = 0.0
end do
!call srand(iseed) !seed to generate random number
wtime = omp_get_wtime ( )

do t=run_start+1,nstep
   call force_calc(t,fx,fy,fz,sigma,epot,kappa_h,kappa_c,ch_sys_h,ch_sys_c) !force calculation
   if (t .gt. nequ_step) then
      close(12)   
      close(8)
   endif
   call integration(t,fx,fy,fz,D1, D2,gama,gama2,part_hot,part_cold,time_step) ! integration 
 !-------------------------------------------------------------------------------------------------------------
   if (t .ge. nsteady) then 
      if (mod(t,tconfig) .eq. 0) then
         call pair_correlation_s(t,bine1,bine2,step) !Pair correlqation function in soret effect
         do k=1,max_bin
	    bin1(k) = bin1(k) + bine1(k)
            bin2(k) = bin2(k) + bine2(k)
         enddo
         !do i=1,natom
            !write(9,*) Rx(i) , Ry(i) , Rz(i) , i , t-nequ_step 
         !enddo
      endif
      dens_hot  = dens_hot  + part_hot
      dens_cold = dens_cold + part_cold
   endif
 ! to save configuration of non Equilibrium calculation
   !if ( t .gt. nequ_step) then
      !if (mod(t,tconfig) .eq. 0) then
         !do i=1,natom
            !write(9,*) Rx(i) , Ry(i) , Rz(i) , i , t-nequ_step 
         !enddo
      !endif
   !endif
 !-------------------------------------------------------------------------------------------------------------
   if (t .ge. nequ .and. t .le. nequ_step) then
      if (mod(t,tconfig) .eq. 0) then
         call pair_correlation(t,bine_eq,dens,step_eq) !Pair correlqation function in equilibrium calculation
         do k=1,max_bin
	    bin_eq(k) = bin_eq(k) + bine_eq(k)
         enddo
         !do i=1,natom
            !write(8,*) Rx(i) , Ry(i) , Rz(i) , i , t 
         !enddo
      endif
   endif
! to save configuration of Equilibrium calculation
   !if ( t  .le. nequ_step) then
      !if (mod(t,tconfig) .eq. 0) then
         !do i=1,natom
             !write(8,*) Rx(i) , Ry(i) , Rz(i) , i , t 
         !enddo
      !endif
   !endif
   
   open(7,file = 'final-coordinate1.dat')
   do i=1,natom
      write(7,*) Rx_npc(i) , Ry_npc(i), Rz_npc(i) , i ,t ! saving final last coordinate for some quick checking
   enddo
   close(7)
enddo
2 FORMAT (5(F8.4,'  '))
do i=1,max_bin
   write(2,*) (i-1)*step, bin1(i) , bin2(i) , (i-1)*step_eq , bin_eq(i)
enddo
!write(5,*) 'charge in dilution: ' ,charge_dilute , 'charge in the system: ', charge_system, 'inverese debye screen length' ,kappa 
write(5,*) 'density of hotter side :' ,dens_hot/(Lx/2.0*L*L)/(nstep-nsteady) 
write(5,*) 'density of colder side :' ,dens_cold/(Lx/2.0*L*L)/(nstep-nsteady) 

!do i=1,natom
   !write(7,*) Rx(i) , Ry(i), Rz(i) , i  ! saving final last coordinate for some quick checking
!enddo
wtime = omp_get_wtime ( ) - wtime
write ( *, '(g14.6 )' ) wtime

end program BD

subroutine force_calc(td,fxd,fyd,fzd,sigmad,epotd,kh,kc,ch_sys_h,ch_sys_c)
common / BLOCK1 / Rx, Ry, Rz, Rx_npc, Ry_npc, Rz_npc
common /BLOCK2/ nstep , nequ , nequ_step , nsteady , tconfig
common /BLOCK3/ alpha, alph, alpc 
integer,parameter :: natom=6400
integer :: nstep,  nequ  , nequ_step, nsteady , tconfig
integer i,j,td
real,dimension(natom) :: Rx,Ry,Rz, Rx_npc, Ry_npc, Rz_npc
real,dimension(natom) :: fxd,fyd,fzd
real sigmad
real,parameter :: epsilon1=1.0
real::f,epotd
real,parameter::L=18.00,Lx=39.5
real::Rxij,Ryij,Rzij,Rxi ,Ryi, Rzi
real :: Rijsq,Rcut,sigsq
real :: disij , disexp , dis , dis3 ,coloumb_h ,coloumb_c, ch_sys_h,ch_sys_c, kh,kc ! k = kappa inverse screen length , c_system = charge in system
real :: s ,delta,ssq,svar,alpha,alph,alpc
delta = 1.27
coloumb_h = (ch_sys_h*ch_sys_h*exp(kh))/(1+kh/2.0)/(1+kh/2.0)
coloumb_c = (ch_sys_c*ch_sys_c*exp(kc))/(1+kc/2.0)/(1+kc/2.0)
coloumb_hc = (ch_sys_h*ch_sys_c*exp(kh))/(1+kh/2.0)/(1+kh/2.0)
coloumb_ch = (ch_sys_h*ch_sys_c*exp(kc))/(1+kc/2.0)/(1+kc/2.0)
Rcut = 5*sigmad
Rcutsq = Rcut*Rcut
sigsq = sigmad*sigmad
epotd = 0.0
do i=1,natom   ! set forces to zero
   fxd(i)=0.0
   fyd(i)=0.0  
   fzd(i)=0.0
end do

!$omp parallel &
!!$omp   shared ( h, n ) &
!$omp shared (Rx,Ry,Rz,Rcutsq,epsilon_h,epsilon_c,epsilon_ch,epsilon_hc) &
!$omp private (i,j,Rxij,Ryij,Rzij,f,disij,dis,dis3,disexp,s,Rijsq,Rxi,Ryi,Rzi)

!$omp do reduction(+:epotd,fxd,fyd,fzd)

do i=1,natom-1
   Rxi = Rx(i)
   Ryi = Ry(i)
   Rzi = Rz(i)
   do j=i+1,natom
      Rxij = Rxi-Rx(j)
      Ryij = Ryi-Ry(j)
      Rzij = Rzi-Rz(j)
      Rxij = Rxij - Lx*anint(Rxij/Lx) ! Periodic Boundary Condition
      Ryij = Ryij - L*anint(Ryij/L)
      Rzij = Rzij - L*anint(Rzij/L)
      Rijsq = Rxij*Rxij+Ryij*Ryij+Rzij*Rzij
      if (Rijsq .le. Rcutsq) then   !cutoff test
         disij = sqrt(Rijsq)
         dis   = 1/disij
         dis3 = dis*dis*dis 
         s = disij - 1
         if (td .gt. nequ_step) then
            if (abs(Rxi) .le. Lx/4.0 .and. abs(Rx(j)) .le. Lx/4.0)then
               disexp = exp(-kh*disij)
             if (s .le. 2*delta-1) then 
                 epotd = epotd+ coloumb_h*disexp/disij + alph*s*s !potential calculation
                 f= coloumb_h*disexp*(kh*disij + 1)*dis3 - alph*2*s*dis   !force calculation
             else 
                 epotd = epotd+ coloumb_h*disexp/disij  !potential calculation
                 f= coloumb_h*disexp*(kh*disij + 1)*dis3    !force calculation
             endif
            else if(abs(Rxi) .ge. Lx/4.0 .and. abs(Rx(j)) .ge. Lx/4.0) then
               disexp = exp(-kc*disij)
             if (s .le. 2*delta-1) then 
                 epotd = epotd+ coloumb_c*disexp/disij + alpc*s*s !potential calculation
                 f= coloumb_c*disexp*(kc*disij + 1)*dis3 - alpc*2*s*dis   !force calculation
             else 
                 epotd = epotd+ coloumb_c*disexp/disij  !potential calculation
                 f= coloumb_c*disexp*(kc*disij + 1)*dis3    !force calculation
             endif
            else if(abs(Rxi) .le. Lx/4.0 .and. abs(Rx(j)) .ge. Lx/4.0) then
               disexp = exp(-kc*disij)
             if (s .le. 2*delta-1) then 
                 epotd = epotd+ coloumb_hc*disexp/disij + alph*s*s !potential calculation
                 f= coloumb_hc*disexp*(kc*disij + 1)*dis3 - alph*2*s*dis   !force calculation
             else 
                 epotd = epotd+ coloumb_hc*disexp/disij  !potential calculation
                 f= coloumb_hc*disexp*(kc*disij + 1)*dis3    !force calculation
             endif
            else if(abs(Rxi) .ge. Lx/4.0 .and. abs(Rx(j)) .le. Lx/4.0) then
                 disexp = exp(-kh*disij)
             if (s .le. 2*delta-1) then 
                 epotd = epotd+ coloumb_ch*disexp/disij + alpc*s*s !potential calculation
                 f= coloumb_ch*disexp*(kh*disij + 1)*dis3 - alpc*2*s*dis   !force calculation
             else 
                 epotd = epotd+ coloumb_ch*disexp/disij  !potential calculation
                 f= coloumb_ch*disexp*(kh*disij + 1)*dis3    !force calculation
             endif
            endif
         else
             disexp = exp(-kh*disij)
             if (s .le. 2*delta-1) then 
                 epotd = epotd+ coloumb_h*disexp/disij + alpha*s*s !potential calculation
                 f= coloumb_h*disexp*(kh*disij + 1)*dis3 - alpha*2*s*dis   !force calculation
             else 
                 epotd = epotd+ coloumb_h*disexp/disij  !potential calculation
                 f= coloumb_h*disexp*(kh*disij + 1)*dis3    !force calculation
             endif
             
         endif
         fxd(j)=fxd(j)-(f*Rxij)   !update force
         fyd(j)=fyd(j)-(f*Ryij)
         fzd(j)=fzd(j)-(f*Rzij) 
         fxd(i)=fxd(i)+(f*Rxij)   !update force
         fyd(i)=fyd(i)+(f*Ryij)
         fzd(i)=fzd(i)+(f*Rzij) 
         
      end if
   enddo
enddo

!$omp end do 
!$omp end parallel 

epotd = epotd/natom
if (td .le. nequ_Step) then
   write(12,*)td,epotd
else
   write(13,*)td-nequ_Step,epotd
endif
return
end 

subroutine integration(td,fxd,fyd,fzd,D1t,D2t,gamat,gama2t,part_h,part_c,dt)
common / BLOCK1 / Rx, Ry, Rz , Rx_npc, Ry_npc, Rz_npc
common /BLOCK2/ nstep , nequ , nequ_step , nsteady , tconfig
integer,parameter :: natom=6400
integer :: nstep,  nequ  , nequ_step, nsteady , tconfig
integer i,j,td 
real cold_flux, hot_flux , part_h, part_c
real,dimension(natom) :: Rx,Ry,Rz, Rx_npc, Ry_npc, Rz_npc
real,dimension(natom) :: fxd,fyd,fzd 
real D1t,gamat,D2t,gama2t , dt ! dt = integration time step
real,parameter ::  L=18.00,Lx=39.5
real Rxi,Ryi,Rzi,random, Rxdif,Rydif,Rzdif ,Rdif,R
hot_flux = 0.0 
cold_flux  = 0.0
part_h = 0.0
part_c  = 0.0
do i=1,natom
   if (td .gt. nequ_step) then
    if (abs(Rx(i)) .le. Lx/4.0) then
      call Gaussrandom(td, random)
      R =   (fxd(i)/gama2t)*dt + sqrt(2*D2t*dt)*random
      Rxi = Rx(i) + R
      Rx_npc(i) = Rx_npc(i) + R
      call Gaussrandom(td, random)
      R =   (fyd(i)/gama2t)*dt + sqrt(2*D2t*dt)*random
      Ryi = Ry(i) + R
      Ry_npc(i) = Ry_npc(i) + R
      call Gaussrandom(td, random)
      R =   (fzd(i)/gama2t)*dt + sqrt(2*D2t*dt)*random
      Rzi = Rz(i) + R
      Rz_npc(i) = Rz_npc(i) + R
      Rxi = Rxi - Lx*anint(Rxi/Lx) !Periodic Boundary Condition
      Ryi = Ryi - L*anint(Ryi/L)
      Rzi = Rzi - L*anint(Rzi/L) 
       
      
    else
      call Gaussrandom(td, random)
      R = (fxd(i)/gamat)*dt + sqrt(2*D1t*dt)*random
      Rxi = Rx(i) + R
      Rx_npc(i) = Rx_npc(i) + R
      call Gaussrandom(td, random)
      R = (fyd(i)/gamat)*dt + sqrt(2*D1t*dt)*random
      Ryi = Ry(i) + R
      Ry_npc(i) = Ry_npc(i) + R
      call Gaussrandom(td, random)
      R = (fzd(i)/gamat)*dt + sqrt(2*D1t*dt)*random
      Rzi = Rz(i) + R
      Rz_npc(i) = Rz_npc(i) + R
      Rxi = Rxi - Lx*anint(Rxi/Lx) !Periodic Boundary Condition
      Ryi = Ryi - L*anint(Ryi/L)
      Rzi = Rzi - L*anint(Rzi/L) 
 
    endif
    ! to save configuration of non Equilibrium calculation
    if (mod(td,tconfig) .eq. 0) then
         write(9,*) Rx_npc(i) , Ry_npc(i) , Rz_npc(i) , i , td-nequ_step
    endif
   
    if (abs(Rx(i)) .le. Lx/4.0 .and. abs(Rxi) .ge. Lx/4.0) then
       hot_flux = hot_flux + 1
    endif
    if (abs(Rx(i)) .ge. Lx/4.0 .and. abs(Rxi) .le. Lx/4.0) then
       cold_flux = cold_flux + 1  
    endif  
    Rx(i) = Rxi
    Ry(i) = Ryi
    Rz(i) = Rzi 
    if (abs(Rxi) .le. Lx/4.0) then
       part_h = part_h + 1
    else
       part_c  = part_c  + 1
    endif
   
   else 
       call Gaussrandom(td, random)
       R = (fxd(i)/gama2t)*dt + sqrt(2*D2t*dt)*random
       Rxi = Rx(i) + R
       Rx_npc(i) = Rx_npc(i) + R
       call Gaussrandom(td, random)
       R = (fyd(i)/gama2t)*dt + sqrt(2*D2t*dt)*random
       Ryi = Ry(i) + R
       Ry_npc(i) = Ry_npc(i) + R
       call Gaussrandom(td, random)
       R = (fzd(i)/gama2t)*dt + sqrt(2*D2t*dt)*random
       Rzi = Rz(i) + R
       Rz_npc(i) = Rz_npc(i) + R
       
   
       Rxi = Rxi - Lx*anint(Rxi/Lx) !Periodic Boundary Condition
       Ryi = Ryi - L*anint(Ryi/L)
       Rzi = Rzi - L*anint(Rzi/L) 
 
       Rx(i) = Rxi
       Ry(i) = Ryi
       Rz(i) = Rzi
       ! to save configuration of Equilibrium calculation
       if (mod(td,tconfig) .eq. 0) then
         write(8,*) Rx_npc(i) , Ry_npc(i) , Rz_npc(i) , i , td
      endif
   endif
enddo  
if (td .gt. nequ_step) then
   write(3,*) td-nequ_step, part_h , part_c
   write(4,*) td-nequ_step , hot_flux, cold_flux
endif
return
end 



subroutine pair_correlation_s(tm,hist1,hist2,step)
common / BLOCK1 / Rx, Ry, Rz , Rx_npc, Ry_npc, Rz_npc
common /BLOCK2/ nstep , nequ , nequ_step , nsteady , tconfig  
integer,parameter::natom=6400!no of particles present
integer :: nstep,  nequ  , nequ_step, nsteady , iseed , tconfig 
integer i,j , cold_atom,hot_atom
real,parameter::L=18.00,Lx=39.5!box length
real,dimension(natom)::Rx,Ry,Rz, Rx_npc, Ry_npc, Rz_npc
integer tm,bin,trial_config!bin number
integer,parameter:: max_bin=120
real,dimension(max_bin)::binc1,binc2,hist1,hist2
real Rxi,Ryi,Rzi,Rxij,Ryij,Rzij,Rijsq,disij,density
real Rmax,Rmin,nideal1,nideal2,step
real,parameter:: pi=3.1415927 
trial_conf =200
cold_atom = 0
hot_atom  = 0
do i =1,max_bin 
   binc1(i) = 0.0
   binc2(i) = 0.0
enddo
do i=1,natom-1
   Rxi = Rx(i)
   Ryi = Ry(i)
   Rzi = Rz(i)
   if (abs(Rxi) .le. Lx/4.0) then
   do j = i+1,natom
   if (abs(Rx(j)) .le. Lx/4.0) then
   Rxij = Rxi - Rx(j)
   Ryij = Ryi - Ry(j)
   Rzij = Rzi - Rz(j)
   Rxij = Rxij - Lx*anint(Rxij/Lx)
   Ryij = Ryij - L*anint(Ryij/L)
   Rzij = Rzij - L*anint(Rzij/L)
   Rijsq = Rxij*Rxij+Ryij*Ryij+Rzij*Rzij
   disij = sqrt(Rijsq)
   bin = int(disij/step) + 1
   !write(*,*) bin
   if (bin .le. max_bin) then
       binc1(bin) = binc1(bin) + 2 ! taking the fact distance ij = ji
       !write(*,*) binc1(bin)
   endif
   endif
   enddo
   hot_atom = hot_atom + 1
   else
      do j = i+1,natom
      if (abs(Rx(j)) .ge. Lx/4.0) then
      Rxij = Rxi - Rx(j)
      Ryij = Ryi - Ry(j)
      Rzij = Rzi - Rz(j)
      Rxij = Rxij - Lx*anint(Rxij/Lx)
      Ryij = Ryij - L*anint(Ryij/L)
      Rzij = Rzij - L*anint(Rzij/L)
      Rijsq = Rxij*Rxij+Ryij*Ryij+Rzij*Rzij
      disij = sqrt(Rijsq)
      bin = int(disij/step) + 1
      !write(*,*) bin
      if (bin .le. max_bin) then
         binc2(bin) = binc2(bin) + 2 ! taking the fact distance ij = ji
      endif
      
      endif
      enddo
      cold_atom = cold_atom + 1
      endif
enddo
if (abs(Rx(natom)) .le. Lx/4.0) then
   hot_atom = hot_atom + 1
else
   cold_atom  = cold_atom  + 1
endif
do j =1,max_bin 
   Rmin = real(j-1)*step
   Rmax = Rmin + step
   density = hot_atom/(Lx/2.0*L*L)
   nideal1 = 4*pi*density *(Rmax**3 - Rmin**3)/3
   density = cold_atom/(Lx/2.0*L*L)
   nideal2 = 4*pi*density*(Rmax**3 - Rmin**3)/3
   hist1(j) = binc1(j)/real(hot_atom)/real(trial_conf)/nideal1
   hist2(j) = binc2(j)/real(cold_atom)/real(trial_conf)/nideal2
enddo
!write(*,*) tm, 'I am reaching here'
return
end

subroutine pair_correlation(tm,hist,dens,step)
common / BLOCK1 / Rx, Ry, Rz , Rx_npc, Ry_npc, Rz_npc
common /BLOCK2/ nstep , nequ , nequ_step , nsteady , tconfig 
integer,parameter::natom=6400!no of particles present
integer :: nstep,  nequ  , nequ_step, nsteady , iseed , tconfig 
integer i,j
real,parameter::L=18.00,Lx=39.5!box length
real,dimension(natom)::Rx,Ry,Rz, Rx_npc, Ry_npc, Rz_npc
integer tm,bin,trial_config!bin number
integer,parameter:: max_bin=120
real,dimension(max_bin)::binc,hist
real Rxi,Ryi,Rzi,Rxij,Ryij,Rzij,Rijsq,disij,dens
real Rmax,Rmin,nideal,step
real,parameter:: pi=3.1415927 
trial_conf =100

do i =1,max_bin 
   binc(i) = 0.0
enddo
do i=1,natom-1
   Rxi = Rx(i)
   Ryi = Ry(i)
   Rzi = Rz(i)
   do j = i+1,natom
   Rxij = Rxi - Rx(j)
   Ryij = Ryi - Ry(j)
   Rzij = Rzi - Rz(j)
   Rxij = Rxij - Lx*anint(Rxij/Lx)
   Ryij = Ryij - L*anint(Ryij/L)
   Rzij = Rzij - L*anint(Rzij/L)
   Rijsq = Rxij*Rxij+Ryij*Ryij+Rzij*Rzij
   disij = sqrt(Rijsq)
   bin = int(disij/step) + 1
   !write(*,*) bin
   if (bin .le. max_bin) then
       binc(bin) = binc(bin) + 2 ! taking the fact distance ij = ji
   endif
   enddo
enddo

do j =1,max_bin 
   Rmin = real(j-1)*step
   Rmax = Rmin + step
   nideal = 4*pi*dens*(Rmax**3 - Rmin**3)/3
   hist(j) = binc(j)/real(natom)/real(trial_conf)/nideal
enddo
!write(*,*) tm, 'I am reaching here'
return
end

!--------------------------------------------------------------------------
subroutine Gaussrandom(tm, randomd)
real :: randomd , u1,u2, ran1 , ran2,toss
integer tm
data pi/3.1415927/ 
u1 = ranf()
u2 = ranf()
ran1 = cos(2*pi*u1)*sqrt(-2.*log(u2));
ran2 = cos(2*pi*u2)*sqrt(-2.*log(u1));

toss = ranf()
if (toss > 0.5) then
   randomd = ran1
else
   randomd = ran2
endif
if (ran1-1 .eq. ran1) then 
   randomd = ran2
   !write(*,*) 'ran1 is a inf' , tm
endif
if (ran2-1 .eq. ran2) then 
   randomd = ran1
   !write(*,*) 'ran2 is a inf' , tm
endif
return
end  



!****************************************************************
      REAL function ranf()
      common /rjran/ i3,i2,i1,i0

!     berkeley random number generator
!     range changed to 0 < 1

      INTEGER  I0, I1, I2, I3, J0, J1, J2, J3, K0, K1, K2, K3
      INTEGER  M0, M1, M2, M3, MM
      REAL     T1, T2, T3, T4     

      parameter (m3=647,m2=1442,m1=3707,m0= 373)
      parameter (t4=2.0**48,t3=2.0**36,t2=2.0**24,t1=2.0**12)
      parameter (mm=4096)

!     mm = 2 ** 12

!     the random number is:
!     (i3*t3+i2*t2+i1*t1+i0)/2.0**48
!     the multiplier is:
!     (m3*t3+m2*t2+m1*t1+m0)

      ranf = float(i3)/t1 + float(i2)/t2 + float(i1)/t3 + float(i0)/t4
      if(ranf.ge.0.9999999)  ranf = 0.0

!     multiply i's and m's

      j0 = m0 * i0
      j1 = m0 * i1 + m1 * i0
      j2 = m0 * i2 + m1 * i1 + m2 * i0
      j3 = m0 * i3 + m1 * i2 + m2 * i1 + m3 * i0
      k0 = j0
      k1 = j1 + k0 / mm
      k2 = j2 + k1 / mm
      k3 = j3 + k2 / mm
      i0 = mod(k0,mm)
      i1 = mod(k1,mm)
      i2 = mod(k2,mm)
      i3 = mod(k3,mm)
      return
      end      


      subroutine ranint(istart)

!     initialize random number generator

      INTEGER  IRAN(4), ISTART

      iran(1)=12345 + istart*1000
      iran(2)=12345 + istart*1000
      iran(3)=12345 + istart*1000
      iran(4)=12345 + istart*1000
      call ranset(iran)
      call ranget(iran)

      return
      end


      SUBROUTINE RANSET(IRAN)
      COMMON /RJRAN/ II3,II2,II1,II0

      INTEGER    IRAN(4), MM, NN
      PARAMETER (MM = 4096, NN = 100000)
      INTEGER    II0, II1, II2, II3
      INTEGER    I0, I1, I2, I3, J0, J1, J2, J3

      i3 = iran(1)
      i2 = iran(2)
      i1 = iran(3)
      i0 = iran(4)
      call divide(i3,i2,i1,i0,j3,j2,j1,j0,nn,mm,ii0)
      call divide(j3,j2,j1,j0,i3,i2,i1,i0,nn,mm,ii1)
      call divide(i3,i2,i1,i0,j3,j2,j1,j0,nn,mm,ii2)
      call divide(j3,j2,j1,j0,i3,i2,i1,i0,nn,mm,ii3)
      return
      end


      subroutine ranget(iran)
      common /rjran/ ii3,ii2,ii1,ii0

      INTEGER    mm, m10, m21, m20, m32, m31, m30

      parameter (mm = 100000)
      parameter (m10 = 4096)
      parameter (m21 =  167, m20 = 77216)
      parameter (m32 =    6, m31 = 87194, m30 = 76736)

      INTEGER    iran(4)
      INTEGER    ii3, ii2, ii1, ii0  
      INTEGER    J0, J1, J2, J3, K0, K1, K2, K3

      j0 = ii0 + m10 * ii1 + m20 * ii2 + m30 * ii3
      j1 =                   m21 * ii2 + m31 * ii3
      j2 =                               m32 * ii3
      j3 =                                       0
      k0 = j0
      k1 = j1 + k0 / mm
      k2 = j2 + k1 / mm
      k3 = j3 + k2 / mm
      iran(4) = mod(k0,mm)
      iran(3) = mod(k1,mm)
      iran(2) = mod(k2,mm)
      iran(1) = mod(k3,mm)
      return
      end


      subroutine divide(i3,i2,i1,i0,j3,j2,j1,j0,nn,id,ir)

      INTEGER     I3,I2,I1,I0,J3,J2,J1,J0,ID,IR,NN,K0,K1,K2,K3

!     given the integer i = i0 + nn * (i1 + nn * (i2 + nn * i3))
!     this routine calculates j = i / id and ir = mod(i, id)
!     j is expressed as i, ir is just an integer

      j3 = i3 / id
      k3 = mod(i3, id)
      k2 = k3 * nn + i2
      j2 = k2 / id
      k2 = mod(k2, id)
      k1 = k2 * nn + i1
      j1 = k1 / id
      k1 = mod(k1, id)
      k0 = k1 * nn + i0
      j0 = k0 / id
      ir = mod(k0, id)
      return
      end
!*********************************************************************************
