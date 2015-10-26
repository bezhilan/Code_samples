**********************************************************************
C Program chemotaxis:
C     Simulates the dynamics in a 3D suspension of chemotactic bacteria
C     in thin films.

C Generates the numerical simulation data in 
C "Chaotic dynamics and oxygen transport in thin films of aerotactic bacteria"
C  B. Ezhilan, A. Alizadeh Pahlavan, D. Saintillan, Physics of Fluids, 24 091701 (2012)

C Please see http://stokeslet.ucsd.edu/publications_files/chemotaxis.pdf for details
C of the underlying model 

C Last used : August 10, 2015
C Download the appropriate MPIF and lfftw libraries before compiling 

C To compile use "mpif90 -O3 filename.f -lfftw_mpi -lfftw -lm"
C To execute use "mpirun -np 8 a.out"
C**********************************************************************



C*********************************************************************
C 				MODULE GLOB
C********************************************************************

      Module Glob

      real,dimension(:,:,:,:,:),allocatable :: c,bEW,cnew !The five dimensional array c is \psi(x,p,t), bEW = \beta E + W term in the rotational flux \dot{p}
      real,dimension(:,:,:,:,:),allocatable :: E,f1,S11,S12 !E is the fluid rate of strain
      real,dimension(:,:,:,:,:),allocatable :: F,Fo ! F and Fo are the flux tendency terms

      real,dimension(:,:,:),allocatable :: cremesh


      real,dimension(:,:,:),    allocatable :: Sigma11i,Sigma12i
      real,dimension(:,:,:),    allocatable :: Sigma13i,Sigma22i
      real,dimension(:,:,:),    allocatable :: Sigma23i,Sigma33i

      real,dimension(:,:,:),    allocatable :: Sigma11fi,Sigma12fi
      real,dimension(:,:,:),    allocatable :: Sigma13fi,Sigma22fi
      real,dimension(:,:,:),    allocatable :: Sigma23fi,Sigma33fi


      real,dimension(:,:,:),    allocatable :: fif,fib,fis
      real,dimension(:,:,:),    allocatable :: fifT,fibT,fisT

      real,dimension(:,:,:),allocatable     :: f1tot,f2tot,f3tot


      real,dimension(:,:,:),allocatable     :: ethapi,kpi,nupi
      real,dimension(:,:,:),allocatable     :: ethapT,kpT,nupT

      real,dimension(:,:),allocatable       :: cl,clt

      real,dimension(:,:),allocatable       :: uxwait,uywait,uzwait


      real,dimension(:,:,:,:),  allocatable :: fcompx,fcompy,cwait
      real,dimension(:,:,:,:),  allocatable :: fcompz,ukx,uky,ukz

      real,dimension(:,:,:,:),  allocatable :: csend
 
      real,dimension(:,:,:),    allocatable :: ci,ux,uy,uz,S11i,S12i
      
      real,dimension(:,:,:),    allocatable :: pxi,pyi,pzi,fx,fy,fz,LEI 

      real,dimension(:),    allocatable     :: theta,phi 
      
      real :: Lx,Ly,Lz,time,alpha,dfx,dfr,gamma,dt,dy2,rt
      real :: pi,dx,dy,dz,dph,dth,Shear,Stime,br,Ext,VolumeFraction

      integer :: nx,ny,nz,nph,nth,itime
      integer :: totnx,totny,totnz,req,req1,req2,req3

      real,dimension(:,:,:),allocatable :: ciphi,LET

      integer :: Scoords(1),Dcoords(1),rootR,rankR,rankL
      integer :: rankF,rankO,coords(1),rank,ierr,reorder,tag,cartcomm
      integer :: periods(1),dims(1),numtasks,rc
      integer :: coordsR(1),coordsL(1) 
      integer :: Imin,Imax,Jmin,Jmax,Kmin,Kmax

      real,dimension(:,:,:),allocatable :: fxT,fyT,fzT
      real,dimension(:,:,:,:),allocatable :: ftotcompx,ftotcompy
      real,dimension(:,:,:,:),allocatable :: ftotcompz,fTcompx,fTcompy
      real,dimension(:,:,:,:),allocatable :: fTcompz,ukxT,ukyT,ukzT
      real,dimension(:,:,:,:),allocatable :: ukxtot,ukytot,ukztot
      real,dimension(:,:,:),allocatable :: ciT
      real,dimension(:,:,:),allocatable :: pxiT
      real,dimension(:,:,:),allocatable :: pyiT,pziT
      real,dimension(:,:,:),allocatable :: uxT,uyT,uzT

      integer,allocatable,dimension(:) :: stat

      integer*8 :: plan,iplan

      integer :: lxs,lnx,lnyt,lyst,lsize

      integer :: fftw_fwd, fftw_bkwd,fftw_est,fftw_norm_order

      double complex,dimension(:),allocatable :: fxft,fyft,fzft,work
      double complex,dimension(:),allocatable :: uxft,uyft,uzft

      double complex,dimension(:),allocatable :: skft,skxft,skyft,skzft
      real,dimension(:,:,:,:),allocatable :: sk,skx,sky,skz,skc,skcb
      real,dimension(:,:,:),allocatable :: sT
      real,dimension(:,:,:,:),allocatable :: fk,cs,fko,cso

      double precision,dimension(:,:),allocatable :: psi0
      real,dimension(:,:),allocatable :: psi00

      complex img
      real D, ocr, beta, deltai, tau, xi, l0

C      real :: MPI_REAL

      end module


C*********************************************************************
C 				MAIN PROGRAM
C********************************************************************



      program continuum
    
      use Glob

      implicit none

      include 'mpif.h'

C      integer stat(MPI_STATUS_SIZE)

      allocate (stat(MPI_STATUS_SIZE))

      call MPI_INIT(ierr)
      if (ierr .ne. MPI_SUCCESS) then
         print *,'Error starting MPI program. Terminating.'
         call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
      end if

      call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)

      periods(1) = 1

      dims(1)    = 64

      reorder    = 0

      call MPI_CART_CREATE(MPI_COMM_WORLD,1,dims,periods,reorder,
     $                     cartcomm,ierr)

      call MPI_COMM_RANK(cartcomm, rank, ierr)

      call MPI_CART_COORDS(cartcomm,rank,1,coords,ierr)

      print*,'rank = ',rank


      nx = 64
      ny = 64
      nz = 1

      nz = nz + 2

         
C Define simulation parameters

      call param

C Start MPI

C      Call Init_MPI

C Define initial configuration field

      call initial
C Perform initial step using Euler for the bacteria field 
      call euler
C Perform initial step using Euler for the oxygen field      
      call eulerscalar
C Force the boundary condition      
      call BC    
C skc changes, skx,sky,skz changes, skcb remains initial
C Time marching
      do itime = 1,80000
C Update concentration field
	 rt =0.d0
         call update
         skcb = skc
C Update oxygen field
	 call updatescalar
	 call BC
C        New skc now
	 call stressmagnitude
       
c        print *, rt
      enddo
	 
      Call MPI_Barrier(cartcomm,ierr)

      call MPI_Finalize(ierr)

      end

C**********************************************************************
C Subroutine param:
C     Define simulation parameters
C**********************************************************************

      subroutine param
      
      use Glob

      implicit none

      Shear = 0.0d0
      Stime = 0.0d0

C Box dimensions in x, y and z
      Lx = 50.0d0
C y is the wall-normal direction      
      Ly = 48.0d0
      Lz = 50.0d0

C Number of points in x, y, z and phi, theta

      totnx = 64
      totny = 64
      totnz = 64
      nph = 16
! ideally use more values in the theta direction
      nth = 8

C Dimensionless stresslet 
      alpha = -10.d0
      
C Center of mass and orientation diffusivities
      dfx = 0.52d0
      dfr = 0.12d0

C Jeffery's parameter
c      A = 10.00
c      gamma = (A**2.00-1.00)/(A**2.00+1.00)
      gamma = 1.d0

C Time marching
      dt = 0.004d0
      time = 0.0d0

      pi = 4.d0*atan(1.d0)

      dx = Lx/(totnx)
      dy = Ly/(totny)
      dz = Lz/(totnz)
      dph = 2.d0*pi/real(nph)
      dth = pi/real(nth)

      dy2 = dy

C Dimensionless Stress Groups
      VolumeFraction = 0.1d0

! br and Ext are coefficients of the Brownian and extensile/flow contributions to the extra stress (See Ezhilan POF 2013). Both are set to zero here.

      br    =  0.d0
      Ext   =  0.d0

      img = (0,1) 

C Oxygen Translational Diffusivity
      D = 8.4d0

C Oxygen Consumption Rate
      ocr = 0.0833d0

C Coefficient for Alignment of Particles with Oxygen Gradient
      beta = 0.1d0

C Mean Run Length Coeff
      tau = 5.d0
      l0 = 1.2d0
      xi = 5.d0 
      
C Correlation of Tumbling Direction Coeff
      deltai = 1.d0

      print*,'tau = ',tau

      return
      end



C**********************************************************************
C Subroutine initial:
C     Define initial configuration field
C**********************************************************************

      subroutine initial
 
      use Glob

      implicit none
      include 'mpif.h'


      integer ix,iy,iz,iph,ith,i,j,k,is
      integer nk,ik,idum,ikx,iky,ikz,ind
      character*32 filename1,filename2,filename3,filename4

      real fftheta(nth),ffph(nph)
      real kx,ky,kz,eps
      real ran,cst
      real eps1(40,40,40),eps2(40,40,40)
      real a(8),mean

      allocate(uxwait(nx,2*ny),uywait(nx,2*ny),uzwait(nx,2*ny))
      
      allocate(psi0(128,128),psi00(nph,nth))

      if (rank.eq.0) then

      print*,'rank 0 allocation started'

C Allocate all variables

      allocate(f1tot(totnx,2*totny,totnz),f2tot(totnx,2*totny,totnz))
      allocate(f3tot(totnx,2*totny,totnz))

      allocate(sT(nx,2*ny,nz))
      allocate(fifT(nx,2*ny,nz),fibT(nx,2*ny,nz),fisT(nx,2*ny,nz))
      allocate(fxT(nx,2*ny,nz),fyT(nx,2*ny,nz))
      allocate(fzT(nx,2*ny,nz))
      allocate(ciT(nx,2*ny,nz))
      allocate(pxiT(nx,2*ny,nz))
      allocate(pyiT(nx,2*ny,nz))
      allocate(pziT(nx,2*ny,nz))
      allocate(uxT(nx,2*ny,nz))
      allocate(uyT(nx,2*ny,nz))
      allocate(uzT(nx,2*ny,nz))

      allocate(LET(nx,2*ny,nz))
      allocate(ethapT(nx,2*ny,nz))
      allocate(kpT(nx,2*ny,nz))
      allocate(nupT(nx,2*ny,nz))

      print*,'rank 0 allocation done'

      end if

      call MPI_Barrier(cartcomm,ierr)



      allocate(fif(nx,2*ny,nz),fib(nx,2*ny,nz),fis(nx,2*ny,nz))
      allocate(cl(nph,nth),clt(nph,nth))

      allocate(LEI(nx,2*ny,nz))
      allocate(ethapi(nx,2*ny,nz),cwait(nx,2*ny,nph,nth))
      allocate(csend(nx,2*ny,nph,nth))


      allocate(nupi(nx,2*ny,nz))

      allocate(kpi(nx,2*ny,nz))

      allocate(E(3,3,nx,2*ny,nz))

      allocate(c(nx,2*ny,nz,nph,nth),ci(nx,2*ny,nz))
      allocate(S11i(nx,2*ny,nz))
      allocate(S12i(nx,2*ny,nz))

      allocate(fx(nx,2*ny,nz),fy(nx,2*ny,nz),fz(nx,2*ny,nz))
      allocate(ux(nx,2*ny,nz),uy(nx,2*ny,nz),uz(nx,2*ny,nz))
      allocate(phi(nph),theta(nth))

      allocate(F(nx,2*ny,nz,nph,nth),Fo(nx,2*ny,nz,nph,nth))

      allocate(bEW(3,3,nx,2*ny,nz))
      allocate(ciphi(nx,2*ny,nz))
      allocate(pxi(nx,2*ny,nz),pyi(nx,2*ny,nz),pzi(nx,2*ny,nz))
      allocate(fcompx(2,nx,2*ny,nz),fcompy(2,nx,2*ny,nz))
      allocate(fcompz(2,nx,2*ny,nz))
      allocate(ukx(2,nx,2*ny,nz),uky(2,nx,2*ny,nz),ukz(2,nx,2*ny,nz))

      allocate(f1(nx,2*ny,nz,nph,nth))

      allocate(sk(2,nx,2*ny,nz),skx(2,nx,2*ny,nz),sky(2,nx,2*ny,nz))
      allocate(skz(2,nx,2*ny,nz),skc(2,nx,2*ny,nz),skcb(2,nx,2*ny,nz))
      allocate(fk(2,nx,2*ny,nz),cs(2,nx,2*ny,nz))
      allocate(fko(2,nx,2*ny,nz),cso(2,nx,2*ny,nz))

      call MPI_Barrier(cartcomm,ierr)

      print*,'all ranks allocation done'


      do iph = 1,nph
         phi(iph) = (iph-1.0)/nph*2*pi
      end do
 
      do ith = 1,nth
         theta(ith) = (ith-0.5d0)/nth*pi
      end do

C  Consider nearly uniform/isotropic suspension
      eps = 0.015d0


C Build random initial distribution
      c = 1.d0

C Initialize c
       
 
C Number of spatial modes in x and y directions
      nk = 8
C Get amplitudes of the various modes
      idum = -854881
      do ikz = 1,2*nk
         do iky = 1,2*nk
            do ikx = 1,2*nk
               eps1(ikx,iky,ikz) = eps*2.0*(ran(idum)-0.50)
               eps2(ikx,iky,ikz) = eps*2.0*(ran(idum)-0.50)
            enddo
         enddo
      enddo
      do ind = 1,8
         a(ind) = 2.0*(ran(idum)-0.5)
      enddo

!3rd Order Polynomial for fi

      do iph = 1,nph
         ffph(iph) = (a(1)+a(2)*cos(phi(iph))+
     $   a(3)*cos(phi(iph))**2.0+a(4)*cos(phi(iph))**3.0)*
     $   (a(5)+a(6)*sin(phi(iph))
     $   +a(7)*sin(phi(iph))**2.0+a(8)*sin(phi(iph))**3.0)
      enddo

!3rd Order Polynomial for Theta

      do ith = 1,nth
         fftheta(ith) = (a(1)+a(2)*cos(theta(ith))
     $   +a(3)*cos(theta(ith))**2.0
     $   +a(4)*cos(theta(ith))**3.0)*
     $   (a(5)+a(6)*sin(theta(ith))
     $   +a(7)*sin(theta(ith))**2.0
     $   +a(8)*sin(theta(ith))**3.0)
      enddo
 
C For each mode, find angle dependence

      call MPI_CART_COORDS(cartcomm,rank,1,Scoords,ierr)

      Kmin = Scoords(1)*(nz-2) + 1
      Kmax = Scoords(1)*(nz-2) + nz-2

      do ikz = 1,2*nk
         if(rank.eq.0) then
            print *,'wave number = ',ikz
         end if
         do iky = 1,2*nk
            do ikx = 1,2*nk
          
            do ith = 1,nth
               do iph = 1,nph
                  do iz = Kmin,Kmax
                     do iy = 1,2*ny
                        do ix = 1,nx




              k = iz + 1 - Scoords(1)*(nz-2)
             
!      if (sqrt((ikx-nx)**2.00+(iky-ny)**2.00+
!     $   (ikz-nz)**2.00).gt.7.00) then
      c(ix,iy,k,iph,ith) = c(ix,iy,k,iph,ith)+
     $(eps1(ikx,iky,ikz)*cos(2.0*pi*(ikx-nk)*(ix-1.0)
     $/real(totnx)+2.0*pi*(iky-nk)*(iy-1.0)/real(2*totny)
     $+2.0*pi*(ikz-nk)*(iz-1.0)/real(totnz))
     $+eps2(ikx,iky,ikz)*sin(2.0*pi*(ikx-nk)*(ix-1.0)
     $/real(totnx)+2.0*pi*(iky-nk)*(iy-1.0)/real(2*totny)
     $+2.0*pi*(ikz-nk)*(iz-1.0)/real(totnz)))
     $*ffph(iph)*fftheta(ith)
!      endif

c              end if

            enddo
            enddo
            enddo
            enddo
            enddo 
       
           enddo
         enddo
      enddo


      allocate(cnew(nx,2*ny,nz,nph,nth))


      Do ith = 1,nth
         Do iph = 1,nph
            Do iz = 2,nz-1
               Do iy = 1,2*ny
                  Do ix = 1,nx

                     i = 2*ny-iy+2
                     if (i.gt.2*ny) i = i - 2*ny
                     j = nph-iph+2
                     if (j.gt.nph) j = j - nph

                    cnew(ix,iy,iz,iph,ith) = c(ix,iy,iz,iph,ith) +
     $                                       c(ix,i,iz,j,ith)


                  End Do
               End Do
            End Do
         End Do
      End Do
     

      c = cnew

      deallocate(cnew)

C Make sure mean concentration is 1
 
      call intphi(c,ci)

      mean = 0.d0
      do iz = 2,nz-1
         do iy = 1,2*ny
            do ix = 1,nx
               mean = mean + ci(ix,iy,iz)
            enddo
         enddo
      enddo


      mean = mean/real(nx*(2*ny)*(nz-2))

      print*,'mean = ',mean

      
      c = c / mean  

c      call reflect

      call intphi(c,ci)


      call MPI_Barrier(cartcomm,ierr)
 

C COPY TO GHOST CELLS FROM NEIGHBORS


C******************************************************************
C	        COORDINATES AND RANKS OF NEIGHBORS
 
         call MPI_CART_COORDS(cartcomm,rank,1,Scoords,ierr)

         coordsR(:) = Scoords(:)
         coordsR(1) = Scoords(1) + 1

         call MPI_CART_RANK(cartcomm, coordsR, rankR, ierr)


         coordsL(:) = Scoords(:)
         coordsL(1) = Scoords(1) - 1

         call MPI_CART_RANK(cartcomm, coordsL, rankL, ierr)


C********************************************************************

         tag = 0
         Call MPI_IRECV(cwait(:,:,:,:),nx*2*ny*nph*nth,MPI_REAL,rankL,
     $                          tag,cartcomm,req,ierr)
 
         csend(:,:,:,:) = c(:,:,nz-1,:,:) 

         tag = 0
         Call MPI_SEND(csend(:,:,:,:),nx*2*ny*nph*nth,MPI_REAL,rankR,
     $                          tag,cartcomm,ierr)


         call MPI_WAIT(req,stat,ierr) 

         c(:,:,1,:,:) = cwait(:,:,:,:)

C     	RECEIVE DATA FROM THE LEFT NEIGHBOR 

C*********************************************************************

C********************************************************************

         tag = 1
         Call MPI_IRECV(cwait(:,:,:,:),nx*2*ny*nph*nth,MPI_REAL,rankR,
     $                          tag,cartcomm,req,ierr)
 
         csend(:,:,:,:) = c(:,:,2,:,:) 


         tag = 1
         Call MPI_SEND(csend(:,:,:,:),nx*2*ny*nph*nth,MPI_REAL,rankL,
     $                          tag,cartcomm,ierr)

         call MPI_WAIT(req,stat,ierr)

         c(:,:,nz,:,:) = cwait(:,:,:,:)


      call intphi(c,ci)


      fftw_fwd  = -1
      fftw_bkwd =  1
      fftw_est  =  0
      fftw_norm_order = 0


    
c      dy2 = dy*sqrt(1.d0+(Shear*Stime)**2.d0)

C*********************************************************************************
C*********************************************************************************
C*********************************************************************************


       Call fftw3d_f77_mpi_create_plan(plan,cartcomm,totnx,2*totny,
     #                                 totnz,fftw_fwd,fftw_est)

       Call fftw3d_f77_mpi_create_plan(iplan,cartcomm,totnx,2*totny,
     $                                 totnz,fftw_bkwd,fftw_est)

 
       Call fftwnd_f77_mpi_local_sizes(plan,lnx,lxs,lnyt,lyst,lsize)


       do is = 1,numtasks

          if (rank.eq.is-1) then
             print*,'rank = ',rank
             print*,'lxs  = ',lxs
             print*,'lnx  = ',lnx
             print*,'lsize =',lsize
          end if
          call MPI_Barrier(cartcomm,ierr)
       end do

      allocate(fxft(lsize),fyft(lsize),fzft(lsize))
      allocate(uxft(lsize),uyft(lsize),uzft(lsize))
      allocate(work(lsize),skft(lsize),skxft(lsize))
      allocate(skyft(lsize),skzft(lsize))


      cst = real(totnx*2*totny*totnz)

C	INITIAL CONDITION FOR THE SCALAR FIELD

      sk             = 0.d0
      sk(1,:,1,:)    = 1.d0
      sk(1,:,ny+1,:) = 1.d0
      skc=sk
      skcb=sk 
      skx=0.d0
      sky=0.d0
      skz=0.d0
   
       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx

                skft(((iz-2)*2*ny+(iy-1))*nx+ix) = sk(1,ix,iy,iz)
         
             end do
          end do
       end do
       

       work = 0
       Call fftwnd_f77_mpi(plan,1,skft,work,0,fftw_norm_order)
      
 
       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx

              sk(1,ix,iy,iz) = 
     $        real(skft(((iz-2)*2*ny+(iy-1))*nx+ix))/cst 

              sk(2,ix,iy,iz) = 
     $        aimag(skft(((iz-2)*2*ny+(iy-1))*nx+ix))/cst
        
             end do
          end do
       end do
 
    
      return
      end


C**********************************************************************
C Subroutine intphi:
C     Calculates the integral over sphere of orientation using trapezoidal rule
C Input = 5D matrix; Output = 3D matrix
C**********************************************************************

      subroutine intphi(cll,clli)

      use Glob

      implicit none

      integer ix,iy,iz,iph,ith

      real cll(nx,2*ny,nz,nph,nth),clli(nx,2*ny,nz)

   
      clli = 0.0

      do ith = 1,nth
         ciphi = 0.0
         do iph = 1,nph
            do iz = 2,nz-1
               do iy = 1,2*ny
                  do ix = 1,nx
                     ciphi(ix,iy,iz) = ciphi(ix,iy,iz)
     $               + cll(ix,iy,iz,iph,ith)*dph
                  enddo
               enddo
            enddo
         enddo
         clli = clli + ciphi*sin(theta(ith))*dth
      enddo
   

      return
      end


C**********************************************************************
C Subroutine intphip:
C     Calculates the integral over phi using trapezoidal rule
C     Input = 2D matrix; Output = scalar
C**********************************************************************

      subroutine intphip(cll,clli)

      use Glob

      implicit none

      integer iph,ith

      real cll(nph,nth),clli,ciphip

   
      clli = 0.0

      do ith = 1,nth
         ciphip = 0.0
         do iph = 1,nph
            ciphip = ciphip + cll(iph,ith)*dph
         enddo
         clli = clli + ciphip*sin(theta(ith))*dth
      enddo
   

      return
      end


C**********************************************************************
C Subroutine reflect:
C Reflects the probability field across the film boundary (y=H)
C**********************************************************************

      subroutine reflect

      use glob

      implicit none

      integer ix,iy,iph,i,j,ith,iz
      real cst

C       REFLECTION (symmetric about ny+1)

      Do ith = 1,nth
         Do iph = 1,nph
            Do iz = 2,nz-1
               Do iy = ny+2,2*ny
                  Do ix = 1,nx

                     i = 2*ny-iy+2
                     j = nph-iph+2
                     if (j.gt.nph) j = j - nph

                    c(ix,iy,iz,iph,ith) = c(ix,i,iz,j,ith)


                  End Do
               End Do
            End Do
         End Do
      End Do



      return
      end

C**********************************************************************
C Subroutine Imposing Boundary Conditions:
C**********************************************************************

      subroutine BC

      use glob

      implicit none

      integer ix,iy,iph,iz,ith
      real mean,ky,cst

C       CONSTANT CONCENTRATION OF OXYGEN AT THE BOUNDARIES OF THIN FILM

       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx

C We convert the 3D array "sk" to a 1D vector "skft" because at the time I was during it, fftwnd_f77_mpi wouldn't handle multi-dimensional arrays
C There is a much easier way of doing it now.

                skft(((iz-2)*2*ny+(iy-1))*nx+ix) =
     $          (sk(1,ix,iy,iz)+img*sk(2,ix,iy,iz))
        
             end do
          end do
       end do

       work = 0

C Inverse fourier transform the oxygen concentration field to get it in real space. 
C (NOTE: IPLAN is for inverse fourier transform : FOURIER -> REAL space)
       Call fftwnd_f77_mpi(iplan,1,skft,work,0,fftw_norm_order)

       cst = real(totnx*2*totny*totnz) !you have to divide the array volume. See the fftwnd_ documentation

       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx

C Transform the real oxygen field from 1D array into 3D array form

                sk(1,ix,iy,iz) = real(skft(((iz-2)*2*ny+(iy-1))*nx+ix))
       
             end do
          end do
       end do

C Set the dirichlet boundary condition for the oxygen concentration at the boundary
       sk(1,:,1,:)     = 1.d0
       sk(1,:,ny+1,:)  = 1.d0
C Force the imaginary part to zero (as it should be)
       sk(2,:,:,:)     = 0.d0
C Make a copy of the oxygen concentration field at the current time step
       skc = sk

       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx

                skft(((iz-2)*2*ny+(iy-1))*nx+ix) = sk(1,ix,iy,iz)
         
             end do
          end do
       end do

C Transform to fourier space      
       work = 0
       Call fftwnd_f77_mpi(plan,1,skft,work,0,fftw_norm_order)
      
 
       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx
C Transform the real oxygen field from 1D array into 3D array form 
C and separate the real and imaginary parts
              sk(1,ix,iy,iz) = 
     $        real(skft(((iz-2)*2*ny+(iy-1))*nx+ix))/cst 

              sk(2,ix,iy,iz) = 
     $        aimag(skft(((iz-2)*2*ny+(iy-1))*nx+ix))/cst
        
             end do
          end do
       end do
 

      return
      end


C**********************************************************************
C Subroutine stressmagnitude:
C     Get shear stress magnitude in suspension
C**********************************************************************

      subroutine stressmagnitude
    
      use Glob

      implicit none
      include 'mpif.h'


      real TotEntropy,TEntropy,TEntropyT,Entropy,EntropyZ,EntropyY
      integer ix,iy,iz,iph,ith,ifile,iorth
      integer i,j,k,is !stat(MPI_STATUS_SIZE),is
      character*32 filename1,filename2,filename3,filename4
      character*32 filename5,filename6,filename7,filename8
      character*32 filename9,filename10,filename11
      character*32 filename12,filename13,filename14
      character*32 filename15,filename16,filename17,fn


      real csum,kpSum,nupSum,ethapSum,xorth,dxm,dxp,p(3),mean,meanT
      real fifSum,fibSum,fisSum,fifY,fibY,fisY,fifZ,fibZ,fisZ
      real OFlux1,OFlux2,OFlux1T,OFlux2T,OFlux1TT,OFlux2TT
      integer sh
      real thetaa,phii!psi00(nph,nth)
c      double precision psi0(128,128)

      rootR = 0


C Get off-diagonal term in stress tensor

      do ith = 1,nth
         do iph = 1,nph
            do iz = 1,nz
               do iy = 1,ny
                  do ix = 1,nx
C p is the unit vector on the sphere of orientation

                     p(1) = sin(theta(ith))*cos(phi(iph))
                     p(2) = sin(theta(ith))*sin(phi(iph))
                     p(3) = cos(theta(ith))

C \psi * p1*p2

                 f1(ix,iy,iz,iph,ith) = c(ix,iy,iz,iph,ith)*p(1)*p(2)

                 enddo
               enddo
            enddo
         enddo
      enddo

C D_{12}  = \int \psi * p1 * p2   dp
      call intphi(f1,S12i)

      do ith = 1,nth
         do iph = 1,nph
            do iz = 1,nz
               do iy = 1,ny
                  do ix = 1,nx

                     p(1) = sin(theta(ith))*cos(phi(iph))
                     p(2) = sin(theta(ith))*sin(phi(iph))
                     p(3) = cos(theta(ith))

C \psi * ( p1^2 -1/3)
                 f1(ix,iy,iz,iph,ith) = c(ix,iy,iz,iph,ith)
     $          *(p(1)*p(1)-1.d0/3.d0)

                 enddo
               enddo
            enddo
         enddo
      enddo

C D_{11}  = \int \psi * ( p1^2 -1/3)   dp

      call intphi(f1,S11i)

      if (rank.eq.0) then
          write(7,*) time,S11i(32,32,5),S12i(32,32,5)
      end if


     
      if (25.0*(itime/25).eq.1.0*itime.OR.itime.eq.1) then
         ifile = itime/25

         if (rank.eq.0) then
            print *,'file',ifile
         end if


C Calculate Flux

         OFlux1T = 0.d0
         OFlux2T = 0.d0
        
         Do iz = 2,nz-1
            OFlux1 = 0.d0
            OFlux2 = 0.d0
            Do ix = 1,nx

            OFlux1 = OFlux1 + (sky(1,ix,1,iz))*dx
            OFlux2 = OFlux2 + (sky(1,ix,ny+1,iz))*dx
    
            End Do
            OFlux1T = OFlux1T + OFlux1*dz
            OFlux2T = OFlux2T + OFlux2*dz
         End Do


         Call MPI_REDUCE (OFlux1T,OFlux1TT,1,MPI_REAL,MPI_SUM,
     $                   rootR,cartcomm,ierr)


         Call MPI_REDUCE (OFlux2T,OFlux2TT,1,MPI_REAL,MPI_SUM,
     $                   rootR,cartcomm,ierr)


         if (rank.eq.0) then

             write(11,*) time, OFlux1TT, OFlux2TT, 
     $       OFlux1TT-OFlux2TT, OFlux1TT+OFlux2TT

          end if


         call intphi(c,ci)


         call meanp

         
C Total Entropy
C See Saintillan and Shelley POF (2008) for the definition

        f1 = 4.d0*pi*c*(log(4.d0*pi*abs(c)))

        call intphi(f1,LEI)

        Do ith = 1,nth
           Do iph = 1,nph
              csum = 0.d0
              Do iz = 2,nz-1
                 Do iy = 1,2*ny
                    Do ix = 1,nx

                       csum = csum + c(ix,iy,iz,iph,ith)

                    End Do
                 End Do
              End Do
              cl(iph,ith) = csum/real(nx*2*ny*(nz-2))
           End Do
        End Do

        Call MPI_REDUCE (cl,clt,nph*nth,MPI_REAL,MPI_SUM,
     $                   rootR,cartcomm,ierr)



      mean = 0.d0
      do iz = 2,nz-1
         do iy = 1,ny+1
            do ix = 1,nx
               mean = mean + ci(ix,iy,iz)
            enddo
         enddo
      enddo


      mean = mean/real(nx*(ny+1)*(nz-2))

 
        Call MPI_REDUCE (mean,meanT,1,MPI_REAL,MPI_SUM,
     $                   rootR,cartcomm,ierr)

      if (rank.eq.0) write(8,*) time, meanT/numtasks
     


C*****************************************************************

C     SEND DATA FROM ALL PROCESSORS TO RANK ZERO 

      rootR = 0


      Do is = 1,numtasks-1

         if (rank.eq.is.OR.rank.eq.0) then


         if (rank.eq.rootR) then
            tag   = 0
            Call MPI_IRECV(ciT(:,:,:),nx*2*ny*nz,MPI_REAL,
     $                         is,tag,cartcomm,req1,ierr)

            tag   = 1
            Call MPI_IRECV(sT(:,:,:),nx*2*ny*nz,MPI_REAL,
     $                         is,tag,cartcomm,req2,ierr)

            tag   = 9
            Call MPI_IRECV(LET(:,:,:),nx*2*ny*nz,MPI_REAL,
     $                         is,tag,cartcomm,req3,ierr)

         end if


 
         if (rank.eq.is) then

            tag   = 0
            Call MPI_SEND(ci(:,:,:),nx*2*ny*nz,MPI_REAL,rootR,
     $                          tag,cartcomm,ierr)

            tag   = 1
            Call MPI_SEND(skc(1,:,:,:),nx*2*ny*nz,MPI_REAL,rootR,
     $                          tag,cartcomm,ierr)

            tag   = 9
            Call MPI_SEND(LEI(:,:,:),nx*2*ny*nz,MPI_REAL,rootR,
     $                          tag,cartcomm,ierr)
c
         end if

         if (rank.eq.rootR) then

            call MPI_WAIT(req1,stat,ierr)
            call MPI_WAIT(req2,stat,ierr)
            call MPI_WAIT(req3,stat,ierr)

         end if


C	RANK ZERO RECEIVES DATA FROM ALL PROCESSORS ONE BY ONE


         if (rank.eq.0) then

            call MPI_CART_COORDS(cartcomm,is,1,Scoords,ierr)

C	RANK ZERO COPIES THE RECEIVED DATA INTO A GLOBAL SIZE ARRAY

            Do iz = 2,nz-1
               Do iy = 1,2*ny
                  Do ix = 1,nx

                     k = Scoords(1)*(nz-2) + iz - 1

                     f1tot(ix,iy,k)    = ciT(ix,iy,iz) 

                     f2tot(ix,iy,k)    = LET(ix,iy,iz)

                     f3tot(ix,iy,k)    = sT(ix,iy,iz)

                 End Do
               End Do
            End Do

         end if

         end if

         Call MPI_Barrier(cartcomm,ierr)

      end do

      Call MPI_Barrier(cartcomm,ierr)


C	INCLUDING THE DATA FROM RANK ZERO ITSELF

      if (rank.eq.rootR) then

         call MPI_CART_COORDS(cartcomm,rootR,1,Scoords,ierr)

         Do iz = 2,nz-1
            Do iy = 1,2*ny
               Do ix = 1,nx

                  k = Scoords(1)*(nz-2) + iz - 1

                  f1tot(ix,iy,k)  = ci(ix,iy,iz)

                  f2tot(ix,iy,k)  = LEI(ix,iy,iz)

                  f3tot(ix,iy,k)  = skc(1,ix,iy,iz)

              End Do
            End Do
         End Do



C       ENTROPY

         TEntropy = 0.d0
         do ix = 1,totnx
            EntropyY = 0.d0
            do iy = 1,totny
               EntropyZ = 0.d0
               do iz = 1,totnz
                  Entropy = f2tot(ix,iy,iz)
                  EntropyZ = EntropyZ + Entropy*dz
               end do
               EntropyY = EntropyY + EntropyZ*dy
            end do
            TEntropy = TEntropy + EntropyY*dx
         end do

         write(56,*) time, TEntropy/real(Lx*Ly*Lz)

C Output con = bacteria concentration / number density field
         write(filename1,'(A3,I4.4,A4)') 'con',ifile,'.txt'
         open(unit=5,file=filename1,status='unknown')
         write(5,*) 'Variables = "x", "y", "z" , "c"'
         write(5,*) 'Zone I= ',totnx, ' J= ', totny+1, ' K= ', totnz

         do iz = 1,totnz
            do iy = 1,totny+1
               do ix = 1,totnx

            write(5,*)  (ix-1)*(Lx/totnx),(Ly/totny)*(iy-1),
     $          (Lz/totnz)*(iz-1),f1tot(ix,iy,iz)

               enddo
            enddo
         enddo
         close(5)

c Output mix = oxygen concentration / number density field

         write(filename1,'(A3,I4.4,A4)') 'mix',ifile,'.txt'
         open(unit=5,file=filename1,status='unknown')
         write(5,*) 'Variables = "x", "y", "z" , "s"'
         write(5,*) 'Zone I= ',totnx, ' J= ', totny+1, ' K= ', totnz

         do iz = 1,totnz
            do iy = 1,totny+1
               do ix = 1,totnx

           write(5,*)  (ix-1)*(Lx/totnx),(Ly/totny)*(iy-1),
     $          (Lz/totnz)*(iz-1),f3tot(ix,iy,iz)

               enddo
            enddo
         enddo
         close(5)

c Output lcon(ith,iph) is the spatially averaged orientation distribution \int_{V} \psi dV. 

         write(filename9,'(A4,I4.4,A4)') 'lcon',ifile,'.txt'
         open(unit=5,file=filename9,status='unknown')
         write(5,*) 'Variables = "phi", "theta", "lcon" '
         write(5,*) 'Zone I = ',nph,' J= ', nth

         do ith = 1,nth

            do iph = 1,nph

               write(5,*) phi(iph),theta(ith),clt(iph,ith)/numtasks

            enddo

         enddo
         close(5)

      end if

      call MPI_Barrier(cartcomm,ierr)




      Do is = 1,numtasks-1

         if (rank.eq.is.OR.rank.eq.0) then


         if (rank.eq.rootR) then

            tag   = 6
            Call MPI_IRECV(uxT(:,:,:),nx*2*ny*nz,MPI_REAL,
     $                         is,tag,cartcomm,req1,ierr)

            tag   = 7
            Call MPI_IRECV(uyT(:,:,:),nx*2*ny*nz,MPI_REAL,
     $                         is,tag,cartcomm,req2,ierr)

            tag   = 8
            Call MPI_IRECV(uzT(:,:,:),nx*2*ny*nz,MPI_REAL,
     $                         is,tag,cartcomm,req3,ierr)

         end if


         if (rank.eq.is) then


            tag   = 6
            Call MPI_SEND(ux(:,:,:),nx*2*ny*nz,MPI_REAL,rootR,
     $                          tag,cartcomm,ierr)

            tag   = 7
            Call MPI_SEND(uy(:,:,:),nx*2*ny*nz,MPI_REAL,rootR,
     $                          tag,cartcomm,ierr)

            tag   = 8
            Call MPI_SEND(uz(:,:,:),nx*2*ny*nz,MPI_REAL,rootR,
     $                          tag,cartcomm,ierr)

c
         end if

         if (rank.eq.rootR) then

            call MPI_WAIT(req1,stat,ierr)
            call MPI_WAIT(req2,stat,ierr)
            call MPI_WAIT(req3,stat,ierr)

         end if

C	RANK ZERO RECEIVES DATA FROM ALL PROCESSORS ONE BY ONE


         if (rank.eq.0) then

            call MPI_CART_COORDS(cartcomm,is,1,Scoords,ierr)

C	RANK ZERO COPIES THE RECEIVED DATA INTO A GLOBAL SIZE ARRAY

            Do iz = 2,nz-1
               Do iy = 1,2*ny
                  Do ix = 1,nx

                     k = Scoords(1)*(nz-2) + iz - 1

                     f1tot(ix,iy,k)   = uxT(ix,iy,iz)
                     f2tot(ix,iy,k)   = uyT(ix,iy,iz)
                     f3tot(ix,iy,k)   = uzT(ix,iy,iz)


                 End Do
               End Do
            End Do

         end if

         end if

         Call MPI_Barrier(cartcomm,ierr)

      end do

      Call MPI_Barrier(cartcomm,ierr)


C	INCLUDING THE DATA FROM RANK ZERO ITSELF

      if (rank.eq.rootR) then

         call MPI_CART_COORDS(cartcomm,rootR,1,Scoords,ierr)

         Do iz = 2,nz-1
            Do iy = 1,2*ny
               Do ix = 1,nx

                  k = Scoords(1)*(nz-2) + iz - 1

                  f1tot(ix,iy,k)   = ux(ix,iy,iz)
                  f2tot(ix,iy,k)   = uy(ix,iy,iz)
                  f3tot(ix,iy,k)   = uz(ix,iy,iz)

 
              End Do
            End Do
         End Do

C Output the velocity field
         write(filename1,'(A3,I4.4,A4)') 'vel',ifile,'.txt'
         open(unit=5,file=filename1,status='unknown')
         write(5,*) 'Variables = "x", "y", "z" , "ux", "uy", "uz"'
         write(5,*) 'Zone I= ',totnx, ' J= ', totny+1, ' K= ', totnz

         do iz = 1,totnz
            do iy = 1,totny+1
               do ix = 1,totnx

            write(5,*)  (ix-1)*(Lx/totnx),
     $          (Ly/totny)*(iy-1),(Lz/totnz)*(iz-1),f1tot(ix,iy,iz)
     $          ,f2tot(ix,iy,iz),f3tot(ix,iy,iz)
               enddo

            enddo
         enddo
         close(5)


      end if

      call MPI_Barrier(cartcomm,ierr)



      Do is = 1,numtasks-1

         if (rank.eq.is.OR.rank.eq.0) then


         if (rank.eq.rootR) then

            tag   = 6
            Call MPI_IRECV(uxT(:,:,:),nx*2*ny*nz,MPI_REAL,
     $                         is,tag,cartcomm,req1,ierr)

            tag   = 7
            Call MPI_IRECV(uyT(:,:,:),nx*2*ny*nz,MPI_REAL,
     $                         is,tag,cartcomm,req2,ierr)

            tag   = 8
            Call MPI_IRECV(uzT(:,:,:),nx*2*ny*nz,MPI_REAL,
     $                         is,tag,cartcomm,req3,ierr)

         end if


         if (rank.eq.is) then


            tag   = 6
            Call MPI_SEND(pxi(:,:,:),nx*2*ny*nz,MPI_REAL,rootR,
     $                          tag,cartcomm,ierr)

            tag   = 7
            Call MPI_SEND(pyi(:,:,:),nx*2*ny*nz,MPI_REAL,rootR,
     $                          tag,cartcomm,ierr)

            tag   = 8
            Call MPI_SEND(pzi(:,:,:),nx*2*ny*nz,MPI_REAL,rootR,
     $                          tag,cartcomm,ierr)

c
         end if

         if (rank.eq.rootR) then

            call MPI_WAIT(req1,stat,ierr)
            call MPI_WAIT(req2,stat,ierr)
            call MPI_WAIT(req3,stat,ierr)

         end if

C	RANK ZERO RECEIVES DATA FROM ALL PROCESSORS ONE BY ONE


         if (rank.eq.0) then

            call MPI_CART_COORDS(cartcomm,is,1,Scoords,ierr)

C	RANK ZERO COPIES THE RECEIVED DATA INTO A GLOBAL SIZE ARRAY

            Do iz = 2,nz-1
               Do iy = 1,2*ny
                  Do ix = 1,nx

                     k = Scoords(1)*(nz-2) + iz - 1

                     f1tot(ix,iy,k)   = uxT(ix,iy,iz)
                     f2tot(ix,iy,k)   = uyT(ix,iy,iz)
                     f3tot(ix,iy,k)   = uzT(ix,iy,iz)


                 End Do
               End Do
            End Do

         end if

         end if

         Call MPI_Barrier(cartcomm,ierr)

      end do

      Call MPI_Barrier(cartcomm,ierr)


C	INCLUDING THE DATA FROM RANK ZERO ITSELF

      if (rank.eq.rootR) then

         call MPI_CART_COORDS(cartcomm,rootR,1,Scoords,ierr)

         Do iz = 2,nz-1
            Do iy = 1,2*ny
               Do ix = 1,nx

                  k = Scoords(1)*(nz-2) + iz - 1

                  f1tot(ix,iy,k)   = pxi(ix,iy,iz)
                  f2tot(ix,iy,k)   = pyi(ix,iy,iz)
                  f3tot(ix,iy,k)   = pzi(ix,iy,iz)

 
              End Do
            End Do
         End Do


C mnp is an abbreviation for "mean p". It is the first orientational moment m = \int \psi p dp

         write(filename1,'(A3,I4.4,A4)') 'mnp',ifile,'.txt'
         open(unit=5,file=filename1,status='unknown')
         write(5,*) 'Variables = "x", "y", "z" , "px", "py", "pz"'
         write(5,*) 'Zone I= ',totnx, ' J= ', totny+1, ' K= ', totnz

         do iz = 1,totnz
            do iy = 1,totny+1
               do ix = 1,totnx

            write(5,*)  (ix-1)*(Lx/totnx),
     $          (Ly/totny)*(iy-1),(Lz/totnz)*(iz-1),f1tot(ix,iy,iz)
     $          ,f2tot(ix,iy,iz),f3tot(ix,iy,iz)
               enddo

            enddo
         enddo
         close(5)


      end if

      call MPI_Barrier(cartcomm,ierr)


      end if


      return
      end


C**********************************************************************
C Subroutine meanp:
C     Get the polarization field
C**********************************************************************

      subroutine meanp

      use Glob

      implicit none

      integer ix,iy,iz,iph,ith
 
  

C Get orientation vector
      do ith = 1,nth
         do iph = 1,nph
            do iz = 2,nz-1
               do iy = 1,2*ny
                  do ix = 1,nx

            f1(ix,iy,iz,iph,ith) = c(ix,iy,iz,iph,ith)
     $                               *sin(theta(ith))*cos(phi(iph))


                  enddo
               enddo
            enddo
         enddo
      enddo
      
      call intphi(f1,pxi)

      
      do ith = 1,nth
         do iph = 1,nph
            do iz = 2,nz-1
               do iy = 1,2*ny
                  do ix = 1,nx

            f1(ix,iy,iz,iph,ith) = c(ix,iy,iz,iph,ith)
     $                               *sin(theta(ith))*sin(phi(iph))

                  enddo
               enddo
            enddo
         enddo
      enddo
      
      
      call intphi(f1,pyi)
     


       do ith = 1,nth
         do iph = 1,nph
            do iz = 2,nz-1
               do iy = 1,2*ny
                  do ix = 1,nx

 
            f1(ix,iy,iz,iph,ith) = c(ix,iy,iz,iph,ith)
     $                               *cos(theta(ith))


                  enddo
               enddo
            enddo
         enddo
      enddo
      
   
      call intphi(f1,pzi)

C Integrate over phi to obtain stress field
      call intphi(c,ci)
      
      pxi = pxi / ci
      pyi = pyi / ci
      pzi = pzi / ci
           

      return
      end

C**********************************************************************
C Subroutine stressfield:
C     Obtain stressfield (RHS in Stokes equations)
C**********************************************************************

      subroutine stressfield
 
      use Glob

      implicit none

      integer ix,iy,iz,iph,ith
      integer ixp,ixm,iyp,iym,izp,izm
C      integer stat(MPI_STATUS_SIZE)
      real p(3),temp(3),dot


c      dy2 = dy*sqrt(1.d0+(Shear*Stime)**2.d0)

      
C Get concentration gradient using centered 2nd order FD

      do ith = 1,nth
         do iph = 1,nph

C Vector p (theta measured from 1 direction)
            p(1) = sin(theta(ith))*cos(phi(iph))
            p(2) = sin(theta(ith))*sin(phi(iph))
            p(3) = cos(theta(ith))

            do iz = 2,nz-1

               izm = iz-1
               izp = iz+1

               do iy = 1,2*ny

                  iym = iy-1
                  if (iym.lt.1)  iym = 2*ny
                  iyp = iy+1
                  if (iyp.gt.2*ny) iyp = 1

                  do ix = 1,nx
                     
                     ixm = ix-1
                     if (ixm.lt.1) ixm = nx
                     ixp = ix+1
                     if (ixp.gt.nx) ixp = 1

C Concentration gradient
               temp(1) = (c(ixp,iy,iz,iph,ith)
     $                   -c(ixm,iy,iz,iph,ith))/2.d0/dx
               temp(2) = (c(ix,iyp,iz,iph,ith)
     $                   -c(ix,iym,iz,iph,ith))/2.d0/dy2
               temp(3) = (c(ix,iy,izp,iph,ith)
     $                   -c(ix,iy,izm,iph,ith))/2.d0/dz

C Do matrix vector multiply
               dot = p(1)*temp(1)+p(2)*temp(2)+p(3)*temp(3)
               f1(ix,iy,iz,iph,ith) = dot*p(1)-temp(1)/3.d0


               enddo
             enddo
           enddo
         enddo
      enddo

      f1 = alpha * f1
      call intphi(f1,fx)



      do ith = 1,nth
         do iph = 1,nph

C Vector p (theta measured from 1 direction)
            p(1) = sin(theta(ith))*cos(phi(iph))
            p(2) = sin(theta(ith))*sin(phi(iph))
            p(3) = cos(theta(ith))

            do iz = 2,nz-1

               izm = iz-1
               izp = iz+1

               do iy = 1,2*ny

                  iym = iy-1
                  if (iym.lt.1)  iym = 2*ny
                  iyp = iy+1
                  if (iyp.gt.2*ny) iyp = 1

                  do ix = 1,nx

                     ixm = ix-1
                     if (ixm.lt.1) ixm = nx
                     ixp = ix+1
                     if (ixp.gt.nx) ixp = 1


C Concentration gradient
               temp(1) = (c(ixp,iy,iz,iph,ith)
     $                   -c(ixm,iy,iz,iph,ith))/2.d0/dx
               temp(2) = (c(ix,iyp,iz,iph,ith)
     $                   -c(ix,iym,iz,iph,ith))/2.d0/dy2
               temp(3) = (c(ix,iy,izp,iph,ith)
     $                   -c(ix,iy,izm,iph,ith))/2.d0/dz

C Do matrix vector multiply
               dot = p(1)*temp(1)+p(2)*temp(2)+p(3)*temp(3)
               f1(ix,iy,iz,iph,ith) = dot*p(2)-temp(2)/3.d0


               enddo
             enddo
           enddo
         enddo
      enddo


      f1 = alpha * f1
      call intphi(f1,fy)


      do ith = 1,nth
         do iph = 1,nph

C Vector p (theta measured from 1 direction)
            p(1) = sin(theta(ith))*cos(phi(iph))
            p(2) = sin(theta(ith))*sin(phi(iph))
            p(3) = cos(theta(ith))

            do iz = 2,nz-1

               izm = iz-1
               izp = iz+1

               do iy = 1,2*ny

                  iym = iy-1
                  if (iym.lt.1)  iym = 2*ny
                  iyp = iy+1
                  if (iyp.gt.2*ny) iyp = 1

                  do ix = 1,nx

                     ixm = ix-1
                     if (ixm.lt.1) ixm = nx
                     ixp = ix+1
                     if (ixp.gt.nx) ixp = 1


C Concentration gradient
               temp(1) = (c(ixp,iy,iz,iph,ith)
     $                   -c(ixm,iy,iz,iph,ith))/2.d0/dx
               temp(2) = (c(ix,iyp,iz,iph,ith)
     $                   -c(ix,iym,iz,iph,ith))/2.d0/dy2
               temp(3) = (c(ix,iy,izp,iph,ith)
     $                   -c(ix,iy,izm,iph,ith))/2.d0/dz

C Do matrix vector multiply
               dot = p(1)*temp(1)+p(2)*temp(2)+p(3)*temp(3)
               f1(ix,iy,iz,iph,ith) = dot*p(3)-temp(3)/3.d0


               enddo
             enddo
           enddo
         enddo
      enddo


      f1 = alpha * f1
      call intphi(f1,fz)



      return
      end


C**********************************************************************
C Subroutine velfield:
C     Obtain velocity field spectrally (Hasimoto solution)
C**********************************************************************

      subroutine velfield
      
      use Glob

      implicit none
      include 'mpif.h'



      integer ix,iy,iz,nn(3),is
      integer ixp,ixm,iyp,iym,izp,izm,i,j,k
C      integer stat(MPI_STATUS_SIZE)


      real kx,ky,kz,ksq,fku(2,3),dot(2)
      real cst
      real W(3,3),dU(3,3)
    
c      dy2 = dy*sqrt(1.d0+(Shear*Stime)**2.d0)

C*********************************************************************************
C*********************************************************************************
C*********************************************************************************

       
       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx

                fxft(((iz-2)*2*ny+(iy-1))*nx+ix) = fx(ix,iy,iz)
                fyft(((iz-2)*2*ny+(iy-1))*nx+ix) = fy(ix,iy,iz)
                fzft(((iz-2)*2*ny+(iy-1))*nx+ix) = fz(ix,iy,iz)
         
             end do
          end do
       end do
       

       work = 0
       Call fftwnd_f77_mpi(plan,1,fxft,work,0,fftw_norm_order)
       work = 0
       Call fftwnd_f77_mpi(plan,1,fyft,work,0,fftw_norm_order)
       work = 0
       Call fftwnd_f77_mpi(plan,1,fzft,work,0,fftw_norm_order)
 
       
 
       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx

             fcompx(1,ix,iy,iz) = real(fxft(((iz-2)*2*ny+(iy-1))*nx+ix))
             fcompy(1,ix,iy,iz) = real(fyft(((iz-2)*2*ny+(iy-1))*nx+ix))
             fcompz(1,ix,iy,iz) = real(fzft(((iz-2)*2*ny+(iy-1))*nx+ix))

            fcompx(2,ix,iy,iz) = aimag(fxft(((iz-2)*2*ny+(iy-1))*nx+ix))
            fcompy(2,ix,iy,iz) = aimag(fyft(((iz-2)*2*ny+(iy-1))*nx+ix))
            fcompz(2,ix,iy,iz) = aimag(fzft(((iz-2)*2*ny+(iy-1))*nx+ix))
        
             end do
          end do
       end do
 
       


C Get Fourier coefficients of the velocity field (Hasimoto)

      call MPI_CART_COORDS(cartcomm,rank,1,Scoords,ierr)

      Kmin = Scoords(1)*(nz-2) + 1
      Kmax = Scoords(1)*(nz-2) + nz-2



      do iz = Kmin,Kmax
         kz = iz-1
         if (iz.gt.(totnz/2)) kz = iz-1-totnz
         kz = kz/Lz

         do iy = 1,2*ny
            ky = iy-1
            if (iy.gt.(totny)) ky = iy-1-2*totny
            ky = ky/(2.d0*Ly)

            do ix = 1,nx
               kx = ix-1
               if (ix.gt.(totnx/2)) kx = ix-1-totnx
               kx = kx/Lx
           

              i = ix
              j = iy
              k = iz + 1 - Scoords(1)*(nz-2)          
 
              ksq = kx*kx+ky*ky+kz*kz

              fku(1,1) = fcompx(1,i,j,k)
              fku(1,2) = fcompy(1,i,j,k)
              fku(1,3) = fcompz(1,i,j,k)
              fku(2,1) = fcompx(2,i,j,k)
              fku(2,2) = fcompy(2,i,j,k)
              fku(2,3) = fcompz(2,i,j,k)
                 
              dot(1) = fku(1,1)*kx+fku(1,2)*ky+fku(1,3)*kz
              dot(2) = fku(2,1)*kx+fku(2,2)*ky+fku(2,3)*kz
                  
              ukx(1,i,j,k) = (fku(1,1)-dot(1)*kx/ksq)/ksq/(4*pi**2.00)
              ukx(2,i,j,k) = (fku(2,1)-dot(2)*kx/ksq)/ksq/(4*pi**2.00)
              uky(1,i,j,k) = (fku(1,2)-dot(1)*ky/ksq)/ksq/(4*pi**2.00)
              uky(2,i,j,k) = (fku(2,2)-dot(2)*ky/ksq)/ksq/(4*pi**2.00)
              ukz(1,i,j,k) = (fku(1,3)-dot(1)*kz/ksq)/ksq/(4*pi**2.00)
              ukz(2,i,j,k) = (fku(2,3)-dot(2)*kz/ksq)/ksq/(4*pi**2.00)
               


            enddo
         enddo
      enddo


      ix = 1
      iy = 1
      iz = 1

      if (iz.ge.Kmin.AND.iz.le.Kmax) then
 
          k = iz + 1 - Scoords(1)*(nz-2)
               
          ukx(1,ix,iy,k) = 0.00
          ukx(2,ix,iy,k) = 0.00
          uky(1,ix,iy,k) = 0.00
          uky(2,ix,iy,k) = 0.00
          ukz(1,ix,iy,k) = 0.00
          ukz(2,ix,iy,k) = 0.00

       end if


      
       Call MPI_Barrier(cartcomm,ierr)


       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx

                uxft(((iz-2)*2*ny+(iy-1))*nx+ix) =
     $          (ukx(1,ix,iy,iz)+img*ukx(2,ix,iy,iz))

                uyft(((iz-2)*2*ny+(iy-1))*nx+ix) =
     $          (uky(1,ix,iy,iz)+img*uky(2,ix,iy,iz))

                uzft(((iz-2)*2*ny+(iy-1))*nx+ix) =
     $          (ukz(1,ix,iy,iz)+img*ukz(2,ix,iy,iz))

        
             end do
          end do
       end do

       work = 0

       Call fftwnd_f77_mpi(iplan,1,uxft,work,0,fftw_norm_order)
       work = 0

       Call fftwnd_f77_mpi(iplan,1,uyft,work,0,fftw_norm_order)
       work = 0

       Call fftwnd_f77_mpi(iplan,1,uzft,work,0,fftw_norm_order)
       work = 0


       cst = real(totnx*2*totny*totnz)

       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx

              ux(ix,iy,iz) = real(uxft(((iz-2)*2*ny+(iy-1))*nx+ix))/cst 
              uy(ix,iy,iz) = real(uyft(((iz-2)*2*ny+(iy-1))*nx+ix))/cst
              uz(ix,iy,iz) = real(uzft(((iz-2)*2*ny+(iy-1))*nx+ix))/cst
       
             end do
          end do
       end do
 
c ux, uy,uz real

C COPY TO GHOST CELLS FROM NEIGHBORS


         call MPI_CART_COORDS(cartcomm,rank,1,Scoords,ierr)

C	SEND DATA TO THE RIGHT NEIGHBOR 
         coordsR(:) = Scoords(:)
         coordsR(1) = Scoords(1) + 1

         call MPI_CART_RANK(cartcomm, coordsR, rankR, ierr)

         coordsL(:) = Scoords(:)
         coordsL(1) = Scoords(1) - 1

         call MPI_CART_RANK(cartcomm, coordsL, rankL, ierr)

        
C     	RECEIVE DATA FROM THE LEFT NEIGHBOR 

         tag = 0
         Call MPI_IRECV(uxwait(:,:),nx*2*ny,MPI_REAL,rankL,
     $                          tag,cartcomm,req1,ierr)
         tag = 1
         Call MPI_IRECV(uywait(:,:),nx*2*ny,MPI_REAL,rankL,
     $                          tag,cartcomm,req2,ierr)
         tag = 2
         Call MPI_IRECV(uzwait(:,:),nx*2*ny,MPI_REAL,rankL,
     $                          tag,cartcomm,req3,ierr)



         tag = 0
         Call MPI_SEND(ux(:,:,nz-1),nx*2*ny,MPI_REAL,rankR,
     $                          tag,cartcomm,ierr)
         tag = 1
         Call MPI_SEND(uy(:,:,nz-1),nx*2*ny,MPI_REAL,rankR,
     $                          tag,cartcomm,ierr)
         tag = 2
         Call MPI_SEND(uz(:,:,nz-1),nx*2*ny,MPI_REAL,rankR,
     $                          tag,cartcomm,ierr)

 
         call MPI_WAIT(req1,stat,ierr)
         call MPI_WAIT(req2,stat,ierr)
         call MPI_WAIT(req3,stat,ierr)

         ux(:,:,1) = uxwait(:,:)
         uy(:,:,1) = uywait(:,:)
         uz(:,:,1) = uzwait(:,:)


C	RECEIVE DATA FROM THE RIGHT NEIGHBOR

         tag = 3
         Call MPI_IRECV(uxwait(:,:),nx*2*ny,MPI_REAL,rankR,
     $                          tag,cartcomm,req1,ierr)
         tag = 4
         Call MPI_IRECV(uywait(:,:),nx*2*ny,MPI_REAL,rankR,
     $                          tag,cartcomm,req2,ierr)
         tag = 5
         Call MPI_IRECV(uzwait(:,:),nx*2*ny,MPI_REAL,rankR,
     $                          tag,cartcomm,req3,ierr)



C	SEND DATA TO THE LEFT NEIGHBOR
         tag = 3
         Call MPI_SEND(ux(:,:,2),nx*2*ny,MPI_REAL,rankL,
     $                          tag,cartcomm,ierr)
         tag = 4
         Call MPI_SEND(uy(:,:,2),nx*2*ny,MPI_REAL,rankL,
     $                          tag,cartcomm,ierr)
         tag = 5
         Call MPI_SEND(uz(:,:,2),nx*2*ny,MPI_REAL,rankL,
     $                          tag,cartcomm,ierr)


         call MPI_WAIT(req1,stat,ierr)
         call MPI_WAIT(req2,stat,ierr)
         call MPI_WAIT(req3,stat,ierr)

         ux(:,:,nz) = uxwait(:,:)
         uy(:,:,nz) = uywait(:,:)
         uz(:,:,nz) = uzwait(:,:)



C	RECEIVING DATA FROM NEIGHBORS AND STORING IT IN THE GHOST CELLS


C Get tensor bEW for angle dynamics

      do iz = 2,nz-1

         izm = iz-1
         izp = iz+1

         do iy = 1,2*ny

            iym = iy-1
            if (iym.lt.1)  iym = 2*ny
            iyp = iy+1
            if (iyp.gt.2*ny) iyp = 1

            do ix = 1,nx
               ixm = ix-1
               if (ixm.lt.1)  ixm = nx
               ixp = ix+1
               if (ixp.gt.nx) ixp = 1

                  dU(1,1) = (ux(ixp,iy,iz)-ux(ixm,iy,iz))/2.00/dx

                  dU(1,2) = (ux(ix,iyp,iz)-ux(ix,iym,iz))/2.00/dy2

                  dU(1,3) = (ux(ix,iy,izp)-ux(ix,iy,izm))/2.00/dz


                  dU(2,1) = (uy(ixp,iy,iz)-uy(ixm,iy,iz))/2.00/dx

                  dU(2,2) = (uy(ix,iyp,iz)-uy(ix,iym,iz))/2.00/dy2

                  dU(2,3) = (uy(ix,iy,izp)-uy(ix,iy,izm))/2.00/dz


                  dU(3,1) = (uz(ixp,iy,iz)-uz(ixm,iy,iz))/2.00/dx

                  dU(3,2) = (uz(ix,iyp,iz)-uz(ix,iym,iz))/2.00/dy2

                  dU(3,3) = (uz(ix,iy,izp)-uz(ix,iy,izm))/2.00/dz

 
            do j = 1,3
            do i = 1,3
               E(i,j,ix,iy,iz) = (dU(i,j)+dU(j,i))/2.00
               W(i,j) = (dU(i,j)-dU(j,i))/2.00
               bEW(i,j,ix,iy,iz) = gamma*E(i,j,ix,iy,iz)+W(i,j)
            enddo
            enddo
            
            enddo
         enddo
      enddo



      return
      end



C**********************************************************************
C Subroutine fluxes:
C     Calculate fluxes in evolution equation for the configuration
C     field
C**********************************************************************

      subroutine flux
 
      use Glob

      implicit none
      include 'mpif.h'

      real skxn,skyn,skzn,sknorm,uskn,Kernelip

      integer ix,iy,iz,iph,ith,iphth,ithpp,iphpp
      integer ixm,iym,izm,iphm,ithm,ixp,iyp,izp,iphp,ithp
      integer i,j,k,is
C      integer stat(MPI_STATUS_SIZE)


      real Vx,Vy,Vz,Vxm,Vym,Vzm,Vxp,Vyp,Vzp,Tumble(nph,nth),pp(3)
      real delta(3,3),pdot(3),p(3),Vymx,Vypx!,Tumblei,Kernel,Tumbleip
      real Tumblei(nx,ny+1,nz),Tumbleip(nx,ny+1,nz),Kernel,Kerneli
c      dy2 = dy*sqrt(1.d0+(Shear*Stime)**2.d0)

      do j = 1,3
         do i = 1,3
            delta(i,j) = 0.00
         enddo
         delta(j,j) = 1.00
      enddo

C Initialize flux

      F = 0.d0

C Calculate convective flux in x-y-z
      do ith = 1,nth
         do iph = 1,nph
            do iz = 2,nz-1

               izm = iz-1
               izp = iz+1

               do iy = 1,ny+1

                  iym = iy-1
                  if (iym.lt.1)  iym = 2*ny
                  iyp = iy+1
c                  if (iyp.gt.2*ny) iyp = 1

                  do ix = 1,nx

                     ixm = ix-1
                     if (ixm.lt.1)  ixm = nx
                     ixp = ix+1
                     if (ixp.gt.nx) ixp = 1

                     Vxm = sin(theta(ith))*cos(phi(iph)) +
     $                     ux(ixm,iy,iz)
                     Vxp = sin(theta(ith))*cos(phi(iph)) +
     $                     ux(ixp,iy,iz)
                     F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) +
     $               (Vxp*c(ixp,iy,iz,iph,ith)-
     $               Vxm*c(ixm,iy,iz,iph,ith))/2.d0/dx


                     Vym = sin(theta(ith))*sin(phi(iph)) +
     $                     uy(ix,iym,iz)
                     Vyp = sin(theta(ith))*sin(phi(iph)) +
     $                     uy(ix,iyp,iz)
                     F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) +
     $               (Vyp*c(ix,iyp,iz,iph,ith)-
     $               Vym*c(ix,iym,iz,iph,ith))/2.d0/dy2 


                     Vzm = cos(theta(ith)) + uz(ix,iy,izm)
                     Vzp = cos(theta(ith)) + uz(ix,iy,izp)
                     F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) +
     $               (Vzp*c(ix,iy,izp,iph,ith)-
     $               Vzm*c(ix,iy,izm,iph,ith))/2.d0/dz

                enddo
              enddo
            enddo
         enddo
      enddo

C Calculate convective flux in theta
      do ith = 1,nth
         do iph = 1,nph
            p(1) = sin(theta(ith))*cos(phi(iph))
            p(2) = sin(theta(ith))*sin(phi(iph))
            p(3) = cos(theta(ith))

            do iz = 2,nz-1
               do iy = 1,ny+1
                  do ix = 1,nx

C Calculate thetadot
               do i = 1,3
               pdot(i) = 0.00
               do k = 1,3
               do j = 1,3
                  pdot(i) = pdot(i) + (delta(i,j)-p(i)*p(j))
     $                 *bEW(j,k,ix,iy,iz)*p(k)
               enddo
               enddo
               enddo


C	f1 = pthetadot, f2 = pphidot

               f1(ix,iy,iz,iph,ith) =
     $         cos(theta(ith))*cos(phi(iph))*pdot(1) +
     $         cos(theta(ith))*sin(phi(iph))*pdot(2) -
     $         sin(theta(ith))*pdot(3)


            enddo
            enddo 
            enddo
         enddo
      enddo



      ith = 1
   
      ithp = ith+1
      do iph = 1,nph

         iphm = iph-1
         if (iphm.lt.1) iphm = nph
         iphp = iph+1
         if (iphp.gt.nph) iphp = 1

         iphth = iph + nph/2
         if (iphth.gt.nph) iphth = iphth - nph

         do iz = 2,nz-1
            do iy =1,ny+1
               do ix =1,nx


               F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) +
     $        (f1(ix,iy,iz,iph,ithp)*c(ix,iy,iz,iph,ithp)*
     $         sin(theta(ithp))
     $        +f1(ix,iy,iz,iphth,ith)*c(ix,iy,iz,iphth,ith)*
     $         sin(theta(ith))) / 2.d0/dth/Sin(theta(ith))
 

 
              end do
           end do
        end do
      end do 


      do ith = 2,nth-1

         ithm = ith-1
         ithp = ith+1
         do iph = 1,nph
            iphm = iph-1
            if (iphm.lt.1) iphm = nph
            iphp = iph+1
            if (iphp.gt.nph) iphp = 1
            do iz = 2,nz-1
               do iy =1,ny+1
                  do ix =1,nx


             F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) +
     $        ((f1(ix,iy,iz,iph,ithp)*c(ix,iy,iz,iph,ithp)
     *         *sin(theta(ithp))
     $        -f1(ix,iy,iz,iph,ithm)*c(ix,iy,iz,iph,ithm)
     $         *sin(theta(ithm)))/sin(theta(ith))/2.d0/dth)
              
                 enddo 
               enddo               
            enddo
         enddo
      enddo

 
      ith = nth
      
      ithm = ith-1
      do iph = 1,nph
         iphm = iph-1
         if (iphm.lt.1) iphm = nph
         iphp = iph+1
         if (iphp.gt.nph) iphp = 1

         iphth = iph + nph/2
         if (iphth.gt.nph) iphth = iphth - nph

         do iz = 2,nz-1
            do iy =1,ny+1
               do ix =1,nx

               F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) +
     $        (-f1(ix,iy,iz,iphth,ith)*c(ix,iy,iz,iphth,ith)*
     $         sin(theta(ith))
     $        -f1(ix,iy,iz,iph,ithm)*c(ix,iy,iz,iph,ithm)*
     $         sin(theta(ithm)))/2.d0/dth/sin(theta(ith))

               end do
            end do
          end do
       end do



      do ith = 1,nth
         do iph = 1,nph
            p(1) = sin(theta(ith))*cos(phi(iph))
            p(2) = sin(theta(ith))*sin(phi(iph))
            p(3) = cos(theta(ith))

            do iz = 2,nz-1
               do iy = 1,ny+1
                  do ix = 1,nx

C Calculate thetadot
               do i = 1,3
               pdot(i) = 0.00
               do k = 1,3
               do j = 1,3
                  pdot(i) = pdot(i) + (delta(i,j)-p(i)*p(j))
     $                 *bEW(j,k,ix,iy,iz)*p(k)
               enddo
               enddo
               enddo

C	f1 = pphidot

               f1(ix,iy,iz,iph,ith) =
     $         -sin(phi(iph))*pdot(1)+cos(phi(iph))*pdot(2)

!               if (sin(theta).ne.0.00) then
!                  thdot(ix,iy,ith) = -pdot(1)/sin(theta)
!               else
!                  thdot(ix,iy,ith) = pdot(2)/cos(theta)
!               endif

            enddo
            enddo 
            enddo
         enddo
      enddo


      do ith = 1,nth

         ithm = ith-1
         ithp = ith+1

         do iph = 1,nph
            iphm = iph-1
            if (iphm.lt.1) iphm = nph
            iphp = iph+1
            if (iphp.gt.nph) iphp = 1

            p(1) = sin(theta(ith))*cos(phi(iph))
            p(2) = sin(theta(ith))*sin(phi(iph))
            p(3) = cos(theta(ith))
 
            do iz = 2,nz-1
               do iy =1,ny+1
                  do ix =1,nx


             F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) +
     $        (f1(ix,iy,iz,iphp,ith)*c(ix,iy,iz,iphp,ith)
     $        -f1(ix,iy,iz,iphm,ith)*c(ix,iy,iz,iphm,ith))
     $        /2.d0/dph/sin(theta(ith))


C       INCLUDE RUN AND TUMBLE EFFECT

             skxn = skx(1,ix,iy,iz)!/sknorm
             skyn = sky(1,ix,iy,iz)!/sknorm
             skzn = skz(1,ix,iy,iz)!/sknorm
	     

C       BACTERIA THAT STOPS MOVING DUE TO TUMBLING


	     F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) +
     $       c(ix,iy,iz,iph,ith)*l0*exp(-1.d0*xi*((skc(1,ix,iy,iz)
     $       -skcb(1,ix,iy,iz))/dt +(p(1)+ux(ix,iy,iz))*skxn
     $       +(p(2)+uy(ix,iy,iz))*skyn+(p(3)+uz(ix,iy,iz))
     $       *skzn))
              
                 enddo 
               enddo               
            enddo
	enddo
      enddo

   

C       INCLUDE RUN AND TUMBLE EFFECT
C       BACTERIA THAT STARTS MOVING HAVING TUMBLED

      do ith = 1,nth
         do iph = 1,nph
            p(1) = sin(theta(ith))*cos(phi(iph))
            p(2) = sin(theta(ith))*sin(phi(iph))
            p(3) = cos(theta(ith))
            Tumblei = 0.d0
	    Kerneli = 0.d0

            Do ithpp = 1,nth
               Tumbleip = 0.d0
	       Kernelip = 0.d0

               Do iphpp = 1,nph
                  pp(1) = sin(theta(ithpp))*cos(phi(iphpp))
                  pp(2) = sin(theta(ithpp))*sin(phi(iphpp))
                  pp(3) = cos(theta(ithpp))

                  Kernel = 1.d0/(4.d0*pi)

                  Kernelip = Kernelip + Kernel*dph
                  Do iz = 2,nz-1
                     Do iy = 1,ny+1
                        Do ix = 1,nx

                           skxn = skx(1,ix,iy,iz)!/sknorm
                           skyn = sky(1,ix,iy,iz)!/sknorm
                           skzn = skz(1,ix,iy,iz)!/sknorm

			   Tumble(iphpp,ithpp) = Kernel*
     $                    c(ix,iy,iz,iphpp,ithpp)*l0*exp(-1.d0*xi*
     $                    ((skc(1,ix,iy,iz)-skcb(1,ix,iy,iz))/dt+
     $                    (pp(1)+ux(ix,iy,iz))*skxn+
     $                    (pp(2)+uy(ix,iy,iz))*skyn+
     $                    (pp(3)+uz(ix,iy,iz))*skzn))


			  
                           Tumbleip(ix,iy,iz) = Tumbleip(ix,iy,iz) + 
     $                                 Tumble(iphpp,ithpp)*dph

                        End Do
                     End Do
                  End Do
 	
               End Do
               Tumblei = Tumblei + Tumbleip*sin(theta(ithpp))*dth
	       Kerneli = Kerneli + Kernelip*sin(theta(ithpp))*dth
            End Do
          F(:,:,:,iph,ith) = F(:,:,:,iph,ith) - Tumblei(:,:,:)/Kerneli

          End Do
      End Do



c               sknorm = sqrt(skx(1,ix,iy,iz)**2.d0+sky(1,ix,iy,iz)**2.d0
c     $                      +skz(1,ix,iy,iz)**2.d0)


c               if (sknorm.ne.0.d0) then

c                skxn = skx(1,ix,iy,iz)!/sknorm
c                skyn = sky(1,ix,iy,iz)!/sknorm
c                skzn = skz(1,ix,iy,iz)!/sknorm

c               else
c                  skxn = 0.d0
c                  skyn = 0.d0
c                  skzn = 0.d0
c               end if

C       BACTERIA THAT STOPS MOVING DUE TO TUMBLING

c               F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) -
c     $         c(ix,iy,iz,iph,ith)*(p(1)*skxn+p(2)*skyn+p(3)*skzn)/tau

C       BACTERIA THAT STARTS MOVING HAVING TUMBLED

c               Tumblei = 0.d0

c               Do ithpp = 1,nth
c                  Tumbleip = 0.d0

c                  Do iphpp = 1,nph
c                  pp(1) = sin(theta(ithpp))*cos(phi(iphpp))
c                  pp(2) = sin(theta(ithpp))*sin(phi(iphpp))
c                  pp(3) = cos(theta(ithpp))

c                  Kernel = (deltai/(4.d0*pi*sinh(deltai)))*
c     $            exp(deltai*(p(1)*pp(1)+p(2)*pp(2)+p(3)*pp(3)))

c                  Tumble(iphpp,ithpp) = Kernel*c(ix,iy,iz,iphpp,ithpp)*
c     $                          (pp(1)*skxn+pp(2)*skyn+pp(3)*skzn)/tau

c                  Tumbleip = Tumbleip + Tumble(iphpp,ithpp)*dph

c                  End Do
c                  Tumblei = Tumblei + Tumbleip*sin(theta(ithpp))*dth

c               End Do

c               F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) + Tumblei
               
c                End Do
c              End Do
c           enddo
c         enddo
c      enddo





C Calculate diffusive flux in x-y, and diffusive flux in theta

      ith = 1
      
      ithp = ith+1
      do iph = 1,nph
         iphm = iph-1
         if (iphm.lt.1) iphm = nph
         iphp = iph+1
         if (iphp.gt.nph) iphp = 1

         iphth = iph + nph/2
         if (iphth.gt.nph) iphth = iphth - nph

            do iz = 2,nz-1

               izm = iz-1
               izp = iz+1

               do iy = 1,ny+1

                  iym = iy-1
                  if (iym.lt.1)  iym = 2*ny
                  iyp = iy+1
c                  if (iyp.gt.2*ny) iyp = 1

                  do ix = 1,nx

                     ixm = ix-1
                     if (ixm.lt.1)  ixm = nx
                     ixp = ix+1
                     if (ixp.gt.nx) ixp = 1


               F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) - dfx*(
     $         (c(ixp,iy,iz,iph,ith)-2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ixm,iy,iz,iph,ith))/dx**2.00
     $         +(c(ix,iyp,iz,iph,ith)-2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ix,iym,iz,iph,ith))/dy2**2.00
     $         +(c(ix,iy,izp,iph,ith)-2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ix,iy,izm,iph,ith))/dz**2.00)



               F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) - dfr*(
     $         (1.d0/((sin(theta(ith)))**2.d0))*
     $         (c(ix,iy,iz,iphp,ith)
     $         -2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ix,iy,iz,iphm,ith))/(dph*dph)
     $         +
     $          (sin(theta(ith)+pi/2.d0/nth)*(c(ix,iy,iz,iph,ithp)-
     $         c(ix,iy,iz,iph,ith))
     $         - sin(theta(ith)-pi/2.d0/nth)*(c(ix,iy,iz,iph,ith)-
     $         c(ix,iy,iz,iphth,ith)))
     $         /(dth*dth*sin(theta(ith))))




              end do
           end do
        end do
      end do 


      do ith = 2,nth-1

         ithm = ith-1
         ithp = ith+1
         do iph = 1,nph
            iphm = iph-1
            if (iphm.lt.1) iphm = nph
            iphp = iph+1
            if (iphp.gt.nph) iphp = 1

            do iz = 2,nz-1

               izm = iz-1
               izp = iz+1

               do iy = 1,ny+1

                  iym = iy-1
                  if (iym.lt.1)  iym = 2*ny
                  iyp = iy+1
c                  if (iyp.gt.2*ny) iyp = 1

                  do ix = 1,nx

                     ixm = ix-1
                     if (ixm.lt.1)  ixm = nx
                     ixp = ix+1
                     if (ixp.gt.nx) ixp = 1



               F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) - dfx*(
     $         (c(ixp,iy,iz,iph,ith)-2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ixm,iy,iz,iph,ith))/dx**2.00
     $         +(c(ix,iyp,iz,iph,ith)-2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ix,iym,iz,iph,ith))/dy2**2.00
     $         +(c(ix,iy,izp,iph,ith)-2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ix,iy,izm,iph,ith))/dz**2.00)



              F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) - dfr*(
     $         (1.d0/(sin(theta(ith))**2.d0))*
     $         (c(ix,iy,iz,iphp,ith)
     $         -2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ix,iy,iz,iphm,ith))/(dph*dph)
     $         +
     $          (sin(theta(ith)+pi/2.d0/nth)*(c(ix,iy,iz,iph,ithp)-
     $         c(ix,iy,iz,iph,ith))
     $         - sin(theta(ith)-pi/2.d0/nth)*(c(ix,iy,iz,iph,ith)-
     $         c(ix,iy,iz,iph,ithm)))
     $         /(dth*dth*sin(theta(ith))))



               
                 enddo 
               enddo               
            enddo
         enddo
      enddo


      ith = nth
      
      ithm = ith-1
      do iph = 1,nph
         iphm = iph-1
         if (iphm.lt.1) iphm = nph
         iphp = iph+1
         if (iphp.gt.nph) iphp = 1

         iphth = iph + nph/2
         if (iphth.gt.nph) iphth = iphth - nph

            do iz = 2,nz-1

               izm = iz-1
               izp = iz+1

               do iy = 1,ny+1

                  iym = iy-1
                  if (iym.lt.1)  iym = 2*ny
                  iyp = iy+1
c                  if (iyp.gt.2*ny) iyp = 1

                  do ix = 1,nx

                     ixm = ix-1
                     if (ixm.lt.1)  ixm = nx
                     ixp = ix+1
                     if (ixp.gt.nx) ixp = 1




               F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) - dfx*(
     $         (c(ixp,iy,iz,iph,ith)-2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ixm,iy,iz,iph,ith))/dx**2.00
     $         +(c(ix,iyp,iz,iph,ith)-2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ix,iym,iz,iph,ith))/dy2**2.00
     $         +(c(ix,iy,izp,iph,ith)-2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ix,iy,izm,iph,ith))/dz**2.00)



              F(ix,iy,iz,iph,ith) = F(ix,iy,iz,iph,ith) - dfr*(
     $         (1.d0/(sin(theta(ith))**2.d0))*
     $         (c(ix,iy,iz,iphp,ith)
     $         -2.d0*c(ix,iy,iz,iph,ith)+
     $         c(ix,iy,iz,iphm,ith))/(dph*dph)
     $         +
     $          (sin(theta(ith)+pi/2.d0/nth)*(c(ix,iy,iz,iphth,ith)-
     $         c(ix,iy,iz,iph,ith))
     $         - sin(theta(ith)-pi/2.d0/nth)*(c(ix,iy,iz,iph,ith)-
     $         c(ix,iy,iz,iph,ithm)))
     $         /(dth*dth*sin(theta(ith))))

               end do
            end do
          end do
       end do

      return
      end


C**********************************************************************
C Subroutine euler:
C     Perform initial step using Euler algorithm
C**********************************************************************

      subroutine euler

      use Glob

      implicit none

      integer ix,iy,iz,iph,ith
  
c      dy2 = dy*sqrt(1.d0+(Shear*Stime)**2.d0)

      
C Get stress field for Stokes equations 
      call stressfield
      call MPI_Barrier(cartcomm,ierr)
C Get velocity field spectrally
      call velfield
      call MPI_Barrier(cartcomm,ierr)
C
C Get fluxes
 
      call flux
      call MPI_Barrier(cartcomm,ierr)

      Do iy = 1,ny+1
     
         c(:,iy,:,:,:) = c(:,iy,:,:,:) - dt * F(:,iy,:,:,:)

      End Do

      Fo = F

      Call reflect

C Update time
      time = time + dt

      return
      end


C**********************************************************************
C Subroutine update:
C     Advance configuration field using 2nd order Adams-Bashforth
C**********************************************************************

      subroutine update
 
      use Glob

      implicit none
      include 'mpif.h'



      integer ix,iy,iz,iph,ith
      integer i,j,k,is,iremesh
      real imp,imp1
  
c      dy2 = dy*sqrt(1.d0+(Shear*Stime)**2.d0)



C COPY TO GHOST CELLS FROM NEIGHBORS


C******************************************************************
C	        COORDINATES AND RANKS OF NEIGHBORS
 
         call MPI_CART_COORDS(cartcomm,rank,1,Scoords,ierr)

         coordsR(:) = Scoords(:)
         coordsR(1) = Scoords(1) + 1

         call MPI_CART_RANK(cartcomm, coordsR, rankR, ierr)


         coordsL(:) = Scoords(:)
         coordsL(1) = Scoords(1) - 1

         call MPI_CART_RANK(cartcomm, coordsL, rankL, ierr)


C********************************************************************

         tag = 0
         Call MPI_IRECV(cwait(:,:,:,:),nx*2*ny*nph*nth,MPI_REAL,rankL,
     $                          tag,cartcomm,req,ierr)
 
         csend(:,:,:,:) = c(:,:,nz-1,:,:)

         tag = 0
         Call MPI_SEND(csend(:,:,:,:),nx*2*ny*nph*nth,MPI_REAL,rankR,
     $                          tag,cartcomm,ierr)

         call MPI_WAIT(req,stat,ierr)

         c(:,:,1,:,:) = cwait(:,:,:,:)

C     	RECEIVE DATA FROM THE LEFT NEIGHBOR 

C*********************************************************************


C********************************************************************

         tag = 1
         Call MPI_IRECV(cwait(:,:,:,:),nx*2*ny*nph*nth,MPI_REAL,rankR,
     $                          tag,cartcomm,req,ierr)
 
 
         csend(:,:,:,:) = c(:,:,2,:,:)


         tag = 1
         Call MPI_SEND(csend(:,:,:,:),nx*2*ny*nph*nth,MPI_REAL,rankL,
     $                          tag,cartcomm,ierr)

         call MPI_WAIT(req,stat,ierr)

         c(:,:,nz,:,:) = cwait(:,:,:,:)

C**********************************************************************
C********************************************************************
C***********************************************************************

      
C Get stress field for Stokes equations 

      call stressfield
      Call MPI_Barrier(cartcomm,ierr)

C Get velocity field spectrally

      call velfield
      Call MPI_Barrier(cartcomm,ierr)

C Get fluxes

      call flux
      Call MPI_Barrier(cartcomm,ierr)

      Do iy = 1,ny+1

         c(:,iy,:,:,:)  = c(:,iy,:,:,:) - 
     $                    0.5d0*dt*(3.d0*F(:,iy,:,:,:) - Fo(:,iy,:,:,:))

      End Do

      Fo = F

      Call reflect

C Update time
      time = time + dt


      return
      end

**********************************************************************
C Subroutine eulerscalar:
C     Perform initial step using Euler algorithm for scalar field
C**********************************************************************

      subroutine eulerscalar

      use glob

      implicit none


      integer i,j,k,ix,iy,iz
      real cst,kx,ky,kz,k2,ss

      Call BC


C Get fluxes
C First, get velocity in real space
C Then, get scalar density gradient in real space
      call MPI_CART_COORDS(cartcomm,rank,1,Scoords,ierr)

      Kmin = Scoords(1)*(nz-2) + 1
      Kmax = Scoords(1)*(nz-2) + nz-2

      do iz = Kmin,Kmax
         kz = 2*pi*(iz-1)
         if (iz.gt.(totnz/2)) kz = 2*pi*(iz-1-totnz)
         kz = kz/Lz

         do iy = 1,2*ny
            ky = 2*pi*(iy-1)
            if (iy.gt.(totny)) ky = 2*pi*(iy-1-2*totny)
            ky = ky/(2.d0*Ly)

            do ix = 1,nx
               kx = 2*pi*(ix-1)
               if (ix.gt.(totnx/2)) kx = 2*pi*(ix-1-totnx)
               kx = kx/Lx
           
               i = ix
               j = iy
               k = iz + 1 - Scoords(1)*(nz-2)          

               skx(1,i,j,k) = -kx*sk(2,i,j,k)
               skx(2,i,j,k) =  kx*sk(1,i,j,k)
               sky(1,i,j,k) = -ky*sk(2,i,j,k)
               sky(2,i,j,k) =  ky*sk(1,i,j,k)
               skz(1,i,j,k) = -kz*sk(2,i,j,k)
               skz(2,i,j,k) =  kz*sk(1,i,j,k)

            enddo
         enddo
      enddo

       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx

                skxft(((iz-2)*2*ny+(iy-1))*nx+ix) =
     $          (skx(1,ix,iy,iz)+img*skx(2,ix,iy,iz))

                skyft(((iz-2)*2*ny+(iy-1))*nx+ix) =
     $          (sky(1,ix,iy,iz)+img*sky(2,ix,iy,iz))

                skzft(((iz-2)*2*ny+(iy-1))*nx+ix) =
     $          (skz(1,ix,iy,iz)+img*skz(2,ix,iy,iz))
      
             end do
          end do
       end do

       work = 0
       Call fftwnd_f77_mpi(iplan,1,skxft,work,0,fftw_norm_order)
       work = 0
       Call fftwnd_f77_mpi(iplan,1,skyft,work,0,fftw_norm_order)
       work = 0
       Call fftwnd_f77_mpi(iplan,1,skzft,work,0,fftw_norm_order)


       cst = real(totnx*2*totny*totnz)

       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx

               skx(1,ix,iy,iz) = real(skxft(((iz-2)*2*ny+(iy-1))*nx+ix))
               sky(1,ix,iy,iz) = real(skyft(((iz-2)*2*ny+(iy-1))*nx+ix))
               skz(1,ix,iy,iz) = real(skzft(((iz-2)*2*ny+(iy-1))*nx+ix))
    
             end do
          end do
       end do

C Take dot product of gradient with velocity in real space
      do iz = 2,nz-1
         do iy = 1,2*ny
            do ix = 1,nx

               fk(1,ix,iy,iz) = ux(ix,iy,iz)*skx(1,ix,iy,iz)
     $                         +uy(ix,iy,iz)*sky(1,ix,iy,iz)
     $                         +uz(ix,iy,iz)*skz(1,ix,iy,iz)

               fk(2,ix,iy,iz) = 0.d0

            enddo
         enddo
      enddo
C Get flux in Fourier by taking transform

       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx

                skft(((iz-2)*2*ny+(iy-1))*nx+ix) = fk(1,ix,iy,iz)
         
             end do
          end do
       end do
       

       work = 0
       Call fftwnd_f77_mpi(plan,1,skft,work,0,fftw_norm_order)
      
 
       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx

              fk(1,ix,iy,iz) = 
     $        real(skft(((iz-2)*2*ny+(iy-1))*nx+ix))/cst 

              fk(2,ix,iy,iz) = 
     $        aimag(skft(((iz-2)*2*ny+(iy-1))*nx+ix))/cst
        
             end do
          end do
       end do
 

C CALCULATING THE SINK TERM IN THE OXYGEN TRANSPORT EQ
      cs        = 0.d0

      cs(1,:,:,:) = ci(:,:,:)*((skc(1,:,:,:)-0.1)/0.9)**(0.2)



       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx

                skft(((iz-2)*2*ny+(iy-1))*nx+ix) = cs(1,ix,iy,iz)
     
		if (skc(1,ix,iy,iz).le.0.1) then
		skft(((iz-2)*2*ny+(iy-1))*nx+ix) = 0.d0
		endif

             end do
          end do
       end do
       

       work = 0
       Call fftwnd_f77_mpi(plan,1,skft,work,0,fftw_norm_order)
      
 
       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx

              cs(1,ix,iy,iz) = 
     $        real(skft(((iz-2)*2*ny+(iy-1))*nx+ix))/cst 

              cs(2,ix,iy,iz) = 
     $        aimag(skft(((iz-2)*2*ny+(iy-1))*nx+ix))/cst
        
             end do
          end do
       end do

C Update scalar field using explicit Euler scheme
C and copy flux      
c      D = 0.0003
      do iz = Kmin,Kmax
         kz = 2*pi*(iz-1)
         if (iz.gt.(totnz/2)) kz = 2*pi*(iz-1-totnz)
         kz = kz/Lz

         do iy = 1,2*ny
            ky = 2*pi*(iy-1)
            if (iy.gt.(totny)) ky = 2*pi*(iy-1-2*totny)
            ky = ky/(2.d0*Ly)

            do ix = 1,nx
               kx = 2*pi*(ix-1)
               if (ix.gt.(totnx/2)) kx = 2*pi*(ix-1-totnx)
               kx = kx/Lx
           
               i = ix
               j = iy
               k = iz + 1 - Scoords(1)*(nz-2)          

               k2 = kx*kx+ky*ky+kz*kz

               sk(1,i,j,k) = (sk(1,i,j,k) - 
     $            dt*(fk(1,i,j,k)+ocr*cs(1,i,j,k))) / (1+dt*D*k2)
               sk(2,i,j,k) = (sk(2,i,j,k) - 
     $            dt*(fk(2,i,j,k)+ocr*cs(2,i,j,k))) / (1+dt*D*k2)

               fko(:,i,j,k) = fk(:,i,j,k)
               cso(:,i,j,k) = cs(:,i,j,k)

            enddo
         enddo
      enddo

C Calculate mix norm using fourier representation
c      snorm = 0.d0
c      do ix = 1,nsx
c         if (ix.lt.nsx/2) then
c            kx = 2*pi*(ix-1)/Lx
c         else
c            kx = -2*pi*(nsx-ix+1)/Lx
c         endif
c         do iy = 1,2*nsy
c            if (iy.lt.nsy) then
c               ky = 2*pi*(iy-1)/(2*Ly)
c            else
c              ky = -2*pi*(2*nsy-iy+1)/(2*Ly)
c           endif
c
c            k2 = kx**2.d0+ky**2.d0
c            ss = sk(1,ix,iy)**2.d0+sk(2,ix,iy)**2.d0
c            snorm = snorm + ss/sqrt(1+k2)
c         enddo
c      enddo
c      snorm = sqrt(snorm)
c      write(10,*) time,snorm

      return
      End Subroutine

C**********************************************************************
C Subroutine updatescalar:
C     Advance scalar field
C**********************************************************************

      subroutine updatescalar

      use glob

      implicit none

      integer i,j,k,ix,iy,iz
      real cst,kx,ky,kz,k2,ss


C Calculate mix norm using fourier representation
c      snorm = 0.d0
c      do iy = 1,2*nsy
c         if (iy.lt.nsy) then
c            ky = 2*pi*(iy-1)/(2*Ly)
c         else
c            ky = -2*pi*(2*nsy-iy+1)/(2*Ly)
c         endif
c
c         do ix = 1,nsx
c           if (ix.lt.nsx/2) then
c               kx = 2*pi*(ix-1)/Lx
c            else
c               kx = -2*pi*(nsx-ix+1)/Lx
c            endif
c
c            k2 = kx**2.d0+ky**2.d0
c
c            ss = sk(1,ix,iy)**2.d0+sk(2,ix,iy)**2.d0
c            snorm = snorm + ss/sqrt(1+k2)
c         enddo
c      enddo
c      snorm = sqrt(snorm)
c      write(10,*) time,snorm



      Call BC

C Then, get scalar density gradient in real space
      call MPI_CART_COORDS(cartcomm,rank,1,Scoords,ierr)

      Kmin = Scoords(1)*(nz-2) + 1
      Kmax = Scoords(1)*(nz-2) + nz-2

      do iz = Kmin,Kmax
         kz = 2*pi*(iz-1)
         if (iz.gt.(totnz/2)) kz = 2*pi*(iz-1-totnz)
         kz = kz/Lz

         do iy = 1,2*ny
            ky = 2*pi*(iy-1)
            if (iy.gt.(totny)) ky = 2*pi*(iy-1-2*totny)
            ky = ky/(2.d0*Ly)

            do ix = 1,nx
               kx = 2*pi*(ix-1)
               if (ix.gt.(totnx/2)) kx = 2*pi*(ix-1-totnx)
               kx = kx/Lx
           
               i = ix
               j = iy
               k = iz + 1 - Scoords(1)*(nz-2)          

               skx(1,i,j,k) = -kx*sk(2,i,j,k)
               skx(2,i,j,k) =  kx*sk(1,i,j,k)
               sky(1,i,j,k) = -ky*sk(2,i,j,k)
               sky(2,i,j,k) =  ky*sk(1,i,j,k)
               skz(1,i,j,k) = -kz*sk(2,i,j,k)
               skz(2,i,j,k) =  kz*sk(1,i,j,k)

            enddo
         enddo
      enddo

       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx

                skxft(((iz-2)*2*ny+(iy-1))*nx+ix) =
     $          (skx(1,ix,iy,iz)+img*skx(2,ix,iy,iz))

                skyft(((iz-2)*2*ny+(iy-1))*nx+ix) =
     $          (sky(1,ix,iy,iz)+img*sky(2,ix,iy,iz))

                skzft(((iz-2)*2*ny+(iy-1))*nx+ix) =
     $          (skz(1,ix,iy,iz)+img*skz(2,ix,iy,iz))
      
             end do
          end do
       end do

       work = 0
       Call fftwnd_f77_mpi(iplan,1,skxft,work,0,fftw_norm_order)
       work = 0
       Call fftwnd_f77_mpi(iplan,1,skyft,work,0,fftw_norm_order)
       work = 0
       Call fftwnd_f77_mpi(iplan,1,skzft,work,0,fftw_norm_order)


       cst = real(totnx*2*totny*totnz)

       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx

               skx(1,ix,iy,iz) = real(skxft(((iz-2)*2*ny+(iy-1))*nx+ix))
               sky(1,ix,iy,iz) = real(skyft(((iz-2)*2*ny+(iy-1))*nx+ix))
               skz(1,ix,iy,iz) = real(skzft(((iz-2)*2*ny+(iy-1))*nx+ix))
    
             end do
          end do
       end do

C Take dot product of gradient with velocity in real space
      do iz = 2,nz-1
         do iy = 1,2*ny
            do ix = 1,nx

               fk(1,ix,iy,iz) = ux(ix,iy,iz)*skx(1,ix,iy,iz)
     $                         +uy(ix,iy,iz)*sky(1,ix,iy,iz)
     $                         +uz(ix,iy,iz)*skz(1,ix,iy,iz)

               fk(2,ix,iy,iz) = 0.d0

            enddo
         enddo
      enddo
C Get flux in Fourier by taking transform

       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx

                skft(((iz-2)*2*ny+(iy-1))*nx+ix) = fk(1,ix,iy,iz)
         
             end do
          end do
       end do
       

       work = 0
       Call fftwnd_f77_mpi(plan,1,skft,work,0,fftw_norm_order)
      
 
       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx

              fk(1,ix,iy,iz) = 
     $        real(skft(((iz-2)*2*ny+(iy-1))*nx+ix))/cst 

              fk(2,ix,iy,iz) = 
     $        aimag(skft(((iz-2)*2*ny+(iy-1))*nx+ix))/cst
        
             end do
          end do
       end do
 

C CALCULATING THE SINK TERM IN THE OXYGEN TRANSPORT EQ
      cs        = 0.d0
       cs(1,:,:,:) = ci(:,:,:)*((skc(1,:,:,:)-0.1)/0.9)**(0.2)

       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx

                skft(((iz-2)*2*ny+(iy-1))*nx+ix) = cs(1,ix,iy,iz)
         
		 if (skc(1,ix,iy,iz).le.0.1) then
	                skft(((iz-2)*2*ny+(iy-1))*nx+ix) = 0.d0
                endif


             end do
          end do
       end do
       

       work = 0
       Call fftwnd_f77_mpi(plan,1,skft,work,0,fftw_norm_order)
      
 
       do iz = 2,nz-1
          do iy = 1,2*ny
             do ix = 1,nx

              cs(1,ix,iy,iz) = 
     $        real(skft(((iz-2)*2*ny+(iy-1))*nx+ix))/cst 

              cs(2,ix,iy,iz) = 
     $        aimag(skft(((iz-2)*2*ny+(iy-1))*nx+ix))/cst
        
             end do
          end do
       end do

C Update scalar field using explicit Euler scheme
C and copy flux      
c      D = 0.0003
      do iz = Kmin,Kmax
         kz = 2*pi*(iz-1)
         if (iz.gt.(totnz/2)) kz = 2*pi*(iz-1-totnz)
         kz = kz/Lz

         do iy = 1,2*ny
            ky = 2*pi*(iy-1)
            if (iy.gt.(totny)) ky = 2*pi*(iy-1-2*totny)
            ky = ky/(2.d0*Ly)

            do ix = 1,nx
               kx = 2*pi*(ix-1)
               if (ix.gt.(totnx/2)) kx = 2*pi*(ix-1-totnx)
               kx = kx/Lx
           
               i = ix
               j = iy
               k = iz + 1 - Scoords(1)*(nz-2)          

               k2 = kx*kx+ky*ky+kz*kz
	       
	       
		
               sk(1,i,j,k) = (sk(1,i,j,k) - 0.5d0*dt*
     $             (3.d0*(fk(1,i,j,k)+ocr*cs(1,i,j,k))-
     $                 (fko(1,i,j,k)+ocr*cso(1,i,j,k))))/(1+dt*D*k2)
               sk(2,i,j,k) = (sk(2,i,j,k) - 0.5d0*dt*
     $             (3.d0*(fk(2,i,j,k)+ocr*cs(2,i,j,k))-
     $                 (fko(2,i,j,k)+ocr*cso(2,i,j,k))))/(1+dt*D*k2)

	       
               fko(:,i,j,k) = fk(:,i,j,k)
               cso(:,i,j,k) = cs(:,i,j,k)

            enddo
         enddo
      enddo


      return
      end






C*******************************************************************
C Function fourn:
C     Computes the fast Fourier transform of a multidimensional
C     array of complex numbers.
C     From Numerical Recipes, Press et al., 2201.
C*******************************************************************

      SUBROUTINE fourn(data,nn,ndim,isign) 

      INTEGER isign,ndim,nn(ndim) 
      real data(*) 
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2 
      INTEGER ip3,k1,k2,n,nprev,nrem,ntot 
      real tempi,tempr 
      REAL theta,wi,wpi,wpr,wr,wtemp 

      ntot=1 
      do idim=1,ndim 
         ntot=ntot*nn(idim) 
      enddo 
      nprev=1 
      do idim=1,ndim 
         n=nn(idim) 
         nrem=ntot/(n*nprev) 
         ip1=2*nprev 
         ip2=ip1*n 
         ip3=ip2*nrem 
         i2rev=1 
         do i2=1,ip2,ip1 
            if(i2.lt.i2rev) then 
               do i1=i2,i2+ip1-2,2 
                  do i3=i1,ip3,ip2 
                     i3rev=i2rev+i3-i2 
                     tempr=data(i3) 
                     tempi=data(i3+1) 
                     data(i3)=data(i3rev) 
                     data(i3+1)=data(i3rev+1) 
                     data(i3rev)=tempr 
                     data(i3rev+1)=tempi 
                  enddo 
               enddo 
            endif 
            ibit=ip2/2 
 1          if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then 
               i2rev=i2rev-ibit 
               ibit=ibit/2 
               goto 1 
            endif 
            i2rev=i2rev+ibit 
         enddo 
         ifp1=ip1  
 2       if(ifp1.lt.ip2) then 
            ifp2=2*ifp1 
            theta=isign*6.2831852071795900/(ifp2/ip1) 
            wpr=-2.00*sin(0.500*theta)**2 
            wpi=sin(theta) 
            wr=1.00 
            wi=0.00 
            do i3=1,ifp1,ip1 
               do i1=i3,i3+ip1-2,2 
                  do i2=i1,ip3,ifp2 
                     k1=i2 
                     k2=k1+ifp1 
                     tempr=(wr)*data(k2)-(wi)*data(k2+1) 
                     tempi=(wr)*data(k2+1)+(wi)*data(k2) 
                     data(k2)=data(k1)-tempr 
                     data(k2+1)=data(k1+1)-tempi
                     data(k1)=data(k1)+tempr 
                     data(k1+1)=data(k1+1)+tempi 
                  enddo 
               enddo 
               wtemp=wr  
               wr=wr*wpr-wi*wpi+wr 
               wi=wi*wpr+wtemp*wpi+wi 
            enddo 
            ifp1=ifp2 
            goto 2 
         endif 
         nprev=n*nprev 
      enddo 
      return 
      END


C*******************************************************************
C Function ran:   
C     Random number generator.
C*******************************************************************

      function ran(idum)

      implicit none
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real ran,am,eps,rnmx
      parameter (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     &     ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     &     ir2=3791,ntab=16,ndiv=1+imm1/ntab,eps=1.2E-7,rnmx=1.-eps)
      integer idum2,j,k,iv(ntab),iy
      save iv,iy,idum2
      data idum2/123456789/, iv/ntab*0/, iy/0/
      
      if (idum.le.0) then
         idum = max(-idum,1)
         idum2 = idum
         do j = ntab+8,1,-1
            k = idum/iq1
            idum = ia1*(idum-k*iq1)-k*ir1
            if (idum.lt.0) idum = idum+im1
            if (j.le.ntab) iv(j) = idum
         enddo
         iy = iv(1)
      endif
      k = idum/iq1
      idum = ia1*(idum-k*iq1)-k*ir1
      if (idum.lt.0) idum = idum+im1
      k = idum2/iq2
      idum2 = ia2*(idum2-k*iq2)-k*ir2
      if ( idum2.lt.0) idum2 = idum2+im2
      j = 1+iy/ndiv
      iy = iv(j)-idum2
      iv(j) = idum
      if (iy.lt.1) iy = iy+imm1
      ran = min(am*iy,rnmx)

      return
      end
