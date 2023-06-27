!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THIS VERSION MAKES A LOOP OVER REALIZATIONS
! This Fortran code is based on the Geophysical High-Order Suite for Turbulence
! (GHOST: http://wp.df.uba.ar/mininni/ghost/) by Pablo Minnini
! 
! The present version solves the 2D incompressible Navier-Stokes equations
! in a rectangular periodic domain of size [2pi/Qx, 2pi/Qy]. 
!
! It outputs energy and enstrophy balance files,
! energy spectra, fluxes, stream function and vorticity fields, and 
! calculatures polarization using two measures: 
!
! m = (<v^2> - <u^2>) / (<u^2> + <v^2>)
! 
! or, alternatively, E_ls_y, the kinetic energy contained in the modes (kx,ky) = (0,1), (0,-1), where
!
! m ~ 0  E_ls_y >> 1 correspond to a strong large-scale vortex, while
! m > 0, E_ls_y << 1 correspond to a pure unidirectional flow in the y direction (e.g. for Qx<<1, Qy=1)
!
! Modified by Adrian van Kan (06/02/2023)

!=================================================================
      PROGRAM ANIS2D
!=================================================================
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!=================================================================

      USE mpivars
      USE fft
      USE ali
      USE var
      USE kes
      USE grid
      USE random
      use iso_fortran_env
      IMPLICIT NONE

!
! Integration parameters
!     ord  : order of the Runge-Kutta method used

      INTEGER, PARAMETER    :: ord = 4
      INTEGER               :: ini
      INTEGER(kind=int64)   :: step
      INTEGER               :: tstep
      INTEGER               :: cstep
      INTEGER               :: sstep

!
! streamfunction, vector potential, z component 
! of the fields and external force matrixes


      DOUBLE COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: ps   ! 2D streamfunction
      DOUBLE COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: fp
 
!
! Temporal data storing matrixes

      DOUBLE COMPLEX,   ALLOCATABLE, DIMENSION (:,:)    :: C1
      DOUBLE COMPLEX,   ALLOCATABLE, DIMENSION (:,:)    :: C2
      DOUBLE COMPLEX,   ALLOCATABLE, DIMENSION (:,:)    :: C3
      DOUBLE COMPLEX,   ALLOCATABLE, DIMENSION (:,:)    :: C4
      DOUBLE COMPLEX,   ALLOCATABLE, DIMENSION (:,:)    :: C5
      DOUBLE COMPLEX,   ALLOCATABLE, DIMENSION (:,:)    :: C6

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:)    :: R1
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:)    :: R2
!
! Some auxiliary matrixes

      DOUBLE PRECISION :: ener,Ekx, Eky, polar1
      DOUBLE PRECISION :: enerv,enst
      DOUBLE PRECISION :: eneru,jnst
      DOUBLE PRECISION :: dt,dt_new,CFL,dt_corr,timecorr
      DOUBLE PRECISION :: kup,kdn,kr
      DOUBLE PRECISION :: kmup,kmdn
      DOUBLE PRECISION :: prm1,prm2
      DOUBLE PRECISION :: dump,tmp
      DOUBLE PRECISION :: tmp1,tmp2,tmp3,tmp4,tmp5
      DOUBLE PRECISION :: fp0,u0
      DOUBLE PRECISION :: time
      DOUBLE PRECISION :: nu,hnu
      DOUBLE PRECISION :: phase1,phase2
      DOUBLE PRECISION :: phase3,phase4
      DOUBLE PRECISION :: tmpx,tmpy
      INTEGER :: Nreal = 1000
      INTEGER :: stat, start, ireal
      INTEGER :: t,o,nn,mm
      INTEGER :: i,j,ir,jr
      INTEGER :: ki,kj
      INTEGER :: ic,id,iu
      INTEGER :: jc,jd,ju,jt
      INTEGER :: timet,timec,times,timec2
      INTEGER :: seed,iflow, seed1

      CHARACTER     :: c,d,u,th
      CHARACTER*3   :: node,ext
      CHARACTER*4   :: ext4
      CHARACTER*100 :: odir
      CHARACTER*100 :: idir

!
! Initializes the MPI library

      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      ic = 48+int(myrank/100)
      id = 48+int(myrank/10)-int(myrank/100)*10
      iu = 48+int(myrank)-int(myrank/10)*10
      c = char(ic)
      d = char(id)
      u = char(iu)
      node = c // d // u

!
! Allocates memory for distributed blocks

      CALL range(1,nx/2+1,nprocs,myrank,ista,iend)
      CALL range(1,ny,nprocs,myrank,jsta,jend)

      ALLOCATE( R1(nx,jsta:jend) )
      ALLOCATE( R2(nx,jsta:jend) )
      ALLOCATE( C1(ny,ista:iend) )
      ALLOCATE( C2(ny,ista:iend) )
      ALLOCATE( C3(ny,ista:iend) )
      ALLOCATE( C4(ny,ista:iend) )
      ALLOCATE( C5(ny,ista:iend) )
      ALLOCATE( C6(ny,ista:iend) )
      ALLOCATE( ps(ny,ista:iend) )
      ALLOCATE( fp(ny,ista:iend) )
      ALLOCATE( kx(nx),ky(ny), kk2(ny,ista:iend), kn2(ny,ista:iend) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Reads from the external file 'status.txt'
! the status of a previous run (if any)
!     stat: last output of a previous run

      IF (myrank.eq.0) THEN
         OPEN(1,file='status.inp',status='unknown')
         READ(1,*) stat
         READ(1,*) time
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(stat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(time,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      IF (myrank.eq.0) THEN
         print*,'READING PARAMETERS'     
         OPEN(1,file='parameter.inp',status='unknown')
         READ(1,*) Qx                       ! 1
         print*,'Qx ',Qx
         READ(1,*) Qy                       ! 2
         print*,'Qy ',Qy
         READ(1,*) DQ                       ! 3
         print*,'DQ ',DQ
         READ(1,*) CFL                      ! 4
         print*,'CFL ',CFL
         READ(1,*) step                     ! 5
         print*,'step ',step
         READ(1,*) cstep                    ! 6
         print*,'cstep ',cstep
         READ(1,*) sstep                    ! 7
         print*,'sstep ',sstep         
         READ(1,*) tstep                    ! 8
         print*,'tstep ',tstep         
         READ(1,*) fp0                      ! 9
         print*,'fp0 ',fp0
         READ(1,*) u0                       ! 10
         print*,'u0 ',u0
         READ(1,*) kdn                      ! 11
         READ(1,*) kup                      ! 12
         print*,'kdn,kup ',kdn,kup
         READ(1,*) nu                       ! 13
         READ(1,*) hnu                      ! 14
         READ(1,*) nn                       ! 15
         READ(1,*) mm                       ! 16
         print*,'nu,hnu,nn,mm ',nu,hnu,nn,mm         
         READ(1,*) seed                     ! 17
         print*,'seed ',seed
         READ(1,*) iflow                    ! 18
         print*,'iflow ',iflow
         READ(1,*) start                    ! 19
         print*,'start ',start              
         READ(1,'(a100)') idir              ! binary input directory
         READ(1,'(a100)') odir              ! output directory
         CLOSE(1)
!         step = step
!         tstep = tstep
!         sstep = sstep
!         cstep = cstep
      ENDIF
      CALL MPI_BCAST(Qx   ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 1
      CALL MPI_BCAST(Qy   ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 2
      CALL MPI_BCAST(DQ   ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 3
      CALL MPI_BCAST(  CFL,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 4
      CALL MPI_BCAST( step,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 5
      CALL MPI_BCAST(cstep,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 6
      CALL MPI_BCAST(sstep,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 7
      CALL MPI_BCAST(tstep,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 8
      CALL MPI_BCAST(  fp0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 9
      CALL MPI_BCAST(   u0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 10
      CALL MPI_BCAST(  kdn,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 11
      CALL MPI_BCAST(  kup,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 12
      CALL MPI_BCAST(   nu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 13
      CALL MPI_BCAST(  hnu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) ! 14
      CALL MPI_BCAST(   nn,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 15
      CALL MPI_BCAST(   mm,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 16
      CALL MPI_BCAST( seed,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 17
      CALL MPI_BCAST(iflow,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 18
      CALL MPI_BCAST(start,1,MPI_INTEGER,         0,MPI_COMM_WORLD,ierr) ! 19
      CALL MPI_BCAST(idir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

!
! Some numerical constants

      ic = 48
      id = 48
      iu = 48
      jt = 48
      jc = 48
      jd = 48
      ju = 48

!
! Some constants for the FFT
!     kmax: maximum truncation for dealiasing
!     tiny: minimum truncation for dealiasing


      !obn  = 1.0d0/9.0d0
      !kmax = (dble(n)/3.d0)**2
      kmax = 1.0d0/9.0d0 !maximum WN when each direction scaled by 1/n_i**2
      nmax = int(max(nx*Qx/DQ,ny*Qy/DQ))
      tiny =  0.000001d0

!
! Builds the wave number and the square wave 
! number matrixes
!!!!!!!!!!!!!!!!!
!  (Part A) Non-dimentional WN for dialiasing and testing convergence
      DO i = 1,nx/2
         kx(i)      = dble(i-1)
         kx(i+nx/2) = dble(i-nx/2-1)
      END DO
      DO j = 1,ny/2
         ky(j)      = dble(j-1)
         ky(j+ny/2) = dble(j-ny/2-1)
      END DO
      tmpx = 1.0d0 / dble(nx)**2
      tmpy = 1.0d0 / dble(ny)**2
      DO i = ista,iend
         DO j = 1,ny
         kn2(j,i) = 0.1d0
         if ((kx(i)**2*tmpx+ky(j)**2*tmpy).ge.kmax) kn2(j,i) = 1.0d0
         if (abs(kx(i))+abs(ky(j)).lt.tiny) kn2(j,i) = 0.0d0
         ENDDO
      ENDDO 
!!!!!!!!!!!!!!
!  (Part B) Dimensional WN for taking derivatives
      DO i = 1,nx/2
         kx(i)      = dble(i-1     )*Qx
         kx(i+nx/2) = dble(i-nx/2-1)*Qx
      END DO
      DO j = 1,ny/2
         ky(j)      = dble(j-1     )*Qy
         ky(j+ny/2) = dble(j-ny/2-1)*Qy
      END DO
      DO i = ista,iend
         DO j = 1,ny
               kk2(j,i) = kx(i)**2     +  ky(j)**2   
         END DO
      END DO
!!!!!!!!!!!!!!!!!!!!!!!

!print*,'DBG: before create_plan'

!
! Initializes the FFT library
! Use FFTW_ESTIMATE in short runs and FFTW_MEASURE 
! in long runs

      CALL fftp2d_create_plan(planrc,nx,ny,FFTW_REAL_TO_COMPLEX, &
                             FFTW_MEASURE)
      CALL fftp2d_create_plan(plancr,nx,ny,FFTW_COMPLEX_TO_REAL, &
                             FFTW_MEASURE)

!
! Sets the initial conditions.

!print*,'DBG: after create_plan'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! START OF REALIZATION LOOP
   DO ireal = 1, Nreal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      seed = 2*ireal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! INITIALIZATION OF STREAMFUNCTION !!!!!!
     IF (stat.eq.0) THEN
         ini = 1
         timet = tstep
         timec = cstep
         timec2 = cstep
         times = sstep
!         timecorr = dt_corr

!STREAM FUNCTION R1 & PHI R2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         DO j = jsta,jend
            DO i = 1,nx
               R1(i,j) = 0.0d0
               R2(i,j) = 0.0d0
            END DO
         END DO
         DO ir = 1,int(kup)+1
         DO jr = 1,int(kup)+1
            kr = sqrt(dble(ir*ir+jr*jr))
         phase1=randu(seed)*2.0d0 *pi
         phase2=randu(seed)*2.0d0 *pi
         phase3=randu(seed)*2.0d0 *pi
         phase4=randu(seed)*2.0d0 *pi
         DO j = jsta,jend
            DO i = 1,nx
               R1(i,j) = R1(i,j) &
                       + cos(2.0d0*(ir-1)*pi*(dble(i)-1)/dble(nx)+phase1) &
                       * cos(2.0d0*(jr-1)*pi*(dble(j)-1)/dble(ny)+phase2)/kr
            END DO
         END DO
         END DO
         END DO


         CALL fftp2d_real_to_complex(planrc,R1,ps,MPI_COMM_WORLD)
         CALL energy(ps,ener,1)
         CALL MPI_BCAST(ener,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
         tmp=u0/sqrt(ener)
         DO i = ista,iend
            DO j = 1,ny
               IF ((kn2(j,i).le.kmax).and.(kn2(j,i).ge.tiny)) THEN
                  ps(j,i) = tmp*ps(j,i)
               ELSE
                  ps(j,i) = 0.0d0
               ENDIF
            END DO
         END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ELSE ! stat
         print*,'READING...',stat
         ini = int((stat-1)*tstep)
         dump = dble(ini)/dble(sstep)+1
         times = 0
         timet = 0
         timec = 0
         timec2 = 0
         timecorr = 0 
        
         jt = 48+int(dump/1000)
         jc = 48+int(dump/100)-int(dump/1000)*10
         jd = 48+int(dump/10)-int(dump/100)*10
         ju = 48+int(dump)-int(dump/10)*10

         ic = 48+int(float(stat)/100)
         id = 48+int(float(stat)/10)-int(float(stat)/100)*10
         iu = 48+int(stat)-int(float(stat)/10)*10
         c = char(ic)
         d = char(id)
         u = char(iu)

         OPEN(1,file=trim(idir) // '/ps.' // node // '.' &
                           // c // d // u //'.out',form='unformatted')
         READ(1) R1
         CLOSE(1) 
         CALL fftp2d_real_to_complex(planrc,R1,ps,MPI_COMM_WORLD)

        !CALL energy(ps,ener,1)
        !  IF (myrank.eq.0) THEN
        !   print*, "DBG:",ener
        !  ENDIF


      ENDIF ! stat

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                         FORCING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !DETERMINISTIC FORCING
         IF (iflow.eq.1) THEN
         kdn = kup
         DO j = jsta,jend
            DO i = 1,nx
               R1(i,j) = cos(2*kup/Qx*pi*(dble(i)-1)/dble(nx)) !&
                       !+ sin(2*kup/Qy*pi*(dble(j)-1)/dble(ny))
            END DO
         END DO
         CALL fftp2d_real_to_complex(planrc,R1,fp,MPI_COMM_WORLD)
         CALL energy(fp,eneru,1)
         CALL MPI_BCAST(eneru,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr) 
         tmp=fp0/sqrt(eneru)
         DO i = ista,iend
            DO j = 1,ny
               IF ((kn2(j,i).le.kmax).and.(kn2(j,i).ge.tiny)) THEN
                  fp(j,i) = tmp*fp(j,i)
               ELSE
                  fp(j,i) = 0.0d0
               ENDIF
            END DO
         END DO
         
         ELSEIF (iflow.eq.2) THEN
         seed1 = seed+1
         CALL const_inj(ps,kdn,kup,fp0,fp,1,seed1)

         ELSEIF (iflow.eq.3) THEN
         seed1 = seed+1
         CALL CFL_condition(CFL,ps,nu,nn,dt)
         CALL rand_force(kdn,kup,fp0,dt,seed,1,fp)
         ELSEIF (iflow.eq.4) THEN
         seed1 = seed+1
         CALL CFL_condition(CFL,ps,nu,nn,dt)
         CALL rand_force_broad(kdn,kup,fp0,dt,seed,1,fp)
         ENDIF ! iflow
         seed = seed + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         CALL energy(ps,ener,1)
!           IF (myrank.eq.0) THEN
        print*, "DBG pre RK:",ener
!          ENDIF
!
! Time integration scheme starts here
! Uses Runge-Kutta of order 'ord'
!#################### MAIN LOOP ######################


 RK : DO t = ini,step

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! CHECK CONDITION FOR BREAK
    !! COMPUTE, <u^2> 
       CALL derivk2(ps,C6,2)
       CALL energy(C6,Ekx,0)

    !! COMPUTE <v^2>
       CALL derivk2(ps,C6,1)
       CALL energy(C6,Eky,0)

    !! COMPUTE M = (<v^2> - <u^2>)/(<u^2> + <v^2>), such that m ~ 1 for Lx >> Ly
       polar1 = (Eky - Ekx)/(Ekx + Eky)
     
       IF (start.eq.1) THEN !IF WE START IN THE LSV
            IF (polar1.gt.0.95) THEN
               IF (myrank.eq.0) THEN
                  OPEN(1,file='DecayLSV_Time.txt',position='append')
                  WRITE(1,*) ireal, time
                  CLOSE(1)           
               ENDIF
               GOTO 150
            ENDIF

       ELSEIF (start.eq.0) THEN  !IF WE START IN THE JET      
            IF (polar1.lt.0.2) THEN
               IF (myrank.eq.0) THEN
                  OPEN(1,file='DecayJet_Time.txt',position='append')
                  WRITE(1,*) ireal, time
                  CLOSE(1)
               ENDIF
               GOTO 150
            ENDIF

       ENDIF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         CALL energy(ps,ener,1)
!        IF (myrank.eq.0) THEN
!        print*,"DBG top RK pre CFL",ener
!         ENDIF
       CALL CFL_condition(CFL,ps,nu,nn,dt)
! Every 'cstep' steps, generates external files 
! to check consistency and convergence. See the 
! cond_check subroutine for details.
!         CALL energy(ps,ener,1)
!           IF (myrank.eq.0) THEN
!        print*,"DBG top RK",ener
!         ENDIF
          IF (timec.eq.cstep) THEN   
              timec = 0
              CALL cond_check(ps,fp,time,nn,nu,mm,hnu,kup)
          ENDIF

!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!
!   FORCING UPDATE IN RK LOOP 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          IF (iflow.eq.2) THEN
          seed1 = seed+1
          CALL const_inj(ps,kdn,kup,fp0,fp,1,seed1)
          ENDIF
          IF (iflow.eq.3) THEN
          CALL rand_force(kdn,kup,fp0,dt,seed,1,fp)
          ENDIF
          IF (iflow.eq.4) THEN
          CALL rand_force_broad(kdn,kup,fp0,dt,seed,1,fp)
          ENDIF ! iflow3 
!         CALL energy(ps,ener,1)
!          IF (myrank.eq.0) THEN
!        print*, "DBG post iflow:",ener
!          ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Every 'sstep' steps, generates external files 
! with the power spectrum

         IF (times.eq.sstep) THEN
            times = 0
            ju = ju+1
            IF (ju.eq.58) THEN
               ju = 48
               jd = jd+1
            ENDIF
            IF (jd.eq.58) THEN
               jd = 48
               jc = jc+1
            ENDIF
            IF (jc.eq.58) THEN
               jc = 48
               jt = jt+1
            ENDIF
            th= char(jt)
            c = char(jc)
            d = char(jd)
            u = char(ju)
            ext4 = th // c // d // u 
           CALL spectrum(ps,ext4,odir)
           CALL transfers(ps,ext4,odir)

           IF (myrank.eq.0) THEN
            OPEN(1,file='time_spec.txt',position='append')
            WRITE(1,13) ext4,time
   13       FORMAT( A4,    F12.3)
            CLOSE(1)
           ENDIF
         ENDIF

! Every 'tstep' steps, stores the results of the integration

         IF (timet.eq.tstep) THEN
            timet = 0
            iu = iu+1
            IF (iu.eq.58) THEN
               iu = 48
               id = id+1
            ENDIF
            IF (id.eq.58) THEN
               id = 48
               ic = ic+1
            ENDIF
            IF (jc.eq.58) THEN
               jc = 48
               jt = jt+1
            ENDIF 
            th= char(jt)           
            c = char(ic)
            d = char(id)
            u = char(iu)
            ext = th // c // d // u
            tmp=1.0d0/dble(nx)/dble(ny)
            DO i = ista,iend
               DO j = 1,ny
                  C1(j,i) = ps(j,i)*tmp
               END DO
            END DO
            CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
            OPEN(1,file=trim(odir) // '/ps.' // node // '.' &
                // c // d // u // '.out',form='unformatted')
            WRITE(1) R1
            CLOSE(1)
            DO i = ista,iend
               DO j = 1,ny
                  C1(j,i) = ps(j,i)*kk2(j,i)*tmp
               END DO
            END DO
            CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
            OPEN(1,file=trim(odir) // '/ww.' // node // '.' &
                // c // d // u // '.out',form='unformatted')
            WRITE(1) R1
            CLOSE(1)
           IF (myrank.eq.0) THEN
           OPEN(1,file='time_field.txt',position='append')
           WRITE(1,12) c//d//u,time
!   10      FORMAT( I10,   A30)
!   11      FORMAT( F12.6, A30)
   12      FORMAT( A4,    F12.3) 
           CLOSE(1)
           ENDIF      

         ENDIF

         timet = timet+1
         times = times+1
         timec = timec+1
         timec2 = timec2+1
!         timecorr = timecorr+dt
         time = time+dt
!         CALL energy(ps,ener,1)
!         CALL energy(phi,enerphi,0)
!           IF (myrank.eq.0) THEN
!        print*, "DBG post sstep:",ener,enerphi
!          ENDIF


! Runge-Kutta step 1
! Copies the streamfunction into the auxiliary matrix C1

         DO i = ista,iend
            DO j = 1,ny
               C1(j,i) = ps(j,i)
            END DO
         END DO


! Runge-Kutta step 2

         DO o = ord,1,-1
!         CALL energy(C1,ener,1)
!         CALL energy(C2,enerphi,0)
!           IF (myrank.eq.0) THEN
!            print*,"DBG top ord",ener,enerphi
!           endif
       seed = seed +1 
       CALL laplak2(C1,C4)               ! make - w_z
       CALL poisson(C1,C4,C5)            ! - ez.curl(u_z x w_z)

         DO i = ista,iend
            DO j = 1,ny
            IF ((kn2(j,i).le.kmax).and.(kn2(j,i).gt.tiny)) THEN
            tmp1 = dt/dble(o) 
            tmp2 = 1.0d0/kk2(j,i) 

            !  ps
            tmp3 = (1.0d0 +(nu*kk2(j,i)**nn + hnu*tmp2**mm)*tmp1) 
            C1(j,i) =  ps(j,i)+(-C5(j,i)*tmp2 +fp(j,i))*tmp1      ! Vorticity equation:
                                                                  ! d/dt(omega)  = - d/dt(Laplacian(psi)) = - ((u,v).nabla)omega + rest
                                                                  ! = - (u.nabla) omega +rest = J(psi,omega) + rest 
                                                                  ! = - J(psi,Laplacian(psi)) + rest
            C1(j,i) =  C1(j,i)/tmp3                               ! Uses semi-implicit scheme, see e.g.
                                                                  ! https://www.staff.uni-oldenburg.de/hannes.uecker/pre/030-mcs-hur.pdf
            ELSE  
            C1(j,i) = 0.0d0
            ENDIF
            END DO
         END DO
!         CALL energy(C1,ener,1)
!         CALL energy(C2,enerphi,0)
!           IF (myrank.eq.0) THEN
!                print*,"DBG end ord",ener,enerphi
!           ENDIF
         END DO  ! ord

! Runge-Kutta step 3
! Copies the result from the auxiliary matrix into ps

         CALL derivk2(C1,ps,0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      END DO RK
!##############  END OF MAIN LOOP ###################
  150 CONTINUE

ENDDO   ! END OF REALIZATION LOOPS
!
! End of Runge-Kutta

      CALL MPI_FINALIZE(ierr)
      CALL fftp2d_destroy_plan(plancr)
      CALL fftp2d_destroy_plan(planrc)
      DEALLOCATE( R1 )
      DEALLOCATE( R2 )
      DEALLOCATE( ps )
      DEALLOCATE( fp )
      DEALLOCATE( C1 )
      DEALLOCATE( C2 )
      DEALLOCATE( C3 )
      DEALLOCATE( C5 )
      DEALLOCATE( C6 )
      
      DEALLOCATE( kx,ky, kk2, kn2 )

      END PROGRAM ANIS2D
