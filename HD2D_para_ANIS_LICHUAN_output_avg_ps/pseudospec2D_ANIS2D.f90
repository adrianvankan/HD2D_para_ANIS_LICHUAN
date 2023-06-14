!=================================================================
! PSEUDOSPECTRAL subroutines
!=================================================================

!*****************************************************************
      SUBROUTINE derivk2(a,b,dir)
!-----------------------------------------------------------------
!
! Two-dimensional derivative of the matrix 'a'
!
! Parameters
!     a  : input matrix
!     b  : at the output contains the derivative da/dk_dir
!     dir: =1 derivative in the x-direction
!          =2 derivative in the y-direction
!          =0 just copies a to b
!
      USE ali
      USE kes
      USE var
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(ny,ista:iend) :: a,b
      INTEGER :: dir
      INTEGER :: i,j

!
! Derivative in the x-direction
!
      IF (dir.eq.1) THEN
         DO i = ista,iend
            DO j = 1,ny
               IF ((kn2(j,i).le.kmax).and.(kn2(j,i).gt.tiny)) THEN
                  b(j,i) = im*kx(i)*a(j,i)
               ELSE
                  b(j,i) = 0.0d0
               ENDIF
            END DO
         END DO
!
! Derivative in the y-direction
!
      ELSE IF (dir.eq.2) THEN
         DO i = ista,iend
            DO j = 1,ny
               IF ((kn2(j,i).le.kmax).and.(kn2(j,i).gt.tiny)) THEN
                  b(j,i) = im*ky(j)*a(j,i)
               ELSE
                  b(j,i) = 0.0d0
               ENDIF
            END DO
         END DO
      ELSE  ! copy
         DO i = ista,iend
            DO j = 1,ny
               IF ((kn2(j,i).le.kmax).and.(kn2(j,i).gt.tiny)) THEN
                  b(j,i) = a(j,i)
               ELSE
                  b(j,i) = 0.0d0
               ENDIF
            END DO
         END DO 
      ENDIF

      RETURN
      END SUBROUTINE derivk2

!*****************************************************************
      SUBROUTINE laplak2(a,b)
!-----------------------------------------------------------------
!
! Two-dimensional Laplacian of the matrix 'a'
!
! Parameters
!     a: input matrix
!     b: at the output contains the Laplacian d2a/dka2
!
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(ny,ista:iend) :: a,b
      INTEGER :: i,j

      DO i = ista,iend
         DO j = 1,ny
            b(j,i) = -kk2(j,i)*a(j,i)
         END DO
      END DO

      RETURN
      END SUBROUTINE laplak2

!*****************************************************************
      SUBROUTINE inv_laplak2(a,b)
!-----------------------------------------------------------------
!
! Two-dimensional Laplacian of the matrix 'a'
!
! Parameters
!     a: input matrix
!     b: at the output contains the Laplacian d2a/dka2
!
      USE ali
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(ny,ista:iend) :: a,b
      INTEGER :: i,j

      DO i = ista,iend
         DO j = 1,ny
            IF (kk2(j,i).gt.tiny) THEN
            b(j,i) = -a(j,i)/kk2(j,i)
            ELSE
            b(j,i) = 0.0d0
            ENDIF
         END DO
      END DO

      RETURN
      END SUBROUTINE inv_laplak2


!*****************************************************************
      SUBROUTINE curl_2D_z(a,b,c)
!-----------------------------------------------------------------
!
! Parameters
!     a:  input matrix
!     b:  input matrix
!     c: output matrix
!
      USE mpivars
      USE kes
      USE ali
      USE grid
      USE fft
      USE var
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(ny,ista:iend) :: a,b,c
      DOUBLE COMPLEX, DIMENSION(ny,ista:iend) :: c1,c2
      INTEGER :: i,j
      DO i = ista,iend
         DO j = 1,ny
            IF ((kn2(j,i).le.kmax)) THEN
               c(j,i) = im*( kx(i)*b(j,i) - ky(j)*a(j,i) )
            ELSE
               c(j,i) = 0.0d0
            ENDIF
         END DO
      END DO
      RETURN
      END SUBROUTINE curl_2D_z



!#################################################################
!#####################    NONLINEARITIES   #######################
!#################################################################

!*****************************************************************
      SUBROUTINE poisson(a,b,c)
!-----------------------------------------------------------------
!
! Poisson bracket of the scalar fields A and B 
! in real space.
!
! Parameters
!     a: input matrix
!     b: input matrix
!     c: Poisson bracket {A,B} [output]
!
      USE mpivars
      USE kes
      USE ali
      USE grid
      USE fft
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(ny,ista:iend) :: a,b,c
      DOUBLE COMPLEX, DIMENSION(ny,ista:iend) :: c1,c2
      DOUBLE PRECISION, DIMENSION(nx,jsta:jend)    :: r1,r2,r3
      INTEGER :: i,j

!
! Computes dA/dx.dB/dy
!
      CALL derivk2(a,c1,1)
      CALL derivk2(b,c2,2)
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      DO j = jsta,jend
         DO i = 1,nx
            r3(i,j) = r1(i,j)*r2(i,j)
         END DO
      END DO

!
! Computes dA/dy.dB/dx
!
      CALL derivk2(a,c1,2)
      CALL derivk2(b,c2,1)
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      DO j = jsta,jend
         DO i = 1,nx
            r3(i,j) = (r3(i,j)-r1(i,j)*r2(i,j))/dble(nx)**2/dble(ny)**2
         END DO
      END DO
      CALL fftp2d_real_to_complex(planrc,r3,c,MPI_COMM_WORLD)
      DO i = ista,iend
         DO j = 1,ny
            IF ((kn2(j,i).ge.kmax).or.(kn2(j,i).le.tiny)) THEN
               c(j,i) = 0.0d0
            ENDIF
         END DO
      END DO
      RETURN
      END SUBROUTINE poisson

!*****************************************************************
      SUBROUTINE pmult(a,b,c)
!-----------------------------------------------------------------
!
! Pointwise multiplication of the scalar fields A and B
! in real space.
!
! Parameters
!     a: input matrix
!     b: input matrix
!     c: pointwise multiplication of a*b [output]
!
      USE mpivars
      USE kes
      USE ali
      USE grid
      USE fft
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(ny,ista:iend) :: a,b,c
      DOUBLE COMPLEX, DIMENSION(ny,ista:iend) :: c1,c2
      DOUBLE PRECISION, DIMENSION(nx,jsta:jend)    :: r1,r2,r3
      INTEGER :: i,j

      CALL derivk2(a,c1,0)
      CALL derivk2(b,c2,0)

      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      DO j = jsta,jend
         DO i = 1,nx
            r3(i,j) = (r1(i,j)*r2(i,j))/dble(nx)/dble(ny)**2
         END DO
      END DO
      CALL fftp2d_real_to_complex(planrc,r3,c,MPI_COMM_WORLD)
      DO i = ista,iend
         DO j = 1,ny
            IF ((kn2(j,i).ge.kmax).or.(kn2(j,i).le.tiny)) THEN
               c(j,i) = 0.0d0
            ENDIF
         END DO
      END DO


      END SUBROUTINE pmult

!#################################################################
!########################   ANALYSIS   ###########################
!#################################################################

!*****************************************************************
      SUBROUTINE mom_calc(a,b,kin)
!-----------------------------------------------------------------
!
! Computes the <|a|^kin >
! The output is valid only in the first node.
!
! Parameters
!     a  : input matrix with the scalar field
!     b  : at the output contains the moment
!     kin: =1 computes the first moment <|a|>
!          =2 computes the second moment <|a|^2>
!          =3 computes the third moment <|a|^3>,
!          etc.
!

      USE kes
      USE grid
      USE mpivars
      USE ali
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(ny,ista:iend) :: a
      DOUBLE PRECISION  :: b
      DOUBLE PRECISION  :: bloc
      DOUBLE PRECISION  :: tmp
      INTEGER :: kin,two
      INTEGER :: i,j

      bloc = 0.0d0
      tmp = 1.0d0/dble(nx)**2/dble(ny)**2

!
! Computes the kin'th moment of the scalar field
!
        DO i = ista,iend
           two = 2
           if (i.eq.1) two = 1
               DO j = 1,ny
                  if (kn2(j,i).ge.tiny) then
                  bloc = bloc+two*abs(a(j,i))**kin*tmp
                  endif
               END DO
        END DO
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(bloc,b,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(  b,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE mom_calc

!###############################


!*****************************************************************
      SUBROUTINE energy(a,b,kin)
!-----------------------------------------------------------------
!
! Computes the mean kinetic or magnetic energy in 2D, or mean vorticity. 
! The output is valid only in the first node.
!
! Parameters
!     a  : input matrix with the scalar field (stream function)
!     b  : at the output contains the energy
!     kin: =0 computes the square of the scalar field
!          =1 computes the energy
!          =2 computes the vorticity
!
      USE kes
      USE grid
      USE mpivars
      USE ali
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(ny,ista:iend) :: a
      DOUBLE PRECISION  :: b
      DOUBLE PRECISION  :: bloc
      DOUBLE PRECISION  :: tmp
      INTEGER :: kin,two
      INTEGER :: i,j

      bloc = 0.0d0
      tmp = 1.0d0/dble(nx)**2/dble(ny)**2

!
! Computes (square of the scalar field) * k^(2kin)
!
        IF (kin.ge.0) THEN
            DO i = ista,iend
               two = 2
               if (i.eq.1) two = 1  
               DO j = 1,ny
                  bloc = bloc+two*abs(a(j,i))**2*kk2(j,i)**kin*tmp
               END DO
            END DO
        ELSE
        DO i = ista,iend
               two = 2
               if (i.eq.1) two = 1
               DO j = 1,ny
                  if (kn2(j,i).ge.tiny) then
                  bloc = bloc+two*abs(a(j,i))**2*kk2(j,i)**kin*tmp
                  endif
               END DO
            END DO
        ENDIF
!
! Computes the reduction between nodes
!
      CALL MPI_REDUCE(bloc,b,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(  b,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      RETURN
      END SUBROUTINE energy

!###############################
!*****************************************************************
      SUBROUTINE inerprod(a,b,kin,rslt)
!-----------------------------------------------------------------
! Parameters
!     a  : first  input matrix
!     b  : second input matrix
!     kin: = multiplies by the laplacian to this power

      USE kes
      USE grid
      USE mpivars
      USE ali

      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(ny,ista:iend) :: a,b
      DOUBLE PRECISION  :: tmp,tmq
      DOUBLE PRECISION  :: rslt
      INTEGER :: kin,two
      INTEGER :: i,j

      tmp = 0.0d0
      tmq = 1./dble(nx)**2/dble(ny)**2
      IF (kin.ge.0) THEN
         DO i = ista,iend
            two = 2
            if (i.eq.1) two = 1 
            DO j = 1,ny
            IF ((kn2(j,i).le.kmax).and.(kn2(j,i).ge.tiny)) THEN
!            IF ((kn2(j,i).le.kmax)) THEN
            tmp = tmp+two*(kk2(j,i)**kin)*dble(b(j,i)*conjg(a(j,i)))*tmq
            ENDIF
            END DO
         END DO
       ELSE
         DO i = ista,iend
            two = 2
            if (i.eq.1) two = 1
            DO j = 1,ny
            IF ((kn2(j,i).le.kmax).and.(kn2(j,i).ge.tiny)) THEN
!            IF ((ka2(j,i).le.kmax)) THEN
            tmp = tmp+two*(kk2(j,i)**kin)*dble(b(j,i)*conjg(a(j,i)))*tmq
            ENDIF
            END DO
         END DO
       ENDIF
      CALL MPI_REDUCE(tmp,rslt,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rslt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)  

      
      RETURN
      END SUBROUTINE inerprod
!###############################

!*****************************************************************
      SUBROUTINE cond_check(ps,fp,time,nn,nu,mm,hnu,kup)
!-----------------------------------------------------------------
!
! Condition check for conservation of energy, enstrophy, polarization in HD2D
!
      USE ali
      USE kes
      USE grid
      USE mpivars
      USE fft

      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(ny,ista:iend) ::ps,fp
      DOUBLE COMPLEX, DIMENSION(ny,ista:iend) :: C1
      DOUBLE PRECISION, DIMENSION(nx,jsta:jend)    :: R1
      DOUBLE PRECISION    :: Ep1,Ep2,Epv,Eph
      DOUBLE PRECISION    :: Ek1,Ek2
      DOUBLE PRECISION    :: Ekx, Eky, polar1, polar2tmp, polar2
      DOUBLE PRECISION    :: Dsp,Hds,Dspw,Hdsw
      DOUBLE PRECISION    :: Dspv,Hdsv,Dsph1,Hdsh1,Dsph2,Hdsh2
      DOUBLE PRECISION    :: injp1,injp2,injh1,injh2
      DOUBLE PRECISION    :: coup,maxw,minw
      DOUBLE PRECISION    :: nu,hnu,nuv,hnuv
      DOUBLE PRECISION    :: Efk, Efp, kup, kmn
      DOUBLE PRECISION    :: tmq,tmp,tmp0,tmp1,tmp2,tmp3,tmp4,two,time
      INTEGER :: i,j,nn,mm


!! ENERGY
      CALL energy(ps,Ep1,1)        ! |u|^2
      CALL inerprod(ps,fp,1,injp1) ! energy injection
      CALL energy(ps,Dsp,nn+1)     ! Dissipation 
      CALL energy(ps,Hds,1-mm)

!! ENSTROPHY
      CALL energy(ps,Ep2,2)  ! |w|^2
      CALL inerprod(ps,fp,2,injp2) ! enstrophy injection
      CALL energy(ps,Dspw,nn+2)
      CALL energy(ps,Hdsw,2-mm)


!! COMPUTE, u^2 
     CALL derivk2(ps,C1,2)
     CALL energy(C1,Ekx,0)
     
!! COMPUTE v^2
     CALL derivk2(ps,C1,1)
     CALL energy(C1,Eky,0)

!! COMPUTE M = (<v^2> - <u^2>)/(<u^2> + <v^2>), such that m ~ 1 for Lx >> Ly
     polar1 = (Eky - Ekx)/(Ekx + Eky)

!! COMPUTE FRACTION OF KINETIC ENERGY IN |kx| = 0, |ky| = 1.
      tmp       = 1.0d0/dble(nx)**2/dble(ny)**2
      polar2    = 0.0d0
      polar2tmp = 0.0d0
     
      
      DO i = ista,iend
         two=2.0d0
         IF (i.eq.1) two=1.0d0
         DO j = 1,ny
            IF (((abs(abs(ky(j)) - 1.0d0)).lt.tiny).and.(abs(kx(i)).lt.tiny)) THEN
               polar2tmp = polar2tmp + two*kk2(j,i)*abs(ps(j,i))**2*tmp/(Ekx + Eky)
            ENDIF
         END DO
      END DO
      CALL MPI_REDUCE(polar2tmp,polar2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

!!!!!!!!!!! Computes the energy at largest scale!!!!!!!!!!!!!!
!      tmp = 1.0d0/dble(nx)**2/dble(ny)**2
!      tmp1=0.0d0
!      DO i = ista,iend
!         two=2.0d0
!         IF (i.eq.1) two=1.0d0
!         DO j = 1,ny
!            kmn = int(sqrt(k2(j,i))+.5d0)
!            IF ((ka2(j,i).le.(2.01)).and.(ka2(j,i).ge.0)) THEN
!            IF ((kmn.gt.0).and.(kmn.le.1)) THEN
!               tmp1 = tmp1+two*ka2(j,i)*abs(ps(j,i))**2*tmp
!            ENDIF
!         END DO
!      END DO
!      CALL MPI_REDUCE(tmp1,Ek1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
!                      MPI_COMM_WORLD,ierr)
!!!!
!!!!!!!!!!! Computes the energy at second-largest scale!!!!!!!!!!!!!!
!      tmp = 1.0d0/dble(n)**4
!      tmp2=0.0d0
!      DO i = ista,iend
!         two=2.0d0
!         IF (i.eq.1) two=1.0d0
!         DO j = 1,n
!            kmn = int(sqrt(ka2(j,i))+.5d0)
!            IF ((ka2(j,i).le.(4.01)).and.(ka2(j,i).ge.(2.01))) THEN
!            IF ((kmn.gt.1).and.(kmn.le.2)) THEN
!               tmp2 = tmp2+two*ka2(j,i)*abs(ps(j,i))**2*tmp
!            ENDIF
!         END DO
!      END DO
!      CALL MPI_REDUCE(tmp2,Ek2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
!                      MPI_COMM_WORLD,ierr)
!!!!

!!!!!!!!!!! Computes the energy at KUP scale!!!!!!!!!!!!!!
!      tmp = 1.0d0/dble(n)**4
!      tmp1=0.0d0
!      tmp2=0.0d0
!      DO i = ista,iend
!         two=2.0d0
!         IF (i.eq.1) two=1.0d0
!         DO j = 1,n
!            IF ((ka2(j,i).le.(2.01)*kup*kup).and.(ka2(j,i).ge.kup*kup )) THEN
!               tmp1 = tmp1+two*ka2(j,i)*abs(ps(j,i))**2*tmp
!            ENDIF
!         END DO
!      END DO
!      CALL MPI_REDUCE(tmp1,Efp,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
!                      MPI_COMM_WORLD,ierr)

!!!!!!!! COMPUTING MAX AND MIN VALUES !!!!!!!!!!!!
!!!!! FINDING MIN/MAX vorticity
!      CALL laplak2(ps,C1)               ! make - W_2D
!      CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
!      tmp=maxval(R1) !max vorticity
!      tmp1=minval(R1) !min vorticity
!      call MPI_REDUCE(tmp,maxw,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
!      call MPI_REDUCE(tmp1,minw,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)


!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='energy_bal.txt',position='append')
         WRITE(1,20) time, Ep1, injp1, nu*Dsp, hnu*Hds
   20    FORMAT(E23.14E3,E23.14E3,E23.14E3,E23.14E3, E23.14E3)
         CLOSE(1)
         OPEN(1,file='enstrophy_bal.txt',position='append')
         WRITE(1,22) time, Ep2, injp2, nu*Dspw, hnu*Hdsw
   22    FORMAT( E23.14E3, E23.14E3, E23.14E3, E23.14E3, E23.14E3 ) 
         CLOSE(1)
         
         OPEN(1,file='polarization.txt',position='append')
         WRITE(1,24) time, polar1, polar2
   24    FORMAT( E23.14E3, E23.14E3, E24.14E3)
         CLOSE(1)
      ENDIF
      
      RETURN
      END SUBROUTINE cond_check

!*****************************************************************
      SUBROUTINE spectrum(ps,ext,odir)
!-----------------------------------------------------------------
!
! Computes the energy power spectrum in 2D. 
! The output is written to a file by the first node.
!
! Parameters
!     a  : streamfunction 
!     b  : vector potential
!     ext: the extension used when writting the file
!     kin: =1 computes the kinetic spectrum
!          =0 computes the magnetic spectrum
!
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nmax/2+1)      :: Ek,Ektot1
      DOUBLE COMPLEX, DIMENSION(ny,ista:iend)    :: ps
      DOUBLE PRECISION        :: tmp,two,tmp1,tmp2,tmp3
      INTEGER     :: kin
      INTEGER     :: kmn
      INTEGER     :: i,j
      CHARACTER*4 :: ext
      CHARACTER*100 :: odir

      tmp = 1.0d0/dble(nx)**2/dble(ny)**2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Computes the VELOCITY energy spectrum!!!!!!!!!!!!!!
      DO i = 1,nx/2+1
         Ek(i) = 0.0d0
      END DO
      DO i = ista,iend
         two=2.0d0
         IF (i.eq.1) two=1.0d0
         DO j = 1,ny
            kmn = int(sqrt(kk2(j,i))/DQ+.5d0)
            IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
               Ek(kmn) = Ek(kmn)+two*kk2(j,i)*abs(ps(j,i))**2*tmp
            ENDIF
         END DO
      END DO
      CALL MPI_REDUCE(Ek,Ektot1,nmax/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

!!!!!!! FOR DEBUGGING: spectrum of forcing !!!!!

!      DO i = 1,n/2+1
!         Fk(i) = 0.0d0
!      END DO
!      DO i = ista,iend
!         two=2.0d0
!         IF (i.eq.1) two=1.0d0
!         DO j = 1,n
!            kmn = int(sqrt(ka2(j,i))+.5d0)
!            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
!               Fk(kmn) = Fk(kmn)+two*ka2(j,i)*abs(fp(j,i))**2*tmp
!            ENDIF
!         END DO
!      END DO
!      CALL MPI_REDUCE(Fk,Fktot1,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
!                      MPI_COMM_WORLD,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(odir) // '/spectrum.' // ext // '.txt')
!         WRITE(1,30) Q,0.0d0,0.5d0*tmp1
         do i=1,nx/2+1
         WRITE(1,30) Ektot1(i)
         enddo
         CLOSE(1)
      ENDIF
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  30    FORMAT( E24.15E3,E24.15E3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RETURN
      END SUBROUTINE spectrum

!*****************************************************************
      SUBROUTINE transfers(ps,ext,odir)
!-----------------------------------------------------------------
!
! Computes the energy power spectrum in 2D. 
! The output is written to a file by the first node.
!
! Parameters
!     a  : streamfunction 
!     b  : vector potential
!     ext: the extension used when writting the file
!     kin: =1 computes the kinetic spectrum
!          =0 computes the magnetic spectrum
!
      USE kes
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nmax/2+1)        :: Fl
      DOUBLE PRECISION, DIMENSION(nmax/2+1)        :: Fl0,Fl1,Fl2
      DOUBLE COMPLEX, DIMENSION(ny,ista:iend)    :: ps,vz
      DOUBLE COMPLEX, DIMENSION(ny,ista:iend)    :: c1,c2
      DOUBLE PRECISION        :: tmp,two,tmp1,fx0,fx1,fx2
      INTEGER     :: kin
      INTEGER     :: kmn
      INTEGER     :: i,j
      CHARACTER*4 :: ext
      CHARACTER*100 :: odir

      tmp = 1.0d0/dble(nx)**2/dble(ny)**2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         DO i = ista,iend
            DO j = 1,ny
               C1(j,i) = 0.0d0
            END DO
         END DO

         DO i = ista,iend
            DO j = 1,ny
               C2(j,i) = 0.0d0
            END DO
         END DO
         CALL laplak2(ps,c1)               ! make - W_z
         CALL poisson(ps,c1,c2)            ! - curl(u x W_z) = [psi,W_z]

!!!!!!!!!!!!!!!!!!!!!!!!    Enstrophy flux 
      DO i = 1,nx/2+1
         Fl(i) = 0.0d0
      END DO
      DO i = ista,iend
         two=2.0d0
         IF (i.eq.1) two=1.0d0
         DO j = 1,ny
            kmn = int(sqrt(kk2(j,i))/DQ+.5d0)
            IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
            Fl(kmn) = Fl(kmn)+two*kk2(j,i)*dble(ps(j,i)*conjg(c2(j,i)))*tmp
            ENDIF
         END DO
      END DO
      CALL MPI_REDUCE(Fl,Fl0,nmax/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!    Energy flux: u_2D
      DO i = 1,nx/2+1
         Fl(i) = 0.0d0
      END DO
      DO i = ista,iend
         two=2.0d0
         IF (i.eq.1) two=1.0d0
         DO j = 1,ny
            kmn = int(sqrt(kk2(j,i))/DQ+.5d0)
            IF ((kmn.gt.0).and.(kmn.le.nmax/2+1)) THEN
            Fl(kmn) = Fl(kmn)+two*dble(ps(j,i)*conjg(c2(j,i)))*tmp
            ENDIF
         END DO
      END DO
      CALL MPI_REDUCE(Fl,Fl1,nmax/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      IF (myrank.eq.0) THEN
         fx0=0.0d0 
         fx1=0.0d0
         fx2=0.0d0 
         OPEN(1,file=trim(odir) // '/transfer.' // ext // '.txt')
         do i=1,nx/2+1
         WRITE(1,40) Fl0(i), Fl1(i), Fl2(i)
         enddo
         CLOSE(1)
         OPEN(1,file=trim(odir) // '/fluxes.' // ext // '.txt')
         do i=1,nx/2+1
         fx0=Fl0(i)+fx0
         fx1=Fl1(i)+fx1
         fx2=Fl2(i)+fx2
         WRITE(1,40) fx0, fx1, fx2
         enddo
         CLOSE(1)
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  40    FORMAT( E24.15E3, E24.15E3, E24.15E3 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RETURN
      END SUBROUTINE transfers

!*****************************************************************
      SUBROUTINE CFL_condition(cfl,ps,nu,nn,dt)
!-----------------------------------------------------------------
!        Parameters
!     cfl :cfl factor

      USE mpivars
      USE kes
      USE ali
      USE grid
      USE fft
      IMPLICIT NONE
      DOUBLE COMPLEX, DIMENSION(ny,ista:iend) :: ps
      DOUBLE COMPLEX, DIMENSION(ny,ista:iend) :: c1,c2,c3
      DOUBLE PRECISION, DIMENSION(nx,jsta:jend)    :: r1,r2,r3
      INTEGER :: i,j,nn,nnv
      DOUBLE PRECISION        :: tmp,dt,cfl,nu
      DOUBLE PRECISION        :: tmp1,tmp2,kx_max,ky_max,kt_max,nrm

      nrm = dble(nx)*dble(ny)

!!!!! FINDING MAX(|nabla psi|**2) 
      CALL derivk2(ps,c1,1)
      CALL derivk2(ps,c2,2)    
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      DO j = jsta,jend
         DO i = 1,nx
            r3(i,j) = r1(i,j)*r1(i,j)+r2(i,j)*r2(i,j)
         END DO
      END DO


      tmp=maxval(r3) !max energy density
      call MPI_REDUCE(tmp,tmp1,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tmp1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
!!!!!

      kx_max=(dble(nx)/3.0d0)*Qx
      ky_max=(dble(ny)/3.0d0)*Qy
      kt_max=sqrt(kx_max**2 + ky_max**2)
      tmp=sqrt(tmp1)/nrm  ! max horizontal velocity 
      dt = cfl/(kt_max*tmp+nu*(kt_max**(2*nn)))

 
      RETURN
      END SUBROUTINE CFL_condition


!NOT ANIS
!*****************************************************************
      SUBROUTINE test_sub(time,ps,nu,nn,dt)
!-----------------------------------------------------------------
      USE mpivars
      USE kes
      USE ali
      USE grid
      USE fft
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(ny,ista:iend) :: ps
      DOUBLE COMPLEX, DIMENSION(ny,ista:iend) :: c1,c2,c3,c4
      DOUBLE PRECISION, DIMENSION(nx,jsta:jend)    :: r1,r2,r3
      INTEGER :: i,j,nn
      DOUBLE PRECISION        :: tmp,nu,dt,div,cfl,time
      DOUBLE PRECISION        :: tmp1,tmp2,kx_max,ky_max,kt_max,nrm

      kx_max = dble(nx)/3.0d0
      ky_max = dble(ny)/3.0d0
      kt_max = sqrt(kx_max**2 + ky_max**2)
      nrm = dble(nx) * dble(ny)
      tmp = 1/nrm

      CALL derivk2(ps,c1,1)     !D_x psi = -vy
      CALL derivk2(ps,c2,2)     !D_y psi =  vx

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Divergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL derivk2(c2,c3,1)     !c3 = D_x vx
      CALL derivk2(-c1,c4,2)     !c2 = D_y vy
      CALL fftp2d_complex_to_real(plancr,c3,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c4,r2,MPI_COMM_WORLD)
      DO j = jsta,jend
         DO i = 1,nx
            tmp1 = (r1(i,j)+r2(i,j))*tmp
         END DO
      END DO
       
      CALL MPI_REDUCE(tmp1,div,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       CFL     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      DO j = jsta,jend
         DO i = 1,nx
                r3(i,j) = r1(i,j)*r1(i,j)+r2(i,j)*r2(i,j)
         END DO
      END DO

      tmp=maxval(r3)
      call MPI_REDUCE(tmp,tmp1,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tmp1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

      !kcut=(dble(n)/3.0d0)
      !tmp=sqrt(tmp1)/nrm
      cfl = dt*(kt_max*tmp+nu*(kt_max**(2*nn)))

!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='test_sub.txt',position='append')
         WRITE(1,20) time, div, cfl
   20 FORMAT(E23.14E3,E23.14E3,E23.14E3 )
         CLOSE(1)

      ENDIF

      RETURN
      END SUBROUTINE test_sub


!NOT ANIS
!*****************************************************************
      SUBROUTINE const_inj(ps,kdn,kup,fp0,fp,kin,seed1)
!-----------------------------------------------------------------
!       This subroutine assures that we inject constant energy.
!       It is called when iflow == 2
!       kin == 0 for vz
!           == 1 for ps

      USE mpivars
      USE kes
      USE ali
      USE grid
      USE fft
      USE var
      USE random
      IMPLICIT NONE
!                                               ps
      DOUBLE COMPLEX, DIMENSION(ny,ista:iend) :: ps,fp
      INTEGER :: i,j,seed1,kin
      DOUBLE PRECISION        :: tmp,kdn,kup,Efp,fp0
      DOUBLE PRECISION        :: tmp1,tmp2,tmp3,two,phase2d


!!!!!!!!!!! 2D Part!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Computes the energy at FORCING scale!!!!!!!!!!!!!!
      tmp = 1.0d0/dble(nx)**2/dble(ny)**2
      tmp1= 0.0d0
      tmp2= 0.0d0

      DO i = ista,iend
         IF (i.ne.1) THEN
         two = 2.0d0
         DO j = 1,ny
            IF ((kk2(j,i).le.kup*kup).and.(kk2(j,i).ge.kdn*kdn )) THEN
               fp(j,i)=  ps(j,i)/(abs(ps(j,i))**2+1.0d0)
            ELSE
               fp(j,i) = 0.0d0
            ENDIF
         END DO
         ELSEIF (i.eq.1) THEN
         fp(    1,i) = 0.0d0
         fp(ny/2+1,i) = 0.0d0
         DO j = 2,ny/2
            IF ((kk2(j,i).le.kup*kup).and.(kk2(j,i).ge.kdn*kdn ) ) THEN
               fp(    j,i) =  ps(j,i)/(abs(ps(j,i))**2+1.0d0)
               fp(ny-j+2,i) = conjg(fp(j,i))
            ELSE
               fp(    j,i) = 0.0d0
               fp(ny-j+2,i) = 0.0d0
            ENDIF
         END DO
         ENDIF
      END DO
        
      CALL inerprod(ps,fp,kin,Efp) ! Finds the dot product: either |nabla psi \nabla fpsi| or |vz fz| 
!!!!!Rescaling of forcing!!!
      seed1=myrank
      DO i = ista,iend
         IF (i.ne.1) THEN
         DO j = 1,ny
            IF ((kk2(j,i).le.kup*kup).and.(kk2(j,i).ge.kdn*kdn )) THEN
               tmp3    = randu(seed1)*sqrt(kk2(j,i)) 
               fp(j,i) = fp(j,i)*fp0/Efp + im*tmp3*ps(j,i)
            ELSE
               fp(j,i) = 0.0d0
            ENDIF
         END DO
         ELSEIF (i.eq.1) THEN
         fp(    1,i) = 0.0d0
         fp(ny/2+1,i) = 0.0d0
         DO j = 2,ny/2
            IF ((kk2(j,i).le.kup*kup).and.(kk2(j,i).ge.kdn*kdn )) THEN
               tmp3    = randu(seed1)*sqrt(kk2(j,i))
               fp(    j,i) = fp(j,i)*fp0/Efp + im*tmp3*ps(j,i)
               fp(ny-j+2,i) = conjg(fp(j,i))
            ELSE
               fp(    j,i) = 0.0
               fp(ny-j+2,i) = 0.0
            ENDIF
         END DO
         ENDIF
      END DO
      
      RETURN
      END SUBROUTINE const_inj

!*****************************************************************
      SUBROUTINE rand_force(kdn,kup,fp0,dt,seed,kin,fp)
!-----------------------------------------------------------------
!       This subroutine creates random forcing.
!       It is called when iflow == 3.
!       kin  == 1   for ps forcing
      USE var
      USE mpivars
      USE kes
      USE ali
      USE grid
      USE fft
      USE random
        IMPLICIT NONE
!                                              
      DOUBLE COMPLEX, DIMENSION(ny,ista:iend) ::fp
      INTEGER :: i,j,seed,kin
      DOUBLE PRECISION        :: tmp,kdn,kup,fp0,energyfp,energyfp2
      DOUBLE PRECISION        :: dt,tmp1,tmp2,two,phase,kxf,kyf,theta

!!!!!!!!!!! Like Chan et al. 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Choose random vector of length kup and random phase !!!
        IF (myrank.eq.0) THEN
                 theta = randu(seed)*pi
                 phase=randu(seed+1)*pi
                 kxf = floor(kup/Qx*cos(theta)+0.5)*Qx
                 kyf = floor(kup/Qy*sin(theta)+0.5)*Qy
!                 print*,"DBG (kx,ky)", kxf,kyf
         ENDIF
         CALL MPI_BCAST(kxf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(kyf,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(phase,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         tmp = sqrt(2*fp0)/(sqrt(dt)*kup**kin)
         tmp1 = tmp*cos(phase)
         tmp2 = tmp*sin(phase)
         
      DO i = ista,iend
         IF (i.ne.1) THEN
         DO j = 1,ny
            !IF (((kx(i).eq.kxf).or.(kx(i).eq.(-kxf))).and.(ky(j).eq.kyf).or.(ky(j).eq.(-kyf))) THEN
            IF ((kx(i).eq.kxf).and.(ky(j).eq.kyf)) THEN
            fp(j,i) = (tmp1+im*tmp2)*dble(nx)*dble(ny)!/2.0d0
            ELSE
               fp(j,i) = 0.0d0
            ENDIF
         END DO
         ELSEIF (i.eq.1) THEN
         fp(    1,i) = 0.0d0
         fp(ny/2+1,i) = 0.0d0
         DO j = 2,ny/2
            IF ((kx(i).eq.kxf).and.(ky(j).eq.kyf)) THEN
               fp(j,i) = (tmp1+im*tmp2)*dble(nx)*dble(ny)!/2.0d0
               fp(ny-j+2,i) = conjg(fp(j,i))!/2.0d0
            ELSE
               fp(j,i) = 0.0d0
               fp(ny-j+2,i) = 0.0d0
            ENDIF
         END DO
         ENDIF
      END DO

      RETURN
      END SUBROUTINE rand_force

!*****************************************************************
      SUBROUTINE rand_force_broad(kdn,kup,fp0,dt,seed,kin,fp)
!-----------------------------------------------------------------
!       This subroutine creates random forcing in arbitrary shell [kdn,kup]
!       It is called if iflow == 4.
!       kin == 1   for ps
      USE var
      USE mpivars
      USE kes
      USE ali
      USE grid
      USE fft
      USE random
      IMPLICIT NONE
!
      DOUBLE COMPLEX, DIMENSION(ny,ista:iend) :: fp
      INTEGER :: i,j,seed,seed1,kin
      DOUBLE PRECISION        :: tmp,kdn,kup,fp0,energyfp,energyfp2,kh
      DOUBLE PRECISION        :: dt,tmp1,tmp2,two,phase,theta,dump

      CALL MPI_BCAST(seed,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      seed1=seed+myrank
     
      DO i = ista,iend
         DO j = 1,ny
               kh = SQRT(kk2(j,i))
               IF ((kh.le.kup).and.(kh.ge.kdn)) THEN
                  dump = 1.0d0
                  !IF (i.eq.1) dump=0.0d0
                  phase = 2*pi*randu(seed1)
                  fp(j,i) = (COS(phase)+im*SIN(phase))*dump
               ELSE
                  fp(j,i) = 0.
               ENDIF
         END DO

         ! MAKE SURE THAT the kx = 0 modes satisfy f(-ky,kx=0) = conjug(f(ky,kx=0))
         IF (ista.eq.1) THEN
         DO j = 1,ny/2+1
            IF ((sqrt(kk2(j,1)).le.kup).and.(sqrt(kk2(j,1)).ge.kdn)) THEN
             seed1   = seed1 + 1
             phase   = 2*pi*randu(seed1)
             fp(j,i) = (COS(phase)+im*SIN(phase))*dump
           ELSE
             fp(j,1) = 0.0d0
            ENDIF
             fp(ny-j+2,1) = CONJG(fp(j,1))
         END DO
         ENDIF
      END DO

      CALL energy(fp,tmp,1)
      CALL MPI_BCAST(tmp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      fp = fp*sqrt(2*fp0)/sqrt(tmp*dt)
      seed=seed1
      CALL MPI_BCAST(seed,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      END SUBROUTINE rand_force_broad




!*****************************************************************
      SUBROUTINE yavg_ps(ps,C,R,ext,odir)
!-----------------------------------------------------------------
!
! Computes the energy power spectrum in 2D.
! The output is written to a file by the first node.
!
! Parameters
!     a  : streamfunction
!     b  : vector potential
!     ext: the extension used when writting the file
!     kin: =1 computes the kinetic spectrum
!          =0 computes the magnetic spectrum
!
      USE kes
      USE grid
      USE mpivars
      USE fft
      USE grid
      USE ali
      
      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(nx)            :: psyavg
      DOUBLE PRECISION, DIMENSION(nmax/2+1)      :: Ek,Ektot1
      DOUBLE COMPLEX, DIMENSION(ny,ista:iend)    :: ps,C
      DOUBLE PRECISION, DIMENSION(nx,jsta:jend)  :: R
      DOUBLE PRECISION        :: tmp,two,tmp1,tmp2,tmp3
      INTEGER     :: kin
      INTEGER     :: kmn
      INTEGER     :: i,j
      CHARACTER*4 :: ext
      CHARACTER*100 :: odir


     CALL derivk2(ps,C,0)

      DO j=1,ny
       DO i=ista,iend
         IF (ky(j).gt.(1.0d0-tiny)) THEN
           C(j,i) = 0.0d0
         ENDIF
       ENDDO
      ENDDO

      tmp=1.0d0/dble(nx)/dble(ny)
      DO i = ista,iend
       DO j = 1,ny
          C(j,i) = C(j,i)*tmp
       END DO
      END DO
      CALL fftp2d_complex_to_real(plancr,C,R,MPI_COMM_WORLD)
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(odir) // '/yavg_ps.' // ext // '.txt')
!         WRITE(1,30) Q,0.0d0,0.5d0*tmp1
         do i=1,nx
         WRITE(1,55) R(i,1)
         enddo
         CLOSE(1)
      ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  55    FORMAT( E24.15E3,E24.15E3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      RETURN
      END SUBROUTINE yavg_ps
