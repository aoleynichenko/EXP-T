        subroutine spectrum(secteur,irrep,mainms,eone) 
        implicit none             
        character (len=80) :: ligne
        character(len=4)   :: secteur
        character(len=40)  :: fmt
        integer*4          :: irrep,mainms,stat,i,j,idummy
        complex*16         :: eone(*)
        real*8             :: dummy,au2cm
        real*8,allocatable,dimension(:)           :: es,work
        integer*4,allocatable,dimension(:)        :: irreps,numeros,
     &                                               iwork,innerno
        character(len=:),allocatable,dimension(:) :: secteurs,cwork
        
        au2cm=219474.631280634d0
        fmt='(1x,i4,f18.12,f15.3,4x,a4,a3,i3,a2,i4)'
        
        open(unit=17,form='formatted',file='FULLSPECTRUM')
c determiner le nombre de lignes
        i=0
        do while (1 .eq. 1)
            read(17,'(a)', iostat=stat) ligne
            if (stat .ne. 0) exit
            i=i+1
        end do
c allocate
        allocate(character(4)::secteurs(i+mainms))
        allocate(character(4)::cwork(i+mainms))
        allocate(es(i+mainms),irreps(i+mainms),innerno(i+mainms),
     &  work(i+mainms),iwork(i+mainms),numeros(i+mainms))        
        rewind 17
        if(i .ne. 0) then
            do j=1,i
              read(17,fmt,iostat=stat)
     &            idummy,es(j),dummy,secteurs(j),ligne(1:3),
     &            innerno(j),ligne(1:2),irreps(j)
            enddo
        endif
        do j=1,mainms
          innerno(i+j)=j
          es(i+j)=real(eone(j))
          irreps(i+j)=irrep
          secteurs(i+j)=secteur
          if( secteur .eq. '0h0p') then
            es(i+j)=0.d0
          endif
        enddo        
c
c   rearrange
c        
        call asc(es,numeros,i+mainms)
        do j =1,i+mainms
          iwork(j)=irreps(numeros(j))
          cwork(j)=secteurs(numeros(j))
        enddo
        do j=1,i+mainms
          irreps(j)=iwork(j)
          secteurs(j)=cwork(j)
        enddo
        do j =1,i+mainms
          iwork(j)=innerno(numeros(j))
        enddo
        do j=1,i+mainms
          innerno(j)=iwork(j)
        enddo
        
        rewind 17
        do j=1,i+mainms
        write(17,fmt)j, es(j),(es(j)-es(1))*au2cm,
     &        secteurs(j),'  (',innerno(j),')',irreps(j)
        enddo
        close(unit=17,status='keep')              
        deallocate(secteurs,cwork,es,irreps,numeros,work,iwork,innerno)      
        write (*,'(2x,a,i4,a,i4,a)') '  FULLSPECTRUM:',
     &        i,' lines read and',i+mainms,' written' 
        return
        end
************************************************************************
*                                                                      *
*                             zinver                                   *
*                                                                      *
************************************************************************
      subroutine zinver(A,M)
c
c  wrapper for ZGETRF & ZGETRI - matrix inversion
c
      implicit real*8 (a-h,o-z)
      complex*16 A(M,M)
      complex*16,allocatable::WORK(:)
      integer,allocatable::IPIV(:)
      integer error
       allocate(WORK(M),IPIV(M),stat=error)
       if (error .ne. 0)then
          write(*,*) '    !!! zinver: allocation error'
          stop
       end if

c      lu factorization

       call ZGETRF(M,M,A,M,IPIV,info)
       if(info .ne. 0) then
         write(*,*) '     !!! zinver: ZGETRF died'
         stop
       end if

       call ZGETRI(M,A,M,IPIV,WORK,M,info)
       if(info .ne. 0) then
         write(*,*)  '      !!! zinver: ZGETRI died'
         stop
       end if
       deallocate(IPIV,WORK)
       return
       end
*      subroutine primam(a,b,c,d,m,n)
*      implicit real*8(a-h,o-z)
*      character*80 ligne
*      character*1 c
*      dimension a(m,n),b(n)
*c
*c  
*      nvec=n
*c
*   12 format(2x,a1,f4.1,' ',f7.1,8f9.1)
*   14 format(f6.1,' ',8f9.5)
*   13 format(2x,77('-'))
*      write(6, 13)
*      j1=1
*      j2=8
*      if(nvec.le.8)j2=nvec
*   17 continue
*      write(6, 12)c,d,(b(j),j=j1,j2)
*      do 4 i=1,n
*      write(6,14)b(i),(a(i,j),j=j1,j2)
*    4 continue
*      if(j2.eq.nvec)go to 99
*      j1=j1+8
*      j2=j2+8
*      if(j2.gt.nvec)j2=nvec
*      go to 17
*   99 write(6, 13)
*      return
*      end
************************************************************************
*                                                                      *
*                              eiggen / zgeev, standard                *
*                                                                      *
************************************************************************
      subroutine eiggen(a,ag,ad,eva,max)
c
c   wrapper for lapack's zgeev
c
      implicit real*8 (a-h,o-z)
c
c   arguments
c      
      complex*16 eva(*)
      complex*16 a(max,max),ag(max,max),ad(max,max)
c
c   internal
c      
      complex*16, allocatable :: cwork(:)
      complex*16, allocatable :: work(:) 
      real*8, allocatable :: rwork(:)
      integer, allocatable :: iwork(:)
      n=max
      mfound=n
      lwork=4*max
c
c     allocation des internes
c
c      
      allocate(work(lwork))
      allocate(rwork(2*max))
      allocate(iwork(2*max))
      allocate(cwork(max))
c      allocate(ifail(n))
c      write(6,*) max
      call zgeev( 'V', 'V', 
c      compute both left and right eigenvectors
     &      max, a, max,
c      matrix size, input matrix, matrix size       
     &   eva, 
c      eigenvalues
     &   ag,max,ad,max,  
c      ev gauches,  ev droits        
     &  work, lwork,
c       complex working array, better of size lwork=4*max 
     &  rwork,
c       real working array 2*max      
     &  info )
      if(info .ne. 0) then
      write(*,*)'   ----- zgeev crashed. info:', info,'  --------' 
      stop
      endif
c
c    arrange: energies (reelles) croissantes
c     
      call ascrea(eva,iwork,max)
c      
      do iii =1,max
        do jjj=1,max
           cwork(jjj)=ad(iii,iwork(jjj))
        enddo
        do jjj=1,max
           ad(iii,jjj)=cwork(jjj)
        enddo
        do jjj=1,max
           cwork(jjj)=ag(iii,iwork(jjj))
        enddo
        do jjj=1,max
           ag(iii,jjj)=cwork(jjj)
        enddo
      enddo
c
c     normalize to get biorthonormal sets
c     right vectors are normalized to 1
c
      do iii=1,max
        cwork(1)=(0.d0,0.d0)
        do jjj=1,max
        cwork(1)=cwork(1)+conjg(ad(jjj,iii))*ad(jjj,iii)
        enddo
c
        do jjj=1,max
        ad(jjj,iii)=ad(jjj,iii)/cwork(1)
        enddo
      enddo 
c           
      do iii=1,max
        cwork(1)=(0.d0,0.d0)
         do jjj=1,max
         cwork(1)=cwork(1)+conjg(ag(jjj,iii))*ad(jjj,iii)
         enddo
c
c     left vectors are normalized to left-right product = 1
c
         do jjj=1,max
         ag(jjj,iii)=ag(jjj,iii)/conjg(cwork(1))
         enddo
      enddo      
c      
      deallocate(cwork,work,rwork,iwork)
      return
      end
************************************************************************
*                                                                      *
*                              asc                                     *
*                                                                      *
************************************************************************
      subroutine asc(x,n,maxmax)
      implicit real*8(a-h,o-z)
c      complex*16 x(*),a
      dimension n(*),x(*)
      m=maxmax
c      
      do  i=1,m
        n(i)=i
      enddo
c    
      l=m-1
      if (l.lt.1) return
      do i=1,l
        ip1=i+1
        do  j=ip1,m
          if( x(i) .le. x(j) ) cycle
          a=x(i)
          x(i)=x(j)
          x(j)=a
          k=n(i)
          n(i)=n(j)
          n(j)=k
        enddo
      enddo   
      return
      end
************************************************************************
*                                                                      *
*                              ascrea                                  *
*                                                                      *
************************************************************************
      subroutine ascrea(x,n,maxmax)
      implicit real*8(a-h,o-z)
      complex*16 x(*),a
      dimension n(*)
      m=maxmax
c      
      do  i=1,m
        n(i)=i
      enddo
c    
      l=m-1
      if (l.lt.1) return
      do i=1,l
        ip1=i+1
        do j=ip1,m
          if( real(x(i)) .le. real(x(j)) ) cycle
          a=x(i)
          x(i)=x(j)
          x(j)=a
          k=n(i)
          n(i)=n(j)
          n(j)=k
        enddo
      enddo
      return
      end
************************************************************************
*                                                                      *
*                            primac                                    *
*                                                                      *
************************************************************************
      subroutine primac(ca,m,n,nvec,ex)
      implicit real*8(a-h,o-z)
      complex*16 ca(m,nvec)
      character*128 ligne
      dimension a(10),aa(5)
c
c  prints rectangular complex matrix n*nvec
c
   12 format(6x,'   ',i9,4i18)
   14 format((i6,'  ',5(f9.5,f8.5,1x)))
   15 format(2x,77(' '))
c   13 format(2x,77(' '))
      write(6,15)
      j1=1
      j2=5
      if(nvec.le.5)j2=nvec
   17 continue
      write(6, 12)(j,j=j1,j2)
      do 4 i=1,n
      s=0.d+00
      irtr=0
          do j=j1,j2
          car=real(ca(i,j))
          irtr=irtr+1
          a(irtr)=car
          cai=dimag(ca(i,j))
          irtr=irtr+1
          a(irtr)=cai
          dabaij=dabs(car)+dabs(cai)
          iaa=irtr/2
          aa(iaa)=dabaij
          if( dabaij.gt.s ) s=dabaij
        enddo
c
        ligne='                                        '
c        ligne(41:80)='                                        '
        write(ligne, 14)i,(a(i9),i9=1,irtr)
c
      irtr2=irtr/2
      lstar=9
      do i9=1,irtr2
      if(aa(i9).lt.ex)then
        ligne(lstar:lstar+17)=
     &                        '   .       .      '
      endif
      lstar=lstar+18
      enddo
      lstar=lstar-1
      write(6,'(a)') ligne(1:lstar)
c
c
    4 continue
      if(j2.eq.nvec)go to 99
      j1=j1+5
      j2=j2+5
      if(j2.gt.nvec)j2=nvec
      go to 17
   99 write(6, 15)
      return
      end
************************************************************************
*                                                                      *
*                            primaclong                                *
*                                                                      *
************************************************************************
      subroutine primaclong(ca,m,n,nvec,ex)
      implicit real*8(a-h,o-z)
      complex*16 ca(m,nvec)
      character*128 ligne
      dimension a(10),aa(5)
c
c  prints rectangular complex matrix n*nvec
c
   12 format(6x,'   ',i14,3i23)
   14 format((i6,'  ',4(f12.2,f9.2,'i',1x)))
   15 format(2x,77(' '))
   13 format(2x,77(' '))
      write(6,15)
      j1=1
      j2=4
      if(nvec.le.4)j2=nvec
   17 continue
      write(6, 12)(j,j=j1,j2)
      do 4 i=1,n
      s=0.d+00
      irtr=0
          do j=j1,j2
          car=real(ca(i,j))
          irtr=irtr+1
          a(irtr)=car
          cai=dimag(ca(i,j))
          irtr=irtr+1
      a(irtr)=cai
      dabaij=dabs(car)+dabs(cai)
          iaa=irtr/2
          aa(iaa)=dabaij
      if( dabaij.gt.s ) s=dabaij
      enddo
c   if( s.lt.ex ) goto 4
c
        ligne=' '
c        ligne(41:80)='                                        '
        write(ligne, 14)i,(a(i9),i9=1,irtr)
c
      irtr2=irtr/2
      lstar=9
      do i9=1,irtr2
      if(aa(i9).lt.ex)then
        ligne(lstar:lstar+22)=
     &                        '         .        .  '
      endif
      lstar=lstar+23
      enddo
      lstar=lstar-1
      write(6,'(a)') ligne(1:lstar)
c
c
    4 continue
      if(j2.eq.nvec)go to 99
      j1=j1+4
      j2=j2+4
      if(j2.gt.nvec)j2=nvec
      go to 17
   99 write(6, 15)
      return
      end
************************************************************************
*                                                                      *
*                            loeworth                                  *
*                                                                      *
************************************************************************
      subroutine loeworth(zvecnon,zvecorth,m,n,label,accu,iprint)
*
*     zvecnon  matrice de vecteurs non orthogonaux,       m x n
*              survit    
*            
*     zvecorth vecteurs orthogonalises,                   m x n
*     
*              si label(1:1)='?', zvecorth peut etre =zvecnon
*                               et seules les valeurs propres de
*                               la matrice de recouvrement sont calcules
*     m        longeur de vecteurs
*     n        nombre de vecteurs
*
      implicit complex*16(z)
      dimension zvecnon(m,*),zvecorth(m,*)
      character*16 label
      real*8 accu
      allocatable :: zs(:,:),zs12(:,:),zwork(:,:)
      real*8, allocatable :: eva(:)
*
      allocate( zs(n,n),zs12(n,n),zwork(m,n) )
      allocate(eva(n))
*
*     evaluer le recouvrement, mettre dans zs
*
      call scapro(zvecnon,zvecnon,zs,m,n)
      if ( iprint.gt.2) then
         write(6,'(a)') '    loewdin orthogonalization: overlap matrix'
         call  primac(zs,n,n,n,1.d-5)
      endif
c      call eiggen(zs,zsl,zsr,zse,n)
c
      call eigher(zs,eva,n,n)
c
c      do i=1,n
c      zse(i)=cmplx(eva(i),0.d+00)
c      do j=1,n
c      zsr(i,j)=zs(i,j)
c      enddo
c      enddo


      if ( iprint.ge.1 .or. label(1:1).eq.'?'
     &      .or. eva(1).le.accu )
     &  write(6,'(5x,2a,f12.8/)')label(2:16),
     &   ' smallest overlap eigenvalue ', eva(1)
      if ( iprint.ge.2) then
        write(6,'(a)')'     overlap matrix eigenvalues '
        write(6,'((3(i6,f15.8)))') (iii,eva(iii),iii=1,n)
        write(6,*)' '
      endif
*
      if ( label(1:1).eq.'?' ) goto 13
*
      do i=1,n
         do j=1,n
            zzz=(0.d0,0.d0)
            do k=1,n
               zzz=zzz+zs(i,k)*conjg(zs(j,k))
     &                 /sqrt( eva(k) )
            enddo
            zs12(i,j)=zzz
         enddo
      enddo   
      if ( iprint.gt.2) then
         write(6,'(2a)') '     loewdin orthogonalization: ',
     &                'overlap matrix eigenvectors   '
         call  primac(zs,n,n,n,1.d-5)
         write(6,'(a)')'     loewdin orthogonalization: overlap^(-1/2)'
         call  primac(zs12,n,n,n,1.d-5)
      endif   

c
c     orthogonalize
c
      do i=1,n
         do j=1,m
            zzz=(0.d0,0.d0)
            do k=1,n
               zzz=zzz+zs12(k,i)*zvecnon(j,k)
            enddo
            zvecorth(j,i)=zzz
         enddo
      enddo   
c
c         call scapro(zvecorth,zvecorth,zs,m,n)
c         write(6,'(a)') '   orthogonality check   '
c         call  primac(zs,n,n,n,1.d-5)
c         call  primac(zvecnon,m,m,n,1.d-5)
c         call  primac(zvecorth,m,m,n,1.d-5)
c         
13    continue
      deallocate( zs,zs12,zwork,eva )
      return
      end
************************************************************************
*                                                                      *
*                             projinto                                 *
*                                                                      *
************************************************************************
      subroutine projinto(zctarg,zcmod,zenergies,zheffectif,ground,
     &  m,n,y_socoupled,iprint)
*      
*     projection sur l espace dit modele (vecteurs de base dans zcmod)
*     des vecteurs dit cible (dans zctarg), probablement non-orthogonaux
* 
*     if (y_socoupled ) then
*         modele, vecteurs so, target, so-free
*     else
*         
*     endif 
*      
      implicit complex*16(z)
      real*8 ground
      logical y_socoupled
      dimension zctarg(m,*),zcmod(m,*),zenergies(*),
     & zheffectif(n,*)
      allocatable :: zs(:,:),zwork(:,:),zorth(:,:)
c      
c
      if ( iprint.ge.2) then
         write(6,'(/,a)') '    model space orthogonal basis vectors'
         call  primac(zcmod,m,m,n,1.d-5)
         write(6,'(a)')   '    so-coupled eigenvalues '
         write(6,'((2x,5(i4,f13.7)))')
     &        (i,real(zenergies(i)),i=1,n)
         write(6,'(a)')   ' '
         write(6,'(a)')   '    target non-orthogonal vectors '
         call  primac(zctarg,m,m,n,1.d-5)
         write(6,'(a)')   ' '
      endif
c
      allocate( zs(n,n),zwork(m,n),zorth(m,n) )
c      
      do i=1,n
         do j=1,n
           zx=(0.d+0,0.d+0)
           do k=1,m
             zx=zx+conjg(zcmod(k,i))*zctarg(k,j)
           enddo
           zs(i,j)=zx
         enddo
      enddo
c 
      if ( iprint.ge.2) then
         write(6,'(a)') '   overlap target - model '
         call  primac(zs,n,n,n,1.d-5)
      endif
c 
      do i=1,m
         do j=1,n
            zx=(0.d+0,0.d+0)
            do k=1,n
               zx=zx+zcmod(i,k)*zs(k,j)
            enddo
            zwork(i,j)=zx
         enddo
      enddo
c      
      if ( iprint.ge.2) then
         write(6,'(a)') '    projected target vectors '
         call  primac(zwork,m,m,n,1.d-5)
      endif
*
*     orthogonaliser
*
      call loeworth(zwork,zorth,m,n,
     &' target2model   ',1.5d+0,iprint)
* 
c      
      if ( iprint.ge.2) then
         write(6,'(a)') '   projected & orthogonalized target vectors'
         call  primac(zorth,m,m,n,1.d-5)
      endif
*         call scapro(zorth,zorth,zs,m,n)
c         write(6,'(a)') '   orthogonality check   '
c         call  primac(zs,n,n,n,1.d-5)
c
c     construire les matrices de transformation
c
      if( y_socoupled ) then
*
*      
        call scapro(zcmod,zorth,zs,m,n)
*        
      else
*      
        call scapro(zorth,zcmod,zs,m,n)
*        
      endif
*
*
      if ( iprint.ge.2 ) then
        write(6,'(a)') '   transformation   '
        call  primac(zs,n,n,n,1.d-5)
      endif
      do i=1,n
         do j=1,n
            zx=(0.d+0,0.d+0)
            do k=1,n
               zx=zx+conjg(zs(k,i))*zs(k,j)*zenergies(k)
            enddo
            zheffectif(i,j)=zx
         enddo
      enddo
c      if ( iprint.ge.1 ) then
         if(y_socoupled) then
            write(6,'(2a)') '     effective hamiltonian ',
     &     'in the subspace spanned by so-coupled eigenvectors'
         else 
            write(6,'(2a)') '     effective hamiltonian ',
     &     'in the subspace spanned by so-free eigenvectors'
         endif        
         call  primac(zheffectif,n,n,n,1.d-5)
         write(6,'(1x,2a,f16.2)')
     &    '        diagonal, au      diagonal, cm-1  '
         do iki=1,n
         if(abs(ground-666.d+0) .le. 1.d-5) then
            write(6,'(i6,f17.7,f17.2)')iki,
     &       real(zheffectif(iki,iki)),
     &       (real(zheffectif(iki,iki))-real(zheffectif(1,1))) 
     &      *219474.631280634d+0
         else
            write(6,'(i6,f17.7,f17.2)')iki,
     &       real(zheffectif(iki,iki)),
     &       (real(zheffectif(iki,iki))-ground) 
     &      *219474.631280634d+0
         endif
         enddo
c      endif   
c
c!      call eiggen(zheffectif,zs,zwork,zenergies,n)
c!      do i=1,n
c!      write(*,'(i5,F15.8,f15.2)')i,real(zenergies(i)),
c!     &   (real(zenergies(i))-ground)*219474.631280634d+0
c!      enddo
c!      stop
c
      deallocate( zs, zwork, zorth )
*
      return
      end
*
*
      subroutine scapro(za,zb,zc,m,n)
*      
*  zc(square)=za^\dag * zb (both rectangular, same size)      
*
      implicit complex*16(z)
      dimension za(m,*),zb(m,*),zc(n,*)
      do i=1,n
         do j=1,n
           zx=(0.d+0,0.d+0)
           do k=1,m
             zx=zx+conjg(za(k,i))*zb(k,j)
           enddo
           zc(i,j)=zx
         enddo
      enddo   
      return
      end
c--------------------------------------------------------
      subroutine eigher(a,eva,max,n)
c      
c eigenvalues eva /eigenvectors a  of complex hermitian matrix a 
c wrapper for zheev
c a is destructed
c eigenvalues in ascending order
c      
      implicit real*8 (a-h,o-z)
      dimension eva(*)
      complex*16 a(max,max)
c      ,zstore 
      complex*16, allocatable :: work(:) 
      real*8, allocatable :: rwork(:)
      integer, allocatable :: iwork(:)
c      integer, allocatable :: ifail(:)
      mfound=n
      lwork=4*n
c      allocate(a(n,n))
      allocate(work(lwork))
      allocate(rwork(7*n))
      allocate(iwork(5*n))
c      allocate(ifail(n))
      call zheev( 'V',
     & 'L', n, a, max, 
     &   eva, 
     &  work, lwork, rwork, 
     &  info )
      if(info .ne. 0) then
      write(*,*)'   ----- zheev crashed. info:', info,'  --------' 
      stop
      endif
      deallocate(work,rwork,iwork)      
      return
      end
     
