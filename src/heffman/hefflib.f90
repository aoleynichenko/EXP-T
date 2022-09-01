        subroutine spectrum(secteur,irrep,mainms,eone) 
        use nombres_mod
!
        implicit none             
        character (len=80) :: ligne
        character(len=4)   :: secteur
        character(len=40)  :: fmt='(1x,i4,f18.12,f15.3,4x,a4,a3,i3,a2,i4)'
        integer            :: irrep,mainms,stat,i,j,idummy
        complex(kind=pd)         :: eone(*)
        real(kind=pd)             :: dummy
        real(kind=pd),allocatable,dimension(:)  :: es,work
        integer,allocatable,dimension(:)        :: irreps,numeros, iwork,innerno
        character(len=:),allocatable,dimension(:) :: secteurs,cwork

!  
        
        open(unit=17,form='formatted',file='FULLSPECTRUM')
! compter les lignes
        i=0
        stat=0
        do while (stat == 0)
            read(17,'(a)', iostat=stat) ligne
            if (stat == 0) i=i+1
        end do
        rewind 17
! allocate
        allocate(character(4)::secteurs(i+mainms))
        allocate(character(4)::cwork(i+mainms))
        allocate(es(i+mainms),irreps(i+mainms),innerno(i+mainms),                   &
                 work(i+mainms),iwork(i+mainms),numeros(i+mainms))        
        if(i /= 0) then
            do j=1,i
              read(17,fmt,iostat=stat) idummy,es(j),dummy,secteurs(j),ligne(1:3),   &
                                       innerno(j),ligne(1:2),irreps(j)
            enddo
        endif
        do j=1,mainms
          innerno(i+j)=j
          es(i+j)=real(eone(j),kind=pd)
          irreps(i+j)=irrep
          secteurs(i+j)=secteur
          if( secteur .eq. '0h0p') then
            es(i+j)=0.d0
          endif
        enddo        
!
!   rearrange
!        
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
        write(17,fmt)j, es(j),(es(j)-es(1))*au2cm, secteurs(j),'  (',innerno(j),')',irreps(j)
        enddo
        close(unit=17,status='keep')              
        deallocate(secteurs,cwork,es,irreps,numeros,work,iwork,innerno)      
        write (*,'(2x,a,i4,a,i4,a)') '  FULLSPECTRUM:', i,' lines read and',i+mainms,' written' 
        return
        end
!***********************************************************************
!                                                                      *
!                             zinver                                   *
!                                                                      *
!***********************************************************************
      subroutine zinver(A,M)
!
!  wrapper for ZGETRF & ZGETRI - matrix inversion
!
      use nombres_mod
      implicit none
!      
      integer, intent(in)                :: M
      complex(kind=pd),intent(in)        :: A(M,M)
      complex(kind=pd),allocatable       :: WORK(:)
      integer,allocatable                :: IPIV(:)
      integer error,info
       allocate(WORK(M),IPIV(M),stat=error)
       if (error /= 0)then
          write(*,*) '    !!! zinver: allocation error'
          stop
       end if

!      lu factorization

       call ZGETRF(M,M,A,M,IPIV,info)
       if(info /= 0) then
         write(*,*) '     !!! zinver: ZGETRF died'
         stop
       end if

       call ZGETRI(M,A,M,IPIV,WORK,M,info)
       if(info /= 0) then
         write(*,*)  '      !!! zinver: ZGETRI died'
         stop
       end if
       deallocate(IPIV,WORK)
       return
       end
!***********************************************************************
!                                                                      *
!                              eiggen / zgeev, standard                *
!                                                                      *
!***********************************************************************
      subroutine eiggen(a,ag,ad,eva,max)
!
!   wrapper for lapack zgeev
!
      use nombres_mod
      implicit none
!
!   arguments
!     
      integer,intent(in)                           :: max 
      complex(kind=pd), dimension(*)               :: eva
      complex(kind=pd), dimension(max,max)         :: a, ag, ad
!
!   internal
!      
      complex(kind=pd), allocatable :: cwork(:)
      complex(kind=pd), allocatable :: work(:) 
      real(kind=pd), allocatable :: rwork(:)
      integer, allocatable :: iwork(:)
      integer :: n, mfound, lwork, info, iii,jjj
      n=max
      mfound=n
      lwork=4*max
!
!     allocation des internes
!
!      
      allocate(work(lwork))
      allocate(rwork(2*max))
      allocate(iwork(2*max))
      allocate(cwork(max))
      call zgeev( 'V', 'V',        & !  compute both left and right eigenvectors
                   max, a, max,    & !  matrix size, input matrix, matrix size       
                   eva,            & !  eigenvalues
                   ag,max,ad,max,  & !  ev gauches,  ev droits        
                   work, lwork,    & !  complex working array, better of size lwork=4*max 
                   rwork,          & !  real working array 2*max      
                   info )
      if(info /= 0) then
      write(*,*)'   ----- zgeev crashed. info:', info,'  --------' 
      stop
      endif
!
!    arrange: energies (reelles) croissantes
!     
      call ascrea(eva,iwork,max)
!      
      do iii =1,max
        do jjj=1,max
           cwork(jjj)=ad(iii,iwork(jjj))
        enddo
        ad(iii,1:max)=cwork(1:max)
        do jjj=1,max
           cwork(jjj)=ag(iii,iwork(jjj))
        enddo
        ag(iii,1:max)=cwork(1:max)
      enddo
!
!     normalize to get biorthonormal sets
!     right vectors are normalized to 1
!
      do iii=1,max
        cwork(1)=(0.d0,0.d0)
        do jjj=1,max
        cwork(1)=cwork(1)+conjg(ad(jjj,iii))*ad(jjj,iii)
        enddo
        ad(1:max,iii)=ad(1:max,iii)/cwork(1)
      enddo 
!           
      do iii=1,max
        cwork(1)=(0.d0,0.d0)
        do jjj=1,max
          cwork(1)=cwork(1)+conjg(ag(jjj,iii))*ad(jjj,iii)
        enddo
!
!     left vectors are normalized to left-right product = 1
!
        ag(1:max,iii)=ag(1:max,iii)/conjg(cwork(1))
      enddo      
!      
      deallocate(cwork,work,rwork,iwork)
      return
      end
!***********************************************************************
!                                                                      *
!                              asc                                     *
!                                                                      *
!***********************************************************************
      subroutine asc(x,n,maxmax)
      use nombres_mod
      implicit none
      real(kind=pd),dimension(*),intent(inout) :: x
      integer,dimension(*),intent(out)         :: n
      integer, intent(in)                      :: maxmax
      real(kind=pd)                            :: a
      integer                                  :: m,l,i,j,k
      m=maxmax
!      
      do  i=1,m
        n(i)=i
      enddo
!    
      l=m-1
      if (l.lt.1) return
      do i=1,l
!        ip1=i+1
        do  j=i+1,m
          if( x(i) <= x(j) ) cycle
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
!***********************************************************************
!                                                                      *
!                              ascrea                                  *
!                                                                      *
!***********************************************************************
      subroutine ascrea(x,n,maxmax)
      use nombres_mod
      implicit none
      complex(kind=pd),dimension(*),intent(inout) :: x
      complex(kind=pd)                            :: a
      integer,dimension(*),intent(out)            :: n
      integer, intent(in)                      :: maxmax
      integer                                  :: m,l,i,j,k
      m=maxmax
!      
      do  i=1,m
        n(i)=i
      enddo
!    
      l=m-1
      if (l.lt.1) return
      do i=1,l
        do j=i+1,m
          if( real(x(i),kind=pd) <= real(x(j),kind=pd) ) cycle
          a=x(i)
          x(i)=x(j)
          x(j)=a
          k=n(i)
          n(i)=n(j)
          n(j)=k
        enddo
      enddo
      return
      end subroutine ascrea
!***********************************************************************
!                                                                      *
!                            loeworth                                  *
!                                                                      *
!***********************************************************************
      subroutine loeworth(zvecnon,zvecorth,m,n,label,accu,iprint)
!
!     zvecnon  matrice de vecteurs non orthogonaux,       m x n
!              survit    
!            
!     zvecorth vecteurs orthogonalises,                   m x n
!     
!              si label(1:1)='?', zvecorth peut etre =zvecnon
!                               et seules les valeurs propres de
!                               la matrice de recouvrement sont calcules
!     m        longeur de vecteurs
!     n        nombre de vecteurs
!

      use nombres_mod
      implicit none

!      implicit complex*16(z)
      complex(kind=pd), intent(in)  :: zvecnon(m,*)
      complex(kind=pd), intent(out) ::zvecorth(m,*)
      character(len=16),intent(in)  :: label
      real(kind=pd),    intent(in)  :: accu
      integer,          intent(in)  :: m,n,iprint 
      complex(kind=pd),allocatable  :: zs(:,:),zs12(:,:),zwork(:,:)
      real(kind=pd),   allocatable  :: eva(:)
      complex(kind=pd)              :: zzz
      integer                       :: i,j,k
!
      allocate( zs(n,n),zs12(n,n),zwork(m,n) )
      allocate(eva(n))
!
!     evaluer le recouvrement, mettre dans zs
!
      call scapro(zvecnon,zvecnon,zs,m,n)
      if ( iprint > 2) then
         write(6,'(a)') '    loewdin orthogonalization: overlap matrix'
         call imprec(zs,n,n,5)
      endif
      call eigher(zs,eva,n,n)
!
!
      if ( iprint > 1 .or. label(1:1).eq.'?' .or. eva(1) <= accu )         &
        write(*,'(5x,2a,f12.8/)')label(2:16),' smallest overlap eigenvalue ', eva(1)
      if ( iprint >= 2) then
        write(*,'(a)')'     overlap matrix eigenvalues '
        write(*,'((3(i6,f15.8)))') (i,eva(i),i=1,n)
        write(*,*)' '
      endif
!
      if ( label(1:1).eq.'?' ) goto 13
!
      do i=1,n
         do j=1,n
            zzz=(0.d0,0.d0)
            do k=1,n
               zzz=zzz+zs(i,k)*conjg(zs(j,k))/sqrt( eva(k) )
            enddo
            zs12(i,j)=zzz
         enddo
      enddo   
      if ( iprint > 2 ) then
         write(6,'(a)')'     loewdin orthogonalization: overlap matrix eigenvectors '
         call  imprec(zs,n,n,5) 
         write(6,'(a)')'     loewdin orthogonalization: overlap^(-1/2)'
         call  imprec(zs12,n,n,5)
      endif   

!
!     orthogonalize
!
      do i=1,n
         do j=1,m
            zzz=(0.d0,0.d0)
            do k=1,n
               zzz=zzz+zs12(k,i)*zvecnon(j,k)
            enddo
            zvecorth(j,i)=zzz
         enddo
      enddo   
13    continue
      deallocate( zs,zs12,zwork,eva )
      return
      end
!***********************************************************************
!                                                                      *
!                             projinto                                 *
!                                                                      *
!***********************************************************************
      subroutine projinto(zctarg,zcmod,zenergies,zheffectif,ground,m,n,y_socoupled,iprint)
!      
!     projection sur l espace dit modele (vecteurs de base dans zcmod)
!     des vecteurs dit cible (dans zctarg), probablement non-orthogonaux
! 
!     if (y_socoupled ) then
!         modele, vecteurs so, target, so-free
!     else
!         
!     endif 
!      
      use nombres_mod
      implicit none
!      implicit complex*16(z)
      real(kind=pd),intent(in)       :: ground
      logical      ,intent(in)       :: y_socoupled
      integer      ,intent(in)       :: m, n,iprint
      complex(kind=pd),intent(inout) :: zctarg(m,*),zcmod(m,*)
      complex(kind=pd),intent(in)    :: zenergies(*)
      complex(kind=pd),intent(out)   :: zheffectif(n,*)
      complex(kind=pd), allocatable  :: zs(:,:),zwork(:,:),zorth(:,:)
      integer                        :: i,j,k
      complex(kind=pd)               :: zx
!      
!
      if ( iprint > 2) then
         write(*,'(/,a)') '    model space orthogonal basis vectors'
         call  imprec(zcmod,m,n,5)
         write(*,'(a)')   '    so-coupled eigenvalues '
         write(*,'((2x,5(i4,f13.7)))') (i,real(zenergies(i),kind=pd),i=1,n)
         write(*,'(a)')   ' '
         write(*,'(a)')   '    target non-orthogonal vectors '
         call  imprec(zctarg,m,n,5)
         write(6,'(a)')   ' '
      endif
!
      allocate( zs(n,n),zwork(m,n),zorth(m,n) )
!      
      do i=1,n
         do j=1,n
           zx=(0.d+0,0.d+0)
           do k=1,m
             zx=zx+conjg(zcmod(k,i))*zctarg(k,j)
           enddo
           zs(i,j)=zx
         enddo
      enddo
! 
      if ( iprint > 2) then
         write(6,'(a)') '   overlap target - model '
         call  imprec(zs,n,n,5)
      endif
! 
      do i=1,m
         do j=1,n
            zx=(0.d+0,0.d+0)
            do k=1,n
               zx=zx+zcmod(i,k)*zs(k,j)
            enddo
            zwork(i,j)=zx
         enddo
      enddo
!      
      if ( iprint >= 2) then
         write(6,'(a)') '    projected target vectors '
         call  imprec(zwork,m,n,5)
      endif
!
!     orthogonaliser
!
      call loeworth(zwork,zorth,m,n,' target2model   ',1.5d+0,iprint)
! 
!      
      if ( iprint >= 2) then
         write(6,'(a)') '   projected & orthogonalized target vectors'
         call  imprec(zorth,m,n,5)
      endif
!         call scapro(zorth,zorth,zs,m,n)
!         write(6,'(a)') '   orthogonality check   '
!         call  primac(zs,n,n,n,1.d-5)
!
!     construire les matrices de transformation
!
      if( y_socoupled ) then
        call scapro(zcmod,zorth,zs,m,n)
      else
        call scapro(zorth,zcmod,zs,m,n)
      endif
!
      if ( iprint >= 2 ) then
        write(6,'(a)') '   transformation   '
        call  imprec(zs,n,n,5)
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
         if(y_socoupled) then
            write(6,'(2a,/)') '     effective hamiltonian ',                   &
                 'in the subspace spanned by so-coupled eigenvectors' 
         else 
            write(6,'(2a,/)') '     effective hamiltonian ',                   &
                 'in the subspace spanned by so-free eigenvectors'
         endif        
         if(iprint >= 1) call  imprec(zheffectif,n,n,5)
         write(6,'(1x,2a,f16.2)')  '        diagonal, au      diagonal, cm-1  '
         do i=1,n
         if(abs(ground-666.d+0) <= 1.d-5) then
         write(6,'(a2,i6,f20.10,f17.2)')'@E',i, real(zheffectif(i,i),kind=pd), &
           (real(zheffectif(i,i))-real(zheffectif(1,1)))*au2cm 
         else
         write(6,'(a2,i6,f20.10,f17.2)')'@E',i, real(zheffectif(i,i),kind=pd), &
                                (real(zheffectif(i,i),kind=pd)-ground)*au2cm
         endif
         enddo
      deallocate( zs, zwork, zorth )
!
      return
      end subroutine projinto
!
!
      subroutine scapro(za,zb,zc,m,n)
!      
!  zc(square)=za^\dag * zb (both rectangular, same size)      
!
      use nombres_mod
      implicit none 
      complex(kind=pd), intent(in) :: za(m,*),zb(m,*)
      complex(kind=pd), intent(out):: zc(n,*)
      integer,          intent(in) :: m,n
      complex(kind=pd)             :: zx
      integer                      :: i,j,k
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
      end subroutine scapro
!--------------------------------------------------------
      subroutine eigher(a,eva,max,n)
!      
! eigenvalues eva /eigenvectors a  of complex hermitian matrix a 
! wrapper for zheev
! a is destructed
! eigenvalues in ascending order
!     
      use nombres_mod
      implicit none
 
      complex(kind=pd),dimension(max,*),intent(inout) :: a
      integer, intent(in)                             :: max,n
      real(kind=pd),dimension(*),intent(out)          :: eva(*)
      complex(kind=pd)                                :: zstore 
      complex(kind=pd), allocatable                   :: work(:) 
      real(kind=pd), allocatable                      :: rwork(:)
      integer, allocatable                            :: iwork(:)
      integer                                         :: mfound, lwork, info
      mfound=n
      lwork=4*n
      allocate(work(lwork))
      allocate(rwork(7*n))
      allocate(iwork(5*n))
      call zheev( 'V','L', n, a, max,eva,work,lwork,rwork,info )
      if( info /= 0 ) then
      write(*,*)'   ----- zheev crashed. info:', info,'  --------' 
      stop
      endif
      deallocate(work,rwork,iwork)      
      return
      end subroutine eigher
     
      
      subroutine imprec(a,n,m,chiffres)
      
      use nombres_mod
      
      implicit none
      integer, intent(in)          ::  n        ! taille de la matrice
      integer, intent(in)          ::  m        ! nombre de colonnes a imprimer
      complex(kind=pd), intent(in) ::  a(n,n)   ! matrice complexe a imprimer
      integer, intent(in)          ::  chiffres ! bon chffres
      integer debut, fin, i,j, trouve, numline
       
      integer,parameter :: l120=120
      real(kind=pd) maxval
      character(len=l120) ligne
      character(len=30) fmat_c,fmat_i
      character(len=20),parameter :: dotzero= '.0000000000000000000'
      character(len=20),parameter :: blanc=   '.                   '
      integer                     :: n4,mant
      
!
!    imprime m premieres colonnes de la matrice carree
!    partie importante
!
        maxval=0.d0
        do i=1,n
          do j=1,m
            if(maxval < abs(a(i,j))) maxval=abs(a(i,j))
          enddo
        enddo
!
!       n4 chiffres avant point
!
        n4=max(int(log10(maxval))+2,2)      
!
!       mant chiffres apres
        mant=chiffres-int(log10(maxval))
        if( mant > 14 ) mant=14
        if( mant < 0 )  mant=0
!        
!       nombre de nimbres per ligne
!        
        numline=(l120-10)/(2*n4+2+2*mant)
        
        write(fmat_c,'(a,i2,a,i2,a,i2,a)')            &
         '(1x,i4,2x,',numline,'((2f',mant+n4+1,'.',mant,',1x)))'
        write(fmat_i,'(a,i2,a,i2,a,i2,a)')            &
         '(',mant+n4+6,'x,',numline,'((i4,',mant*2+n4*2-1,'x)))'
        
        debut=1
        fin=debut+numline-1

        do while (debut <= m )
          fin=min(fin,m)
          write(*,*)' '
            write(*,fmat_i)(i,i=debut,fin)
          do j=1,n
            ligne=' '
            write(ligne,fmat_c)j,( real(a(j,i),kind=pd),aimag(a(j,i)) ,i=debut,fin)
            do while (index(ligne,dotzero(1:mant+1)) > 0 )
              trouve=index(ligne,dotzero(1:mant+1))
              if(trouve > 0) then
                 ligne(trouve:trouve+mant)=blanc(1:mant+1)
                 if( trouve > 2 .and.       &
                   ( ligne(trouve-2:trouve-1) == ' 0' .or. ligne(trouve-2:trouve-1) == '-0' ) ) then
                   ligne(trouve-2:trouve-1)='  '
                 endif  
              endif   
            enddo
            write(*,'(a)') trim(ligne)
          enddo 
          fin=fin+numline
          debut=debut+numline
        enddo
      write(*,*) ' '      
      return
      end subroutine imprec

