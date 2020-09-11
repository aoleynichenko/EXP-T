************************************************************************
*                                                                      *
*                             puma                                     *
*                                                                      *
************************************************************************
      subroutine puma(mlar,msma,zlar,zsma,ii,jj,wtd)
c     
      implicit real*8(a-h,o-x)
      implicit complex*16(z)
      implicit logical(y)
      character*3 wtd
c
c     puts square matrix in a cell (ii,jj) of large supermatrix
c     or gets square matrix from a cell of large supermatrix
c     all is complex
c
      dimension zlar(mlar,*),zsma(msma,msma)
      ncells=mlar/msma
      if(mlar .ne. ncells*msma) then
        write(6,'(1x,a)')' subroutine puma:'
        write(6,'(1x,a,2i8)')' wrong supermatrix/matrix sizes',mlar,msma
        stop 
      endif
      if(wtd .ne. 'get' .and. wtd .ne. 'put') then
        write(6,'(1x,a)')' subroutine puma:'
        write(6,'(1x,a,a)')' do not know what to do:',wtd
        stop 
      endif
      istart=(ii-1)*msma
      jstart=(jj-1)*msma
      do i=1,msma
        do j=1,msma
          if(wtd .eq. 'put') then
          zlar(istart+i,jstart+j)=zsma(i,j)
          else
          zsma(i,j)=zlar(istart+i,jstart+j)
          endif
        enddo
      enddo
      return
      end

*      subroutine comaop(a,b,what,c,n)
*      implicit complex*16(a-h,o-z) 
*      dimension a(n,n),b(n,n),c(n,n)
*      complex*16, allocatable :: z(:,:),zz(:,:) 
*      character*1 what
*      nn=n*n
*      if(what.eq.'*')then
*      allocate( z(n,n) )
*      do i=1,n
*      do j=1,n
*      s=(0.d0,0.d0)
*      do k=1,n
*      s=s+b(i,k)*c(k,j)
*      enddo
*      z(i,j)=s
*      enddo
*      enddo
*      call zcopy(nn,z,1,a,1)
*      deallocate( z )
*      else if (what.eq.'+')then
*      allocate( z(n,n) )
*      do i=1,n
*      do j=1,n
*      z(i,j)=b(i,j)+c(i,j)
*      enddo
*      enddo
*      call zcopy(nn,z,1,a,1)
*      deallocate( z )
*      endif
*      return
*      end    
          
************************************************************************
*                                                                      *
*                             matpade                                  *
*                                                                      *
************************************************************************
      subroutine matpade(max,ll,mm,zc,zpade,impression)
c
c  [ll/mm] (right) rational pade approximant 
c          for a complex (max x max) matrix series zc
c        
      implicit real*8(a-h,o-x)
      implicit complex*16(z)
      implicit logical(y)
      dimension zc(max*max,*),zpade(max,max)
      complex*16, allocatable :: z(:,:),zb(:,:),za(:,:),zright(:,:) 
     & , zbt(:,:)
*czai(((((((     
*czai)))))))
*czai)))))))
*      complex*16, allocatable :: zcr(:,:),zb1(:,:),zb2(:,:),za1(:,:),
*     &            zw1(:,:),zw2(:,:)
*czai)))))))
*czai)))))))
*czai)))))))
      if(impression .gt. 0) then
          write(6,'(/4x,a,i1,a,i1,a/)') '----- [',ll,'/',mm,
     &     '] matrix rational pade approximant -----'
      endif 
c
c   input data
c
      if(impression .ge. 10) then
        do i=1,ll+mm+1
          write(6,'(1x,a,i2)') ' series: term',i-1
          call primac(zc(1,i),max,max,max,1.d-5)
        enddo
      endif
c
      if(mm.le.0) then
c
c   trivial [ll/0] pade approximant = taylor partial sum
c
        ii=0
        do j=1,max
        do i=1,max
          ii=ii+1   
          zpade(i,j)=(0.d0,0.d0)
          do k=1,ll+mm+1
            zpade(i,j)=zpade(i,j)+zc(ii,k)
          enddo
        enddo
        enddo
        return
c      
      endif      
c      
czai)))))))
*czai)))))))
*czai)))))))
*      if( ll.eq.0 .and. mm.eq.2 ) then
*         maxmax=max*max
*         allocate( zcr(max,max),zw1(max,max),zw2(max,max) )
*         allocate( zb1(max,max) )
*              write(*,*) ' un'
*         call zcopy(maxmax,zc(1,1),1,zcr,1)
*         call zinver(zcr,max)
*         call comaop(zw1,zcr,'*',zc(1,2),max)
*c  zw1=c0^{-1}*c1 
*         call comaop(zw2,zw1,'*',zw1,max)
*c  zw2=( c0^{-1}*c1 )^2      
*         do i=1,max
*         do j=1,max
*            zb1(i,j)=-zw1(i,j)+zw2(i,j)
*            if(i .eq. j)zb1(i,j)=zb1(i,j) +(1.d0,0.d0) 
*         enddo
*         enddo
*c         zb1= 1 - c0^{-1}*c1 + ( c0^{-1}*c1 )^2
*         call comaop(zw1,zcr,'*',zc(1,3),max)
*c        zw1=c0^{-1}*c2
*         do i=1,max
*         do j=1,max
*            zb1(i,j)=zb1(i,j)-zw1(i,j)
*         enddo
*         enddo
*c         zb1= 1 - c0^{-1}*c1 + ( c0^{-1}*c1 )^2 - c0^{-1}*c2
*         call zinver(zb1,max)
*         call comaop(zpade,zc(1,1),'*',zb1,max)                             
*         return     
*      endif
*czai)))))))
*czai)))))))
*czai)))))))
c
c   allo
c
      allocate ( z(max*mm,max*mm),zbt(max*mm,max),zright(max*mm,max) )
c
c   clean z
c
      do i=1,max*mm
        do j=1,max*mm
          z(i,j)=(0.d0,0d0)
        enddo
      enddo
c
c   fill 
c
      mlar=max*mm
      do i=1, mm
        call puma(mlar,max,zright,zc(1,ll+i+1),i,1,'put')
        do j=1,mm
          if( ll+i-j+1 .ge. 1 .and. ll+i-j+1 .le. ll+mm+1 ) then
            call puma(mlar,max,z,zc(1,ll+i-j+1),i,j,'put')
          endif
        enddo
      enddo
c
c   rhs sign
c  
      do i=1,mlar
        do j = 1,max
          zright(i,j)=-zright(i,j)
        enddo
      enddo
      if(impression.gt.15)then
        write(6,'(1x,a)')' leq system matrix'
        call primac(z,mlar,mlar,mlar,1.d-5)
        write(6,'(1x,a)')' right hand side'
        call primac(zright,mlar,mlar,max,1.d-5)
      endif
c
c   invert z          
c
      call zinver(z,mlar)
      if(impression.gt.15)then
        write(6,'(1x,a)')' inverted leq system matrix'
        call primac(z,mlar,mlar,mlar,1.d-5)
      endif
c
c   compute denominator matrix coefs
c
      do i=1,mlar
        do j=1,max
          zs=(0.d0,0.d0)
          do k=1,mlar
            zs=zs+z(i,k)*zright(k,j)
          enddo
          zbt(i,j)=zs
        enddo
      enddo 
c
c    plus besoin de z,zright
c
      deallocate (z,zright)
      allocate ( za(max*max,ll+1) ) 
      allocate ( zb(max*max,mm) ) 
c
c   extract denominator matrix coefs to zb      
c
      do i=1,mm
        call puma(mlar,max,zbt,zb(1,i),i,1,'get')
        if(impression .gt. 10) then
          write(6,'(1x,a,i1)')' b',i
          call primac(zb(1,i),max,max,max,1.d-5)
        endif        
      enddo
c      
      deallocate (zbt)
c
c   compute numerator and put it into za
c
      maxmax=max*max
      do i=0,ll
c  
        call zcopy(maxmax,zc(1,i+1),1,za(1,i+1),1)
        if(i .gt. 0) then
          min_mm_i=i
          if(mm .lt. min_mm_i) min_mm_i=mm
          if(min_mm_i .ge. 1) then
            do k=1,min_mm_i
            call summat( za(1,i+1), zc(1,i+1-k), zb(1,k),max )
            enddo
          endif 
        endif       
        if(impression .gt. 10) then
          write(6,'(1x,a,i1)')' a',i
          call primac(za(1,i+1),max,max,max,1.d-5)
        endif 
      enddo
c
c     compute total numerator for z=1 in za( ,1)
c   
      if(ll.gt.0) then 
        do i=2,ll+1 
          do j=1,maxmax 
            za(j,1)=za(j,1)+za(j,i)
          enddo
        enddo 
      endif   
c
c     compute total denominator for z=1 in zb( ,1)
c     
      if (mm .gt. 1) then
        do i=2,mm
          do j=1,maxmax
          zb(j,1)=zb(j,1)+zb(j,i)
          enddo
        enddo
      endif  
c add unit
      iii=0
      do i=1,max
        do j=1,max
          iii=iii+1
          if(i.eq.j) zb(iii,1)=zb(iii,1)+(1.d0,0.d0)
        enddo
      enddo
        if(impression .gt. 10) then
          write(6,'(1x,a)')' total numerator'
          call primac(za(1,1),max,max,max,1.d-5)
          write(6,'(1x,a)')' total denominator'
          call primac(zb(1,1),max,max,max,1.d-5)
        endif        
c  invert
      call zinver(zb,max)
c
c     compute [ll/mm] value at z=0 and put it to zpade
c
      do i=1,max
        do j=1,max
          zpade(i,j)=(0.d0,0.d0)
        enddo
      enddo
      call summat( zpade, za, zb, max ) 
      if( impression.ge.10)then
        write(6,'(1x,a)')' approximant at z=1'                    
        call primac(zpade,max,max,max,1.d-5)
      endif  
c    
c   deallo
c
      deallocate (zb,za) 
c           
      return
      end
************************************************************************
*                                                                      *
*                             summat                                   *
*                                                                      *
************************************************************************
      subroutine summat(za,zb,zc,max)
c     
      implicit real*8(a-h,o-x)
      implicit complex*16(z)
      implicit logical(y)
c
c     
c     a=a+b*c
c
      dimension za(max,max),zb(max,max),zc(max,max)
c
      do i=1,max
        do j=1,max
          do k=1,max
            za(i,j)=za(i,j)+zb(i,k)*zc(k,j)
          enddo
        enddo
      enddo
      return
      end
            
