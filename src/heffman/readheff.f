************************************************************************
*                                                                      *
*                             readheff                                 *
*                                                                      *
************************************************************************
      subroutine readheff(ycarith,iunit,maxmax,mainms,nskip,
     &  sector,irrep,e,vg,vd,llmmpade,y_cut,y_dc,arrange,iprint)
c
c     input 
c     ycarith        logical, is complex arithmetic used
c     iunit          logical unit for heff
c     maxmax         heff size
c     mainms         number of main states
c     nskip          number of skipped low-energy states
c     sector         0h1p or 1h0p or 0h2p etc
c     irrep          representation
c     llmmpade       encoded [ll/mm]
c     y_cut          true if cropped heff are used 
c     y_dc           true if heff is transformed to des cloiseaux form
c     iprint         printing regime
c
c     output
c     e              1d array of eigenvalues (at least formally complex)
c     vg, vd         left and right eigenvectors, possibly maxmax*maxmax
c
      implicit real*8(a-h,o-x,z)
      implicit logical(y)
      character*4 sector
      character*128 ligne
      character*80 arrange
      complex*16 vg(maxmax,maxmax),vd(maxmax,maxmax),e(*)
c
      complex*16, allocatable :: h(:,:)      
      complex*16, allocatable :: smn(:,:)
      complex*16, allocatable :: work(:)
      complex*16, allocatable :: hseries(:,:,:)
      complex*16, allocatable :: cuttedhseries(:,:,:),vdmain(:,:)
     &                           ,vgmain(:,:)
      real*8, allocatable :: hr(:,:)
c      real*8 arr(9)
      character*4 sector2
      character*80 bla
      bla(1:40)='                                        '
      bla(41:80)='                                        '
c
c   allo
c
      allocate (h(maxmax,maxmax))
      allocate (hr(maxmax,maxmax))
      allocate (smn(mainms,mainms))
c
c   count the effective hamiltonials glued in a single file 
c
      i_file=0
c      
      if(llmmpade.ne.0)then
        if(llmmpade.gt.0) then
          llpade=llmmpade/10
          mmpade=llmmpade-llpade*10
          lseries=llpade+mmpade+2
        else
          llpade=(-llmmpade)/10
          mmpade=(-llmmpade)-llpade*10
          lseries=llpade+mmpade+1
        endif          
c
c   space for pade approximant calculation
c      
        allocate(hseries(maxmax,maxmax,lseries))
        if(y_cut)then
          allocate(cuttedhseries(mainms,mainms,lseries))
          allocate(vdmain(mainms,mainms))
          allocate(vgmain(mainms,mainms))
        endif        
      endif
c     
      ysector=.false.
      mm=0
  13  read(iunit,'(a)',err=14,end=14) ligne(1:80)
         iw=index(ligne(1:80),'sector')
         if(iw.ge.1) then
           read (ligne(1:80),'(a)') sector2
           if (sector.eq.sector2.and. (.not.ysector)) ysector=.true.
           if (sector.ne.sector2.and. ysector) ysector=.false.
           goto 13
         endif  
         if(.not. ysector)goto 13
         iw=index(ligne(1:80),'heff size') 
         if(iw.ge.1 .and. ysector) then
           read (ligne(1:80),*) mdummy,mm
           if(mdummy.eq.irrep) goto 14
         endif  
         goto 13
  14  continue
      if( mm. ne. maxmax) then
        write(6,*)' wrong heff size', mm
        stop
      endif
c
c   required effective hamiltonian found
c
c   new format (luuk visscher, dirac 17+)
c
      if(ycarith) then
        read(iunit,'(4e21.12)') h
      else
        read(iunit,'(4e21.12)') hr
        h=cmplx(hr,0.d0)
      endif
c
c   end new format of eff hamiltonian reading
c
c   if pade approximation used, move heff to hseries and go on reading
c
      if(llmmpade .ne. 0) then
        i_file = i_file + 1
        do iii=1,maxmax
          do jjj=1,maxmax
            hseries(iii,jjj,i_file)=h(iii,jjj)
          enddo
        enddo
        if(i_file .lt. lseries) then
c
c    na dvorou byl kol, na kolou motchalo
c        
          ysector=.false.
          mm=0
          goto 13
        else
          write(6,'(4x,i2,a)')lseries,
     &                        ' required effective hamiltonians found'
c
c   the effective hamiltonian sequence is long enough and read ok
c
c    rearrange if required
c
         if(arrange(1:40).ne.bla(1:40))  then
c
c
           write(6,'(/,4x,2a)')'--- heff series rearrangement: ',arrange
c      
c   convert sequence into series
c
          do kkk=lseries,2,-1
            do iii=1,maxmax
              do jjj=1,maxmax
                hseries(iii,jjj,kkk)=hseries(iii,jjj,kkk)
     &                              -hseries(iii,jjj,kkk-1)
              enddo
            enddo  
          enddo
c  
c   converted      
c
           iarr=0
           do while(iarr.lt.80)
             iarr=iarr+1
             if(arrange(iarr:iarr).eq.'h') then
c  new formula starts
               iarr=iarr+1
               if (iarr.gt.80) exit
               read(arrange(iarr:iarr),'(i1)') inew
               write(*,'(4x,a,i1,a)')'--- reconstruction of h(',inew,
     &                                 '): ------'
               arrange(iarr-1:iarr)='  '
c   find "="
               do while(arrange(iarr:iarr).ne.'=')
                 arrange(iarr:iarr)=' '
                 iarr=iarr+1
               enddo
               arrange(iarr:iarr)=' '
c  decode composition
               sign_c=1.d+0
               do while(iarr .lt. 80)
                 iarr=iarr+1
c                 if(iarr.ge.80)exit
                 if(arrange(iarr:iarr).eq.' ') then
                   cycle
                 elseif(arrange(iarr:iarr).eq.'+') then
                   sign_c=1.d+0
                   arrange(iarr:iarr)=' '
                   cycle
                 elseif(arrange(iarr:iarr).eq.'-') then
                   sign_c=-1.d+0
                   arrange(iarr:iarr)=' '
                   cycle
                 elseif(arrange(iarr:iarr).eq.',') then
c  go to a new equation
                   arrange(iarr:iarr)=' '                  
                   exit 
                 elseif(arrange(iarr:iarr).eq.'h') then
c  next effective hamiltonian  will be added                  
c                   write(*,*)'arrange', arrange 
                   read(arrange(iarr+1:iarr+1),'(i1)')iold
                   if(arrange(1:iarr-1).eq.bla(1:iarr-1))then
                     arr=sign_c
                   else
                     read(arrange(1:iarr-1),*)arr
                     arrange(1:iarr-1)=bla(1:iarr-1)
                     arr=arr*sign_c
                   endif
c
                   do iii=1,maxmax
                     do jjj=1,maxmax
                       if(iold .eq. inew) then
                         hseries(iii,jjj,inew)=hseries(iii,jjj,inew)
     &                              * arr
                       else
                         hseries(iii,jjj,inew)=hseries(iii,jjj,inew)
     &                              + hseries(iii,jjj,iold) * arr
                       endif
                     enddo
                   enddo  
c
                   write(6,'(4x,a,i1,a,f8.4)')'--- h(',iold,
     &               ') added with weight',arr
                   arrange(iarr:iarr+1)='  '
                   iarr=iarr+1
                 endif
               enddo
             endif
           enddo
c      
c   convert back - series into sequence 
c
          do kkk=2, lseries,1
            do iii=1,maxmax
              do jjj=1,maxmax
                hseries(iii,jjj,kkk)=hseries(iii,jjj,kkk)
     &                              +hseries(iii,jjj,kkk-1)
              enddo
            enddo  
          enddo
c  
c   converted      
c
         endif
c
c    end of rearrangement
c
         if(y_cut) then
c
c   projection if required
c
c
c   compute eigenvectors of the best term (= last one)
c
           call eiggen(h,vg,vd,e,maxmax) 
           if ( nskip.ge.1 ) then
c
c       vectors / energies to be skipped go to the end of the list 
c
             allocate (work(maxmax))
*
             do iii=1,maxmax
               if ( iii .le. nskip ) then
                 work(maxmax-nskip+iii)=e(iii)
               else
                 work(iii-nskip)=e(iii)
               endif  
             enddo
             do iii=1,maxmax
               e(iii)=work(iii)
             enddo
*        
             do kkk=1,maxmax
*        
               do iii=1,maxmax
                 if ( iii .le. nskip ) then
                   work(maxmax-nskip+iii)=vg(kkk,iii)
                 else
                   work(iii-nskip)=vg(kkk,iii)
                 endif  
               enddo
               do iii=1,maxmax
                 vg(kkk,iii)=work(iii)
               enddo
*        
               do iii=1,maxmax
                 if ( iii .le. nskip ) then
                   work(maxmax-nskip+iii)=vd(kkk,iii)
                 else
                   work(iii-nskip)=vd(kkk,iii)
                 endif  
               enddo
               do iii=1,maxmax
                 vd(kkk,iii)=work(iii)
               enddo
*
             enddo
             deallocate (work)
c
c    vecteurs sont rearranges   
c
           endif
c            
c    main eigenvectors span the space to project onto
c
           call loeworth(vd,vg,maxmax,mainms,
     &                   ' proj for pade  ',1.d-1,iprint)
c
c   orthogonalized main RIGHT eigenvectors are in <<vg>>
c
           do kkk=1,lseries
*
             do iii=1,maxmax
               do jjj=1,mainms
                 vd(iii,jjj)=(0.d+0,0.d+0)
                 do lll=1,maxmax
                   vd(iii,jjj)=vd(iii,jjj)
     &                        +hseries(iii,lll,kkk)*vg(lll,jjj)
                 enddo
               enddo
             enddo
c               
             do iii=1,mainms
               do jjj=1,mainms
                 cuttedhseries(iii,jjj,kkk)=(0.d+0,0.d+0)
                 do lll=1,maxmax
                   cuttedhseries(iii,jjj,kkk)=cuttedhseries(iii,jjj,kkk)
     &                        +conjg(vg(lll,iii))*vd(lll,jjj)
                 enddo
               enddo
             enddo
           enddo
c
c   note: vd are destroyed and orthogonalised right vectors in vg, kept
c 
c
c   cropped series is converted into difference series 
c          
           do kkk=lseries,2,-1
             do iii=1,mainms
               do jjj=1,mainms
                 cuttedhseries(iii,jjj,kkk)=
     &           cuttedhseries(iii,jjj,kkk)
     &          -cuttedhseries(iii,jjj,kkk-1)
               enddo
             enddo  
           enddo        
c
          do iii=1,mainms
            do jjj=1,mainms
              smn(iii,jjj)=(0.d0,0.d0)
            enddo
          enddo  
          if(llmmpade.gt.0)then
c
c    build the approximant for the sum 2 through lseries in h 
c
            call matpade(mainms,llpade,mmpade,cuttedhseries(1,1,2),
     &                   smn,iprint)
c
c    now add the first term in series
c
              do iii=1,mainms
                do jjj=1,mainms
                  smn(iii,jjj)=smn(iii,jjj)+cuttedhseries(iii,jjj,1)
                enddo
              enddo 
c
c
c
            else
c
c    build the approximant for the whole trimmed series 
c
              call matpade(mainms,llpade,mmpade,cuttedhseries,
     &                   smn,iprint)
            endif 
c
c    smn contains the trimmed effective hamiltonian  
c                           
            if(iprint.gt.2) then
              write(6,*) ' trimmed & extrapolated ',
     &                   'effective hamiltonian ',iunit-10 
              call primac(smn,mainms,mainms,mainms,5.d-5)
            endif
c
c    trimmed effective Hamiltonian diagonalization
c            
            call eiggen(smn,vgmain,vdmain,e,mainms)
c
c    orthogonalize the minivectors
c            
            call loeworth(vdmain,vgmain,mainms,mainms,
     &           '  p(heff)p evecs',1.1d+0,iprint)

c
c    restore the whole vectors. dc form of small hamiltonian is supposed
c
            do i3=1,maxmax
              do j3=1,mainms
                vd(i3,j3)=(0.d0,0.d0)
                do k3=1,mainms
                  vd(i3,j3)=vd(i3,j3)+vg(i3,k3)*vgmain(k3,j3)
                enddo
              enddo
            enddo
            vg=vd  
      if(iprint.ge.3)then
      write(6,*)'  eigenvalues of the cropped hamiltonian',iunit-10 
      do iii=1,mainms
      write(6,'(i6,f15.8,f13.8)') iii, dreal(e(iii)),dimag(e(iii))
      enddo
      write(6,*)' ' 
      write(6,*)'  eigenvectors of the cropped hamiltonian',iunit-10 
      call primac(vd,maxmax,maxmax,mainms,5.d-5)
      endif
            
            goto 1626

            
            write(6,'(/12x,a,a/)')
     &             '  eigenenergies, au  relative energies, cm-1 ', 
     &             '                   imag part, au'             
            do iii=1,mainms
              if( abs(aimag(e(iii))) .lt. 1.d-7) then
              write(6,'(10x,i3,f14.7,2f20.2)')iii,
     &        real(e(iii)),
     &        real(e(iii)-e(1))*219474.631280634d+0               
     &        , real(e(iii))*219474.631280634d+0   
              else
              write(6,'(10x,i3,f14.7,2f20.2,f19.7)')iii,
     &        real(e(iii)),
     &        real(e(iii)-e(1))*219474.631280634d+0               
     &        , real(e(iii))*219474.631280634d+0,aimag(e(iii))   
              endif            
            enddo   
            write(6,'(a)')' '
            if(iprint.gt.0) then
              write(6,*) '   eigenvectors  ',iunit-10 
              call primac(vdmain,mainms,mainms,mainms,5.d-5)
            endif

1626        continue
            deallocate(cuttedhseries,vdmain,vgmain)        
            deallocate (h,hr,smn,hseries)
            return
c
c
c  pour l'instant, ici finit la version heff projete
c
c          
         endif
c
c
c   pade approximant calculations
c      
c   convert sequence into series
c
          do kkk=lseries,2,-1
            do iii=1,maxmax
              do jjj=1,maxmax
                hseries(iii,jjj,kkk)=hseries(iii,jjj,kkk)
     &                              -hseries(iii,jjj,kkk-1)
              enddo
            enddo  
          enddo        
c
c    clean h 
c
          do iii=1,maxmax
            do jjj=1,maxmax
              h(iii,jjj)=(0.d0,0.d0)
            enddo
          enddo  
          if(llmmpade.gt.0)then
c
c    build the approximant for the sum 2 through lseries in h 
c
            call matpade(maxmax,llpade,mmpade,hseries(1,1,2),
     &                   h,iprint)
c
c    now add the first term in series
c
            do iii=1,maxmax
              do jjj=1,maxmax
                h(iii,jjj)=h(iii,jjj)+hseries(iii,jjj,1)
              enddo
            enddo 
          else
c
c    build the approximant for the whole lseries in h directly
c
            call matpade(maxmax,llpade,mmpade,hseries,
     &                   h,iprint)
          endif                  
c
c   no more need of hseries space      
c
          deallocate(hseries)
c          
c   end of pade approximant calculation
c
        endif
      endif
c
      if(iprint.gt.2) then
        write(6,*) ' effective hamiltonian ',iunit-10 
        call primac(h,maxmax,maxmax,maxmax,5.d-5)
      endif
      call eiggen(h,vg,vd,e,maxmax) 
c
c
c
            if(y_dc)then 
c            
               call loeworth(vd,vg,maxmax,maxmax,
     &           '  heff eigenvecs',1.1d+0,iprint)
               do iii=1,maxmax
                 do jjj=1,maxmax
                   vd(iii,jjj)=vg(iii,jjj)
                 enddo
               enddo
            else
c            
               call loeworth(vd,vd,maxmax,maxmax,
     &           '? heff eigenvecs',1.1d+0,iprint)
c
            endif
            
      if ( nskip.ge.1 ) then
c       vectors / energies to be skipped go to the end of the list 
        allocate (work(maxmax))
*
        do iii=1,maxmax
          if ( iii .le. nskip ) then
            work(maxmax-nskip+iii)=e(iii)
          else
            work(iii-nskip)=e(iii)
          endif  
        enddo
        do iii=1,maxmax
            e(iii)=work(iii)
        enddo
*        
        do kkk=1,maxmax
*        
           do iii=1,maxmax
             if ( iii .le. nskip ) then
               work(maxmax-nskip+iii)=vg(kkk,iii)
             else
               work(iii-nskip)=vg(kkk,iii)
             endif  
           enddo
           do iii=1,maxmax
               vg(kkk,iii)=work(iii)
           enddo
*        
           do iii=1,maxmax
             if ( iii .le. nskip ) then
               work(maxmax-nskip+iii)=vd(kkk,iii)
             else
               work(iii-nskip)=vd(kkk,iii)
             endif  
           enddo
           do iii=1,maxmax
               vd(kkk,iii)=work(iii)
           enddo
*
        enddo
        
        deallocate (work)
c       vecteurs sont rearranges   
      endif
c
c
c
      if(iprint.ge.3)then
      write(6,*)'  eigenvalues of the hamiltonian',iunit-10 
      do iii=1,maxmax
      write(6,'(i6,f15.8,f13.8)') iii, dreal(e(iii)),dimag(e(iii))
      enddo
      write(6,*)' ' 
      write(6,*)'  right eigenvectors of the hamiltonian',iunit-10 
      call primac(vd,maxmax,maxmax,mainms,5.d-5)
      write(6,*)'  left eigenvectors of the hamiltonian',iunit-10 
      call primac(vg,maxmax,maxmax,mainms,5.d-5)
      endif
c
c      deallocate (h,hr,hi,smn)
c      hi is outdated, deleted
c
      deallocate (h,hr,smn)
c      
      return
      end
