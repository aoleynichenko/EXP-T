!***********************************************************************
!                                                                      *
!                             readheff                                 *
!                                                                      *
!***********************************************************************
      subroutine readheff(ycarith,iunit,maxmax,mainms,nskip,  &
                          sector,irrep,e,vg,vd,llmmpade,y_cut,y_dc,arra,iprint)
!
!     input 
!     ycarith        logical, is complex arithmetic used
!     iunit          logical unit for heff
!     maxmax         heff size
!     mainms         number of main states
!     nskip          number of skipped low-energy states
!     sector         0h1p or 1h0p or 0h2p etc
!     irrep          representation
!     llmmpade       encoded [ll/mm]
!     y_cut          true if cropped heff are used 
!     y_dc           true if heff is transformed to des cloiseaux form
!     arra           rule to compose heff or to transform the heff series
!     iprint         printing regime
!
!     output
!     e              1d array of eigenvalues (at least formally complex)
!     vg, vd         left and right eigenvectors, possibly maxmax*maxmax
!
      use nombres_mod
!  
      implicit none
      logical, intent(in):: ycarith,y_cut,y_dc
      integer, intent(in):: iunit, maxmax, mainms, nskip,irrep,llmmpade,iprint
      character(len=4)   :: sector, sector2
      character(len=128) :: ligne
      character(len=80)  :: arrange, arra, bla
      logical y_combine, ysector
      integer            :: llpade,mmpade,lseries,iw,i_file, mm, mdummy,i,j,k, l, inu
      complex(kind=pd),intent(out) :: vg(maxmax,maxmax),vd(maxmax,maxmax),e(*)
!
      complex(kind=pd), allocatable :: h(:,:)      
      complex(kind=pd), allocatable :: smn(:,:)
      complex(kind=pd), allocatable :: work(:)
      complex(kind=pd), allocatable :: hseries(:,:,:)
      complex(kind=pd), allocatable :: cuttedhseries(:,:,:),vdmain(:,:),vgmain(:,:)
      real(kind=pd), allocatable    :: hr(:,:)
!
      bla ='                                                                                '
      y_combine=.false.
      llpade=0
      mmpade=0
      lseries=1
      arrange=arra
!      
      iw=index(arrange(1:80),'H')
      if(iw.ge.1) then
        y_combine=.true.
        if(llmmpade.ne.0)then
          write(*,'(/,a)')'    helas! '
          write(*,'(a,a)')'    combined scheme ',arrange(1:80)
          write(*,'(a)')'    is not compatible with the [L/M] option'
          write(*,'(a,/)')'    i am forced to stop'
          stop
        endif
      endif  
!
!   allo
!
      allocate (h(maxmax,maxmax))
      allocate (hr(maxmax,maxmax))
      allocate (smn(mainms,mainms))
!
!   if a pade approximant is to be computed
!   or composite scheme used 
!   count the effective hamiltonials expected in a single file 
!      
      nombre_heff: if (llmmpade.ne.0 .or. y_combine) then
        if(llmmpade.gt.0) then
          llpade=llmmpade/10
          mmpade=llmmpade-llpade*10
          lseries=llpade+mmpade+2
        endif
        if(llmmpade.lt.0) then
          llpade=(-llmmpade)/10
          mmpade=(-llmmpade)-llpade*10
          lseries=llpade+mmpade+1
        endif 
        if(y_combine)then
          do i=3,80
            if(arrange(i:i) .eq. 'H') then
              read(arrange(i+1:i+1),'(i1)')inu
              if(inu .gt. lseries)lseries=inu
            endif
          enddo
        endif         
!
!   space for pade approximant calculation
!      
        allocate(hseries(maxmax,maxmax,lseries))
        if(y_cut)then
          allocate(cuttedhseries(mainms,mainms,lseries))
          allocate(vdmain(mainms,mainms))
          allocate(vgmain(mainms,mainms))
        endif        
      endif nombre_heff
!   
!    trouver un bon secteur / rep
!  
!
      i_file=0
!
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
           if(mdummy == irrep) goto 14
         endif  
         goto 13
  14  continue
      if( mm /= maxmax) then
        write(6,*)' wrong heff size', mm
        stop
      endif
!
!   required effective hamiltonian found
!
!   new format (luuk visscher, dirac 17+)
!
      if(ycarith) then
        read(iunit,'(4e21.12)') h
      else
        read(iunit,'(4e21.12)') hr
        h=cmplx(hr,0.d0,kind=pd)
      endif
!
!   end new format of eff hamiltonian reading
!
!   if pade approximation or combined scheme used, move heff to hseries and go on reading
!
!
      plusieurs_heff: if( (llmmpade /= 0) .or. y_combine) then
        i_file = i_file + 1
        hseries(:,:,i_file)=h(:,:)
        if(i_file < lseries) then
!
!    na dvorou byl kol, na kolou motchalo
!        
          ysector=.false.
          mm=0
          goto 13
        else
          write(6,'(4x,i2,a)')lseries,' required effective hamiltonians found'
!
!   the effective hamiltonian sequence is long enough and looks ok
!
        endif
!
!    rearrange if required
!
        if(arrange(1:40) /= bla(1:40))  call rear(hseries,h,arrange,maxmax,lseries)
!
      endif plusieurs_heff
      
      if(y_combine) then
        deallocate(hseries)
      endif
!
!    end of rearrangement
!
      padecalc: if(llmmpade /= 0) then
!
!     pade
!      
        h_coupe: if(y_cut) then
!
!   projection if required
!
!
!   compute eigenvectors of the best term (= last one)
!
          call eiggen(h,vg,vd,e,maxmax) 
          skipping: if ( nskip >= 1 ) then
!
!       vectors / energies to be skipped go to the end of the list 
!
            allocate (work(maxmax))
!
            do i=1,maxmax
              if ( i <= nskip ) then
                work(maxmax-nskip+i)=e(i)
              else
                work(i-nskip)=e(i)
              endif  
            enddo
            e(1:maxmax)=work(1:maxmax)
!        
            do k=1,maxmax
!        
              do i=1,maxmax
                if ( i <= nskip ) then
                   work(maxmax-nskip+i)=vg(k,i)
                else
                   work(i-nskip)=vg(k,i)
                endif  
              enddo
              vg(k,:)=work
!        
              do i=1,maxmax
                if ( i <= nskip ) then
                  work(maxmax-nskip+i)=vd(k,i)
                else
                  work(i-nskip)=vd(k,i)
                endif  
              enddo
              vd(k,:)=work
!
            enddo
            deallocate (work)
!
!    vecteurs sont rearranges   
!
          endif skipping
!            
!    main eigenvectors span the space to project onto
!
          call loeworth(vd,vg,maxmax,mainms,' proj for pade  ',1.d-1,iprint)
!
!   orthogonalized main RIGHT eigenvectors are in <<vg>>
!
          do k=1,lseries
!
            do i=1,maxmax
              do j=1,mainms
                vd(i,j)=(0.d+0,0.d+0)
                do l=1,maxmax
                  vd(i,j) = vd(i,j) + hseries(i,l,k)*vg(l,j)
                enddo
              enddo
            enddo
!               
            do i=1,mainms
              do j=1,mainms
                cuttedhseries(i,j,k)=(0.d+0,0.d+0)
                do l=1,maxmax
                  cuttedhseries(i,j,k)=cuttedhseries(i,j,k)+conjg(vg(l,i))*vd(l,j)
                enddo
              enddo
            enddo
          enddo
!
!   note: vd are destroyed and orthogonalised right vectors in vg, kept
! 
!
!   cropped series is converted into difference series 
!          
         do k=lseries,2,-1
           cuttedhseries(:,:,k) = cuttedhseries(:,:,k)-cuttedhseries(:,:,k-1) 
         enddo        
!
         smn=(0.d0,0.d0)

         if(llmmpade > 0)then
!
!    build the approximant for the sum 2 through lseries in h 
!
           call matpade(mainms,llpade,mmpade,cuttedhseries(:,:,2),smn,iprint)
!
!    now add the first term in series
!
            smn=smn+cuttedhseries(:,:,1)
!
          else
!
!    build the approximant for the whole trimmed series 
!
            call matpade(mainms,llpade,mmpade,cuttedhseries,smn,iprint)
          endif 
!
!    smn contains the trimmed effective hamiltonian  
!                           
          if(iprint > 2) then
            write(6,*) ' trimmed & extrapolated effective hamiltonian ',iunit-10 
            call imprec(smn,mainms,mainms,5)
          endif
!
!    trimmed effective Hamiltonian diagonalization
!            
          call eiggen(smn,vgmain,vdmain,e,mainms)
!
!    orthogonalize the minivectors
!            
          call loeworth(vdmain,vgmain,mainms,mainms,'  p(heff)p evecs',1.1d+0,iprint)

!
!    restore the whole vectors. dc form of small hamiltonian is supposed
!
          do i=1,maxmax
            do j=1,mainms
              vd(i,j)=(0.d0,0.d0)
              do k=1,mainms
                vd(i,j)=vd(i,j)+vg(i,k)*vgmain(k,j)
              enddo
            enddo
          enddo
          vg=vd  
          impression: if(iprint >= 3)then
            write(6,*)'  eigenvalues: cropped hamiltonian',iunit-10 
            do i=1,mainms
              write(6,'(i6,f15.8,f13.8)')i,real(e(i)),aimag(e(i))
            enddo
            write(6,*)' ' 
            write(6,*)'  eigenvectors: cropped hamiltonian',iunit-10 
            call imprec(vd,maxmax,maxmax,5)
            endif impression
          deallocate(cuttedhseries,vdmain,vgmain)        
          deallocate (h,hr,smn,hseries)
          return
!
!        ici finit la version heff projete
!          
        endif h_coupe 
!
!
!   full-sized pade approximant calculations
!      
!   convert sequence into series
!
        do k=lseries,2,-1
          hseries(:,:,k)=hseries(:,:,k)-hseries(:,:,k-1)
        enddo        
!
!    clean h 
!
        h(1:maxmax,1:maxmax)= (0.d0,0.d0) 
        if(llmmpade.gt.0)then
!
!    build the approximant for the sum 2 through lseries in h 
!
          call matpade(maxmax,llpade,mmpade,hseries(1,1,2), h,iprint)
!
!    now add the first term in series
!
          h(:,:)=h(:,:) + hseries(:,:,1)              
        else
!
!    build the approximant for the whole lseries in h directly
!
          call matpade(maxmax,llpade,mmpade,hseries,h,iprint)
        endif                  
!
!   no more need of hseries space      
!
        deallocate(hseries)
!          
!   end of pade approximant calculation
!
      endif padecalc
!      
!
      if(iprint > .2) then
        write(6,*) ' effective hamiltonian ',iunit-10 
        call imprec(h,maxmax,maxmax,5)
      endif
!      
!   diagonaliser l'hamiltonien total
!
      call eiggen(h,vg,vd,e,maxmax) 
!
!
!
            if(y_dc)then 
!            
              call loeworth(vd,vg,maxmax,maxmax,'  heff eigenvecs',1.1d+0,iprint)
              vd=vg
            else
!            
              call loeworth(vd,vd,maxmax,maxmax,'? heff eigenvecs',1.1d+0,iprint)
!
            endif
!            
      if ( nskip >= 1 ) then
!       vectors / energies to be skipped go to the end of the list 
        allocate (work(maxmax))
!
        do i=1,maxmax
          if ( i <= nskip ) then
            work(maxmax-nskip+i)=e(i)
          else
            work(i-nskip)=e(i)
          endif  
        enddo
        e(1:maxmax)=work(1:maxmax)
!        
        do k=1,maxmax
!        
           do i=1,maxmax
             if ( i <= nskip ) then
               work(maxmax-nskip+i)=vg(k,i)
             else
               work(i-nskip)=vg(k,i)
             endif  
           enddo
           vg(k,1:maxmax)=work(1:maxmax)
!        
           do i=1,maxmax
             if ( i <= nskip ) then
               work(maxmax-nskip+i)=vd(k,i)
             else
               work(i-nskip)=vd(k,i)
             endif  
           enddo
           vd(k,1:maxmax)=work(1:maxmax)
!
        enddo
        
        deallocate (work)
!       vecteurs sont rearranges   
      endif
!
!
!
      if(iprint >= 3)then
      write(6,*)'  eigenvalues of the hamiltonian',iunit-10 
      do i=1,maxmax
      write(6,'(i6,f15.8,f13.8)') i, real(e(i),kind=pd),aimag(e(i))
      enddo
      write(6,*)' ' 
      write(6,*)'  right eigenvectors of the hamiltonian',iunit-10 
      call imprec(vd,maxmax,mainms,5)
      write(6,*)'  left eigenvectors of the hamiltonian',iunit-10 
      call imprec(vg,maxmax,mainms,5)
      endif
!
!
      deallocate (h,hr,smn)
!      
      return
      end
!--------------------------------------------------------------------
!               rear
!--------------------------------------------------------------------
      subroutine rear(hseries,h,arrange,maxmax,lseries)
      implicit none
      complex*16 :: h(maxmax,maxmax),hseries(maxmax,maxmax,lseries)
      real*8     :: sign_coef,coef
      integer    :: maxmax,lseries
      integer    :: iarr, i,j,k, iold,inew, bad, lentr, iw
      logical    :: yyy
      integer, parameter :: l80=80
      character(len=l80) :: arrange
      character(len=l80),parameter :: bla=' '
!
      lentr=bad(arrange,l80)
! 
!   voyons si on fait addition ou series transformation
!    
      i=index(arrange,'H')
      if(i.ge.1) then
         yyy=.true.
      else
         yyy=.false.
      endif
      
!
      if (yyy) then
!
!     composition
!
        write(*,'(/,4x,2a,/)')'--- composing heff           : ',arrange
!       initialisation
        h=(0.d0,0.d0)
        iw=index(arrange,'=')
        if(iw > 0) then
          arrange(1:iw)=bla(1:iw)
        else
          write(*,'(/,4x,2a)')'--- suspicious composite scheme ',arrange
          stop 
        endif
        do iarr = 1,lentr
          if(arrange(iarr:iarr).ne.'H') cycle
          read(arrange(iarr+1:iarr+1),'(i1)') k
          coef=1.d0
          select case (arrange(iarr-1:iarr-1))
            case('+',' ')
              coef=1.d0
            case('-')
              coef=-1.d0
            case default
              read(arrange(1:iarr-1),*) coef  
          end select
          arrange(1:iarr+1)=' '
          h(1:maxmax,1:maxmax)=h(1:maxmax,1:maxmax)+hseries(1:maxmax,1:maxmax,k)*coef
            write(*,'(4x,a,i1,a,f9.6,a)')'--- heff: H(',k,') times ',coef,' added'
        enddo 
!      
!
      else
!
!     transformation de la serie
!
        write(*,'(/,4x,2a,/)')'--- heff series rearrangement: ',arrange
!
!    sequence -> serie
!
          do k=lseries,2,-1
            hseries(1:maxmax,1:maxmax,k)=hseries(1:maxmax,1:maxmax,k)-hseries(1:maxmax,1:maxmax,k-1)
          enddo
!  
!   converted      
!      
        iarr=0
        lire: do while(iarr .lt. lentr)
          iarr=iarr+1
          neq: if(arrange(iarr:iarr).eq.'h') then
!
!  nouvelle formule
!
            iarr=iarr+1
            if (iarr .gt. lentr) exit
            read(arrange(iarr:iarr),'(i1)') inew
            write(*,'(4x,a,i1,a)')'--- reconstruction of h(',inew,'): ------'
            arrange(iarr-1:iarr)='  '
!
!   cherche "="
!
            do while(arrange(iarr:iarr).ne.'=')
              arrange(iarr:iarr)=' '
              iarr=iarr+1
            enddo
            arrange(iarr:iarr)=' '
!
!  decode composition
!
            sign_coef=1.d+0
            decode: do while(iarr .lt. lentr)
              iarr=iarr+1
              select case (arrange(iarr:iarr))
                case (' ')
                  cycle
                case ('+')  
                  sign_coef=1.d+0
                  arrange(iarr:iarr)=' '
                  cycle
                case ('-')  
                  sign_coef=-1.d+0
                  arrange(iarr:iarr)=' '
                  cycle
                case (',',';')
!   equation nouvelle
                  arrange(iarr:iarr)=' '                  
                  exit
                case ('h') 
!
!  ajouter l'hamiltonien effectif suivant
!
                  read(arrange(iarr+1:iarr+1),'(i1)')iold
                  if(arrange(1:iarr-1).eq.bla(1:iarr-1))then
                    coef=sign_coef
                  else
                    read(arrange(1:iarr-1),*)coef
                    arrange(1:iarr-1)=bla(1:iarr-1)
                    coef=coef*sign_coef
                  endif
!
                  do i=1,maxmax
                    do j=1,maxmax
                      if(iold .eq. inew) then
                        hseries(i,j,inew)=hseries(i,j,inew)*coef
                      else
                        hseries(i,j,inew)=hseries(i,j,inew) + hseries(i,j,iold) * coef
                      endif
                    enddo
                  enddo  
!
                  write(*,'(4x,a,i1,a,f8.4)')'--- h(',iold,') added with weight',coef
                  arrange(iarr:iarr+1)='  '
                  iarr=iarr+1
              end select
!               
            enddo decode
!               
          endif neq
!
        enddo lire
!      
!    serie -> sequence 
!
        do k=2, lseries,1
          hseries(1:maxmax,1:maxmax,k)=hseries(1:maxmax,1:maxmax,k)+hseries(1:maxmax,1:maxmax,k-1)
        enddo
!
      endif   
!      
      write(*,'(a)')' '
      return
      end subroutine rear
!--------------------------------------------------------------------
!               bad
!--------------------------------------------------------------------
!          blancs a droite, resultat=nombre de non-blancs
!
      integer function bad(a,l)
      implicit none
      integer          :: i,j,l,lt
      character(len=l) :: a
      j = 0
      do i = 1,len(trim(a))
        if(a(i:i) .eq. ' ') cycle
        j=j+1
        a(j:j)=a(i:i)
      enddo
      if (j .lt. l) then
        a(j+1:l)=' '
      endif
      bad=j
      return
      end function bad   
!----------------------------------------------------------------------     
