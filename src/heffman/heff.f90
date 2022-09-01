!     
      program heff
! 
      use nombres_mod
!  
      implicit none
      character(len=128)fone, ftwo, ligne      
      character(len=80) arrange
      character(len=4) sector,sector2
      complex(kind=pd), allocatable :: eone(:),etwo(:)
      complex(kind=pd), allocatable :: voneg(:,:),vtwog(:,:), &
                                             voned(:,:),vtwod(:,:), tdm1(:,:),tdm2(:,:),heffectif(:,:)
      complex(kind=pd) tmoment
      logical carith,ysector,y_ff,y_socoupled,y_cut,y_dc,yyy
      character(len=50), parameter :: fmt9999 =  '(1x,i3,a4,i3,f14.3,f12.6,f11.6,f12.6,f12.8,2f12.6)'
      character(len=50), parameter :: fmt9999s = '(1x,i3,a4,i3,f14.3,23x,        f12.6,f12.8,2f12.6)'
      character(len=10), parameter :: fmt9998 = '(/,a,a,/)' 
      integer,  parameter :: maxirreps=32 
      integer :: nombre_irreps,irrepset(maxirreps),nskipset(maxirreps), mainmsset(maxirreps)
      integer :: numero_irrep,maxmax, mainms, nskip, irrep       

      integer :: iprint, llmmpade, llpade, mmpade,mdummy,mmmm
      integer :: not_end, iw, iww, iwww, iwwww, iwwwww, i, iconf, iii,jjj, icompo, jcompo
      real(kind=pd) ethr,step,scale, ground,energie0h0pone,energie0h0ptwo, wei,thelevel,tenergy,stmoment

!
!     default (universal) 
!
      iprint=1
      ethr=30
      fone='                     '
      ftwo='                     '
      sector='0h0p'
      carith=.false.
      y_ff=.false.
      y_dc=.true.
      y_cut=.false.
      step=1.d-4
      scale=0.00d+0
      ground=666.d+00
      do i=1,maxirreps
        irrepset(i)=-99 
        mainmsset(i)=99999
        nskipset(i)=0
      enddo
      nombre_irreps=0


      llmmpade=0
      arrange='                                                                                '
      energie0h0pone=zero  
      energie0h0ptwo=zero  
!
!     chapeau
!
      write(*,'(/3a)')'    ---- effective hamiltonian manip program ----'
      write(*,'(3a/)')'    ----         le 06 janvier  2021, az     ----'
!
!     parse input file
!
         not_end=0
lecture: do while (not_end == 0)
           read(*,'(a)',iostat=not_end) ligne
           if(not_end .lt. 0) then
             exit
           elseif (not_end .gt. 0) then
             write(*,'(/3a)')'    ---- i cannot read the input stream      ----'
             stop
           endif
          ligne=adjustl(ligne)     
!
!        commentaire
!
           if( ligne(1:1) .eq. '#' ) cycle
           iw=index(ligne,'heff:')
           if(iw.ge.1) then
             y_ff=.false.
             cycle
           endif  
           iw=index(ligne,']')
           iww=index(ligne,'[')
           iwww=index(ligne,'/')
           if(iw.ge.1 .and. iww.ge.1 .and. iwww.ge.1) then
!
!   pade
!
           ligne(iw:iw)=' '
           ligne(iww:iww)=' '
           ligne(iwww:iwww)=' '
           iwwww=index(ligne,'f')
           if(iwwww >= 1) then
             ligne(iwwww:iwwww)=' '
           endif
           iwwwww=index(ligne,'cut')
           if(iwwwww.ge.1)then
             ligne(iwwwww:iwwwww+2)='   '
             y_cut=.true.
           endif
           read (ligne(1:128),*) llpade,mmpade
           llmmpade=10*llpade+mmpade
           if(iwwww >= 1)llmmpade=-llmmpade
             cycle
           endif 
!
!           
           select case (ligne(1:4)) 
             case ('file')
               ligne(1:5)='     '
               ligne=adjustl(ligne)     
               if(fone(1:20) .eq. '                    ') then
!  premier fichier
                 read (ligne,'(a)') fone
               else
!  second fichier
               read (ligne,'(a)') ftwo
               endif   
             case ('step')
!   longueur de pas - champs exterieur
               read (ligne(6:128),*) step
               y_ff=.true.
             case('scal')
!   facteur spinorbite
               read (ligne(7:128),*) scale
             case('grou')
!   energie etat fondamental
               read (ligne(8:128),*) ground
             case('prin')
!   impression
               read (ligne(7:128),*) iprint
             case('sect')
!   secteur espace de fock
               ligne(1:7)='       '
               ligne=adjustl(ligne)     
               read (ligne,'(a)') sector
             case('rep:')
               nombre_irreps=nombre_irreps+1
               read (ligne(5:128),*) irrepset(nombre_irreps)
               iw=index(ligne,'main:')
               if(iw >= 1) then
!   sous-espace principal
                 iww=index(ligne,'-')
                 if(iww >= 1) then
                   ligne(iww:iww)=' '
                   read (ligne(iw+6:128),*) nskipset(nombre_irreps),mainmsset(nombre_irreps)
                   nskipset(nombre_irreps)=nskipset(nombre_irreps)-1
                   mainmsset(nombre_irreps)=mainmsset(nombre_irreps)-nskipset(nombre_irreps)
                 else
                   read (ligne(iw+6:128),*) mainmsset(nombre_irreps)
                 endif
               endif
             case('main')
!
!    archaique mais soit, tuerons plus tard   
!
               iww=index(ligne,'-')
               if(iww >= 1) then
                 ligne(iww:iww)=' '
                 if(nombre_irreps > 0)then
                   read (ligne(iw+6:128),*) nskipset(nombre_irreps),mainmsset(nombre_irreps)
                   nskipset(nombre_irreps)=nskipset(nombre_irreps)-1
                   mainmsset(nombre_irreps)=mainmsset(nombre_irreps)-nskipset(nombre_irreps)
                 else
                   read (ligne(6:128),*) nskipset(1),mainmsset(1)
                     nskipset(1)=nskipset(1)-1
                     mainmsset(1)=mainmsset(1)-nskipset(1)
                 endif
               else
                 if(nombre_irreps > 0)then
                   read (ligne(iw+6:128),*) mainmsset(nombre_irreps)
                 else
                   read (ligne(6:128),*) mainmsset(1)
                 endif
               endif
!
!    fin du morceau caduc
!
             case ('arra')
!   rearranger la serie des hamiltoniens
               arrange=ligne(iw+9:iw+88)
             case('Bloc')
!   formulation strictement non-hermitique
               y_dc=.false.
             case default
               write(*,*)' '
               write(*,*)'    warning: i do not understand the record -> ',ligne 
           end select 
                     
         enddo lecture
!
!
!   fin lecture parametres
!
!   consistency check
!
      if( (.not.y_ff) .and. ( y_dc .and. (.not. y_cut) ) ) then
        if(ftwo(1:4).ne.'    ') then
          y_dc=.false.
          write(*,'(a,a)')'    warning: no need to use heff des cloizeaux, switched off '
        endif
      endif
!      
      if( (.not.y_dc) .and. y_cut) then
      y_dc=.true.
        write(*,'(a,a)')'    warning: hermitian heffs assumed for trimmed series extrapolation '
      endif
!      
!      
      if(y_ff)then
!
!     finite-field calculations declared
!      
         if(ftwo(1:4) /= '    ')write(*,'(a)')'    finite-field transition dipole calculations'
         write(*,'(a,a)')'    effective hamiltonian 1 in ',trim(fone),   &
                         '    effective hamiltonian 2 in ',trim(ftwo)
      else
!
!     quasidiabatization via projection declared
!      
         if(ftwo(1:4).ne.'    ') then
         write(*,'(a)')'    effective hamiltonian projection'
         write(*,'(a,a)')'    spin-orbit-free heff    in ',trim(fone),    & 
                         '    spin-orbit-coupled heff in ',trim(ftwo)
         else
         write(*,'(a,a)')'    eff hamiltonian(s) file    ',trim(fone)
         endif
      endif
!      
      if(y_ff) then
         if(ftwo(1:4).ne.'    ') write(*,'(/a,f10.7)')'    finite-field step         ', step
      else
         if(ftwo(1:4).ne.'    ') write(*,'(/a,f10.7)')'    so fraction in so-free    ', scale
      endif   
      write(*,'(a,a)')  '    sector                          ', sector 
 
 
!     grande boucle irreps ----------------------------------------            
!     grande boucle irreps ---------------------------------------- 
!     grande boucle irreps ----------------------------------------            
      
reps: do numero_irrep=1,nombre_irreps  
      
        irrep= irrepset(numero_irrep)         
        mainms=mainmsset(numero_irrep)
        nskip= nskipset(numero_irrep)
!       
!
!
      inquire(unit=11,opened=yyy)      
      if(yyy) then
        rewind(unit=11)
      else
        open(unit=11,file=trim(fone),form='FORMATTED',status='OLD')
      endif
!      
!      
!      
!      
      if(ftwo(1:20) .ne. '                    ') then
        inquire(unit=12,opened=yyy)      
        if(yyy) then
          rewind(unit=12)
        else
          open(unit=12,file=trim(ftwo),form='FORMATTED',status='OLD')
        endif
      endif
!
!   HEFFF file prescanning, determine arithmetic type, find Heff size
!
      ysector=.false.
      not_end=0
lire1:do while (not_end == 0)
        read(11,'(a)',iostat=not_end) ligne(1:80)
        if (not_end /= 0) then
! fin du fichier
          exit
        endif
        iw=index(ligne(1:80),'arithmetic')
        if(iw >= 1) then
          iww=index(ligne(1:80),'comp')
          if (iww.ge.1) carith=.true.
          cycle
        endif  
        iw=index(ligne(1:80),'sector')
        if(iw >= 1) then
          read (ligne(1:80),'(a)') sector2
!
          if ( (sector2 == '0h0p') .and.(ftwo(1:20) == '                    ') ) then
            read(11,*)
            read(11,*)energie0h0pone
            cycle
          endif
!
          if (sector == sector2 .and. (.not.ysector)) ysector=.true.
          if (sector /= sector2 .and. ysector) ysector=.false.
          cycle
        endif  
        if ( .not. ysector) then
          cycle
        endif
        iw=index(ligne(1:80),'heff size') 
        if(iw.ge.1 .and. ysector) then
          read (ligne(1:80),*) mdummy,mmmm
          if(mdummy.eq.irrep) maxmax=mmmm
        endif  
      enddo lire1 
!
!  
      if( mainms > maxmax ) mainms=maxmax
      if(numero_irrep == 1) then
      write(*,'(a,l10)')'    complex arithmetic        ', carith 
      write(*,'(a,l10)')'    heff hermitization        ', y_dc 
      if(llmmpade > 0) then
        write(*,'(a,i1,a1,i1,a1)')'    pade extrapolation  (diff)     [',llpade,'/',mmpade,']'
      endif 
      if(llmmpade < 0) then
        write(*,'(a,i1,a1,i1,a1)')'    pade extrapolation  (full)     [',llpade,'/',mmpade,']'
      endif 
      if(llmmpade.ne.0) write(*,'(a,l8)')'    heff pre-projection         ', y_cut
      if(abs(energie0h0pone) > 1.d-9) write(*,'(a,f16.8)')'    0h0p energy         ', & 
                                                              energie0h0pone
      endif
      write(*,'(a)')'    '      
      write(*,'(a,i10)')'    rep no.                   ', irrep
      write(*,'(a,i10)')'    heff size                 ', maxmax  
      write(*,'(a,i2,a,i4)')'    target states               ', &
                             nskip+1,' -', mainms+nskip  
      write(*,'(a)')'    '      
      rewind(11)
!
!    space allocation
!      
      allocate(eone(maxmax))  
      allocate(etwo(maxmax))  
      allocate(voneg(maxmax,maxmax))  
      allocate(voned(maxmax,maxmax))  
      allocate(vtwog(maxmax,maxmax))  
      allocate(vtwod(maxmax,maxmax)) 
      allocate(heffectif(mainms,mainms))
!
!      get heffs and compute left and right eigenvectors
!        
      call readheff(carith,11,maxmax,mainms,nskip,sector,      &
                    irrep,eone,voneg,voned,llmmpade,y_cut,y_dc,arrange,iprint)
      if(ftwo(1:20) .eq. '                    ') then
!
!     single point energy extraction only
!
        write(6,'(/12x,a,a/)') '  eigenenergies, au   relative energies, cm-1 ',   &
                               ' imag part (if any), au'
        if(abs(ground-666.d+0) > 1.d-5) then
              thelevel=ground
        else
              thelevel=real(eone(1),kind=pd)
        endif        
        do iii=1,mainms
          if( abs(aimag(eone(iii))) .lt. 1.d-7) then
            write(6,'(10x,i3,f18.8,f16.2)') iii, real(eone(iii),kind=pd)+energie0h0pone,  &
                                            real((eone(iii)-thelevel)*au2cm,kind=pd)        
          else
            write(6,'(10x,i3,f18.8,f16.2,f22.8)')iii,real(eone(iii),kind=pd)+energie0h0pone,  &
                                                 real((eone(iii)-thelevel)*au2cm,kind=pd),aimag(eone(iii))              
          endif       
        enddo   
        write(6,'(a)')' '
!
!      on emprime et enregistre sur un fichier formatte
!      les vecteurs propres de l'hamiltonien effectif
!
        write(6,*) '    effective hamiltonian eigenvectors'
        call imprec(voned,maxmax,mainms,4)
!
        open(16,form='formatted',file='HEFFEVF')
        write(16,'(a,a)')sector, '             # sector' 
        write(16,'(i4,i6,i5,a)')irrep,maxmax,mainms,   &
                                '  # rep No     det basis size    number of states'
        write(16,'(a,l1,a)')'dC=',y_dc,'             # symmetrized Heff ?'
        write(16,'(a)')'eigenvectors'
        write(16,'(4e21.12)')((voned(iii,jjj),iii=1,maxmax),jjj=1,mainms)
        write(16,'(a)')'eigenvalues'
        write(16,'(4e21.12)') (eone(jjj),jjj=1,mainms)
        close(16,status='keep')
        call spectrum(sector,irrep,mainms,eone)              
        goto 13013      
!
!     end of single effective hamiltonian calculations
!
      endif
!
!     read second effective hamiltonian
!      
      call readheff(carith,12,maxmax,mainms,nskip,sector,                 &
                    irrep,etwo,vtwog,vtwod,llmmpade,y_cut,y_dc,arrange,iprint)
      if(.not. y_ff) then
!
!  <<<<<<<<<<<< debut calcul so-versus-sf <<<<<<<<<<<<<<<<<<<<<<<<
!
        write(6,'(/15x,a/15x,a/)') '  eigenenergies, au      relative energies, cm-1 ', &
                                 '  so-free    so-coupled    so-free    so-coupled '
        if(abs(ground-666.d+0) > 1.d-5) then
          thelevel=ground
        else
          thelevel=real(etwo(1),kind=pd)
      endif        
      do iii=1,mainms
      write(6,'(10x,i3,2f12.7,2f12.2)')iii, real(eone(iii),kind=pd),real(etwo(iii),kind=pd),  &
            real((eone(iii)-eone(1))*au2cm,kind=pd),(real(etwo(iii),kind=pd)-thelevel)*au2cm               
      enddo   
      write(6,'(a)')' '
!
!      ---- aller ----
!
      y_socoupled=.true.
!
!     orthogonalisation loewdin de vecteurs avec so
!       
      call loeworth(vtwod,vtwog,maxmax,mainms,' so-coupled vecs',1.1d+0,iprint)
!      
!     les vecteurs so orthogonaux sont dans vtwog
!
      if(iprint > 0) write(6,*)'   projecting so-free vectors onto so-coupled vector subspace'     
!
!    hamiltonien effectif espace so-coupled
!
      call projinto(voned,vtwog,etwo,heffectif,ground,maxmax,mainms,y_socoupled,iprint)
      do iii=1,mainms
      heffectif(iii,iii)=( heffectif(iii,iii)-eone(iii) )*au2cm
        if (iii > 1) then
          do jjj=1,iii-1
          heffectif(iii,jjj)=heffectif(iii,jjj)*au2cm/(1.d0-scale)
          heffectif(jjj,iii)=heffectif(jjj,iii)*au2cm/(1.d0-scale)
          enddo
        endif
      enddo
      write(6,'(/15x,a/15x,a)') '         effective spin-orbit interactions, cm-1          ', &
                                ' basis: so-free states projected onto so-coupled subspace '                  

      do iii = 1, mainms
          do jjj = iii + 1, mainms
              write(6,'(a3,2i3,f16.6)') '@SO', iii, jjj, real(heffectif(iii,jjj))
          end do
      end do
      do iii = 1, mainms
          write(6,'(a7,2i3,f16.6)') '@SODIAG', iii, iii, real(heffectif(iii,iii))
      end do

      call imprec(heffectif,mainms,mainms,6)
!      
!      
!      ---- retour ----
!
      y_socoupled=.false.
!
!     orthogonalisation loewdin de vecteurs so-free
!       
      call loeworth(voned,voneg,maxmax,mainms,' so-free vecs   ',1.1d+0,iprint)
!      
!     les vecteurs so-free orthogonaux sont dans voneg
!
      write(6,'(a)')'        composition of so-coupled states   '
      do icompo=1,mainms
        write(6,'(a,i4,a,f14.8,a)')'     ------- so-coupled state',icompo,  &
                                   '  e= ',real(etwo(icompo),kind=pd),' ---'
        do jcompo=1,maxmax
          tmoment=(0.d+0,0.d+0)
          do iconf=1,maxmax
            tmoment=tmoment+conjg(vtwog(iconf,icompo))*voneg(iconf,jcompo)      
          enddo
          wei=abs(tmoment)**2
          if(wei.gt.0.01d+0)write(6,'(a,i4,a,f14.8,a,f7.4)')          &
          '     so-free state',jcompo,' (e=',real(eone(jcompo),kind=pd),      &
          '): weight ',wei
        enddo
      enddo
!
!    hamiltonien effectif espace so-coupled
!
      call projinto(vtwod,voneg,etwo,heffectif,ground,maxmax,mainms,y_socoupled,iprint)
      do iii=1,mainms
      heffectif(iii,iii)=( heffectif(iii,iii)-eone(iii) )*au2cm
      if (iii.gt.1) then
        do jjj=1,iii-1
          heffectif(iii,jjj)=heffectif(iii,jjj)*au2cm/(1.d0-scale)
          heffectif(jjj,iii)=heffectif(jjj,iii)*au2cm/(1.d0-scale)
        enddo
      endif
      enddo
      write(6,'(/15x,a/15x,a)') '         effective spin-orbit interactions, cm-1          ', &
                                '               basis: so-free states  '                  
      call imprec(heffectif,mainms,mainms,6)
!      
      close(12)
      goto 13013
!
!   >>>>>>>>>>>>> fin calcul so-versus-sf >>>>>>>>>>>>>>>>>>>>>> 
!
      endif
!
!   <<<<<<<<<<<<<< ff tdm  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      
!      
      allocate(tdm1(mainms,mainms))  
      allocate(tdm2(mainms,mainms)) 
!
      write(6,'(/15x,a/15x,a/)') '  eigenenergies, au      relative energies, cm-1 ',  &
                                 '  left file right file    left file  right file '
      if(abs(ground-666.d+0) .gt. 1.d-5) then
              thelevel=ground
      else
              thelevel=0.5d0*real(etwo(1)+eone(1),kind=pd)
      endif        
      do iii=1,mainms
      write(6,'(10x,i3,2f12.7,2f12.2)')iii, real(eone(iii),kind=pd),real(etwo(iii),kind=pd),    &
                                       real((eone(iii)-thelevel)*au2cm,kind=pd),                &
                                       real((etwo(iii)-thelevel)*au2cm,kind=pd)               
      enddo   
      write(6,'(a)')' '
!
!      
      call  scapro(voneg,vtwod,tdm1,maxmax,mainms)
      if ( iprint >= 2) then
         write(6,'(a)') '   overlap matrix  '
         call  imprec(tdm1,mainms,mainms,5)
      endif
! 
      write(6,fmt9998) ' transition  energy, cm^-1    Re(d)      Im(d)   ',     &
                       '     |d|^2       |d|       osc str   A,10^6 s-1'
      do iii = 1, mainms
        do jjj=iii+1,mainms
        tenergy=(real(eone(jjj),kind=pd)-real(eone(iii),kind=pd)+real(etwo(jjj),kind=pd)  &
                -real(etwo(iii),kind=pd))*half
        if(abs(tenergy*au2cm) > ethr) then
          tmoment=(tenergy/step)*tdm1(iii,jjj)
          stmoment=abs(tmoment)*abs(tmoment)
          write(6,fmt9999)iii,'  ->',jjj,tenergy*au2cm, real(tmoment,kind=pd),aimag(tmoment),  &
                          stmoment,abs(tmoment), twothirds*tenergy*stmoment,  &
                          a_dau_eau * stmoment * tenergy**3
          if(.not.y_dc) then     
            tmoment=(tenergy/step)*tdm1(jjj,iii)
            stmoment=abs(tmoment)*abs(tmoment)
            write(6,fmt9999)jjj,'  ->',iii,tenergy*au2cm,           &  
                            real(tmoment,kind=pd),aimag(tmoment),   &  
                            stmoment,abs(tmoment),                  &  
                            twothirds*tenergy*stmoment,             &  
                            a_dau_eau * stmoment * tenergy**3
            tmoment=(tenergy/step)* sqrt(abs(tdm1(iii,jjj))*abs(tdm1(jjj,iii)))
            stmoment=abs(tmoment)*abs(tmoment)
            write(6,fmt9999s)iii,' <=>',jjj,tenergy*au2cm, stmoment,abs(tmoment),    &
            twothirds*tenergy*stmoment,  a_dau_eau * stmoment * tenergy**3
          endif     
        endif
        enddo
      enddo
      write(*,'(a)')'    '
!
!     transitions de droite a gauche
!     que pour tester la stabilite numerique 
!
      if(iprint > 1) then
!
      call  scapro(vtwog,voned,tdm2,maxmax,mainms)
      if ( iprint.ge.2) then
         write(6,'(a)') '   overlap matrix  '
         call  imprec(tdm2,mainms,mainms,5)
      endif
      write(6,fmt9998) ' transition  energy, cm^-1    Re(d)      Im(d)   ',      &
                       '     |d|^2       |d|       osc str   A,10^6 s-1'
      do iii = 1, mainms
        do jjj=iii+1,mainms
          tenergy=(real(eone(jjj))-real(eone(iii),kind=pd)         &
                 +  real(etwo(jjj),kind=pd)-real(etwo(iii),kind=pd))*half
        if(abs(tenergy*au2cm) .gt. ethr) then
          tmoment=(tenergy/step)*tdm2(iii,jjj)
          stmoment=abs(tmoment)*abs(tmoment)
          write(6,fmt9999)jjj,' -> ',iii,tenergy*au2cm,                  &
                          real(tmoment,kind=pd),aimag(tmoment),          &
                          stmoment,abs(tmoment),                         &
                          twothirds*tenergy*stmoment,                    &
                          a_dau_eau * stmoment * tenergy**3
          if(.not.y_dc) then     
            tmoment=(tenergy/step)*tdm2(jjj,iii)
            stmoment=abs(tmoment)*abs(tmoment)
            write(6,fmt9999)iii,' -> ',jjj,tenergy*au2cm,                &
                            real(tmoment,kind=pd),aimag(tmoment),        &
                            stmoment,abs(tmoment),                       &
                            twothirds*tenergy*stmoment,                  &
                            a_dau_eau * stmoment * tenergy**3
          endif     
        endif
        enddo
      enddo
      write(*,'(a)')'    '
!
      endif
!
!     fin recalc gauche - droite            
!
!     bye
!     
      deallocate(tdm1,tdm2)
!      close(12) 
13013 continue   
      deallocate(eone,etwo,voneg,voned,vtwog,vtwod,heffectif)
!     fin grande boucle reps
      enddo reps     
      inquire(unit=12,opened=yyy)      
      if(yyy) close(12)
      close(11)
      stop 
      end 
