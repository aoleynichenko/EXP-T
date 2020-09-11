c     
      program heff
c 
c     
c  
      implicit real*8(a-h,o-z)
      character(len=128)fone, ftwo, ligne      
      character*80 arrange
      character*4 sector,sector2
      complex*16, allocatable :: eone(:),etwo(:)
      complex*16, allocatable :: voneg(:,:),vtwog(:,:),
     &      voned(:,:),vtwod(:,:),
     &      tdm1(:,:),tdm2(:,:),heffectif(:,:)
      complex*16 tmoment
      logical carith,ysector,y_ff,y_socoupled,y_cut,y_dc
      character(len=50), parameter :: fmt9999 = 
     &      '(1x,i3,a4,i3,f14.3,f12.6,f11.6,f12.6,f12.8,2f12.6)'
      character(len=50), parameter :: fmt9999s = 
     &      '(1x,i3,a4,i3,f14.3,23x,        f12.6,f12.8,2f12.6)'
c9998  format(/,a,a,/) 
      character(len=10), parameter :: fmt9998 = '(/,a,a,/)'   
c
c     constants
c
      a_dau_eau=21420.007784d0 
c         (luuk visscher's value)        
      au2cm=219474.631280634d0
c
c     default 
c
      iprint=1
      mainms=99999
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
      nskip=0
      llmmpade=0
      arrange(1:40)='                                        '
      arrange(41:80)='                                        '
      energie0h0pone=0.d+00
      energie0h0ptwo=0.d+00
c
c     chapeau
c
      write(*,'(/3a)')'    ----',
     &' effective hamiltonian manip program',
     &' ----'
      write(*,'(3a/)')'    ----',
     &'         le 28 juin 2020, az        ',
     &' ----'
c
c     parse input file
c
         not_end=0
         do while (not_end == 0)
         read(*,'(a)',iostat=not_end) ligne
         if(not_end .lt. 0) then
           exit
         elseif (not_end .gt. 0) then
           write(*,'(/3a)')'    ----',
     &      ' i cannot read the input stream     ',
     &      ' ----'
           stop
         endif
c
c commentaire
c
         if( ligne(1:1) .eq. '#' ) cycle
         iw=index(ligne,'heff:')
         if(iw.ge.1) then
           y_ff=.false.
           cycle
         endif  
         iw=index(ligne,'file:') 
         if(iw.ge.1) then
c 
c  fichier heff 
c       
           iwp6=iw+5
c
c  supprimer les blancs
c
           do while ( iwp6.lt.80 .and. ligne(iwp6:iwp6) .eq. ' ')
             iwp6=iwp6+1
           enddo
           if(fone(1:20) .eq. '                    ') then
c
c  premier fichier
c
              read (ligne(iwp6:128),'(a)') fone
           else
c
c  seconde fichier
c
              read (ligne(iwp6:128),'(a)') ftwo
           endif   
           cycle
         endif  
         iw=index(ligne,'step:')
         if(iw.ge.1) then
c
c   pas champs exterieur
c
           read (ligne(iw+6:128),*) step
           y_ff=.true.
           cycle
         endif  
         iw=index(ligne,'scale:')
         if(iw.ge.1) then
c
c   facteur spinorbite
c
           read (ligne(iw+7:128),*) scale
           cycle
         endif  
         iw=index(ligne,'ground:')
         if(iw.ge.1) then
c
c   energie etat fondamental
c
           read (ligne(iw+8:128),*) ground
           cycle
         endif  
         iw=index(ligne,'print:')
         if(iw.ge.1) then
c
c   impression
c
           read (ligne(iw+7:128),*) iprint
           cycle
         endif  
         iw=index(ligne,'sector:')
         if(iw.ge.1) then
c
c   secteur espace de fock
c
           do while ( iw.lt.80 .and. ligne(iw:iw) .eq. ' ')
             iw=iw+1
           enddo
           read (ligne(iw+8:128),'(a)') sector
           cycle
         endif           
         iw=index(ligne,'rep:')
         if(iw.ge.1) then
c
c   (co)rep
c
           read (ligne(iw+5:128),*) irrep
           cycle
         endif           
         iw=index(ligne,'main:')
         if(iw.ge.1) then
c
c   sous-espace principal
c
           iww=index(ligne,'-')
           if(iww.ge.1) then
              ligne(iww:iww)=' '
              read (ligne(iw+6:128),*) nskip, mainms
              nskip=nskip-1
              mainms=mainms-nskip
           else
              read (ligne(iw+6:128),*) mainms
           endif
           cycle
         endif
         iw=index(ligne,'arrange:')
         if(iw.ge.1) then
c
c   rearranger la serie des hamiltoniens
c
           arrange=ligne(iw+9:iw+88)
           cycle
         endif
         iw=index(ligne,'Bloch')
         if(iw.ge.1) then
c
c   formulation strictement non-hermitique
c
           y_dc=.false.
           cycle
         endif           
         iw=index(ligne,']')
         iww=index(ligne,'[')
         iwww=index(ligne,'/')
         iwwww=index(ligne,'f')
         iwwwww=index(ligne,'cut')
         if(iw.ge.1 .and. iww.ge.1 .and. iwww.ge.1) then
c
c   pade
c
           ligne(iw:iw)=' '
           ligne(iww:iww)=' '
           ligne(iwww:iwww)=' '
           if(iwwww.ge.1)ligne(iwwww:iwwww)=' '
           if(iwwwww.ge.1)then
             ligne(iwwwww:iwwwww+2)='   '
             y_cut=.true.
           endif
           read (ligne(1:128),*) llpade,mmpade
           llmmpade=10*llpade+mmpade
           if(iwwww.ge.1)llmmpade=-llmmpade
c
c           
           cycle
         endif           
         enddo
c
c
c   fin lecture parametres
c
c   consistency check
c
      if( (.not.y_ff) .and. ( y_dc .and. (.not. y_cut) ) ) then
        if(ftwo(1:4).ne.'    ') then
          y_dc=.false.
          write(*,'(a,a)')
     *    '    warning: no need to use heff des cloizeaux,',
     *    ' switched off '
        endif
      endif
c      
      if( (.not.y_dc) .and. y_cut) then
      y_dc=.true.
        write(*,'(a,a)')'    warning: heff des cloizeaux assumed for',
     *  ' trimmed series extrapolation '
      endif
c      
c      
      if(y_ff)then
c
c     finite-field calculations declared
c      
         if(ftwo(1:4).ne.'    ')
     &   write(*,'(a)')'    finite-field transition dipole calculations'
         write(*,'(a,a)')'    effective hamiltonian 1 in ',fone(1:60),
     &                '    effective hamiltonian 2 in ',ftwo(1:60)
      else
c
c     quasidiabatization via projection declared
c      
         if(ftwo(1:4).ne.'    ') then
         write(*,'(a)')'    effective hamiltonian projection'
         write(*,'(a,a)')'    spin-orbit-free heff    in ',fone(1:60),
     &                '    spin-orbit-coupled heff in ',ftwo(1:60)
    
         else
         write(*,'(a,a)')'    eff hamiltonian(s) file    ',fone(1:60)
         endif
      endif
      if(y_ff) then
         if(ftwo(1:4).ne.'    ')
     &    write(*,'(/a,f10.7)')'    finite-field step         ', step
      else
         if(ftwo(1:4).ne.'    ')
     &   write(*,'(/a,f10.7)')'    so fraction in so-free    ', scale
      endif    
      write(*,'(a,a)')  '    sector                          ', sector 
      write(*,'(a,i10)')'    rep no.                   ', irrep
      open(unit=11,file=fone,form='FORMATTED',status='OLD')

c      
c      
c      
c      open(unit=66,file='kar_so',form='FORMATTED')
c      
c      
      

      if(ftwo(1:20) .ne. '                    ')
     &    open(unit=12,file=ftwo,form='FORMATTED',status='OLD')
c
c   HEFFF file prescanning, determine arithmetic type, find Heff size
c
      ysector=.false.
      not_end=0
      do while (not_end == 0)
         read(11,'(a)',iostat=not_end) ligne(1:80)
         if (not_end .ne. 0) then
c fin du fichier
           exit
         endif
         iw=index(ligne(1:80),'arithmetic')
         if(iw.ge.1) then
           iww=index(ligne(1:80),'comp')
           if (iww.ge.1) carith=.true.
           cycle
         endif  
         iw=index(ligne(1:80),'sector')
         if(iw.ge.1) then
           read (ligne(1:80),'(a)') sector2
c
           if ( (sector2 .eq. '0h0p') .and.
     &      (ftwo(1:20) .eq. '                    ') ) then
           read(11,*)
           read(11,*)energie0h0pone
           cycle
           endif
c
           if (sector.eq.sector2.and. (.not.ysector)) ysector=.true.
           if (sector.ne.sector2.and. ysector) ysector=.false.
           cycle
         endif  
         if(.not. ysector)then
           cycle
         endif
         iw=index(ligne(1:80),'heff size') 
         if(iw.ge.1 .and. ysector) then
           read (ligne(1:80),*) mdummy,mmmm
           if(mdummy.eq.irrep) maxmax=mmmm
         endif  
      enddo 
c
c  
      if( mainms. gt. maxmax ) mainms=maxmax
      write(*,'(a,i10)')'    heff size                 ', maxmax  
      write(*,'(a,i2,a,i4)')'    target states               ',
     &                        nskip+1,' -', mainms+nskip  
      write(*,'(a,l10)')'    complex arithmetic        ', carith 
      write(*,'(a,l10)')'    heff des cloizeaux        ', y_dc 
      if(llmmpade.gt.0) then
        write(*,'(a,i1,a1,i1,a1)')
     &   '    pade extrapolation  (diff)     [',llpade,'/',mmpade,']'
      endif 
      if(llmmpade.lt.0) then
        write(*,'(a,i1,a1,i1,a1)')
     &   '    pade extrapolation  (full)     [',llpade,'/',mmpade,']'
      endif 
      if(llmmpade.ne.0)
     &      write(*,'(a,l8)')'    heff pre-projection         ', y_cut
      if(abs(energie0h0pone) .gt. 1.d-9) 
     &      write(*,'(a,f16.8)')'    0h0p energy         ', 
     &      energie0h0pone
      write(*,'(a)')'    '      
      rewind(11)
c
c    space allocation
c      
      allocate(eone(maxmax))  
      allocate(etwo(maxmax))  
      allocate(voneg(maxmax,maxmax))  
      allocate(voned(maxmax,maxmax))  
      allocate(vtwog(maxmax,maxmax))  
      allocate(vtwod(maxmax,maxmax)) 
      allocate(heffectif(mainms,mainms))
c
c      get heffs and compute left and right eigenvectors
c        
      call readheff(carith,11,maxmax,mainms,nskip,sector,
     &  irrep,eone,voneg,voned,llmmpade,y_cut,y_dc,arrange,iprint)
      if(ftwo(1:20) .eq. '                    ') then
c
c     single point energy extraction only
c
        write(6,'(/12x,a,a/)')
     &     '  eigenenergies, au   relative energies, cm-1 ',
     &     ' imag part (if any), au'
        if(abs(ground-666.d+0) .gt. 1.d-5) then
              thelevel=ground
        else
              thelevel=real(eone(1))
        endif        
        do iii=1,mainms
          if( abs(aimag(eone(iii))) .lt. 1.d-7) then
            write(6,'(10x,i3,f18.8,f16.2)')iii,
     &      real(eone(iii))+energie0h0pone,
     &      real((eone(iii)-thelevel)*au2cm)        
          else
            write(6,'(10x,i3,f18.8,f16.2,f22.8)')iii,
     &      real(eone(iii))+energie0h0pone,
     &      real((eone(iii)-thelevel)*au2cm),aimag(eone(iii))              
          endif       
        enddo   
        write(6,'(a)')' '
c
c      on emprime et enregistre sur un fichier formatte
c      les vecteurs propres de l'hamiltonien effectif
c
        write(6,*) ' effective hamiltonian eigenvectors'
        call primac(voned,maxmax,maxmax,mainms,5.d-5) 
c
        open(16,form='formatted',file='HEFFEVF')
        write(16,'(a,a)')sector, '             # sector' 
        write(16,'(i4,i6,i5,a)')irrep,maxmax,mainms,
     &  '  # rep No     det basis size    number of states'
        write(16,'(a,l1,a)')'dC=',y_dc,
     &     '             # symmetrized Heff ?'
        write(16,'(a)')'eigenvectors'
        write(16,'(4e21.12)')
     &   ((voned(iwrite,jwrite),iwrite=1,maxmax),jwrite=1,mainms)
        write(16,'(a)')'eigenvalues'
        write(16,'(4e21.12)') (eone(jwrite),jwrite=1,mainms)
        close(16,status='keep')
        call spectrum(sector,irrep,mainms,eone)              
        goto 13013      
c
c     end of single effective hamiltonian calculations
c
      endif
c
c     read second effective hamiltonian
c      
      call readheff(carith,12,maxmax,mainms,nskip,sector,
     &  irrep,etwo,vtwog,vtwod,llmmpade,y_cut,y_dc,arrange,iprint)
      if(.not. y_ff) then
*
*  <<<<<<<<<<<< debut calcul so-versus-sf <<<<<<<<<<<<<<<<<<<<<<<<
*
      write(6,'(/15x,a/15x,a/)')
     &   '  eigenenergies, au      relative energies, cm-1 ',
     &   '  so-free    so-coupled    so-free    so-coupled '
      if(abs(ground-666.d+0) .gt. 1.d-5) then
              thelevel=ground
      else
              thelevel=real(etwo(1))
      endif        
      do iii=1,mainms
      write(6,'(10x,i3,2f12.7,2f12.2)')iii,
     &  real(eone(iii)),real(etwo(iii)),
     &  real((eone(iii)-eone(1))*au2cm),
     &  real((etwo(iii)-thelevel)*au2cm)               
      enddo   
      write(6,'(a)')' '
c
c      ---- aller ----
c
      y_socoupled=.true.
c
c     orthogonalisation loewdin de vecteurs avec so
c       
      call loeworth(vtwod,vtwog,maxmax,mainms,
     &' so-coupled vecs',1.1d+0,iprint)
c      
c     les vecteurs so orthogonaux sont dans vtwog
c
      if(iprint.gt.0) write(6,*)'   projecting so-free vectors ',
     &          'onto so-coupled vector subspace'     
c
c    hamiltonien effectif espace so-coupled
c
      call projinto(voned,vtwog,etwo,heffectif,ground,
     &           maxmax,mainms,y_socoupled,iprint)
      do iii=1,mainms
      heffectif(iii,iii)=( heffectif(iii,iii)-etwo(1  ) )*au2cm
      if (iii.gt.1) then
        do jjj=1,iii-1
        heffectif(iii,jjj)=heffectif(iii,jjj)*au2cm/(1.d0-scale)
        heffectif(jjj,iii)=heffectif(jjj,iii)*au2cm/(1.d0-scale)
        enddo
      endif
      enddo
      write(6,'(/15x,a/15x,a)')
     &   '         effective spin-orbit interactions, cm-1          ',
     &   ' basis: so-free states projected onto so-coupled subspace '                  
      call primaclong(heffectif,mainms,mainms,mainms,0.05d0)


c      
c      
c      write(66,'(6f12.3)')real(heffectif(2,1)),
c     &  real(heffectif(3,1)),real(heffectif(3,2)),
c     &  real(heffectif(1,1)),real(heffectif(2,2)),
c     &  real(heffectif(3,3))                     
c      close(66,status='keep')
c      

*      
c
c      ---- retour ----
c
      y_socoupled=.false.
c
c     orthogonalisation loewdin de vecteurs so-free
c       
      call loeworth(voned,voneg,maxmax,mainms,
     &' so-free vecs   ',1.1d+0,iprint)
c      
c     les vecteurs so-free orthogonaux sont dans voneg
c
      if(iprint.gt.0) then
         write(6,*)'   projecting everything ',
     &          'onto so-free vector subspace'  
      endif   
c
c
      write(6,'(a)')'        composition of so-coupled states   '
      do icompo=1,mainms
        write(6,'(a,i4,a,f14.8,a)')'     ------- so-coupled state',
     &      icompo, '  e= ',real(etwo(icompo)),' ---'
        do jcompo=1,maxmax
          tmoment=(0.d+0,0.d+0)
          do iconf=1,maxmax
c            write(6,*)vtwod(iconf,icompo),voneg(iconf,jcompo)
            tmoment=tmoment+conjg(vtwog(iconf,icompo))*
     &            voneg(iconf,jcompo)      
          enddo
          wei=abs(tmoment)**2
          if(wei.gt.0.01d+0)write(6,'(a,i4,a,f14.8,a,f7.4)')
     &    '     so-free state',jcompo,' (e=',real(eone(jcompo)),
     &    '): weight ',wei
        enddo
      enddo
c
c    hamiltonien effectif espace so-coupled
c
      call projinto(vtwod,voneg,etwo,heffectif,ground,
     &           maxmax,mainms,y_socoupled,iprint)
      do iii=1,mainms
      heffectif(iii,iii)=( heffectif(iii,iii)-eone(iii) )*au2cm
      if (iii.gt.1) then
        do jjj=1,iii-1
          heffectif(iii,jjj)=heffectif(iii,jjj)*au2cm/(1.d0-scale)
          heffectif(jjj,iii)=heffectif(jjj,iii)*au2cm/(1.d0-scale)
        enddo
      endif
      enddo
      write(6,'(/15x,a/15x,a)')
     &   '         effective spin-orbit interactions, cm-1          ',
     &   '               basis: so-free states  '                  
      call primaclong(heffectif,mainms,mainms,mainms,0.05d0)
*      
      close(12)
      goto 13013
*
*   >>>>>>>>>>>>> fin calcul so-versus-sf >>>>>>>>>>>>>>>>>>>>>> 
*
      endif
*
*   <<<<<<<<<<<<<< ff tdm  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      
*      
      allocate(tdm1(mainms,mainms))  
      allocate(tdm2(mainms,mainms)) 
*
      write(6,'(/15x,a/15x,a/)')
     &   '  eigenenergies, au      relative energies, cm-1 ',
     &   '  left file right file    left file  right file '
      if(abs(ground-666.d+0) .gt. 1.d-5) then
              thelevel=ground
      else
              thelevel=0.5d0*real(etwo(1)+eone(1))
      endif        
      do iii=1,mainms
      write(6,'(10x,i3,2f12.7,2f12.2)')iii,
     &  real(eone(iii)),real(etwo(iii)),
     &  real((eone(iii)-thelevel)*au2cm),
     &  real((etwo(iii)-thelevel)*au2cm)               
      enddo   
      write(6,'(a)')' '
c
c      
      call  scapro(voneg,vtwod,tdm1,maxmax,mainms)
      if ( iprint.ge.2) then
         write(6,'(a)') '   overlap matrix  '
         call  primac(tdm1,mainms,mainms,mainms,5.d-5)
      endif
c 
      write(6,fmt9998)
     &  ' transition  energy, cm^-1    Re(d)      Im(d)   ',
     &  '     |d|^2       |d|       osc str   A,10^6 s-1'
      do iii = 1, mainms
        do jjj=iii+1,mainms
        tenergy=(real(eone(jjj))-real(eone(iii))
     &       +  real(etwo(jjj))-real(etwo(iii)))*0.5d0
        if(abs(tenergy*au2cm) .gt. ethr) then
          tmoment=(tenergy/step)*tdm1(iii,jjj)
          stmoment=abs(tmoment)*abs(tmoment)
          write(6,fmt9999)iii,'  ->',jjj,tenergy*au2cm,
     &    real(tmoment),dimag(tmoment),
     &    stmoment,abs(tmoment),
     &    0.6666666666667d+0*tenergy*stmoment,
     &    a_dau_eau * stmoment * tenergy**3
          if(.not.y_dc) then     
            tmoment=(tenergy/step)*tdm1(jjj,iii)
            stmoment=abs(tmoment)*abs(tmoment)
            write(6,fmt9999)jjj,'  ->',iii,tenergy*au2cm,
     &      real(tmoment),dimag(tmoment),
     &      stmoment,abs(tmoment),
     &      0.6666666666667d+0*tenergy*stmoment,
     &      a_dau_eau * stmoment * tenergy**3
            tmoment=(tenergy/step)*
     &        sqrt(abs(tdm1(iii,jjj))*abs(tdm1(jjj,iii)))
            stmoment=abs(tmoment)*abs(tmoment)
            write(6,fmt9999s)iii,' <=>',jjj,tenergy*au2cm,
     &      stmoment,abs(tmoment),
     &      0.6666666666667d+0*tenergy*stmoment,
     &      a_dau_eau * stmoment * tenergy**3
          endif     
        endif
        enddo
      enddo
      write(*,'(a)')'    '
c
c     transitions de droite a gauche
c     que pour tester la stabilite numerique 
c
      if(iprint .ge. 1) then
c
      call  scapro(vtwog,voned,tdm2,maxmax,mainms)
      if ( iprint.ge.2) then
         write(6,'(a)') '   overlap matrix  '
         call  primac(tdm2,mainms,mainms,mainms,5.d-5)
      endif
      write(6,fmt9998)
     &  ' transition  energy, cm^-1    Re(d)      Im(d)   ',
     &  '     |d|^2       |d|       osc str   A,10^6 s-1'
      do iii = 1, mainms
        do jjj=iii+1,mainms
        tenergy=(real(eone(jjj))-real(eone(iii))
     &       +  real(etwo(jjj))-real(etwo(iii)))*0.5d0
        if(abs(tenergy*au2cm) .gt. ethr) then
          tmoment=(tenergy/step)*tdm2(iii,jjj)
          stmoment=abs(tmoment)*abs(tmoment)
          write(6,fmt9999)jjj,' -> ',iii,tenergy*au2cm,
     &    real(tmoment),dimag(tmoment),
     &    stmoment,abs(tmoment),
     &    0.6666666666667d+0*tenergy*stmoment,
     &    a_dau_eau * stmoment * tenergy**3
          if(.not.y_dc) then     
            tmoment=(tenergy/step)*tdm2(jjj,iii)
            stmoment=abs(tmoment)*abs(tmoment)
            write(6,fmt9999)iii,' -> ',jjj,tenergy*au2cm,
     &      real(tmoment),dimag(tmoment),
     &      stmoment,abs(tmoment),
     &      0.6666666666667d+0*tenergy*stmoment,
     &      a_dau_eau * stmoment * tenergy**3
          endif     
        endif
        enddo
      enddo
      write(*,'(a)')'    '
c
      endif
c
c     fin recalc gauche - droite            
c
c     bye
c     
      deallocate(tdm1,tdm2)
      close(12) 
13013 continue   
      deallocate(eone,etwo,voneg,voned,vtwog,vtwod,
     &           heffectif)
      close(11)
      stop 
      end 
