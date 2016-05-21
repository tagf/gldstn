! mains.F main program for gldstn model
! variable depth
! extra outputs for global version

!	program  gldstn
	subroutine  gldstn(first_gldstn, have_TS_atm, sdedy_oc,ts_oc_for_atm,k2)
	USE IFPORT	!module to use TIMEF() function (= elapsed_time)
	include 'var.f90'
!	include 'varAGCM.f90' 
	     real  TS_atm_for_oc(0:72+1,0:72+1), QS_atm_for_oc(0:72+1,0:72+1),&
             PREC_atm_for_oc(0:72+1,0:72+1)  
         common /varsAGCM/  TS_atm_for_oc, QS_atm_for_oc, PREC_atm_for_oc
	real elapsed_time, elapsed_time0
	real avn, avs, sum, sums(8*maxl), tv(maxl,maxi,maxj,maxk), rms
	real sum1, sum2, sum3,sum4,cor, pacsf
	integer k2(0:maxi+1,0:maxj+1)              !GLD GGGGGGGGGG
! to average osc runs
! real ats(maxl,maxi,maxj,maxk), au(3,maxi,maxj,maxk)
! 1   , afn(maxl,maxi,maxj,maxk)

      real opsi(0:maxj,0:maxk), ou(maxj,maxk)
      real opsia(0:maxj,0:maxk), omina, omaxa
      real opsip(0:maxj,0:maxk), ominp, omaxp

	real mean_v, mean_u, windconst, windconst1
	real fricU(maxj),fricV(maxj), windu(0:maxi,0:maxj),&
		windv(0:maxi,0:maxj),height

      real zpsi(0:maxi,0:maxk), zu(maxi,maxk)

      integer nsteps, npstp, iwstp, itstp, iw, icount,&
           i, j, k, l, istep, iterun, isol, isl, natm,iice,sdeyr0,sdedy_oc
      character ans,lout*3,lin*6,ext*3,conv*3,name*10

      logical flat, osc

      real hft(3), hfp(3), hfa(3), phfmax, tv2, tv3
      real psisl(0:maxi,0:maxj,isles), ubisl(2,0:maxi+1,0:maxj,isles)&
          ,erisl(isles+1,isles+1), psibc(2)
!ctest1    ,erisl(isles,isles+1), psibc(2)
!############  GLD GGGGGGGG
	real  ts_oc_for_atm(0:maxi+1,0:maxj+1), arg  !for AGCM
	logical end_of_day, first_gldstn, have_TS_atm
	
	real yatm(1:74,1:46)!for AGCM
	real SCOSZ(maxj),solfor1(maxj),plot(maxi,maxj) !seasonal run
	character id_mnth*2 !=Jn or Jl
      integer istep1
!############

! for repeated runs forward or backward

! open(20,file='output.553R')
! do iterun=10,250,10
! do iterun=240,0,-10
! open(5,file='input.553')
!******** file write tq1 for diagn and surfer
!******         open(555,file=trim(path_results)//'tq1')
!********
!GLD GGGGGGGGG
      sdedy=sdedy_oc
      if (.not. first_gldstn) goto 123
      istep=1    
!GLD GGGGGGGGG
        elapsed_time0 = TIMEF( )
       path_source='..\..\gldstn_F90_1980\source\'
!       path_source='c:\gldstn_F90_1980\source\'
!       path_results='c:\gldstn_F90_1980\results\'
        path_results='..\..\gldstn_F90_1980\results\'
!********  open run parameters file *************
         open(5,file=trim(path_source)//'goinS') !seasonal basic exper.
!         open(5,file=trim(path_source)//'goinMelt') !ice melting exper.

      print*,'nsteps npstp iwstp itstp'
      read(5,*)nsteps,npstp,iwstp,itstp
      ! nsteps=400 !1000 !500000 !stop after nsteps steps
!At multiples of npstp the subroutine diag is called which writes some
!information to unit 6, normally the screen. diag is called again at the
!next step so you can spot any oscillatory instability easily. The second
!time, ie one step after multiples of npstp, mains calculates an average
!rate of change of all dynamic variables (if npstp .ge. 2)
!	npstp=1000
!At multiples of iwstp the subroutine outm writes the out files
!	iwstp=6250 !1000 !500000
!At multiples of itstp the subroutine diag2 writes the out files
!	itstp=50 !1000 !500
      print*,nsteps,npstp,iwstp,itstp
      print*,'new or continuing run ?'
      read(5,'(a1)')ans
!	ans='c'   !continuing run
!	ans='n'   !new run
      print*,ans
!gseto sets up geometry etc variable depth
      call gseto
! EMBM
! subroutine gseta, sets up  atmosphere and sea-ice
      call gseta
      sdedy_oc=sdedy
	  solfor1=solfor  !fix mean year insolation

      print*,'file extension for output (a3) ?'
      read(5,'(a3)')lout
      lout='w29'   !GLD GGGGGGGGG
      print*,lout

      open(4,file=trim(path_results)//lout//'.'//'t')
      write(4,'(11(a11,3x))')'%time      ',' Pac_T_d   ',' Atl_T_d   ',' Ind_T_d   ',' Sou_T_d   ',&
      ' Pac_T_u   ',' Atl_T_u   ',' Ind_T_u   ',' Sou_T_u   ',' drho/dz   ',' speed     '

      open(14,file=trim(path_results)//lout//'.'//'s')
      write(14,'(11(a11,3x))')'%time      ',' Pac_S_d   ',' Atl_S_d   ',' Ind_S_d   ',' Sou_S_d   ',&
      ' Pac_S_u   ',' Atl_S_u   ',' Ind_S_u   ',' Sou_S_u   ',' drho/dz   ',' speed     '

      if(ans.eq.'n'.or.ans.eq.'N')then  !new run
! define an i.c. (initial conditions) done in gset
         read(5,'(a6)')lin
      else
         print*,'input file extension for input (a6)'
         read(5,'(a6)')lin
         lin='w29.3'  !GLD GGGGGGGGG
      print*,lin
         open(1,file=trim(path_results)//lin)

! subroutine inm.f reads in data gldstn
! expanded to read in atmos and sea ice data

         call inm(1)
         close(1)
         sdedy_oc=sdedy
         t0 = 0. !Init time
         t = t0 !Init time
! perturb the salinity
! print*,'perturbing the i.c.'
! do k=10 ,16
! do i=9,11
! do j=1,1
! ts(1, i,j ,k ) = ts(1, i,j ,k ) + 0.2*(1-0.5*abs(
! 1           real(i-10)))
! enddo
! enddo
! enddo
         do k=1,kmax
            do j=1,jmax
               do i=1,imax
                  do l=1,lmax
                     ts1(l,i,j,k) = ts(l,i,j,k)
                  enddo  !rho - water density
                  rho(i,j,k) = ec(1)*ts(1,i,j,k) + ec(2)*ts(2,i,j,k)&
                         + ec(3)*ts(1,i,j,k)**2 + ec(4)*ts(1,i,j,k)**3
               enddo
            enddo
         enddo
! EMBM atm
         do j=1,jmax
            do i=1,imax !temp (1) and specific humidity (2) in atm
               tq1(1,i,j) = tq(1,i,j)
               tq1(2,i,j) = tq(2,i,j)
            enddo
         enddo
! EMBM sea-ice
         do j=1,jmax
            do i=1,imax
	!sea ice variables: average height (1) and fractional area (2)
               varice1(1,i,j) = varice(1,i,j)
               varice1(2,i,j) = varice(2,i,j)
            enddo
         enddo
      endif  !continuing run

      print*,'seasonal or mean solar run ?'
      read(5,'(a1)')ans
      ans='s' !GLD GGGGGGGGG
      print*,ans
      if (ans.eq.'s'.or.ans.eq.'S') then
       print*,'seasonal solar run '
	endif
      if (ans.eq.'m'.or.ans.eq.'M') then
       print*,'mean solar run '
	endif
! periodic b.c. (required for implicit code)

      do k=1,kmax
         do j=1,jmax
            rho(0,j,k) = rho(imax,j,k)
            rho(imax+1,j,k) = rho(1,j,k)
            do l=1,lmax
               ts(l,0,j,k) = ts(l,imax,j,k)
               ts(l,imax+1,j,k) = ts(l,1,j,k)
! for cimp.ne.1 need
               ts1(l,0,j,k) = ts(l,imax,j,k)
               ts1(l,imax+1,j,k) = ts(l,1,j,k)
            enddo
         enddo
      enddo

! oscillating forcing

! if(abs(sda1).gt. 1e-5)then
! osc = .true.
! print*,'oscillatory forcing amplitude',sda1
! do k=1,kmax
! do j=1,jmax
! do i=1,imax
! do l=1,lmax
! ats(l,i,j,k)= 0
! afn(l,i,j,k) = 0
! enddo
! do l=1,3
! au(l,i,j,k) = 0
! enddo
! enddo
! enddo
! enddo
! else
         osc = .false.
! endif

      flat = .true.
      do i=1,imax
         do j=1,jmax
            if(k1(i,j).gt.1.and.k1(i,j).le.kmax)flat = .false.
         enddo
      enddo
      if(flat)then
         print*,'flat bottom'
      else
         print*,'topography present'
      endif
! invert.f subroutine to invert matrix for barotropic streamfunction
      call invert

!********#ifdef disle2
!*******c code for multiple islands
!*********      if(isles.ne.2)stop 'Makefile and var.fi have different no. isles'
!********#endif

      do isol=1,isles

! set source term to 1 on the ith island (i+1th landmass) only

         do j=0,jmax
            do i=1,imax
               k=i + j*imax !gbold(maxi*maxj)- initial value of gb
               if(gbold(k).eq.isol + 1)then
                  gb(k) = 1.0
               else
                  gb(k) = 0.0
               endif
            enddo
         enddo
! ubarsolv.f subroutine to calculate barotropic velocity on c grid
         call ubarsolv(ubisl(1,0,0,isol),psisl(0,0,isol))

! find island path integral due to unit source on boundary

         do isl=1,isles
! island.f subroutine to calculate path integral around island.
            call island(ubisl(1,0,0,isol),erisl(isl,isol),isl,0)
         enddo
      enddo

      print*,'island path integrals due to unit sources',&
               ((erisl(isl,isol),isl=1,isles),isol=1,isles)
! wind.f sets wind stress forcing for barotropic streamfunction
      call wind_mean

      iw = 1
      icount = 0
!##########   for year means
	    i_avr=0
          tq_avr=0.
	    ice_avr=0.
!#################
 123     k2=k1    !globe map
   if (have_TS_atm) then   !GLD GGGGGGG
  ! Must be after interpolation: 
         do j=1,jmax
           do i=1,imax  
            write (144,*) i,j,tq(1,i,j)-TS_atm_for_oc(i,j),k2(i,j) !GGGGGGGGGG
    !        write (145,*) i,j,k2(i,j) !GGGGGGGGGG
           tq(1,i,j)=TS_atm_for_oc(i,j) 
           tq(2,i,j)=QS_atm_for_oc(i,j) 
           pptn(i,j)=PREC_atm_for_oc(i,j) 
           enddo
         enddo
      ! stop 'tq-ts'  !test difference
   endif
 ! time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  do istep1=1,1 !****nsteps
      istep=istep+1
      istepCO2=istepCO2+1

! calculation of atmospheric advective velocities
! wind #####
!	open (39,file='c:\gldstn model\results\result',status='old')
!	write (39,*), 'day', sdedy, 'year', sdeyr
	height = 1000.  !height of atmosphere, m
	do j=1,2
		fricU(j) = 1.e-4      !coefficient of friction
	    fricU(jmax+1-j) = 1.e-4
		fricV(j) = 1.e-4
	    fricV(jmax+1-j) = 1.e-5
	enddo
	do j=3,jmax/2
		fricU(j) = 3.e-4
		fricU(jmax+1-j) = 3.e-4		
		fricV(j) = 1.e-5
		fricV(jmax+1-j) = 1.e-5
	enddo
	do j=jmax/2+1,jmax/2+4
		fricU(j) = 0.05
		fricV(j) = 0.001
	enddo

	do j=1,jmax-1
	  do i=1,imax
		windconst = ((tq(1,i,j)+273.)-gsc/cpa*height)/(tq(1,i,j)+273.)&
		/(tq(1,i,j)+273.)
		windconst1=gsc*height*windconst/rsc/(4.*fsc*fsc*s(j)*s(j)+fricU(j)&
		*fricV(j))
		if (i.eq.1) then
		  windu(i,j)=-windconst1*(fricV(j)*(tq(1,i,j)-tq(1,imax,j))/c(j)&
		  /dphi+2.*c(j)*fsc*s(j)*(tq(1,i,j+1)-tq(1,i,j))/ds)
		  windv(i,j)=windconst1*(2.*fsc*s(j)*(tq(1,i,j)-tq(1,imax,j))/c(j)&
		  /dphi-0.3*fricU(j)*c(j)*(tq(1,i,j+1)-tq(1,i,j))/ds)
		else
		  windu(i,j)=-windconst1*(fricV(j)*(tq(1,i,j)-tq(1,i-1,j))/c(j)&
		  /dphi+2.*c(j)*fsc*s(j)*(tq(1,i,j+1)-tq(1,i,j))/ds)
		  windv(i,j)=windconst1*(2.*fsc*s(j)*(tq(1,i,j)-tq(1,i-1,j))/c(j)&
		  /dphi-0.3*fricU(j)*c(j)*(tq(1,i,j+1)-tq(1,i,j))/ds)
		endif
!c		write(39, 113),'j=',j,'i=',i ,'sin=',s(j),'cos=',c(j)
! 1		,'v=',windv(i,j),'u=',windu(i,j),'T=',tq(1,i,j)
! 2		,'T(j-1)=',tq(1,i,j-1),'T(i-1)=',tq(1,i-1,j)
!		write (39,113),'j=',j,'i=',i,'u=',windu(i,j),'v=',windv(i,j)
113	    format(a2,i2,x,a2,i2,x,a2,e10.4,x,a2,e8.2)
		uatm(1,i,j)=windu(i,j)
		uatm(2,i,j)=windv(i,j)
	  end do
	end do

	
! write out atmospheric advective velocities

! do j=1,jmax
! do i=1,imax
! plot(i,j)=windu(i,j)
! enddo
! enddo

! call DataForSurAtm(lout,'Jl',plot,'uwind')

! write out atmospheric advective velocities

! do j=1,jmax
! do i=1,imax
! plot(i,j)=windv(i,j)
! enddo
! enddo

! call DataForSurAtm(lout,'Jl',plot,'vwind')
	close (39)
! end calculation of atmospheric advective velocities


! if(osc)call oscer
         if(.not.flat.or.osc.or.istep.eq.1)then
! jbar.f subroutine to calculate jbar forcing for streamfunction
            call jbar
! ubarsolv.f subroutine to calculate barotropic velocity on c grid
            call ubarsolv(ub,psi) !ub-barotropic velocity components
                                  !psi(0:maxi,0:maxj)- barotropic streamfunction
! find island path integral due to wind and jbar terms

            do isl=1,isles
! island.f subroutine to calculate path integral around island.
              call island(ub,erisl(isl,isles+1),isl,1)
            enddo

!****#ifdef disle2
! solve system of simultaneous equations. Zero division here might
! suggest not enough islands in the .psiles file

!****            psibc(1) = (- erisl(1,3)*erisl(2,2) + erisl(2,3)*erisl(1,2))
!****     1                  /(erisl(2,2)*erisl(1,1) - erisl(2,1)*erisl(1,2))
!****            psibc(2) =   (erisl(1,3)*erisl(2,1) - erisl(2,3)*erisl(1,1))
!****     1                  /(erisl(2,2)*erisl(1,1) - erisl(2,1)*erisl(1,2))
           if (isles.eq.0)  then !- only one continent
! to artificially set flows around islands to zero set psibc() to zero
            psibc(1) = 0.
            psibc(2) = 0.
           else
!****#else
            do isl=1,isles
             psibc(isl) = - erisl(isl,isl+1)/erisl(isl,isl)
!            psibc(1) = - erisl(1,2)/erisl(1,1)
            enddo
           endif
!****#endif

            do j=1,jmax
               do i=0,imax+1
                  do isl=1,isles  !ub - barotropic velocity components
                     ub(1,i,j) = ub(1,i,j) + ubisl(1,i,j,isl)*psibc(isl)
                     ub(2,i,j) = ub(2,i,j) + ubisl(2,i,j,isl)*psibc(isl)
                  enddo
               enddo
            enddo

! update diagnostic psi, not always necessary

            do j=0,jmax
               do i=0,imax
                  do isl=1,isles
                     psi(i,j) = psi(i,j) + psisl(i,j,isl)*psibc(isl)
                  enddo
               enddo
            enddo
         endif

! update velocities

         call velc

! EMBM change E-P forcing
! extra0 = extra0 + ?? check def'n same

! HYSTERESIS option
! gradual increase/decrease P-E in N.Atlantic
! plus/minus range0 Sv @ 0.05 Sv/1000 yr [after Rahmstorf 1996]
         extra0 = extra0 + (range0/real(nsteps))

! MELTWATER PULSE option
! reset extra0 to zero after 25/50/100 years:
         if(istep.ge.nsteps_extra0) extra0 = 0.0

! EMBM surface fluxes
! EMBM sea ice now updated in surflux to allow adjustment of surface fluxes
         call surflux

!****#ifdef dimpo
! timestep subroutine for program gldstn
         call tstipo(istep)
!****#else
!****         call tstepo(istep)
!****#endif

! EMBM update 1-layer atmosphere
         do natm = 1,ndta
!****#ifdef dimpa
! subroutine tstipa.f atmospheric timestep for gldstn
            call tstipa
!****#else
!****            call tstepa
!****#endif
! if(mod(istep,npstp).lt.1) call diaga

         enddo

         t = istep*dt(kmax) + t0
	   istepT=istep
 !############# seasonal solar
!         UPDATE DAY, EARTH-SUN DISTANCE, SOLAR DECLINATION
!      end_of_day= (MOD(istep,24).EQ.0)
         end_of_day= .true. ! - if dt = 1day!!!!

         if (end_of_day) then
            do j=0,jmax+1
               do i=0,imax+1
                ts_oc_for_atm(i,j)=ts(1,i,j,kmax) !GLD GGGGGGGGG
               enddo
            enddo
	 ! sdedy - day of year
	 ! SDEYR - year
            call SDET_OC(SCOSZ) !seasonal solar
            sdedy_oc=sdedy      !GLD GGGGGGGGG
            if (ans.eq.'s'.or.ans.eq.'S') then
	         solfor=SCOSZ
	        endif
	     endif  !end_of_day
 !#############

         if(mod(istep,npstp).lt.1)then  !npstp=1000
          elapsed_time =  TIMEF( )   !MCLOCK( )

          print*
          print 111 ,istep,t*tsc/86400.,dt(kmax),(elapsed_time-elapsed_time0)*3650./npstp
 111    format (1x,'## step ',I6,'## t(days) ',F10.0,'## dt ',E9.3,'## 10 year run time ', F6.2,'##')
          elapsed_time0=elapsed_time
  !        print*,'step ',istep,' t(days) ',t*tsc/86400.,'  dt ',dt(kmax)
          print*,'sdedy ',sdedy,' SDEYR ',SDEYR
          print*,'psi on islands ',(psibc(isl),isl=1,isles)
            call diag   !print diagnostic values on screen
            if(mod(istep,npstp).eq.0)then
               do k=1,kmax
                  do j=1,jmax
                     do i=1,imax
                        do l=1,lmax
                           tv(l,i,j,k) = ts(l,i,j,k)
                        enddo
                     enddo
                  enddo
               enddo
            else if(mod(istep,npstp).eq.1.and.istep.gt.1)then
! open(7,file=trim(path_results)//'tmp.1')
               rms = 0
               do j=1,jmax
                  do i=1,imax
                     do k=1,kmax
                        do l=1,lmax
                           rms = rms + (tv(l,i,j,k) - ts(l,i,j,k))**2
! change = (tv(l,i,j,k) - ts(l,i,j,k))
! if(j.eq.3.and.k.eq.kmax-1)write(6,*)i,j,k,l,change
! write(7,*)(tv(l,i,j,k) - ts(l,i,j,k))/dt(kmax)
                        enddo
                     enddo
                  enddo
               enddo
! close (7)
               rms = sqrt(rms/lmax/ntot/dt(kmax)/dt(kmax))
               print*,'r.m.s. r.o.c.',rms
            endif 
	      call diaga
         endif  !npstp
         if(mod(istep,iwstp).eq.0)then
	goto 999
            ext=conv(mod(iw,10)) !output (after iwstp steps) file number
            open(2,file=trim(path_results)//lout//'.'//ext)
            rewind 2
	      ! data for restart
            ! write out: ts(l,i,j,k), u(l,i,j,k),tq(l,i,j),varice(l,i,j),
                                  !tice(i,j),t after any iwstp steps
            call outm(2)
            close(2)

            !#########################
            !Interpolation from ts(1,i,j,kmax) to yatm (74,46)
            ts_oc_for_atm(0:imax+1,0:jmax+1)=ts(1,0:imax+1,0:jmax+1,kmax)
            ! call oc_atm(ts_oc_for_atm,yatm,imax,jmax) possibly not necessary !!

            open(57,file=trim(path_results)//lout//'.'//'yoc')

	      do 20 j=1,jmax !write ts as in AGCM
	       arg=amin1(2./36.*(j-0.5)-1.,1.)
            do 22 i=8,imax
             if(k1(i,j).le.kmax)then !ocean and # interior wet points
              write(57,1 ) i-7,5.*asin(arg),ts(1,i,j,kmax)
	         else
              write(57,1 ) i-7,5.*asin(arg)
	         endif
   22       continue  !i
            do 20 i=1,8
             if(k1(i,j).le.kmax)then !ocean and # interior wet points
              write(57,1 ) 29+i,5.*asin(arg),ts(1,i,j,kmax)
	         else
              write(57,1 ) 29+i,5.*asin(arg)
	         endif
   20       continue !i
            close (57) !yoc
            open(58,file=trim(path_results)//lout//'.'//'yatm')

   	      do 21 j=1,46 !write interpolated temp for AGCM
            do 21 i=2,74
             if(abs(yatm(i,j)).lt.1.e-7)then
                     write(58,2 ) i-1,j  !land
		     else     !ocean and # interior wet points
                     write(58,2 ) i-1,j,yatm(i,j)
	         endif

   21     continue  !j,i
           close (58) !yatm

   1       format (1x, i3,',',e12.4,',',e12.4)
   2       format (1x, i3,',',i3,',',e12.4)

!       stop 'interpol'
          !##############################
          ! write oscillating streamfunction
999		continue
            if(.not.flat.or.osc)then
               open(2,file=trim(path_results)//lout//'.psi.'//ext)
               do j=0,jmax
                  do i=0,imax
                     write(2,*)psi(i,j)
                  enddo
               enddo
               close(2)
            endif
            iw = iw + 1
         endif  !iwstp

         if(mod(istep,itstp).eq.0)then
! open(4,file=trim(path_results)//lout//'.'//'t'
! 1          ,access='append')
! for ibm no access parameter and don't close file faster on sg too
            open(4,file=trim(path_results)//lout//'.'//'t')
            open(14,file=trim(path_results)//lout//'.'//'s')
            open(40,file=trim(path_results)//lout//'.'//'opsit')
            call diag2(sums,avn,avs)
            write(4,110)t*tsc/86400.,(sums(i),i=1,8),avn,avs
            write(14,110)t*tsc/86400.,(sums(i),i=9,16),avn,avs

! close(4)

! if(osc)then
! open(8,file=trim(path_results)//lout//'.osi')
! write(8,100)(ts(2,i,20,7   ),i=1,imax)
! write(8,100)( u(2,i, 1,7   ),i=1,imax)
! open(9,file=trim(path_results)//lout//'.osj')
! write(9,100)(ts(2,1 ,j, 7),j=1,jmax)
! open(12,file=trim(path_results)//lout//'.osk'
! 1             ,access='append')
! write(12,100)(ts(1,10,10,k ),k=1,kmax)
! close(12)
! to average osc runs
! if(osc)call oscav(ats,au,afn,icount)
! endif

! Calculate meridional overturning streamfunction opsi on C grid only

            do j=0,jmax
             do k=0,kmax
              opsi(j,k) = 0.
              opsia(j,k) = 0.
              opsip(j,k) = 0.
             enddo
            enddo

            do 35 j=1,jmax-1
             do 40 k=1,kmax-1
              ou(j,k) = 0.
               do 45 i=1,imax
                ou(j,k) = ou(j,k) + cv(j)*u(2,i,j,k)*dphi
   45          continue
               opsi(j,k) = opsi(j,k-1) - dz(k)*ou(j,k)
   40        continue
   35       continue

! Pacific and Atlantic overturning streamfunctions

            ominp = 0.
            omaxp = 0.
            do j=jsf+1,jmax-1
             do k=1,kmax-1
              ou(j,k) = 0.
		      i=ips(j)  !####
		      do while (i.ne.(ipf(j)+1))
				 ou(j,k) = ou(j,k) + cv(j)*u(2,i,j,k)*dphi
				 i = mod(i,imax)+1    !####
              enddo
              opsip(j,k) = opsip(j,k-1) - dz(k)*ou(j,k)
              if(opsip(j,k).lt.ominp)ominp = opsip(j,k)
              if(opsip(j,k).gt.omaxp)omaxp = opsip(j,k)
             enddo
            enddo

            omina = 0.
            omaxa = 0.
            do j=jsf+1,jmax-1
! if continuing run, then just N.Atlantic
! do j=10,jmax-1
! do k=1,kmax-1
! excluding surface wind-driven cells (below 500m)
             do k=1,kmax-3
              ou(j,k) = 0.
			  i=ias(j) !####
		       do while (i.ne.(iaf(j)+1))
				 ou(j,k) = ou(j,k) + cv(j)*u(2,i,j,k)*dphi
				 i = mod(i,imax)+1 !####
               enddo
               opsia(j,k) = opsia(j,k-1) - dz(k)*ou(j,k)
               if(opsia(j,k).lt.omina)omina = opsia(j,k)
               if(opsia(j,k).gt.omaxa)omaxa = opsia(j,k)
             enddo
            enddo

            write(40,'(6e15.5)')t,extra0,ominp,omaxp,omina,omaxa

         endif !itstp
         
         
	!######## every month output
            open(42,file=trim(path_results)//lout//'.'//'MonthGlob')

         if(sdedy==15 .or.sdedy==197.or.sdedy==106 .or.sdedy==288)then
!         if(MNTHDY.eq.15)then
            call diag4(sum1,sum2,sum3,sum4,iice)
            write(42,119)istep,sum1,sum2,sum3,sum4,iice
         endif
      !###########
!############## for seasonal run, mean year values
	   tq_avr=tq_avr+tq
	   ice_avr=ice_avr+varice1
	   i_avr=i_avr+1

         if(mod(istep,itstp).eq.0)then
	           sdeyr0=sdeyr
         endif !itstp

	   if (sdeyr==sdeyr0.and.sdedy==365) then
	    tq_avr=tq_avr/i_avr
	    ice_avr=ice_avr/i_avr
! EMBM
            open(41,file=trim(path_results)//lout//'.'//'TimeGlob')
!########            open(42,file=trim(path_results)'//lout//'.'//'q')
!#######            open(43,file=trim(path_results)//lout//'.'//'hiceT')
!	 diag3: lout.airt  averaged values of air temperature
!      diag3: lout.q	   averaged values of specific humidity
            call diag3(sum1,sum2,sum3,sum4,iice)
            write(41,119)istep,sum1,sum2,sum3,sum4,iice
  119        format(i7,4e14.6,i4)
!            write(41,110)t*tsc/86400.,sum1
!            write(42,110)t*tsc/86400.,sum2
	   endif

	   if (sdedy==365) then
	    i_avr=0
          tq_avr=0.
	    ice_avr=0.
	   endif
!####################
!         if((istep.gt.24820).and.sdedy==15)then
! write out atmos air temp Jan

!      open(20,file=trim(path_results)//lout//'.tJan')
!      do  j=1,jmax
!         do  i=1,imax
!                     write(20,199 ) i,j,tq(1,i,j)
!	   enddo
!	enddo
  199    format (1x, i3,',',i3,',',e12.4)
!      close(20)
!	   endif

! Seasonal run: Jan and Jul last year output ##########

         if (istep.gt.(nsteps-365).and.(sdedy==15 .or.sdedy==197)) then
	    if (sdedy==15)  then
		   id_mnth='Jn'
	    endif
	    if (sdedy==197) then
		   id_mnth='Jl'
	    endif

          print*, ' ### Last year### '
          print*,'step ',istep,' t(days) ',t*tsc/86400.,'  dt ',dt(kmax)
          print*,'sdedy ',sdedy,' SDEYR ',SDEYR

! diag.f end-of-run diagnostics for gldstn
		call diagend(lout,id_mnth)

! put all the following junk into diagend when get round to it.??

! write out convective frequency array. Divide by 2*nsteps if call co twice

          open(3,file=trim(path_results)//lout//'.cost'//id_mnth)
! cost(maxi,maxj) - array of convection depth (in gridpoints) averaged over a run
          if(nsteps.gt.0)write(3,'(e15.8)')((cost(i,j)/nsteps,i=1,imax),j=1,jmax)
! if(nsteps.gt.0)write(3,'(e15.8)')((0.5*cost(i,j)/nsteps
          close(3)

! write out barotropic streamfunction

          open(3,file=trim(path_results)//lout//'.psi'//id_mnth)
          do 60 j=0,jmax
          do 60 i=0,imax
            write(3,*)psi(i,j)
  60      continue
          close(3)

! Calculate meridional overturning streamfunction opsi on C grid only

          do j=0,jmax
           do k=0,kmax
            opsi(j,k) = 0.
            opsia(j,k) = 0.
            opsip(j,k) = 0.
           enddo
          enddo

          do 70 j=1,jmax-1
           do 80 k=1,kmax-1
            ou(j,k) = 0
            do 90 i=1,imax
               ou(j,k) = ou(j,k) + cv(j)*u(2,i,j,k)*dphi
   90       continue
            opsi(j,k) = opsi(j,k-1) - dz(k)*ou(j,k)
   80      continue
   70     continue

          open(10,file=trim(path_results)//lout//'.opsi'//id_mnth)
          write(10,100)((opsi(j,k),j=0,jmax),k=0,kmax)
          close(10)

! Pacific and Atlantic overturning streamfunctions

          ominp = 0.
          omaxp = 0.
          do j=jsf+1,jmax-1
           do k=1,kmax-1
            ou(j,k) = 0.
            i=ips(j)  !####
		    do while (i.ne.(ipf(j)+1))
			 ou(j,k) = ou(j,k) + cv(j)*u(2,i,j,k)*dphi
			 i = mod(i,imax)+1  !####
            enddo
            opsip(j,k) = opsip(j,k-1) - dz(k)*ou(j,k)
            if(opsip(j,k).lt.ominp)ominp = opsip(j,k)
            if(opsip(j,k).gt.omaxp)omaxp = opsip(j,k)
           enddo
          enddo

          open(10,file=trim(path_results)//lout//'.opsip'//id_mnth)
          write(10,100)((opsip(j,k),j=0,jmax),k=0,kmax)
          close(10)
          open(101,file=trim(path_results)//lout//'.opsipS'//id_mnth)
          write(101,'(1x,2i4,e12.4)')((j,k,opsip(j,k)*1592.5,j=0,jmax),k=0,kmax)
          close(101)

          omina = 0.
          omaxa = 0.
          do j=jsf+1,jmax-1
           do k=1,kmax-1
            ou(j,k) = 0.
            i=ias(j)   !####
		    do while (i.ne.(iaf(j)+1))
			 ou(j,k) = ou(j,k) + cv(j)*u(2,i,j,k)*dphi
			 i = mod(i,imax)+1   !####
            enddo
            opsia(j,k) = opsia(j,k-1) - dz(k)*ou(j,k)
            if(opsia(j,k).lt.omina)omina = opsia(j,k)
            if(opsia(j,k).gt.omaxa)omaxa = opsia(j,k)
           enddo
          enddo

          open(10,file=trim(path_results)//lout//'.opsia'//id_mnth)
          write(10,100)((opsia(j,k),j=0,jmax),k=0,kmax)
          close(10)
          open(102,file=trim(path_results)//lout//'.opsiaS'//id_mnth)
          write(102,'(1x,2i4,e12.4)')((j,k,opsia(j,k)*1592.5,j=0,jmax),k=0,kmax)
          close(102)

! zonal overturning streamfunction

          do i=0,imax
           do k=0,kmax
            zpsi(i,k) = 0.
           enddo
          enddo

          do i=1,imax-1
           do k=1,kmax-1
            zu(i,k) = 0
            do j=1,jmax
               zu(i,k) = zu(i,k) + u(1,i,j,k)/c(j)*ds
            enddo
            zpsi(i,k) = zpsi(i,k-1) - dz(k)*zu(i,k)
           enddo
          enddo

          open(10,file=trim(path_results)//lout//'.zpsi'//id_mnth)
          write(10,100)((zpsi(i,k),i=0,imax),k=0,kmax)
          close(10)

  100     format(e14.7)
  110     format(11e14.6)

! write poleward heat flux in Atlantic and Pacific and total

          pi=4*atan(1.0)
          open(15,file=trim(path_results)//lout//'.fofy'//id_mnth)
          write(15,'(10(a11,3x))')' latitude  ',' tot_tot   ',' Pac_tot   ',' Atl_tot   ',&
          ' tot_adv   ',' Pac_adv   ',' Atl_adv   ',' tot_dif   ',' Pac_dif   ',' Atl_dif   '
          phfmax = 0
          do j=1,jmax-1
           do l=1,3
            hft(l) = 0
            hfp(l) = 0
            hfa(l) = 0
           enddo
           do i=1,imax
            if(k1(i,j).le.kmax.and.k1(i,j+1).le.kmax)then
               tv2 = 0
               tv3 = 0
               do k=k1(i,j),kmax
                  tv2 = tv2 + 0.5*cv(j)*u(2,i,j,k)*(ts(1,i,j+1,k) +ts(1,i,j,k))*dz(k)*dphi
                  tv3 = tv3 - cv(j)*cv(j)*(ts(1,i,j+1,k) -ts(1,i,j,k))/ds*diff(1)*dz(k)*dphi
               enddo
               hft(1) = hft(1) + tv2 + tv3
               hft(2) = hft(2) + tv2
               hft(3) = hft(3) + tv3
               if(i.ge.ips(j).and.i.le.ipf(j))then
                  hfp(1) = hfp(1) + tv2 + tv3
                  hfp(2) = hfp(2) + tv2
                  hfp(3) = hfp(3) + tv3
               elseif(i.ge.ias(j).and.i.le.iaf(j))then
                  hfa(1) = hfa(1) + tv2 + tv3
                  hfa(2) = hfa(2) + tv2
                  hfa(3) = hfa(3) + tv3
               endif
            endif
           enddo
           write(15,110)180.0/pi*asin(s(j)),(hft(l),hfp(l),hfa(l),l=1,3)
           if(abs(hft(3)).gt.phfmax)phfmax = abs(hft(3))
          enddo

          write(6,*)'max poleward heat flux ',phfmax

          close(15)

! write out potential vorticity or density

! open(11,file=trim(path_results)//lout//'.pv')
          open(11,file=trim(path_results)//lout//'.rho'//id_mnth)

! print*,'input k for Sz plot'
! read (5,*)k
          do j=1,jmax
           do i=1,imax
! do k=1,kmax-1
            do k=1   ,kmax
               if(k.ge.k1(i,j))then
! tmp = s(j)*(rho(i,j,k+1)-rho(i,j,k))/dza(k)
! tmp = (ts(2,i,j,k+1)-ts(2,i,j,k))/dza(k)
! write(11,*)tmp
                  write(11,*)rho(i,j,k)
               else
                  write(11,*)0.0
               endif
            enddo
           enddo
          enddo

          close(11)

! call routine to write terms in psi equation, gb is subsequently wrong

! call goose(lout)

! ropsi calculates overturning on steady density surfaces

! call ropsi

! open(11,file=trim(path_results)//'tmp.w')
! do k=kmax,1,-1
! tmp = 0
! do i=1,imax
! do j=1,jmax/2
! tmp =  u(3,5,8,k)
! enddo
! enddo
! write(11,*)tmp
! enddo
! close(11)

! write average osc data

! if(osc.and.icount.gt.0)then
! open(10,file=trim(path_results)//lout//'.av')
! write(10,100)((((ats(l,i,j,k)/icount
! 1             ,l=1,lmax),k=1,kmax),i=1,imax),j=1,jmax)
! close(10)
! open(10,file=trim(path_results)//lout//'.avu')
! write(10,100)((((au(l,i,j,k)/icount,l=1,3)
! 1             ,k=1,kmax),i=1,imax),j=1,jmax)
! close(10)
! open(10,file=trim(path_results)//lout//'.avf')
! write(10,100)(((afn(1,i,j,k)/icount,k=1,kmax),i=1,imax)
! 1                ,j=1,jmax)
! close(10)
! endif

! for repeated runs

! write(20,'(8e15.5)')0.01*iterun,ominp,omaxp,omina,omaxa
! 1                   ,phfmax,avn,rms
          write(6,'(a)')'ominp,omaxp,omina,omaxa,avn,rms'
          write(6,'(6e15.5)')ominp,omaxp,omina,omaxa,avn,rms
! close(5)
! enddo
! close(20)
	   endif ! Seasonal run: Jan and Jul last year output #######
      if (istep==nsteps) goto 567 
      enddo   ! time loop end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        return   ! to AGCM

 567     close(4)
      close(5)
      close(11)
      close(14)
	  close(39)
      close(40)
      close(41)
      close(42)
	  close(52) !####
!#########      close(43)
! close(8)
! close(9)
      print *,' sdedy_oc = ',sdedy_oc, istep
      !pause ' pause GLDSTN STEP END'
      stop 'normal end'

      end !program end

      character*3 function conv(i)
      character*1 a,b,c
      if(i.lt.10)then
        a=char(i+48)
        conv=a//'  '
      else if(i.lt.100)then
        i1=i/10
        i2=i-i1*10
        a=char(i1+48)
        b=char(i2+48)
        conv=a//b//' '
      else
        i1=i/100
        itemp=i-100*i1
        i2=itemp/10
        i3=itemp-10*i2
        a=char(i1+48)
        b=char(i2+48)
        c=char(i3+48)
        conv=a//b//c
      endif
      end
