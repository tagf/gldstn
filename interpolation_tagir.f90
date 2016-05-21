
!// subroutine interpolate.f for program gldstn
!// interpolates temperature of ocean on land linearly along altitudes
!// to include in mains:  call interpolate(ts(1,:,:,kmax), k1) before stop 'ts'


    subroutine interpolate(toc, k_norm)
    
    include 'var.f90'
     
    real toc(0:maxi+1,0:maxj+1), toc_shift(0:maxi,0:maxj+1) ! only 0:maxi dimension for i required in shifted array
	integer k_norm(0:maxi+1,0:maxj+1), k_shift(0:maxi, 0:maxj+1)
    integer i, j
    
    ! Tagir required
    integer i1, i2, normal_shift, strange_shift ! start and end of land on interpolation
    logical boolLand
    normal_shift = 18
    strange_shift = 0 ! if everything is correct, not required
    
    ! Experimantal for interpolation on AGCM grid
    real t_AGCM(0:maxi+1, 1:46)
    
    ! oc(1:56,:)=ts_oc_for_atm(17:72,:)             !Shift ts_oc_for_atm grid
    ! oc(57:72,:)=ts_oc_for_atm(1:16,:)
    ! k1(1:56,:)=k2(17:72,:)                        !Shift k1
    ! k1(57:72,:)=k2(1:16,:)
    
    ! code by Tagir
    ! diagnostics write (1447, *) maxi-normal_shift, maxi, maxi-normal_shift+1, normal_shift-2
    k_shift(0:maxi-normal_shift,:)=k_norm(normal_shift:maxi,:) ! Tagir
    k_shift(maxi-normal_shift+1: maxi-1,:)=k_norm(0:normal_shift-2, : ) ! Tagir
    k_shift(maxi, : )=k_norm(normal_shift, : ) ! only 0..72 array actually required
    
    ! strange_shift = 0
    toc_shift(0:maxi-normal_shift-strange_shift,:)=toc(normal_shift+strange_shift: maxi, : ) ! Tagir
    toc_shift(maxi-normal_shift-strange_shift+1: maxi-1, : )=toc(0:normal_shift+strange_shift-2, : ) ! Tagir
    toc_shift(maxi, : )=toc(normal_shift+strange_shift, : ) ! Tagir
    
    ! do i=1,72
    !    write (1287, *) i, k1(0,i), k1(72,i)
    ! enddo -- all correct
    
    ! without interpolation
    do j=72,1,-1 ! 
      do i=1,72 ! altitude 
        write (1446, *) i, j, toc_shift(i,j)
      enddo
    enddo
    
    ! ocean surf type 1 - 8  (ocean to atm interp begin, Tagir)
    do j=1,maxj ! 
      do i=0,maxi ! altitude 
        if (k_shift(i,j).gt.8) then ! not ocean
          if (boolLand) then ! we've stepped and know bounds i1 - i2 of Land
            ! T(i) = interpolate(i; i1, i2, T1, T2)
            toc_shift(i,j)=toc_shift(i1-1,j)+float(i-i1+1)*(toc_shift(i2,j)-toc_shift(i1-1,j))/float(i2-i1+1)
          else
            i1 = i ! start of land
            boolLand = .true. ! we've stepped onLand flag, for interpolating further
            do i2=i,72 ! we search for the end of land (ocean)
              if ((k_shift(i2,j).ge.1).and.(k_shift(i2,j).le.8)) then ! ocean, found end of land
                ! first cell on land interpolation
                toc_shift(i1,j)=toc_shift(i1-1,j) + 1*(toc_shift(i2,j)-toc_shift(i1-1,j))/float(i2-i1+1)
                ! write (1448, *), i2, j, toc_shift(i2, j) ! purely diagnostics of land-altitude-line end
                exit
              endif
            enddo
          endif ! boolLand
        else ! ocean
          boolLand = .false.
          ! temperature on ocean without changes
        endif ! ocean / not ocean
      enddo
    enddo
    ! ocean to atm interp procedure end (Tagir)

    ! with interpolation    
    do j=72,1,-1 ! 
      do i=1,72 ! altitude 
        write (1447, *) i, j, toc_shift(i,j)
      enddo
    enddo

  !!! Experimental block Begin: interpolation to AGCM grid
  do i=1,maxi
    do j=1,46 ! only along longtitude
       t_AGCM(i,j)=floor(sin(float(i * coeff)))
    enddo
  enddo
      
  !!! Experimantal block End
        
  return
  end subroutine	

