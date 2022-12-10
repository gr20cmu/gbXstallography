C*****************************************
       program graph_discrete_dist
C*****************************************
       implicit none
	   include 'common.fi'

       write(6,1) '================================================='
	   write(6,1) 'PROGRAM TO GRAPH 2D OR 5D GRAIN BOUNDARY PLANE '
       write(6,1) 'DISTRIBUTIONS FROM DISCRETIZED DATA '
	   write(6,1) ' '
	   write(6,1) 'rohrer@cmu.edu'
	   write(6,1) 'version 12/10/2022'
       write(6,1) '================================================='
	   
	   !version = 'version 04/09/2022'
       !For single 2D plots, use DS_1
       !for 5D plots, provide options for different symmetries.
       version = 'version 12/10/2022' !fixed dimensions of the gbd array to
                                      !compute for resolutions .gt. 9

       !some constants that might be useful
       pi = 4.0*atan(1.0)
       eps = 1.e-6
       degrad = 180.0/pi
       out = 0.0
       AxAng = 0.0
       AA = 0.0

       !some formating statements used to format typed output
 1     format(A)
 2     format(A,A)
 4     format(A,A,A,I10,A)
 5     format(A,I8,A)
 6     format(f8.5,3x,f8.5,3x,f8.3,2(i4))
 7     format(A,f5.2,A)
 8     format(A,f7.2,A)
 9     format(A,A,A,f4.2,1x,f4.2,1x,f4.2,A)
 10    format(A,A,A,f4.2,1x,f4.2,1x,f4.2,A)
 11    format(A,f5.2,2X,f5.2,2X,f5.2,2X,f6.2,2X)
 12    format(A,I1,1x,A,A,f6.2,1x,f7.2,1x,f6.2,A)


       !The program requires a file, 'input.txt' that specifies all the parameters
       !Read the lower part of the input file for the definition of these parameters
       open(21, file='input.txt',status='old')
       read(21,*)gbpd_FileName           !name of the file with the 2D GB distribution
       read(21,*)gbcd_FileName           !name of the file with the 5D GB distribution
       read(21,*)msym                    !crystal symmetry
       read(21,*)CD, CD2                 !bin size, 90°/CD
       read(21,*)out(1),out(2),out(3),out(4),out(5)!specify output type: 0/1 = off/on
       read(21,*)normal_dir
       read(21,1)label
       read(21,*)Num5DPlots
       read(21,*)NumAAPlots
       if (Num5DPlots.gt.9) then
        write (6,1)'9 is the maximum number of 5D plots I can draw.  Bye'
        close(21)
        goto 9000
       endif
       if (NumAAPlots.gt.9) then
        write (6,1)'9 is the maximum number of axis-angle plots I can draw.  Bye'
        close(21)
        goto 9000
       endif
       do i1=1,Num5DPlots
        read(21,*)AxAng(i1,1),AxAng(i1,2),AxAng(i1,3),AxAng(i1,4)
       enddo
       do i1=1,NumAAPlots
        read(21,*)AA(i1)
       enddo

       close(21)                         !close the parameter file
 
	   if (msym.eq.1) then
        write(6,1)'you have specified the symmetry as tetragonal'
	   endif
	   if (msym.eq.2) then
        write(6,1)'you have specified the symmetry as hexagonal'
	   endif
	   if (msym.eq.3) then
        write(6,1)'you have specified the symmetry as cubic'
	   endif
	   if (msym.eq.4) then
        write(6,1)'you have specified the symmetry as trigonal'
	   endif
	   if (msym.eq.5) then
        write(6,1)'you have specified the symmetry as orthorhombic'
	   endif


       !------------------------------------------------------------------
       ! Begin the section of code that graphs the 2D grain boundary
       ! plane distribution.
       !------------------------------------------------------------------
       if (out(1).ne.1.AND.out(2).ne.1) goto 2000
	   write(6,2)'I am graphing data in the file labeled: ',trim(gbpd_FileName)
       write(6,1)'The [001] direction will be normal to the plane of the drawing.'
	   if (out(1).eq.1) then
	    write(6,1)'The stereogram will be plotted in the fundamental zone'
	   endif
	   if (out(2).eq.1) then
	    write(6,1)'The stereogram will be plotted in the complete hemisphere'
	   endif
       pd = 0.0
       call get_symop (msym, O, nsymm)
       !call get_symop (msym, O, nsymm) calls the subroutine that contains the symmetry operators.
       !The variable msym, read from the input file, specifies the symmetry (1 = tetragonal,
       !2 = hexagonal, 3 = cubic, 4 = trigonal, 5 = orthorhombic), O is a matrix containing the
       !operators, and nsymm is the number of symmetry operators.  The routine returns O and nsymm
       point = index(gbpd_FileName,'.') - 1  !number of characters in filename, before the extension
       open (30, file=trim(gbpd_FileName),status='unknown') !Open the data file
       nnline = 0  !counter
 100   continue
       read(30,*,end=101)
       nnline=nnline+1
       goto 100
 101   continue
       close(30)
       !Check each line for the string '#'. If found, we know this line is part
       !of the header.  If not, we have reached the data section.
       open(30, file=trim(gbpd_FileName),status='unknown')
       header = 0
       hash = 0
       do i1=1,nnline
        read(30, "(a)") inline
        hash = index(inline, '#')
        if (hash.ne.0) then
         header=header+1
         hash=0
        else
         goto 120
        endif
       enddo ! closes the i1=1,nnline loop
 120   close(30)
       !Here, inform the user of what we found
       write(6,4)'File ',trim(gbpd_FileName),' has ',nnline,' lines.'
       write(6,5)'There are ',header,' lines in the header'
       !read the distribution into matrix pd
       open (30, file=gbpd_FileName, status='old')
       do i1=1,header  ! read through the header
        read(30, "(a)") inline
       enddo
       do i1=1,CD2
        do i2=1,4*CD2
         read (30,*) pd(i1,i2)
        enddo
       enddo
       close (30)
       avg=0.0 ! determine the average value
       do i1=1,CD2
        do i2=1,4*CD2
         avg=avg+pd(i1,i2)
        enddo
       enddo
       avg=avg/float(4*CD2*CD2)

       if (out(1).eq.1) then
        open (40, file=gbpd_FileName(1:point)//'_2d_gmt1.gpf',status='unknown')
        write(40,"(4(1x,f4.1))")  0.0, 0.0, 0.0, 0.0
       endif

       if (out(2).eq.1) then
        open (41, file=gbpd_FileName(1:point)//'_2d_s_gmt.dat',status='unknown')
        write(41,"(4(1x,f4.1))")  0.0, 0.0, 0.0, 0.0
       endif

	   summrd=0.0
	   maxmrd=0.0
	   minmrd=1000.0
	   do i_plot=1,31
		do j_plot=1,121
		 theta_plot=3.0*float(i_plot-1)*(pi/180.0)
         phi_plot=3.0*float(j_plot-1)*(pi/180.0)
         !this call takes angle theta and phi and convert them to a vector, n_plot
         call AnglesToV (theta_plot, phi_plot, n)
         sum=0.0
         ct=0
         do i1=1,nsymm
          call symop (O, i1, so_1) !returns the first symmetry operator
          call MToV (so_1,n,nf) ! rotate n by the first, non-transposed operator
          call VToAngles (nf, sa) !returns azimuthal (0-90) and in-plane (0-360) angles
          c4=int(float(CD2)*cos(sa(1)))+1
          c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
          if (c4.eq.CD2+1) c4=CD2
          if (c4.eq.0) c4=1
          if (c5.eq.4*CD2+1) c5=4*CD2
          sum=sum+pd(c4,c5)
          ct=ct+1
         enddo
         !this projects the angles for the stereographic projection
	     xp=tan(theta_plot*0.5)*cos(phi_plot)
	     yp=tan(theta_plot*0.5)*sin(phi_plot)
	     amrd=sum/float(ct)
         if (out(1).eq.1) then
          write (40,*) phi_plot*degrad, (90.-(theta_plot*degrad)), amrd
         endif
         if (out(2).eq.1) then
          write (41,*) phi_plot*degrad, (90.-(theta_plot*degrad)), amrd
         endif
		 summrd=summrd+amrd
		 if (maxmrd.lt.amrd) maxmrd=amrd
		 if (minmrd.gt.amrd) minmrd=amrd
        enddo ! j_plot=1,121
	   enddo ! i_plot=1,31
       if (out(1).eq.1) then
        close(40)
       endif
       if (out(2).eq.1) then
        close(41)
       endif

	   write(6,7) 'The average value is ',summrd/3751,' MRD.'
	   write(6,8) 'The minimum value is ',minmrd,' MRD.'
	   write(6,8) 'The maximum value is ',maxmrd,' MRD.'
       interval=(maxmrd-minmrd)/10

       if (out(1).eq.1) then
        write(6,1) ' '
        if (msym.eq.3) then
         write(6,1)'to initiate the gmt script, enter ./Draw_stereograms 1 [filename]_2d_gmt IF rainbow [min] [max] [interval] stereo CUBIC'
         write(6,10)'For example: ./Draw_stereograms 1 ',gbpd_FileName(1:point),'_2d_gmt IF rainbow ',minmrd,maxmrd,interval,' stereo CUBIC'
        endif
        if (msym.eq.2) then
         write(6,1)'to initiate the gmt script, enter ./Draw_stereograms 1 [filename]_2d_gmt IF rainbow [min] [max] [interval] stereo HEX'
         write(6,10)'For example: ./Draw_stereograms 1 ',gbpd_FileName(1:point),'_2d_gmt IF rainbow ',minmrd,maxmrd,interval,' stereo HEX'
        endif
       endif

       if (out(2).eq.1) then
        write(6,1) ' '
        if (msym.eq.3) then
         write(6,1)'to initiate the gmt script, enter ./Draw_stereograms 1 [filename]_2d_gmt 2d rainbow [min] [max] [interval]'
         write(6,9)'For example: ./Draw_stereograms_1 1 ',gbpd_FileName(1:point),'_2d_s_gmt 2d rainbow ',minmrd,maxmrd,interval,' stereo CUBIC'
        endif
        if (msym.eq.2.OR.msym.eq.4) then
         write(6,1)'to initiate the gmt script, enter ./Draw_stereograms 1 [filename]_2d_gmt 2d rainbow [min] [max] [interval]'
         write(6,9)'For example: ./Draw_stereograms_1 1 ',gbpd_FileName(1:point),'_2d_s_gmt 2d rainbow ',minmrd,maxmrd,interval,' stereo HEX'
        endif
        if (msym.eq.1.OR.msym.eq.5) then
         write(6,1)'to initiate the gmt script, enter ./Draw_stereograms 1 [filename]_2d_gmt 2d rainbow [min] [max] [interval]'
         write(6,9)'For example: ./Draw_stereograms_1 1 ',gbpd_FileName(1:point),'_2d_s_gmt 2d rainbow ',minmrd,maxmrd,interval,' stereo ORT'
        endif
       endif

 2000  continue

       !------------------------------------------------------------------
       ! Begin the section of code that graphs the axis-angle
       ! distribution.
       !------------------------------------------------------------------
       if (out(3).ne.1) goto 3000
       gbd = 0.0
       aad = 0.0
	   write(6,2)'I am graphing data in the file labeled: ',trim(gbcd_FileName)
       call get_symop (msym, O, nsymm)
       !call get_symop (msym, O, nsymm) calls the subroutine that contains the symmetry operators.
       !The variable msym, read from the input file, specifies the symmetry (1 = tetragonal,
       !2 = hexagonal, 3 = cubic, 4 = trigonal, 5 = orthorhombic), O is a matrix containing the
       !operators, and nsymm is the number of symmetry operators.  The routine returns O and nsymm
       point = index(gbcd_FileName,'.') - 1  !number of characters in filename, before the extension
       open (32, file=trim(gbcd_FileName),status='unknown') !Open the data file
       nnline = 0  !counter
 2100  continue
       read(32,*,end=2101)
       nnline=nnline+1
       goto 2100
 2101  continue
       close(32)
       !Check each line for the string '#'. If found, we know this line is part
       !of the header.  If not, we have reached the data section.
       open(32, file=trim(gbcd_FileName),status='unknown')
       header = 0
       hash = 0
       do i1=1,nnline
        read(32, "(a)") inline
        hash = index(inline, '#')
        if (hash.ne.0) then
         header=header+1
         hash=0
        else
         goto 2120
        endif
       enddo ! closes the i1=1,nnline loop
 2120  close(32)
       !Here, inform the user of what we found
       write(6,4)'File ',trim(gbcd_FileName),' has ',nnline,' lines.'
       write(6,5)'There are ',header,' lines in the header'
       write (6,"(A,I1,A)") 'I am planning to make ',NumAAPlots,' plots'
       !read the distribution into matrix gbd
       if(msym.eq.3) then  ! for cubic
        open (32, file=trim(gbcd_FileName), status='old')
		do i1=1,header
		 read(32, "(a)") inline
		enddo
        do i1=1,CD
         do i2=1,CD
          do i3=1,CD
           do i4=1,CD2
            do i5=1,4*CD2
             read(32,*) gbd(i1,i2,i3,i4,i5)
            enddo
           enddo
          enddo
         enddo
        enddo
        close (32)
       else           ! for lower symmetry
        open (32, file=trim(gbcd_FileName), status='old')
		do i1=1,header
		 read(32, "(a)") inline
		enddo
        do i1=1,CD*4
         do i2=1,CD*2
          do i3=1,CD*4
           do i4=1,CD2
            do i5=1,4*CD2
             read(32,*) gbd(i1,i2,i3,i4,i5)
            enddo
           enddo
          enddo
         enddo
        enddo
        close (32)
       endif
       !Next we project the 5D distribution to 3D
       if(msym.eq.3) then  ! for cubic
        do i1=1,CD
         do i2=1,CD
          do i3=1,CD
           do i4=1,CD2
            do i5=1,4*CD2
             tem=gbd(i1,i2,i3,i4,i5)
             aad(i1,i2,i3)=aad(i1,i2,i3)+tem
            enddo
           enddo
          enddo
         enddo
        enddo
       else
        do i1=1,CD*4
         do i2=1,CD*2
          do i3=1,CD*4
           do i4=1,CD2
            do i5=1,4*CD2
             tem=gbd(i1,i2,i3,i4,i5)
             aad(i1,i2,i3)=aad(i1,i2,i3)+tem
            enddo
           enddo
          enddo
         enddo
        enddo
       endif
       !next, normalized the distribution
       if(msym.eq.3) then  ! for cubic
        avg = 0.0
        do i1=1,CD
         do i2=1,CD
          do i3=1,CD
           avg=avg+aad(i1,i2,i3)
          enddo
         enddo
        enddo
        avg=avg/float(CD*CD*CD)
        do i1=1,CD
         do i2=1,CD
          do i3=1,CD
           aad(i1,i2,i3)=aad(i1,i2,i3)/avg
          enddo
         enddo
        enddo
       else
        avg = 0.0
        do i1=1,CD*4
         do i2=1,CD*2
          do i3=1,CD*4
           avg=avg+aad(i1,i2,i3)
          enddo
         enddo
        enddo
        avg=avg/float(4*2*4*CD*CD*CD)
        do i1=1,CD*4
         do i2=1,CD*2
          do i3=1,CD*4
           aad(i1,i2,i3)=aad(i1,i2,i3)/avg
          enddo
         enddo
        enddo
       endif

       do plot=1,NumAAPlots
        first = .true.
        write(6,"(A,I1)") 'plot number= ',plot
        ang=AA(plot)
        write (6,"(A,f6.2)") 'now plotting: ', ang
        ang=ang*(pi/180.0)
        sum=0.0
        ct=0
        one=char(48+plot)
        if (msym.eq.2.OR.msym.eq.3) then
         open(43,file=gbcd_FileName(1:point)//'_'//trim(label)//'_AA_gmt_'//one//'.gpf',status='unknown')
        else
         open(43,file=gbcd_FileName(1:point)//'_'//trim(label)//'_AA_gmt_'//one//'.dat',status='unknown')
        endif
        write(43,"(1x,f5.1)")ang*degrad
        summrd=0.0
        maxmrd = 0.0
        minmrd = 1.e+36  ! large value
        do i_plot=1,31
         do j_plot=1,121
          theta_plot=3.0*float(i_plot-1)*(pi/180.0)
          phi_plot=3.0*float(j_plot-1)*(pi/180.0)
          call AnglesToV (theta_plot, phi_plot, ax)
          sum=0.0
          ct=0
          len=sqrt(ax(1)*ax(1)+ax(2)*ax(2)+ax(3)*ax(3))
          ax(1)=ax(1)/len
          ax(2)=ax(2)/len
          ax(3)=ax(3)/len
          !this converts the axis and angle values to a g matrix
          call AAToG (ax, ang, dg)
          !in the next loops, dg is hit with all the symmetry operators
          !then the bin corresponding to each misorientation is queried.
          !Values read from all of indistinguishable bins are then
          !averaged to yeild a value for this specific misorientation.
          do ii_plot=1,nsymm
           do jj_plot=1,nsymm
            call symop (O, ii_plot, so_1)
            call symop (O, jj_plot, so_2)
            call trans (so_2,so_2t)!transpose of the second symmetry operator
            call MToM (dg, so_2t, gg_1)!determine misorientation and  Euler angles
            call MToM (so_1, gg_1, dg_1)
            call GToE (dg_1, mort)
            !find the bin that corresponds to three parameters, mort(3)
            if (msym.eq.3) then
             if (mort(1).lt.(pi*0.5+eps).and.mort(2).lt.(pi*0.5+eps).and.mort(3).lt.(pi*0.5+eps)) then
     		  c1=int(float(CD)*2.0*mort(1)/pi)+1
			  c2=int(float(CD)*cos(mort(2)))+1
			  c3=int(float(CD)*2.0*mort(3)/pi)+1
              if (c1.eq.CD+1) c1=CD
              if (c2.eq.CD+1) c2=CD
              if (c3.eq.CD+1) c3=CD
              sum=sum+aad(c1,c2,c3)!add the value in this bin
              ct=ct+1!increment the counter
             endif
            else
             c1=int(float(2*CD)*mort(1)/pi)+1
             if (mort(2).lt.pi/2.0) then
              c2=int(float(CD)*cos(mort(2)))+1
             else
              c2=9+int(float(CD)*abs(cos(mort(2))))+1
             endif
             c3=int(float(2*CD)*mort(3)/pi)+1
             if (c1.eq.4*CD+1) c1=4*CD
             if (c2.eq.2*CD+1) c2=2*CD
             if (c3.eq.4*CD+1) c3=4*CD
             sum=sum+aad(c1,c2,c3)!add the value in this bin
             ct=ct+1!increment the counter
            endif
            !apply crystal exchange symmetry
            call trans (dg, dgt)
            call MToM (dgt, so_2, gg_2)
            call MToM (so_1, gg_2, dg_2)
            call GToE (dg_2, mort)
            if (msym.eq.3) then
             if (mort(1).lt.(pi*0.5+eps).and.mort(2).lt.(pi*0.5+eps).and.mort(3).lt.(pi*0.5+eps)) then
     		  c1=int(float(CD)*2.0*mort(1)/pi)+1
			  c2=int(float(CD)*cos(mort(2)))+1
			  c3=int(float(CD)*2.0*mort(3)/pi)+1
              if (c1.eq.CD+1) c1=CD
              if (c2.eq.CD+1) c2=CD
              if (c3.eq.CD+1) c3=CD
              sum=sum+aad(c1,c2,c3)!add the value in this bin
              ct=ct+1!increment the counter
             endif
            else
             c1=int(float(2*CD)*mort(1)/pi)+1
             if (mort(2).lt.pi/2.0) then
              c2=int(float(CD)*cos(mort(2)))+1
             else
              c2=9+int(float(CD)*abs(cos(mort(2))))+1
             endif
             c3=int(float(2*CD)*mort(3)/pi)+1
             if (c1.eq.4*CD+1) c1=4*CD
             if (c2.eq.2*CD+1) c2=2*CD
             if (c3.eq.4*CD+1) c3=4*CD
             sum=sum+aad(c1,c2,c3)
             ct=ct+1
            endif
           enddo   !jj_plot loop
          enddo   !ii_plot loop
          xp=tan(theta_plot*0.5)*cos(phi_plot)
          yp=tan(theta_plot*0.5)*sin(phi_plot)
          amrd=sum/float(ct)
		  write (43,*) phi_plot*degrad, (90.-(theta_plot*degrad)), amrd
		  summrd = summrd + amrd
		  if (maxmrd.lt.amrd) maxmrd = amrd
		  if (minmrd.gt.amrd) minmrd = amrd
          first = .false.
         enddo  ! j_plot loop
        enddo  ! i_plot loop

        write(6,7) 'For this misorientation, the average value is ',summrd/3600,' MRD.'
        write(6,8) 'For this misorientation, the maximum value is ',maxmrd,' MRD.'
        write(6,8) 'For this misorientation, the minimum value is ',minmrd,' MRD.'
        close(43)

       enddo !   ends the plot=1,NumAAPlots loop

       interval=(maxmrd-minmrd)/10
       fname=gbcd_FileName(1:point)//'_'//trim(label)

        write(6,1) ' '

        if (msym.eq.3) then
         write(6,1)'to initiate the gmt script, enter ./Draw_stereograms 1 [filename]_AA_gmt IF rainbow [min] [max] [interval] stereo CUBIC'
         write(6,12)'For example: ./Draw_stereograms ',NumAAPlots,trim(fname),'_AA_gmt_ IF rainbow ',minmrd,maxmrd,interval,' stereo CUBIC'
        endif

        if (msym.eq.2) then
         write(6,1)'to initiate the gmt script, enter ./Draw_stereograms 1 [filename]_AA_gmt IF rainbow [min] [max] [interval] stereo HEX'
         write(6,12)'For example: ./Draw_stereograms ',NumAAPlots,trim(fname),'_AA_gmt_ IF rainbow ',minmrd,maxmrd,interval,' stereo HEX'
        endif


       if (msym.ne.2.AND.msym.ne.3) then
        write(6,1)'to initiate the gmt script, enter ./Draw_stereograms [Number of plots] [filename]_gmt_ 5d rainbow [min] [max] [interval] stereo ORT'
        write(6,12)'For example: ./Draw_stereograms ',NumAAPlots,trim(fname),'_AA_gmt_ 5d rainbow ',minmrd,maxmrd,interval,' stereo ORT'
       endif

 3000  continue

       !------------------------------------------------------------------
       ! Begin the section of code that graphs the 5D grain boundary
       ! plane distribution at a fixed misorientation.
       !------------------------------------------------------------------
       if (out(4).ne.1.AND.out(5).ne.1) goto 4000
	   write(6,2)'I am graphing data in the file labeled: ',trim(gbcd_FileName)
	   if (out(4).eq.1) then
	    write(6,1)'a single stereographic projection will be plotted'
	   endif
	   if (out(5).eq.1) then
	    write(6,1)'multiple stereographic projections will be plotted'
	   endif
       !set up the reference direction
       !for a [001] projection, there is no rotation of the data
	   if (normal_dir.eq.1) then
        ax_ref(1)=0.0
        ax_ref(2)=0.0
        ax_ref(3)=1.0
        ang_ref= 0.0
        write(6,1)'The [001] direction will be normal to the plane of the drawing.'
	   endif
       !for a [110] projection, rotate 90 deg around [1-10]
	   if (normal_dir.eq.2) then
        ax_ref(1)=0.7071
        ax_ref(2)=-0.7071
        ax_ref(3)=0.0
        ang_ref= pi*(90.0/180.0)
        write(6,1)'The [110] direction will be normal to the plane of the drawing.'
	   endif
       !for a [111] projection, rotate 54.7 deg around [1-10]
	   if (normal_dir.eq.3) then
        ax_ref(1)=0.7071
        ax_ref(2)=-0.7071
        ax_ref(3)=0.0
        ang_ref= pi*(54.7/180.0)
        write(6,1)'The [111] direction will be normal to the plane of the drawing.'
	   endif
      !for a [100] projection, rotate 90 deg around [0-10]
	   if (normal_dir.eq.4) then
        ax_ref(1)=0.0
        ax_ref(2)=-1.0
        ax_ref(3)=0.0
        ang_ref= pi*(90.0/180.0)
        write(6,1)'The [100] direction will be normal to the plane of the drawing.'
	   endif
	   write(6,1) ' '
       !this is the matrix that will rotate the vectors on the plot.
	   call AAToG (ax_ref, ang_ref, g_ref)
	   call trans (g_ref, g_ref_t)
       gbd = 0.0
       call get_symop (msym, O, nsymm)
       !call get_symop (msym, O, nsymm) calls the subroutine that contains the symmetry operators.
       !The variable msym, read from the input file, specifies the symmetry (1 = tetragonal,
       !2 = hexagonal, 3 = cubic, 4 = trigonal, 5 = orthorhombic), O is a matrix containing the
       !operators, and nsymm is the number of symmetry operators.  The routine returns O and nsymm
       point = index(gbcd_FileName,'.') - 1  !number of characters in filename, before the extension
       open (31, file=trim(gbcd_FileName),status='unknown') !Open the data file
       nnline = 0  !counter
 3100  continue
       read(31,*,end=3101)
       nnline=nnline+1
       goto 3100
 3101  continue
       close(31)
       !Check each line for the string '#'. If found, we know this line is part
       !of the header.  If not, we have reached the data section.
       open(31, file=trim(gbcd_FileName),status='unknown')
       header = 0
       hash = 0
       do i1=1,nnline
        read(31, "(a)") inline
        hash = index(inline, '#')
        if (hash.ne.0) then
         header=header+1
         hash=0
        else
         goto 3120
        endif
       enddo ! closes the i1=1,nnline loop
 3120  close(31)
       !Here, inform the user of what we found
       write(6,4)'File ',trim(gbcd_FileName),' has ',nnline,' lines.'
       write(6,5)'There are ',header,' lines in the header'

       !read the distribution into matrix gbd
       if(msym.eq.3) then  ! for cubic
        open (31, file=trim(gbcd_FileName), status='old')
		do i1=1,header
		 read(31, "(a)") inline
		enddo
        do i1=1,CD
         do i2=1,CD
          do i3=1,CD
           do i4=1,CD2
            do i5=1,4*CD2
             read(31,*) gbd(i1,i2,i3,i4,i5)
            enddo
           enddo
          enddo
         enddo
        enddo
        close (31)
       else           ! for lower symmetry
        open (31, file=trim(gbcd_FileName), status='old')
		do i1=1,header
		 read(31, "(a)") inline
		enddo
        do i1=1,CD*4
         do i2=1,CD*2
          do i3=1,CD*4
           do i4=1,CD2
            do i5=1,4*CD2
             read(31,*) gbd(i1,i2,i3,i4,i5)
            enddo
           enddo
          enddo
         enddo
        enddo
        close (31)
       endif

       if(out(5).eq.1) then
        write(6,"(A,i1,A)")'I am planning to make ',Num5DPlots,' plots'
       endif

       do plot=1,Num5DPlots  ! loop over multiple plots
        first = .true.
        write(6,"(A,I1)") 'plot number= ',plot
        write (6,11) 'now plotting: ',AxAng(plot,1),AxAng(plot,2),AxAng(plot,3),AxAng(plot,4)
        sum=0.0
        ct=0
        one=char(48+plot)
		!open the file for the results
        open (42, file=gbcd_FileName(1:point)//'_'//trim(label)//'_gmt_'//one//'.dat',status='unknown')
        write(42,"(4(1x,f5.1))")  AxAng(plot,1),AxAng(plot,2),AxAng(plot,3),AxAng(plot,4)
        !normalize the axis angle description of the boundary
        len=sqrt(AxAng(plot,1)*AxAng(plot,1)+AxAng(plot,2)*AxAng(plot,2)+AxAng(plot,3)*AxAng(plot,3))
        ax(1)=AxAng(plot,1)/len
        ax(2)=AxAng(plot,2)/len
        ax(3)=AxAng(plot,3)/len
        ang=AxAng(plot,4)*(pi/180.0)
        !convert axis-angle values to a g matrix
        call AAToG (ax, ang, dg)
        summrd=0.0
        maxmrd = 0.0
        minmrd = 1.e+36  ! large value
        !i_plot is the azimuthal angle.  It takes 30 discrete values between
        !0 and 90°, assigned to theta_plot.  j_plot is the in plane angle and
        !it takes 120 values in range between 0 and 360, assigned to phi_plot.
        do i_plot=1,31
         do j_plot=1,121
          theta_plot=3.0*float(i_plot-1)*(pi/180.0)
          phi_plot=3.0*float(j_plot-1)*(pi/180.0)
          call AnglesToV (theta_plot, phi_plot, n_plot_1)!use theta and phi to make a vector
          call MToV (g_ref_t, n_plot_1, n_1) ! rotate to reference direction
          call trans (dg, dgt)  !transpose of dg is dgt
          call MToV (dgt, n_1, n_2) !rotate n_1 by dgt to get n_2
          ct=0
          sum=0.0
          !in this next set of loops, n1 and n2 are hit with all the symmetry
          !operators, then the bin corresponding to each vector is queried.
          !Values read from all of indistinguishable bins are then averaged
          !to yeild a value for this specific orientation.
          do ii_plot=1,nsymm
           do jj_plot=1,nsymm
            call symop (O, ii_plot, so_1)!get the first symmetry operator
            call MToV (so_1, n_1, nf) !rotate n_1 by the first, non-transposed operator
            call symop (O, jj_plot, so_2) !get the second symmetry operator
            call trans (so_2,so_2t) !transpose the second symmetry operator
            call MToM (dg, so_2t, gg_1) ! determine misorientation
            call MToM (so_1, gg_1, dg_1)
            call GToE (dg_1, mort) !and Euler angles
            if (msym.eq.3) then
            if (mort(1).lt.(pi*0.5+eps).and.mort(2).lt.(pi*0.5+eps).and.mort(3).lt.(pi*0.5+eps)) then
              call VToAngles (nf, sa) !get spherical angles for the vector
     		  c1=int(float(CD)*2.0*mort(1)/pi)+1
			  c2=int(float(CD)*cos(mort(2)))+1
			  c3=int(float(CD)*2.0*mort(3)/pi)+1
              c4=int(float(CD2)*cos(sa(1)))+1
              c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
              if (c1.eq.CD+1) c1=CD
              if (c2.eq.CD+1) c2=CD
              if (c3.eq.CD+1) c3=CD
              if (c4.eq.CD2+1) c4=CD2
              if (c4.eq.0) c4=1
              if (c5.eq.4*CD2+1) c5=4*CD2
              sum=sum+gbd(c1,c2,c3,c4,c5)
              ct=ct+1
             endif
            else
             call VToAngles (nf, sa) !get spherical angles for the vector
             c1=int(float(2*CD)*mort(1)/pi)+1
             if (mort(2).lt.pi/2.0) then
              c2=int(float(CD)*cos(mort(2)))+1
             else
              c2=9+int(float(CD)*abs(cos(mort(2))))+1
             endif
             c3=int(float(2*CD)*mort(3)/pi)+1
			 c4=int(float(CD2)*cos(sa(1)))+1
             c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
             if (c1.eq.4*CD+1) c1=4*CD
             if (c2.eq.2*CD+1) c2=2*CD
             if (c3.eq.4*CD+1) c3=4*CD
             if (c4.eq.CD2+1) c4=CD2
             if (c4.eq.0) c4=1
             if (c5.eq.4*CD2+1) c5=4*CD2
             sum=sum+gbd(c1,c2,c3,c4,c5)
             ct=ct+1
            endif
            !Next, we apply crystal exchange symmetry
            call MToV (so_1, n_2, nf)
            call trans (dg, dgt)
            call MToM (dgt, so_2, gg_2)
            call MToM (so_1, gg_2, dg_2)
            call GToE (dg_2, mort)
            if (msym.eq.3) then
            if (mort(1).lt.(pi*0.5+eps).and.mort(2).lt.(pi*0.5+eps).and.mort(3).lt.(pi*0.5+eps)) then
              call VToAngles (nf, sa) !get spherical angles for the vector
     		  c1=int(float(CD)*2.0*mort(1)/pi)+1
			  c2=int(float(CD)*cos(mort(2)))+1
			  c3=int(float(CD)*2.0*mort(3)/pi)+1
              c4=int(float(CD2)*cos(sa(1)))+1
              c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
              if (c1.eq.CD+1) c1=CD
              if (c2.eq.CD+1) c2=CD
              if (c3.eq.CD+1) c3=CD
              if (c4.eq.CD2+1) c4=CD2
              if (c4.eq.0) c4=1
              if (c5.eq.4*CD2+1) c5=4*CD2
              sum=sum+gbd(c1,c2,c3,c4,c5)
              ct=ct+1
             endif
            else
             call VToAngles (nf, sa) !get spherical angles for the vector
             c1=int(float(2*CD)*mort(1)/pi)+1
             if (mort(2).lt.pi/2.0) then
              c2=int(float(CD)*cos(mort(2)))+1
             else
              c2=9+int(float(CD)*abs(cos(mort(2))))+1
             endif
             c3=int(float(2*CD)*mort(3)/pi)+1
			 c4=int(float(CD2)*cos(sa(1)))+1
             c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
             if (c1.eq.4*CD+1) c1=4*CD
             if (c2.eq.2*CD+1) c2=2*CD
             if (c3.eq.4*CD+1) c3=4*CD
             if (c4.eq.CD2+1) c4=CD2
             if (c4.eq.0) c4=1
             if (c5.eq.4*CD2+1) c5=4*CD2
             sum=sum+gbd(c1,c2,c3,c4,c5)
             ct=ct+1
            endif
           enddo   !jj_plot loop
          enddo   !ii_plot loop
          !determine position on a stereographic projection
          xp=tan(theta_plot*0.5)*cos(phi_plot)
          yp=tan(theta_plot*0.5)*sin(phi_plot)
          amrd=sum/float(ct) !compute  average from the indistinuishable bins
          write (42,*) phi_plot*degrad, (90.-(theta_plot*degrad)), amrd
          summrd = summrd + amrd
          if (maxmrd.lt.amrd) maxmrd = amrd
          if (minmrd.gt.amrd) minmrd = amrd
          first = .false.
         enddo !    ends the j_plot=1,121 loop
        enddo !   ends the i_plot=1,31 loop

        write(6,7) 'For this misorientation, the average value is ',summrd/3600,' MRD.'
        write(6,8) 'For this misorientation, the maximum value is ',maxmrd,' MRD.'
        write(6,8) 'For this misorientation, the minimum value is ',minmrd,' MRD.'
        close(42)
	   
       enddo !   ends the plot=1,Num5DPlots loop

       interval=(maxmrd-minmrd)/10
       fname=gbcd_FileName(1:point)//'_'//trim(label)
       if (out(4).eq.1) then
        if (msym.eq.3) then
         write(6,1) ' '
         write(6,1)'to initiate the gmt script, enter ./Draw_stereograms_1 [Number of plots] [filename]_gmt_ 5d rainbow [min] [max] [interval]'
         write(6,12)'For example: ./Draw_stereograms_1 ',num5dplots,trim(fname),'_gmt_ 5d rainbow ',minmrd,maxmrd,interval,' stereo CUBIC'
        endif
        if (msym.eq.2.OR.msym.eq.4) then
         write(6,1) ' '
         write(6,1)'to initiate the gmt script, enter ./Draw_stereograms_1 [Number of plots] [filename]_gmt_ 5d rainbow [min] [max] [interval]'
         write(6,12)'For example: ./Draw_stereograms_1 ',num5dplots,trim(fname),'_gmt_ 5d rainbow ',minmrd,maxmrd,interval,' stereo HEX'
        endif
        if (msym.eq.1.OR.msym.eq.5) then
         write(6,1) ' '
         write(6,1)'to initiate the gmt script, enter ./Draw_stereograms_1 [Number of plots] [filename]_gmt_ 5d rainbow [min] [max] [interval]'
         write(6,12)'For example: ./Draw_stereograms_1 ',num5dplots,trim(fname),'_gmt_ 5d rainbow ',minmrd,maxmrd,interval,' stereo ORT'
        endif
       endif

       if (out(5).eq.1) then
        if (msym.eq.3) then
         write(6,1) ' '
         write(6,1)'to initiate the gmt script, enter ./Draw_stereograms [Number of plots] [filename]_gmt_ 5d rainbow [min] [max] [interval]'
         write(6,12)'For example: ./Draw_stereograms ',num5dplots,trim(fname),'_gmt_ 5d rainbow ',minmrd,maxmrd,interval,' stereo CUBIC'
        endif
        if (msym.eq.2.OR.msym.eq.4) then
         write(6,1) ' '
         write(6,1)'to initiate the gmt script, enter ./Draw_stereograms [Number of plots] [filename]_gmt_ 5d rainbow [min] [max] [interval]'
         write(6,12)'For example: ./Draw_stereograms ',num5dplots,trim(fname),'_gmt_ 5d rainbow ',minmrd,maxmrd,interval,' stereo HEX'
        endif
       if (msym.eq.1.OR.msym.eq.5) then
         write(6,1) ' '
         write(6,1)'to initiate the gmt script, enter ./Draw_stereograms [Number of plots] [filename]_gmt_ 5d rainbow [min] [max] [interval]'
         write(6,12)'For example: ./Draw_stereograms ',num5dplots,trim(fname),'_gmt_ 5d rainbow ',minmrd,maxmrd,interval,' stereo ORT'
        endif
       endif

 4000  continue


       write(6,1) 'Program complete'
 9000  continue

       end






	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
