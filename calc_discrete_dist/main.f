C*****************************************
       program calculate_discrete_distributions
C*****************************************
       implicit none
	   include 'common.fi'

       write(6,1) '================================================='
	   write(6,1) 'PROGRAM FOR COMPUTING DISCRETIZED GRAIN BOUNDARY '
       write(6,1) 'DISTRIBUTIONS. IT CAN COMPUTE 1D DISORIENTATION'
	   write(6,1) 'DISTRIBUTIONS, 2D GRAIN BOUNDARY PLANE DISTRIBUTIONS'
	   write(6,1) '(GBPDs), AND 5D GRAIN BOUNDARY CHARACTER'
       write(6,*) 'DISTRIBUTIONS (GBCDs).  3D TRIANGLE DATA IS SIMPLY'
	   write(6,1) 'DISCRETIZED INTO BINS.  DATA CONSISTING OF GRAIN'
	   write(6,1) 'BOUNDARY TRACES OBSERVED ON A SECTION PLANE ARE'
	   write(6,1) 'INTERPRETED STEREOLOGICALLY ACCORDING TO '
	   write(6,1) 'SAYLOR ET AL., MET TRANS. 35A (2004) 1981. LINE'
       write(6,1) 'SEGMENT DATA CAN ALSO BE BINNED UNDER THE ASSUMPTION'
       write(6,1) 'THAT THE MICROSTRUCTURE IS COLUMNAR'
       write(6,1) ' '
	   write(6,1) 'rohrer@cmu.edu'
	   write(6,1) 'version 04/09/2022'
       write(6,1) '================================================='
	   
       !This program combines all prior programs for computing discrete
       !distributions.  Grain boundary line segment (trace) data is
       !interpreted stereologically or under the columnar assumption.
       !Grain boundary triangle data, from 3D studies, is interpreted
       !as 3D data.  This version uses the full domain of misorientations
       !for non-cubic structures and improves variable definitions.

	   version = 'version 04/09/2022' !corrected trigonal sym operators, and
                                      !problem with DisG routine in subs2

       !some constants that might be useful
       pi = 4.0*atan(1.0)
       eps = 1.e-6
       HCD = 90.0
       out=0

       !some formating statements used to format typed output
 1     format(A)
 2     format(A,A)
 3     format(A,I2,2X,A,I2,2X,A,I2,2X,A,F5.2)
 4     format(A,A,A,I10,A)
 5     format(A,I8,A)
 6     format(A,F5.1,A)
 7     format(f5.1,6x,F6.4,5x,F6.4)
 8     format(A,I2,A,I1)
 9     format(A,I1,A,I1,A,f7.3)
 10    format(A,I2,A,I2)
 11    format(A,7I2)

       !The program requires a file, 'input.txt' that specifies all the parameters
       !Read the lower part of the input file for the definition of these parameters
       open(21, file='input.txt',status='old')
       read(21,*)keyword                 !file name containing the GB line segment data
       read(21,*)NCol                    !number of columns (10 or 12) in the data file
       read(21,*)msym, rotation, radian  !crystal symmetry, reference frame, units of angles
       read(21,*)CD, CD2                 !bin size, 90°/CD
       read(21,*)StepSize                !pixel spacing
       read(21,*)out(1),out(2),out(3),out(4),out(5),out(6),out(7)
                                         !flags to cumpute(1) or not(0) distributions
       read(21,1)comment                 !comment added to output files
       close(21)                         !close the parameter file

       !Begin by testing if the parameters have allowable values
       if (NCol.ne.10.AND.NCol.ne.12) then
        write(6,1)'NCol must be 10 or 12. Read notes in input.txt; good bye.'
        goto 9000
       endif
       if (msym.ne.1.AND.msym.ne.2.AND.msym.ne.3.AND.msym.ne.4.AND.msym.ne.5) then
        write(6,1)'msym must be 1, 2, 3, 4, or 5. Read notes in input.txt; good bye.'
        goto 9000
       endif
       if (rotation.ne.0.AND.rotation.ne.1.AND.rotation.ne.2) then
        write(6,1)'rotation must be 0, 1, or 2. Read notes in input.txt; good bye.'
        goto 9000
       endif
       if (CD.gt.18.OR.CD2.gt.18) then
        write(6,1)'Neither CD nor CD2 can exceed 18; good bye.'
        goto 9000
       endif
       if (CD.le.3.OR.CD2.le.3) then
        write(6,1)'Neither CD nor CD2 can be beloe 4; good bye.'
        goto 9000
       endif
       if (StepSize.le.0.0) then
        write(6,1)'StepSize must be greater than zero. Read notes in input.txt; good bye.'
        goto 9000
       endif

       write(6,2)'Using data in file: ', trim(keyword)
       write(6,3)'msym =',msym,'rotation =',rotation,'radian =',radian,'step size =',StepSize

       point = index(keyword,'.') - 1    !Determines number of characters in the filename, before the extension
       open (30, file=trim(keyword),status='unknown') !Open the data file
       !The next section of code determines the length of the file (nnline)
       !and the length of the header (header)
       nnline = 0  !counter
 100   continue
       read(30,*,end=101)
       nnline=nnline+1
       goto 100
 101   continue
       close(30)
       !Check each line for the string '#'. If found, we know this line is part
       !of the header.  If not, we have reached the data section.
       open(30, file=trim(keyword),status='unknown')
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
       write(6,4)'File ',trim(keyword),' has ',nnline,' lines.'
       write(6,5)'There are ',header,' lines in the header'

       !Initialize matrices
       gbd=0.0

       !to keep it real, determine the average number of pixels per segment
       !this procedure is only used for line segment data
       if(out(2).eq.1.OR.out(3).eq.1.OR.out(6).eq.1.OR.out(7).eq.1) then
        open(30, file=trim(keyword),status='old')
        do i1=1,header
         read(30, "(a)") inline
        enddo
        l_tot=0.0
        t_out=0
 210    continue
        if (NCol.eq.10) then
         read(30,*,end=220)e1(1),e1(2),e1(3),e2(1),e2(2),e2(3),x1,y1,x2,y2
        else
         read(30,*,end=220)e1(1),e1(2),e1(3),e2(1),e2(2),e2(3),temp,temp,x1,y1,x2,y2
        endif
        nl(1)=x2-x1
        nl(2)=y2-y1
        length=sqrt((nl(1)*nl(1))+(nl(2)*nl(2)))
        l_step=length/StepSize
        l_tot=l_tot+l_step
        t_out=t_out+1
        goto 210
 220    continue
        close(30)
        l_ave=l_tot/float(t_out)
        write(6,6)'on average, there are ',l_ave,' steps per segment'
       endif

       call get_symop (msym, O, nsymm)
       !call get_symop (msym, O, nsymm) calls the subroutine that contains the symmetry operators.
       !The variable msym, read from the input file, specifies the symmetry (1 = tetragonal,
       !2 = hexagonal, 3 = cubic, 4 = trigonal, 5 = orthorhombic), O is a matrix containing the
       !operators, and nsymm is the number of symmetry operators.  The routine returns O and nsymm
       
       !------------------------------------------------------------------
       ! Begin the section of code that computes the disorientation
       ! distribution.
       !------------------------------------------------------------------
       if (out(1).ne.1) goto 2000
       Write(6,1)'calculating the disorientation angle distribution ...'
       ! This defines the maximum angle for a disorientation distribution
       if (msym.eq.1) limit=99  !max eq 98.42,  D4, 422, 4/m mm
	   if (msym.eq.2) limit=94  !max eq 93.84,  D6, 622, 6/m mm
	   if (msym.eq.3) limit=63  !max eq 62.7994, O, 432, m-3m
	   if (msym.eq.4) limit=105 !max eq 104.5,  D3,  32, -3m
	   if (msym.eq.5) limit=120 !max eq 120.0,  D2, 222, mmm
       bin = 1.0 ! this is the bin width in degrees.  Perhaps make it a variable?
	   bin_num = int(limit/bin) !computes the number of bins in the distribution
       !this just ensures the initial values of the distribution are zero
       do i1=1,2000
        d(i1,1)=0.0
        d(i1,2)=0.0
        d_norm(i1,1)=0.0
        d_norm(i1,2)=0.0
       enddo ! ends the i1=1,2000 loop
       l_tot=0.0
       !open the file for the results
       open(40, file=keyword(1:point)//'_disor_dist.txt',status='unknown')
       !open the file with the data
       open(30, file=trim(keyword), status='old')
       write(40,2)'# written by calc_discrete_dist, ',version
       write(40,2)'# based on data in file: ',trim(keyword)
       write(40,8)'# number of columns = ',NCol,', symmetry = ',msym
       write(40,9)'# ref. frame = ',rotation,', radians = ',radian,', StepSize = ',StepSize
       write(40,1)'# the discrization of the distribution is 1°'
       write(40,11)'# output flags: ',out(1),out(2),out(3), out(4),out(5),out(6),out(7)
       write(40,2)'# note: ',trim(comment)
       write(40,1)'# the following lines are the header of the source file'
       do i1=1,header
        read(30,1) inline ! read a line of the header
        write(40,1)inline ! write that line into the output file
       enddo! closes i1=1,header
       blank = 0  ! Initialize counter for lines with zero Euler angles
       !this is the main loop over all lines of data
       do i1=header+1,nnline
        if (NCol.eq.10) then
         if (out(4).eq.1.OR.out(5).eq.1) then
          read(30,*)e1(1),e1(2),e1(3),e2(1),e2(2),e2(3),nl(1),nl(2),nl(3),area
         else
          read(30,*)e1(1),e1(2),e1(3),e2(1),e2(2),e2(3),x1,y1,x2,y2
         endif
        else
         read(30,*)e1(1),e1(2),e1(3),e2(1),e2(2),e2(3),tem,tem,x1,y1,x2,y2
        endif
        if (out(4).eq.1.OR.out(5).eq.1) then
         length=area
        else
         nl(1)=x2-x1
         nl(2)=y2-y1
         length=sqrt((nl(1)*nl(1))+(nl(2)*nl(2)))
        endif
        !outputs from some programs produce zero Euler angles when
        !one or both orientiations are unknown.  These lines are
        !ignored because they will cause problems later.
        if (e1(1).eq.0.00.and.e1(2).eq.0.00.and.e1(3).eq.0.00) then
         blank=blank+1 ! increment counter
         goto 1000 ! Sends it to the next line
        else
         continue
        endif
        if (e2(1).eq.0.00.and.e2(2).eq.0.00.and.e2(3).eq.0.00) then
         blank=blank+1 ! increment counter
         goto 1000 ! Sends it to the next line
        else
         continue
        endif
        ! If the Euler angles are provided in units of degrees, these
        ! lines convert them to radians
        if (radian.eq.0) then
         call DToRad3(e1, e1r)
         e1=e1r
         call DToRad3(e2, e2r)
         e2=e2r
        else
         continue
        endif
        call DisG(O, nsymm, e1, e2, gg, axis, disor2)
        if (disor2.lt.0.01) then
         blank=blank+1
         goto 1000 ! occasional false boundaries removed
        endif
        ! The next section bins the observation in the matix d.
        ! The first column of d accumulates the lengths of the
        ! lines and the second the number accumulates the
        ! number of lines.
        disor2=disor2*(180.0/pi)     ! convert the value to degrees
        ind = floor(disor2/bin)+1    ! finds index of d (smallest integer > 0 for bin number)
        if (ind.eq.0) ind=1          ! make sure the index is non-zero
        if (ind.gt.bin_num) then     ! make sure it is not > the uper limit
         ind=bin_num
         write(6,*)'misorientation out of bounds'
        endif
        d(ind,1)=d(ind,1)+length     ! add the area to the first column of d
        d(ind,2)=d(ind,2)+1.0        ! increment the second column of d
        l_tot=l_tot+length           ! keep track ot the total area of all triangles
        ! when you reach this point, the orientation pair has been classified
 1000   continue
       enddo !this closes the i1=header+1,nnline loop
       close(30)

      !here we compute and write the distribution
      ct=nnline-header-blank         ! ct is the lines of useful data
      ! The next four lines normalize the data
      do i1=1,bin_num                ! loop over each bin
       d_norm(i1,1)=d(i1,1)/l_tot    ! divide area by the total area
       d_norm(i1,2)=d(i1,2)/float(ct)! divide number by the total number
      enddo
      !these statements write header information
      !write(6,1)' '
      write(40,1)'# '
      !write(6,1)'Disorientation distribution: rotation angle, area fraction, and number fraction'
      write(40,1)'# Disorientation distribution: rotation angle, area fraction, and number fraction'
      !write(6,1)' '
      write(40,1)'# '
      !write(6,1)'degrees    area       number'
      write(40,1)'  degrees  area       number'
      !these statements write each line of the distribution
      do i1=1,bin_num
       write(40,7)i1*bin, d_norm(i1,1),d_norm(i1,2)
       !write(6,7)i1*bin, d_norm(i1,1),d_norm(i1,2)
      enddo ! closes the i1=0,bin_num loop
      close(40) ! closes the output file
      write(6,1)'disorientation angle distribution complete'

      !------------------------------------------------------------------
      ! Begin the section of code that computes the two-dimensional
      ! grain boundary character distribution using stereology.
      !------------------------------------------------------------------

      !Note: this is based on calc_pd, version 12/19/2008.
      !On 03/28/22, this program exactly reproduced the GBPD calculated
      !by calc_pd version 12/19/2008 using all_IN100_segments.txt as input.
 2000 if (out(2).ne.1) goto 3000
      write(6,1)'calculating the 2D grain boundary plane distribution ...'
      pd=0.0
      open(41, file=keyword(1:point)//'_gbpd.txt',status='unknown')
      !open the file with the data
      open(30, file=trim(keyword), status='old')
       write(41,2)'# written by calc_discrete_dist, ',version
       write(41,2)'# based on data in file: ',trim(keyword)
       write(41,8)'# number of columns = ',NCol,', symmetry = ',msym
       write(41,9)'# ref. frame = ',rotation,', radians = ',radian,', StepSize = ',StepSize
       write(41,10)'# Discretization = ',CD,', ',CD2
       write(41,11)'# output flags: ',out(1),out(2),out(3), out(4),out(5),out(6),out(7)
       write(41,2)'# note: ',trim(comment)
       write(41,1)'# the following lines are the header of the source file'
      pd=0.0
      do i1=1,header
       read(30,1) inline ! read a line of the header
       write(41,1)inline ! write that line into the output file
      enddo! closes i1=1,header
      do i1=header+1,nnline
       if (NCol.eq.10) then
        read(30,*)e1(1),e1(2),e1(3),e2(1),e2(2),e2(3),x1,y1,x2,y2
       else
        read(30,*)e1(1),e1(2),e1(3),e2(1),e2(2),e2(3),tem,tem,x1,y1,x2,y2
       endif
       !outputs from some programs produce zero Euler angles when
       !one or both orientiations are unknown.  These lines are
       !ignored because they will cause problems later.
       if (e1(1).eq.0.00.and.e1(2).eq.0.00.and.e1(3).eq.0.00) then
        goto 2050 ! Sends it to the next line
       endif
       if (e2(1).eq.0.00.and.e2(2).eq.0.00.and.e2(3).eq.0.00) then
        goto 2050 ! Sends it to the next line
       endif
       ! If the Euler angles are provided in units of degrees, these
       ! lines convert them to radians
       if (radian.eq.0) then
        call DToRad3(e1, e1r)
        e1=e1r
        call DToRad3(e2, e2r)
        e2=e2r
       else
        continue
       endif
       !If the data are in the HKL reference frame, and the system is hexagonal,
       !then the data are rotated by -30° so that [2-1-10] is parallel to +x in
       !the output (this direction is to the right in the plane of the paper)
       if (rotation.eq.0.and.msym.eq.2) then
        if (e1(3).gt.(30*pi/180)) then
         e1(3) = e1(3)-(30*pi/180)
        else
         e1(3) = e1(3)+(330*pi/180)
        endif
        if (e2(3).gt.(30*pi/180)) then
         e2(3) = e2(3)-(30*pi/180)
        else
         e2(3) = e2(3)+(330*pi/180)
        endif
       endif
       !This section determines the in plane components of the grain boundary segment,
       !based on the endpoints in the data file.  The first calculation of nl(i)
       !transforms line segments into the same reference frame as the Euler angles,
       !based on and assumption of default settings in TSL or HKL software.  The second
       !calculation determines the normal to the grain boundary trace.
       if (rotation.eq.1) then       !TSL
        nl(1)=y1-y2
        nl(2)=x1-x2
       endif
       if (rotation.eq.0) then       !HKL
        nl(1)=x2-x1
        nl(2)=-(y2-y1)
       endif
       if (rotation.eq.2) then       !D3D/natural
        nl(1)=x2-x1
        nl(2)=y2-y1
       endif
       templ=nl(2)                   !normal to the GB trace
       nl(2)=nl(1)
       nl(1)=-templ
       nl(3)=0.0
       length=sqrt(nl(1)*nl(1)+nl(2)*nl(2))
       if (nl(1).ne.0.0) then
        p=atan(nl(2)/nl(1))
       else
        p=pi*0.5
       endif
       p=atan2(nl(2),nl(1))
       if (p.lt.0.0) p=pi+p
       call EToG(e1,g_1)
       call EToG(e2,g_2)
       do i2=1,nsymm
		  call symop(O,i2,so_1)!call the i2-th 3x3 symmetry operator
		  call MToM (so_1,g_1,gg_1)!hit the first g  with the i2-th symmetry operator
          call MToM (so_1,g_2,gg_2)
          do i4=1,int(HCD)
           t=acos2((2.0*(float(i4)-0.5)/HCD)-1.0)
           call AnglesToV (t, p, n)
           call MToV (gg_1, n, nf)
           call VToAngles (nf, sa)
           c4=int(float(CD2)*cos(sa(1)))+1
           c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
           if (c4.eq.CD2+1) c4=CD2
           if (c4.eq.0) c4=1
           if (c5.eq.4*CD2+1) c5=4*CD2
           pd(c4,c5)=pd(c4,c5)+length
          enddo ! ends the i4=1,int(HCD) loop
          do i4=1,int(HCD)
           t=acos2((2.0*(float(i4)-0.5)/HCD)-1.0)
           call AnglesToV (t, p, n)
           call MToV (gg_2, n, nf)
           call VToAngles (nf, sa)
           c4=int(float(CD2)*cos(sa(1)))+1
           c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
           if (c4.eq.CD2+1) c4=CD2
           if (c4.eq.0) c4=1
           if (c5.eq.4*CD2+1) c5=4*CD2
           pd(c4,c5)=pd(c4,c5)+length
          enddo ! ends the i4=1,int(HCD) loop
       enddo ! end the i2=1,nsymm loop
       !notify user of the progress of the calculation
       if (mod(i1,10000).eq.0) then
        write(6,"(I9,A)") i1,' line segments classified'
       endif
 2050  continue
      enddo ! this closes the i1=header+1,nnline loop
      close(30)
      !determine the average line length in each bin
      avg=0.0
      do i1=1,CD2
       do i2=1,4*CD2
        avg=avg+pd(i1,i2)
       enddo
      enddo
      avg=avg/float(4*CD2*CD2)
      !subtract the background using the mean field approximation
      pnc=(2.0*float(CD2)-1.0)/(2.0*float(CD2))
      f=(1.0/float(CD2))*(pnc/(1.0-pnc))
      do i1=1,CD2
       do i2=1,4*CD2
        tem=pd(i1,i2)
        tem=(tem-avg*(pnc-f+f*pnc))/(1.0+f)
        if (tem.gt.0.0) then
         pd(i1,i2)=tem
        else
         pd(i1,i2)=0.0
        endif
       enddo
      enddo
      !determine the average value of the distribution
      avg=0.0
      do i1=1,CD2
       do i2=1,4*CD2
        avg=avg+pd(i1,i2)
       enddo
      enddo
      avg=avg/float(4*CD2*CD2)
      !write the values in MRD units
      do i1=1,CD2
       do i2=1,4*CD2
        pd(i1,i2)=pd(i1,i2)/avg
       write (41,*) pd(i1,i2)
       enddo
      enddo
      close(41)
      write(6,1)'2D grain boundary plane distribution complete'




      !------------------------------------------------------------------
      ! Begin the section of code that computes the five-dimensional
      ! grain boundary character distribution using stereology.
      !------------------------------------------------------------------

      !Note: this is based on calc_gbcd_stereo_fd, version 07/27/2013.
      !On 03/28/22, this program exactly reproduced the GBCD calculated
      !by calc_gbcd_stereo_fd using all_IN100_segments.txt as input.
3000  if (out(3).ne.1) goto 4000
      write(6,1)'calculating the grain boundary character distribution ...'
      ref(1)=1.0
      ref(2)=0.0
      mag_ref=1.0
      t_out=0.0
      open(42, file=keyword(1:point)//'_gbcd.txt',status='unknown')
      !open the file with the data
      open(30, file=trim(keyword), status='old')
       write(42,2)'# written by calc_discrete_dist, ',version
       write(42,2)'# based on data in file: ',trim(keyword)
       write(42,8)'# number of columns = ',NCol,', symmetry = ',msym
       write(42,9)'# ref. frame = ',rotation,', radians = ',radian,', StepSize = ',StepSize
       write(42,10)'# Discretization = ',CD,', ',CD2
       write(42,11)'# output flags: ',out(1),out(2),out(3), out(4),out(5),out(6),out(7)
       write(42,2)'# note: ',trim(comment)
       write(42,1)'# the following lines are the header of the source file'
      do i1=1,header
       read(30,1) inline ! read a line of the header
       write(42,1)inline ! write that line into the output file
      enddo! closes i1=1,header
      gbd=0.0
      do i1=header+1,nnline
       if (NCol.eq.10) then
        read(30,*)e1(1),e1(2),e1(3),e2(1),e2(2),e2(3),x1,y1,x2,y2
       else
        read(30,*)e1(1),e1(2),e1(3),e2(1),e2(2),e2(3),tem,tem,x1,y1,x2,y2
       endif
       if (e1(1).eq.0.00.and.e1(2).eq.0.00.and.e1(3).eq.0.00) then
        goto 3050 ! Sends it to the next line
       endif
       if (e2(1).eq.0.00.and.e2(2).eq.0.00.and.e2(3).eq.0.00) then
        goto 3050 ! Sends it to the next line
       endif
       ! If the Euler angles are provided in units of degrees, these
       ! lines convert them to radians
       if (radian.eq.0) then
        call DToRad3(e1, e1r)
        e1=e1r
        call DToRad3(e2, e2r)
        e2=e2r
       else
        continue
       endif
       !If the data are in the HKL reference frame, and the system is hexagonal,
       !then the data are rotated by -30° so that [2-1-10] is parallel to +x in
       !the output (this direction is to the right in the plane of the paper)
       if (rotation.eq.0.and.msym.eq.2) then
        if (e1(3).gt.(30*pi/180)) then
         e1(3) = e1(3)-(30*pi/180)
        else
         e1(3) = e1(3)+(330*pi/180)
        endif
        if (e2(3).gt.(30*pi/180)) then
         e2(3) = e2(3)-(30*pi/180)
        else
         e2(3) = e2(3)+(330*pi/180)
        endif
       endif
       !This section determines the in plane components of the grain boundary segment,
       !based on the endpoints in the data file.  The first calculation of nl(i)
       !transforms line segments into the same reference frame as the Euler angles,
       !based on and assumption of default settings in TSL or HKL software.  The second
       !calculation determines the normal to the grain boundary trace.
       if (rotation.eq.1) then       !TSL
        nl(1)=y1-y2
        nl(2)=x1-x2
       endif
       if (rotation.eq.0) then       !HKL
        nl(1)=x2-x1
        nl(2)=-(y2-y1)
       endif
       if (rotation.eq.2) then       !D3D/natural
        nl(1)=x2-x1
        nl(2)=y2-y1
       endif
       templ=nl(2)                   !normal to the GB trace
       nl(2)=nl(1)
       nl(1)=-templ
       nl(3)=0.0
       length=sqrt(nl(1)*nl(1)+nl(2)*nl(2))

       !This section calculates phi for the GB normal, with respect to [100]
       mag_dir=((nl(1)**2)+(nl(2)**2))**0.5
       dot = (nl(1)*ref(1))+(nl(2)*ref(2))
       p = acos2(dot/(mag_dir*mag_ref))
       if (nl(2).lt.0.0) p=-p  !vectors pointing below the x-z plane assigned
                               ! negative angles (acos make them positive) so
                               ! the range of phi is 0 - 180 for a positive y
                               ! coorinate and 0 - -180 for a negative y coord
       if (p.lt.0.0) p=pi+p    ! For all directions with negative y (negative angle)
							   ! we consider the opposite direction
       !convert the Euler angles from the input file to a g matrix
       call EToG(e1,g_1)
       call EToG(e2,g_2)
       do i2=1,nsymm
        do i3=1,nsymm
		  call symop(O,i2,so_1)!call the i2-th 3x3 symmetry operator
		  call MToM (so_1,g_1,gg_1)!hit the first g  with the i2-th symmetry operator
		  call symop(O,i3,so_2)!call the i3-th 3x3 symmetry operator
		  call MToM(so_2,g_2,gg_2)!hit the second g  with the i3-th symmetry operator
		  call trans(gg_2,gg_2t)!transpose the second g matrix
		  call MToM(gg_1,gg_2t,dg)!misorientation between g_1 and g_2
		  call GToE(dg,mort) !get  Euler angles for this misorientation
          !determine the indices the GBCD matrix (c1,c2,c3)
          if (msym.eq.3) then
		   if (mort(1).lt.(pi*0.5+eps).and.mort(2).lt.(pi*0.5+eps).and.mort(3).lt.(pi*0.5+eps) ) then
     		c1=int(float(CD)*2.0*mort(1)/pi)+1
			c2=int(float(CD)*cos(mort(2)))+1
			c3=int(float(CD)*2.0*mort(3)/pi)+1
			if (c1.eq.CD+1) c1=CD
			if (c2.eq.CD+1) c2=CD
			if (c3.eq.CD+1) c3=CD
            !puts line segments into bins cooresponding to the orientation of the
            !vector normal to the trace.  Consider orientations in 1° increments, from 1 to 90.
            do i4=1,int(HCD)
             t=acos2((2.0*(float(i4)-0.5)/HCD)-1.0)
             call AnglesToV (t, p, n)
             call MToV (gg_1, n, nf)
             call VToAngles (nf, sa)
             c4=int(float(CD2)*cos(sa(1)))+1
             c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
             if (c4.eq.CD2+1) c4=CD2
             if (c4.eq.0) c4=1
             if (c5.eq.4*CD2+1) c5=4*CD2
             gbd(c1,c2,c3,c4,c5)=gbd(c1,c2,c3,c4,c5)+length
            enddo ! end the i4=1,90 loop
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
           !puts line segments into bins cooresponding to the orientation of the
           !vector normal to the trace.  Consider orientations in 1° increments, from 1 to 90.
           do i4=1,int(HCD)
            t=acos2((2.0*(float(i4)-0.5)/HCD)-1.0)
            call AnglesToV (t, p, n)
            call MToV (gg_1, n, nf)
            call VToAngles (nf, sa)
            c4=int(float(CD2)*cos(sa(1)))+1
            c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
            if (c4.eq.CD2+1) c4=CD2
            if (c4.eq.0) c4=1
            if (c5.eq.4*CD2+1) c5=4*CD2
            gbd(c1,c2,c3,c4,c5)=gbd(c1,c2,c3,c4,c5)+length
           enddo ! end the i4=1,90 loop
          endif
          !implement crystal exchange symmetry, and repeat the process above
		  call trans (gg_1, gg_1t)!transpose the fist g matrix
		  call MToM (gg_2, gg_1t, dg)!misorientation between g2 and g1
		  call GToE (dg, mort) !get  Euler angles for this misorientation
          !determine the indices the GBCD matrix (c1,c2,c3)
          if (msym.eq.3) then
		   if (mort(1).lt.(pi*0.5+eps).and.mort(2).lt.(pi*0.5+eps).and.mort(3).lt.(pi*0.5+eps)) then
     		c1=int(float(CD)*2.0*mort(1)/pi)+1
			c2=int(float(CD)*cos(mort(2)))+1
			c3=int(float(CD)*2.0*mort(3)/pi)+1
			if (c1.eq.CD+1) c1=CD
			if (c2.eq.CD+1) c2=CD
			if (c3.eq.CD+1) c3=CD
            do i4=1,int(HCD)
             t=acos2((2.0*(float(i4)-0.5)/HCD)-1.0)
             call AnglesToV (t, p, n)
             call MToV (gg_2, n, nf)
             call VToAngles (nf, sa)
             c4=int(float(CD2)*cos(sa(1)))+1
             c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
             if (c4.eq.CD2+1) c4=CD2
             if (c4.eq.0) c4=1
             if (c5.eq.4*CD2+1) c5=4*CD2
             gbd(c1,c2,c3,c4,c5)=gbd(c1,c2,c3,c4,c5)+length
            enddo ! end the i4=1,90 loop
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
           do i4=1,int(HCD)
            t=acos2((2.0*(float(i4)-0.5)/HCD)-1.0)
            call AnglesToV (t, p, n)
            call MToV (gg_2, n, nf)
            call VToAngles (nf, sa)
            c4=int(float(CD2)*cos(sa(1)))+1
            c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
            if (c4.eq.CD2+1) c4=CD2
            if (c4.eq.0) c4=1
            if (c5.eq.4*CD2+1) c5=4*CD2
            gbd(c1,c2,c3,c4,c5)=gbd(c1,c2,c3,c4,c5)+length
           enddo ! end the i4=1,90 loop
          endif
         enddo ! end the i3=1,nsymm loop
        enddo ! end the i2=1,nsymm loop
        !notify user of the progress of the calculation
		if (mod(i1,20000).eq.0) then
	     write(6,"(I9,A)") i1,' line segments classified'
	    endif
 3050   continue
      enddo ! this closes the i1=header+1,nnline loop
      close(30)

      !The size of the cubic distribution is smaller than for
      !other symmetries
      if (msym.eq.3) then
      !determine the average line length in each bin
      avg=0.0
      do i1=1,CD
       do i2=1,CD
        do i3=1,CD
         do i4=1,CD2
          do i5=1,4*CD2
           avg=avg+gbd(i1,i2,i3,i4,i5)
          enddo
         enddo
        enddo
       enddo
      enddo
      avg=avg/float(4*CD*CD*CD*CD2*CD2)

      !background subtraction by the mean field aprpoximation
      pnc=(2.0*float(CD2)-1.0)/(2.0*float(CD2))
      f=(1.0/float(CD2))*(pnc/(1.0-pnc))
      do i1=1,CD
       do i2=1,CD
        do i3=1,CD
         avg=0.0
         do i4=1,CD2
          do i5=1,4*CD2
           avg=avg+gbd(i1,i2,i3,i4,i5)
          enddo
         enddo
         avg=avg/float(4*CD2*CD2)
         do i4=1,CD2
          do i5=1,4*CD2
           tem=gbd(i1,i2,i3,i4,i5)
           tem=(tem-avg*(pnc-f+f*pnc))/(1.0+f)
           if (tem.gt.0.0) then
            gbd(i1,i2,i3,i4,i5)=tem
           else
            gbd(i1,i2,i3,i4,i5)=0.00
           endif
          enddo
         enddo
        enddo
       enddo
      enddo

      !find the new average value of the distribution
      avg=0.0
      do i1=1,CD
       do i2=1,CD
        do i3=1,CD
         do i4=1,CD2
          do i5=1,4*CD2
           avg=avg+gbd(i1,i2,i3,i4,i5)
          enddo
         enddo
        enddo
       enddo
      enddo
      avg=avg/float(4*CD*CD*CD*CD2*CD2)

      !write the MRD values of the distribution to the gbcd file
      do i1=1,CD
       do i2=1,CD
        do i3=1,CD
         do i4=1,CD2
          do i5=1,4*CD2
           gbd(i1,i2,i3,i4,i5)=gbd(i1,i2,i3,i4,i5)/avg
           write (42,*) gbd(i1,i2,i3,i4,i5)
          enddo
         enddo
        enddo
       enddo
      enddo

      else

      !determine the average line length in each bin
      avg=0.0
      do i1=1,CD*4
       do i2=1,CD*2
        do i3=1,CD*4
         do i4=1,CD2
          do i5=1,4*CD2
           avg=avg+gbd(i1,i2,i3,i4,i5)
          enddo
         enddo
        enddo
       enddo
      enddo
      avg=avg/float(4*4*4*2*CD*CD*CD*CD2*CD2)

      !background subtraction by the mean field aprpoximation
      pnc=(2.0*float(CD2)-1.0)/(2.0*float(CD2))
      f=(1.0/float(CD2))*(pnc/(1.0-pnc))
      do i1=1,CD*4
       do i2=1,CD*2
        do i3=1,CD*4
         avg=0.0
         do i4=1,CD2
          do i5=1,4*CD2
           avg=avg+gbd(i1,i2,i3,i4,i5)
          enddo
         enddo
         avg=avg/float(4*CD2*CD2)
         do i4=1,CD2
          do i5=1,4*CD2
           tem=gbd(i1,i2,i3,i4,i5)
           tem=(tem-avg*(pnc-f+f*pnc))/(1.0+f)
           if (tem.gt.0.0) then
            gbd(i1,i2,i3,i4,i5)=tem
           else
            gbd(i1,i2,i3,i4,i5)=0.00
           endif
          enddo
         enddo
        enddo
       enddo
      enddo

      !find the new average value of the distribution
      avg=0.0
      do i1=1,CD*4
       do i2=1,CD*2
        do i3=1,CD*4
         do i4=1,CD2
          do i5=1,4*CD2
           avg=avg+gbd(i1,i2,i3,i4,i5)
          enddo
         enddo
        enddo
       enddo
      enddo
      avg=avg/float(4*4*4*2*CD*CD*CD*CD2*CD2)

      !write the MRD values of the distribution to the gbcd file
      do i1=1,CD*4
       do i2=1,CD*2
        do i3=1,CD*4
         do i4=1,CD2
          do i5=1,4*CD2
           gbd(i1,i2,i3,i4,i5)=gbd(i1,i2,i3,i4,i5)/avg
           write (42,*) gbd(i1,i2,i3,i4,i5)
          enddo
         enddo
        enddo
       enddo
      enddo
      endif
      close (42)




      !------------------------------------------------------------------
      ! Begin the section of code that bins three-dimensional data into
      ! a discrete grain boundary plane distribution (2D).
      !------------------------------------------------------------------



 4000 if (out(4).ne.1) goto 5000
      if (NCol.ne.10) then
       write(6,1)'to bin 3D data, the number of columns (NCol) must equal 1.  bye'
       goto 9000
      endif
      write(6,1)'calculating the 2D grain boundary plane distribution ...'
      open(43, file=keyword(1:point)//'_3D_gbpd.txt',status='unknown')
      !open the file with the data
      open(30, file=trim(keyword), status='old')
       write(43,2)'# written by calc_discrete_dist, ',version
       write(43,2)'# based on data in file: ',trim(keyword)
       write(43,8)'# number of columns = ',NCol,', symmetry = ',msym
       write(43,9)'# ref. frame = ',rotation,', radians = ',radian,', StepSize = ',StepSize
       write(43,10)'# Discretization = ',CD,', ',CD2
       write(43,11)'# output flags: ',out(1),out(2),out(3), out(4),out(5),out(6),out(7)
       write(43,2)'# note: ',trim(comment)
       write(43,1)'# the following lines are the header of the source file'
      do i1=1,header
       read(30,1) inline ! read a line of the header
       write(43,1)inline ! write that line into the output file
      enddo! closes i1=1,header
      pd=0.0
      do i1=header+1,nnline
       read(30,*)e1(1),e1(2),e1(3),e2(1),e2(2),e2(3),nl(1),nl(2),nl(3),area
       if (e1(1).eq.0.00.and.e1(2).eq.0.00.and.e1(3).eq.0.00) then
        goto 4050 ! Sends it to the next line
       endif
       if (e2(1).eq.0.00.and.e2(2).eq.0.00.and.e2(3).eq.0.00) then
        goto 4050 ! Sends it to the next line
       endif

       ! If the Euler angles are provided in units of degrees, these
       ! lines convert them to radians
       if (radian.eq.0) then
        call DToRad3(e1, e1r)
        e1=e1r
        call DToRad3(e2, e2r)
        e2=e2r
       else
        continue
       endif
       !If the data are in the HKL reference frame, and the system is hexagonal,
       !then the data are rotated by -30° so that [2-1-10] is parallel to +x in
       !the output (this direction is to the right in the plane of the paper)
       if (rotation.eq.0.and.msym.eq.2) then
        if (e1(3).gt.(30*pi/180)) then
         e1(3) = e1(3)-(30*pi/180)
        else
         e1(3) = e1(3)+(330*pi/180)
        endif
        if (e2(3).gt.(30*pi/180)) then
         e2(3) = e2(3)-(30*pi/180)
        else
         e2(3) = e2(3)+(330*pi/180)
        endif
       endif

       !convert the Euler angles from the input file to a g matrix
       call EToG(e1,g_1)
       call EToG(e2,g_2)
       do i2=1,nsymm
		  call symop(O,i2,so_1)!call the i2-th 3x3 symmetry operator
		  call MToM (so_1,g_1,gg_1)!hit the first g  with the i2-th symmetry operator
          call MToM (so_1,g_2,gg_2)
          call MToV (gg_1, nl, nf)
          call VToAngles (nf, sa)
          c4=int(float(CD2)*cos(sa(1)))+1
          c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
          if (c4.eq.CD2+1) c4=CD2
          if (c4.eq.0) c4=1
          if (c5.eq.4*CD2+1) c5=4*CD2
          pd(c4,c5)=pd(c4,c5)+area
          call MToV (gg_2, nl, nf)
          call VToAngles (nf, sa)
          c4=int(float(CD2)*cos(sa(1)))+1
          c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
          if (c4.eq.CD2+1) c4=CD2
          if (c4.eq.0) c4=1
          if (c5.eq.4*CD2+1) c5=4*CD2
          pd(c4,c5)=pd(c4,c5)+area
       enddo ! end the i2=1,nsymm loop
       !notify user of the progress of the calculation
       if (mod(i1,50000).eq.0) then
        write(6,"(I9,A)") i1,' line segments classified'
       endif
 4050  continue
      enddo ! this closes the i1=header+1,nnline loop
      close(30)
      !determine the average value of the distribution
      avg=0.0
      do i1=1,CD2
       do i2=1,4*CD2
        avg=avg+pd(i1,i2)
       enddo
      enddo
      avg=avg/float(4*CD2*CD2)
      !write the values in MRD units
      do i1=1,CD2
       do i2=1,4*CD2
        pd(i1,i2)=pd(i1,i2)/avg
       write (43,*) pd(i1,i2)
       enddo
      enddo
      close(43)














      !------------------------------------------------------------------
      ! Begin the section of code that bins three-dimensional data into
      ! a discrete grain boundary character distribution.
      !------------------------------------------------------------------

      !Note: this is based on calc_gbcd_fd (Calc GBCD 3D full domain),
      !version 07/18/2018.

5000  if (out(5).ne.1) goto 6000
      if (NCol.ne.10) then
       write(6,1)'to bin 3D data, the number of columns (NCol) must equal 1.  bye'
       goto 9000
      endif
      write(6,1)'calculating the grain boundary character distribution from 3D data ...'
      open(44, file=keyword(1:point)//'_3D_gbcd.txt',status='unknown')
      !open the file with the data
      open(30, file=trim(keyword), status='old')
       write(44,2)'# written by calc_discrete_dist, ',version
       write(44,2)'# based on data in file: ',trim(keyword)
       write(44,8)'# number of columns = ',NCol,', symmetry = ',msym
       write(44,9)'# ref. frame = ',rotation,', radians = ',radian,', StepSize = ',StepSize
       write(44,10)'# Discretization = ',CD,', ',CD2
       write(44,11)'# output flags: ',out(1),out(2),out(3), out(4),out(5),out(6),out(7)
       write(44,2)'# note: ',trim(comment)
       write(44,1)'# the following lines are the header of the source file'
      gbd=0.0
      do i1=1,header
       read(30,1) inline ! read a line of the header
       write(44,1)inline ! write that line into the output file
      enddo! closes i1=1,header
      do i1=header+1,nnline
       read(30,*)e1(1),e1(2),e1(3),e2(1),e2(2),e2(3),nl(1),nl(2),nl(3),area
       if (e1(1).eq.0.00.and.e1(2).eq.0.00.and.e1(3).eq.0.00) then
        goto 5050 ! Sends it to the next line
       endif
       if (e2(1).eq.0.00.and.e2(2).eq.0.00.and.e2(3).eq.0.00) then
        goto 5050 ! Sends it to the next line
       endif
       ! If the Euler angles are provided in units of degrees, these
       ! lines convert them to radians
       if (radian.eq.0) then
        call DToRad3(e1, e1r)
        e1=e1r
        call DToRad3(e2, e2r)
        e2=e2r
       else
        continue
       endif
       !If the data are in the HKL reference frame, and the system is hexagonal,
       !then the data are rotated by -30° so that [2-1-10] is parallel to +x in
       !the output (this direction is to the right in the plane of the paper)
       if (rotation.eq.0.and.msym.eq.2) then
        if (e1(3).gt.(30*pi/180)) then
         e1(3) = e1(3)-(30*pi/180)
        else
         e1(3) = e1(3)+(330*pi/180)
        endif
        if (e2(3).gt.(30*pi/180)) then
         e2(3) = e2(3)-(30*pi/180)
        else
         e2(3) = e2(3)+(330*pi/180)
        endif
       endif

       !convert the Euler angles from the input file to a g matrix
       call EToG(e1,g_1)
       call EToG(e2,g_2)
       do i2=1,nsymm
        do i3=1,nsymm
		  call symop(O,i2,so_1)!call the i2-th 3x3 symmetry operator
		  call MToM (so_1,g_1,gg_1)!hit the first g  with the i2-th symmetry operator
		  call symop(O,i3,so_2)!call the i3-th 3x3 symmetry operator
		  call MToM(so_2,g_2,gg_2)!hit the second g  with the i3-th symmetry operator
		  call trans(gg_2,gg_2t)!transpose the second g matrix
		  call MToM(gg_1,gg_2t,dg)!misorientation between g_1 and g_2
		  call GToE(dg,mort) !get  Euler angles for this misorientation
          if (msym.eq.3) then
 		   if (mort(1).lt.(pi*0.5+eps).and.mort(2).lt.(pi*0.5+eps).and.mort(3).lt.(pi*0.5+eps)) then
     		c1=int(float(CD)*2.0*mort(1)/pi)+1
			c2=int(float(CD)*cos(mort(2)))+1
			c3=int(float(CD)*2.0*mort(3)/pi)+1
			if (c1.eq.CD+1) c1=CD
			if (c2.eq.CD+1) c2=CD
			if (c3.eq.CD+1) c3=CD
            !puts line segments into bin cooresponding to the orientation of the
            !vector normal.
            call MToV (gg_1, nl, nf)
            call VToAngles (nf, sa)
             c4=int(float(CD2)*cos(sa(1)))+1
            c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
            if (c4.eq.CD2+1) c4=CD2
            if (c4.eq.0) c4=1
            if (c5.eq.4*CD2+1) c5=4*CD2
            gbd(c1,c2,c3,c4,c5)=gbd(c1,c2,c3,c4,c5)+area
           endif
          else
           !determine the indices the GBCD matrix (c1,c2,c3)
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
           !puts line segments into bin cooresponding to the orientation of the
           !vector normal.
           call MToV (gg_1, nl, nf)
           call VToAngles (nf, sa)
           c4=int(float(CD2)*cos(sa(1)))+1
           c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
           if (c4.eq.CD2+1) c4=CD2
           if (c4.eq.0) c4=1
           if (c5.eq.4*CD2+1) c5=4*CD2
           gbd(c1,c2,c3,c4,c5)=gbd(c1,c2,c3,c4,c5)+area
          endif
          !implement crystal exchange symmetry, and repeat the process above
		  call trans (gg_1, gg_1t)!transpose the fist g matrix
		  call MToM (gg_2, gg_1t, dg)!misorientation between g2 and g1
		  call GToE (dg, mort) !get  Euler angles for this misorientation
          if (msym.eq.3) then
		   if (mort(1).lt.(pi*0.5+eps).and.mort(2).lt.(pi*0.5+eps).and.mort(3).lt.(pi*0.5+eps)) then
     		c1=int(float(CD)*2.0*mort(1)/pi)+1
			c2=int(float(CD)*cos(mort(2)))+1
			c3=int(float(CD)*2.0*mort(3)/pi)+1
			if (c1.eq.CD+1) c1=CD
			if (c2.eq.CD+1) c2=CD
			if (c3.eq.CD+1) c3=CD
            call MToV (gg_2, nl, nf)
		    call VToAngles (nf, sa)
            c4=int(float(CD2)*cos(sa(1)))+1
            c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
            if (c4.eq.CD2+1) c4=CD2
            if (c4.eq.0) c4=1
            if (c5.eq.4*CD2+1) c5=4*CD2
            gbd(c1,c2,c3,c4,c5)=gbd(c1,c2,c3,c4,c5)+area
           endif
          else
           !determine the indices the GBCD matrix (c1,c2,c3)
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
           call MToV (gg_2, nl, nf)
		   call VToAngles (nf, sa)
           c4=int(float(CD2)*cos(sa(1)))+1
           c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
           if (c4.eq.CD2+1) c4=CD2
           if (c4.eq.0) c4=1
           if (c5.eq.4*CD2+1) c5=4*CD2
           gbd(c1,c2,c3,c4,c5)=gbd(c1,c2,c3,c4,c5)+area
          endif
         enddo ! end the i3=1,nsymm loop
        enddo ! end the i2=1,nsymm loop
        !notify user of the progress of the calculation
		if (mod(i1,50000).eq.0) then
	     write(6,"(I9,A)") i1,' line segments classified'
	    endif
 5050   continue
      enddo ! this closes the i1=header+1,nnline loop
      close(30)

      if (msym.eq.3) then
      !determine the average line length in each bin
      avg=0.0
      do i1=1,CD
       do i2=1,CD
        do i3=1,CD
         do i4=1,CD2
          do i5=1,4*CD2
           avg=avg+gbd(i1,i2,i3,i4,i5)
          enddo
         enddo
        enddo
       enddo
      enddo
      avg=avg/float(4*CD*CD*CD*CD2*CD2)

      !write the MRD values of the distribution to the gbpd file
      do i1=1,CD
       do i2=1,CD
        do i3=1,CD
         do i4=1,CD2
          do i5=1,4*CD2
           gbd(i1,i2,i3,i4,i5)=gbd(i1,i2,i3,i4,i5)/avg
           write (44,*) gbd(i1,i2,i3,i4,i5)
          enddo
         enddo
        enddo
       enddo
      enddo

      else

      !determine the average line length in each bin
      avg=0.0
      do i1=1,CD*4
       do i2=1,CD*2
        do i3=1,CD*4
         do i4=1,CD2
          do i5=1,4*CD2
           avg=avg+gbd(i1,i2,i3,i4,i5)
          enddo
         enddo
        enddo
       enddo
      enddo
      avg=avg/float(4*4*4*2*CD*CD*CD*CD2*CD2)

      !write the MRD values of the distribution to the gbpd file
      do i1=1,CD*4
       do i2=1,CD*2
        do i3=1,CD*4
         do i4=1,CD2
          do i5=1,4*CD2
           gbd(i1,i2,i3,i4,i5)=gbd(i1,i2,i3,i4,i5)/avg
           write (44,*) gbd(i1,i2,i3,i4,i5)
          enddo
         enddo
        enddo
       enddo
      enddo


      endif

      close (44)
      goto 6000


      !------------------------------------------------------------------
      ! Begin the section of code that bins line segment data into a
      ! discrete grain boundary plane distribution, as if it is
      ! grain boundary plane data, under the assumption that the
      ! microstructure is purely columnar.
      !------------------------------------------------------------------

      !Note: this is based on calc_gbcd_col (gbcd_columnar),
      !version 05/18/2013.

6000  if (out(6).ne.1) goto 7000
      write(6,1)'calculating the 2D GBPD assuming a columnar microstructure ...'
      pd=0.0
      open(45, file=keyword(1:point)//'_col_gbpd.txt',status='unknown')
      !open the file with the data
      open(30, file=trim(keyword), status='old')
       write(45,2)'# written by calc_discrete_dist, ',version
       write(45,2)'# based on data in file: ',trim(keyword)
       write(45,8)'# number of columns = ',NCol,', symmetry = ',msym
       write(45,9)'# ref. frame = ',rotation,', radians = ',radian,', StepSize = ',StepSize
       write(45,10)'# Discretization = ',CD,', ',CD2
       write(45,11)'# output flags: ',out(1),out(2),out(3), out(4),out(5),out(6),out(7)
       write(45,2)'# note: ',trim(comment)
       write(45,1)'# the following lines are the header of the source file'
      pd=0.0
      do i1=1,header
       read(30,1) inline ! read a line of the header
       write(45,1)inline ! write that line into the output file
      enddo! closes i1=1,header
      do i1=header+1,nnline
       if (NCol.eq.10) then
        read(30,*)e1(1),e1(2),e1(3),e2(1),e2(2),e2(3),x1,y1,x2,y2
       else
        read(30,*)e1(1),e1(2),e1(3),e2(1),e2(2),e2(3),tem,tem,x1,y1,x2,y2
       endif
       if (e1(1).eq.0.00.and.e1(2).eq.0.00.and.e1(3).eq.0.00) then
        goto 6050 ! Sends it to the next line
       endif
       if (e2(1).eq.0.00.and.e2(2).eq.0.00.and.e2(3).eq.0.00) then
        goto 6050 ! Sends it to the next line
       endif

       ! If the Euler angles are provided in units of degrees, these
       ! lines convert them to radians
       if (radian.eq.0) then
        call DToRad3(e1, e1r)
        e1=e1r
        call DToRad3(e2, e2r)
        e2=e2r
       else
        continue
       endif
       !If the data are in the HKL reference frame, and the system is hexagonal,
       !then the data are rotated by -30° so that [2-1-10] is parallel to +x in
       !the output (this direction is to the right in the plane of the paper)
       if (rotation.eq.0.and.msym.eq.2) then
        if (e1(3).gt.(30*pi/180)) then
         e1(3) = e1(3)-(30*pi/180)
        else
         e1(3) = e1(3)+(330*pi/180)
        endif
        if (e2(3).gt.(30*pi/180)) then
         e2(3) = e2(3)-(30*pi/180)
        else
         e2(3) = e2(3)+(330*pi/180)
        endif
       endif
       !This section determines the in plane components of the grain boundary segment,
       !based on the endpoints in the data file.  The first calculation of nl(i)
       !transforms line segments into the same reference frame as the Euler angles,
       !based on and assumption of default settings in TSL or HKL software.  The second
       !calculation determines the normal to the grain boundary trace.
       if (rotation.eq.1) then       !TSL
        nl(1)=y1-y2
        nl(2)=x1-x2
       endif
       if (rotation.eq.0) then       !HKL
        nl(1)=x2-x1
        nl(2)=-(y2-y1)
       endif
       if (rotation.eq.2) then       !D3D/natural
        nl(1)=x2-x1
        nl(2)=y2-y1
       endif
       templ=nl(2)                   !normal to the GB trace
       nl(2)=nl(1)
       nl(1)=-templ
       nl(3)=0.0
       length=sqrt(nl(1)*nl(1)+nl(2)*nl(2))
       if (nl(1).ne.0.0) then
        p=atan(nl(2)/nl(1))
       else
        p=pi*0.5
       endif
       p=atan2(nl(2),nl(1))
       if (p.lt.0.0) p=pi+p
       call EToG(e1,g_1)
       call EToG(e2,g_2)
       do i2=1,nsymm
		  call symop(O,i2,so_1)!call the i2-th 3x3 symmetry operator
		  call MToM (so_1,g_1,gg_1)!hit the first g  with the i2-th symmetry operator
          call MToM (so_1,g_2,gg_2)
          t=pi/2.0 !for the columnar approximation, theta is fixed
          call AnglesToV (t, p, n)
          call MToV (gg_1, n, nf)
          call VToAngles (nf, sa)
          c4=int(float(CD2)*cos(sa(1)))+1
          c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
          if (c4.eq.CD2+1) c4=CD2
          if (c4.eq.0) c4=1
          if (c5.eq.4*CD2+1) c5=4*CD2
          pd(c4,c5)=pd(c4,c5)+length
          t=pi/2.0
          call AnglesToV (t, p, n)
          call MToV (gg_2, n, nf)
          call VToAngles (nf, sa)
          c4=int(float(CD2)*cos(sa(1)))+1
          c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
          if (c4.eq.CD2+1) c4=CD2
          if (c4.eq.0) c4=1
          if (c5.eq.4*CD2+1) c5=4*CD2
          pd(c4,c5)=pd(c4,c5)+length
       enddo ! end the i2=1,nsymm loop
       !notify user of the progress of the calculation
       if (mod(i1,5000).eq.0) then
        write(6,"(I9,A)") i1,' line segments classified'
       endif
 6050  continue
      enddo ! this closes the i1=header+1,nnline loop
      close(30)
      !determine the average line length in each bin
      avg=0.0
      do i1=1,CD2
       do i2=1,4*CD2
        avg=avg+pd(i1,i2)
       enddo
      enddo
      avg=avg/float(4*CD2*CD2)
      !write the values in MRD units
      do i1=1,CD2
       do i2=1,4*CD2
        pd(i1,i2)=pd(i1,i2)/avg
       write (45,*) pd(i1,i2)
       enddo
      enddo
      close(45)
      write(6,1)'grain boundary plane distribution complete'















      !------------------------------------------------------------------
      ! Begin the section of code that bins line segment data into a
      ! discrete grain boundary character distribution, as if it is
      ! grain boundary plane data, under the assumption that the
      ! microstructure is purely columnar.
      !------------------------------------------------------------------

      !Note: this is based on calc_gbcd_col (gbcd_columnar),
      !version 05/18/2013.

7000  if (out(7).ne.1) goto 9000
      write(6,1)'calculating the GBCD using the columnar assumption ...'
      open(46, file=keyword(1:point)//'_col_gbcd.txt',status='unknown')
      !open the file with the data
      open(30, file=trim(keyword), status='old')
       write(46,2)'# written by calc_discrete_dist, ',version
       write(46,2)'# based on data in file: ',trim(keyword)
       write(46,8)'# number of columns = ',NCol,', symmetry = ',msym
       write(46,9)'# ref. frame = ',rotation,', radians = ',radian,', StepSize = ',StepSize
       write(46,10)'# Discretization = ',CD,', ',CD2
       write(46,11)'# output flags: ',out(1),out(2),out(3), out(4),out(5),out(6),out(7)
       write(46,2)'# note: ',trim(comment)
       write(46,1)'# the following lines are the header of the source file'
      gbd=0.0
      do i1=1,header
       read(30,1) inline ! read a line of the header
       write(46,1)inline ! write that line into the output file
      enddo! closes i1=1,header
      do i1=header+1,nnline
       if (NCol.eq.10) then
        read(30,*)e1(1),e1(2),e1(3),e2(1),e2(2),e2(3),x1,y1,x2,y2
       else
        read(30,*)e1(1),e1(2),e1(3),e2(1),e2(2),e2(3),tem,tem,x1,y1,x2,y2
       endif
       if (e1(1).eq.0.00.and.e1(2).eq.0.00.and.e1(3).eq.0.00) then
        goto 7050 ! Sends it to the next line
       endif
       if (e2(1).eq.0.00.and.e2(2).eq.0.00.and.e2(3).eq.0.00) then
        goto 7050 ! Sends it to the next line
       endif

       ! If the Euler angles are provided in units of degrees, these
       ! lines convert them to radians
       if (radian.eq.0) then
        call DToRad3(e1, e1r)
        e1=e1r
        call DToRad3(e2, e2r)
        e2=e2r
       else
        continue
       endif
       !If the data are in the HKL reference frame, and the system is hexagonal,
       !then the data are rotated by -30° so that [2-1-10] is parallel to +x in
       !the output (this direction is to the right in the plane of the paper)
       if (rotation.eq.0.and.msym.eq.2) then
        if (e1(3).gt.(30*pi/180)) then
         e1(3) = e1(3)-(30*pi/180)
        else
         e1(3) = e1(3)+(330*pi/180)
        endif
        if (e2(3).gt.(30*pi/180)) then
         e2(3) = e2(3)-(30*pi/180)
        else
         e2(3) = e2(3)+(330*pi/180)
        endif
       endif
       !This section determines the in plane components of the grain boundary segment,
       !based on the endpoints in the data file.  The first calculation of nl(i)
       !transforms line segments into the same reference frame as the Euler angles,
       !based on and assumption of default settings in TSL or HKL software.  The second
       !calculation determines the normal to the grain boundary trace.
       if (rotation.eq.1) then       !TSL
        nl(1)=y1-y2
        nl(2)=x1-x2
       endif
       if (rotation.eq.0) then       !HKL
        nl(1)=x2-x1
        nl(2)=-(y2-y1)
       endif
       if (rotation.eq.2) then       !D3D/natural
        nl(1)=x2-x1
        nl(2)=y2-y1
       endif
       templ=nl(2)                   !normal to the GB trace
       nl(2)=nl(1)
       nl(1)=-templ
       nl(3)=0.0
       length=sqrt(nl(1)*nl(1)+nl(2)*nl(2))

       !This section calculates phi for the GB normal, with respect to [100]
       mag_dir=((nl(1)**2)+(nl(2)**2))**0.5
       dot = (nl(1)*ref(1))+(nl(2)*ref(2))
       p = acos2(dot/(mag_dir*mag_ref))
       if (nl(2).lt.0.0) p=-p  !vectors pointing below the x-z plane assigned
                               ! negative angles (acos make them positive) so
                               ! the range of phi is 0 - 180 for a positive y
                               ! coorinate and 0 - -180 for a negative y coord
       if (p.lt.0.0) p=pi+p    ! For all directions with negative y (negative angle)
							   ! we consider the opposite direction
       !convert the Euler angles from the input file to a g matrix
       call EToG(e1,g_1)
       call EToG(e2,g_2)
       do i2=1,nsymm
        do i3=1,nsymm
		  call symop(O,i2,so_1)!call the i2-th 3x3 symmetry operator
		  call MToM (so_1,g_1,gg_1)!hit the first g  with the i2-th symmetry operator
		  call symop(O,i3,so_2)!call the i3-th 3x3 symmetry operator
		  call MToM(so_2,g_2,gg_2)!hit the second g  with the i3-th symmetry operator
		  call trans(gg_2,gg_2t)!transpose the second g matrix
		  call MToM(gg_1,gg_2t,dg)!misorientation between g_1 and g_2
		  call GToE(dg,mort) !get  Euler angles for this misorientation
          if (msym.eq.3) then
		   if (mort(1).lt.(pi*0.5+eps).and.mort(2).lt.(pi*0.5+eps).and.mort(3).lt.(pi*0.5+eps)) then
     		c1=int(float(CD)*2.0*mort(1)/pi)+1
			c2=int(float(CD)*cos(mort(2)))+1
			c3=int(float(CD)*2.0*mort(3)/pi)+1
			if (c1.eq.CD+1) c1=CD
			if (c2.eq.CD+1) c2=CD
			if (c3.eq.CD+1) c3=CD
            t=pi/2.0 !for the columnar approximation, theta is fixed
            call AnglesToV (t, p, n)
            call MToV (gg_1, n, nf)
            call VToAngles (nf, sa)
            c4=int(float(CD2)*cos(sa(1)))+1
            c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
            if (c4.eq.CD2+1) c4=CD2
            if (c4.eq.0) c4=1
            if (c5.eq.4*CD2+1) c5=4*CD2
            gbd(c1,c2,c3,c4,c5)=gbd(c1,c2,c3,c4,c5)+length
           endif
          else
           !determine the indices the GBCD matrix (c1,c2,c3)
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
           t=pi/2.0 !for the columnar approximation, theta is fixed
           call AnglesToV (t, p, n)
           call MToV (gg_1, n, nf)
           call VToAngles (nf, sa)
           c4=int(float(CD2)*cos(sa(1)))+1
           c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
           if (c4.eq.CD2+1) c4=CD2
           if (c4.eq.0) c4=1
           if (c5.eq.4*CD2+1) c5=4*CD2
           gbd(c1,c2,c3,c4,c5)=gbd(c1,c2,c3,c4,c5)+length
          endif
          !implement crystal exchange symmetry, and repeat the process above
		  call trans (gg_1, gg_1t)!transpose the fist g matrix
		  call MToM (gg_2, gg_1t, dg)!misorientation between g2 and g1
		  call GToE (dg, mort) !get  Euler angles for this misorientation
          if (msym.eq.3) then
		   if (mort(1).lt.(pi*0.5+eps).and.mort(2).lt.(pi*0.5+eps).and.mort(3).lt.(pi*0.5+eps)) then
     		c1=int(float(CD)*2.0*mort(1)/pi)+1
			c2=int(float(CD)*cos(mort(2)))+1
			c3=int(float(CD)*2.0*mort(3)/pi)+1
			if (c1.eq.CD+1) c1=CD
			if (c2.eq.CD+1) c2=CD
			if (c3.eq.CD+1) c3=CD
            t=pi/2.0
            call AnglesToV (t, p, n)
            call MToV (gg_2, n, nf)
            call VToAngles (nf, sa)
            c4=int(float(CD2)*cos(sa(1)))+1
            c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
            if (c4.eq.CD2+1) c4=CD2
            if (c4.eq.0) c4=1
            if (c5.eq.4*CD2+1) c5=4*CD2
            gbd(c1,c2,c3,c4,c5)=gbd(c1,c2,c3,c4,c5)+length
           endif
          else
           !determine the indices the GBCD matrix (c1,c2,c3)
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
           t=pi/2.0
           call AnglesToV (t, p, n)
           call MToV (gg_2, n, nf)
           call VToAngles (nf, sa)
           c4=int(float(CD2)*cos(sa(1)))+1
           c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
           if (c4.eq.CD2+1) c4=CD2
           if (c4.eq.0) c4=1
           if (c5.eq.4*CD2+1) c5=4*CD2
           gbd(c1,c2,c3,c4,c5)=gbd(c1,c2,c3,c4,c5)+length
          endif
         enddo ! end the i3=1,nsymm loop
        enddo ! end the i2=1,nsymm loop
        !notify user of the progress of the calculation
		if (mod(i1,5000).eq.0) then
	     write(6,"(I9,A)") i1,' line segments classified'
	    endif
 7050   continue
      enddo ! this closes the i1=header+1,nnline loop
      close(30)
      if (msym.eq.3) then
      avg=0.0
      do i1=1,CD
       do i2=1,CD
        do i3=1,CD
         do i4=1,CD2
          do i5=1,4*CD2
           avg=avg+gbd(i1,i2,i3,i4,i5)
          enddo
         enddo
        enddo
       enddo
      enddo
      avg=avg/float(4*CD*CD*CD*CD2*CD2)

      !write the MRD values of the distribution to the gbcd file
      do i1=1,CD
       do i2=1,CD
        do i3=1,CD
         do i4=1,CD2
          do i5=1,4*CD2
           gbd(i1,i2,i3,i4,i5)=gbd(i1,i2,i3,i4,i5)/avg
           write (46,*) gbd(i1,i2,i3,i4,i5)
          enddo
         enddo
        enddo
       enddo
      enddo

      else

      !determine the average line length in each bin
      avg=0.0
      do i1=1,CD*4
       do i2=1,CD*2
        do i3=1,CD*4
         do i4=1,CD2
          do i5=1,4*CD2
           avg=avg+gbd(i1,i2,i3,i4,i5)
          enddo
         enddo
        enddo
       enddo
      enddo
      avg=avg/float(4*4*4*2*CD*CD*CD*CD2*CD2)

      !write the MRD values of the distribution to the gbcd file
      do i1=1,CD*4
       do i2=1,CD*2
        do i3=1,CD*4
         do i4=1,CD2
          do i5=1,4*CD2
           gbd(i1,i2,i3,i4,i5)=gbd(i1,i2,i3,i4,i5)/avg
           write (46,*) gbd(i1,i2,i3,i4,i5)
          enddo
         enddo
        enddo
       enddo
      enddo
      endif
      close (46)






 9000 continue
      write(6,1) 'Program complete'

      end






	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
