C*****************************************
       program dist_graph
C*****************************************
       implicit none
	   include 'common.fi'

       write(6,4) '================================================='
	   write(6,4) ' USES GRAIN BOUNDARY TRIANGLES FROM 3D DATA TO '
	   write(6,4) ' TO COMPUTE DISTRIBUTIONS OF AREAS OR PROPERTIES '
	   write(6,4) ' AT SPECIFIC SETS OF POINTS IN THE FIVE-DIMENSIONAL '
	   write(6,4) ' (5D) SPACE.  THE PROGRAM WILL ALSO DETERMINE '
       write(6,4) ' TWO-DIMENSIONAL PROJECTIONS OF THE GRAIN BOUNDARY '
       write(6,4) ' PLANE DISTRIBUTION OR THE ONE-DIMENSIONAL  '
	   write(6,4) ' DISORIENTATION ANGLE DISTRIBUTION'
	   write(6,4) ' '
	   write(6,4) ' rohrer@cmu.edu'
	   write(6,4) ' version 05/27/2022'
       write(6,4) '================================================='
	   
	   version = 'version 05/27/2022'
       !Note: if there are more and 10 million representations
       !of triangles within the aperture, the program will crash.
       !if this is the problem, the fix is to increase the dimensions
       !of tris(10000000,13) in the common, and recompile with make.
       !This program synthesizes many of the functions previously available in
       !programs: disorientation_calculator, KDA_2d_graph, KDA_graph, and 
       !KDA_points.
       !The 1D disorientation distribution is set to work only for numbers
       !and areas of boundaries.
       !The 2D grain boundary plane distibution is set up to work only for cubic
       !and hexagonal crystals.
       !The error estimation in 2D cases is incorrect and should not be used.

       !useful constants defined here
       pi = 4.0*atan(1.0)

       !some formating statements used to format typed output
  4    format(A)
  5    format(A,A)
  6    format(A,I3,2X,A,I3)
  7    format(A,I3,A)
  8    format(A,I12,A)
  9    format(A,f4.1)
  10   format(f5.1,6x,F6.4,5x,F6.4)
  11   format(A,A,A,A,f8.3,1x,f8.3,1x,f8.3)
  12   format(A,A,A,A,f8.3,1x,f8.3,1x,f8.3)
  13   format(A,I4,A,I4)
  14   format(A,3I3,A,f5.2)
  15   format(A,I9,A,I9,A)
  20   format(A,A,A,A,f5.2,1x,f8.2,1x,f7.2,A)
  22   format(A,A,A,A,f6.3,1x,f8.3,1x,f8.3)
  23   format(f6.2,1x,f10.4,1x,I12,1x,f10.4)
  24   format(I5,f8.2,1x,f10.2,1x,I9,1x,I8,f10.2)

       !this section reads the user specified parameters in the file, input.txt
       open(21, file='input.txt', status='old') ! Open the file
	   read(21,*) keyword                       ! This reads the filename with the data
       read(21,*) mode                          ! Specifies the calculation mode
	   read(21,*) tempword, tempword            ! Reads a line containing headers
       read(21,*) msym, radian, dist_type       ! Specifies the symmetry, units for angles, 5D plot type
       read(21,*) tempword                      ! Read a header line
	   read(21,*) rho_m, rho_p, BallVolumeIndex ! Read the dimensions of the aperture
       read(21,*) ax(1), ax(2), ax(3), ang      ! This is the fixed misorientation for the 5D mode
       read(21,*) projection                    ! projection axis for output in 5D mode
       read(21,*) name                          ! added to the output filename
       read(21,*) keyword2                      ! read filename containing list of fixed points for mode 6
       read(21,*) errors                        ! determines whether or not error estimates are written
       close(21)                                ! Close the input file
       point = index(keyword,'.') - 1           ! Determines the number of characters in the filename, before the extension
       point2 = index(keyword2,'.') - 1
       !this updates the user about the parameters on the consol.
       write(6,5)'I am going to look for data in the file labeled: ', keyword
       write(6,6)'msym =',msym,'radian =',radian
       if (mode.eq.1) then
        write(6,4)'calculating the disorientation distribution'
       endif
       if (mode.eq.2) then
        write(6,4)'calculating the distribution as function of grain boundary plane orientation, ignoring disorientation'
       endif
       if (mode.eq.5) then
        write(6,4)'calculating the distribution as function of grain boundary plane orientation, at a fixed disorientation'
       endif
       if (mode.eq.6) then
        write(6,4)'calculating the value at a user defined list of points in the 5D space'
        write(6,5)'Finding the values at the points listed in: ',trim(keyword2)
       endif
       if (mode.ne.1.and.mode.ne.2.and.mode.ne.5.and.mode.ne.6) then
        write(6,7)'selected mode (',mode,') not valid - end program'
        goto 9000
       endif

       if (dist_type.eq.1) then
        write(6,4)'the input is being interpreted as grain boundary area data'
       endif
       if (dist_type.eq.2) then
        write(6,4)'the input is being interpreted as grain boundary energy data'
       endif
       if (dist_type.eq.3) then
        write(6,4)'the input is being interpreted as grain boundary curvature data'
       endif
       if (dist_type.eq.4) then
        write(6,4)'the input is being interpreted as grain boundary velocity data'
       endif
       if (dist_type.eq.5) then
        write(6,4)'the input is being interpreted as grain boundary dihedral angle data'
       endif
       if (dist_type.eq.6) then
        write(6,4)'the input is being interpreted as grain boundary velocity x curvature data'
       endif

       ! This part of the code determines the length of the file and the length of the header
       ! and has to be determined for all modes
       ! This part determines the total number of lines
       ! in the data file and assigns it to the variable, nnline
		nnline = 0 ! Initialize counter
        open(41, file=keyword, status='old')
 60     continue
		read(41,*,end=61) ! read line, if end of file, goto 61
		nnline=nnline+1  ! Increment counter
		goto 60
 61     continue
        close(41)  ! we are finished, close the file, report result.
		write(6,8) 'There are  ',nnline,' lines in the data file.'
		
        ! Check each line for the string '#'. If found, we know this line
        ! is part of the header.  If not, we have reached the data section.
        open(41, file=keyword, status='old') ! open the file
        header = 0 ! Initialize counter
        hash = 0   ! Initialize flag
	    do i1=1,nnline  ! loop over each line in the file
	     read(41, "(a)") inline    ! reads line into a string
	     hash = index(inline, '#') ! if the character "#" in the string, hash>1.
	     if (hash.ne.0) then       ! check the flag
	      header=header+1          ! when "#" found, it is a header line
		  hash=0                   ! reset the flag
	     else
		  goto 80                  ! when hash is = 0, we have reached the end of the header
	     endif
	    enddo  ! closes i1=1,nnline
 80    close(41)				   ! close the file, report the results
	   write(6,8)'There are ',header,' lines in the header'
       write(6,8) 'There are  ',nnline-header,' triangles in the data file.'

       call get_symop (msym, O, nsymm)
       !call get_symop (msym, O, nsymm) calls the subroutine that contains the symmetry operators.
       !The variable msym, read from the input file, specifies the symmetry (1 = tetragonal,
       !2 = hexagonal, 3 = cubic, 4 = trigonal, 5 = orthorhombic), O is a matrix containing the
       !operators, and nsymm is the number of symmetry operators.  The routine returns O and nsymm

       !The value of mode determines which part of the program is used.
       !After this point, each calculation mode proceeds independently
       if (mode.eq.1) goto 1000
       if (mode.eq.2) goto 2000
       if (mode.eq.5) goto 5000
       if (mode.eq.6) goto 6000

!======================================================================
!      This is the mode = 1 calculation of the disorientation angle distribution
!======================================================================

 1000  continue
       ! This defines the maximum angle for a disorientation distribution
       if (msym.eq.1) limit=99  !max eq 98.42,  D4, 422, 4/m mm
	   if (msym.eq.2) limit=94  !max eq 93.84,  D6, 622, 6/m mm	   
	   if (msym.eq.3) limit=63  !max eq 62.7994, O, 432, m-3m
	   if (msym.eq.4) limit=105 !max eq 104.5,  D3,  32, -3m
	   if (msym.eq.5) limit=120 !max eq 120.0,  D2, 222, mmm
       bin = 1.0 ! this is the bin width in degrees.  Perhaps make it a variable?
	   bin_num = int(limit/bin) !computes the number of bins in the distribution
       !repeats back to you the parameters you specified in the input.
       write(6,9)'degrees per bin =',bin
       !this just ensures the initial values of the distribution are zero
       do i1=1,2000
        d(i1,1)=0.0
        d(i1,2)=0.0
       enddo
       tris=0.0
       !open the file for the results
       open(30, file=keyword(1:point)//trim(name)//'_dist.txt',status='unknown')
       !open the file with the data
       open(41, file=keyword, status='old')
       write(30,5)'# written by dist_graph, version: ',version
       write(30,5)'# based on data in file: ',trim(keyword)
       write(30,4)'# the following lines are the header of the source file'
       do i1=1,header
        read(41, 4) inline ! read a line of the header
        write(30, 4)inline ! write that line into the output file
       enddo! closes i1=1,header
       blank = 0  ! Initialize counter for lines with zero Euler angles
       !this is the main loop over all lines of data
       do i1=header+1,nnline
        if (dist_type.eq.1) then !This case is for a GB plane distribution
         read(41,*) e1(1), e1(2), e1(3), e2(1), e2(2), e2(3), normal_lab(1), normal_lab(2), normal_lab(3), area_real
         area=area_real
        endif
        if (dist_type.eq.2) then !This case is for a GB energy distribution
         read(41,*) e1(1), e1(2), e1(3), e2(1), e2(2), e2(3), normal_lab(1), normal_lab(2), normal_lab(3), area_real, energy
         area=abs(energy)
        endif
        if (dist_type.eq.3) then !This case is for a GB curvature distribution
         read(41,*) e1(1), e1(2), e1(3), e2(1), e2(2), e2(3), normal_lab(1), normal_lab(2), normal_lab(3), area_real, energy, curv
         area=abs(curv)
        endif
        if (dist_type.eq.4) then !This case is for a GB velocity distribution
         read(41,*) e1(1), e1(2), e1(3), e2(1), e2(2), e2(3), normal_lab(1), normal_lab(2), normal_lab(3), area_real, energy, curv, vel
        area=abs(vel)
        endif
        if (dist_type.eq.5) then !This case is for the GB dihedral angle distribution
         read(41,*) e1(1), e1(2), e1(3), e2(1), e2(2), e2(3), normal_lab(1), normal_lab(2), normal_lab(3), area_real, energy, curv, vel, dihedral
         area=abs(dihedral)
        endif
        if (dist_type.eq.6) then !This case is for a GB velocity * curvature distribution
         read(41,*) e1(1), e1(2), e1(3), e2(1), e2(2), e2(3), normal_lab(1), normal_lab(2), normal_lab(3), area_real, energy, curv, vel, dihedral, velxcurve
         area=abs(velxcurve)
        endif

        !outputs from some programs produce zero Euler angles when
        !one or both orientiations are unknown.  These lines are
        !ignored because they will cause problems later.
		if (e1(1).eq.0.00.and.e1(2).eq.0.00.and.e1(3).eq.0.00) then
          blank=blank+1 ! increment counter
		  goto 1100 ! Sends it to the next line
        else
          continue
		endif
        if (e2(1).eq.0.00.and.e2(2).eq.0.00.and.e2(3).eq.0.00) then
          blank=blank+1 ! increment counter
          goto 1100 ! Sends it to the next line
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
        !if (disor2.lt.0.01) then !I don't understand it, but alot of GBs come through
        ! goto 1100               !D3D that have misorientations less than 1°.  Given
        !endif                    !the reconstruction, they can't be physical, so skip
        ! The next section bins the observation in the matix d.
        ! The first column of d accumulates the areas of the
        ! trianges and the second the number accumulates the
        ! numnber of triangles.
        disor2=disor2*(180.0/pi)     ! convert the value to degrees
        ind = floor(disor2/bin)+1    ! finds index of d (smallest integer > 0 for bin number)
        if (ind.eq.0) ind=1          ! make sure the index is non-zero
        if (ind.gt.bin_num) then     ! make sure it is not > the uper limit
         ind=bin_num
         write(6,*)'misorientation out of bounds'
        endif
        d(ind,1)=d(ind,1)+area         ! add the area to the first column of d
        d(ind,2)=d(ind,2)+1.0          ! increment the second column of d
        TotalArea=TotalArea+area ! keep track ot the total area of all triangles
        ! when you reach this point, the orientation pair has been classified
 1100   continue
        ! notify user of progress
        if (mod(i1,200000).eq.0) then
         write(6,15) 'classified ',i1, ' of ',nnline,' triangles'
        endif

       enddo !this closes the i1=header+1,nnline loop
       close(41)

       !here we compute and write the distribution
       ct=nnline-header-blank             ! ct is the lines of useful data
       ! The next four lines normalize the data
       do i1=1,bin_num                    ! loop over each bin
        d_norm(i1,1)=d(i1,1)/TotalArea ! divide area by the total area
        d_norm(i1,2)=d(i1,2)/float(ct)    ! divide number by the total number
       enddo
       !these statements write header information
       write(6,4)' '
       write(30,4)'# '
       write(6,4)'Disorientation distribution: rotation angle, area fraction, and number fraction'
       write(30,4)'# Disorientation distribution: rotation angle, area fraction, and number fraction'
       write(6,4)' '
       write(30,4)'# '
       write(6,4)'degrees    area       number'
       write(30,4)'  degrees  area       number'
       !these statements write each line of the distribution
       do i1=1,bin_num
		write(30,10)i1*bin, d_norm(i1,1),d_norm(i1,2)
        write(6,10)i1*bin, d_norm(i1,1),d_norm(i1,2)
       enddo ! closes the i1=0,bin_num loop
       close(30) ! closes the output file

       goto 9000

!================================================================================
!      This is the mode = 2 calculation of the 2D Grain Boundary Plane Distribution
!================================================================================
 2000  continue

       !This opens the output file, where results will be written
       open(32, file=keyword(1:point)//trim(name)//'_gmt1.gpf',status='unknown')
       !if data on errors was requested, a file is opened for this distribution
       if (errors.eq.1) then
        open(33, file=keyword(1:point)//trim(name)//'_errors_gmt1.gpf',status='unknown')
       else
        continue
       endif

       !This is the factor needed for the normalization
       BallVolume2D = 2.0*nsymm*(1.0 - cos(rho_p*(pi/180.0)))
       rho_p=rho_p*(pi/180.0)  !We have to convert the aperture values to radians

       !read in the sampling points
       if (msym.eq.3) then
        open(22, file='2d_grid_points_cubic.txt',status='old')
       else  !currently, hexagonal is the only possible choice other than cubic
        open(22, file='2d_grid_points_hex.txt',status='old')
       endif
       read(22,*)num_samplPts
       do i1=1,num_samplPts
        read(22,*)samplPts(1,i1),samplPts(2,i1),samplPts(3,i1)
       enddo
       close(22)

       do i1=1,num_samplPts
        dist(i1)=0.0
        err(i1)=0.0
        num(i1)=0.0
       enddo

       open(41, file=keyword, status='old') ! open the file and read the header
       do i1=1,header
        read(41, "(a)") inline ! read a line of the header
       enddo! ends i1=1,header loop

       write(6,4)'reading data'
       totalArea = 0.0
       do i1=header+1,nnline !read all data to the array tris
        if (dist_type.eq.1) then !This case is for a GB plane distribution
         read(41,*) tris(i1,1),tris(i1,2),tris(i1,3),tris(i1,4),tris(i1,5),tris(i1,6),tris(i1,7),tris(i1,8),tris(i1,9),tris(i1,10)
         area=tris(i1,10)
        endif
        if (dist_type.eq.2) then !This case is for a GB energy distribution
         read(41,*) tris(i1,1),tris(i1,2),tris(i1,3),tris(i1,4),tris(i1,5),tris(i1,6),tris(i1,7),tris(i1,8),tris(i1,9),tris(i1,10),tris(i1,11)
         area=abs(tris(i1,11))
        endif
        if (dist_type.eq.3) then !This case is for a GB curvature distribution
         read(41,*) tris(i1,1),tris(i1,2),tris(i1,3),tris(i1,4),tris(i1,5),tris(i1,6),tris(i1,7),tris(i1,8),tris(i1,9),tris(i1,10),tris(i1,11),tris(i1,12)
         area=abs(tris(i1,12))
        endif
        if (dist_type.eq.4) then !This case is for a GB velocity distribution
         read(41,*) tris(i1,1),tris(i1,2),tris(i1,3),tris(i1,4),tris(i1,5),tris(i1,6),tris(i1,7),tris(i1,8),tris(i1,9),tris(i1,10),tris(i1,11),tris(i1,12),tris(i1,13)
        area=abs(tris(i1,13))
        endif

        if (tris(i1,1).eq.0.00.and.tris(i1,2).eq.0.00.and.tris(i1,3).eq.0.00) then
         goto 2500 ! Sends it to the next line
        else
         continue
        endif
        if (tris(i1,4).eq.0.00.and.tris(i1,5).eq.0.00.and.tris(i1,6).eq.0.00) then
         goto 2500 ! Sends it to the next line
        else
         continue
        endif
        totalArea=totalArea+area

2500    continue
       enddo
       close(41)

       write(6,4)'computing distribution'
       do j1=1,num_samplPts  ! this is a loop over all grid points
        sum = 0.0 !zero the counters for each grid point
        c = 0.0
        numTris = 0.0
        do i1=header+1,nnline
         e1(1)=tris(i1,1)    ! here we assign the variables to their conventional names
         e1(2)=tris(i1,2)    ! to be consistent with borrowed code below.
         e1(3)=tris(i1,3)
         e2(1)=tris(i1,4)
         e2(2)=tris(i1,5)
         e2(3)=tris(i1,6)
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
         normal_lab(1)=tris(i1,7)
         normal_lab(2)=tris(i1,8)
         normal_lab(3)=tris(i1,9)
         if (dist_type.eq.1) then !This case is for a GB plane distribution
          area=tris(i1,10)
         endif
         if (dist_type.eq.2) then !This case is for a GB energy distribution
          area=abs(tris(i1,11))
         endif
         if (dist_type.eq.3) then !This case is for a GB curvature distribution
          area=abs(tris(i1,12))
         endif
         if (dist_type.eq.4) then !This case is for a GB velocity distribution
          area=abs(tris(i1,13))
         endif

        !The Euler angles are converted to 3x3 orientation matrices, g_o1 and g_o2.
         call EToG(e1,g_o1)
         call EToG(e2,g_o2)
         fixedNormal1(1)=samplPts(1,j1)
         fixedNormal1(2)=samplPts(2,j1)
         fixedNormal1(3)=samplPts(3,j1)


         do i_sy=1,nsymm             ! loop over all symmetry elements
          call MToV(g_o1, normal_lab, normal_grain1)
          call MToV(g_o2, normal_lab, normal_grain2)
          call symop(O,i_sy,so_1)     ! get first symmetry element
          call MToV(so_1, normal_grain1, sym_normal1)
          call MToV(so_1, normal_grain2, sym_normal2)
          sign=1.0
          chi1=acos2(sign*( (sym_normal1(1)*fixedNormal1(1)) + (sym_normal1(2)*fixedNormal1(2)) + (sym_normal1(3)*fixedNormal1(3)) ))
          chi2=acos2(sign*( (sym_normal2(1)*fixedNormal1(1)) + (sym_normal2(2)*fixedNormal1(2)) + (sym_normal2(3)*fixedNormal1(3)) ))
          if (chi1.lt.rho_p) then
           numTris=numTris+1.0
           y = area - c
           t = sum + y
           c = (t - sum) - y
           sum = t
          endif
          if (chi2.lt.rho_p) then
           numTris=numTris+1.0
           sum = sum + area
          endif

          sign=-1.0
          chi1=acos2(sign*( (sym_normal1(1)*fixedNormal1(1)) + (sym_normal1(2)*fixedNormal1(2)) + (sym_normal1(3)*fixedNormal1(3)) ))
          chi2=acos2(sign*( (sym_normal2(1)*fixedNormal1(1)) + (sym_normal2(2)*fixedNormal1(2)) + (sym_normal2(3)*fixedNormal1(3)) ))
          if (chi1.lt.rho_p) then
           numTris=numTris+1.0
           y = area - c
           t = sum + y
           c = (t - sum) - y
           sum = t
          endif
          if (chi2.lt.rho_p) then
           numTris=numTris+1.0
           sum = sum + area
          endif

         enddo  !This closes the i_sy=1,nsymm loop

        enddo    ! This closes i1=header+1,nnline loop

        dist(j1)=sum ! assigns the total area at this sample point
        num(j1)=numTris ! assigns number of energy values added at this sample point

        if (mod(j1,50).eq.0) then
         write(6,13) 'At grid point ',j1,' of ',num_samplPts
        endif

       enddo !j1 loop over sample points

       dist_max=0.0
       dist_min=10000000.0

       !normalize the distributions, and find the extreme values
       if (dist_type.eq.1) then !This part is for grain boundary plane distributions
        do i1=1, num_samplPts
         dist(i1)=dist(i1)/totalArea
         dist(i1) = dist(i1)/BallVolume2D
         if (dist(i1).gt.dist_max) then
          dist_max=dist(i1)
         endif
         if (dist(i1).lt.dist_min) then
          dist_min=dist(i1)
         endif
        enddo !ends the i1=1, num_samplPts loop
       else   !This part is for grain boundary property distributions
        do i1=1, num_samplPts
         dist(i1)=dist(i1)/num(i1)
         if (dist(i1).gt.dist_max) then
          dist_max=dist(i1)
         endif
         if (dist(i1).lt.dist_min) then
          dist_min=dist(i1)
         endif
        enddo !ends the i1=1, num_samplPts loop
       endif

       do i1=1,num_samplPts !here we write the output file
        vec(1)=samplPts(1,i1)
        vec(2)=samplPts(2,i1)
        vec(3)=samplPts(3,i1)
        call VToAngles(vec,angles)
        theta_plot=angles(1)
        phi_plot=angles(2)
        write (32,*) phi_plot*(180.0/3.14159), (90.0-(theta_plot*(180.0/3.14159))), dist(i1)
       enddo !closes the i1=1,num_samplPts loop
       close(32) ! closes the output file

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !this section of code is incorrect and needs to be fixed
       numDistinctGBs=100
       if (errors.eq.1) then  ! Write a file with errors, when requested (needs to be checked)
        do i1=1, num_samplPts
         err(i1)=sqrt((dist(i1)/totalArea/float(numDistinctGBs)))/BallVolume2D
        enddo
        do i1=1,num_samplPts
         vec(1)=samplPts(1,i1)
         vec(2)=samplPts(2,i1)
         vec(3)=samplPts(3,i1)
         call VToAngles(vec,angles)
         write (33,*) angles(2)*(180.0/pi), (90.0-(angles(1)*(180.0/pi))), dist(i1)
        enddo
        close(33) ! closes the output file
       else
        continue
       endif
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!

       interval = (dist_max-dist_min)/10.0 !Provide some guidance for setting the limits on the plot

       !this part outputs instructions for initiating Draw_stereograms.  One case for
       !hexagonal and another for cubic.
       if (msym.eq.2) then
        write(6,*)' '
        write(6,4)'to initiate the gmt script, enter ./Draw_stereograms_1 [Number of plots] [filename]_gmt_ IF rainbow [min] [max] [interval] stereo HEX'
        write(6,11)'For example: ./Draw_stereograms_1 1 ',keyword(1:point),trim(name),'_gmt IF rainbow ',dist_min,dist_max,interval,' stereo HEX'
       else
        write(6,*)' '
        write(6,4)'to initiate the gmt script, enter ./Draw_stereograms_1 [Number of plots] [filename]_gmt_ IF rainbow [min] [max] [interval]'
        write(6,12)'For example: ./Draw_stereograms_1 1 ',keyword(1:point),trim(name),'_gmt IF rainbow ',dist_min,dist_max,interval
       endif
       If (dist_type.eq.2) then
        write(6,4)'use ./Draw_stereograms_e'
       endif
       If (dist_type.eq.3) then
        write(6,4)'use ./Draw_stereograms_h'
       endif
       If (dist_type.eq.5) then
        write(6,4)'use ./Draw_stereograms_d'
       endif

       goto 9000

!=================================================================
!      This is the mode = 5 calculation of MRD value at all points
!      on a hemisphere at a fixed disorientation
!=================================================================
 5000  continue

       !this is the section where we set the viewpoint for the projection
       !for a [001] projection, there is no rotation of the data
	   if (projection.eq.1) then
		   ax_ref(1)=0.0
	       ax_ref(2)=0.0
	       ax_ref(3)=1.0
	       ang_ref= 0.0
	       write(6,4)'The [001] direction will be normal to the plane of the drawing.'
	   endif
       !for a [110] projection, rotate 90 deg around [1-10]
	   if (projection.eq.2) then
	       ax_ref(1)=0.7071
	       ax_ref(2)=-0.7071
	       ax_ref(3)=0.0
           ang_ref= pi*(90.0/180.0)
	       write(6,4)'The [110] direction will be normal to the plane of the drawing.'
	   endif
       !for a [111] projection, rotate 54.7 deg around [1-10]
	   if (projection.eq.3) then
	       ax_ref(1)=0.7071
	       ax_ref(2)=-0.7071
	       ax_ref(3)=0.0
           ang_ref= pi*(54.7/180.0)
	       write(6,4)'The [111] direction will be normal to the plane of the drawing.'
	   endif
       !for a [100] projection, rotate 90 deg around [0-10]
	   if (projection.eq.4) then
	       ax_ref(1)=0.0
	       ax_ref(2)=-1.0
	       ax_ref(3)=0.0
           ang_ref= pi*(90.0/180.0)
	       write(6,4)'The [100] direction will be normal to the plane of the drawing.'
	   endif
	   write(6,4) ' '
       !this is the matrix that will rotate the vectors on the plot.
	   call AAToG (ax_ref, ang_ref, g_ref)
	   call trans (g_ref, g_ref_t)

       !the following are normalization constants for grain boundary plane distribution data
       BallVolume(1) = 0.0000641361 !3° mis, 7° plane
       BallVolume(2) = 0.000139158  !5° mis, 5° plane
       BallVolume(3) = 0.000287439  !5° mis, 7° plane
       BallVolume(4) = 0.00038019   !5° mis, 8° plane
       BallVolume(5) = 0.000484151  !6° mis, 7° plane
       BallVolume(6) = 0.000747069  !7° mis, 7° plane
       BallVolume(7) = 0.00145491   !8° mis, 8° plane

       !This opens the output file, where results will be written
       open(34, file=keyword(1:point)//trim(name)//'_gmt1.dat',status='unknown')
       write(34,"(4(1x,f4.1))")  ax(1), ax(2), ax(3), ang
       write(6,14) 'misorientation = [',int(ax(1)), int(ax(2)), int(ax(3)),'] ', ang
       !if data on errors was requested, a file is opened for this distribution
       if (errors.eq.1) then
        open(35, file=keyword(1:point)//trim(name)//'_errors_gmt1.dat',status='unknown')
        write(35,"(4(1x,f4.1))")  ax(1), ax(2), ax(3), ang
       else
        continue
       endif

       !here we have to adjust the normalization factors to account for symmetry
       do i1=1,7
        BallVolume(i1)=BallVolume(i1)*(float(nsymm)/24.0)
       enddo

       rho_m=rho_m*(pi/180.0)  !We have to convert the aperture values to radians
       rho_p=rho_p*(pi/180.0)
       rho_pSq=rho_p*rho_p     !compute rho_p squared - a constant used later

       !we convert the axis angle description of the misorientation to a 3x3 orientation matrix, gFixed.
       !we start by normalizing the axis.
       call normalize(ax, axn)
       ax=axn
       ang=ang*(pi/180.0)           !convert the misorientation angle to radians
       write(6,"(A,4(1x,f6.3))") 'misorientation = ',ax(1), ax(2), ax(3), ang
       call AAToG (ax, ang, gFixed) !convert the axis and angle to an orientation matrix, g.
       call trans (gFixed, gFixedT) !get the transpose

       !read in the sampling points
       open(22, file='5d_grid_points.txt',status='old')
       read(22,*)num_samplPts
       do i1=1,num_samplPts
        read(22,*)samplPts(1,i1),samplPts(2,i1),samplPts(3,i1)
       enddo
       close(22)

       !intialize the arrays where we tabulate quantities
       do i1=1,cnt
        dist(i1)=0.0
        err(i1)=0.0
        num(i1)=0.0
       enddo

       !This is the main part of the code, where we read the triangles in one
       !at a time.  If they are within the aperture, we save their information.

       open(42, file=keyword, status='old') ! open the file and read the header
       do i1=1,header
        read(42, "(a)") inline ! read a line of the header
       enddo! ends i1=1,header loop
       totalArea = 0 !introduce area accumulator
       ct = 0 !introduce a counter for triangles
       num_unique =0 !introduce a counter for unique triangles
       NumTri = 0
       num=0.0
       dist=0.0

       do i1=header+1,nnline !start on the line after the header, continue to end.
        area_real = 0.0      !provide zeros for all values that might be read
        energy = 0.0
        curv = 0.0
        vel = 0.0
        velxcurve = 0.0
        dihedral = 0.0
        if (dist_type.eq.1) then !This case is for a GB plane distribution
         read(42,*) e1(1), e1(2), e1(3), e2(1), e2(2), e2(3), normal_lab(1), normal_lab(2), normal_lab(3), area_real
         area=area_real
        endif
        if (dist_type.eq.2) then !This case is for a GB energy distribution
         read(42,*) e1(1), e1(2), e1(3), e2(1), e2(2), e2(3), normal_lab(1), normal_lab(2), normal_lab(3), area_real, energy
         area=abs(energy)
        endif
        if (dist_type.eq.3) then !This case is for a GB curvature distribution
         read(42,*) e1(1), e1(2), e1(3), e2(1), e2(2), e2(3), normal_lab(1), normal_lab(2), normal_lab(3), area_real, energy, curv
         area=abs(curv)
        endif
        if (dist_type.eq.4) then !This case is for a GB velocity distribution
         read(42,*) e1(1), e1(2), e1(3), e2(1), e2(2), e2(3), normal_lab(1), normal_lab(2), normal_lab(3), area_real, energy, curv, vel
        area=abs(vel)
        endif
        if (dist_type.eq.5) then !This case is for the GB dihedral angle distribution
         read(42,*) e1(1), e1(2), e1(3), e2(1), e2(2), e2(3), normal_lab(1), normal_lab(2), normal_lab(3), area_real, energy, curv, vel, dihedral
         area=abs(dihedral)
        endif
        if (dist_type.eq.6) then !This case is for a GB velocity * curvature distribution
         read(42,*) e1(1), e1(2), e1(3), e2(1), e2(2), e2(3), normal_lab(1), normal_lab(2), normal_lab(3), area_real, energy, curv, vel, dihedral, velxcurve
         area=abs(velxcurve)
        endif

        totalArea=totalArea+area !add the area

        !If the Euler angles are provided in units of degrees, these
        !lines convert them to radians
	     if (radian.eq.0) then
          call DToRad3(e1, e1r)
          e1=e1r
          call DToRad3(e2, e2r)
          e2=e2r
         else
          continue
	     endif

        ! The Euler angles are convert to 3x3 orientation matrices, g_o1 and g_o2.
        call EToG(e1,g_o1)
        call EToG(e2,g_o2)

        ! Next, we compute delta_m and for all symm equivalents of this
        ! boundary and compare it to the aperture, rho_m.  If it is smaller than
        ! the aperture, we find the normal vectors and record them, with the
        ! area, in the array tris.  This is repeated for every triangle until
        ! the end of the file.
          in_aperture=0
		  do i_sy=1,nsymm             ! loop over all symmetry elements
           call symop(O,i_sy,so_1)    ! get first symmetry element for inner loop (so_1)
           call MToM (so_1,g_o1,gg_1) ! hit first orientation with a symmetry element
           call MToV(gg_1, normal_lab, normal_grain1)
           do j_sy=1,nsymm            ! loop over all symmetry elements
		    call symop(O,j_sy,so_2)   ! get the second symmetry element for outer loop (so_2)
		    call MToM(so_2,g_o2,gg_2) ! hit second orientation with a symmetry element
		    call trans(gg_2,g2sT)     ! take the transpose of the second orientation
		    call MToM(gg_1,g2sT,dg)   ! calculate the misorientation (dg)
            call trans(dg,dgT)
            call MToM(dg,gFixedT,diffFromFixed)
            tem=(diffFromFixed(1,1)+diffFromFixed(2,2)+diffFromFixed(3,3)-1)/2
            if (tem.gt.1.0) tem=1.0
            if (tem.lt.-1.0) tem=-1.0
            delta_m=acos2(tem)
            if (delta_m.lt.rho_m) then  ! check to see if it is within the aperture
             ct=ct+1 !increment ct everytime a triangle is within the aperture
             in_aperture=in_aperture+1
             call MToV(dgT, normal_grain1, normal_grain2)
             tris(ct,1)=normal_grain1(1)
             tris(ct,2)=normal_grain1(2)
             tris(ct,3)=normal_grain1(3)
             tris(ct,4)=-normal_grain2(1)
             tris(ct,5)=-normal_grain2(2)
             tris(ct,6)=-normal_grain2(3)
             tris(ct,7)=area
             tris(ct,8)=e1(1)
             tris(ct,9)=e1(2)
             tris(ct,10)=e1(3)
             tris(ct,11)=e2(1)
             tris(ct,12)=e2(2)
             tris(ct,13)=e2(3)
            else
             continue
            endif
            call MToM(dgT,gFixedT,diffFromFixed)
            tem=(diffFromFixed(1,1)+diffFromFixed(2,2)+diffFromFixed(3,3)-1)/2
            if (tem.gt.1.0) tem=1.0
            if (tem.lt.-1.0) tem=-1.0
            delta_m=acos2(tem)
            if (delta_m.lt.rho_m) then  ! check to see if it is within the aperture
             ct=ct+1 !increment ct everytime a triangle is within the aperture
             in_aperture=in_aperture+1
             call MToV(dgT, normal_grain1, normal_grain2)
             tris(ct,1)=-normal_grain2(1)
             tris(ct,2)=-normal_grain2(2)
             tris(ct,3)=-normal_grain2(3)
             tris(ct,4)=normal_grain1(1)
             tris(ct,5)=normal_grain1(2)
             tris(ct,6)=normal_grain1(3)
             tris(ct,7)=area
             tris(ct,8)=e1(1)
             tris(ct,9)=e1(2)
             tris(ct,10)=e1(3)
             tris(ct,11)=e2(1)
             tris(ct,12)=e2(2)
             tris(ct,13)=e2(3)
            else
             continue
            endif
		   enddo ! closes j_sy=1,nsymm
		  enddo ! closes i_sy=1,nsymm
          !this is a progress update
		  if (mod(i1,25000).eq.0) then
	       write(6,15) 'classified ',i1, ' of ',nnline,' triangles'
	      endif

          if (in_aperture.gt.0.) then
           num_unique=num_unique+1
           eus(num_unique,1)=e1(1)
           eus(num_unique,2)=e1(2)
           eus(num_unique,3)=e1(3)
           eus(num_unique,4)=e2(1)
           eus(num_unique,5)=e2(2)
           eus(num_unique,6)=e2(3)
          else
           continue
          endif
 5100     continue ! so we go to the next line of the data file and repeat
         enddo    ! This closes i1=header+1,nnline, the loop over each line of data.
         close(42)

         !Before the next step, we have to find out how many distinct
         !grain faces are in the data to get the error estimate correct
         NumberFaces=1
         Ueus(NumberFaces,1)=eus(1,1)
         Ueus(NumberFaces,2)=eus(1,2)
         Ueus(NumberFaces,3)=eus(1,3)
         Ueus(NumberFaces,4)=eus(1,4)
         Ueus(NumberFaces,5)=eus(1,5)
         Ueus(NumberFaces,6)=eus(1,6)
         do i3=2,num_unique
          do i4=1,NumberFaces
           diff=eus(i3,1)-Ueus(i4,1)+eus(i3,2)-Ueus(i4,2)+eus(i3,3)-Ueus(i4,3)+eus(i3,4)-Ueus(i4,4)+eus(i4,5)-Ueus(i4,5)+eus(i3,6)-Ueus(i4,6)
           diff=abs(diff)
           if (diff.lt.0.01) then
            goto 5200
           endif
          enddo
          NumberFaces=NumberFaces+1 !count it, it must be distinct
          Ueus(NumberFaces,1)=eus(i3,1)
          Ueus(NumberFaces,2)=eus(i3,2)
          Ueus(NumberFaces,3)=eus(i3,3)
          Ueus(NumberFaces,4)=eus(i3,4)
          Ueus(NumberFaces,5)=eus(i3,5)
          Ueus(NumberFaces,6)=eus(i3,6)
 5200     continue
         enddo

         write(6,8)'there are ',NumberFaces,' grain faces within the misorientation aperture'
         write(6,8)'there are ',ct,' triangles within the misorientation aperture'
         write(6,8)'there are ',num_unique,' unique triangles within the misorientation aperture'

         do i1=1,num_samplPts
          fixedNormal1_pre(1)=samplPts(1,i1)
          fixedNormal1_pre(2)=samplPts(2,i1)
          fixedNormal1_pre(3)=samplPts(3,i1)
          call MToV (g_ref_t, fixedNormal1_pre, fixedNormal1)!change projection
          call MToV(gFixedT, fixedNormal1, fixedNormal2)
          sum = 0.0
          numTri = 0
          do i2=1,ct
           normal_grain1(1)=tris(i2,1)
           normal_grain1(2)=tris(i2,2)
           normal_grain1(3)=tris(i2,3)
           normal_grain2(1)=tris(i2,4)
           normal_grain2(2)=tris(i2,5)
           normal_grain2(3)=tris(i2,6)
           area=tris(i2,7)
           sign=1.0
           chi1=acos2(sign*( (normal_grain1(1)*fixedNormal1(1)) + (normal_grain1(2)*fixedNormal1(2)) + (normal_grain1(3)*fixedNormal1(3)) ))
           chi2=acos2(-sign*( (normal_grain2(1)*fixedNormal2(1)) + (normal_grain2(2)*fixedNormal2(2)) + (normal_grain2(3)*fixedNormal2(3)) ))
           delta_pSq=0.5*((chi1*chi1)+(chi2*chi2))
           if (delta_pSq.lt.rho_pSq) then
            numTri=numTri+1
            sum = sum + area
            eus(numTri,1)=tris(i2,8) !record euler angles of accepted triangle
            eus(numTri,2)=tris(i2,9)
            eus(numTri,3)=tris(i2,10)
            eus(numTri,4)=tris(i2,11)
            eus(numTri,5)=tris(i2,12)
            eus(numTri,6)=tris(i2,13)
           else
            continue
           endif
           sign=-1.0
           chi1=acos(sign*( (normal_grain1(1)*fixedNormal1(1)) + (normal_grain1(2)*fixedNormal1(2)) + (normal_grain1(3)*fixedNormal1(3)) ))
           chi2=acos(-sign*( (normal_grain2(1)*fixedNormal2(1)) + (normal_grain2(2)*fixedNormal2(2)) + (normal_grain2(3)*fixedNormal2(3)) ))
           delta_pSq=0.5*((chi1*chi1)+(chi2*chi2))
           if (delta_pSq.lt.rho_pSq) then
            numTri=numTri+1
            sum = sum + area
            eus(numTri,1)=tris(i2,8) !record euler angles of accepted triangle
            eus(numTri,2)=tris(i2,9)
            eus(numTri,3)=tris(i2,10)
            eus(numTri,4)=tris(i2,11)
            eus(numTri,5)=tris(i2,12)
            eus(numTri,6)=tris(i2,13)
           else
            continue
           endif
          enddo !i2 loop over triangles
          dist(i1)=sum ! assigns the total area at this sample point
          num(i1)=float(numTri) ! assigns number of energy values added at this sample point
          if (mod(i1,500).eq.0) then
           write(6,"(A,I5,A,I5)") 'At grid point ',i1,' of ',num_samplPts
          endif

          !Before going to the next sample point, we have to find out how many distinct
          !grain faces were used at this point to get the error estimate
           NumberFaces=1
           Ueus(NumberFaces,1)=eus(1,1) !Take the first one as unique
           Ueus(NumberFaces,2)=eus(1,2)
           Ueus(NumberFaces,3)=eus(1,3)
           Ueus(NumberFaces,4)=eus(1,4)
           Ueus(NumberFaces,5)=eus(1,5)
           Ueus(NumberFaces,6)=eus(1,6)
           do i5=2,numTri
            do i6=1,NumberFaces
             diff=eus(i5,1)-Ueus(i6,1)+eus(i5,2)-Ueus(i6,2)+eus(i5,3)-Ueus(i6,3)+eus(i5,4)-Ueus(i6,4)+eus(i5,5)-Ueus(i6,5)+eus(i5,6)-Ueus(i6,6)
             diff=abs(diff)
            if (diff.lt.0.01) then
             goto 5300
            endif
           enddo !ends the i6=1,NumberFaces loop
           NumberFaces=NumberFaces+1
           Ueus(NumberFaces,1)=eus(i5,1)
           Ueus(NumberFaces,2)=eus(i5,2)
           Ueus(NumberFaces,3)=eus(i5,3)
           Ueus(NumberFaces,4)=eus(i5,4)
           Ueus(NumberFaces,5)=eus(i5,5)
           Ueus(NumberFaces,6)=eus(i5,6)
 5300      continue
          enddo ! ends the i5=2,numTri loop
          NumF(i)=NumberFaces

         enddo !i1 loop over sample points

         dist_max=0.0
         dist_min=1000.0

         !normalize the distributions, and find find the extreme values
         if (dist_type.eq.1) then !This part is for grain boundary plane distributions
          do i1=1, num_samplPts
           dist(i1)=dist(i1)/totalArea
           dist(i1) = dist(i1)/BallVolume(BallVolumeIndex)
            if (dist(i1).gt.dist_max) then
             dist_max=dist(i1)
            endif
            if (dist(i1).lt.dist_min) then
             dist_min=dist(i1)
            endif
          enddo ! ends the i1=1, num_samplPts loop
         else   !This part is for all other distributions
          do i1=1, num_samplPts
           dist(i1)=dist(i1)/num(i1)
            if (dist(i1).gt.dist_max) then
             dist_max=dist(i1)
            endif
            if (dist(i1).lt.dist_min) then
             dist_min=dist(i1)
            endif
          enddo !ends the i1=1, num_samplPts loop
         endif

         !write result at each sampling point
         do i1=1,num_samplPts
          vec(1)=samplPts(1,i1)
          vec(2)=samplPts(2,i1)
          vec(3)=samplPts(3,i1)
          call VToAngles(vec,angles)
          theta_plot=angles(1)
          phi_plot=angles(2)
          write (34,*) phi_plot*(180.0/3.14159), (90.0-(theta_plot*(180.0/3.14159))), dist(i1)
         enddo
         close(34) ! closes the output file

         if (errors.eq.1) then  ! Write a file with errors, when requested
          do i1=1, num_samplPts
           err(i1)=1/sqrt(NumF(i1)*BallVolume(BallVolumeIndex)*dist(i1))
          enddo
          do i1=1,num_samplPts
           vec(1)=samplPts(1,i1)
           vec(2)=samplPts(2,i1)
           vec(3)=samplPts(3,i1)
           call VToAngles(vec,angles)
           phi_plot=angles(1)
           theta_plot=angles(2)
           write (35,*) phi_plot*(180.0/3.14159), (90.0-(theta_plot*(180.0/3.14159))), err(i1)
          enddo
          close(35) ! closes the output file
         else
          continue
         endif

         interval = (dist_max-dist_min)/10.0 !Provide some guidance for setting the limits on the plot

         !this part outputs instructions for initiating Draw_stereograms.  One case for
         !hexagonal and another for cubic.
         if (msym.eq.2) then
          write(6,4)' '
          write(6,4)'to initiate the gmt script, enter ./Draw_stereograms [Number of plots] [filename]_gmt_ 5d rainbow [min] [max] [interval] stereo HEX'
          write(6,20)'For example: ./Draw_stereograms 1 ',keyword(1:point),trim(name),'_gmt 5d rainbow ',dist_min,dist_max,interval,' stereo HEX'
         else
          if (dist_type.eq.1) then
          write(6,4)' '
           write(6,4)'to initiate the gmt script, enter '
           write(6,22)'For example: ./Draw_stereograms 1 ',keyword(1:point),trim(name),'_gmt 5d rainbow ',dist_min,dist_max,interval
          endif
          if (dist_type.eq.2) then
          write(6,4)' '
           write(6,4)'to initiate the gmt script, enter ./Draw_stereograms [Number of plots] [filename]_gmt_ 5d rainbow [min] [max] [interval]'
           write(6,22)'For example: ./Draw_stereograms_e 1 ',keyword(1:point),trim(name),'_gmt 5d rainbow ',dist_min,dist_max,interval
          endif
          if (dist_type.eq.3) then
          write(6,4)' '
           write(6,4)'to initiate the gmt script, enter ./Draw_stereograms [Number of plots] [filename]_gmt_ 5d rainbow [min] [max] [interval]'
           write(6,22)'For example: ./Draw_stereograms_h 1 ',keyword(1:point),trim(name),'_gmt 5d rainbow ',dist_min,dist_max,interval
          endif
          if (dist_type.eq.4) then
          write(6,4)' '
           write(6,4)'to initiate the gmt script, enter ./Draw_stereograms_v [Number of plots] [filename]_gmt_ 5d rainbow [min] [max] [interval]'
           write(6,22)'For example: ./Draw_stereograms_v 1 ',keyword(1:point),trim(name),'_gmt 5d rainbow ',dist_min,dist_max,interval
          endif
          if (dist_type.eq.5) then
          write(6,4)' '
           write(6,4)'to initiate the gmt script, enter ./Draw_stereograms_d [Number of plots] [filename]_gmt_ 5d rainbow [min] [max] [interval]'
           write(6,22)'For example: ./Draw_stereograms_d 1 ',keyword(1:point),trim(name),'_gmt 5d rainbow ',dist_min,dist_max,interval
          endif
         endif


       goto 9000

!=================================================================
!      This is the mode = 6 calculation of the MRD value at an
!      arbitrary list of points read from a file
!=================================================================

 6000  continue

       !read in the sampling points
       open(43, file=trim(keyword2), status='old')
       read(43,*)num_misors
       if (num_misors.gt.200) then
        write(6,4)'The number of points exceeds the maximum, 200.  To proceed,'
        write(6,4)'you must adjust the dimensions of AxList, angs, and normals'
        write(6,4)'in common, adjust this conditional, and recompile with make.'
        write(6,4)'goodbye.'
        goto 9000
       endif
       do i1=1,num_misors
        read(43,*)AxList(i1,1),AxList(i1,2),AxList(i1,3),angs(i1),normals(i1,1),normals(i1,2),normals(i1,3)
       enddo
       close(43)

       !the following are normalization constants for grain boundary plane distribution data
       BallVolume(1) = 0.0000641361 !3° mis, 7° plane
       BallVolume(2) = 0.000139158  !5° mis, 5° plane
       BallVolume(3) = 0.000287439  !5° mis, 7° plane
       BallVolume(4) = 0.00038019   !5° mis, 8° plane
       BallVolume(5) = 0.000484151  !6° mis, 7° plane
       BallVolume(6) = 0.000747069  !7° mis, 7° plane
       BallVolume(7) = 0.00145491   !8° mis, 8° plane

       !here we have to adjust the normalization factors to account for symmetry
       do i1=1,7
        BallVolume(i1)=BallVolume(i1)*(float(nsymm)/24.0)
       enddo

       rho_m=rho_m*(pi/180.0)  !We have to convert the aperture values to radians
       rho_p=rho_p*(pi/180.0)
       rho_pSq=rho_p*rho_p     !compute rho_p squared - a constant used later

       do i1=1,num_misors
        dist(i1)=0.0
        err(i1)=0.0
        num(i1)=0.0
       enddo

       !This opens the output file, where results will be written
       open(44, file=keyword(1:point)//trim(name)//keyword2(1:point2)//'.txt',status='unknown')
       write(44,5)'# written by dist_graph: ',version
       write(44,5)'# based on data in: ',keyword
       write(44,5)'# using list of points labeled: ',keyword2
       if (dist_type.eq.1) then
        write(44,6)'# this is grain boundary relative area'
       endif
       if (dist_type.eq.2) then
        write(44,6)'# this is grain boundary relative energy'
       endif
       if (dist_type.eq.3) then
        write(44,6)'# this is grain boundary curvature'
       endif
       if (dist_type.eq.4) then
        write(44,6)'# this is grain boundary velocity'
       endif
       if (dist_type.eq.5) then
        write(44,6)'# this is the grain boundary dihedral angle'
       endif
       write(44,"(A,I4)")  '# number of points: ',num_misors
       write(6,"(A,I4)")  'number of points: ',num_misors
       write(44,6)'point   angle      value    number    faces    std error'

       do i1=1,num_misors  !this is the main loop over the points in the list
        write(6,4)'working on a new point ...'
        !we start by normalizing the axis and normal vectors.
        ax(1)=AxList(i1,1)
        ax(2)=AxList(i1,2)
        ax(3)=AxList(i1,3)
        call normalize(ax, axis)
        normal_lab(1)=normals(i1,1)
        normal_lab(2)=normals(i1,2)
        normal_lab(3)=normals(i1,3)
        call normalize(normal_lab, fixedNormal1)
        angs(i1)=angs(i1)*(pi/180.0)!convert the misorientation angle to radians
        !convert the axis and angle to 3x3 orientation matrix.
        call AAToG (axis, angs(i1), gFixed)
        call trans (gFixed, gFixedT) !get the transpose

        !In the next section of code, we read the triangles in one at a time.
        !If they are within the aperture, we save their information.
        open(41, file=keyword, status='old') ! open the file and read the header
        do i2=1,header
         read(41, "(a)") inline ! read a line of the header
        enddo! ends i2=1,header loop
        totalArea = 0 !introduce area accumulator
        ct = 0 !introduce a counter for triangles
        num_unique =0 !introduce a counter for unique triangles
        NumTri = 0
        num=0.0
        dist=0.0
        do i2=header+1,nnline !start on the line after the header, continue to end.
         area_real = 0.0      !provide zeros for all values that might be read
         energy = 0.0
         curv = 0.0
         vel = 0.0
         velxcurve = 0.0
         dihedral = 0.0
        if (dist_type.eq.1) then !This case is for a GB plane distribution
         read(41,*) e1(1), e1(2), e1(3), e2(1), e2(2), e2(3), normal_lab(1), normal_lab(2), normal_lab(3), area_real
         area=area_real
        endif
        if (dist_type.eq.2) then !This case is for a GB energy distribution
         read(41,*) e1(1), e1(2), e1(3), e2(1), e2(2), e2(3), normal_lab(1), normal_lab(2), normal_lab(3), area_real, energy
         area=abs(energy)
        endif
        if (dist_type.eq.3) then !This case is for a GB curvature distribution
         read(41,*) e1(1), e1(2), e1(3), e2(1), e2(2), e2(3), normal_lab(1), normal_lab(2), normal_lab(3), area_real, energy, curv
         area=abs(curv)
        endif
        if (dist_type.eq.4) then !This case is for a GB velocity distribution
         read(41,*) e1(1), e1(2), e1(3), e2(1), e2(2), e2(3), normal_lab(1), normal_lab(2), normal_lab(3), area_real, energy, curv, vel
        area=abs(vel)
        endif
        if (dist_type.eq.5) then !This case is for the GB dihedral angle distribution
         read(41,*) e1(1), e1(2), e1(3), e2(1), e2(2), e2(3), normal_lab(1), normal_lab(2), normal_lab(3), area_real, energy, curv, vel, dihedral
         area=abs(dihedral)
        endif
        if (dist_type.eq.6) then !This case is for a GB velocity * curvature distribution
         read(41,*) e1(1), e1(2), e1(3), e2(1), e2(2), e2(3), normal_lab(1), normal_lab(2), normal_lab(3), area_real, energy, curv, vel, dihedral, velxcurve
         area=abs(velxcurve)
        endif

         totalArea=totalArea+area !add the value to the accumulator
         !If the Euler angles are provided in units of degrees, convert to rad
	     if (radian.eq.0) then
          call DToRad3(e1, e1r)
          e1=e1r
          call DToRad3(e2, e2r)
          e2=e2r
         else
          continue
	     endif

         !The Euler angles are convert to 3x3 orientation matrices, g_o1 and g_o2.
         call EToG(e1,g_o1)
         call EToG(e2,g_o2)
         ! Next, we compute delta_m and for all symm equivalents of this
         ! boundary and compare it to the aperture, rho_m.  If it is smaller than
         ! the aperture, we find the normal vectors and record them, with the
         ! area, in the array tris.  This is repeated for every triangle until
         ! the end of the file.
          in_aperture=0
		  do i_sy=1,nsymm             ! loop over all symmetry elements
           call symop(O,i_sy,so_1)    ! get first symmetry element for inner loop (so_1)
           call MToM (so_1,g_o1,gg_1) ! hit first orientation with a symmetry element
           call MToV(gg_1, normal_lab, normal_grain1)
           do j_sy=1,nsymm            ! loop over all symmetry elements
		    call symop(O,j_sy,so_2)   ! get the second symmetry element for outer loop (so_2)
		    call MToM(so_2,g_o2,gg_2) ! hit second orientation with a symmetry element
		    call trans(gg_2,g2sT)     ! take the transpose of the second orientation
		    call MToM(gg_1,g2sT,dg)   ! calculate the misorientation (dg)
            call trans(dg,dgT)
            call MToM(dg,gFixedT,diffFromFixed)
            tem=(diffFromFixed(1,1)+diffFromFixed(2,2)+diffFromFixed(3,3)-1)/2
            if (tem.gt.1.0) tem=1.0
            if (tem.lt.-1.0) tem=-1.0
            delta_m=acos2(tem)
            if (delta_m.lt.rho_m) then  ! check to see if it is within the aperture
             ct=ct+1 !increment ct everytime a triangle is within the aperture
             in_aperture=in_aperture+1
             call MToV(dgT, normal_grain1, normal_grain2)
             tris(ct,1)=normal_grain1(1)
             tris(ct,2)=normal_grain1(2)
             tris(ct,3)=normal_grain1(3)
             tris(ct,4)=-normal_grain2(1)
             tris(ct,5)=-normal_grain2(2)
             tris(ct,6)=-normal_grain2(3)
             tris(ct,7)=area
             tris(ct,8)=e1(1)
             tris(ct,9)=e1(2)
             tris(ct,10)=e1(3)
             tris(ct,11)=e2(1)
             tris(ct,12)=e2(2)
             tris(ct,13)=e2(3)
            else
             continue
            endif
            call MToM(dgT,gFixedT,diffFromFixed)
            tem=(diffFromFixed(1,1)+diffFromFixed(2,2)+diffFromFixed(3,3)-1)/2
            if (tem.gt.1.0) tem=1.0
            if (tem.lt.-1.0) tem=-1.0
            delta_m=acos2(tem)
            if (delta_m.lt.rho_m) then  ! check to see if it is within the aperture
             ct=ct+1 !increment ct everytime a triangle is within the aperture
             in_aperture=in_aperture+1
             call MToV(dgT, normal_grain1, normal_grain2)
             tris(ct,1)=-normal_grain2(1)
             tris(ct,2)=-normal_grain2(2)
             tris(ct,3)=-normal_grain2(3)
             tris(ct,4)=normal_grain1(1)
             tris(ct,5)=normal_grain1(2)
             tris(ct,6)=normal_grain1(3)
             tris(ct,7)=area
             tris(ct,8)=e1(1)
             tris(ct,9)=e1(2)
             tris(ct,10)=e1(3)
             tris(ct,11)=e2(1)
             tris(ct,12)=e2(2)
             tris(ct,13)=e2(3)
            else
             continue
            endif
		   enddo ! closes j_sy=1,nsymm
		  enddo ! closes i_sy=1,nsymm

          !this is a progress update
		  if (mod(i2,250000).eq.0) then
	       write(6,15) 'classified ',i2, ' of ',nnline,' triangles'
	      endif

          continue ! so we go to the next line of the data file and repeat
         enddo    ! This closes i2=header+1,nnline, the loop over each line of data.

         write(6,15)'there are ',ct,' triangles within the misorientation aperture'
         close(41)



         call MToV(gFixedT, fixedNormal1, fixedNormal2)
         sum = 0.0
         numTri = 0
         do i2=1,ct
          normal_grain1(1)=tris(i2,1)
          normal_grain1(2)=tris(i2,2)
          normal_grain1(3)=tris(i2,3)
          normal_grain2(1)=tris(i2,4)
          normal_grain2(2)=tris(i2,5)
          normal_grain2(3)=tris(i2,6)
          area=tris(i2,7)
          sign=1.0
          chi1=acos2(sign*( (normal_grain1(1)*fixedNormal1(1)) + (normal_grain1(2)*fixedNormal1(2)) + (normal_grain1(3)*fixedNormal1(3)) ))
          chi2=acos2(-sign*( (normal_grain2(1)*fixedNormal2(1)) + (normal_grain2(2)*fixedNormal2(2)) + (normal_grain2(3)*fixedNormal2(3)) ))
          delta_pSq=0.5*((chi1*chi1)+(chi2*chi2))
          if (delta_pSq.lt.rho_pSq) then
           numTri=numTri+1
           sum = sum + area
           eus(numTri,1)=tris(i2,8) !record euler angles of accepted triangle
           eus(numTri,2)=tris(i2,9)
           eus(numTri,3)=tris(i2,10)
           eus(numTri,4)=tris(i2,11)
           eus(numTri,5)=tris(i2,12)
           eus(numTri,6)=tris(i2,13)
          else
           continue
          endif
          sign=-1.0
          chi1=acos(sign*( (normal_grain1(1)*fixedNormal1(1)) + (normal_grain1(2)*fixedNormal1(2)) + (normal_grain1(3)*fixedNormal1(3)) ))
          chi2=acos(-sign*( (normal_grain2(1)*fixedNormal2(1)) + (normal_grain2(2)*fixedNormal2(2)) + (normal_grain2(3)*fixedNormal2(3)) ))
          delta_pSq=0.5*((chi1*chi1)+(chi2*chi2))
          if (delta_pSq.lt.rho_pSq) then
           numTri=numTri+1
           sum = sum + area
           eus(numTri,1)=tris(i2,8) !record euler angles of accepted triangle
           eus(numTri,2)=tris(i2,9)
           eus(numTri,3)=tris(i2,10)
           eus(numTri,4)=tris(i2,11)
           eus(numTri,5)=tris(i2,12)
           eus(numTri,6)=tris(i2,13)
          else
           continue
          endif
         enddo !i2 loop over triangles
         dist(i1)=sum ! assigns the total area at this sample point
         num(i1)=float(numTri) ! assigns number of energy values added at this sample point
         
         !We have to find out how many distinct
         !grain faces were used at this point to get the error estimate
         NumberFaces=1
         Ueus(NumberFaces,1)=eus(1,1) !Take the first one as unique
         Ueus(NumberFaces,2)=eus(1,2)
         Ueus(NumberFaces,3)=eus(1,3)
         Ueus(NumberFaces,4)=eus(1,4)
         Ueus(NumberFaces,5)=eus(1,5)
         Ueus(NumberFaces,6)=eus(1,6)
         do i5=2,numTri
          do i6=1,NumberFaces
           diff=eus(i5,1)-Ueus(i6,1)+eus(i5,2)-Ueus(i6,2)+eus(i5,3)-Ueus(i6,3)+eus(i5,4)-Ueus(i6,4)+eus(i5,5)-Ueus(i6,5)+eus(i5,6)-Ueus(i6,6)
           diff=abs(diff)
           if (diff.lt.0.01) then
            goto 6300
           endif
          enddo !ends the i6=1,NumberFaces loop
          NumberFaces=NumberFaces+1
          Ueus(NumberFaces,1)=eus(i5,1)
          Ueus(NumberFaces,2)=eus(i5,2)
          Ueus(NumberFaces,3)=eus(i5,3)
          Ueus(NumberFaces,4)=eus(i5,4)
          Ueus(NumberFaces,5)=eus(i5,5)
          Ueus(NumberFaces,6)=eus(i5,6)
 6300     continue
         enddo ! ends the i5=2,numTri loop
         NumF(i1)=NumberFaces

         write(6,6) 'finished point ',i1,' of ',num_misors

         if (dist_type.eq.1) then
          dist(i1)=dist(i1)/totalArea
          dist(i1) = dist(i1)/BallVolume(BallVolumeIndex)
         else
          dist(i1)=dist(i1)/num(i1)
         endif
          err(i1)=1/sqrt(NumF(i1)*BallVolume(BallVolumeIndex)*dist(i1))
         write(6,24)i1,angs(i1)*(180.0/pi),dist(i1),int(num(i1)),NumberFaces,err(i1)
         write(44,24)i1,angs(i1)*(180.0/pi),dist(i1),int(num(i1)),NumberFaces,err(i1)
        enddo !this is the i1 loop
        close(44) ! closes the output file

       goto 9000



 9000    continue
         write(6,4)'program complete'

         end

