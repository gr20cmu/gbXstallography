C*****************************************
       program condition triangles
C*****************************************
       implicit none
	   include 'common.fi'

       write(6,5) '================================================='
	   write(6,5) ' PROGRAM TO READ TRIANGLE DATA, REMOVE LINES '
       write(6,5) ' WHERE EULER ANGLES WITH ZERO VALUES ARE FOUND '
       write(6,5) ' '
	   write(6,5) ' rohrer@cmu.edu'
	   write(6,5) ' version 05/07/2022'
       write(6,5) '================================================='
	   
	   version = 'version 05/07/2022'
       !05/07/2022: This version just removes lines with zero Euler angles
       !            The parts used for filtering based on thresholds are
       !            commented out.

       !some constants that might be useful
       pi = 4.0*atan(1.0)

       !some formating statements used to format typed output
  5    format(A,A)
  6    format(A)
  15   format(A,I7,A)
  28   format(10f8.4)

       open(21, file='input.txt', status='old') ! Open the file with the program
                                                ! parameters.
	   read(21,*) keyword                       ! This reads the filename with the data
       !read(21,*) ThreshCurv
       !read(21,*) ThreshVel

       point = index(keyword,'.') - 1

       write(6,5)'I am going to look for data in the file labeled: ', keyword
       close(21)

       ! This part is to determines the total number of lines
       ! in the data file and assigns it to the variable, nnline
		nnline = 0 ! Initialize counter
		open(41, file=keyword, status='old') ! open the file
 60     continue
		read(41,*,end=61) ! read line, if end of file, goto 61
		nnline=nnline+1  ! Increment counter
		goto 60
 61     continue
        close(41)  ! we are finished, close the file
		write(6,15) 'There are  ',nnline,' lines in the data file.'
		
        ! Check each line for the string '#'. If found, we know this line
        ! is part of the header.  If not, we have reached the data section.
	    open(41, file=keyword, status='old') ! open the file
        header = 0 ! Initialize counter
        hash = 0   ! Initialize flag
	    do i=1,nnline  ! loop over each line in the file
	     read(41, "(a)") inline    ! reads line into a string
	     hash = index(inline, '#') ! if the character "#" in the string, hash>1.
	     if (hash.ne.0) then       ! check the flag
	      header=header+1          ! when "#" found, it is a header line
		  hash=0                   ! reset the flag
	     else
		  goto 80                  ! when hash is = 0, we have reached the end of the header
	     endif
	    enddo  ! closes i=1,nnline
 80    close(41)				   ! close the file
	   write(6,*)header,' lines in the header'

	   ! This opens the output file, where results are written

       open(31, file=keyword(1:point)//'_cnd.txt', status='unknown')
       write(31, 5)'# written by condition faces: ',version
       write(31, 5)'# initial file: ',keyword
       ! Now we read the data
       open(41, file=keyword, status='old') ! Open the input file
       nline = 0 ! Before reading the data, we have to read the neader
       do i=1,header
        read(41, "(a)") inline ! read a line of the header
        write(31, 6)inline  ! write that line into the output file
       enddo! closes i=1,header

       !write(31, 22)'# eliminated lines with curvature greater than',ThreshCurv
       !write(31, 22)'# eliminated lines with velocity greater than',ThreshVel

       blank = 0  ! Initialize counter for zero Euler angle lines

       ! Next is the start of the main loop, reading reading and
       ! analyzing one line of data at a time.
       do kk=header+1,nnline !start on the line after the header, continue to end.
        read(41,*) e1(1), e1(2), e1(3), e2(1), e2(2), e2(3),
     &  normal(1),normal(2),normal(3),area!, curv, vel
         area=abs(area)
         !curv=abs(curv)
         !vel=abs(vel)

        ! Outputs from some programs produce zero Euler angles when the
        ! when one or both orientiations are not known.  These lines are
        ! ignored because they will cause problems later.
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

        !if (curv.gt.ThreshCurv) then
        ! high_val=high_val+1
        ! goto 1000
        !else
        ! if (vel.gt.ThreshVel) then
        !  high_vel=high_vel+1
        !  goto 1000
        ! else
          write(31,28) e1(1), e1(2), e1(3), e2(1), e2(2), e2(3), normal(1),normal(2),normal(3), area!, curv, vel
         !endif
        !endif

         ! when you reach this point, the line has been written or ignored
 1000    continue ! so we go to the next line of the data file and repeat
       enddo    ! This closes kk=header+1,nnline, the loop over each line of data.


 9000    continue
         !write(6,*)'high_val=',high_val
         !write(6,*)'high_vel=',high_vel
         write(6,*)'zero Eulers:',blank
         write(6,*)'final number of lines =',nnline-header-blank
         close(31)
       end






	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
