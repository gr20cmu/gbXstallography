C*****************************************
       program match
C*****************************************
       implicit none
	   include 'common.fi'

       write(6,*) '================================================='
	   write(6,*) ' PROGRAM FINDS TRIPLE POINTS ON ADJACENT '
	   write(6,*) ' LAYERS THAT ARE CONNECTED.  THE INPUT '
	   write(6,*) ' ARE TEXT FILES PRODUCED BY FIND_TJS. '
	   write(6,*) ' THESE FILES CONTAIN TRIPLETS OF LINE  '
	   write(6,*) ' SEGMENTS THAT MEET AT A POINT. '
	   write(6,*) ' THE OUTPUT OF THIS PROGRAM IS A LIST '
	   write(6,*) ' OF MATCHED TRIPLETS IN MATCH_NNN.TXT'
	   write(6,*) ' '
	   write(6,*) ' rohrer@cmu.edu'
	   write(6,*) ' version 06/17/2022'
       write(6,*) '================================================='
	   
	   version = 'version 06/17/2022'

       !some constants that might be useful
       pi = 4.0*atan(1.0)
       eps = 1.e-6

       !some formating statements used to format typed output
 1     format(A)
 2     format(A,A)
 3     format(A,I7,A,I3)
 4     format(7f9.3,f8.2,4f11.2)
 5     format(6f9.3,4f11.2)
 6     format(f5.2,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.0,1x,f5.0)
 7     format(A,I6,A)
 8     format(A,A,A)

       !here we read the program paramters from the input file
       open(21, file='input.txt', status='old')
       read(21,*) keyword
       read(21,*) slice
       read(21,*) first
	   read(21,*) last
	   read(21,*) skip
	   read(21,*) msym
       read(21,*) cols
       close(21)

       !some checks to make sure parameters are in acceptable ranges
       if (first.gt.998.or.last.gt.999) then
	    write(6,*) 'first and last exceed limits - edit input.txt - bye'
		goto 9000
	   endif
       if (cols.ne.12.and.cols.ne.10) then
	    write(6,*) 'only coded for 10 or 12 columns - edit input.txt - bye'
		goto 9000
	   endif

       !This starts the main loop over all the data files
       do file_n=first,last,skip
	    if (file_n+skip.gt.last) goto 8900
        !These statements resolve the names of the files with the data
	    msf=int(file_n/100) ! first file name
        nsf=int((file_n-(msf*100))/10)
		lsf=file_n-(100*msf)-(10*nsf)
		one=char(48+msf)
		two=char(48+nsf)
		three=char(48+lsf)
        fname1=trim(keyword)//'_'//one//two//three//'.txt'
        open (32, file=fname1,status='old')
	    msf=int((file_n+skip)/100) !second file name
        nsf=int(((file_n+skip)-(msf*100))/10)
		lsf=(file_n+skip)-(100*msf)-(10*nsf)
        one=char(48+msf)
        two=char(48+nsf)
        three=char(48+lsf)
        fname2=trim(keyword)//'_'//one//two//three//'.txt'
        open (33, file=fname2,status='old')

        !This part is to determines the total number of triple lines in the two data files
		nnline_1 = 0
 100    continue
		read(32,*,end=101)
		nnline_1=nnline_1+1
		goto 100
 101    continue
        close(32)

		nnline_2 = 0
 120    continue
		read(33,*,end=121)
		nnline_2=nnline_2+1
		goto 120
 121    continue
        close(33)

        !Now we figure out how long the header is by checking for '#' in each line.
        !there is an assumption here that each file has a header of the same length
        open (32, file=fname1,status='old')
        header = 0
        hash = 0
	    do i1=1,nnline_1
	     read(32, "(a)") inline
	     hash = index(inline, '#')
	     if (hash.ne.0) then
	      header=header+1
		  hash=0
	     else
		  goto 125
	     endif
	    enddo
 125    close(32)
        ! The number of triple lines is computed here
        n_1=int((nnline_1-header)/3)
		write(6,3) 'There are',n_1,' triple junctions in file ',file_n
		n_2=int((nnline_2-header)/3)
		write(6,3) 'There are',n_2,' triple junctions in file ',file_n+skip
      
        !this completes the preliminaries.  Now we reopen the files and
        !read the data.
        open (32, file=fname1,status='old')     ! Open first file
        open (33, file=fname2,status='old')     ! Open second file
        !These statements resolve the names of the files with the data
	    msf=int(file_n/100) ! first file name
        nsf=int((file_n-(msf*100))/10)
		lsf=file_n-(100*msf)-(10*nsf)
		one=char(48+msf)
		two=char(48+nsf)
		three=char(48+lsf)
        fname3=trim(keyword)//'_matches_'//one//two//three//'.txt'
        open (39, file=fname3,status='unknown') ! Open the output file
        write(39,8)'# written by match: ',trim(version),' previous header(s) follows.'
		counter = 0
        call get_symop (msym, O, nsymm)
        !first we copy the hear of one of the files into the output
        do i1=1,header
	     read(32, "(a)") inline
		 write(39,"(a)") inline
	    enddo

        write(39,1)'#  column labels for this file:'
        write(39,1)'#   phi1     PHI      phi2     phi1     PHI      phi2       x1          y1        x2          y2'
        !next we read the data
        if (cols.eq.12) then
		 do i1=1,n_1 ! read the first file
          read(32,4)tri_1(1,1,i1),tri_1(2,1,i1),tri_1(3,1,i1),tri_1(4,1,i1),
     &    tri_1(5,1,i1),tri_1(6,1,i1),tri_1(7,1,i1),tri_1(8,1,i1),
     &	  tri_1(9,1,i1),tri_1(10,1,i1),tri_1(11,1,i1),tri_1(12,1,i1)
		  read(32,4)tri_1(1,2,i1),tri_1(2,2,i1),tri_1(3,2,i1),tri_1(4,2,i1),
     &    tri_1(5,2,i1),tri_1(6,2,i1),tri_1(7,2,i1),tri_1(8,2,i1),
     &	  tri_1(9,2,i1),tri_1(10,2,i1),tri_1(11,2,i1),tri_1(12,2,i1)
		  read(32,4)tri_1(1,3,i1),tri_1(2,3,i1),tri_1(3,3,i1),tri_1(4,3,i1),
     &    tri_1(5,3,i1),tri_1(6,3,i1),tri_1(7,3,i1),tri_1(8,3,i1),
     &	  tri_1(9,3,i1),tri_1(10,3,i1),tri_1(11,3,i1),tri_1(12,3,i1)
	     enddo
	     close(32) ! close the first file
         do i1=1,header !read through the header of the second file
	      read(33, "(a)") inline
	     enddo
		 do i1=1,n_2 ! read the second file
          read(33,4)tri_2(1,1,i1),tri_2(2,1,i1),tri_2(3,1,i1),tri_2(4,1,i1),
     &    tri_2(5,1,i1),tri_2(6,1,i1),tri_2(7,1,i1),tri_2(8,1,i1),
     &	  tri_2(9,1,i1),tri_2(10,1,i1),tri_2(11,1,i1),tri_2(12,1,i1)
		  read(33,4)tri_2(1,2,i1),tri_2(2,2,i1),tri_2(3,2,i1),tri_2(4,2,i1),
     &    tri_2(5,2,i1),tri_2(6,2,i1),tri_2(7,2,i1),tri_2(8,2,i1),
     &	  tri_2(9,2,i1),tri_2(10,2,i1),tri_2(11,2,i1),tri_2(12,2,i1)
		  read(33,4)tri_2(1,3,i1),tri_2(2,3,i1),tri_2(3,3,i1),tri_2(4,3,i1),
     &    tri_2(5,3,i1),tri_2(6,3,i1),tri_2(7,3,i1),tri_2(8,3,i1),
     &	  tri_2(9,3,i1),tri_2(10,3,i1),tri_2(11,3,i1),tri_2(12,3,i1)
	     enddo
	     close(33) ! close the second file
        else
		 do i1=1,n_1 ! read the first file
          read(32,5)tri_1(1,1,i1),tri_1(2,1,i1),tri_1(3,1,i1),tri_1(4,1,i1),
     &    tri_1(5,1,i1),tri_1(6,1,i1),
     &	  tri_1(9,1,i1),tri_1(10,1,i1),tri_1(11,1,i1),tri_1(12,1,i1)
		  read(32,5)tri_1(1,2,i1),tri_1(2,2,i1),tri_1(3,2,i1),tri_1(4,2,i1),
     &    tri_1(5,2,i1),tri_1(6,2,i1),
     &	  tri_1(9,2,i1),tri_1(10,2,i1),tri_1(11,2,i1),tri_1(12,2,i1)
		  read(32,5)tri_1(1,3,i1),tri_1(2,3,i1),tri_1(3,3,i1),tri_1(4,3,i1),
     &    tri_1(5,3,i1),tri_1(6,3,i1),
     &	  tri_1(9,3,i1),tri_1(10,3,i1),tri_1(11,3,i1),tri_1(12,3,i1)
	     enddo
	     close(32) ! close the first file
         do i1=1,header !read through the header of the second file
	      read(33, "(a)") inline
	     enddo
		 do i1=1,n_2 ! read the second file
          read(33,5)tri_2(1,1,i1),tri_2(2,1,i1),tri_2(3,1,i1),tri_2(4,1,i1),
     &    tri_2(5,1,i1),tri_2(6,1,i1),
     &	  tri_2(9,1,i1),tri_2(10,1,i1),tri_2(11,1,i1),tri_2(12,1,i1)
		  read(33,5)tri_2(1,2,i1),tri_2(2,2,i1),tri_2(3,2,i1),tri_2(4,2,i1),
     &    tri_2(5,2,i1),tri_2(6,2,i1),
     &	  tri_2(9,2,i1),tri_2(10,2,i1),tri_2(11,2,i1),tri_2(12,2,i1)
		  read(33,5)tri_2(1,3,i1),tri_2(2,3,i1),tri_2(3,3,i1),tri_2(4,3,i1),
     &    tri_2(5,3,i1),tri_2(6,3,i1),
     &	  tri_2(9,3,i1),tri_2(10,3,i1),tri_2(11,3,i1),tri_2(12,3,i1)
	     enddo
	     close(33) ! close the second file
        endif
        !zero the values in the triple line matrix
        do i1=1,15000
         do i2=1,6
          tl(i2,i1)=0.0
         enddo
        enddo

        !this is the first main loop, where triple points on layer 1
        !are matched with points on the second layer.
        !We start by finding the coordinate of the triple junction.  We know that
        !one of the endpoints of the first segment must match one of the end points
        !of the second segment.  The match is the tiple point coordinate.
		counter=0
        do i1=1,n_1
		 if(tri_1(9,1,i1).eq.tri_1(9,2,i1).and.tri_1(10,1,i1)
     &	  .eq.tri_1(10,2,i1)) then
	       tl(1,i1)=tri_1(9,1,i1)
		   tl(2,i1)=tri_1(10,1,i1)
		   goto 135
		 endif
         if(tri_1(9,1,i1).eq.tri_1(11,2,i1).and.tri_1(10,1,i1)
     &	  .eq.tri_1(12,2,i1)) then
	      tl(1,i1)=tri_1(9,1,i1)
		  tl(2,i1)=tri_1(10,1,i1)
		  goto 135
		 endif
         !if the first coodinate doesn't match either of those from the second segment,
         !then the second coordinate must be the triple point.
		 tl(1,i1)=tri_1(11,1,i1)
		 tl(2,i1)=tri_1(12,1,i1)
 135     continue
         !In this loop, we identify the five triple junctions in the second layer
         !that are closest to the one in first layer. This is done by locating
         !each junction coordinate in the second layer (as above) and then
         !calculating the distance between it and the one on the first layer.  The
         !coordinates with the shortest distances are saved in the matrix clos.
         !Initialize clos
		 do i2=1,5
		  clos(i2,1)=200.0
		  clos(i2,2)=0.0
		 enddo
		 do i3=1,n_2
 		 if(tri_2(9,1,i3).eq.tri_2(9,2,i3).and.tri_2(10,1,i3)
     &	  .eq.tri_2(10,2,i3)) then
          tl(3,i3)=tri_2(9,1,i3)
		  tl(4,i3)=tri_2(10,1,i3)
		  goto 140
		 endif
         if(tri_2(9,1,i3).eq.tri_2(11,2,i3).and.tri_2(10,1,i3)
     &	  .eq.tri_2(12,2,i3)) then
	      tl(3,i3)=tri_2(9,1,i3)
		  tl(4,i3)=tri_2(10,1,i3)
		  goto 140
		 endif
		 tl(3,i3)=tri_2(11,1,i3)
		 tl(4,i3)=tri_2(12,1,i3)
 140    continue
        !calculate the distance between the TJ on layer 1 and the one on layer 2
		dx=tl(3,i3)-tl(1,i1)
        dy=tl(4,i3)-tl(2,i1)
		length=((dx*dx)+(dy*dy))**0.5
        !The following cascade remembers the TJ if it is close
		if (length.lt.clos(1,1)) then
		 clos(1,1)=length
		 clos(1,2)=float(i3)
		 clos(1,3)=tl(3,i3)
		 clos(1,4)=tl(4,i3)
		 goto 150
		endif
 		if (length.lt.clos(2,1)) then
		 clos(2,1)=length
		 clos(2,2)=float(i3)
		 clos(2,3)=tl(3,i3)
		 clos(2,4)=tl(4,i3)
		 goto 150
		endif
		if (length.lt.clos(3,1)) then
		 clos(3,1)=length
		 clos(3,2)=float(i3)
		 clos(3,3)=tl(3,i3)
		 clos(3,4)=tl(4,i3)
		 goto 150
		endif
		if (length.lt.clos(4,1)) then
		 clos(4,1)=length
		 clos(4,2)=float(i3)
		 clos(4,3)=tl(3,i3)
		 clos(4,4)=tl(4,i3)
		 goto 150
		endif
		if (length.lt.clos(5,1)) then
		 clos(5,1)=length
		 clos(5,2)=float(i3)
		 clos(5,3)=tl(3,i3)
		 clos(5,4)=tl(4,i3)
		 goto 150
		endif
 150	continue
       enddo  !  this is for the i3 loop, that finds the five closest

        !if (tri_1(1,1,i1).eq.4.229.and.tri_1(2,1,i1).eq.0.224.and.tri_1(3,1,i1).eq.1.061.and.tri_1(4,1,i1).eq.5.518) then
         !do i3=1,5
         !write(6,*)clos(i3,1),clos(i3,2),clos(i3,3),clos(i3,4)

         !enddo
         !write(6,*)' '
        !endif



       !Next, we have to decide if any of these close ones correspond to the
       !the triple junction on layer 1.
       do i4=1,5 !was k
        indx=int(clos(i4,2))
        !we start by finding the three unique euler angles for each junction
        !We know the first two have to be different
        e1_1(1)=tri_1(1,1,i1)
        e1_1(2)=tri_1(2,1,i1)
        e1_1(3)=tri_1(3,1,i1)
        e2_1(1)=tri_1(4,1,i1)
        e2_1(2)=tri_1(5,1,i1)
        e2_1(3)=tri_1(6,1,i1)
        !If the third set has even one different compentent, then this is the one we want.
        if(tri_1(1,2,i1).ne.e1_1(1).and.tri_1(1,2,i1).ne.e2_1(1)) then
         e3_1(1)=tri_1(1,2,i1)
         e3_1(2)=tri_1(2,2,i1)
         e3_1(3)=tri_1(3,2,i1)
         goto 160
        endif
        if(tri_1(2,2,i1).ne.e1_1(2).and.tri_1(2,2,i1).ne.e2_1(2)) then
         e3_1(1)=tri_1(1,2,i1)
         e3_1(2)=tri_1(2,2,i1)
         e3_1(3)=tri_1(3,2,i1)
         goto 160
        endif
        if(tri_1(3,2,i1).ne.e1_1(3).and.tri_1(3,2,i1).ne.e2_1(3)) then
         e3_1(1)=tri_1(1,2,i1)
         e3_1(2)=tri_1(2,2,i1)
         e3_1(3)=tri_1(3,2,i1)
         goto 160
        endif
        !If the components of the third all matched components in the first two, then
        !the fourth one must be the one we want.
        e3_1(1)=tri_1(4,2,i1)
        e3_1(2)=tri_1(5,2,i1)
        e3_1(3)=tri_1(6,2,i1)
 160    continue
        !now we repeat everything that was done above, but for the second layer
        e1_2(1)=tri_2(1,1,indx)
        e1_2(2)=tri_2(2,1,indx)
        e1_2(3)=tri_2(3,1,indx)
        e2_2(1)=tri_2(4,1,indx)
        e2_2(2)=tri_2(5,1,indx)
        e2_2(3)=tri_2(6,1,indx)
        if(tri_2(1,2,indx).ne.e1_2(1).and.tri_2(1,2,indx).ne.e2_2(1)) then
         e3_2(1)=tri_2(1,2,indx)
         e3_2(2)=tri_2(2,2,indx)
         e3_2(3)=tri_2(3,2,indx)
         goto 170
        endif
        if(tri_2(2,2,indx).ne.e1_2(2).and.tri_2(2,2,indx).ne.e2_2(2)) then
         e3_2(1)=tri_2(1,2,indx)
         e3_2(2)=tri_2(2,2,indx)
         e3_2(3)=tri_2(3,2,indx)
         goto 170
        endif
        if(tri_2(3,2,indx).ne.e1_2(3).and.tri_2(3,2,indx).ne.e2_2(3)) then
         e3_2(1)=tri_2(1,2,indx)
         e3_2(2)=tri_2(2,2,indx)
         e3_2(3)=tri_2(3,2,indx)
         goto 170
        endif
        e3_2(1)=tri_2(4,2,indx)
        e3_2(2)=tri_2(5,2,indx)
        e3_2(3)=tri_2(6,2,indx)
170     continue
        !next, we use a disorientation threshold to determine if the
        !grains on the first layer are the same as those on the second layer
        min_disor(1) = pi ! assign a large starting disorientation
        min_disor(2) = pi
        min_disor(3) = pi
        !first, we find the minimum disorientation between the first grain on
        !layer one and the three grains on layer two
        call DisG(O, nsymm, e1_1, e1_2, dg, ax, angle) !first grain
        if (angle.lt.min_disor(1)) then
         min_disor(1)=angle
        endif
        call DisG(O, nsymm, e1_1, e2_2, dg, ax, angle) !second grain
        if (angle.lt.min_disor(1)) then
         min_disor(1)=angle
        endif
        call DisG(O, nsymm, e1_1, e3_2, dg, ax, angle) !third grain
        if (angle.lt.min_disor(1)) then
         min_disor(1)=angle
        endif
        !next, we find the minimum disorientation between the second grain on
        !layer one and the three grains on layer two
        call DisG(O, nsymm, e2_1, e1_2, dg, ax, angle) !first grain
        if (angle.lt.min_disor(2)) then
         min_disor(2)=angle
        endif
        call DisG(O, nsymm, e2_1, e2_2, dg, ax, angle) !second grain
        if (angle.lt.min_disor(2)) then
         min_disor(2)=angle
        endif
        call DisG(O, nsymm, e2_1, e3_2, dg, ax, angle) !third grain
        if (angle.lt.min_disor(2)) then
         min_disor(2)=angle
        endif
        !finally, we find the minimum disorientation between the third grain on
        !layer one and the three grains on layer two
        call DisG(O, nsymm, e3_1, e1_2, dg, ax, angle) !first grain
        if (angle.lt.min_disor(3)) then
         min_disor(3)=angle
        endif
        call DisG(O, nsymm, e3_1, e2_2, dg, ax, angle) !second grain
        if (angle.lt.min_disor(3)) then
         min_disor(3)=angle
        endif
        call DisG(O, nsymm, e3_1, e3_2, dg, ax, angle) !third grain
        if (angle.lt.min_disor(3)) then
         min_disor(3)=angle
        endif

        !if (tri_1(1,1,i1).eq.4.229.and.tri_1(2,1,i1).eq.0.224.and.tri_1(3,1,i1).eq.1.061.and.tri_1(4,1,i1).eq.5.518) then
         !do i3=1,5
         !write(6,*)min_disor(1)*(180./pi),min_disor(2)*(180./pi),min_disor(3)*(180./pi)

         !enddo
         !write(6,*)' '
        !endif


        !Here we implement the criterion, if all three minimum disorientations
        !are less than 5°, the we assume that it is a match
        if (min_disor(1).lt.0.0873.and.min_disor(2).lt.0.0873
     &  .and.min_disor(3).lt.0.0873) then
         counter=counter+1
         tl(3,i1)=clos(i4,3)
         tl(4,i1)=clos(i4,4)
         tl(5,i1)=float(i1)
         tl(6,i1)=float(indx)
         if (cols.eq.12) then
		  write(39,4)tri_1(1,1,i1),tri_1(2,1,i1),tri_1(3,1,i1),tri_1(4,1,i1),
     &    tri_1(5,1,i1),tri_1(6,1,i1),tri_1(7,1,i1),tri_1(8,1,i1),
     &	  tri_1(9,1,i1),tri_1(10,1,i1),tri_1(11,1,i1),tri_1(12,1,i1)
		  write(39,4)tri_1(1,2,i1),tri_1(2,2,i1),tri_1(3,2,i1),tri_1(4,2,i1),
     &    tri_1(5,2,i1),tri_1(6,2,i1),tri_1(7,2,i1),tri_1(8,2,i1),
     &	  tri_1(9,2,i1),tri_1(10,2,i1),tri_1(11,2,i1),tri_1(12,2,i1)
		  write(39,4)tri_1(1,3,i1),tri_1(2,3,i1),tri_1(3,3,i1),tri_1(4,3,i1),
     &    tri_1(5,3,i1),tri_1(6,3,i1),tri_1(7,3,i1),tri_1(8,3,i1),
     &	  tri_1(9,3,i1),tri_1(10,3,i1),tri_1(11,3,i1),tri_1(12,3,i1)
		  write(39,4)tri_2(1,1,indx),tri_2(2,1,indx),tri_2(3,1,indx),tri_2(4,1,indx),
     &    tri_2(5,1,indx),tri_2(6,1,indx),tri_2(7,1,indx),tri_2(8,1,indx),
     &	  tri_2(9,1,indx),tri_2(10,1,indx),tri_2(11,1,indx),tri_2(12,1,indx)
		  write(39,4)tri_2(1,2,indx),tri_2(2,2,indx),tri_2(3,2,indx),tri_2(4,2,indx),
     &    tri_2(5,2,indx),tri_2(6,2,indx),tri_2(7,2,indx),tri_2(8,2,indx),
     &	  tri_2(9,2,indx),tri_2(10,2,indx),tri_2(11,2,indx),tri_2(12,2,indx)
		  write(39,4)tri_2(1,3,indx),tri_2(2,3,indx),tri_2(3,3,indx),tri_2(4,3,indx),
     &    tri_2(5,3,indx),tri_2(6,3,indx),tri_2(7,3,indx),tri_2(8,3,indx),
     &	  tri_2(9,3,indx),tri_2(10,3,indx),tri_2(11,3,indx),tri_2(12,3,indx)
         else
		  write(39,5)tri_1(1,1,i1),tri_1(2,1,i1),tri_1(3,1,i1),tri_1(4,1,i1),
     &    tri_1(5,1,i1),tri_1(6,1,i1),
     &	  tri_1(9,1,i1),tri_1(10,1,i1),tri_1(11,1,i1),tri_1(12,1,i1)
		  write(39,5)tri_1(1,2,i1),tri_1(2,2,i1),tri_1(3,2,i1),tri_1(4,2,i1),
     &    tri_1(5,2,i1),tri_1(6,2,i1),
     &	  tri_1(9,2,i1),tri_1(10,2,i1),tri_1(11,2,i1),tri_1(12,2,i1)
		  write(39,5)tri_1(1,3,i1),tri_1(2,3,i1),tri_1(3,3,i1),tri_1(4,3,i1),
     &    tri_1(5,3,i1),tri_1(6,3,i1),
     &	  tri_1(9,3,i1),tri_1(10,3,i1),tri_1(11,3,i1),tri_1(12,3,i1)
		  write(39,5)tri_2(1,1,indx),tri_2(2,1,indx),tri_2(3,1,indx),tri_2(4,1,indx),
     &    tri_2(5,1,indx),tri_2(6,1,indx),
     &	  tri_2(9,1,indx),tri_2(10,1,indx),tri_2(11,1,indx),tri_2(12,1,indx)
		  write(39,5)tri_2(1,2,indx),tri_2(2,2,indx),tri_2(3,2,indx),tri_2(4,2,indx),
     &    tri_2(5,2,indx),tri_2(6,2,indx),
     &	  tri_2(9,2,indx),tri_2(10,2,indx),tri_2(11,2,indx),tri_2(12,2,indx)
		  write(39,5)tri_2(1,3,indx),tri_2(2,3,indx),tri_2(3,3,indx),tri_2(4,3,indx),
     &    tri_2(5,3,indx),tri_2(6,3,indx),
     &	  tri_2(9,3,indx),tri_2(10,3,indx),tri_2(11,3,indx),tri_2(12,3,indx)
         endif
         goto 1000
		endif
       enddo       ! this is the i4 loop, over the five closest
 1000  continue
       enddo  ! ends the i1 loop
       write(6,7)'I found and wrote ',counter,' matched triple junctions'
 8900  close(39)

      enddo  ! This is the end of the file_n looop

 9000 continue
      write(6,1) 'Program complete'

      end






	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
