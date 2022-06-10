C*****************************************
       program find triple junctions
C*****************************************
       implicit none
	   include 'common.fi'

       write(6,*) '=============================================='
	   write(6,*) ' PROGRAM SEARCHES LISTS OF LINE SEGMENTS'
	   write(6,*) ' TO FIND TRIPLETS THAT MEET AT A POINT '
	   write(6,*) ' AND THEN WRITES AN ORDERED LIST TO A FILE '
	   write(6,*)
	   write(6,*) ' rohrer@cmu.edu'
	   write(6,*) ' version 04/12/22'
       write(6,*) '=============================================='
	     
       !gfortran find_tjs.f -o find_tjs
       !compile with gfortran -O3 -c -ffixed-line-length-none find_tjs.f -o find_tjs
	   ! modified 11/13/13 to work with more than 100 files
       ! modified 11/30/16 to improve filenaming infrastructure and deal with variable column numbers.
       ! modified 03/04/21 added a check for segments with zero length.
       ! modified 03/07/22: used program template so there is a make file
       ! added the correct column labels
       ! modified 03/09/22: added an option for 10 column output
       ! modified 04/12/22: updated common, sub2, symmetry with corrections, allow any filenames

       version='version 04/12/22'

       !some constants that might be useful
       pi = 4.0*atan(1.0)
       eps = 1.e-6

       !some formating statements used to format typed output
 1     format(A)
 2     format(A,A)
 3     format(A,I7,A,I3,A)
 4     format(7f9.3,f8.2,4f11.2)
 5     format(A,I5,A)
 6     format(6f9.3,4f11.2)

       !Open the input file, read the parameters, close the input file
       open(21, file='input.txt', status='old')
       read(21,*) keyword
       read(21,*) first
	   read(21,*) last
       read(21,*) col_in
       read(21,*) col_out
       close(21)

       !It will only work if the number of columns are 10, 12, 14, or 21
       if (col_in.ne.12.and.col_in.ne.14.and.col_in.ne.21.and.col_in.ne.10) then
        write(6,*)'col_in must be 10, 12, 14, or 21; good bye.'
        goto 8990
       endif

       !These statements resolve the name of the file with data
       do file_n=first,last
        msf=int(file_n/100)
        nsf=int((file_n-(msf*100))/10)
        lsf=file_n-(100*msf)-(10*nsf)
        one=char(48+msf)
        two=char(48+nsf)
        three=char(48+lsf)
        fname=trim(keyword)//'_'//one//two//three//'.txt'
        open (32, file=trim(fname),status='old') !Open the file
        !These few statements figure out how many lines are in the file (nnline)
        nnline = 0
 10     continue
        read(32,*,end=11) !goto 11 at the end of file
        nnline=nnline+1
        goto 10
 11     continue
        close(32)
        !Open the file again and determine the length of the header
        msf=int(file_n/100) !here, the first file is opened
        nsf=int((file_n-(msf*100))/10)
        lsf=file_n-(100*msf)-(10*nsf)
        one=char(48+msf)
        two=char(48+nsf)
        three=char(48+lsf)
        fname=trim(keyword)//'_'//one//two//three//'.txt'
        open (32, file=trim(fname),status='old')
        !Check each line for the string '#'. If found, we know this line is part
        !of the header.  If not, we have reached the data section.
        header = 0
        hash = 0
	    do i1=1,nnline
	     read(32, "(a)") inline
	     hash = index(inline, '#')
	     if (hash.ne.0) then
	      header=header+1
		  hash=0
	     else
		  goto 14
	     endif
	    enddo
 14    close(32)
       !Inform user of the result
	   write(6,3)'There are',nnline,' lines in file',file_n,'.'

       if (nnline.gt.70000) then  !check to see if the array is large enough to store the data
        write(6,*)'you need to increase the size of the s array.  Bye'
        goto 8990
       endif

       msf=int(file_n/100) !here, the first file is opened
       nsf=int((file_n-(msf*100))/10)
       lsf=file_n-(100*msf)-(10*nsf)
       one=char(48+msf)
       two=char(48+nsf)
       three=char(48+lsf)
       fname=trim(keyword)//'_'//one//two//three//'.txt'
       open (32, file=trim(fname),status='old')
       fname2=trim(keyword)//'_triples_'//one//two//three//'.txt'
       open (33, file=trim(fname2),status='unknown')
       !Read the data into the segment matrix
       !First, read the header and write it the output
       write(33,2)'# written by find_tjs: ',version
       write(33,1)'# text below from the header of the data file'
       do i1=1,header
	    read(32, "(a)") inline
		write(33,"(a)") inline
	   enddo

       !This is where we read the data
       do i1=1,nnline-header
        !The section is used when the input file has 21 columns
        !note that the ordering is selected so that the coordinates
        !always go to positions 9, 10, 11, and 12.
        if (col_in.eq.21) then
         read(32,*,end=21) s(1,i1), s(2,i1), s(3,i1), s(4,i1), s(5,i1),
     &	 s(6,i1), s(7,i1), s(8,i1), s(16,i1), s(17,i1), s(18,i1), s(19,i1),
     &   s(13,i1), s(14,i1), s(15,i1), s(9,i1), s(10,i1), s(11,i1),
     &   s(12,i1), s(20,i1), s(21,i1)
        endif
        !The section is used when the input file has 14 columns
        if (col_in.eq.14) then
         read(32,*,end=21) s(1,i1), s(2,i1), s(3,i1), s(4,i1), s(5,i1),
     &	 s(6,i1), s(7,i1), s(8,i1), s(9,i1), s(10,i1), s(11,i1), s(12,i1),
     &   s(13,i1), s(14,i1)
        endif
        !The section is used when the input file has 12 columns
        if (col_in.eq.12) then
		 read(32,*,end=21) s(1,i1), s(2,i1), s(3,i1), s(4,i1), s(5,i1),
     &	 s(6,i1), s(7,i1), s(8,i1), s(9,i1), s(10,i1), s(11,i1), s(12,i1)
        endif
        !The section is used when the input file has 10 columns
        if (col_in.eq.10) then
		 read(32,*,end=21) s(1,i1), s(2,i1), s(3,i1), s(4,i1), s(5,i1),
     &	 s(6,i1), s(9,i1), s(10,i1), s(11,i1), s(12,i1)
        endif
        
        enddo !This end the loop that reads the data
 21     continue
		close (32)

        !There are occasionally line segments that have zero length.
        !In this section, they are identifies.  They must be removed.
        segs=0
        ZeroLen=0
        do i1=1,nnline-header
         length=sqrt(((s(9,i1)-s(11,i1))**2)+((s(10,i1)-s(12,i1))**2))
         if (length.lt.0.0001) then
          ZeroLen=ZeroLen+1
          write(6,*)i1,s(9,i1),s(10,i1),s(11,i1),s(12,i1)
         else
          segs=segs+1
         endif
        enddo
        if (ZeroLen.gt.0) then
         write(6,*)'there are: ',ZeroLen,' segments with zero length.'
         write(6,*)'at the lines listed above.  They should be deleted before running.  Bye'
         goto 8990
        endif

        !The next big section of code finds all the triple junctions, but it finds
        !most of them more then once, so they have to be sorted.  At this stage,
        !the number of line segments involved and the x,y coordinates are saved in
        !the matrix tj_1

		trips=0
		do i1=1,nnline
         !This is here so that you don't detect surface junctions
         if (s(1,i1).eq.0.0.AND.s(2,i1).eq.0.0) then
          goto 70
         endif
		 match=0
		 start=0
		 do i2=1,nnline
		  if (i2.eq.i1) goto 40
          !This is here so that you don't detect surface junctions
          if (s(1,i2).eq.0.0.AND.s(2,i2).eq.0.0) then
           goto 40
          endif
		  if(match.eq.1) goto 30
          !this is the section where you find the first match
          if(s(9,i1).eq.s(9,i2).and.s(10,i1).eq.s(10,i2)) then
		   match=match+1
		   start=1
		   if(match.eq.1) tem=i2
		   goto 40
		  endif
          if(s(9,i1).eq.s(11,i2).and.s(10,i1).eq.s(12,i2)) then
		   match=match+1
		   start=1
		   if(match.eq.1) tem=i2
		   goto 40
		  endif
          if(s(11,i1).eq.s(9,i2).and.s(12,i1).eq.s(10,i2)) then
		   match=match+1
		   start=2
		   if(match.eq.1) tem=i2
		   goto 40
		  endif
          if(s(11,i1).eq.s(11,i2).and.s(12,i1).eq.s(12,i2)) then
		   match=match+1
		   start=2
		   if(match.eq.1) tem=i2
		   goto 40
		  endif
		  goto 40
          !if you already have the first match, you come here
 30       continue
          !this is if the first match was on the first coordinate
		  if (start.eq.1) then
		   if(s(9,i1).eq.s(9,i2).and.s(10,i1).eq.s(10,i2)) then
		   	trips=trips+1
			tj_1(trips,1)=i1
			tj_1(trips,2)=tem
		    tj_1(trips,3)=i2
		    tj_1(trips,4)=s(9,i1)
		    tj_1(trips,5)=s(10,i1)
			goto 60
		   endif
		   if(s(9,i1).eq.s(11,i2).and.s(10,i1).eq.s(12,i2)) then
		   	trips=trips+1
			tj_1(trips,1)=i1
			tj_1(trips,2)=tem
		    tj_1(trips,3)=i2
		    tj_1(trips,4)=s(9,i1)
		    tj_1(trips,5)=s(10,i1)
			goto 60
		   endif
		  endif
          !this is if the first match was on the second coordinate
		  if (start.eq.2) then
		   if(s(11,i1).eq.s(9,i2).and.s(12,i1).eq.s(10,i2)) then
			trips=trips+1
			tj_1(trips,1)=i1
			tj_1(trips,2)=tem
		    tj_1(trips,3)=i2
		    tj_1(trips,4)=s(11,i1)
		    tj_1(trips,5)=s(12,i1)
			goto 60
		   endif
		   if(s(11,i1).eq.s(11,i2).and.s(12,i1).eq.s(12,i2)) then
			trips=trips+1
			tj_1(trips,1)=i1
			tj_1(trips,2)=tem
		    tj_1(trips,3)=i2
		    tj_1(trips,4)=s(11,i1)
		    tj_1(trips,5)=s(12,i1)
			goto 60
		   endif
		  endif
		  
 40       continue
        enddo	!  This closes the j=1,nnline loop
        if(match.ne.2) goto 70

        !if you found the second match, you end up here
 60     continue
			
 70    continue
       enddo   !  This closes the i1=1,nnline loop

       !this take all the tjs in tj_1, eliminated duplicates, and puts the new
       !list in to tj_2; kk=i3,ll=i4

       tj_2(1,1)=tj_1(1,1)
       tj_2(1,2)=tj_1(1,2)
       tj_2(1,3)=tj_1(1,3)
       tj_2(1,4)=tj_1(1,4)
       tj_2(1,5)=tj_1(1,5)
       i5=1

       do i3=2,trips
	    do i4=1,i5
		 if (tj_1(i3,4).eq.tj_2(i4,4).and.tj_1(i3,5).eq.tj_2(i4,5)) then
           goto 200
	     endif
		enddo ! ends the i4=1,i5 loop
		 i5=i5+1
		 tj_2(i5,1)=tj_1(i3,1)
         tj_2(i5,2)=tj_1(i3,2)
         tj_2(i5,3)=tj_1(i3,3)
         tj_2(i5,4)=tj_1(i3,4)
         tj_2(i5,5)=tj_1(i3,5)
 200    continue
	   enddo ! ends the i3=2,trips loop
       write(33,5)'# I found:',i5,' triple junctions.'
       write(33,1)'# column labels for this file:'
       if (col_out.eq.12) then
        write(33,1)'#   phi1     PHI      phi2     phi1     PHI      phi2     dummy    dummy     x1          y1        x2          y2'
       endif
       if (col_out.eq.10) then
        write(33,1)'#   phi1     PHI      phi2     phi1     PHI      phi2       x1          y1        x2          y2'
       endif

       !The next step is to write out all of the data

       if (col_out.eq.12) then
       do i3=1,i5
        indx=int(tj_2(i3,1))
        !conditional added to prevents possibility of 'zero' lines
        if (s(1,indx).eq.0.0.and.s(2,indx).eq.0.0) then
         goto 300
        endif
        write(33,4)s(1,indx),s(2,indx),s(3,indx),s(4,indx),
     &	s(5,indx),s(6,indx),s(7,indx),s(8,indx),s(9,indx),
     &	s(10,indx),s(11,indx),s(12,indx)
        indx=int(tj_2(i3,2))
        write(33,4)s(1,indx),s(2,indx),s(3,indx),s(4,indx),
     &	s(5,indx),s(6,indx),s(7,indx),s(8,indx),s(9,indx),
     &	s(10,indx),s(11,indx),s(12,indx)
        indx=int(tj_2(i3,3))
        write(33,4)s(1,indx),s(2,indx),s(3,indx),s(4,indx),
     &	s(5,indx),s(6,indx),s(7,indx),s(8,indx),s(9,indx),
     &	s(10,indx),s(11,indx),s(12,indx)
 300    continue
       enddo !  ends the i3=1,i5 loop
       endif

       if (col_out.eq.10) then
       do i3=1,i5
        indx=int(tj_2(i3,1))
        !conditional added to prevents possibility of 'zero' lines
        if (s(1,indx).eq.0.0.and.s(2,indx).eq.0.0) then
         goto 400
        endif
        write(33,6)s(1,indx),s(2,indx),s(3,indx),s(4,indx),
     &	s(5,indx),s(6,indx),s(9,indx),
     &	s(10,indx),s(11,indx),s(12,indx)
        indx=int(tj_2(i3,2))
        write(33,6)s(1,indx),s(2,indx),s(3,indx),s(4,indx),
     &	s(5,indx),s(6,indx),s(9,indx),
     &	s(10,indx),s(11,indx),s(12,indx)
        indx=int(tj_2(i3,3))
        write(33,6)s(1,indx),s(2,indx),s(3,indx),s(4,indx),
     &	s(5,indx),s(6,indx),s(9,indx),
     &	s(10,indx),s(11,indx),s(12,indx)
 400    continue
       enddo !  ends the i3=1,i5 loop
       endif

       if (col_out.ne.10.and.col_out.ne.12) then
        write(6,1)'col_out must equal 10 or 12.  It does not, so no output written.'
        goto 8990
       endif


       close(32)
       close(33)
       write(6,5)'I found:',i5,' triple junctions.'
		
      enddo ! this ends the file_n=first,last loop

 8990 continue
      write(6,1) 'Program complete'
 9000 continue

      end






	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
