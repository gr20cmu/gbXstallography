C*****************************************
       program condition segments
C*****************************************
       implicit none
       include 'common.fi'
	   
       write(6,5) '=============================================='
	   write(6,5) ' PROGRAM PREPARES GRAIN BOUNDARY SEGMENT FILES'
	   write(6,5) ' FOR ANALYSIS BY WRITING THEM IN A 10 OR 12  '
	   write(6,5) ' COLUMN FORMAT.  MULTIPLE FILES CAN BE '
	   write(6,5) ' CONCATENATED INTO A SINGLE FILE.'
	   write(6,5)
	   write(6,5) ' rohrer@cmu.edu'
	   write(6,5) ' version 05/07/22'
       write(6,5) '=============================================='

       !compile with the make file.  You can also use gfortran main.f -O3 -o condition_segs
       !9/9/14: fix issue with 14 column files
       !12/01/16: can handle more than 100 files
       !02/10/19: Added make file, comments, and instructions for single
       !file processing.  Also changed name from combine_segs to condition_segs
       !04/03/19: Added a filter to remove segments of zero length, a new feature found in
       !output from TSLv8
       !03/09/22: added feature to accept 10 column input from extract gb traces
       !03/26/22: added feature to eliminate 0 0 0 Euler angles
       !          and to condition files without concatenation
       !05/07/22: improved parameter checks and fixed concatenation option
       version = 'version 05/07/22'

       !statements to format the input or output and referred to later
 5     format(A)
 12    format(A,I7,A,I3,A)
 81    format(A,A)
 82    format(7f9.3,f8.2,4f11.3)
 83    format(6f9.3,4f11.3)

       !The following statements are needed to read the input parameters
       open(21, file='input.txt', status='old') ! open file with input parameters
       read(21,*) keyword   !read the base name of the data file(s)
       read(21,*) out_fname !read the name you want to call the output file
       read(21,*) first     !read the first file number
	   read(21,*) last      !read the last file number
       read(21,*) col_in    !read the number of columns in the data file
       read(21,*) col_out   !the number of columns to be written in the output data
       read(21,*) concat    !flag to concatenate files (1) or not (0)
       read(21,5) inline2   !read user added comment
       close(21)            !cloes the file with input parameters

       !the next section of code just checks the validity of the parameters
       !end program if the number of columns in the input is not allowed.
       if (col_in.ne.10.AND.col_in.ne.12.AND.col_in.ne.14.AND.col_in.ne.21) then
        write(6,5)'col_in must be 10, 12, 14, or 21 - correct input.txt; good bye.'
        goto 9020
       endif
       !end program if the number of columns in the output is not allowed.
       if (col_out.ne.10.AND.col_out.ne.12) then
        write(6,5)'col_out must be 10 or 12 - correct input.txt; good bye.'
        goto 9020
       endif
       !end program if col_out is not less than or equal to col_in.
       if (col_out.gt.col_in) then
        write(6,5)'col_out must be less than or equal to col_in - correct input.txt; good bye.'
        goto 9020
       endif
       !end program if the first file number is not allowed.
       if (first.lt.0.or.first.gt.998) then
        write(6,5)'first must be between 0 and 998 - correct input.txt; good bye.'
        goto 9020
       endif
       !end program if the last file number is not allowed.
       if (last.lt.0.or.last.gt.999) then
        write(6,5)'last must be between 0 and 999 - correct input.txt; good bye.'
        goto 9020
       endif
       !end program if last is not greater than or equal to first.
       if (last.lt.first) then
        write(6,5)'last must be greater than or equal to first - correct input.txt; good bye.'
        goto 9020
       endif

       !This is the main loop over the files
       do file_n=first,last
        !Begin with the first file; first we have to create a string for the file name
        msf=int(file_n/100) !determine the most significant digit in the three digit file number
        nsf=int((file_n-(msf*100))/10) !determine the next (second) most significant digit
        lsf=file_n-(100*msf)-(10*nsf)  !determine the least (third) most significant digit
        one=char(48+msf)               !make it a string
        two=char(48+nsf)               !make it a string
        three=char(48+lsf)             !make it a string
        !create the filename by concatenating the base name, _, and the three digits of the file number.
        fname=trim(keyword)//'_'//one//two//three//'.txt'
        open (32, file=fname,status='old')  !Open this data file

        !We need to know how many lines are in the file, so we count them.
		nnline = 0           !Zero the counter
 100    continue
		read(32,*,end=101)   !read a line; if end of file, go to 101
		nnline=nnline+1      !if there was something on the line, increment the counter
		goto 100             !check to see if there is another line
 101    continue             !you land here when the end of file has been reached
        close(32)            !close the data file
		write(6,12) 'There are',nnline,' lines in file',file_n,'.'! report the result

        !If there are more lines than allowed in the data structure, 
        !stop the program before it fails.
        if (nnline.gt.3000000) then
         write(6,5)'you must increase the size of the s array in common, edit line above, and recompile.  Bye.'
         goto 9020
		endif

        if (concat.ne.0) then
         !open the file where concatenated output is written
	     open (34, file=out_fname, status='unknown')
        else
         !or open the file where unconcatenated output is written
	     open (34, file=trim(keyword)//'_cnd_'//one//two//three//'.txt', status='unknown')
        endif
        if (file_n.eq.first.OR.concat.eq.0) then
	     write(34,81)'# segment file written by condition_segs: ',version
        endif

        !Now that we know how many lines are in the file, we open it again
        open (32, file=fname,status='old')

        !Here we determine the number of lines in the Header.
        !Check each line for the string '#'. If found, we know this line is part
        !f the header.  If not, we have reached the data section.
        header = 0                      !Zero a counter
        hash = 0                        !Zero a indicator
        do i=1,nnline                   !loop over every line of the file
	     read(32, "(a)") inline         !read a line of the file
	     hash = index(inline, '#')      !sets hash to 1 if # (mark for header) is in the line
	     if (hash.ne.0) then            !when hash is 1
	      header=header+1               !# was in line, so count as a line of header
		  hash=0                                   !reset the indicator
		  if (file_n.eq.first.OR.concat.eq.0) then !We write the header to the output
		   write(34, "(a)")inline                  !write the line to the file
		  endif
	     else
		  goto 128                      !If not, we skip
	     endif			
	    enddo
 128   continue
       close(32)                        !close the data file
       write(6,*)header,' lines in the header'    !report the result

       !Add some useful information to the header
       if (file_n.eq.first.OR.concat.eq.0) then
	    write(34,81)'# comment: ',inline2
        write(34,5)'#   phi1     PHI      phi2     phi1     PHI      phi2        x1         y1         x2         y2 '
       endif

       !Now that we know how many lines are in the file, 
       !and we know the header length, we open it again
       open (32, file=fname,status='old')
       !this loop reads the header
       do i=1,header
        read(32, "(a)") inline
       enddo

       !This section reads the data into the matrix, s
       !The loop goes from the first line after the header until
       !the end of the file.
       do i=header+1,nnline
         !The section is used when the input file has 21 columns
         if (col_in.eq.21) then
         read(32,*,end=200) s(1,i), s(2,i), s(3,i), s(4,i), s(5,i),
     &	 s(6,i), s(7,i), s(8,i), s(9,i), s(10,i), s(11,i), s(12,i),
     &   s(13,i), s(14,i), s(15,i), s(16,i), s(17,i), s(18,i),
     &   s(19,i), s(20,i), s(21,i)
         endif
         !The section is used when the input file has 14 columns
         if (col_in.eq.14) then
         read(32,*,end=200) s(1,i), s(2,i), s(3,i), s(4,i), s(5,i),
     &	 s(6,i), s(7,i), s(8,i), s(16,i), s(17,i), s(18,i), s(19,i),
     &   s(13,i), s(14,i)
         endif
         !The section is used when the input file has 12 columns
         if (col_in.eq.12) then
		 read(32,*,end=200) s(1,i), s(2,i), s(3,i), s(4,i), s(5,i),
     &	 s(6,i), s(7,i), s(8,i), s(16,i), s(17,i), s(18,i), s(19,i)
	     endif
        !The section is used when the input file has 10 columns
        if (col_in.eq.10) then
		 read(32,*,end=200) s(1,i), s(2,i), s(3,i), s(4,i), s(5,i),
     &	 s(6,i), s(16,i), s(17,i), s(18,i), s(19,i)
        endif
       enddo ! end the i=header+1,nnline loop

 200   continue
       close (32) !The data has been read, so close the file

       !Here we write out the line segments data
       !The conditionals are set so that it only writes out
       !line segments of non-zero length and Euler angles that are not 0 0 0.
       !Segments with these characteristics sometime show up for a variety of reasons.
       if (col_out.eq.12) then
        do i=header+1,nnline
         if ((s(16,i).eq.s(18,i)).and.(s(17,i).eq.s(19,i))) then
          goto 300
         endif
         if ((s(1,i).eq.0.0).and.(s(4,i).eq.0.0)) then
          goto 300
         endif
         write(34,82) s(1,i), s(2,i), s(3,i), s(4,i), s(5,i),
     &	  s(6,i), s(14,i), s(15,i), s(16,i), s(17,i), s(18,i), s(19,i)
 300     continue
        enddo
       else
        do i=header+1,nnline
         if ((s(16,i).eq.s(18,i)).and.(s(17,i).eq.s(19,i))) then
          goto 400
         endif
         if ((s(1,i).eq.0.0).and.(s(4,i).eq.0.0)) then
          goto 400
         endif
         write(34,83) s(1,i), s(2,i), s(3,i), s(4,i), s(5,i),
     &	 s(6,i), s(16,i), s(17,i), s(18,i), s(19,i)
 400     continue
        enddo
       endif
       if (concat.eq.0) then
        close(34)
       else
        continue
       endif

      enddo ! This is the file_n loop

      close(34) !all files have been read, so close the output file
 9020 continue
      end
