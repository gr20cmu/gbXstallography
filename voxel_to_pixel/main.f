C*****************************************
       program vox_to_pix
C*****************************************
       implicit none
	   include 'common.fi'

       write(6,*) '================================================='
	   write(6,*) ' PROGRAM TAKES 3D ORIENTATION MAPS AND WRITES '
	   write(6,*) ' THE DATA AS PARALLEL 2D MAPS IN .ANG STYLE'
	   write(6,*) ' '
	   write(6,*) ' rohrer@cmu.edu'
	   write(6,*) ' version 03/25/2022'
       write(6,*) '================================================='
	   
       !version = 'version 02/23/2022' !started based on align 3D
       !version = 'version 02/26/2022' !modified for D3D output, read parameters from input.txt
       version = 'version 03/25/2022'  !modified to have keywords other than 8 characters

       !constant used later
       pi = 4.0*atan(1.0)
       eps = 1.e-6

       !used to format typed output
  1    format(5f13.5)
  2    format(A,A)
  3    format(A,2f9.4,2I5)
  22   format(f6.2)
  23   format(A)
  24   format(9f7.3)
  25   format(3I5,5f7.3)

       !read the data
       open(21, file='input.txt', status='old')
       read(21,*) in_fname
       read(21,*) keyword
       read(21,23) comment
       read(21,*) XLimit
       read(21,*) YLimit
       read(21,*) ZLimit
       read(21,*) XSTEP
       read(21,*) YSTEP
       close (21)
       grain=0.0

       !classifying grains
       open(25, file=in_fname, status='old')
       write(6,23)'reading data ...'
       do k=1,ZLimit
        do j=1,YLimit
         do i=1,XLimit
          read(25,*)e1(1),e1(2),e1(3),xp,yp,zp,tem
          if (e1(2).eq.0.0.AND.e1(3).eq.0.0)then
           e1(1)=0.0
          endif
          grain(int(tem),1)=e1(1)
          grain(int(tem),2)=e1(2)
          grain(int(tem),3)=e1(3)
         enddo
        enddo
       enddo
       close(25)


       open(25, file=in_fname, status='old')
       do k=1,ZLimit
        msf=int(k/100) !determine the most significant digit in the three digit file number
        nsf=int((k-(msf*100))/10) !determine the next (second) most significant digit
        lsf=k-(100*msf)-(10*nsf)  !determine the least (third) most significant digit
        one=char(48+msf)              !make it a string
        two=char(48+nsf)              !make it a string
        three=char(48+lsf)            !make it a string
        !create the filename by concatenating the base name, _, and the three digits of the file number.
        fname=trim(keyword)//'_'//one//two//three//'.ang'
        write(6,2)'working on: ',fname
        open (30, file=fname,status='unknown') !This opens the first data file
        write(30,2)'# file written by voxel_to_pixel: ',version
        write(30,2)'# This is a pseudo *.ang file.  It was written based on data in: ',in_fname
        write(30,23)'# The output is suitable for extract_gb_traces'
        write(30,3)'# the XSTEP, YSTEP, X pixels, and Y pixels are: ',XSTEP, YSTEP, XLimit, YLimit
        write(30,2)'# ',trim(comment)
        write(30,23)'#        phi1          PHI         phi2            X            Y'
        
        do j=1,YLimit
         do i=1,XLimit
          read(25,*)e1(1),e1(2),e1(3),xp,yp,zp,tem
          write(30,1)grain(int(tem),1),grain(int(tem),2),grain(int(tem),3),xp*XSTEP,yp*YSTEP
         enddo
        enddo
        close(30)
       enddo
       close(25)

       write(6,23)'program complete'
 9000  continue

       end






	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
