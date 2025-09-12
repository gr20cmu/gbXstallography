C*****************************************
       program calculate_triple_line_distribution
C*****************************************
       implicit none
	   include 'common.fi'

       write(6,1) '================================================='
	   write(6,1) 'PROGRAM TO CALCULATE THE ORIENTATION DISTRIBUTION '
       write(6,1) 'OF TRIPLE LINES AND THE GRAIN BOUNDARY PLANE'
       write(6,1) 'DISTRIBUTION FROM LISTS OF TRIPLE LINES'
	   write(6,1) ' '
	   write(6,1) 'rohrer@cmu.edu'
	   write(6,1) 'version 06/03/2025'
       write(6,1) '================================================='
	   
	   version = 'version 06/03/2025'
       !The input file is the same as for the energy reconstruction
       !I found a *.dat file with one extra line inserted - this leads to null result

       !some constants that might be useful
       pi = 4.0*atan(1.0)
       eps = 1.e-6

       !some formating statements used to format typed output
 1     format(A)
 2     format(A,A)
 3     format(A,I8,A)
 4     format(I8)
 5     format(A,I2,A,I2)
 6     format(2I4,2x,f10.6)
 7     format(f10.6)

       !program here
       
       open(21, file='input.txt',status='old') !Open parameter file
       read(21,*)keyword       !name of the file with the triple lines
       read(21,*) msym         !designates crystal symmetry
       read(21,*) CD,CD2       !discretization of orientations, per 90°
       !read(21,*) LessThanTJNum
       close(21)               !close the input file
       write(6,1) keyword      !write parameters to terminal
       write(6,*)msym,CD, CD2
       call GetFileSize (trim(keyword),header,nnline)
       write(6,3)'The data file has ',nnline,' lines.'
       point = index(keyword,'.') - 1

       !start by determining how many triple lines are in the file
       open(30, file=trim(keyword), status='old')! open the file
       do i1 = 1,header
        read(30,*)inline
       enddo
       do i1=1,(nnline/16)  !this loop reads the data
        read(30,4)TJNum
        do i2=1,15
         read(30,1)inline
        enddo
       enddo
       write(6,3)'The data file has ',TJNum,' triple lines.'
       close(30)

       !This is a check that file structure is correct
       if ((nnline/TJNum).ne.int(float(nnline)/float(TJNum))) then
        write(6,1)'There is something wrong with the structure of this file.'
        write(6,1)'The number of lines and number of triple junctions is not consistent.'
        write(6,1)'Ratio of lines to triple junctions must be 16. The program will quit.'
        goto 9000
       else
        continue
       endif

       !now that we have the number of triple lines, we can read
       !each one, find its orientation in the crystal reference
       !frame, and bin the orientation in the distribution

       call get_symop (msym, O, nsymm) !load symmetry operators needed for the calculation
       pd=0.0      !initialize array to store the distributions data
       gbpd=0
       NumObs2=0   !initialize array to store the number of observations
       NumGBs2=0
    
       open(30, file=trim(keyword), status='old')! open the file
       do i1 = 1,header
        read(30,*)inline
       enddo
       do i1=1, TJNum !this loop reads the data
       !do i1=1, LessThanTJNum !for the case of analyzing a number of TLs < TJNum
        read(30,*)TJNum
        read(30,*)tl_s(1), tl_s(2), tl_s(3)
        read(30,*)e1(1), e1(2), e1(3)
        read(30,*)norm1(1), norm1(2), norm1(3)
        read(30,*)e2(1), e2(2), e2(3)
        read(30,*)norm2(1), norm2(2), norm2(3)
        read(30,*)e3(1), e3(2), e3(3)
        read(30,*)norm3(1), norm3(2), norm3(3)
        !The triple line file is in degrees, so convert angles to radians
        call DToRad3 (e1, e1_r)
        call DToRad3 (e2, e2_r)
        call DToRad3 (e3, e3_r)
        !Convert the Euler angles to 3X3 Orientation matrices
        call EToG(e1_r,g_1)
        call EToG(e2_r,g_2)
        call EToG(e3_r,g_3)
        !Calculate the length of the triple line
        call Vlen(tl_s,length)

        do i2=1,nsymm !Each equivalent representation of the TL is considered
		  call symop(O,i2,so_1)!call the i2-th 3x3 symmetry operator
		  call MToM (so_1,g_1,gg_1)!hit the first g with the i2-th symmetry operator
          call MToM (so_1,g_2,gg_2)!hit the second ...
          call MToM (so_1,g_3,gg_3)!hit the third ...
          !classify the triple line in the crystal reference frame of the first grain
          call MToV (gg_1, tl_s, nf)!convert triple line to frame of first crystal
          call VToAngles (nf, sa)!convert triple line direction to spherical angles
          c4=int(float(CD2)*cos(sa(1)))+1 !discretize the angle
          c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1 !discretize the angle
          if (c4.eq.CD2+1) c4=CD2 !Make sure the indices are within bounds
          if (c4.eq.0) c4=1
          if (c5.eq.4*CD2+1) c5=4*CD2
          pd(c4,c5)=pd(c4,c5)+length
          NumObs2(c4,c5)=NumObs2(c4,c5)+1
          !classify the triple line in the crystal reference frame of the second grain
          call MToV (gg_2, tl_s, nf) !this report the process for the first grain
          call VToAngles (nf, sa)
          c4=int(float(CD2)*cos(sa(1)))+1
          c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
          if (c4.eq.CD2+1) c4=CD2
          if (c4.eq.0) c4=1
          if (c5.eq.4*CD2+1) c5=4*CD2
          pd(c4,c5)=pd(c4,c5)+length
          NumObs2(c4,c5)=NumObs2(c4,c5)+1
          !classify the triple line in the crystal reference frame of the third grain
          call MToV (gg_3, tl_s, nf) !this report the process for the first grain
          call VToAngles (nf, sa)
          c4=int(float(CD2)*cos(sa(1)))+1
          c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
          if (c4.eq.CD2+1) c4=CD2
          if (c4.eq.0) c4=1
          if (c5.eq.4*CD2+1) c5=4*CD2
          pd(c4,c5)=pd(c4,c5)+length
          NumObs2(c4,c5)=NumObs2(c4,c5)+1
        enddo ! end the i2=1,nsymm loop



        do i2=1,nsymm !Each equivalent representation of each GB plane is considered
		  call symop(O,i2,so_1)!call the i2-th 3x3 symmetry operator
		  call MToM (so_1,g_1,gg_1)!hit the first g with the i2-th symmetry operator
          call MToM (so_1,g_2,gg_2)!hit the second ...
          call MToM (so_1,g_3,gg_3)!hit the third ...
          !consider the first GB
          call MToV (gg_1, norm1, nf)
          call VToAngles (nf, sa)
          c4=int(float(CD2)*cos(sa(1)))+1
          c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
          if (c4.eq.CD2+1) c4=CD2
          if (c4.eq.0) c4=1
          if (c5.eq.4*CD2+1) c5=4*CD2
          gbpd(c4,c5)=gbpd(c4,c5)+1.0
          NumGBs2(c4,c5)=NumGBs2(c4,c5)+1
          call MToV (gg_2, norm1, nf)
          call VToAngles (nf, sa)
          c4=int(float(CD2)*cos(sa(1)))+1
          c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
          if (c4.eq.CD2+1) c4=CD2
          if (c4.eq.0) c4=1
          if (c5.eq.4*CD2+1) c5=4*CD2
          gbpd(c4,c5)=gbpd(c4,c5)+1.0
          NumGBs2(c4,c5)=NumGBs2(c4,c5)+1
          !consider the second GB
          call MToV (gg_2, norm2, nf)
          call VToAngles (nf, sa)
          c4=int(float(CD2)*cos(sa(1)))+1
          c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
          if (c4.eq.CD2+1) c4=CD2
          if (c4.eq.0) c4=1
          if (c5.eq.4*CD2+1) c5=4*CD2
          gbpd(c4,c5)=gbpd(c4,c5)+1.0
          NumGBs2(c4,c5)=NumGBs2(c4,c5)+1
          call MToV (gg_3, norm2, nf)
          call VToAngles (nf, sa)
          c4=int(float(CD2)*cos(sa(1)))+1
          c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
          if (c4.eq.CD2+1) c4=CD2
          if (c4.eq.0) c4=1
          if (c5.eq.4*CD2+1) c5=4*CD2
          gbpd(c4,c5)=gbpd(c4,c5)+1.0
          NumGBs2(c4,c5)=NumGBs2(c4,c5)+1
          !consider the third GB
          call MToV (gg_3, norm3, nf)
          call VToAngles (nf, sa)
          c4=int(float(CD2)*cos(sa(1)))+1
          c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
          if (c4.eq.CD2+1) c4=CD2
          if (c4.eq.0) c4=1
          if (c5.eq.4*CD2+1) c5=4*CD2
          gbpd(c4,c5)=gbpd(c4,c5)+1.0
          NumGBs2(c4,c5)=NumGBs2(c4,c5)+1
          call MToV (gg_1, norm3, nf)
          call VToAngles (nf, sa)
          c4=int(float(CD2)*cos(sa(1)))+1
          c5=int(float(4*CD2)*sa(2)/(pi*2.0))+1
          if (c4.eq.CD2+1) c4=CD2
          if (c4.eq.0) c4=1
          if (c5.eq.4*CD2+1) c5=4*CD2
          gbpd(c4,c5)=gbpd(c4,c5)+1.0
          NumGBs2(c4,c5)=NumGBs2(c4,c5)+1

        enddo ! end the i2=1,nsymm loop

       enddo ! end the i1=1, TJNum loop
       close(30) ! close the data file

       !Here we write the data to the output file - open file, write header
       open(43, file=keyword(1:point)//'_tld.txt',status='unknown')
       write(43,2)'# written by calc_TL_GBP_dist, ',version
       write(43,2)'# based on data in file: ',trim(keyword)
       write(43,5)'# Discretization = ',CD,', ',CD2
       open(44, file=keyword(1:point)//'_gbpd.txt',status='unknown')
       write(44,2)'# written by calc_TL_GBP_dist, ',version
       write(44,2)'# based on data in file: ',trim(keyword)
       write(44,5)'# Discretization = ',CD,', ',CD2

       !determine the average value of the TL distribution
       avg=0.0
       do i1=1,CD2
        do i2=1,4*CD2
         avg=avg+pd(i1,i2)
        enddo
       enddo
       avg=avg/float(4*CD2*CD2)
       write(6,*)avg
       !normalize the distribution
       do i1=1,CD2
        do i2=1,4*CD2
         pd(i1,i2)=pd(i1,i2)/avg
         write(43,7)pd(i1,i2)
        enddo
       enddo
       close(43)


      !determine the average value of the GB plane distribution
      gbavg=0.0
      do i1=1,CD2
       do i2=1,4*CD2
        gbavg=gbavg+gbpd(i1,i2)
       enddo
      enddo
      gbavg=gbavg/float(4*CD2*CD2)
      !write the values in MRD units
      do i1=1,CD2
       do i2=1,4*CD2
        gbpd(i1,i2)=gbpd(i1,i2)/gbavg
       write (44,*) gbpd(i1,i2)
       enddo
      enddo
      close(44)


 9000  continue
       write(6,1) 'Program complete'

       end






	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
	   
