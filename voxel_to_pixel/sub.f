c=====================================================================
	subroutine EToG (euler, g)
c=====================================================================
	include 'common.fi'
       
	g(1,1)=cos(euler(1))*cos(euler(3))
     &	                 -sin(euler(1))*sin(euler(3))*cos(euler(2))
	g(1,2)=sin(euler(1))*cos(euler(3))
     &	                 +cos(euler(1))*sin(euler(3))*cos(euler(2))
	g(1,3)=sin(euler(3))*sin(euler(2))
	g(2,1)=-cos(euler(1))*sin(euler(3))
     &	                 -sin(euler(1))*cos(euler(3))*cos(euler(2))
	g(2,2)=-sin(euler(1))*sin(euler(3))
     &	                 +cos(euler(1))*cos(euler(3))*cos(euler(2))
	g(2,3)=cos(euler(3))*sin(euler(2))
	g(3,1)=sin(euler(1))*sin(euler(2))
	g(3,2)=-cos(euler(1))*sin(euler(2))
	g(3,3)=cos(euler(2))


	return
	end




c=====================================================================
	subroutine GToE (delg, euler)
c=====================================================================
	include 'common.fi'
	
	xz=delg(3,3)
	
	if (xz*xz.lt.1.0) then
	sf=sqrt(1.0-xz*xz)
	else
	sf=0.0
	endif

	if (sf.gt.eps) then
	euler(1)=-delg(3,2)/sf
	euler(1)=acos2(euler(1))
	if (delg(3,1).lt.0.0) euler(1)=2.0*pi-euler(1)
	euler(2)=acos2(xz)
	euler(3)=delg(2,3)/sf
	euler(3)=acos2(euler(3))
	if (delg(1,3).lt.0.0) euler(3)=2.0*pi-euler(3)
	else
	euler(1)=delg(2,2)
	euler(1)=acos2(euler(1))
	if (delg(1,2).lt.0.0) euler(1)=2.0*pi-euler(1)
	euler(2)=0.0
	euler(3)=0.0
	endif


	return
	end




c=====================================================================
	subroutine MtoM (g1, g2, g3)
c=====================================================================
	include 'common.fi'

	do ii=1,3
	do jj=1,3
	g3(ii,jj)=0.0
	do kk=1,3
	g3(ii,jj)=g3(ii,jj)+g1(ii,kk)*g2(kk,jj)
	enddo
	enddo
	enddo

	return
	end




c=====================================================================
	subroutine MToV (g, n1, n2)
c=====================================================================
    include 'common.fi'

	do ii=1,3
	n2(ii)=0.0
	do jj=1,3
	n2(ii)=n2(ii)+g(ii,jj)*n1(jj)
	enddo
	enddo

	return
	end




c=====================================================================
	subroutine AnglesToV (theta, psi, vector_out)
c=====================================================================
	include 'common.fi'
	
	vector_out(1)=0.0
	vector_out(2)=0.0
	vector_out(3)=0.0
	
	vector_out(1)=sin(theta)*cos(psi)
	vector_out(2)=sin(theta)*sin(psi)
	vector_out(3)=cos(theta)

	return
	end




c=====================================================================
	subroutine VToAngles (v, a)
c=====================================================================
	include 'common.fi'

	if (v(3).lt.0.0) then
	do ii=1,3
	v(ii)=-v(ii)
	enddo
	endif

	a(1)=acos2(v(3))
c	if (v(1).ne.0.0) then
c	a(2)=atan(v(2)/v(1))
c	else
c	a(2)=pi*0.5
c	endif
	a(2)=atan2(v(2),v(1))
	if (a(2).lt.0.0) a(2)=2.0*pi+a(2)

	return
	end




c=====================================================================
	subroutine AAToG (x, phi, g)
c=====================================================================
	include 'common.fi'
	
	g(1,1)=cos(phi)+(1.0-cos(phi))*(x(1)**2)
	g(1,2)=x(1)*x(2)*(1.0-cos(phi))-x(3)*sin(phi)
	g(1,3)=x(1)*x(3)*(1.0-cos(phi))+x(2)*sin(phi)
	g(2,1)=x(1)*x(2)*(1.0-cos(phi))+x(3)*sin(phi)
	g(2,2)=cos(phi)+(1.0-cos(phi))*(x(2)**2)
	g(2,3)=x(3)*x(2)*(1.0-cos(phi))-x(1)*sin(phi)
	g(3,1)=x(1)*x(3)*(1.0-cos(phi))-x(2)*sin(phi)
	g(3,2)=x(2)*x(3)*(1.0-cos(phi))+x(1)*sin(phi)
	g(3,3)=cos(phi)+(1.0-cos(phi))*(x(3)**2)


	return
	end




c=====================================================================
	subroutine trans (g, gt)
c=====================================================================
	include 'common.fi'

	do ii=1,3
	do jj=1,3
	gt(ii,jj)=g(jj,ii)
	enddo
	enddo

	return
	end




c=====================================================================
	subroutine symop (op, i_s, so)
c=====================================================================
	include 'common.fi'
	

	do ii=1,3
	do jj=1,3
	so(ii,jj)=op(ii,jj,i_s)
	enddo
	enddo

	return
	end




c=====================================================================
      function acos2(value)
c=====================================================================
	include 'common.fi'
	
      temp=value
	if (temp.gt.1.0) temp=1.0
	if (temp.lt.-1.0) temp=-1.0
	acos2=acos(temp)
	
	return
	end
	
	
