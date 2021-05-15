module data_common
	interface display
		subroutine int_display(text,i)
			integer, intent(out) :: i
			character, intent(in) :: text
		end subroutine int_display
		subroutine real_display(text,vi)
			real, intent(out) :: vi
			character, intent(in) :: text
		end subroutine real_display
		subroutine dp_display(text,vi)
			real(kind=8), intent(out) :: vi
			character, intent(in) :: text
		end subroutine dp_display
		subroutine text_display(text,text_out)
		character (len=*) :: text
		character (len=*) :: text_out
		end subroutine text_display
	end interface

	integer, parameter	:: maxlen=256, len_display=50, max_constit=30, max_nr=10, max_conv=100
	character (len=1)	::	test_stop, DIR_SEP_S
	character (len=5)	::	OS_S
	integer(kind=2)	:: int2
	integer	:: &
		c_bwhite=15, c_yellow=14, c_lmagenta=13, c_lred=12, c_lgreen=10, c_lcyan=11, &
		c_lblue=9, c_white=7, c_brown=6, c_red=4, c_cyan=3, c_green=2
	real(kind=8)	:: pk
	real(kind=8)	:: epsilon=1.D-8
end module data_common


SUBROUTINE actp (molal, gact, aw, fi, temp)
	use data_common
	implicit none
	real(kind=8), dimension(0:*)	:: molal, gact
	real(kind=8)	:: aw, fi, temp

	double precision, external	:: temperature, j0, j1, g0, g1 ! functions

	integer, save	:: nc,na,nn,ndepact
	integer, allocatable, save, dimension (:)	:: nzc, nza
	real(kind=8), save	:: ap0,bp0,mh2o
	real(kind=8), allocatable, save, dimension (:,:) :: b0, b1, b2, c0, tc, ta, lc, la
	real(kind=8), allocatable, save, dimension (:,:,:) :: sc, sa, xi

	real(kind=8), allocatable, dimension (:)	:: gc, ga, gn
	real(kind=8), allocatable, dimension (:,:)	:: ec,fc,xc,ea,fa,xa,pp,p,pf,qp,q,qf,cc,bf,b,bp

	integer	:: i,j,k,ii,jj,alloc_res
	character (len=maxlen)	:: x_S
	real(kind=8)	:: at, bt, ct, dt, et, u, z, fj, w, f, v, s, co
	real(kind=8)	:: c(0:9), a(0:11), h(0:3)

  	c(1) = molal(1); c(2) = molal(2); c(3) = molal(3)
  	c(4) = molal(4); c(5) = molal(22); c(6) = molal(5)
  	c(7) = molal(18); c(8) = molal(23); c(9) = molal(11)
  	a(1) = molal(6); a(2) = molal(7); a(3) = molal(25)
  	a(4) = molal(12); a(5) = molal(14); a(6) = molal(13)
  	a(7) = molal(24); a(8) = molal(19); a(9) = molal(20)
  	a(10) = molal(21); a(11) = molal(8)
  	h(1) = molal(10); h(2) = molal(9); h(3) = molal(15)

	IF (ndepact == 0) THEN
    	OPEN (10,file="coefft4",action='READ')
      	read(10,*)  nc, na, nn

		allocate (nzc(0:nc), nza(0:na), stat=alloc_res)
		allocate (b0(0:nc, 0:na), b1(0:nc, 0:na), b2(0:nc, 0:na), c0(0:nc, 0:na), stat=alloc_res)
		allocate (sc(0:nc, 0:nc, 0:na), sa(0:na, 0:na, 0:nc), stat=alloc_res)
		allocate (tc(0:nc, 0:nc), ta(0:na, 0:na), stat=alloc_res)
		allocate (lc(0:nn, 0:nc), la(0:nn, 0:na), xi(0:nn, 0:nc, 0:na), stat=alloc_res)
		
		nzc=0; nza=0
		b0=0.; b1=0.; b2=0.; c0=0.; sc=0.; sa=0.; tc=0.; ta=0.; lc=0.; la=0.; xi=0.
		
      	DO i = 1 , nc; read(10,*)  x_S, nzc(i); ENDDO
      	DO i = 1 , na; read(10,*)  x_S, nza(i); ENDDO
		read(10,*)  x_S, at, bt, ct, dt, et
        ap0 = temperature(at, bt, ct, dt, et, temp)
      	DO i = 1 , nc
        	DO j = 1 , na
          		read(10,*)  x_S, at, bt, ct, dt, et
          		b0(i, j) = temperature(at, bt, ct, dt, et, temp)
          		read(10,*)  x_S, at, bt, ct, dt, et
          		b1(i, j) = temperature(at, bt, ct, dt, et, temp)
          		read(10,*)  x_S, at, bt, ct, dt, et
         		 b2(i, j) = temperature(at, bt, ct, dt, et, temp)
          		read(10,*)  x_S, at, bt, ct, dt, et
          		c0(i, j) = temperature(at, bt, ct, dt, et, temp)
        	ENDDO
      	ENDDO
      	DO i = 1 , nc - 1
        	DO j = i + 1 , nc
          		read(10,*)  x_S, at, bt, ct, dt, et
          		tc(i, j) = temperature(at, bt, ct, dt, et, temp)
          		tc(j, i) = tc(i, j)
        	ENDDO
      	ENDDO
      	DO i = 1 , na - 1
        	DO j = i + 1 , na
          		read(10,*)  x_S, at, bt, ct, dt, et
          		ta(i, j) = temperature(at, bt, ct, dt, et, temp)
          		ta(j, i) = ta(i, j)
        	ENDDO
      	ENDDO
      	DO k = 1 , nc - 1
        	DO i = k + 1 , nc
         	 DO j = 1 , na
				read(10,*)  x_S, at, bt, ct, dt, et
           		sc(k, i, j) = temperature(at, bt, ct, dt, et, temp)
            	sc(i, k, j) = sc(k, i, j)
          	 ENDDO
        	ENDDO
      	ENDDO
      	DO k = 1 , na - 1
        	DO i = k + 1 , na
          		DO j = 1 , nc
            		read(10,*)  x_S, at, bt, ct, dt, et
            		sa(k, i, j) = temperature(at, bt, ct, dt, et, temp)
            		sa(i, k, j) = sa(k, i, j)
          		ENDDO
        	ENDDO
      	ENDDO
      	DO i = 1 , nn
        	DO j = 1 , nc
            	read(10,*)  x_S, at, bt, ct, dt, et
            	lc(i, j) = temperature(at, bt, ct, dt, et, temp)
        	ENDDO
      	ENDDO
      	DO i = 1 , nn
       	 DO j = 1 , na
            read(10,*)  x_S, at, bt, ct, dt, et
            la(i, j) = temperature(at, bt, ct, dt, et, temp)
      	 ENDDO
      	ENDDO
    	CLOSE (10)
    	DO k = 1 , nn
      		DO i = 1 , nc
       			DO j = 1 , na
          			xi(k, i, j) = 0.
        		ENDDO
      		ENDDO
    	ENDDO
    	xi(2, 9, 1) = -.0102D0
    	xi(2, 1, 2) = .046D0
  	END IF

	allocate (ec(0:nc, 0:nc), stat=alloc_res)
	allocate (fc(0:nc, 0:nc), xc(0:nc, 0:nc), ea(0:na, 0:na), fa(0:na, 0:na), xa(0:na, 0:na), stat=alloc_res)
	allocate (pp(0:nc, 0:nc), p(0:nc, 0:nc), pf(0:nc, 0:nc), qp(0:na, 0:na), q(0:na, 0:na), qf(0:na, 0:na), stat=alloc_res)
	allocate (cc(0:nc, 0:na), bf(0:nc, 0:na), b(0:nc, 0:na), bp(0:nc, 0:na), stat=alloc_res)
	allocate (gc(0:nc), ga(0:na), gn(0:nn), stat=alloc_res)
	
	ec=0.; fc=0.; xc=0.; ea=0.; fa=0.; xa=0.; pp=0.; p=0.; pf=0.; qp=0.; q=0.; qf=0.
	cc=0.; bf=0.; b=0.; bp=0.; gc=0.; ga=0.; gn=0.
	
  	bp0 = 1.2D0; mh2o = 55.51D0
  	u = 0.; z = 0.
  	DO i = 1 , nc; u = u + c(i) * nzc(i) ** 2; z = z + c(i) * nzc(i); ENDDO
  	DO j = 1 , na; u = u + a(j) * nza(j) ** 2; z = z + a(j) * nza(j); ENDDO
  	fi = u / 2.; fj = SQRT(fi)
  	u = 6. * ap0 * fj
  	DO i = 1 , nc - 1
    	DO j = i + 1 , nc
      		IF (nzc(i) == nzc(j)) THEN
        		ec(i, j) = 0.; fc(i, j) = 0.
     		 ELSE
       			 xc(i, j) = 2. * u; xc(i, i) = nzc(i) ** 2 * u; xc(j, j) = nzc(j) ** 2 * u
       			 ec(i, j) = (j0(xc(i, j)) - j0(xc(i, i)) / 2. - j0(xc(j, j)) / 2.) / fi / 2.
        		 fc(i, j) = (xc(i, j) * j1(xc(i, j)) - xc(i, i) * j1(xc(i, i)) / 2. - xc(j, j) * j1(xc(j, j)) / 2.) &
        		                    / fi ** 2 / 4. - ec(i, j) / fi
        		 ec(j, i) = ec(i, j); fc(j, i) = fc(i, j)
				 if(i == 1 .and. j == 4) then
				 endif
      		END IF
  	  	ENDDO
  	ENDDO
  	DO i = 1 , na - 1
    	DO j = i + 1 , na
      		IF (nza(i) == nza(j)) THEN
       		 	ea(i, j) = 0.; fa(i, j) = 0.
      		ELSE
        		xa(i, j) = 2. * u; xa(i, i) = nza(i) ** 2 * u; xa(j, j) = nza(j) ** 2 * u
        		ea(i, j) = (j0(xa(i, j)) - j0(xa(i, i)) / 2. - j0(xa(j, j)) / 2.) / fi / 2.
        		fa(i, j) = &
				 (xa(i, j) * j1(xa(i, j)) - xa(i, i) * j1(xa(i, i)) / 2. - xa(j, j) * j1(xa(j, j)) / 2.) &
				  / fi ** 2 / 4. - ea(i, j) / fi
      	  		ea(j, i) = ea(i, j); fa(j, i) = fa(i, j)
      		END IF
    	ENDDO
  	ENDDO
  	DO i = 1 , nc - 1
    	DO j = i + 1 , nc
      		pp(i, j) = fc(i, j); p(i, j) = tc(i, j) + ec(i, j); pf(i, j) = p(i, j) + pp(i, j) * fi
      		pp(j, i) = pp(i, j); p(j, i) = p(i, j); pf(j, i) = pf(i, j)
    	ENDDO
  	ENDDO
  	DO i = 1 , na - 1
    	DO j = i + 1 , na
      		qp(i, j) = fa(i, j); q(i, j) = ta(i, j) + ea(i, j); qf(i, j) = q(i, j) + qp(i, j) * fi
      		qp(j, i) = qp(i, j); q(j, i) = q(i, j); qf(j, i) = qf(i, j)
    	ENDDO
  	ENDDO
  	w = fj * 12.
  	DO i = 1 , nc
    	DO j = 1 , na
      		cc(i, j) = c0(i, j) / SQRT(real(nzc(i) * nza(j))) / 2.
      		IF (nzc(i) == 2 .AND. nza(j) == 2)  v = fj * 1.4D0
      		IF (nzc(i) == 1 .OR. nza(j) == 1)  v = fj * 2.
      		bf(i, j) = b0(i, j) + b1(i, j) * EXP(-v) + b2(i, j) * EXP(-w)
      		b(i, j) = b0(i, j) + b1(i, j) * (g0(v)) + b2(i, j) * (g0(w))
      		bp(i, j) = b1(i, j) * (g1(v)) / fi + b2(i, j) * (g1(w)) / fi
    	ENDDO
  	ENDDO
  	f = -ap0 * (fj / (1 + bp0 * fj) + 2. / bp0 * LOG(1 + bp0 * fj))
  	DO i = 1 , nc; DO j = 1 , na; f = f + c(i) * a(j) * bp(i, j); ENDDO; ENDDO
  	DO i = 1 , nc - 1; DO j = i + 1 , nc; f = f + c(i) * c(j) * pp(i, j); ENDDO; ENDDO
  	DO i = 1 , na - 1; DO j = i + 1 , na; f = f + a(i) * a(j) * qp(i, j); ENDDO; ENDDO
  	DO ii = 1 , nc
    	u = nzc(ii) ** 2 * f
    	DO j = 1 , na; u = u + a(j) * (b(ii, j) * 2. + z * cc(ii, j)); ENDDO
      	DO i = 1 , nc
        	IF (i /= ii) THEN
         	 v = 0.
          	 DO j = 1 , na; v = v + a(j) * sc(ii, i, j); ENDDO
          	 u = u + c(i) * (p(ii, i) * 2. + v)
        	END IF
      	ENDDO
      	DO i = 1 , na - 1; DO j = i + 1 , na; u = u + a(i) * a(j) * sa(i, j, ii); ENDDO; ENDDO
      	DO i = 1 , nc; DO j = 1 , na; u = u + c(i) * a(j) * cc(i, j) * nzc(ii); ENDDO; ENDDO
      	DO i = 1 , nn
        	u = u + h(i) * lc(i, ii) * 2.
      	ENDDO
      	DO k = 1 , nn
        	DO j = 1 , na
          		u = u + h(k) * a(j) * xi(k, ii, j)
        	ENDDO
      	ENDDO
     	gc(ii) = EXP(u)
  	ENDDO 
  	DO jj = 1 , na
    	u = nza(jj) ** 2 * f
    	DO i = 1 , nc; u = u + c(i) * (b(i, jj) * 2. + z * cc(i, jj)); ENDDO
      	DO i = 1 , na
        	IF (i /= jj) THEN
          		v = 0.
          		DO j = 1 , nc; v = v + c(j) * sa(jj, i, j); ENDDO
          		u = u + a(i) * (q(jj, i) * 2. + v)
        	END IF
      	ENDDO
      	DO i = 1 , nc - 1; DO j = i + 1 , nc; u = u + c(i) * c(j) * sc(i, j, jj); ENDDO; ENDDO
      	DO i = 1 , nc; DO j = 1 , na; u = u + c(i) * a(j) * cc(i, j) * nza(jj); ENDDO; ENDDO
      	DO j = 1 , nn
        	u = u + h(j) * la(j, jj)
      	ENDDO
      	DO k = 1 , nn
        	DO i = 1 , nc
          		u = u + h(k) * c(i) * xi(k, i, jj)
        	ENDDO
      	ENDDO
    	ga(jj) = EXP(u)
  	ENDDO 
  	DO k = 1 , nn
    	u = 0.
    	DO i = 1 , nc
      		u = u + c(i) * lc(k, i) * 2.
    	ENDDO
    	DO j = 1 , na
     		u = u + a(j) * la(k, j) * 2.
    	ENDDO
    	DO i = 1 , nc
      		DO j = 1 , na
        		u = u + c(i) * a(j) * xi(k, i, j)
      		ENDDO
    	ENDDO
    	gn(k) = EXP(u)
  	ENDDO
  	u = -ap0 * fi ** 1.5D0 / (1. + bp0 * fj)
  	DO i = 1 , nc; DO j = 1 , na; u = u + c(i) * a(j) * (bf(i, j) + z * cc(i, j)); ENDDO; ENDDO
  	DO i = 1 , nc - 1
    	DO j = i + 1 , nc
      		v = 0.
      		DO k = 1 , na; v = v + a(k) * sc(i, j, k); ENDDO
     		u = u + c(i) * c(j) * (pf(i, j) + v)
    	ENDDO
  	ENDDO
  	DO i = 1 , na - 1
    	DO j = i + 1 , na
      		v = 0.
      		DO k = 1 , nc; v = v + c(k) * sa(i, j, k); ENDDO
     		 u = u + a(i) * a(j) * (qf(i, j) + v)
    	ENDDO
  	ENDDO
  	DO k = 1 , nn
    	DO i = 1 , nc
      		u = u + h(k) * c(i) * lc(k, i)
    	ENDDO
  	ENDDO
  	DO k = 1 , nn
    	DO j = 1 , na
      		u = u + h(k) * a(j) * la(k, j)
    	ENDDO
  	ENDDO
  	DO k = 1 , nn
    	DO i = 1 , nc
      		DO j = 1 , na
        		u = u + h(k) * c(i) * a(j) * xi(k, i, j)
      		ENDDO
    	ENDDO
  	ENDDO
  	s = 0.
  	DO i = 1 , nc; s = s + c(i); ENDDO
  	DO j = 1 , na; s = s + a(j); ENDDO
  	co = 1. + 2. * u / s; aw = EXP(-s * co / mh2o)
  	gact(1) = gc(1); gact(2) = gc(2); gact(3) = gc(3)
  	gact(4) = gc(4); gact(22) = gc(5); gact(5) = gc(6)
  	gact(18) = gc(7); gact(23) = gc(8); gact(11) = gc(9)
  	gact(6) = ga(1); gact(7) = ga(2); gact(25) = ga(3)
  	gact(12) = ga(4); gact(14) = ga(5); gact(13) = ga(6)
  	gact(24) = ga(7); gact(19) = ga(8); gact(20) = ga(9)
  	gact(21) = ga(10); gact(8) = ga(11)
  	gact(10) = aw*aw*gn(1)**log(10.D0); gact(9) = gn(2); gact(15) = gn(3)
  	gact(16) = 1.; gact(17) = 1.
  	ndepact = 1
  
  	deallocate (ec, fc, xc, ea, fa, xa, stat=alloc_res)
	deallocate (pp, p, pf, qp, q, qf, stat=alloc_res)
	deallocate (cc, bf, b, bp, stat=alloc_res)
	deallocate (gc, ga, gn, stat=alloc_res)
  
END SUBROUTINE actp
