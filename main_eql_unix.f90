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


program eql

!EQL
!System: Na-K-Li-Ca-Mg-Cl-SO4-NO3-CO2-HCO3-CO3-H4SiO4-B(OH)3-H-OH-H2O
!Initial equilibrium for the simulation of evaporation EVP
!Ionic interaction (Pitzer)
!Solution of the system of equations by Newton-Raphson method
!Temperature dependence included for minerals
!  and for some Pitzer's coefficients (0-50C)
!Options to change PCO2 and to dilute the initial solution
!Connected files: MATRICE1-MURTF2-MURTF3-COMPLEX3-COEFFT4-DENSITE
!Program chained to EVP through transfer file: *.tra
!Last revision: march 1999
!author: F. RISACHER
!conversion to FORTRAN 90: A. CLEMENT
!
!	use dfwin ! WIN32 only
!       use dfport  ! win32 only
	use data_common
	implicit none


	integer, parameter	:: n = 25, ntot = 12, n0 = 10, max_cat=5, max_an=5
	real(kind=8), allocatable, dimension (:)	:: tot, tot0, totinit, molal, act, gact, atom
	real(kind=8), allocatable, dimension (:)	:: zz, xx, psc,cat,ani
	real(kind=8), allocatable, dimension (:)	:: mu, mum, psol, pai
	real(kind=8), allocatable, dimension (:,:)	:: z, wmin
	integer, allocatable, dimension (:,:)	:: kmat
	integer, allocatable, dimension (:)	:: nchcat, nchani, nch, ica, lmin, nwmin, kinvar
	character (len=maxlen), allocatable, dimension (:)	:: constit_S, aq_S, mineral_s, mineral0_S
	character (len=maxlen), allocatable, dimension (:)	:: nom_ion_S
	real(kind=8), allocatable, dimension (:)	:: c_ion

	integer	:: i,j,k,l,ii,nminer,nm,nt,nc,na,ncomp,icat,iani,ialc,nconstit,alloc_res,indice,result_evap
	integer	:: nm0,nu,ni,nj,nconv,num,nwmp,nwm,kinvariant,nw,nmneuf,npasfich,npasimp,ncompt,npasecran
	character (len=maxlen)	:: num_S,unit_S,x_S,e_S,repinit_S,min_S,pco2_S,pco_S,syst_S,rep_S,min_col1,min_col2
	character (len=maxlen)	:: zone_S,m_S,p_S,y_S,transf_S,transfneuf_S,fichmin_S,molemin_S,fich_S,prnt_S,constituant_S
	character (len=maxlen)	:: at_S,bt_S,ct_S,dt_S,et_S
	real(kind=8)	:: temp,dens,ph,at, bt, ct, dt, et,pk0,pkf,pkstd,pkmol,pkeq,dph,dil,diltot,mh2o
	real(kind=8)	:: a,b,c,u,sc,cmax,sa,amax,dca,delta,xu,eq,s,ee,stdi,aw,fi,std,ef,tinit,stdmax
	real(kind=8)	:: po,po0,poinit,phinit,ph0,pc,poa,alcar,albor,alsil,aloh,alh,altest,d,v,xi,eps
	integer, external	:: system !unix only
	character (len=30)	:: d_and_t_S
	integer	:: hconsole						!WIN32 only
	character(len=100)	:: lpszFileName		!WIN32 only


	interface
	 function f_lcase(c) result (minus)
	  character(len=*) :: c 
	  character (len=1000) :: minus
	 end function f_lcase
	end interface

	interface
	 function f_ucase(c) result (minus)
		character(len=*) :: c 
		character (len=1000) :: minus
	 end function f_ucase
	end interface

! From here WIN32 only
!	lpszFileName="CONOUT$"C
!	hConsole = CreateFile(lpszFileName, IOR(GENERIC_WRITE , GENERIC_READ),                   &
!                         IOR(FILE_SHARE_READ , FILE_SHARE_WRITE),                           &
!                         NULL_SECURITY_ATTRIBUTES, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0)
!	if( hConsole == -1 ) then
!      i =  GetLastError ()
!      if (i .eq. ERROR_ALREADY_EXISTS) then
!        i = 1
!      end if
!      write(*,*)
!      write (*,*) ("Error: Unable to open console.")
!      write(*,*)
!      stop 1
!	end if
!until there Win32 only

	call OS_VERIF ! what's the operating system?
	call color(hconsole,c_lred)
	write(*,*)
	write(*,*) "This is EQL.............."
	write(*,*)
	call color(hconsole,c_white)
	call fdate(d_and_t_S)
	write(*,*) d_and_t_S
	write(*,*)



	repinit_S = "."; fich_S="."; molemin_S="."; fichmin_S="."
	pk = .1D0; pk0 = pk; pkf = .0001D0; pkstd = .1D0
	dph = .2D0; dil = 0.; diltot = 1.
	mh2o = 55.51D0; pco2_S=""
	nconv=2; eps=1.D-12

	allocate (constit_S(0:max_constit), psc(0:14), stat=alloc_res)
	allocate (tot(0:ntot), tot0(0:ntot), totinit(0:n0), stat=alloc_res)
	allocate (nch(0:n), molal(0:n), act(0:n), gact(0:n), aq_S(0:n), atom(0:n), stat=alloc_res)
	allocate (kmat(0:n, 0:n), ica(0:n), stat=alloc_res)
	allocate (z(0:n, 0:n), zz(0:n), xx(0:n), stat=alloc_res)
	allocate (cat(0:max_cat), ani(0:max_an), nchcat(0:max_cat), nchani(0:max_an), stat=alloc_res)
	psc=0.; tot=0.; tot0=0.; totinit=0.; molal=0.; act=0.; gact=0.; atom=0.; z=0.; zz=0.; xx=0.; cat=0.; ani=0.
	nch=0; kmat=0; ica=0; nchcat=0; nchani=0 

	open(10,file="aqu.dat",action="READ")
	DO i = 0 , n
	    READ(10,*) aq_S(i), atom(i), nch(i)
	ENDDO
	close (10)

	DO i = 1 , n;  ica(i) = 1; ENDDO

	OPEN (10,file="matrice1",action="READ")
	    DO i = 1 , n
			read(10,*)  (kmat(i, j),j=1,n)
	    ENDDO
	CLOSE (10)
	kmat(13, 0) = -1; kmat(15, 0) = -1; kmat(20, 0) = 3; kmat(21, 0) = 5

	num_S = ""; temp = 25D0; dens = 1.D0; ph = 7.D0
	tot(1) = 0.; tot(2) = 0.; tot(3) = 0.
	tot(4) = 0.; tot(5) = 0.
	tot(6) = 0.; tot(7) = 0.; tot(8) = 0.
	tot(9) = 0.; tot(10) = 0.
	tot(11) = 10 ** (-ph); tot(12) = 0.

	call input_data ! read the user's data

	prnt_S='.' // DIR_SEP_S // num_S(1:len_trim(num_S)) // ".log"


	nminer = 1
	DO i = 1 , 10
	    IF (tot(i) > 0.)  nminer = nminer + 1
	ENDDO
	IF (dens == 0.)  dens = 1.D0
	IF (dens == 1.D0) THEN 
		unit_S = "molal" 
	ELSE 
		unit_S = "molar"
	END IF

	OPEN (10,file="complex3",action="READ")
	    DO i = 1 , 14
	        read(10,*)  x_S, at, bt, ct, dt, et
	        psc(i) = 10 ** (at + bt / 300D0 * temp + ct / 30000D0 * temp ** 2 + &
	                 dt / 3000000D0 * temp ** 3 + et / 300000000D0 * temp ** 4)
	    ENDDO
	CLOSE (10)

	OPEN (10,file="murtf2",action="READ")
	    read(10,*)  nc, na, nm
	CLOSE (10)
	nt = nc + na

	allocate (nom_ion_S(0:nt), c_ion(0:nt), stat=alloc_res)
	allocate (wmin(0:nm, 0:nt), lmin(0:nm), nwmin(0:nm), stat=alloc_res)
	allocate (mineral_S(0:nm), mineral0_S(0:nm), mu(0:nt), mum(0:nm), stat=alloc_res)
	allocate (psol(0:nm), pai(0:nm), kinvar(0:nt), stat=alloc_res)

	wmin=0.; mu=0.; mum=0.; psol=0.; pai=0.
	lmin=0; nwmin=0; kinvar=0

	OPEN (10,file="murtf2",action="READ")
	    read(10,*)  a, b, c
	    DO i = 0 , 14
	        read(10,*)  x_S, at, bt, ct, dt, et
	        x_S =f_lcase(x_S)
	        DO j = 0 , nt
	            IF (x_S == aq_S(j)) THEN
	                mu(j) = at + bt / 300.D0 * temp + ct / 30000.D0 * temp ** 2 + &
	                        dt / 3000000.D0 * temp ** 3 + et / 300000000.D0 * temp ** 4
	            END IF
	        ENDDO
	    ENDDO
	    DO k = 1 , nm
	        read(10,*)  mineral_S(k),ncomp, (c_ion(i),nom_ion_S(i),i=1,ncomp), at, bt, ct, dt, et
	        DO i = 1 , ncomp
	            x_S = f_lcase(nom_ion_S(i))
	            DO j = 0 , nt
	                IF (x_S == aq_S(j))  wmin(k, j) = c_ion(i)
	            ENDDO
	        ENDDO
	        mum(k) = at + bt / 300.D0 * temp + ct / 30000.D0 * temp ** 2 + &
	                 dt / 3000000.D0 * temp ** 3 + et / 300000000.D0 * temp ** 4
	    ENDDO
	CLOSE (10)

	DO k = 1 , nm
	    u = mum(k)
	    DO i = 0 , nt
	        u = u - wmin(k, i) * mu(i)
	    ENDDO
	    psol(k) = EXP(u)
	ENDDO

	sc = 0.; cmax = 0.
	DO i = 1 , ntot
	    IF (nch(i) > 0 .AND. i /= 11) THEN
	        sc = sc + tot(i) * nch(i)
	        IF (tot(i) * nch(i) >= cmax) THEN
	            cmax = tot(i) * nch(i)
	            icat = i
	        END IF
	    END IF
	ENDDO
	sa = 0.; amax = 0.
	DO i = 1 , ntot
	    IF (nch(i) < 0) THEN
	        sa = sa + tot(i) * (-nch(i))
	        IF (tot(i) * (-nch(i)) >= amax) THEN
	            amax = tot(i) * (-nch(i))
	            iani = i
	        END IF
	    END IF
	ENDDO
	IF (sc+sa /= 0.) dca = 200. * ABS(sc - sa) / (sc + sa)
	delta = sc - sa
	write(*,*) "sum of cations = ", real(sc)
	write(*,*) "sum of anions  = ", real(sa)
	write(*,*) "Electrical balance = ", real(dca * 100 + .5) / 100, "%"
	tot(icat) = tot(icat) - delta / 2. / nch(icat)
	tot(iani) = tot(iani) + delta / 2. / (-nch(iani))
	DO i = 1 , 12
	    tot0(i) = tot(i)
	ENDDO

	OPEN (10,file="murtf3",action="READ")
	    read(10,*)  nc, na, nm0
	    DO i = 1 , nc + na + 1
	        read(10,*)  x_S
	    ENDDO
	    DO k = 1 , nm0
	        read(10,*)  mineral0_S(k)
	    ENDDO
	CLOSE (10)
	DO k = 1 , nm
	    nwmin(k) = 0
	ENDDO
	DO k = 1 , nm
	    DO l = 1 , nm0
	        IF (mineral0_S(l) == mineral_S(k))  nwmin(k) = 1
	    ENDDO
	ENDDO
	min_S = "murtf3"


	500 continue !dilution
	molal(1) = tot(1); molal(2) = tot(2); molal(3) = tot(3)
	molal(6) = tot(6); molal(8) = tot(8)
	tot(11) = 10. ** (-ph)
	molal(11) = tot(11)
	molal(13) = psc(1) / molal(11)
	IF (tot(9) > 0.) THEN
	    a = 1. + molal(13) / psc(7)
	    b = 3. * psc(8) * molal(13)
	    c = 4. * psc(9) * molal(13) ** 2
	    xu = tot(9) / 2.
	    u = xu
	    DO
	        eq = a * xu + b * xu ** 3 + c * xu ** 4
	        IF (200. * ABS(eq - tot(9)) / (eq + tot(9)) < pk)  EXIT 
	        u = u / 2.
	        IF (eq > tot(9)) THEN 
	        	xu = xu - u 
	        ELSE 
	        	xu = xu + u
	        END IF
	    ENDDO
	    molal(9) = xu
	    molal(19) = molal(9) * molal(13) / psc(7)
	    molal(20) = molal(13) * molal(9) ** 3 * psc(8)
	    molal(21) = molal(13) ** 2 * molal(9) ** 4 * psc(9)
	END IF
	molal(14) = (tot(12) + molal(11) - molal(13) - molal(19) - molal(20) - 2. * molal(21)) / (2. + molal(11) / psc(2))
	molal(12) = (tot(12) + molal(11) - molal(13) - molal(19) - molal(20) - 2. * molal(21)) / (1. + 2. * psc(2) / molal(11))
	molal(15) = molal(12) * molal(11) / psc(3)
	molal(4) = tot(4) / (1. + molal(14) / psc(4) + molal(19) / psc(10))
	molal(16) = molal(4) * molal(14) / psc(4)
	molal(22) = molal(4) * molal(19) / psc(10)
	molal(5) = tot(5) / (1. + molal(14) / psc(5) + molal(13) / psc(6) + molal(19) / psc(11))
	molal(17) = molal(5) * molal(14) / psc(5)
	molal(18) = molal(5) * molal(13) / psc(6)
	molal(23) = molal(5) * molal(19) / psc(11)
	molal(10) = tot(10) / (1. + psc(12) / molal(11))
	molal(24) = tot(10) / (1. + molal(11) / psc(12))
	molal(7) = tot(7) / (1. + molal(11) / psc(13))
	molal(25) = molal(7) * molal(11) / psc(13)

	sc = 0.; cmax = 0.
	DO i = 1 , n
	    IF (nch(i) > 0) THEN
	        sc = sc + molal(i) * nch(i)
	        IF (molal(i) * nch(i) > cmax) THEN
	            cmax = molal(i) * nch(i)
	            icat = i
	        END IF
	    END IF
	ENDDO
	sa = 0.; amax = 0.
	DO i = 1 , n
	    IF (nch(i) < 0) THEN
	        sa = sa + molal(i) * (-nch(i))
	        IF (molal(i) * (-nch(i)) > amax) THEN
	            amax = molal(i) * (-nch(i))
	            iani = i
	        END IF
	    END IF
	ENDDO
	delta = sc - sa
	molal(icat) = molal(icat) - delta / 2. / nch(icat)
	molal(iani) = molal(iani) + delta / 2. / (-nch(iani))

	sc = molal(1) + molal(2) + molal(3) + molal(4) * 2. + molal(5) * 2. + molal(11) + molal(18) + molal(22) + molal(23)
	sa = molal(6) + molal(7) * 2. + molal(8) + molal(12) + molal(13) + molal(14) * 2. &
	     + molal(19) + molal(20) + molal(21) * 2. + molal(24) + molal(25)

	write(*,*)
	write(*,*) "Sum of cations = ", real(sc), "  corrected for  ", aq_S(icat)(1:len_trim(aq_S(icat)))
	write(*,*) "Sum of anions  = ", real(sa), "  corrected for  ", aq_S(iani)(1:len_trim(aq_S(iani)))
	write(*,*)

	s = 0.
	DO i = 1 , n
	    s = s + molal(i) * atom(i)
	ENDDO
	IF (unit_S == "molar") THEN
	    ee = 1000. / (1000. * dens - s)
	    DO i = 1 , ntot
	        IF (i /= 11)  tot(i) = tot(i) * ee
	    ENDDO
	    DO i = 1 , n
	        IF (i /= 11) molal(i) = molal(i) * ee
	    ENDDO
	ELSEIF (unit_S == "molal") THEN
	    ee = 1.
	END IF
	stdi = s * ee

	200 continue !reprise
	DO
		nu=1
		ncompt=0
	    DO while(nu/=0)
	        CALL actp(molal, gact, aw, fi, temp)
	        DO i = 1 , n
	            act(i) = molal(i) * gact(i)
	        ENDDO
	        act(0) = aw
	        tot(11) = (10 ** (-ph)) / gact(11)
	        act(11) = 10 ** (-ph)
	        DO i = 1 , 12
	            DO j = 1 , n
	                IF (molal(j) /= 0.)  z(i, j) = kmat(i, j)
	            ENDDO
	            u = 0.
	            DO j = 1 , n
	                u = u + kmat(i, j) * molal(j)
	            ENDDO
	            zz(i) = tot(i) - u
	        ENDDO
      
	        DO i = 13 , n
	            DO j = 1 , n
	                IF (molal(j) /= 0.) THEN
	                    z(i, j) = kmat(i, j) / molal(j)
	                ELSEIF (molal(j) == 0.) THEN
	                    z(i, j) = 0.
	                END IF
	            ENDDO
	            u = 0.
	            DO j = 0 , n
	                IF (act(j) > 0.)  u = u + kmat(i, j) * LOG(act(j))
	            ENDDO
	            zz(i) = LOG(psc(i - 12)) - u
	        ENDDO

	        DO k = 1 , n0
	            IF (tot(k) == 0. .AND. k /= 12) THEN
	                ica(k) = 0
	                DO i = k + 1 , n
	                    IF (kmat(i, k) /= 0)  ica(i) = 0
	                ENDDO
	                DO j = k + 1 , n
	                    IF (kmat(k, j) /= 0)  ica(j) = 0
	                ENDDO
	            END IF
	        ENDDO

	        ni = n; nj = n
	        DO k = n , 1 , -1
	            IF (ica(k) == 0) THEN
	                DO i = k , ni - 1
	                    DO j = 1 , nj
	                        z(i, j) = z(i + 1, j)
	                    ENDDO
	                    zz(i) = zz(i + 1)
	                ENDDO
	                ni = ni - 1
	                DO j = k , nj - 1
	                    DO i = 1 , ni
	                        z(i, j) = z(i, j + 1)
	                    ENDDO
	                ENDDO
	                nj = nj - 1
	            END IF
	        ENDDO

	        DO k = 2 , ni
	            DO i = k , ni
	                IF (z(i, k - 1) /= 0.) THEN
	                    u = z(k - 1, k - 1) / z(i, k - 1)
	                    DO j = k , ni
	                        z(i, j) = z(k - 1, j) - z(i, j) * u
	                    ENDDO
	                    zz(i) = zz(k - 1) - zz(i) * u
	                END IF
	            ENDDO
	        ENDDO

	        xx(ni) = zz(ni) / z(ni, ni)
	        DO i = ni - 1 , 1 , -1
	            s = 0.
	            DO j = i + 1 , ni
	                s = s + z(i, j) * xx(j)
	            ENDDO
	            xx(i) = (zz(i) - s) / z(i, i)
	        ENDDO

	        DO k = 1 , n
	            IF (ica(k) == 0) THEN
	                DO i = ni , k , -1
	                    xx(i + 1) = xx(i)
	                ENDDO
	                xx(k) = 0.
	                ni = ni + 1
	            END IF
	        ENDDO	        
	        ncompt = ncompt + 1
	        write(*,*) "iteration molalities ",ncompt
	        IF (ncompt >= 100) THEN
	           DO i= 1, n
	             IF (molal(i) + xx(i) / nconv < 0.) THEN
	                write(*,*) "the equation set diverges: end of program"
	                stop
	             END IF
	           ENDDO
	        END IF
	        DO i= 1, n
	           IF (molal(i) + xx(i) / nconv < 0.) THEN
	              molal(i)=eps
	           ELSE
	              molal(i) = molal(i) + xx(i) / nconv
	           END IF
	        ENDDO
	        
	        nu = 0	        
	        DO i = 1 , n
	            IF (ica(i) == 1) THEN
	                IF (200. * ABS(xx(i) / nconv / (2. * molal(i) - xx(i) / nconv)) > pk)  nu = 1
	            END IF
	        ENDDO
	    ENDDO
  
	    std = 0.
	    DO i = 0 , n
	        std = std + molal(i) * atom(i)
	    ENDDO
	    call color(hconsole,c_cyan)
	    write(*,*) "tdsi = ", real(stdi)
	    write(*,*) "tds  = ", real(std)
	    call color(hconsole,c_white)
	    IF (ABS(std - stdi) / (std + stdi) * 200. < pkstd) THEN
	        EXIT 
	    ELSE
	        IF (unit_S == "molar") THEN
	            ef = (1000. + std) / dens / 1000.
	            DO i = 1 , ntot
	                IF (i /= 11)  tot(i) = tot(i) / ee * ef
	            ENDDO
	            DO i = 0 , n
	                IF (i /= 11)  molal(i) = molal(i) / ee * ef
	            ENDDO
	            ee = ef
	        END IF
	        call color(hconsole,c_lgreen); write(*,*) "iteration TDS"; call color(hconsole,c_white)
	        stdi = std
	    END IF
	ENDDO


	IF (unit_S == "molal" .AND. dil == 0.) THEN
	    DO i = 1 , 5; cat(i) = tot(i); nchcat(i) = nch(i); ENDDO
	    DO j = 1 , 3; ani(j) = tot(j + 5); nchani(j) = -nch(j + 5); ENDDO
	    ani(4) = molal(12); ani(5) = molal(14) + molal(16) + molal(17)
	    nchani(4) = -nch(12); nchani(5) = -nch(14)
	    CALL densite(cat, ani, nchcat, nchani, unit_S, dens)
	END IF

	po = LOG(act(15) / psc(14)) / LOG(10D0)
	IF (pco2_S == "") THEN
	    IF (diltot == 1.) THEN
	        poinit = po
	        phinit = ph
	    END IF
	    call color(hconsole,c_lred)
	    write(*,*) "LOG PCO2 = ", real(po)
	    call display ("Other LOG PCO2 (else press <enter>): ", pco_S)
	    IF (pco_S /= "") THEN
	    	read(pco_S,*) pc
	        po0 = po; ph0 = ph
	        pco2_S = "y"
	    ELSEIF (pco_S == "") THEN
	        pco2_S = "n"
	    END IF
	    call color(hconsole,c_white)
	END IF
	IF (pco2_S == "y") THEN
	    write(*,*)
	    write(*,*) "Log(PCO2) selected   = ", real(pc)
	    write(*,*) "Log(PCO2) calculated = ", real(po)
	    write(*,*)
	    IF (ABS(po - pc) > .01) THEN
	        IF (po < pc .AND. poa < pc)  ph = ph - dph
	        IF (po < pc .AND. poa > pc)  THEN
	        	dph = dph / 2; ph = ph - dph
	        END IF
	        IF (po > pc .AND. poa > pc)  ph = ph + dph
	        IF (po > pc .AND. poa < pc) THEN 
	        	dph = dph / 2; ph = ph + dph
	        END IF 
	        poa = po
	        GOTO 200 !reprise
	    END IF
	END IF
	IF (pk > pkf) THEN
	    pk = pkf
	    call color(hconsole,c_bwhite)
	    write(*,*) "last iteration"
	    call color(hconsole,c_white)
	    GOTO 200 !reprise
	END IF

	400 continue  !ecran
	DO i = 1 , n0
	    totinit(i) = 0.
	    DO j = 1 , n
	        totinit(i) = totinit(i) + molal(j) * kmat(i, j)
	    ENDDO
	ENDDO
	write(*,*)
	call color(hconsole,c_yellow)
	write(*,*) num_S(1:len_trim(num_S))
	call color(hconsole,c_bwhite)
	write(*,'(t12,a,t28,a,t44,a,t60,a)') "  MOLALITY", "ACT COEFF", "  ACTIVITY", " MOLAL TOT"
	call color(hconsole,c_white)
	write(*,*)
	DO i = 1 , n0
	    IF (molal(i) /= 0.) THEN
	        write(*,'(a,t12,g13.6,t28,g13.6,t44,g13.6,t60,g13.6)') &
	                aq_S(i)(1:len_trim(aq_S(i))),real(molal(i)),real(gact(i)),real(act(i)),real(totinit(i))
	    END IF
	ENDDO
	DO i = n0 + 1 , n
	    IF (molal(i) /= 0.) THEN
	        write(*,'(a,t12,g13.6,t28,g13.6,t44,g13.6)') aq_S(i)(1:len_trim(aq_S(i))),real(molal(i)),real(gact(i)),real(act(i))
	    END IF
	ENDDO
	write(*,*)

	alcar = molal(12) + 2. * (molal(14) + molal(16) + molal(17))
	albor = molal(19) + molal(20) + 2. * molal(21) + molal(22) + molal(23)
	alsil = molal(24)
	aloh = molal(13) + molal(18)
	alh = -molal(11) - molal(25)
	altest = alcar + albor + alsil + aloh + alh

	write(*,*) "ELECTRICAL BALANCE     = ", real(dca), "% corrected on  ", &
	           aq_S(icat)(1:len_trim(aq_S(icat))), " and ", aq_S(iani)(1:len_trim(aq_S(iani)))
	write(*,*) "TOTAL DISSOLVED SOLIDS = ", real(std), "g/kg(H2O)"
	write(*,*) "MOLAL/MOLAR FACTOR     = ", real(1000. * dens / (1000. + std))
	write(*,*) "DENSITY                = ", real(dens)
	write(*,*) "IONIC STRENGTH         = ", real(fi)
	write(*,*) "WATER ACTIVITY         = ", real(aw)
	IF (diltot > 1.) THEN
		write(*,*) "DILUTION               = ", real(diltot)
	END IF
	IF (pco2_S == "" .OR. pco2_S == "n") THEN
	    write(*,*) "pH                     = ", real(ph)
	    write(*,*) "LOG PCO2               = ", real(po)
	ELSEIF (pco2_S == "y") THEN
	    write(*,*) "INITIAL LOG PCO2       = ", real(poinit)
	    write(*,*) "INITIAL pH             = ", real(phinit)
	    write(*,*) "CALCULATED LOG PCO2    = ", real(po)
	    write(*,*) "CALCULATED pH          = ", real(ph)
	END IF
	write(*,*) "CARBONATE ALKALINITY   = ", real(alcar)
	write(*,*) "BORATE ALKALINITY      = ", real(albor)
	write(*,*) "SILICATE ALKALINITY    = ", real(alsil)
	write(*,*) "OH ALKALINITY          = ", real(aloh)
	write(*,*) "H ALKALINITY           = ", real(alh)
	write(*,'(a,g13.6,t50,a,g13.6)') " TOTAL ALKALINITY       = ", real(altest), "init. alk. = ", real(tot(12))
	call display(" Press <enter> to continue:", x_S)

	write(*,*)
	call color(hconsole,c_lcyan)
	write(*,'(a,t15,a,t32,a,t50,a)') "TESTS", "SUM OF SPECIES"," INIT. CONC.","BALANCE %"
	write(*,*)
	call color(hconsole,c_white)
	DO i = 1 , 12
	    IF (ica(i) == 1) THEN
	        u = 0.
	        DO j = 1 , n
	            u = u + molal(j) * kmat(i, j)
	        ENDDO
	        IF (u + tot(i) /= 0.)  d = 200. * ABS(u - tot(i)) / (u + tot(i))
	        IF (i == 12) THEN
	            IF (u + tot(i) /= 0. .AND. tot(i) /= 0.) THEN
	            	write(*,'(a,t15,g13.6,t32,g13.6,t50,g13.6)') "ALK",real(u),real(tot(i)),real(d)
	            ELSE 
	            	write(*,'(a,t15,g13.6,t32,g13.6)') "ALK",real(u),real(tot(i))
	            END IF
	        ELSE
	            zone_S=f_ucase(aq_S(i))
	            write(*,'(a,t15,g13.6,t32,g13.6,t50,g13.6)') zone_S(1:len_trim(zone_S)),real(u),real(tot(i)),real(d)
	        END IF
	    END IF
	ENDDO
	write(*,*) ; write(*,*)
	call color(hconsole,c_lcyan)
	write(*,'(t15,a,t32,a,t50,a)') " LOG(IAP)","  LOG(K)","  BALANCE %"
	write(*,*)
	call color(hconsole,c_white)
	DO i = 13 , n
	    u = 0.
	    IF (ica(i) == 1) THEN
	        DO j = 0 , n
	            IF (act(j) /= 0.)  u = u + kmat(i, j) * LOG(act(j)) / LOG(10D0)
	        ENDDO
	        v = LOG(psc(i - 12)) / LOG(10D0)
	        d = 200. * ABS(u - v) / (u + v)
	        IF (i == 13) THEN
	           write(*,'(a,t15,g13.6,t32,g13.6,t50,g13.6)') "H2O",real(u),real(v),real(d)
	        ELSE
	            zone_S=f_ucase(aq_S(i))
	            write(*,'(a,t15,g13.6,t32,g13.6,t50,g13.6)') zone_S(1:len_trim(zone_S)),real(u), real(v),real(d)
	        END IF
	    END IF
	ENDDO
	call display(" Press <enter> to continue:", x_S)

	write(*,*)
	call color(hconsole,c_lmagenta)
	write(*,'(t20,a,t40,a,t60,a)') " SOLUB. PROD.","ION. ACT. PROD."," SATUR. RATIO"
	call color(hconsole,c_white)
	write(*,*)
	nwm = 0; nwmp = 0
	DO k = 1 , nm
	    pai(k) = 1.
	    DO i = 0 , nt
	        pai(k) = pai(k) * act(i) ** wmin(k, i)
	    ENDDO
	    IF (pai(k) /= 0.) THEN
	        IF (pai(k) / psol(k) >= 1.) call color(hconsole,c_lred)
	        IF (pai(k) / psol(k) >= .5 .AND. pai(k) / psol(k) < 1.) call color(hconsole,c_green)
	        IF (pai(k) / psol(k) >= .1 .AND. pai(k) / psol(k) < .5) call color(hconsole,c_cyan)
	        IF (nwmin(k) == 0)  zone_S=" " // mineral_S(k)(1:15)
	        IF (nwmin(k) == 1)  zone_S="*" // mineral_S(k)(1:15)
	        x_S=" "
	        IF (pai(k) / psol(k) >= 1. .AND. nwmin(k) == 1) THEN
	            nwm = nwm + 1
	            x_S="*"
	        ELSEIF (pai(k) / psol(k) >= .9 .AND. pai(k) / psol(k) < 1. .AND. nwmin(k) == 1) THEN
	            nwmp = nwmp + 1
	        END IF
	        write(*,'(a,t20,g13.6,t40,g13.6,t60,g13.6,t78,a)') &
	              zone_S(1:len_trim(zone_S)),real(psol(k)),real(pai(k)),real(pai(k) / psol(k)),x_S(1:1)      
	        call color(hconsole,c_white)
	    END IF
	    IF (k == 42)  call display(" Press <enter> to continue:", x_S)
	ENDDO

	write(*,*)
	call display("Print the initial solution in a file ? (y/n=<enter>): ", x_S)
	IF (x_S == "y") THEN
		open(9,file=prnt_S)  ! log file  i.e. old printer file
	    write(*,*) "LOG FILE IS ", prnt_S(1:len_trim(prnt_S))
	    write(9,*) num_S(1:len_trim(num_S)); write(9,*)
	    write(9,'(t15,"  MOLALITY",t35,"ACT. COEFF.",t55,"  ACTIVITY")')
	    write(9,*)
	    DO i = 1 , n
	        IF (molal(i) /= 0.) THEN
	            write(9,'(a,t15,g13.6,t35,g13.6,t55,g13.6)') &
	                    aq_S(i)(1:len_trim(aq_S(i))),real(molal(i)),real(gact(i)),real(act(i))
	        END IF
	    ENDDO
	    write(9,*)
	    write(9,*) "ELECTRICAL BALANCE     = ", real(dca), "% corrected on ", &
	               aq_S(icat)(1:len_trim(aq_S(icat))), " and ", aq_S(iani)(1:len_trim(aq_S(iani)))
	    write(9,*) "TOTAL DISSOLVED SOLIDS = ", real(std), "g/kg(H2O)"
	    write(9,*) "MOLAL/MOLAR FACTOR     = ", real(1000. * dens / (1000. + std))
	    write(9,*) "DENSITY                = ", real(dens)
	    write(9,*) "IONIC STRENGTH         = ", real(fi)
	    write(9,*) "WATER ACTIVITY         = ", real(aw)
	    IF (diltot > 1.) THEN
			write(9,*) "DILUTION               = ", real(diltot)
	    END IF
	    IF (pco2_S == "" .OR. pco2_S == "n") THEN
	        write(9,*) "pH                     = ", real(ph)
	        write(9,*) "LOG PCO2               = ", real(po)
	    ELSEIF (pco2_S == "y") THEN
	        write(9,*) "INITIAL LOG PCO2       = ", real(poinit)
	        write(9,*) "INITIAL pH             = ", real(phinit)
	        write(9,*) "CALCULATED LOG PCO2    = ", real(po)
	        write(9,*) "CALCULATED pH          = ", real(ph)
	    END IF
	    write(9,*) "CARBONATE ALKALINITY   = ", real(alcar)
	    write(9,*) "BORATE ALKALINITY      = ", real(albor)
	    write(9,*) "SILICATE ALKALINITY    = ", real(alsil)
	    write(9,*) "OH ALKALINITY          = ", real(aloh)
	    write(9,*) "H ALKALINITY           = ", real(alh)
	    write(9,'(a,g13.6,t50,a,g13.6)') " TOTAL ALKALINITY       = ", real(altest), "init. alk. = ", real(tot(12))
	    write(9,*)
	    write(9,'(t20," SOLUB. PROD.",t40,"ION. ACT. PROD.",t60," SATUR. RATIO")') 
	    write(9,*)
	    DO k = 1 , nm
	        IF (pai(k)>0.) THEN
				write(9,'(a,t20,g13.6,t40,g13.6,t60,g13.6)') mineral_S(k)(1:16),real(psol(k)), real(pai(k)),real(pai(k)/psol(k))
			ENDIF
		ENDDO
		close(9)
	END IF

	DO k = 1 , nm
	    lmin(k) = 0
	ENDDO
	write(*,*)
	write(*,*) "The initial solution is oversaturated in ", nwm, " mineral(s) of the data base MURTF3:"
	call color(hconsole,c_bwhite)
	DO k = 1 , nm
	    IF (pai(k) / psol(k) >= 1. .AND. nwmin(k) == 1) THEN
	        write(*,'(a,t20,g13.6)') mineral_S(k)(1:len_trim(mineral_S(k))),real(pai(k) / psol(k))
	        lmin(k) = 1
	    END IF
	ENDDO
	call color(hconsole,c_white)
	IF (nwm > nminer) THEN
	    write(*,*)
	    call color(hconsole,c_lred)
	    write(*,*) "VIOLATION OF THE PHASE RULE :"
	    write(*,*) "The maximum number of minerals allowed is: ", nminer
	    write(*,*) "The evaporation program cannot start with this paragenesis"
	    call color(hconsole,c_white)
	ELSEIF (nwm <= nminer) THEN
	    write(*,*)
	    CALL invar(nm, nt, wmin, mineral_S, lmin, psol, psc, kinvariant, kinvar)
	    IF (kinvariant > 0) THEN
	        write(*,*) "The activity of water is constrained by : "
	        call color(hconsole,c_red)
	        DO k = 1 , kinvariant
	            write(*,'(a," ",$)') mineral_S(kinvar(k))(1:len_trim(mineral_S(kinvar(k)))) 
	        ENDDO
	        call color(hconsole,c_white)
	        write(*,*)
	    ELSEIF (kinvariant == -1) THEN
	        call color(hconsole,c_red)
	        write(*,*) "System in thermodynamic desequilibrium"
	        write(*,*) "The activity of water is constrained at different values"
	        write(*,*) "by more than one mineral assemblage"
	        call color(hconsole,c_white)
	    ELSEIF (kinvariant == -2) THEN
	        call color(hconsole,c_red)
	        write(*,*) "System in thermodynamic desequilibrium: inconsistent mineral assemblage"
	        call color(hconsole,c_white)
	    ELSEIF (kinvariant == 0) THEN
	        write(*,*) "No invariant paragenesis detected"
	    END IF
	    IF (kinvariant /= 0) THEN
	        call color(hconsole,c_lred)
	        write(*,*) "The evaporation program cannot start with this paragenesis"
	        call color(hconsole,c_white)
	    END IF

	    IF (kinvariant == 0 .AND. nwmp > 0) THEN
	        write(*,*)
	        write(*,*) "The solution is close to saturation in ", nwmp, " mineral(s) of the data base MURTF3:"
	        call color(hconsole,c_brown)
	        DO k = 1 , nm
	            IF (pai(k) / psol(k) >= .9 .AND. pai(k) / psol(k) < 1. .AND. nwmin(k) == 1) THEN
	                write(*,'(a,t20,g13.6)') mineral_S(k)(1:len_trim(mineral_S(k))),real(pai(k) / psol(k))
	                lmin(k) = 1
	            END IF
	        ENDDO
	        call color(hconsole,c_white)

	        CALL invar(nm, nt, wmin, mineral_S, lmin, psol, psc, kinvariant, kinvar)
	        write(*,*)
	        IF (kinvariant > 0) THEN
	            write(*,*) "At the start of evaporation, the activity of water may be constrained by : "
	            call color(hconsole,c_lred)
	            DO k = 1 , kinvariant
	            	write(*,'(a," ",$)') mineral_S(kinvar(k))(1:len_trim(mineral_S(kinvar(k))))
	            ENDDO
	            call color(hconsole,c_white)
	            write(*,*)
	        ELSEIF (kinvariant == -1) THEN
	            write(*,*) "System in thermodynamic desequilibrium"
	            write(*,*) "The activity of water is constrained at different values"
	            write(*,*) "by more than one mineral assemblage"
	        ELSEIF (kinvariant == -2) THEN
	            write(*,*) "System in thermodynamic desequilibrium: inconsistent mineral assemblage"
	        ELSEIF (kinvariant == 0) THEN
	            write(*,*) "No invariant paragenesis detected"
	        END IF
	        IF (kinvariant /= 0) THEN
	            write(*,*) "If the evaporation program does not start"
	            write(*,*) "dilute again slightly the solution"
	        END IF
	    END IF
	END IF


	write(*,*)
	write(*,*) "Dilute the solution ? (else press <enter>:"
	call display("  dilution = 1/", zone_S)
	dil=0.
	if (zone_S /= "") read(zone_S,*) dil
	IF (dil > 1.) THEN
	    diltot = diltot * dil
	    pco2_S = ""; pk = pk0; dph = .2
	    DO i = 1 , 12
	        IF (i /= 11)  tot(i) = tot(i) / dil
	    ENDDO
	    DO i = 1 , 5; cat(i) = tot(i); nchcat(i) = nch(i); ENDDO
	    DO j = 1 , 3; ani(j) = tot(j + 5); nchani(j) = -nch(j + 5); ENDDO
	    ani(4) = molal(12) / dil; ani(5) = (molal(14) + molal(16) + molal(17)) / dil
	    nchani(4) = -nch(12); nchani(5) = -nch(14)
	    unit_S = "molal"
	    CALL densite(cat, ani, nchcat, nchani, unit_S, dens)
	    GOTO 500 !dilution
	END IF
	IF (diltot > 1.) THEN
	    write(*,*) "The initial solution has been diluted ", diltot, " times"
	END IF
	write(*,*)
	call display("Modify the mineral data base ? (y/n=<enter>):", x_S)
	write(*,*)
	IF (x_S == "y") THEN
		write(*,*) 'Minerals already in the file EVP are marked by *'
	    nw = INT(nm / 2 + .5)
	    DO k = 1 , nw
	        IF (nwmin(k) == 1) THEN 
			    min_col1='*'//mineral_S(k)
			ELSE 
				min_col1=' '//mineral_S(k)
			END IF
	        
	        IF (k + nw <= nm) THEN
				IF (nwmin(k + nw) == 1) THEN 
					min_col2='*'//mineral_S(k + nw)
				ELSE 
					min_col2=' '//mineral_S(k + nw)
				END IF
	            write(*,'(i3,t5,a,t40,i3,t45,a)') k,min_col1(1:len_trim(min_col1)),k + nw,min_col2(1:len_trim(min_col2))
			ELSE
				write(*,'(i3,t5,a)') k,min_col1(1:len_trim(min_col1))
	        END IF
	    ENDDO
		call color(hconsole,c_white)
	    write(*,*)
	    write(*,*) "Enter mineral names (lower or upper case) or the number at the left."
	    nmneuf = 0
	    DO
	        call display(" Mineral added (<enter> to stop): ", m_S)
	        IF (m_S == "" ) EXIT 
	        m_S=adjustl(m_S)
			IF (m_S(1:1) >= '0' .and. m_S(1:1) <= '9') then
				read(m_S,*) indice
	        	m_S = mineral_S(indice)
	        END IF
	        m_S = f_ucase(m_S)
	        write(*,*) m_S(1:len_trim(m_s))
	        DO k = 1 , nm
	            IF (mineral_S(k) == m_S) THEN
	                IF (nwmin(k) == 1) THEN
	                    write(*,*) m_S(1:len_trim(m_S)), " already belongs to the mineral evaporation data base"
	                ELSEIF (nwmin(k) == 0) THEN
	                    nmneuf = nmneuf + 1
	                    nwmin(k) = 1
	                END IF
	                EXIT 
	            END IF
	        ENDDO
	        IF (k == nm + 1)  write(*,*) "Wrong name"
	    ENDDO

	    write(*,*)
	    DO
	        call display(" Mineral removed (<enter> to stop): ", m_S)
	        IF (m_S == "")  EXIT 
	        m_S=adjustl(m_S)
			IF (m_S(1:1) >= '0' .and. m_S(1:1) <= '9') then
				read(m_S,*) indice
	        	m_S = mineral_S(indice)
	        END IF
	        m_S = f_ucase(m_S)
	        write(*,*) m_S(1:len_trim(m_s))
	        DO k = 1 , nm
	            IF (mineral_S(k) == m_S) THEN
	                IF (nwmin(k) == 1) THEN
	                    nwmin(k) = 0
	                    nmneuf = nmneuf - 1
	                ELSEIF (nwmin(k) == 0) THEN
	                    write(*,*) m_S(1:len_trim(m_S))," does not belong to the mineral evaporation data base"
	                END IF
	                EXIT 
	            END IF
	        ENDDO
	        IF (k == nm + 1)  write(*,*) "wrong name"
	    ENDDO
  
	    OPEN (10,file="murtf2",action="READ")
	    OPEN (20,file="murtf0",action="WRITE")
	        read(10,*)  nc, na, a
	        write(20,*) nc, ",", na, ",", nm0 + nmneuf
	        DO i = 1 , nc + na + 1
	            read(10,'(a)')  x_S
	            write (20,'(a)') x_S
			ENDDO
			DO k = 1 , nm
				IF (nwmin(k) == 1) THEN
					read(10,'(a)') x_S
					write(20,'(a)') x_S
				ELSEIF (nwmin(k) == 0) THEN
					read(10,'(a)') x_S
				END IF
	        ENDDO
	    CLOSE (20)																
	    CLOSE (10)
	    min_S = "murtf0"
	    GOTO 400   !ecran
	END IF

	write(*,*)
	call display( "Open or closed system ? (o/c) : ", syst_S)
	IF (syst_S == "")  STOP
	write(*,*)
	write(*,*) "For an automatic increment, press <enter>"
	call display( "increment (in %) = ", x_S)
	if (x_S == "") then
		xi=0.
	else
		read(x_S,*) xi
		xi = xi / 100.
	endif
	write(*,*)
	write(*,*) "For no endpoint, press <enter> "
	call display("Endpoint at a total salinity in mg/l of: ", x_S)
	if (x_S == "") then
		stdmax=0.
	else
		read(x_S,*) stdmax
	endif
	write(*,*) 
	call display("Screen output step (1=<enter>) ?", p_S)
	IF (p_S == "") THEN
		npasecran=1
	ELSE
		read(p_S,*) npasecran
	ENDIF
	write(*,*)
	call display("Log the output screen in file ? (y/n=<enter>):", p_S)
	IF (p_S == "y")  THEN
		call display(" Output step (1=<enter>) ? ", x_S)
		p_S=prnt_S 
		IF (x_S == "") THEN
			npasimp=1
		ELSE
			read(x_S,*) npasimp
		ENDIF
	ELSE
	    npasimp=0
		p_S="n"
	ENDIF
	write(*,*)
	call display("Store the simulation in files ? (y/n=<enter>): ", y_S)
	IF (y_S == "y") THEN
	    call display("  storage step : ", npasfich)
	    call display("  unit: molarity (1) or molality (2) ?  ", k)
	    IF (k == 1) THEN 
	    	unit_S = "molar" 
	    ELSE 
	    	unit_S = "molal"
	    END IF
	    write(*,*)
	    write(*,*) "  Three files will be created:"
	    write(*,*) "    file one stores the composition of the evaporated solutions"
	    write(*,'(a,$)') "    file name is " 
	    call color(hconsole,c_lred)
	    write(*,*) num_S(1:len_trim(num_S)) // "." // "j" // syst_S(1:len_trim(syst_S)) // "&" 
	    call color(hconsole,c_white)
	    call display("    other name ? (else press <enter>) :", x_S)
		if (x_S /= "") fich_S=x_S
	    write(*,*)
	    write(*,*) "    file two stores all mineral events"
	    write(*,'(a,$)') "    file name is "
	    call color(hconsole,c_green)
	    write(*,*) num_S(1:len_trim(num_S)) // "." // "j" // syst_S(1:len_trim(syst_S)) // "@" 
	    call color(hconsole,c_white)
	    call display("    other name ? (else press <enter>) : ", x_S)
		if (x_S /= "") fichmin_S=x_S
	    write(*,*)
	    write(*,*) "    file three stores the mole numbers of precipitated minerals"
	    write(*,'(a,$)') "    file name is "
	    call color(hconsole,c_yellow)
	    write(*,*) num_S(1:len_trim(num_S)) // "." // "j" // syst_S(1:len_trim(syst_S)) // "%" 
	    call color(hconsole,c_white)
	    call display("    other name ? (else press <enter>) : ", x_S)
		if (x_S /= "") molemin_S=x_S
	END IF
	write(*,*)
	transf_S = num_S(1:len_trim(num_S)) // ".tra"
	IF (transf_S == ".tra")  transf_S = "default.tra"
	write(*,*) "Transfer file name is ",transf_S(1:len_trim(transf_S))
	call display ("  other name ? (else press <enter>) : ",transfneuf_S)
	IF (transfneuf_S /= "") transf_S = transfneuf_S
	write(*,*)
	call display ("Modify limits of convergence ? (y/n=<enter>):  ", x_S)
	IF (x_S == "y") THEN
	    write(*,*) "mole number iterations limit is .001"
	    call display ("other limit ? (else press <enter>):                 ", x_S)
		pkmol=0.
		if(x_S /= "") read(x_S,*) pkmol
	    IF (pkmol == 0.) pkmol = .001D0
	    write(*,*) "Newton-Raphson iterations limit is .0000000000001"
	    call display ("other limit ? (else press <enter>):                     ", x_S)
		pkeq=0.
		if(x_S /= "") read(x_S,*) pkeq
	    IF (pkeq == 0.)  pkeq = .0000000000001D0
	ELSEIF (x_S == "n" .OR. x_S == "") THEN
	    pkmol = .001D0
	    pkeq = .0000000000001D0
	END IF

	write(*,*) 
	call color(hconsole,c_lred)
	call display( "final agreement ? (n/y=<enter>) ", x_S)
	call color(hconsole,c_white)
	write(*,*) ' '
	call fdate(d_and_t_S)
	write(*,*) d_and_t_S

	IF (x_S == "n") STOP
	IF (y_S == "y") THEN
		IF (fich_S == ".")  fich_S = num_S(1:len_trim(num_S)) // "." // "j" // syst_S(1:len_trim(syst_S)) // "&"
    	IF (constit_S(ialc) == "alc" .OR. e_S == "c") THEN 
			constituant_S = "numero,fc,eva,ds,ph,alc"
		ELSE IF (constit_S(ialc) == "alk" .OR. e_S == "k") THEN 
			constituant_S = "label,fc,eva,ds,ph,alk"
		ENDIF
    	DO i = 1 , 8
        	IF (tot(i) > 0.)  then
        		zone_S=f_lcase(aq_S(i))
        		constituant_S = constituant_S(1:len_trim(constituant_S)) // "," // zone_S(1:len_trim(zone_S))
        	END IF
    	ENDDO
    	IF (tot(9) /= 0.)  constituant_S = constituant_S(1:len_trim(constituant_S)) // ",b"
    	IF (tot(10) /= 0.)  constituant_S = constituant_S(1:len_trim(constituant_S)) // ",si"
    	IF (constit_S(ialc) == "alc" .OR. e_S == "c") THEN 
			constituant_S = constituant_S(1:len_trim(constituant_S)) // ",std"
		ELSE IF (constit_S(ialc) == "alk" .OR. e_S == "k") THEN 
			constituant_S = constituant_S(1:len_trim(constituant_S)) // ",tds"
		ENDIF
		IF (fichmin_S == ".")  fichmin_S = num_S(1:len_trim(num_S)) // "." // "j" // syst_S(1:len_trim(syst_S)) // "@"
		IF (molemin_S == ".")  molemin_S = num_S(1:len_trim(num_S)) // "." // "j" // syst_S(1:len_trim(syst_S)) // "%"
	END IF

	OPEN (10,file="stockage",action="WRITE")
    write(10,*) transf_S(1:len_trim(transf_S))
	CLOSE (10)

	OPEN(10,file=transf_S,action="WRITE")
    write(10,*)  temp
	write(10,*)  tinit
    write(10,*)  ph
	write(10,*)  phinit
    write(10,*)  po
	write(10,*)  poinit
	write(10,*)  diltot
	write(10,'(a)')  constituant_S(1:len_trim(constituant_S))
    DO i = 1 , 10
        write(10,*)  tot(i)
    ENDDO
    write(10,*)  molal(15)
    DO i = 1 , 10
        write(10,*)  molal(i)
    ENDDO
    write(10,*)  mh2o
    write(10,*)  molal(13)
    write(10,*)  molal(11)
    write(10,*)  molal(14)
    write(10,*)  molal(12)
    DO i = 16 , 25
        write(10,*)  molal(i)
    ENDDO
    write(10,'(a)')  syst_S(1:len_trim(syst_S))
    write(10,*)  xi
	write(10,*)  npasecran
    write(10,'(a)')  p_S(1:len_trim(p_S))
    write(10,*)  npasimp
    write(10,'(a)')  y_S(1:len_trim(y_S))
    write(10,*)  npasfich
    write(10,'(a)')  unit_S(1:len_trim(unit_S))
    write(10,'(a)')  fich_S(1:len_trim(fich_S))
    write(10,'(a)')  fichmin_S(1:len_trim(fichmin_S))
    write(10,'(a)')  molemin_S(1:len_trim(molemin_S))
    write(10,'(a)')  min_S(1:len_trim(min_S))
    write(10,*)  stdmax
    write(10,*)  pkmol
    write(10,*)  pkeq
	CLOSE (10)

result_evap=system("evp.exe") ! run the evaporation module

  contains

	subroutine input_data
		character(len=maxlen)	:: f_S,zone,nu_S,tdif_S,x_S,fname_S,xx_S
		character(len=maxlen), allocatable, dimension (:)	:: v_S
		integer	:: i,j,k,nal,max_len_var,nb_col,itemp,compt_al,compt_fn,nf,max_len_fn
		real(kind=8), allocatable, dimension (:)	:: v

		allocate (v_S(0:max_constit), v(0:max_constit), stat=i)

		call display ("Input file/keyboard (f/k) : ", e_S)
		IF (e_S == "")  e_S = "f"
		IF (e_S == "k" .OR. e_S == "c") THEN
			write(*,*) "Enter concentrations in MILLIMOLES/liter of solution (molarities)"
			write(*,*) "                  or in MILLIMOLES/kg/H2O (molalities)"
			write(*,*) "To enter MOLALITIES ==> press <enter> for density"; write(*,*)
			call display ( "Label        : ", num_S)
			call display ( "Temp(deg. C) = ", temp)
			call display ( "Density      = ", x_S)
			dens=0.
			IF (x_S /= " ") read(x_S,*) dens
			call display ( "pH   = ", ph)
			call display ( "Na   = ", tot(1))
			call display ( "K    = ", tot(2))
			call display ( "Li   = ", tot(3))
			call display ( "Ca   = ", tot(4))
			call display ( "Mg   = ", tot(5))
			call display ( "Cl   = ", tot(6))
			call display ( "SO4  = ", tot(7))
			call display ( "Alk  = ", tot(12))
			call display ( "NO3  = ", tot(8))
			call display ( "Si   = ", tot(10))
			call display ( "B    = ", tot(9))
			DO i = 1 , ntot
				tot(i) = tot(i) / 1000.
			ENDDO
			tinit = temp
		ELSEIF (e_S == "f") THEN
			write(*,*) ' '
			call color(hconsole,c_lred)
			call display("Directory (current directory=<enter>) : ",rep_S)
			if (rep_S == "") rep_S=repinit_S
			write(*,*) 'Working Directory : ', rep_S(1:len_trim(rep_S))
			call color(hconsole,c_cyan)
			IF (OS_S == "winNT" .OR. OS_S == "win9x") THEN
				x_S="dir " // rep_S(1:len_trim(rep_S)) // " /b  > file.lst"
			ELSE
				x_S="ls " // rep_S(1:len_trim(rep_S)) // "  > file.lst"
			ENDIF
			result_evap=system(x_S) 
! here we display on screen the names of the files in the choosen directory
			nf=0
			max_len_fn=0
			nconstit=1
			open (10,file='file.lst',action='READ')
				do 
					read(10,*,end=100) fname_S
					nf=nf+1 ! Number of files
					max_len_fn=max(max_len_fn,len_trim(fname_S)) ! max size of file names
				enddo
				100 continue
				rewind (10)				
				max_len_fn=max_len_fn+4 ! more space for the index number
				compt_fn=0
				! now we display the names on the screen
				nb_col=len_display/(max_len_fn+1)
				do i=1,nf,nb_col
				  if(i+nb_col-1.le.nf) then
					zone=' '
					do j=1,nb_col
						read(10,*) fname_S
						compt_fn=compt_fn+1
						write(zone((j-1)*(max_len_fn+1)+1:),'(i3)') compt_fn
						zone((j-1)*(max_len_fn+1)+5:)=fname_S
					enddo
					write(*,*) zone(1:len_trim(zone))
				  endif
				enddo
				zone=' '
				do i=1,mod(nf,nb_col)
					read(10,*) fname_S
					compt_fn=compt_fn+1
					write(zone((i-1)*(max_len_fn+1)+1:),'(i3)') compt_fn
					zone((i-1)*(max_len_fn+1)+5:)=fname_S
				enddo
				write(*,*) zone(1:len_trim(zone))
			close (10)
			call color(hconsole,c_white)
			call display("File Name: ",f_S)
			f_S=rep_S(1:len_trim(rep_S)) // DIR_SEP_S // f_S
	
! Here we read the sample to be used for the simulation		
			nal=0
			max_len_var=0
			nconstit=1
			open (10,file=f_S,action='READ')
				read(10,'(a)') x_S   ! list of components
				IF (index(x_S,',') /= 0) THEN ! , is the separator
					do i=1,len_trim(x_S)
						IF (x_S(i:i) == ',') nconstit=nconstit+1
					enddo
				ELSE
					xx_S=x_S					
					do while(len_trim(xx_S) > 0)
						xx_S=adjustl(xx_S); xx_S=xx_S(1:len_trim(xx_S))
						k=index(xx_S,' ')
						IF (k>1) THEN
							nconstit=nconstit+1
						ENDIF
						xx_S=xx_S(k:len_trim(xx_S))
					enddo
					nconstit=nconstit-1
				END IF
				read(x_S,*) (constit_S(i), i=1,nconstit) ! store the names of components

				do 
					read(10,*,end=50) num_S
					nal=nal+1 ! N of samples
					max_len_var=max(max_len_var,len_trim(num_S)) ! max size of samples names
				enddo
				50 continue
				rewind (10)
				
				max_len_var=max_len_var+4 ! more space for the index number
				compt_al=0
				! now we display the names on the screen
				read(10,*) num_S ! skip first record = components names
				nb_col=len_display/(max_len_var+1)
				do i=1,nal,nb_col
				  if(i+nb_col-1.le.nal) then
					zone=' '
					do j=1,nb_col
						read(10,*) num_S
						compt_al=compt_al+1
						write(zone((j-1)*(max_len_var+1)+1:),'(i3)') compt_al
						zone((j-1)*(max_len_var+1)+5:)=num_S
					enddo
					write(*,*) zone(1:len_trim(zone))
				  endif
				enddo
				zone=' '
				do i=1,mod(nal,nb_col)
					read(10,*) num_S
					compt_al=compt_al+1
					write(zone((i-1)*(max_len_var+1)+1:),'(i3)') compt_al
					zone((i-1)*(max_len_var+1)+5:)=num_S
				enddo
				write(*,*) zone(1:len_trim(zone))

				! now we read the name of the sample to be used
				rewind (10)
				read(10,*) zone ! skip first record = components names
				IF (f_S(len_trim(f_S)-1:len_trim(f_S))=="&") THEN
					call display("Number of analysis (left number) : ", nu_S)
				ELSE
					call display("Number or name of analysis: ", nu_S)
				ENDIF
				nu_S=adjustl(nu_S)
				IF ( (nu_S(1:1) >= '0' .and. nu_S(1:1) <= '9') .OR. (f_S(len_trim(f_S)-1:len_trim(f_S))=="&") ) THEN
					read(nu_S,*) k
					do i=1,nal
						read(10,*) x_S
						IF(k == i) THEN
							num_S=f_lcase(x_S)
							call color(hconsole,c_bwhite)
							write(*,*) " ",num_S(1:len_trim(num_S))  ! found the sample to be used
							call color(hconsole,c_white)
							exit
						ENDIF
					enddo
				ELSE
					num_S=nu_S
				ENDIF

				rewind(10)
				read(10,*) zone ! skip first record = components names
				v_S(1)='0.'  ! first component = sample label
				do k=1,nal
					read(10,*) x_S,(v_S(j),j=2,nconstit)
					x_S=f_lcase(x_S)
					IF (x_S == num_S) THEN
						call verif(v_S,nconstit,v) ! nd na replaced by 0
						do i=2,nconstit
							IF (constit_S(i) == "t")  THEN
								temp = v(i); itemp = i
							ENDIF
							IF (constit_S(i) == "ds")  dens = v(i)
							IF (constit_S(i) == "ph") ph = v(i)
							IF (constit_S(i) == "na")  tot(1) = v(i) / 1000.
							IF (constit_S(i) == "k")  tot(2) = v(i) / 1000.
							IF (constit_S(i) == "li")  tot(3) = v(i) / 1000.
							IF (constit_S(i) == "ca")  tot(4) = v(i) / 1000.
							IF (constit_S(i) == "mg")  tot(5) = v(i) / 1000.
							IF (constit_S(i) == "cl")  tot(6) = v(i) / 1000.
							IF (constit_S(i) == "so4")  tot(7) = v(i) / 1000.
							IF (constit_S(i) == "no3")  tot(8) = v(i) / 1000.
							IF (constit_S(i) == "alk")  THEN
								tot(12) = v(i) / 1000.;ialc=i
							ELSEIF (constit_S(i) == "alc")  THEN
								tot(12) = v(i) / 1000.;ialc=i
							ENDIF
							IF (constit_S(i) == "si")  tot(10) = v(i) / 1000.
							IF (constit_S(i) == "b")  tot(9) = v(i) / 1000.
						enddo
						exit
					ENDIF
				enddo
			close(10)
			write(*,*) ' '
			IF (x_S /= num_S) THEN
				write(*,*) "Analysis not found"
				stop
			ENDIF
    		tinit = temp
			IF (itemp > 1) THEN
				write(*,*) "temperature = ", temp, " deg. C"
				call display(" other temperature (else press <enter>): ", tdif_S)
				IF (tdif_S /= "")  read(tdif_S,*) temp 
			ELSEIF (itemp == 0) THEN
				write(*,*) "temperature is not specified"
				call display (" enter a temperature: temp = ", temp)
				tinit = temp
			END IF
		END IF
	end subroutine input_data

end program eql



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



SUBROUTINE densite (cat, ani, nchcat, nchani, unit_S, dens)
	use data_common
	implicit none

	integer, dimension(0:*)	:: nchcat, nchani
	real(kind=8)	::	dens
	real(kind=8), dimension(0:*)	::	cat, ani
	character(len=*)	:: unit_S

	character(len=maxlen) :: x_S
	real(kind=8), allocatable, dimension (:,:)	:: s,ao,bo,au,bu
	real(kind=8)	::	u
	integer	:: nc, na, i,j,k,alloc_res

  	nc = 5; na = 5

  	allocate (s(0:nc, 0:na), stat=alloc_res)
  	allocate (ao(0:nc, 0:na), bo(0:nc, 0:na), au(0:nc, 0:na), bu(0:nc, 0:na), stat=alloc_res)
	s=0.; ao=0.; au=0.; bo=0.; bu=0.

  	OPEN (10,file="densite",action='READ')
    DO i = 1 , 5
        DO j = 1 , 5
            read(10,*) x_S, ao(i, j), bo(i, j), au(i, j), bu(i, j)
        ENDDO
    ENDDO
  	CLOSE (10)
  	IF (unit_S == "molar") THEN
     dens = 1.; u = 0.
     DO j = 1 , na; u = u + ani(j); ENDDO
     DO i = 1 , nc
      DO j = 1 , na
        s(i, j) = (nchcat(i) + nchani(j)) / 2 * cat(i) * ani(j) / nchcat(i) / nchani(j) / u
        dens = dens + ao(i, j) * s(i, j) + bo(i, j) * s(i, j) ** 2
      ENDDO
     ENDDO
  	ELSEIF (unit_S == "molal") THEN
     dens = 1.; u = 0.
     DO j = 1 , na; u = u + ani(j) * nchani(j); ENDDO
     DO i = 1 , nc
      DO j = 1 , na
        s(i, j) = (nchcat(i) + nchani(j)) / 2 * cat(i) * ani(j) / u
        dens = dens + au(i, j) * s(i, j) + bu(i, j) * s(i, j) ** 2
      ENDDO
   	 ENDDO
  	END IF
	deallocate (s,ao,au,bo,bu,stat=alloc_res)
END SUBROUTINE densite


SUBROUTINE invar (nm, nt, wmin, mineral_S, lmin, psol, psc, kinvariant, kinvar)
	use data_common
	implicit none

	integer	:: nm,kinvariant,nt
	integer, dimension(0:*)	:: lmin,kinvar
	real(kind=8), dimension(0:*)	::	psol,psc
	real(kind=8), dimension(0:nm,0:nt)	:: wmin
	character(len=*), dimension(0:*)	:: mineral_S

	character(len=maxlen), allocatable, dimension (:)	:: minv_S,minvar_S
	real(kind=8), allocatable, dimension (:)	:: psminv,psminvar,tt4
	real(kind=8), allocatable, dimension (:,:)	:: winv,t0,t1,t2,t3,t4
	real(kind=8)	::	swap,u,alea1,alea2,ah2o
	integer, allocatable, dimension (:)	:: kinv
	integer	:: i,j,k,n1,n2,n3,n4,ii,jj,kk,ncond,nbmin,ninvar,alloc_res,det,det1,z,ncm,i1,i2

    kinvariant = 0; ncm = 14
    nbmin = 10
    ninvar = nbmin + 3

    allocate ( kinv(0:ninvar), minv_S(0:ninvar), psminv(0:ninvar), winv(0:ninvar, 0:ncm), stat=alloc_res)
    allocate ( minvar_S(0:ninvar), psminvar(0:ninvar), stat=alloc_res)
    allocate ( t0(0:ninvar, 0:ninvar), t1(0:ninvar, 0:ninvar), stat=alloc_res)
    allocate ( t2(0:ninvar, 0:ninvar), t3(0:ninvar, 0:ncm), t4(0:ncm, 0:ncm), tt4(0:ncm), stat=alloc_res)

	kinv=0; psminv=0.; winv=0.; psminvar=0.; t0=0.; t1=0.; t2=0.; t3=0.; t4=0.; tt4=0.

    DO k = 1 , 3
        psminv(k) = LOG(psc(k)) / LOG(10D0)
    ENDDO
    winv(1, 11) = 1.; winv(1, 13) = 1.; winv(1, 0) = -1.
    winv(2, 11) = 1.; winv(2, 12) = -1.; winv(2, 14) = 1.
    winv(3, 11) = 1.; winv(3, 12) = 1.; winv(3, 0) = -1.

    n1 = 3
    DO k = 1 , nm
      IF (lmin(k) == 1) THEN
        n1 = n1 + 1
        kinv(n1) = k
        minv_S(n1) = mineral_S(k)
        psminv(n1) = LOG(psol(k)) / LOG(10D0)
        DO j = 0 , ncm
          winv(n1, j) = wmin(k, j)
        ENDDO
      END IF
    ENDDO
    DO i = 1 , n1
        swap= winv(i, 0)
		winv(i, 0) = winv(i, 14)
		winv(i, 14) = swap
    ENDDO

    DO i = 1 , n1
      DO j = i , n1
        t1(i, j) = 0
        DO k = 0 , ncm - 1
          t1(i, j) = t1(i, j) + winv(i, k) * winv(j, k)
          t1(j, i) = t1(i, j)
          t0(i, j) = t1(i, j)
          t0(j, i) = t0(i, j)
        ENDDO
      ENDDO
    ENDDO
    DO k = 2 , n1
      DO i = k , n1
        IF (ABS(t1(i, k - 1)) > epsilon) THEN
          u = t1(k - 1, k - 1) / t1(i, k - 1)
          DO j = k , n1
            t1(i, j) = t1(k - 1, j) - t1(i, j) * u
            IF (ABS(t1(i, j)) < epsilon)  t1(i, j) = 0
          ENDDO
        END IF
      ENDDO
    ENDDO
    det = 1
    DO i = 1 , n1
      IF (ABS(t1(i, i)) < epsilon) THEN
        det = 0
        EXIT 
      END IF
    ENDDO

    IF (det == 0) THEN
      n3 = 0; n2 = n1 - 1
      DO kk = 1 , n1
        ii = 0
        DO i = 1 , n1
          IF (i .NE. kk) THEN
            ii = ii + 1
            jj = 0
            DO j = 1 , n1
              IF (j .NE. kk) THEN
                jj = jj + 1
                t2(ii, jj) = t0(i, j)
              END IF
            ENDDO
          END IF
        ENDDO

        DO k = 2 , n2
          DO i = k , n2
            IF (ABS(t2(i, k - 1)) > epsilon) THEN
              u = t2(k - 1, k - 1) / t2(i, k - 1)
              DO j = k , n2
                t2(i, j) = t2(k - 1, j) - t2(i, j) * u
                IF (ABS(t2(i, j)) < epsilon)  t2(i, j) = 0
              ENDDO
            END IF
          ENDDO
        ENDDO

        det1 = 1
        DO i = 1 , n2
          IF (ABS(t2(i, i)) < epsilon) THEN
            det1 = 0
            EXIT 
          END IF
        ENDDO
        IF (det1 == 1) THEN
          n3 = n3 + 1
          kinvar(n3) = kinv(kk)
          minvar_S(n3) = minv_S(kk)
          psminvar(n3) = psminv(kk)
          DO j = 0 , ncm
            t3(n3, j) = winv(kk, j)
          ENDDO
        END IF
      ENDDO
      
	    IF (n3 == 0) THEN
            kinvariant = -1
        ELSEIF (n3 > 0) THEN
            n4 = ncm
            DO j = ncm, 1, -1
				u = 0
                DO i = 1, n3
					u = u + t3(i, j) ** 2
                ENDDO
                IF (u < epsilon) THEN
                    DO k = j+1, n4
						DO i = 1, n3
                            t3(i, k - 1) = t3(i, k)
                        ENDDO
                    ENDDO
                    n4 = n4 - 1
                END IF
            ENDDO
                            
            DO i = 1, n4
				DO j = 1, n4
                    t4(i, j) = 0
                    DO k = 1, n3
					    t4(i, j) = t4(i, j) + t3(k, i) * t3(k, j)
                        t4(j, i) = t4(i, j)
                    ENDDO
                ENDDO
            ENDDO
            DO i = 1, n4
			    tt4(i) = 0
                DO k = 1, n3
				    tt4(i) = tt4(i) + t3(k, i) * psminvar(k)
                ENDDO
            ENDDO

            DO k = 2, n4
				DO i = k, n4
                        IF (ABS(t4(i, k - 1)) > epsilon) THEN
                        u = t4(k - 1, k - 1) / t4(i, k - 1)
                        DO j = k, n4
						    t4(i, j) = t4(k - 1, j) - t4(i, j) * u
                            IF (ABS(t4(i, j)) < epsilon) THEN 
								t4(i, j) = 0
							ENDIF
                        ENDDO
                        tt4(i) = tt4(k - 1) - tt4(i) * u
                    ENDIF
                ENDDO
            ENDDO
                    
            IF (ABS(t4(n4, n4)) > epsilon) THEN
                ah2o = 10 ** (tt4(n4) / t4(n4, n4))
                IF (ah2o > 1 .OR. ah2o <= .01) THEN
                    kinvariant = -2
                ELSE
                    kinvariant = n3
                    DO i = kinvariant, 1, -1
					    IF (kinvar(i) == 0) THEN
                            DO k = 1, kinvariant-1
							    kinvar(k) = kinvar(k + 1)
                            ENDDO
                            kinvariant = kinvariant - 1
                        ENDIF
                    ENDDO
                ENDIF
            ELSEIF (ABS(t4(n4, n4)) <= epsilon) THEN
                kinvariant = -2
            END IF
        END IF
    ENDIF 
          
    deallocate ( kinv, minv_S, psminv, winv, stat=alloc_res)
    deallocate ( minvar_S, psminvar, stat=alloc_res)
    deallocate ( t0,  t1,  stat=alloc_res)
    deallocate ( t2,  t3,  t4, tt4, stat=alloc_res)

END SUBROUTINE invar




DOUBLE PRECISION FUNCTION g0 (x)
	real(kind=8)	:: x
    g0 = 2 * (1 - (1 + x) * EXP(-x)) / x ** 2
END FUNCTION g0

DOUBLE PRECISION FUNCTION g1 (x)
	real(kind=8)	:: x
    g1 = -2 * (1 - (1 + x + x ** 2 / 2) * EXP(-x)) / x ** 2
END FUNCTION g1


DOUBLE PRECISION FUNCTION j0 (x)
    implicit none
    real(kind=8)	:: ya,yb,yc,yd,x
    ya = 4.581; yb = -.7237; yc = -.012; yd = .528
    j0 = x / (4 + ya * x ** yb * EXP(yc * x ** yd))
END FUNCTION j0

DOUBLE PRECISION FUNCTION j1 (x)
    implicit none
    real(kind=8)	:: ya,yb,yc,yd,x
    ya = 4.581; yb = -.7237; yc = -.012; yd = .528
    j1 = (4 + ya * x ** yb * (1 - yb - yc * yd * x ** yd) * EXP(yc * x ** yd)) / (4 + ya * x ** yb * EXP(yc * x ** yd)) ** 2
END FUNCTION j1


double precision FUNCTION temperature (at, bt, ct, dt, et, temp)
	implicit none
	real (kind=8) :: at,bt,ct,dt,et,temp
    temperature = at + bt * temp + ct * temp ** 2 + dt * temp ** 3 + et * temp ** 4
END FUNCTION temperature




SUBROUTINE color(hconsole,n)
!use dfwin !WIN32 only
use data_common
implicit none
integer res,n,BackColor,hconsole
integer(kind=2) entier2
BackColor=0
if(OS_S.eq.'winNT') then
!	entier2=IOR(n, ishft(BackColor, 4))
!	res=SetConsoleTextAttribute (hconsole, entier2) !WIN32 only
endif
END SUBROUTINE color

function f_lcase(c) result (minus)
	implicit none
	character (len=*) c
	character (len=1000) minus
	integer i
	minus=' '
	do i=1,len_trim(c)
		if(c(i:i).ge.'A'.and.c(i:i).le.'Z') then
			minus(i:i)=char(ichar(c(i:i))+32)
		else
			minus(i:i)=c(i:i)
		endif
	enddo
end function f_lcase

function f_ucase(c) result (minus)
	implicit none
	character (len=*) c
	character (len=1000) minus
	integer i
	minus=' '
	do i=1,len_trim(c)
		if(c(i:i).ge.'a'.and.c(i:i).le.'z') then
			minus(i:i)=char(ichar(c(i:i))-32)
		else
			minus(i:i)=c(i:i)
		endif
	enddo
end function f_ucase


subroutine int_display(text,i)
	character (len=*) :: text
	integer :: i
	write(*,'(a," ",$)') text(1:len_trim(text))
	read(*,*) i
end subroutine int_display

subroutine real_display(text,vi)
	character (len=*) :: text
	real :: vi
	write(*,'(a," ",$)') text(1:len_trim(text))
	read(*,*) vi
end subroutine real_display

subroutine dp_display(text,vi)
	character (len=*) :: text
	real(kind=8) :: vi
	write(*,'(a," ",$)') text(1:len_trim(text))
	read(*,*) vi
end subroutine dp_display

subroutine text_display(text,text_out)
	character (len=*) :: text
	character (len=*) :: text_out
	write(*,'(a," ",$)') text(1:len_trim(text))
	read(*,'(a)') text_out
end subroutine text_display


subroutine verif(v_car,n,v_val)
	character (len=*)	:: v_car(0:*)
	integer		:: n,i
	real(kind=8)	:: v_val(0:*)
	do i=1,n
		if(v_car(i).ne.'na'.and.v_car(i).ne.'nd') then
			read(v_car(i),*) v_val(i)
		else
			v_val(i)=0.
		endif
	enddo
end subroutine verif


subroutine OS_VERIF
	use data_common
	implicit none
	logical ios
	inquire(file='c:\command.com',exist=ios) 
	if(ios) then
		OS_S="win9x"
		DIR_SEP_S="\"
	endif
	inquire(file='c:\boot.ini',exist=ios)
	if(ios) then
		OS_S="winNT"
		DIR_SEP_S="\"
	endif
	inquire(file='/bin/ls',exist=ios)
	if(ios) then
		OS_S="unix"
		DIR_SEP_S="/"
	endif

end subroutine OS_VERIF
