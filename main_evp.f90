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

	integer, parameter	:: maxlen=256, len_display=50,  max_conv=100  !,max_constit=30, max_nr=10
	character (len=1)	::	test_stop, DIR_SEP_S
	character (len=5)	::	OS_S
	character (len=maxlen)	:: stockfich_S, molemin_S, fichmin_S
	integer	:: n, ntot, ncm, nm, nminer, nbmin, kinvariant, ncpt
	
	integer(kind=2)	:: int2
	integer	:: &
		c_bwhite=15, c_yellow=14, c_lmagenta=13, c_lred=12, c_lgreen=10, c_lcyan=11, &
		c_lblue=9, c_white=7, c_brown=6, c_red=4, c_cyan=3, c_green=2
	real(kind=8)	:: mh2o
	real(kind=8), parameter	:: epsilon=1.D-08
end module data_common

program evp
!EVP
!Simulation of evaporation
!System: Na-K-Li-Ca-Mg-Cl-SO4-NO3-CO2-HCO3-CO3-H4SiO4-B(OH)3-H-OH-H2O
!Ionic interaction (Pitzer)
!Solution of the system of equations by Newton-Raphson method
!Temperature dependence included for minerals
!  and for some Pitzer's coefficients (0-50C)
!Connected files: MATRICE2-MURTF3/MURTF0-COMPLEX3-COEFFT4-DENSITE
!Program chained to EQL through transfer file: *.tra
!Last revision: march 1999
!author: F. RISACHER
!translated in FORTRAN 90 by A. CLEMENT

	use dfwin		!WIN32 only
	use data_common
	implicit none

	real(kind=8)	:: fc, temp, tinit, ph, phinit, po, poinit, diltot, xi, stdmax, pkmol, pkeq, xi0
	real(kind=8)	:: at, bt, ct, dt, et, molec, u, sc, cmax, sa, amax, dca, delta, ctot0, debp, deb, fin
	real(kind=8)	:: xinv,gmin,m,p,g,psc3,fi,fi0,dens,std,ee,hco3,co3,ctot,s,pco,alc,mwev,mwev0
	integer	:: alloc_res,i,j,k,nc,na,ncomp,ic,ia,ncmpt,ix,iy,initdeseq,nminer0,nw,ksupprim,nu,nperitec
	integer	:: ncomplex, npasimp, npasfich, kneuf, npasi, npasf, ii, jj, l, ncmptinv, ki, npasecran
	character(len=maxlen)	:: q0_S, transf_S, constituant_S, syst_S, prnt_S, unit_S, fich_S, a_S, x_S
	character(len=maxlen)	:: constit_S,zone_S, m0_S, mineraux_S, my_S, my0_S, q_S, miner_S, inc_S, q1_S, q2_S, y_S
	real(kind=8), allocatable, dimension(:)	:: tot, totinit, tot0, totest, psc, gact, gact0, gact1, atom
	real(kind=8), allocatable, dimension(:)	:: mol, mol0, mol1, molal, molal0, act, act0, mu, mum, psol, psol0, pai, pai0
	real(kind=8), allocatable, dimension(:)	:: min, min0, minp, minp0 
	real(kind=8), allocatable, dimension(:,:)	:: wmin 
	real(kind=8), allocatable, dimension(:,:)	:: kmat
	integer, allocatable, dimension(:)	:: nch, ica, kinvar, linvar, lmin, lmin0, lmin1
	character(len=maxlen), allocatable, dimension(:)	:: aq_S, mineral_S, sauve_S
	character (len=maxlen), allocatable, dimension (:)	:: nom_ion_S
	real(kind=8), allocatable, dimension (:)	:: c_ion

	character(len=30)	:: d_and_t_S
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


! From here WIN 32 only
	lpszFileName="CONOUT$"C
	hConsole = CreateFile(lpszFileName, IOR(GENERIC_WRITE , GENERIC_READ),                   &
                         IOR(FILE_SHARE_READ , FILE_SHARE_WRITE),                           &
                         NULL_SECURITY_ATTRIBUTES, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0)
	if( hConsole == -1 ) then
      i =  GetLastError ()
      if (i .eq. ERROR_ALREADY_EXISTS) then
        i = 1
      end if
      write(*,*)
      write (*,*) ("Error: Unable to open console.")
      write(*,*)
      stop 1
	end if
! until there WIN 32 only

	call OS_VERIF

	call color(hconsole,c_lred)
	write(*,*)
	write(*,*) "This is EVP.............."
	write(*,*)

	write(*,*) "STARTING THE EVAPORATION PROGRAM"
	write(*,*) ; call color(hconsole,c_white)
	call fdate(d_and_t_S)
	write(*,*) d_and_t_S
	write(*,*)

	ncpt = 0; mwev = 0; fc = 1; q0_S = ""
	n = 25; ntot = 12; ncomplex = 14
	mh2o = 55.51D0

	allocate (tot(0:ntot), totinit(0:ntot), tot0(0:ntot), totest(0:ntot), psc(0:ncomplex), stat=alloc_res)
	allocate (gact(0:n), gact0(0:n), gact1(0:n),aq_S(0:n), atom(0:n), kmat(0:n, 0:n), nch(0:n), stat=alloc_res)
	allocate (mol(0:n), mol0(0:n), mol1(0:n), molal(0:n), molal0(0:n), act(0:n), act0(0:n), stat=alloc_res)
	tot=0.;  totinit=0.;  tot0=0.;  totest=0.;  psc=0.; 
	gact=0.; gact0=0.; gact1=0.; atom=0.; kmat=0.; nch=0 
	mol=0.; mol0=0.; mol1=0.; molal=0.; molal0=0.; act=0.; act0=0.

	open(10,file='aquv.dat',action='READ')
	DO i = 0 , n
		read(10,*) aq_S(i), atom(i), nch(i)
	ENDDO

	OPEN (10,file="matrice2",action='READ')
		DO i = 1 , n
			read(10,*)  (kmat(i, j),j=1,n)
		ENDDO
	CLOSE (10)

	OPEN (10,file="stockage",action='READ')
		read(10,*)  transf_S
	CLOSE (10)
	OPEN (10,file=transf_S,action='READ')
		read(10,*)  temp
		read(10,*)  tinit
		read(10,*)  ph
		read(10,*)  phinit
		read(10,*)  po
		read(10,*)  poinit
		read(10,*)  diltot
		read(10,'(a)')  constituant_S
		nbmin = 0
		DO i = 1 , 10
			read(10,*)  totinit(i)
			IF (totinit(i) /= 0)  nbmin = nbmin + 1
		ENDDO
		allocate (ica(0:n + nbmin), kinvar(0:nbmin + 3), stat=alloc_res)
		DO i = 0 , n
			read(10,*)  mol(i)
			mol0(i) = mol(i)
		ENDDO
		read(10,'(a)')  syst_S
		read(10,*)  xi
		read(10,*)  npasecran
		read(10,'(a)')  prnt_S
		read(10,*)  npasimp
		read(10,'(a)')  stockfich_S
		read(10,*)  npasfich
		read(10,'(a)')  unit_S
		read(10,'(a)')  fich_S
		read(10,'(a)')  fichmin_S
		read(10,'(a)')  molemin_S
		read(10,'(a)')  miner_S
		read(10,*)  stdmax
		read(10,*)  pkmol
		read(10,*)  pkeq
	CLOSE (10)
	
	IF (npasimp > 0) THEN
		open (9,file=prnt_S,position='APPEND')  ! log file
	END IF

	IF (xi == 0) THEN 
		inc_S = "auto" 
	ELSE 
		inc_S = "manu"
	END IF
	xi0 = xi
	DO i = 1 , n;  ica(i) = 1; ENDDO

	IF (stockfich_S == "y") THEN
		OPEN (10,file=fich_S, action='WRITE')
			write(10,*) constituant_S(1:len_trim(constituant_S))
		CLOSE (10)
		OPEN (10,file=fichmin_S, action='WRITE')
			write(10,*) fichmin_S(1:len_trim(fichmin_S))
			write(10,'(a,g13.6,a,g13.6,a)') "Temperature of solution = ", real(tinit), &
			            "Deg C.   Temperature of simulation = ", real(temp), "Deg C."
			IF (diltot > 1) THEN
				write(10,*) "The initial solution has been diluted ", real(diltot), " times"
			END IF
			IF (ph /= phinit) THEN
				write(10,'(a,g13.6,a,g13.6)') "Initial Log(pco2) = ", real(poinit), "    Selected Log(pco2) = ", real(po)
				write(10,'(a,g13.6,a,g13.6)') "Initial ph = ", real(phinit), "     Calculated ph = ", real(ph)
			END IF
		CLOSE (10)
		OPEN (10, file=molemin_S, action='WRITE')
		CLOSE (10)
	END IF

	OPEN (10,file="complex3",action='READ')
		DO i = 1 , ncomplex
			read(10,*)  a_S, at, bt, ct, dt, et
			psc(i) = 10 ** (at + bt / 300.D0 * temp + ct / 30000.D0 * temp ** 2 + &
			          dt / 3000000.D0 * temp ** 3 + et / 300000000.D0 * temp ** 4)
		ENDDO
	CLOSE (10)
	psc3 = psc(3)
	psc(3) = psc(3) * psc(14) * 10 ** po

	OPEN (10,file=miner_S,action='READ')
		read(10,*)  nc, na, nm
		ncm = nc + na + 1
		allocate (wmin(0:nm, 0:ncm), mu(0:ncm), linvar(0:nm), stat=alloc_res)
		allocate (mineral_S(0:nm), mum(0:nm), psol(0:nm), psol0(0:nm), pai(0:nm), pai0(0:nm), stat=alloc_res)
		allocate (lmin(0:nm), lmin0(0:nm), lmin1(0:nm), stat=alloc_res)
		allocate (min(0:nm), min0(0:nm), minp(0:nm), minp0(0:nm), stat=alloc_res)
		allocate (nom_ion_S(0:ncm), c_ion(0:ncm), stat=alloc_res)
   
   		wmin=0.; mu=0.; linvar=0
   		mum=0.; psol=0.; psol0=0.; pai=0.; pai0=0.
   		lmin=0; lmin0=0; lmin1=0
   		min=0.; min0=0.; minp=0.; minp0=0.
   
		DO i = 1 , ncm
			read(10,*)  a_S, at, bt, ct, dt, et
			a_S = f_lcase(a_S)
			DO j = 1 , ncm
				IF (a_S == aq_S(j)) THEN
					mu(j) = at + bt / 300.D0 * temp + ct / 30000.D0 * temp ** 2 + &
					        dt / 3000000.D0 * temp ** 3 + et / 300000000.D0 * temp ** 4
				END IF
			ENDDO
		ENDDO
		DO k = 1 , nm
			read(10,*)  mineral_S(k),ncomp, (c_ion(i),nom_ion_S(i),i=1,ncomp), at, bt, ct, dt, et
			DO i = 1 , ncomp
				x_S = f_lcase(nom_ion_S(i))
				DO j = 0 , ncm
					IF (x_S == aq_S(j))  wmin(k, j) = c_ion(i)
				ENDDO
			ENDDO
			mum(k) = at + bt / 300.D0 * temp + ct / 30000.D0 * temp ** 2 + &
			         dt / 3000000.D0 * temp ** 3 + et / 300000000.D0 * temp ** 4
        ENDDO
	CLOSE (10)

	DO k = 1 , nm
		u = mum(k)
		DO i = 0 , ncm
			u = u - wmin(k, i) * mu(i)
		ENDDO
		psol(k) = EXP(u)
		psol0(k) = psol(k)
	ENDDO

	IF (stockfich_S == "y") THEN
		OPEN (10,file=molemin_S,action='WRITE')
		  write(10,'(a,$)') "fc,"
		  DO k = 1 , nm - 1
			zone_S=mineral_S(k)(1:len_trim(mineral_S(k)))//","
			write(10,'(a,$)') zone_S(1:len_trim(zone_S))
		  ENDDO
                  zone_S=mineral_S(nm)(1:len_trim(mineral_S(nm)))
	          write(10,'(a)') zone_S(1:len_trim(zone_S))
		CLOSE (10)
	END IF
	sc = 0.; cmax = 0.
	DO i = 1 , n
		IF (nch(i) > 0) THEN
			sc = sc + mol(i) * nch(i)
			IF (mol(i) * nch(i) > cmax) THEN
				cmax = mol(i) * nch(i)
				ic = i
			END IF
		END IF
	ENDDO
	sa = 0.; amax = 0.
	DO i = 1 , n
		IF (nch(i) < 0 ) THEN
			sa = sa + mol(i) * (-nch(i))
			IF (mol(i) * (-nch(i)) > amax) THEN
				amax = mol(i) * (-nch(i))
				ia = i
			END IF
		END IF
	ENDDO
	dca = 200. * ABS(sc - sa) / (sc + sa)
	delta = sc - sa
	mol(ic) = mol(ic) - delta / 2. / nch(ic)
	mol(ia) = mol(ia) + delta / 2. / (-nch(ia))
	sc = 0.; sa = 0.
	DO i = 1 , n
		IF (nch(i) > 0)  sc = sc + mol(i) * nch(i)
		IF (nch(i) < 0)  sa = sa + mol(i) * (-nch(i))
	ENDDO
	write(*,*) "Sum of cations = ", real(sc)
	write(*,*) "Sum of anions  = ", real(sa)
	write(*,*)

	DO i = 1 , 11
		totinit(i) = 0.
		DO j = 1 , n
			totinit(i) = totinit(i) + kmat(i, j) * mol(j)
			tot(i) = totinit(i)
			tot0(i) = totinit(i)
		ENDDO
	ENDDO
	tot(12) = 0.
	ctot0 = mol(0) + mol(14) + mol(15) + mol(16) + mol(17)
	CALL actp(mol, gact, fi, temp)
	DO i = 0 , n
		molal(i) = mol(i) * mh2o / mol(11)
		act(i) = molal(i) * gact(i)
	ENDDO
	DO k = 1 , nm
		pai(k) = 1.
		DO i = 1 , ncm
			pai(k) = pai(k) * act(i) ** wmin(k, i)
		ENDDO
		IF (pai(k) >= psol(k))  lmin(k) = 1
	ENDDO

	IF (stockfich_S == "y") THEN
		OPEN (10,file=fichmin_S,position='APPEND')
			DO k = 1 , nm
				IF (lmin(k) == 1) THEN
					write(10,*) "Initial solution oversaturated in " // mineral_S(k)(1:len_trim(mineral_S(k)))
				END IF
			ENDDO
		CLOSE (10)
	END IF
       
	DO k = 1 , nm
		IF (lmin(k) == 1) THEN
			psol(k) = pai(k)
		END IF
	ENDDO


	500 continue   !debut:
	DO
		ncmpt = 0
		ix = 1; iy = 2
		initdeseq = 0
		IF (mwev == 0.) THEN
			DO k = 1 , nm
				IF (psol(k) > psol0(k))  initdeseq = 1
			ENDDO
			IF (initdeseq == 1) THEN
				DO k = 1 , nm
					IF (psol(k) * .95D0 > psol0(k)) THEN
						psol(k) = psol(k) * .95D0
						write(*,'(a,t20,g13.6,t35,g13.6)') &
						       mineral_S(k)(1:len_trim(mineral_S(k))),real(psol(k)),real(psol0(k))
					ELSEIF (psol(k) * .95D0 <= psol0(k)) THEN
						psol(k) = psol0(k)
					END IF
				ENDDO
			END IF
		ENDIF
   
		nw=1
		DO while (nw /= 0) 

			ncmpt = ncmpt + 1
			m0_S = ""
			DO k = 1 , nm
				IF (lmin(k) == 1) THEN
					m0_S = m0_S(1:len_trim(m0_S)) // "_" // mineral_S(k)(1:len_trim(mineral_S(k)))
				END IF
			ENDDO
			DO i = 0 , n; gact1(i) = gact(i); ENDDO
			CALL actp(mol, gact, fi, temp)
			IF (kinvariant == 0) then
			   DO i = 0 , n; gact(i) = (gact(i) + gact1(i) * ix) / iy; ENDDO
			END IF
			DO i = 0 , n
				molal(i) = mol(i) * mh2o / mol(11)
				act(i) = molal(i) * gact(i)
			ENDDO
			DO k = 1 , nm
				pai(k) = 1.
				DO i = 1 , ncm
					pai(k) = pai(k) * act(i) ** wmin(k, i)
				ENDDO
				IF (pai(k) >= psol(k)) THEN
					IF (min(k) >= 0.) THEN
						lmin(k) = 1
					ELSEIF (min(k) < 0.) THEN
						lmin(k) = 0
						min(k) = 0.
					END IF
				ELSEIF (pai(k) < psol(k)) THEN
					IF (min(k) <= 0) THEN
						lmin(k) = 0
						min(k) = 0.
					ELSEIF (min(k) > 0.) THEN
						lmin(k) = 1
					END IF
				END IF
			ENDDO
       
			DO k = 1 , nm
				IF (psol(k) == 1.D+50) THEN
					IF (pai(k) < psol0(k)*.9) THEN
						linvar(k) = 0
					ELSEIF (pai(k) >= psol0(k)) THEN
						linvar(k) = 1
					END IF
				END IF
			ENDDO
      
			mineraux_S = ""
			nminer = 0
			DO k = 1 , nm
				IF (lmin(k) == 1) THEN
					nminer = nminer + 1
					mineraux_S = mineraux_S(1:len_trim(mineraux_S)) // "_" &
					            // mineral_S(k)(1:len_trim(mineral_S(k)))
				END IF
			ENDDO
			IF (ncpt == 1 .OR. MOD(ncpt,npasecran) == 0) THEN
				IF (nminer == 0) THEN
					write(*,'(i4,t6,a)') ncmpt,"No_minerals"
				ELSE
					write(*,'(i4,1x,a)') ncmpt,mineraux_S(2:len_trim(mineraux_S))
				END IF
			END IF
              
			IF (mwev > 0. .AND. fc /= 1. .AND. nminer - nminer0 >= 2) THEN
				xi = xi / 2.
				IF (xi < epsilon) THEN
					write(*,*)
					write(*,*) "Program unstable"
					write(*,*) "Restart the initialisation program (EQP...)"
					write(*,*) "and lower the limits of convergence"
					IF (stockfich_S == "y") THEN
						open(10,file=fichmin_S,position='APPEND')
						write(10,*) "Program unstable"
						write(10,*) "Restart the initialisation program (EQP...)"
						write(10,*) "and lower the limits of convergence"
						close(10)
						write(*,*)
						call color(hconsole,c_yellow)
						call display("Compact mineral file ? (y/n=<enter):  ", x_S)
						call color(hconsole,c_white)
						IF (x_S == "y")  CALL compact(lmin)
					END IF
					call stop_simulation
				END IF
				write(*,*) "reduction of increment at ", xi
				DO i = 0 , n
					mol(i) = mol0(i)
				ENDDO
				DO i = 1 , ntot
					tot(i) = tot0(i)
				ENDDO
				DO k = 1 , nm
					lmin(k) = lmin0(k)
					min(k) = min0(k)
					minp(k) = minp0(k)
				ENDDO
				mwev = mwev0
				nminer = nminer0
				mwev = mwev + mol(11) * xi
				tot(11) = tot(11) - 2 * mol(11) * xi
				GOTO 500 !debut
			END IF
       
			IF (nminer > 1 .AND. mineraux_S /= m0_S) THEN
				ix = 2; iy = 3
				CALL invar(act(11), wmin, mineral_S, lmin, psol, psc, kinvar)
				IF (kinvariant > 0) THEN
					IF (syst_S == "o") THEN
						DO i = 1 , kinvariant
							IF (minp(kinvar(i)) == 0 .AND. min(kinvar(i)) == 0)  kneuf = kinvar(i)
						ENDDO
						GOTO 2000  !energie
					ELSEIF (syst_S == "c") THEN
						DO i = 0 , n
							mol(i) = mol0(i)
							molal(i) = molal0(i)
							gact(i) = gact0(i)
							act(i) = act0(i)
						ENDDO
						DO i = 1 , ntot
							tot(i) = tot0(i)
						ENDDO
						DO k = 1 , nm
							pai(k) = pai0(k)
							min(k) = min0(k)
						ENDDO
						mwev = mwev0
						nminer = nminer0
						fi=fi0
					END IF
				END IF
			END IF
       
			DO i = 0 , n; mol1(i) = mol(i); ENDDO
			IF (kinvariant == 0) THEN
				CALL reseq(nch, tot, mol, psol, psc, kmat, wmin, gact, act(11), ica, lmin, min, mineral_S, pkeq)
			ELSEIF (kinvariant > 0) THEN
				CALL reseqinv(lmin, wmin, kmat, min, kinvar, xinv)
				mwev = mwev + xinv/2.
				tot(11) = tot(11) - xinv
			END IF
			mol(0) = mol(15) * gact(15) * mol(13) * gact(13) / mol(11) / gact(11) / psc3 / gact(0)
			nw = 0
			DO i = 1 , n
				IF (mol(i) > 0) THEN
					IF (200. * ABS(mol(i) - mol1(i)) / (mol(i) + mol1(i)) > pkmol)  nw = 1
				 END IF
			ENDDO
			ki=kinvariant
			IF (kinvariant > 0) THEN
				DO k = 1 , kinvariant
					IF (min(kinvar(k)) <= 0.) THEN
						lmin(kinvar(k)) = 0
						min(kinvar(k)) = 0.
						psol(kinvar(k)) = 1.D+50
						mwev = mwev + mol(11) * xi
			                        tot(11) = tot(11) - 2. * mol(11) * xi
						ki = 0; nw = 1
					END IF
				ENDDO
			END IF
			kinvariant=ki
			IF (nw == 1) THEN 
        		DO i = 0 , n; mol(i) = (mol(i) + mol1(i)) / 2.; ENDDO
			END IF
			IF (ncmpt == 500) THEN
				write(*,*)
				write(*,*) "Program unstable"
				write(*,*) "Restart the initialisation program (EQP...)"
				write(*,*) "and lower the limits of convergence."
				write(*,*) "Set the increment in manual mode at a value lower than .5"
				IF (stockfich_S == "y") THEN
					open(10,file=fichmin_S,position='APPEND')
					write(10,*)
					write(10,*) "Program unstable"
					write(10,*) "Restart the initialisation program (EQP...)"
					write(10,*) "and lower the limits of convergence."
					write(10,*) "Set the increment in manual mode at a value lower than .5"
					close(10)
					write(*,*)
					call color(hconsole,c_yellow)
					call display("Compact mineral file ? (y/n=<enter):  ", x_S)
					call color(hconsole,c_white)
					IF (x_S == "y")  CALL compact(lmin)
				END IF
				call stop_simulation
			END IF
		ENDDO 

		DO k = 1 , nm
			IF (psol(k) == 1.D+50 .AND. linvar(k) == 0) THEN
				psol(k) = psol0(k)
				call color(hconsole,c_bwhite)
				write(*,*) "resetting: ", mineral_S(k)(1:len_trim(mineral_S(k))); call color(hconsole,c_white)
			END IF
			linvar(k) = 0
		ENDDO

		npasi = npasi + 1
		npasf = npasf + 1
		ncpt = ncpt + 1
   
		IF (syst_S == "o") THEN
			DO k = 1 , nm
				minp(k) = minp(k) + min(k)
			ENDDO
			DO i = 1 , 10
				DO k = 1 , nm
					tot(i) = tot(i) - wmin(k, i) * min(k)
				ENDDO
			ENDDO
			DO j = 1 , ncm
				DO k = 1 , nm
					tot(11) = tot(11) - wmin(k, j) * kmat(11, j) * min(k)
				ENDDO
			ENDDO
		END IF
		DO i = 1 , ntot - 1
			totest(i) = 0
			DO j = 1 , n
				totest(i) = totest(i) + kmat(i, j) * mol(j)
			ENDDO
		ENDDO
		totinit(12) = 0.; totest(12) = 0.
		DO j = 1 , n
			IF (kmat(12, j) > 0.) THEN
				totest(12) = totest(12) + kmat(12, j) * mol(j)
			ELSEIF (kmat(12, j) < 0.) THEN
				totinit(12) = totinit(12) + kmat(12, j) * mol(j)
			END IF
		ENDDO
		totinit(12) = -totinit(12)
		DO i = 1 , 10
			DO k = 1 , nm
				IF (syst_S == "c")  totest(i) = totest(i) + min(k) * wmin(k, i)
				IF (syst_S == "o")  totest(i) = totest(i) + minp(k) * wmin(k, i)
			ENDDO
		ENDDO
		DO j = 1 , ncm
		  DO k = 1 , nm
			  IF (syst_S == "c")  totest(11) = totest(11) + wmin(k, j) * kmat(11, j) * min(k)
			  IF (syst_S == "o")  totest(11) = totest(11) + wmin(k, j) * kmat(11, j) * minp(k)
		  ENDDO
		ENDDO
		totest(11) = totest(11) + mwev * 2.
		CALL densite(mol, kmat, nch, dens)
   
		fc = mh2o / mol(11)
		alc = molal(12) - molal(13) + molal(14) * 2 + molal(15) + (molal(16) + molal(17)) * 2
		alc = alc + molal(18) + molal(19) + molal(20) + molal(21) * 2
		alc = alc + molal(22) + molal(23) + molal(24) - molal(25)
		std = -molal(11) * atom(11)
		DO i = 1 , n
			std = std + molal(i) * atom(i)
		ENDDO
		ee = 1000000D0 * dens / (1000D0 + std)
		hco3 = 0.; co3 = 0.
		DO k = 1 , nm
			IF (syst_S == "c") THEN
				hco3 = hco3 + wmin(k, 15) * min(k)
				co3 = co3 + wmin(k, 14) * min(k)
			ELSEIF (syst_S == "o") THEN
				hco3 = hco3 + wmin(k, 15) * minp(k)
				co3 = co3 + wmin(k, 14) * minp(k)
			END IF
		ENDDO
		ctot = mol(0) + mol(14) + mol(15) + mol(16) + mol(17) + hco3 + co3
		my_S = ""
		DO i = 1 , nm
			IF (lmin(i) == 1 .OR. min(i) /= 0.) THEN
				my_S = my_S(1:len_trim(my_S)) // "_" // mineral_S(i)(1:len_trim(mineral_S(i)))
			END IF
		ENDDO

	600 continue !ecran:
		IF (ncpt == 1 .OR. MOD(ncpt,npasecran) == 0) THEN
			write(*,*)
			q_S = ""
			call color(hconsole,c_lmagenta)
			IF (nminer > 0) THEN
				IF (syst_S == "c") THEN
					write(*,'(t25,a,t67,a)') "MOLES PREC","TESTS"
				ELSEIF (syst_S == "o") THEN
					write(*,'(t25,a,t39,a,t67,a)') "MOLES 1 STEP","MOLES TOT","TESTS"
				END IF
			END IF
			call color(hconsole,c_white)
			write(*,*)
			DO i = 1 , nm
				IF (lmin(i) == 1 .OR. min(i) /= 0.) THEN
					IF (syst_S == "c") THEN
						IF (min(i) > min0(i)) THEN
							q_S = 'P'
						ELSEIF (min(i) < min0(i)) THEN
							q_S = 'D'
						ELSEIF (min(i) == min0(i)) THEN
							q_S = "="
						END IF
					END IF
					x_S= mineral_S(i)(1:len_trim(mineral_S(i))) // " " // q_S(1:1)
					u = 200. * ABS(pai(i) - psol(i)) / (pai(i) + psol(i))
            
					IF (syst_S == "o")  then
						write(*,'(a,t24,G13.6,t38,g13.6,t66,g13.6)') &
								x_S(1:len_trim(x_S)), real(min(i)),real(minp(i)),real(u)
					ELSE
						write(*,'(a,t24,G13.6,t66,g13.6)') x_S(1:len_trim(x_S)), real(min(i)),real(u)
					ENDIF
				END IF
			ENDDO
			IF (my_S == "")  write(*,*) "No_minerals"
			write(*,*)
  
			call color(hconsole,c_bwhite)
			write(*,'(t11,"MOLES",t25,"MOLALITIES",t39,"ACT COEFF",t53,"MOLAL TOT",t67,"TESTS")') 
			write(*,*)
			call color(hconsole,c_white)
			DO i = 1 , ntot
				IF (tot(i) > 0 .OR. i==12) THEN
					u = 200 * ABS(totest(i) - totinit(i)) / (totest(i) + totinit(i))
					IF (i <= 10) THEN
						s = 0.
						DO j = 1 , n
							s = s + molal(j) * kmat(i, j)
						ENDDO
						write(*,'(a,t10,g13.6,t24,g13.6,t38,g13.6,t52,g13.6,t66,g13.6)') &
							aq_S(i)(1:len_trim(aq_S(i))),real(mol(i)),real(molal(i)),real(gact(i)),real(s),real(u)
					ELSEIF (i>10) THEN
						write(*,'(a,t10,g13.6,t24,g13.6,t38,g13.6,t66,g13.6)') &
							aq_S(i)(1:len_trim(aq_S(i))),real(mol(i)),real(molal(i)),real(gact(i)),real(u)
					END IF
				END IF
			ENDDO
			DO i = ntot+1 , n
				IF (mol(i) > 0.) THEN
					p = 1.
					DO j = 1 , n
						p = p * act(j) ** kmat(i, j)
					ENDDO
					u = 200. * ABS(p - psc(i - 12)) / (p + psc(i - 12))
					write(*,'(a,t10,g13.6,t24,g13.6,t38,g13.6,t66,g13.6)') &
						   aq_S(i)(1:len_trim(aq_S(i))),real(mol(i)),real(molal(i)),real(gact(i)),real(u)
				END IF
			ENDDO
			pco = LOG(act(15) * act(13) / act(11) / psc3 / psc(14)) / LOG(10.D0)
			u = 200. * ABS(pco - po) / (pco + po)
			write(*,'(a,t10,g13.6,t24,g13.6,t38,g13.6,t66,g13.6)') &
					aq_S(0)(1:len_trim(aq_S(0))),real(mol(0)),real(molal(0)),real(gact(0)),real(u)
			
			write(*,*)
			write(*,'(a,t40,a,g13.6)') transf_S(1:len_trim(transf_S)),"concentration factor = ", real(fc)
			write(*,'(a,g13.6,t40,a,g13.6)') "ionic strength      = ", real(fi),"salinity (g/kg)      = ", real(std)
			write(*,'(a,g13.6,t40,a,g13.6)') "activity of water   = ", real(act(11)),"water evapor. (mol)  = ", real(mwev)
			write(*,'(a,g13.6,t40,a,g13.6)') "pH                  = ", real(-LOG(act(13)) / LOG(10.D0)),&
											 "CO2 exchanged (mol)  = ", real(ctot - ctot0)
			write(*,'(a,g13.6,t40,a,g13.6)') "alkalinity (eq/kg)  = ", real(alc),"Log PCO2             = ", real(pco)
			write(*,'(a,i5,t40,a,g13.6)') "number of steps     = ", ncpt,"molal/molar factor   = ", real(ee / 1000.)
			IF (kinvariant == 0) THEN
				write(*,'(a,g13.6,t40,a,g13.6)') "increment (%)       = ", real(xi) * 100,"density              = ", real(dens)
			ELSE
				write(*,'(a,g13.6,t40,a,g13.6)') "increment (moles)   = ", real(xinv),    "density              = ", real(dens)
			END IF
			write(*,*)
		ENDIF

  
		IF (prnt_S /= "n" .AND. (ncpt == 1 .OR. my_S /= my0_S .OR. MOD(ncpt,npasimp) == 0)) THEN
			my_S = ""; q_S = ""
			write(9,*)
			IF (syst_S == "c") THEN
					write(9,'(t25,a,t67,a)') "MOLES PREC","TESTS"
				ELSEIF (syst_S == "o") THEN
					write(9,'(t25,a,t39,a,t67,a)') "MOLES 1 STEP","MOLES TOT","TESTS"
				END IF
			write(9,*)
			
			DO i = 1 , nm
				IF (lmin(i) == 1 .OR. min(i) /= 0.) THEN
					IF (syst_S == "c") THEN
						IF (min(i) > min0(i)) THEN
							q_S = "(p)"
						ELSEIF (min(i) < min0(i)) THEN
							q_S = "(d)"
						ELSEIF (min(i) == min0(i)) THEN
							q_S = "(=)"
						END IF
					END IF
										
					x_S= mineral_S(i)(1:len_trim(mineral_S(i))) // " " // q_S(1:3)
					u = 200. * ABS(pai(i) - psol(i)) / (pai(i) + psol(i))
            
					IF (syst_S == "o")  then
						write(9,'(a,t24,G13.6,t38,g13.6,t66,g13.6)') &
								x_S(1:len_trim(x_S)), real(min(i)),real(minp(i)),real(u)
					ELSE
						write(9,'(a,t24,G13.6,t66,g13.6)') x_S(1:len_trim(x_S)), real(min(i)),real(u)
					ENDIF
					my_S = my_S(1:len_trim(my_S)) // "_" // mineral_S(i)(1:len_trim(mineral_S(i)))
				END IF
			ENDDO
			IF (my_S == "")  write(9,*) "No_minerals"
			write(9,*)
              
			IF (ncpt == 1) THEN
				write(9,'(t14, "MOLES", t26, "MOLALITIES", t39, "ACT COEFF", t53, "MOLAL TOT", t67, "TESTS")')
				write(9,*)
			END IF
			
			DO i = 1 , ntot
				IF (tot(i) > 0 .OR. i==12) THEN
					u = 200 * ABS(totest(i) - totinit(i)) / (totest(i) + totinit(i))
					IF (i <= 10) THEN
						s = 0.
						DO j = 1 , n
							s = s + molal(j) * kmat(i, j)
						ENDDO
						write(9,'(a,t10,g13.6,t24,g13.6,t38,g13.6,t52,g13.6,t66,g13.6)') &
							aq_S(i)(1:len_trim(aq_S(i))),real(mol(i)),real(molal(i)),real(gact(i)),real(s),real(u)
					ELSEIF (i>10) THEN
						write(9,'(a,t10,g13.6,t24,g13.6,t38,g13.6,t66,g13.6)') &
							aq_S(i)(1:len_trim(aq_S(i))),real(mol(i)),real(molal(i)),real(gact(i)),real(u)
					END IF
				END IF
			ENDDO
		
			DO i = ntot+1 , n
				IF (mol(i) > 0.) THEN
					p = 1.
					DO j = 1 , n
						p = p * act(j) ** kmat(i, j)
					ENDDO
					u = 200. * ABS(p - psc(i - 12)) / (p + psc(i - 12))
					write(9,'(a,t10,g13.6,t24,g13.6,t38,g13.6,t66,g13.6)') &
						   aq_S(i)(1:len_trim(aq_S(i))),real(mol(i)),real(molal(i)),real(gact(i)),real(u)
				END IF
			ENDDO
			pco = LOG(act(15) * act(13) / act(11) / psc3 / psc(14)) / LOG(10.D0)
			u = 200. * ABS(pco - po) / (pco + po)
			write(9,'(a,t10,g13.6,t24,g13.6,t38,g13.6,t66,g13.6)') &
					aq_S(0)(1:len_trim(aq_S(0))),real(mol(0)),real(molal(0)),real(gact(0)),real(u)
			
			write(9,*)
			write(9,*) transf_S(1:len_trim(transf_S))
			write(9,*) "concentration factor     = ", real(fc)
			write(9,*) "ionic strength           = ", real(fi)
			write(9,*) "salinity                 = ", real(std / 1000.), " g/kg"
			write(9,*) "activity of water        = ", real(act(11))
			write(9,*) "pH                       = ", real(-LOG(act(13)) / LOG(10.D0))
			write(9,*) "alkalinity               = ", real(alc), " eq/kg"
			write(9,*) "water evaporated         = ", real(mwev), " moles"
			write(9,*) "CO2 exchanged            = ", real(ctot - ctot0), " moles"
			write(9,*) "Log PCO2                 = ", real(pco)
			write(9,*) "density                  = ", real(dens)
			write(9,*) "molal/molar factor       = ", real(ee/1000.)
			IF (kinvariant == 0) THEN
				write(9,*) "increment (%)            = ", real(xi) * 100
			ELSE
				write(9,*) "increment (moles)        = ", real(xinv)
			END IF
			write(9,*) "number of steps          = ", ncpt
			write(9,*)
		END IF

		IF (stockfich_S == "y") THEN
			q_S = ""
			write(y_S,*) real(fc)
			DO k = 1 , nm
				write(x_S,*) real(minp(k))
				IF (syst_S == "c") THEN
					IF (lmin(k) == 1 .AND. lmin0(k) == 0) THEN
						q_S = "start of precipitation of " // &
						mineral_S(k)(1:len_trim(mineral_S(k))) // &
						" at fc =" // y_S(1:len_trim(y_S))
					ELSEIF (lmin(k) == 1 .AND. lmin1(k) == 0) THEN
						IF (min(k) < min0(k)) THEN
							lmin1(k) = 1
							q_S = "end of precipitation and start of dissolution of " // &
							      mineral_S(k)(1:len_trim(mineral_S(k))) // " at fc =" // y_S(1:len_trim(y_S))
						END IF
					ELSEIF (lmin(k) == 1 .AND. lmin1(k) == 1) THEN
						IF (min(k) > min0(k)) THEN
							lmin1(k) = 0
							q_S = "end of dissolution and start of precipitation of " // &
							      mineral_S(k)(1:len_trim(mineral_S(k))) // " at fc =" // y_S(1:len_trim(y_S))
						END IF
					ELSEIF (lmin(k) == 0 .AND. lmin1(k) == 1 .AND. lmin0(k) == 1) THEN
						lmin1(k) = 0
						q_S = "end of dissolution and of saturation of " // &
						      mineral_S(k)(1:len_trim(mineral_S(k))) // " at fc =" // y_S(1:len_trim(y_S))
					ELSEIF (lmin(k) == 0 .AND. lmin1(k) == 0 .AND. lmin0(k) == 1) THEN
						lmin1(k) = 0
						q_S = "end of saturation of " // &
						      mineral_S(k)(1:len_trim(mineral_S(k))) // " at fc =" // y_S(1:len_trim(y_S))
					END IF
				ELSEIF (syst_S == "o") THEN
					IF (lmin(k) == 1 .AND. lmin0(k) == 0) THEN
						q_S = "start of precipitation of " // &
						       mineral_S(k)(1:len_trim(mineral_S(k))) // &
						       " at fc =" // y_S(1:len_trim(y_S))
					ELSEIF (lmin(k) == 0 .AND. lmin0(k) == 1) THEN
						q_S = "end of precipitation of " // &
						      mineral_S(k)(1:len_trim(mineral_S(k))) // " at fc =" // &
					              y_S(1:len_trim(y_S)) // ": moles =" // x_S(1:len_trim(x_S))
					END IF
				END IF
				IF (q_S /= "") THEN
					OPEN (10,file=fichmin_S,position='APPEND')
						write(10,'(a)')  q_S(1:len_trim(q_S))
					CLOSE (10)
					q_S = ""
				END IF
			ENDDO
      
			IF (ncpt == 1 .OR. my_S(1:len_trim(my_S)) /= my0_S(1:len_trim(my0_S)) .OR. npasf == npasfich) THEN
				IF (unit_S == "molal")  ee = 1000.D0
				IF (unit_S == "molar")  ee = 1000000.D0 * dens / (1000.D0 + std)
				OPEN (10,file=fich_S,position= 'APPEND')
					IF (my_S == "") THEN
						write(10,'(a,a,5(g13.6,a),$)')  &
						"No_minerals",',',real(fc),",",real(mwev),",", &
						real(dens),",",real(-LOG(act(13)) / LOG(10.D0)),",",real(alc * ee),","
					ELSE
						write(10,'(a,a,5(g13.6,a),$)') my_S(2:len_trim(my_S)),',',real(fc),",",&
						real(mwev),",",real(dens),",",real(-LOG(act(13)) / LOG(10.D0)),",",real(alc * ee),","
					END IF
					DO i = 1 , 10
						IF (tot(i) > 0.) THEN
							s = 0.
							DO j = 1 , n
								s = s + molal(j) * kmat(i, j)
							ENDDO
							write(10,'(g13.6,a,$)') real(s * ee), ","
						END IF
					ENDDO
					write(10,*) real(std * ee)
				CLOSE (10)
           
				OPEN (10,file=molemin_S,position='APPEND')
					write(10,'(g13.6,a,$)') real(fc), ","
					IF (syst_S == "c") THEN
						DO i = 1 , nm - 1
							IF (min(i) >= 0.) THEN
								write(10,'(g13.6,a,$)') real(min(i)), ","
							ELSEIF (min(i) < 0.) THEN
								write(10,'("0., ",$)')
							END IF
						ENDDO
						IF (min(nm) >= 0.) THEN
							write(10,'(g13.6)') real(min(nm))
						ELSEIF (min(nm) < 0.) THEN
							write(10,'("0. ")')
						END IF
					ELSEIF (syst_S == "o") THEN
						DO i = 1 , nm - 1
							IF (minp(i) >= 0.) THEN
								write(10,'(g13.6,a,$)') real(minp(i)), ","
							ELSEIF (minp(i) < 0.) THEN
								write(10,'("0., ",$)')
							END IF
						ENDDO
						IF (minp(nm) >= 0.) THEN
							write(10,'(g13.6)') real(minp(nm))
						ELSEIF (minp(nm) < 0.) THEN
							write(10,'("0. ")')
						END IF
					END IF
				CLOSE (10)
			END IF
		END IF
		IF (stdmax > 0 .AND. std * ee >= stdmax) THEN
			IF (stockfich_S == "y")  CALL compact(lmin)
			call stop_simulation
		END IF

		IF (mwev > 0. .AND. kinvariant == 0 .AND. ABS(act(11) - act0(11)) / (act(11) + act0(11)) < .0000000001D0) THEN
			IF (syst_S == "c") THEN
				nu = 0
				DO k = 1 , nm
					IF (min(k) > 0. .AND. min(k) < min0(k))  nu = 1
				ENDDO
				IF (nu == 0) THEN
					IF (nminer == nbmin)  q_S = "invariant system / eutectic point / end"
					IF (nminer < nbmin)  q_S = "invariant system / pseudo-eutectic point / end"
				ELSEIF (nu == 1) THEN
					IF (nperitec == 0)  CALL peritec(lmin, kmat, wmin, tot, nperitec)
					IF (nperitec == 0) THEN
						IF (nminer == nbmin)  q_S = "invariant system / peritectic point / end"
						IF (nminer < nbmin)  q_S = "invariant system / pseudo-peritectic point / end"
					ELSEIF (nperitec == 1) THEN
						IF (nminer == nbmin)  q_S = "invariant system / peritectic / passing over"
						IF (nminer < nbmin)  q_S = "invariant system / pseudo-peritectic / passing over"
					END IF
				 END IF
			ELSEIF (syst_S == "o") THEN
				IF (nminer == nbmin)  q_S = "invariant system / eutectic / end"
				IF (nminer < nbmin)  q_S = "invariant system / pseudo-eutectic / end"
			END IF
			write(*,*) ; write(*,*) q_S(1:len_trim(q_S))
			IF (INDEX(q_S, "pseudo") /= 0) THEN
				q1_S = "maximum number of minerals allowed by the phase rule = "
				q2_S = "number of minerals in equilibrium with the invariant system = "
				write(*,*) q1_S(1:len_trim(q1_S)), nbmin
				write(*,*) q2_S(1:len_trim(q2_S)), nminer
			END IF
			IF (prnt_S /= "n" .AND. q_S /= q0_S) THEN
				write(9,*) ; write(9,*) q_S(1:len_trim(q_S))
				IF (INDEX(q_S, "pseudo") /= 0) THEN
					write(9,*) q1_S(1:len_trim(q1_S)), nbmin
					write(9,*) q2_S(1:len_trim(q2_S)), nminer
				END IF
			END IF
			IF (stockfich_S == "y" .AND. q_S /= q0_S) THEN
			  OPEN (10,file=fichmin_S,position='APPEND')
			  write(10,*) q_S(1:len_trim(q_S)) // " at fc = ",real(fc)
				IF (INDEX(q_S, "pseudo") /= 0) THEN					
					write(10,*) q1_S(1:len_trim(q1_S)), nbmin
					write(10,*) q2_S(1:len_trim(q2_S)), nminer
				END IF
				q0_S = q_S
			  CLOSE (10)
			END IF
			IF (INDEX(q_S, "end") /= 0) THEN
				IF (stockfich_S == "y")  CALL compact(lmin)
				call stop_simulation
			END IF
		ELSEIF (kinvariant == 0 .AND. ABS(act(11) - act0(11)) / (act(11) + act0(11)) > .0000000001D0) THEN
			nperitec = 0
			q_S = "";q0_S = "";q1_S = "";q2_S = ""
		END IF
  
		IF (npasi == npasimp)  npasi = 0
		IF (npasf == npasfich)  npasf = 0
		IF (my_S /= my0_S)  my0_S = my_S
		DO i = 0 , n
			mol0(i) = mol(i)
			molal0(i) = molal(i)
			gact0(i) = gact(i)
			act0(i) = act(i)
		ENDDO
		DO i = 1 , ntot
			tot0(i) = tot(i)
		ENDDO
		DO k = 1 , nm
			lmin0(k) = lmin(k)
			pai0(k) = pai(k)
			min0(k) = min(k)
			minp0(k) = minp(k)
		ENDDO
		fi0 = fi
		mwev0 = mwev
		nminer0 = nminer
		IF (kinvariant == 0 .AND. initdeseq == 0) THEN
			IF (inc_S == "auto") THEN
				xi = (51. - 8. * LOG(std * 1000.D0) / LOG(10.D0)) / 700.
				xi0 = xi
			ELSEIF (inc_S == "manu") THEN
				xi = xi0
			END IF
			mwev = mwev + mol(11) * xi
			tot(11) = tot(11) - 2. * mol(11) * xi
		END IF
	ENDDO


	2000 continue   !energie:
	write(*,*) "Free energy minimization"
	IF (kinvariant == 2) THEN
		IF (wmin(kinvar(1), ncm) > wmin(kinvar(2), ncm)) THEN
			ksupprim = kinvar(1)
		ELSE
			ksupprim = kinvar(2)
		END IF
	ELSEIF (kinvariant > 2) THEN
		gmin = 0
		DO ii = 1 , kinvariant
			DO i = 0 , n
				mol(i) = mol0(i)
			ENDDO
			DO i = 1 , ntot
				tot(i) = tot0(i)
			ENDDO
			mwev = mwev + mol(11) * xi
			tot(11) = tot(11) - 2 * mol(11) * xi
			DO k = 1 , nm
				lmin(k) = lmin0(k)
				min(k) = min0(k)
				minp(k) = minp0(k)
			ENDDO
			IF (kinvar(ii) > 0 .AND. kinvar(ii) /= kneuf .AND. min(kinvar(ii)) >= 0) THEN
				l = lmin(kinvar(ii)); m = min(kinvar(ii)); p = psol(kinvar(ii))
				lmin(kinvar(ii)) = 0; min(kinvar(ii)) = 0.;
				psol(kinvar(ii)) = 1.D+50
				write(*,*) "mineral removed : ", mineral_S(kinvar(ii))(1:len_trim(mineral_S(kinvar(ii))))
				ncmptinv = 0
				nw=1
				DO while(nw /= 0)
					ncmptinv = ncmptinv + 1
					DO i = 0 , n; gact0(i) = gact(i); ENDDO
					CALL actp(mol, gact, fi, temp)
					DO i = 0 , n; gact(i) = (gact0(i) + gact(i)) / 2.; ENDDO
					DO i = 0 , n
						molal(i) = mol(i) * mh2o / mol(11)
						act(i) = molal(i) * gact(i)
					ENDDO
					DO k = 1 , nm
						pai(k) = 1.
						DO i = 1 , ncm
							pai(k) = pai(k) * act(i) ** wmin(k, i)
						ENDDO
						IF (pai(k) >= psol(k)) THEN
							IF (min(k) >= 0.) THEN
								lmin(k) = 1
							ELSEIF (min(k) < 0.) THEN
								lmin(k) = 0
								min(k) = 0.
							END IF
						ELSEIF (pai(k) < psol(k)) THEN
							IF (min(k) <= 0.) THEN
								lmin(k) = 0
								min(k) = 0.
							ELSEIF (min(k) > 0.) THEN
								lmin(k) = 1
							END IF
						END IF
					ENDDO
					mineraux_S = ""
					nminer = 0
					DO k = 1 , nm
						IF (lmin(k) == 1) THEN
							nminer = nminer + 1
							mineraux_S = mineraux_S(1:len_trim(mineraux_S)) // "_" // mineral_S(k)(1:len_trim(mineral_S(k)))
						END IF
					ENDDO
					write(*,'(i4,t6,a)') ncmptinv,mineraux_S(1:len_trim(mineraux_S))
                                             
					DO i = 0 , n; mol1(i) = mol(i); ENDDO
					CALL reseq(nch, tot, mol, psol, psc, kmat, wmin, gact, act(11), ica, lmin, min, mineral_S, pkeq)
					DO i = 0 , n; molal(i) = mol(i) * mh2o / mol(11); ENDDO
					nw = 0
					DO i = 0 , n
						IF (ABS(mol1(i) - mol(i)) > pkmol)  nw = 1
					ENDDO
				ENDDO
                            
				g = 0.
				DO i = 0 , ncm
					IF (mol(i) > 0.)  g = g + mol(i) * (mu(i) + LOG(act(i)))
				ENDDO
				DO k = 1 , nm
					g = g + min(k) * mum(k)
				ENDDO
				write(*,'(a,g13.6,t25,$)') "g = ", g
				DO i = 1 , kinvariant
					IF (i /= ii) write(*,'(a,1x,$)') mineral_S(kinvar(i))(1:len_trim(mineral_S(kinvar(i))))
				ENDDO
				write(*,*) ""
				IF (g < gmin) THEN
					gmin = g
					ksupprim = kinvar(ii)
				END IF
				lmin(kinvar(ii)) = l; min(kinvar(ii)) = m; psol(kinvar(ii)) = p
			END IF
		ENDDO 
	END IF
	DO i = 0 , n
		mol(i) = mol0(i)
	ENDDO
	DO i = 1 , ntot
		tot(i) = tot0(i)
	ENDDO
	DO k = 1 , nm
		lmin(k) = lmin0(k)
		min(k) = min0(k)
		minp(k) = minp0(k)
	ENDDO
	mwev = mwev0
	nminer = nminer0
	mwev = mwev + mol(11) * xi
	tot(11) = tot(11) - 2. * mol(11) * xi
	kinvariant = 0
	lmin(ksupprim) = 0; min(ksupprim) = 0.; psol(ksupprim) = 1.D+50
	min0(ksupprim) = 0.
	write(*,*) "mineral definitely removed : ", mineral_S(ksupprim)(1:len_trim(mineral_S(ksupprim)))
	GOTO 500 !debut


end program evp




SUBROUTINE actp (mol, gact, fi, temp)
  use data_common
  implicit none
  real(kind=8), dimension (0:*)	:: mol, gact
  real(kind=8)	:: fi, temp
  
  double precision, external	:: temperature, j0, j1, g0, g1 ! functions
  
  integer, save	:: nc, na, nn, ndepact
  integer, save, allocatable, dimension (:)	:: nzc, nza
  real(kind=8), save	:: ap0, bp0
  real(kind=8), save, allocatable, dimension (:)	:: cat, ani, h
  real(kind=8), save, allocatable, dimension (:,:)	:: b0, b1, b2, c0, c, ta, tc, lc ,la
  real(kind=8), save, allocatable, dimension (:,:,:)	:: sc, sa, xi
  
  real(kind=8), allocatable, dimension (:,:)	:: ec, fc, xc, ea, fa, xa, pp, p, pf, qp, q, qf, cc, bf, b, bp
  real(kind=8), allocatable, dimension (:)	:: gc, ga, gn
  
  
  integer	:: alloc_res,i,j,k,ii,jj
  real (kind=8)	:: u, z, fj, w, v, f, s, co, aw, at, bt, ct, dt, et
  character (len=maxlen)	:: x_S
  
  allocate (cat(0:9), ani(0:11), h(0:3), stat=alloc_res)
  cat=0.; ani=0.; h=0.
  cat(1) = mol(1); cat(2) = mol(2); cat(3) = mol(3)
  cat(4) = mol(4); cat(5) = mol(22); cat(6) = mol(5)
  cat(7) = mol(18); cat(8) = mol(23); cat(9) = mol(13)
  ani(1) = mol(6); ani(2) = mol(7); ani(3) = mol(25)
  ani(4) = mol(15); ani(5) = mol(14); ani(6) = mol(12)
  ani(7) = mol(24); ani(8) = mol(19); ani(9) = mol(20)
  ani(10) = mol(21); ani(11) = mol(8)
  h(1) = mol(10); h(2) = mol(9); h(3) = mol(0)
  DO i = 1 , 9
    cat(i) = cat(i) * mh2o / mol(11)
  ENDDO
  DO i = 1 , 11
    ani(i) = ani(i) * mh2o / mol(11)
  ENDDO
  DO i = 1 , 3
    h(i) = h(i) * mh2o / mol(11)
  ENDDO
  IF (ndepact == 0) THEN
    OPEN (10,file="coefft4",action='READ') 
      read(10,*)  nc, na, nn
      allocate (nzc(0:nc), nza(0:na), stat=alloc_res)
      allocate (b0(0:nc, 0:na), b1(0:nc, 0:na), b2(0:nc, 0:na), c0(0:nc, 0:na), stat=alloc_res) 
      allocate (sc(0:nc, 0:nc, 0:na), sa(0:na, 0:na, 0:nc), stat=alloc_res)
      allocate (tc(0:nc, 0:nc), ta(0:na, 0:na), stat=alloc_res)
      allocate (lc(0:nn, 0:nc), la(0:nn, 0:na), xi(0:nn, 0:nc, 0:na), stat=alloc_res)
	  nzc=0; nza=0; b0=0.; b1=0.; b2=0.; c0=0.
	  sc=0.; sa=0.; tc=0.; ta=0.; lc=0.; la=0.; xi=0.
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
  bp0 = 1.2D0
  allocate (ec(0:nc, 0:nc), fc(0:nc, 0:nc), xc(0:nc, 0:nc), ea(0:na, 0:na), fa(0:na, 0:na), xa(0:na, 0:na), stat=alloc_res)
  allocate (pp(0:nc, 0:nc), p(0:nc, 0:nc), pf(0:nc, 0:nc), qp(0:na, 0:na), q(0:na, 0:na), qf(0:na, 0:na), stat=alloc_res)
  allocate (cc(0:nc, 0:na), bf(0:nc, 0:na), b(0:nc, 0:na), bp(0:nc, 0:na), stat=alloc_res)
  allocate (gc(0:nc), ga(0:na), gn(0:nn), stat=alloc_res)
  ec=0.; fc=0.; xc=0.; ea=0.; fa=0.; xa=0.; pp=0.; p=0.; pf=0.; qp=0.; q=0.; qf=0.; gc=0.; ga=0.; gn=0.
  
  ndepact = 1

  u = 0.; z = 0.
  DO i = 1 , nc; u = u + cat(i) * nzc(i) ** 2; z = z + cat(i) * nzc(i); ENDDO
  DO j = 1 , na; u = u + ani(j) * nza(j) ** 2; z = z + ani(j) * nza(j); ENDDO
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
                    / fi ** 2 / 4 - ec(i, j) / fi
        ec(j, i) = ec(i, j); fc(j, i) = fc(i, j)
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
        fa(i, j) = (xa(i, j) * j1(xa(i, j)) - xa(i, i) * j1(xa(i, i)) / 2. - xa(j, j) * j1(xa(j, j)) / 2.) &
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
      IF (nzc(i) == 1 .OR. nza(j) == 1) v = fj * 2.
      bf(i, j) = b0(i, j) + b1(i, j) * EXP(-v) + b2(i, j) * EXP(-w)
      b(i, j) = b0(i, j) + b1(i, j) * (g0(v)) + b2(i, j) * (g0(w))
      bp(i, j) = b1(i, j) * (g1(v)) / fi + b2(i, j) * (g1(w)) / fi
    ENDDO
  ENDDO
  f = -ap0 * (fj / (1. + bp0 * fj) + 2. / bp0 * LOG(1. + bp0 * fj))
  DO i = 1 , nc; DO j = 1 , na; f = f + cat(i) * ani(j) * bp(i, j); ENDDO; ENDDO
  DO i = 1 , nc - 1; DO j = i + 1 , nc; f = f + cat(i) * cat(j) * pp(i, j); ENDDO; ENDDO
  DO i = 1 , na - 1; DO j = i + 1 , na; f = f + ani(i) * ani(j) * qp(i, j); ENDDO; ENDDO
  DO ii = 1 , nc
    u = nzc(ii) ** 2 * f
    DO j = 1 , na; u = u + ani(j) * (b(ii, j) * 2. + z * cc(ii, j)); ENDDO
      DO i = 1 , nc
        IF (i /= ii) THEN
          v = 0.
          DO j = 1 , na; v = v + ani(j) * sc(ii, i, j); ENDDO
          u = u + cat(i) * (p(ii, i) * 2. + v)
        END IF
      ENDDO
      DO i = 1 , na - 1; DO j = i + 1 , na; u = u + ani(i) * ani(j) * sa(i, j, ii); ENDDO; ENDDO
      DO i = 1 , nc; DO j = 1 , na; u = u + cat(i) * ani(j) * cc(i, j) * nzc(ii); ENDDO; ENDDO
      DO i = 1 , nn
        u = u + h(i) * lc(i, ii) * 2.
      ENDDO
      DO k = 1 , nn
        DO j = 1 , na
          u = u + h(k) * ani(j) * xi(k, ii, j)
        ENDDO
      ENDDO
    gc(ii) = EXP(u)
  ENDDO 
  DO jj = 1 , na
    u = nza(jj) ** 2 * f
    DO i = 1 , nc; u = u + cat(i) * (b(i, jj) * 2. + z * cc(i, jj)); ENDDO
      DO i = 1 , na
        IF (i /= jj) THEN
          v = 0.
          DO j = 1 , nc; v = v + cat(j) * sa(jj, i, j); ENDDO
          u = u + ani(i) * (q(jj, i) * 2. + v)
        END IF
      ENDDO
      DO i = 1 , nc - 1; DO j = i + 1 , nc; u = u + cat(i) * cat(j) * sc(i, j, jj); ENDDO; ENDDO
      DO i = 1 , nc; DO j = 1 , na; u = u + cat(i) * ani(j) * cc(i, j) * nza(jj); ENDDO; ENDDO
      DO j = 1 , nn
        u = u + h(j) * la(j, jj)
      ENDDO
      DO k = 1 , nn
        DO i = 1 , nc
          u = u + h(k) * cat(i) * xi(k, i, jj)
        ENDDO
      ENDDO
    ga(jj) = EXP(u)
  ENDDO 
  DO k = 1 , nn
    u = 0.
    DO i = 1 , nc
      u = u + cat(i) * lc(k, i) * 2.
    ENDDO
    DO j = 1 , na
      u = u + ani(j) * la(k, j) * 2
    ENDDO
    DO i = 1 , nc
      DO j = 1 , na
        u = u + cat(i) * ani(j) * xi(k, i, j)
      ENDDO
    ENDDO
    gn(k) = EXP(u)
  ENDDO
  u = -ap0 * fi ** 1.5D0 / (1. + bp0 * fj)
  DO i = 1 , nc; DO j = 1 , na; u = u + cat(i) * ani(j) * (bf(i, j) + z * cc(i, j)); ENDDO; ENDDO
  DO i = 1 , nc - 1
    DO j = i + 1 , nc
      v = 0.
      DO k = 1 , na; v = v + ani(k) * sc(i, j, k); ENDDO
      u = u + cat(i) * cat(j) * (pf(i, j) + v)
    ENDDO
  ENDDO
  DO i = 1 , na - 1
    DO j = i + 1 , na
      v = 0.
      DO k = 1 , nc; v = v + cat(k) * sa(i, j, k); ENDDO
      u = u + ani(i) * ani(j) * (qf(i, j) + v)
    ENDDO
  ENDDO
  DO k = 1 , nn
    DO i = 1 , nc
      u = u + h(k) * cat(i) * lc(k, i)
    ENDDO
  ENDDO
  DO k = 1 , nn
    DO j = 1 , na
      u = u + h(k) * ani(j) * la(k, j)
    ENDDO
  ENDDO
  DO k = 1 , nn
    DO i = 1 , nc
      DO j = 1 , na
        u = u + h(k) * cat(i) * ani(j) * xi(k, i, j)
      ENDDO
    ENDDO
  ENDDO
  s = 0.
  DO i = 1 , nc; s = s + cat(i); ENDDO
  DO j = 1 , na; s = s + ani(j); ENDDO
  co = 1. + 2. * u / s; aw = EXP(-s * co / mh2o)
  gact(0) = gn(3); gact(1) = gc(1); gact(2) = gc(2); gact(3) = gc(3)
  gact(4) = gc(4); gact(22) = gc(5); gact(5) = gc(6)
  gact(18) = gc(7); gact(23) = gc(8); gact(13) = gc(9)
  gact(6) = ga(1); gact(7) = ga(2); gact(25) = ga(3)
  gact(15) = ga(4); gact(14) = ga(5); gact(12) = ga(6)
  gact(24) = ga(7); gact(19) = ga(8); gact(20) = ga(9)
  gact(21) = ga(10); gact(8) = ga(11)
  gact(10) = aw*aw*gn(1)**log(10.D0); gact(9) = gn(2)
  gact(16) = 1D0; gact(17) = 1D0; gact(11) = aw / mh2o
  ndepact = 1
  
  deallocate (ec, fc, xc, ea, fa, xa, stat=alloc_res)
  deallocate (pp, p, pf, qp, q, qf, stat=alloc_res)
  deallocate (cc, bf, b, bp, stat=alloc_res)
  deallocate (gc, ga, gn, stat=alloc_res)

END SUBROUTINE actp



SUBROUTINE compact (lmin)
    use data_common
	implicit none
	integer, dimension (0:*)	:: lmin
	
	real(kind=8)	:: fc
	real(kind=8), allocatable, dimension(:)	:: q_min
	integer	:: i,alloc_res, nbmin_comp, nlast
	integer, parameter	:: maxrec=4096
	character (len=maxlen)	:: x_S
	character (len=maxlen), allocatable, dimension (:)	:: nom_min
	
	allocate(q_min(0:nm),nom_min(0:nm), stat=alloc_res)
		
    PRINT *, ' '
    PRINT *, "Compacting mineral file"
    PRINT *, ' '
    DO i = 1 , nm
        lmin(i) = 0
    ENDDO
    nbmin_comp=0
    
    open(10,file=molemin_S,action='READ',recl=maxrec)
    open(11,status='SCRATCH',recl=maxrec)
    read(10,*) x_S,(nom_min(i),i=1,nm)
    do 
    	read(10,*,end=50) fc,(q_min(i),i=1,nm)
    	do i=1,nm
    		if(q_min(i).gt.0.) then
    			if(lmin(i).eq.0) then
    				lmin(i)=1
    				nlast=max(nlast,i)
    				nbmin_comp=nbmin_comp+1
    			endif
    		endif
    	enddo
    enddo
	50 continue
    
    rewind(10)
    read(10,*) x_S,(nom_min(i),i=1,nm)
    
    write(11,'(A,'',''$)') (X_S(1:len_trim(x_S)))
    do i=1,nm
    	if(lmin(i).eq.1) then
    		if(i.lt.nlast) then
    			write(11,'(A,'',''$)') nom_min(i)(1:len_trim(nom_min(i)))
    		else
    			write(11,'(A)') nom_min(i)(1:len_trim(nom_min(i)))
    		endif
    	endif	
    enddo    
    
    do 
    	read(10,*,end=60) fc,(q_min(i),i=1,nm)
    	write(11,'(g14.8,'',''$)') fc
    	do i=1,nm
    		if(lmin(i).eq.1) then
    			if(i.lt.nlast) then
    				write(11,'(g14.8,'',''$)') q_min(i)
    			else
    				write(11,'(g14.8)') q_min(i)
    			endif
    		endif		
    	enddo
    enddo
	60 continue
   
    close(10)
    open(10,file=molemin_S,action='WRITE',recl=maxrec)
    rewind (11)
    read(11,*) x_S,(nom_min(i),i=1,nbmin_comp)
    write(10,'(A,'',''$)') (X_S(1:len_trim(x_S)))
    
    do i=1,nbmin_comp
    		if(i.lt.nbmin_comp) then
    			write(10,'(A,'',''$)') nom_min(i)(1:len_trim(nom_min(i)))
    		else
    			write(10,'(A)') nom_min(i)(1:len_trim(nom_min(i)))
    		endif
    enddo
    
    do 
    	read(11,*,end=70) fc,(q_min(i),i=1,nbmin_comp)
    	write(10,'(g14.8,'',''$)') fc
    	do i=1,nbmin_comp
    			if(i.lt.nbmin_comp) then
    				write(10,'(g14.8,'',''$)') q_min(i)
    			else
    				write(10,'(g14.8)') q_min(i)
    			endif
    	enddo
    enddo
	70 continue
    
	close (10)
	close (11)
	deallocate(q_min, nom_min, stat=alloc_res)

END SUBROUTINE compact



SUBROUTINE densite (mol, kmat, nch, dens)
	use data_common
	implicit none
	
	integer, dimension (0:*)	:: nch
	real(kind=8), dimension (0:n,0:n)	:: kmat
	real(kind=8), dimension (0:*)	:: mol
	real(kind=8)	:: dens
	
	integer	:: ncdens,nadens,alloc_res,i,j
	integer, allocatable, dimension(:)	:: ic, ia
	real(kind=8), allocatable, dimension(:,:), save	:: au,bu
	real(kind=8), allocatable, dimension(:,:)	:: s
	real(kind=8), allocatable, dimension(:)	:: cat, ani
	real(kind=8)	:: u,v
	character(len=maxlen)	:: x_S

    ncdens = 5; nadens = 5
    allocate (s(0:ncdens, 0:nadens), cat(0:ncdens), ani(0:nadens), ic(0:ncdens), ia(0:nadens), stat=alloc_res)
    s=0.; cat=0.; ani=0.; ic=0; ia=0

    DO i = 1 , 8
        IF (i <= ncdens)  cat(i) = mol(i) / mol(11) * mh2o
        IF (i > ncdens)  ani(i - ncdens) = mol(i) / mol(11) * mh2o
    ENDDO
    ani(4) = mol(15) / mol(11) * mh2o
    ani(5) = mol(14) / mol(11) * mh2o
 
    DO i = 1 , ncdens; ic(i) = nch(i); ENDDO
    DO i = 1 , 3; ia(i) = -nch(i + 5); ENDDO
    ia(4) = -nch(15); ia(5) = -nch(14)
    IF (ncpt == 1) THEN
        allocate (au(0:ncdens, 0:nadens), bu(0:ncdens, 0:nadens), stat=alloc_res)
        OPEN (10,file="densite",action='READ')
          DO i = 1 , 5
              DO j = 1 , 5
                  read(10,*)  x_S, u, v, au(i, j), bu(i, j)
              ENDDO
          ENDDO
        CLOSE (10)
    END IF
    dens = 1D0; u = 0
    DO j = 1 , nadens; u = u + ani(j) * ia(j); ENDDO
    DO i = 1 , ncdens
      DO j = 1 , nadens
        s(i, j) = INT((ic(i) + ia(j)) / 2) * cat(i) * ani(j) / u
        dens = dens + au(i, j) * s(i, j) + bu(i, j) * s(i, j) ** 2
      ENDDO
    ENDDO
    deallocate (s, cat, ani, ic, ia, stat=alloc_res)
END SUBROUTINE densite

DOUBLE PRECISION FUNCTION g0 (x)
	real(kind=8)	:: x
    g0 = 2. * (1. - (1. + x) * EXP(-x)) / x ** 2
END FUNCTION g0

DOUBLE PRECISION FUNCTION g1 (x)
	real(kind=8)	:: x
    g1 = -2. * (1. - (1. + x + x ** 2 / 2.) * EXP(-x)) / x ** 2
END FUNCTION g1


SUBROUTINE invar (aw, wmin, mineral_S, lmin, psol, psc, kinvar)
    use data_common
	implicit none
	real(kind=8)	:: aw
	real(kind=8), dimension (0:*)	:: psol,psc
	real(kind=8), dimension (0:nm,0:ncm)	:: wmin
	integer, dimension (0:*)	:: lmin
	integer, dimension (0:*) :: kinvar
	character (len=*), dimension(0:*)	:: mineral_S


    real (kind=8)	:: alea1,alea2,swap,u,ah2o
	integer	:: alloc_res, ninvar, l,i,j,k,n1,n2,det,nz,ii,jj,kk,det1,n3,n4,i1,i2, ncond  
	real (kind=8), allocatable, dimension (:)	:: psminv, psminvar, tt4
	real (kind=8), allocatable, dimension (:,:)	:: t0, t1, t2, t3, t4, winv
	integer, allocatable, dimension (:)	:: kinv
	character(len=maxlen), allocatable, dimension (:)	:: minv_S, minvar_S

    kinvariant = 0; ncm = 15
    ninvar = nbmin + 3
    allocate (kinv(0:ninvar), minv_S(0:ninvar), psminv(0:ninvar), winv(0:ninvar, 0:ncm), stat=alloc_res)
    allocate (minvar_S(0:ninvar), psminvar(0:ninvar), stat=alloc_res)
    allocate (t0(0:ninvar, 0:ninvar), t1(0:ninvar, 0:ninvar), stat=alloc_res)
    allocate (t2(0:ninvar, 0:ninvar), t3(0:ninvar, 0:ncm), t4(0:ncm, 0:ncm), tt4(0:ncm), stat=alloc_res)
    kinv=0; psminv=0.; winv=0.; psminvar=0.; t0=0;; t1=0.; t2=0.; t3=0; t4=0.; tt4=0.

    DO k = 1 , 3
        psminv(k) = LOG(psc(k)) / LOG(10.D0)
    ENDDO
    winv(1, 13) = 1.; winv(1, 12) = 1.; winv(1, 15) = -1.
    winv(2, 13) = 1.; winv(2, 11) = -1.; winv(2, 14) = 1.
    winv(3, 13) = 1.; winv(3, 11) = 1.; winv(3, 15) = -1.

    n1 = 3
    DO k = 1 , nm
      IF (lmin(k) == 1) THEN
        n1 = n1 + 1
        kinv(n1) = k
        minv_S(n1) = mineral_S(k)
        psminv(n1) = LOG(psol(k)) / LOG(10.D0)
        DO j = 1 , ncm
          winv(n1, j) = wmin(k, j)
        ENDDO
        swap= winv(n1, 11)
        winv(n1, 11)=winv(n1, 15)
        winv(n1, 15)=swap
      END IF
    ENDDO

    DO i = 1 , n1
      DO j = i , n1
        t1(i, j) = 0.
        DO k = 1 , ncm - 1
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
            IF (ABS(t1(i, j)) < epsilon)  t1(i, j) = 0.
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
          IF (i /= kk) THEN
            ii = ii + 1
            jj = 0
            DO j = 1 , n1
              IF (j /= kk) THEN
                jj = jj + 1
                t2(ii, jj) = t0(i, j)
              END IF
            ENDDO
          END IF
        ENDDO

        DO k = 2 , n2
          DO i = k , n2
            IF (t2(i, k - 1) /= 0) THEN
              u = t2(k - 1, k - 1) / t2(i, k - 1)
              DO j = k , n2
                t2(i, j) = t2(k - 1, j) - t2(i, j) * u
                IF (ABS(t2(i, j)) < epsilon) t2(i, j) = 0.
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
          DO j = 1 , ncm
            t3(n3, j) = winv(kk, j)
          ENDDO
        END IF
      ENDDO

      n4 = ncm
      DO j = ncm , 1 , -1
        u = 0
        DO i = 1 , n3
          u = u + t3(i, j) ** 2
        ENDDO
        IF (u == 0) THEN
          DO k = j + 1 , n4
            DO i = 1 , n3
              t3(i, k - 1) = t3(i, k)
            ENDDO
          ENDDO
          n4 = n4 - 1
        END IF
      ENDDO

      DO i = 1 , n4
        DO j = i , n4
           t4(i,j)=0.
           DO k = 1, n4
              t4(i,j) = t4(i,j) + t3(k,i) * t3(k,j)
              t4(j,i) = t4(i,j)
           ENDDO
        ENDDO
      ENDDO
      DO i = 1, n4
         tt4(i) = 0.
         DO k =1 , n3
            tt4(i) = tt4(i) + t3(k,i) * psminvar(k)
         ENDDO
      ENDDO
  

      DO k = 2 , n4
        DO i = k , n4
          IF (ABS(t4(i, k - 1)) > epsilon) THEN
            u = t4(k - 1, k - 1) / t4(i, k - 1)
            DO j = k , n4
              t4(i, j) = t4(k - 1, j) - t4(i, j) * u
              IF (ABS(t4(i, j)) < epsilon)  t4(i, j) = 0.
            ENDDO
            tt4(i) = tt4(k - 1) - tt4(i) * u
          END IF
        ENDDO
      ENDDO



      ah2o = 10 ** (tt4(n4) / t4(n4, n4))
      kinvariant = n3
      DO i = kinvariant , 1 , -1
          IF (kinvar(i) == 0) THEN
              DO k = 1 , kinvariant - 1
                  kinvar(k) = kinvar(k + 1)
              ENDDO
              kinvariant = kinvariant - 1
          END IF
      ENDDO
      write(*,'(a,$)') "invariant system constrained by: "
      DO k = 1 , kinvariant
        write(*,'(a,1X,$)') " "//mineral_S(kinvar(k))(1:len_trim(mineral_S(kinvar(k))))
      ENDDO
      write(*,*) " "
      write(*,*) "invariant aH2O = ", ah2o
      write(*,*) "simulation aH2O    = ", aw
      IF (stockfich_S == "y") THEN
	OPEN (10,file=fichmin_S,position='APPEND')
	    write(10,'(a,$)') "invariant system constrained by : "
        DO k = 1 , kinvariant
          write(10,'(a,1X,$)') " "//mineral_S(kinvar(k))(1:len_trim(mineral_S(kinvar(k))))
        ENDDO
        write(10,*) " "
        write(10,'(a,g13.6,a,g13.6)') "invariant aH2O = ", real(ah2o),"      simulation aH2O    = ", real(aw)
	CLOSE (10)
      END IF
    END IF
    deallocate (kinv, minv_S, psminv, winv, stat=alloc_res)
    deallocate (minvar_S, psminvar, stat=alloc_res)
    deallocate (t0, t1, stat=alloc_res)
    deallocate (t2, t3, t4, tt4, stat=alloc_res)
END SUBROUTINE invar


DOUBLE PRECISION FUNCTION j0 (x)
    implicit none
    real(kind=8)	:: ya,yb,yc,yd,x
    ya = 4.581; yb = -.7237; yc = -.012; yd = .528
    j0 = x / (4. + ya * x ** yb * EXP(yc * x ** yd))
END FUNCTION j0

DOUBLE PRECISION FUNCTION j1 (x)
    implicit none
    real(kind=8)	:: ya,yb,yc,yd,x
    ya = 4.581; yb = -.7237; yc = -.012; yd = .528
    j1 = (4. + ya * x ** yb * (1. - yb - yc * yd * x ** yd) * EXP(yc * x ** yd)) / (4. + ya * x ** yb * EXP(yc * x ** yd)) ** 2
END FUNCTION j1




SUBROUTINE peritec (lmin, kmat, wmin, tot, nperitec)
    use data_common
	implicit none
	real(kind=8), dimension (0:nm,0:ncm)	:: wmin
	real(kind=8), dimension (0:*)	:: tot
	integer, dimension (0:*)	:: lmin
	real(kind=8), dimension (0:n,0:n)	:: kmat
	integer	:: nperitec


	real (kind=8)	:: u,alea1,alea2,s,swap
	integer	:: nt,k,i,j,l,ni,nj,alloc_res,ncomp,i1,i2,ii,nresol
	real (kind=8), allocatable, dimension (:)	:: tt, tt0, xx
	real (kind=8), allocatable, dimension (:,:)	:: t, t0
	character(len=maxlen)	:: x_S

    nt = ntot - 1
    allocate (t(0:nt, 0:nt), tt(0:nt), t0(0:nt, 0:nt), tt0(0:nt), xx(0:nt), stat=alloc_res)
    t=0.; tt=0.; t0=0.; tt0=0.; xx=0

    j = 0
    DO k = 1 , nm
      IF (lmin(k) == 1) THEN
        j = j + 1
        DO i = 1 , nt - 1
          t0(i, j) = wmin(k, i)
        ENDDO
        s = 0.
        DO i = 1 , ncm
            s = s + wmin(k, i) * kmat(11, i)
        ENDDO
        t0(11, j) = s
      END IF
    ENDDO
    j = j + 1
    t0(nt, j) = 2
    nj = j
    DO i = 1 , nt
      tt0(i) = tot(i)
    ENDDO
    
    DO i = 1 , nj
        DO j = i , nj
           t(i,j)=0.
           DO k = 1, nt
              t(i,j) = t(i,j) + t0(k,i) * t0(k,j)
              t(j,i) = t(i,j)
           ENDDO
        ENDDO
      ENDDO
      DO i = 1, nj
         tt(i) = 0.
         DO k =1 , nt
            tt(i) = tt(i) + t0(k,i) * tt0(k)
         ENDDO
    ENDDO
      
    DO k = 2 , nj
      DO i = k , nj
        IF (ABS(t(i, k - 1)) > epsilon) THEN
          u = t(k - 1, k - 1) / t(i, k - 1)
          DO j = k , nj
            t(i, j) = t(k - 1, j) - t(i, j) * u
          ENDDO
          tt(i) = tt(k - 1) - tt(i) * u
        END IF
      ENDDO
    ENDDO
    xx(nj) = tt(nj) / t(nj, nj)
    DO i = nj - 1 , 1 , -1
      s = 0.
      DO j = i + 1 , nj
        s = s + t(i, j) * xx(j)
      ENDDO
      xx(i) = (tt(i) - s) / t(i, i)
    ENDDO
    nperitec = 0
    DO i = 1 , nj - 1
      IF (xx(i) < 0.)  nperitec = 1
    ENDDO
    deallocate (t, tt, t0, tt0, xx, stat=alloc_res)
END SUBROUTINE peritec





SUBROUTINE reseq (nch, tot, mol, psol, psc, kmat, wmin, gact, aw, ica, lmin, min, mineral_S, pkeq)
    use data_common
	implicit none
	real(kind=8)	:: aw,pkeq
	real(kind=8), dimension (0:*)	::tot,mol,gact,min,psol,psc
	real(kind=8), dimension (0:nm,0:ncm)	:: wmin
	integer, dimension (0:*)	:: nch,lmin,ica
	real(kind=8), dimension (0:n,0:n)	:: kmat
	character (len=*), dimension (0:*)	:: mineral_S
	
	character(len=30)	:: d_and_t_S
	real (kind=8)	:: s,sh,p,u
	integer	:: nconv,nt,k,i,j,l,ni,nj,alloc_res,kk,nu
	real (kind=8), allocatable, dimension (:)	:: zz,xx
	real (kind=8), allocatable, dimension (:,:)	:: z

    allocate (z(0:n + nminer, 0:n + nminer), zz(0:n + nminer), xx(0:n + nminer), stat=alloc_res)
    z=0.; zz=0.; xx=0.
    
    nconv = 1
    ncm = 15; nt = n
    nu=1
    DO while (nu /=0)
        DO i = 1 , nt
            DO j = 1 , nt
                z(i, j) = 0.
            ENDDO
            zz(i) = 0.
        ENDDO
        DO i = 1 , 12
            DO j = 1 , n
                IF (mol(j) /= 0.)  z(i, j) = kmat(i, j)
            ENDDO
        ENDDO
              
        DO i = 13 , n
            kmat(i, 0) = 0.
            DO j = 1 , n
                IF (j /= 11)  kmat(i, 0) = kmat(i, 0) + kmat(i, j)
            ENDDO
        ENDDO
        DO i = 13 , n
            p = 1.; u = 0.
            DO j = 1 , n
                IF (mol(j) /= 0. .AND. j /= 11) THEN
                    z(i, j) = kmat(i, j) / mol(j)
                    p = p * gact(j) ** kmat(i, j)
                    u = u + kmat(i, j) * LOG(mol(j))
                ELSEIF (j == 11) THEN
                    z(i, j) = -kmat(i, 0) / mol(j)
                    p = p * aw ** kmat(i, j)
                    u = u - kmat(i, 0) * LOG(mol(j))
                ELSEIF (mol(j) == 0.) THEN
                    z(i, j) = 0
                END IF
            ENDDO
            p = p * mh2o ** kmat(i, 0)
            zz(i) = LOG(psc(i - 12)) - LOG(p) - u
        ENDDO
             
        l = 0
        DO k = 1 , nm
            IF (lmin(k) == 1) THEN
                l = l + 1
                ica(n + l) = 1
                wmin(k, 0) = 0.
                p = 1.; u = 0.
                DO j = 1 , ncm
                    IF (j /= 11)  wmin(k, 0) = wmin(k, 0) + wmin(k, j)
                ENDDO
                DO j = 1 , ncm
                    IF (j /= 11 .AND. mol(j) > 0.) THEN
                        z(n + l, j) = wmin(k, j) / mol(j)
                        p = p * gact(j) ** wmin(k, j)
                        u = u + wmin(k, j) * LOG(mol(j))
                    ELSEIF (j == 11) THEN
                        z(n + l, j) = -wmin(k, 0) / mol(j)
                        p = p * aw ** wmin(k, j)
                        u = u - wmin(k, 0) * LOG(mol(j))
                    END IF
                ENDDO
                p = p * mh2o ** wmin(k, 0)
                zz(n + l) = LOG(psol(k)) - LOG(p) - u
                    
                sh = 0.
                DO j = 1 , ncm
                    sh = sh + wmin(k, j) * kmat(11, j)
                ENDDO
                DO i = 1 , 10
                    z(i, n + l) = wmin(k, i)
                ENDDO
                z(11, n + l) = sh
            END IF
        ENDDO
        nt = n + l
            
        DO i = 1 , 10
            u = 0.
            DO j = 1 , n
                u = u + kmat(i, j) * mol(j)
            ENDDO
            DO k = 1 , nm
                IF (lmin(k) == 1) THEN
                    u = u + min(k) * wmin(k, i)
                END IF
            ENDDO
            zz(i) = tot(i) - u
        ENDDO
        u = 0.
        DO j = 1 , n
            u = u + kmat(11, j) * mol(j)
        ENDDO
        DO j = 1 , ncm
          DO k = 1 , nm
              u = u + wmin(k, j) * kmat(11, j) * min(k)
          ENDDO
        ENDDO
     
        zz(11) = tot(11) - u
        u = 0.
        DO j = 1 , n
            u = u + z(12, j) * mol(j)
        ENDDO
        zz(12) = tot(12) - u
            
        DO k = 1 , 10
            IF (tot(k) == 0.) THEN
                ica(k) = 0
                DO j = k + 1 , n
                    IF (kmat(k, j) /= 0.)  ica(j) = 0
                ENDDO
            END IF
        ENDDO
             
        ni = nt; nj = nt
        DO k = nt , 1 , -1
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
                IF (ABS(z(i, k - 1)) > epsilon) THEN
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

        DO k = 1 , nt
            IF (ica(k) == 0) THEN
                DO i = ni , k , -1
                    xx(i + 1) = xx(i)
                ENDDO
                xx(k) = 0.
                ni = ni + 1
            END IF
        ENDDO
        nu=1
        DO while(nu /= 0)
            DO i = 1 , n
                IF (ica(i) == 1) THEN
                    nu = 0
                    IF (mol(i) + xx(i) / nconv < 0.) THEN
                        nconv = nconv + 1
                        nu = 1
                        EXIT
                    END IF
                END IF
            ENDDO
            IF (nconv >=  max_conv) THEN
                write(*,*)
                write(*,*) "The equation system diverges: end of simulation"
                IF (stockfich_S == "y")  THEN
					open(10,file=fichmin_S,position='APPEND')
					write(10,*)
					write(10,*) "The equation system diverges: end of simulation"
					close(10)
					CALL compact(lmin)
				END IF
				call stop_simulation
            END IF
        ENDDO
        DO i = 1 , n
            mol(i) = mol(i) + xx(i) / nconv
        ENDDO
        i = n
        DO k = 1 , nm
            IF (lmin(k) == 1) THEN
                i = i + 1
                min(k) = min(k) + xx(i) / nconv
            END IF
        ENDDO
        nu = 0
        DO i = 1 , nt
            IF (ica(i) == 1) THEN
                IF (ABS(xx(i)) > pkeq)  nu = 1
            END IF
        ENDDO
    ENDDO
    deallocate (z, zz, xx, stat=alloc_res)
END SUBROUTINE reseq



SUBROUTINE reseqinv (lmin, wmin, kmat, min, kinvar, xinv)
    use data_common
	implicit none
	real(kind=8)	:: xinv
	real(kind=8), dimension (0:*)	:: min
	real(kind=8), dimension (0:nm,0:ncm)	:: wmin
	integer, dimension (0:*)	:: lmin,kinvar
	real(kind=8), dimension (0:n,0:n)	:: kmat

    integer, save	:: ninv
    integer	:: nt,i,j,k,alloc_res,ncomp,ii,ncond,i1,i2,kk,nmin
    real(kind=8)	:: hmin, s, swap, u, alea1, alea2
    real(kind=8), allocatable, dimension (:)	:: tt, tt0, xx
    real(kind=8), allocatable, dimension (:,:)	:: t, t0
	character(len=maxlen)	:: x_S
	character(len=30)	:: d_and_t_S
    
    nt = ntot - 1
    allocate (t(0:nt, 0:nt), tt(0:nt), t0(0:nt, 0:nt), tt0(0:nt), xx(0:nt), stat=alloc_res)
    t=0.; tt=0.; t0=0.; tt0=0.; xx=0.
    
    IF (ninv == 0) THEN
        hmin = 1000.D0
        DO k = 1 , nm
            DO kk = 1 , kinvariant
                IF (k == kinvar(kk) .AND. min(k) > 0.) THEN
                    s = 0.
                    DO i = 1 , 15
                        s = s + wmin(k, i) * kmat(nt, i)
                    ENDDO
                    IF (s > 0.) THEN
                        IF (s * min(k) < hmin)  hmin = s * min(k)
                    END IF
                END IF
            ENDDO
        ENDDO
        xinv = hmin / 100.
        IF (xinv <= 0.) THEN
			call stop_simulation
		END IF
    END IF
    ninv = 1
    j = 0
    DO k = 1 , nm
      IF (lmin(k) == 1) THEN
        DO kk = 1 , kinvariant
            IF (k == kinvar(kk)) THEN
                j = j + 1
                DO i = 1 , nt - 1
                  t0(i, j) = wmin(k, i)
                ENDDO
                s = 0.
                DO i = 1 , 15
                    s = s + wmin(k, i) * kmat(nt, i)
                ENDDO
                t0(nt, j) = s
            END IF
        ENDDO
      END IF
    ENDDO
    nmin = j

    DO i = 1 , nt - 1
        tt(i) = 0.
        DO k = 1 , nm
            DO kk = 1 , kinvariant
                IF (k == kinvar(kk)) THEN
                    tt0(i) = tt0(i) + min(k) * wmin(k, i)
                END IF
            ENDDO
        ENDDO
    ENDDO
    tt0(nt) = 0.
    DO k = 1 , nm
        DO kk = 1 , kinvariant
            IF (k == kinvar(kk)) THEN
                DO i = 1 , 15
                    tt0(nt) = tt0(nt) + min(k) * wmin(k, i) * kmat(11, i)
                ENDDO
            END IF
        ENDDO
    ENDDO
    tt0(nt) = tt0(nt) - xinv

    DO i = 1 , nmin
        DO j = i , nmin
           t(i,j)=0.
           DO k = 1, nt
              t(i,j) = t(i,j) + t0(k,i) * t0(k,j)
              t(j,i) = t(i,j)
           ENDDO
        ENDDO
    ENDDO
    DO i = 1, nmin
         tt(i) = 0.
         DO k =1 , nt
            tt(i) = tt(i) + t0(k,i) * tt0(k)
         ENDDO
    ENDDO


    DO k = 2 , nmin
      DO i = k , nmin
        IF (ABS(t(i, k - 1)) > epsilon ) THEN
          u = t(k - 1, k - 1) / t(i, k - 1)
          DO j = k , nmin
            t(i, j) = t(k - 1, j) - t(i, j) * u
          ENDDO
          tt(i) = tt(k - 1) - tt(i) * u
        END IF
      ENDDO
    ENDDO

    xx(nmin) = tt(nmin) / t(nmin, nmin)
    DO i = nmin - 1 , 1 , -1
      s = 0.
      DO j = i + 1 , nmin
        s = s + t(i, j) * xx(j)
      ENDDO
      xx(i) = (tt(i) - s) / t(i, i)
      if (abs(xx(i)) < epsilon) xx(i)=0.
    ENDDO
    i = 0
    DO k = 1 , kinvariant
      i = i + 1
      min(kinvar(k)) = xx(i)
      IF (min(kinvar(k)) <= 0.)  ninv = 0
    ENDDO
    deallocate (t, tt, t0, tt0, xx, stat=alloc_res)
END SUBROUTINE reseqinv



double precision FUNCTION temperature (at, bt, ct, dt, et, temp)
	implicit none
	real (kind=8) :: at,bt,ct,dt,et,temp
    temperature = at + bt * temp + ct * temp ** 2 + dt * temp ** 3 + et * temp ** 4
END FUNCTION temperature



SUBROUTINE color(hconsole,ncolor)
use dfwin !WIN32 only
use data_common
implicit none
integer res,ncolor,BackColor,hconsole
integer(kind=2) entier2
BackColor=0
if(OS_S.eq.'winNT') then
	entier2=IOR(ncolor, ishft(BackColor, 4))
	res=SetConsoleTextAttribute (hconsole, entier2) !WIN32 only
endif
END SUBROUTINE color

!SUBROUTINE inkey(test_stop)
!	character (len=*)	:: test_stop
!	open(100,file='stop.dat',err=200)
!	read(100,*) test_stop
!	close (100)
!200 continue
!END SUBROUTINE inkey

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

subroutine stop_simulation
	character (len=30) :: d_and_t_S
	call fdate(d_and_t_S)
	write(*,*) d_and_t_S
	write(*,*) "Press <enter> to exit "
	read(*,'(a)') d_and_t_S
	STOP
end

