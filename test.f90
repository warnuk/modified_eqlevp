module data_common
	
end module data_common

function f_ucase(strIn) result(strOut)
    implicit none
    character(len=*), intent(in) :: strIn
    character(len=len(strIn)) :: strOut
    integer :: i,j

    do i = 1, len(strIn)
        j = iachar(strIn(i:i))
        if (j>= iachar("a") .and. j<=iachar("z") ) then
            strOut(i:i) = achar(iachar(strIn(i:i))-32)
        else
            strOut(i:i) = strIn(i:i)
        end if
    end do

end function f_ucase

function f_lcase(strIn) result(strOut)
    implicit none
    character(len=*), intent(in) :: strIn
    character(len=len(strIn)) :: strOut
    integer :: i,j

    do i = 1, len(strIn)
        j = iachar(strIn(i:i))
        if (j>= iachar("A") .and. j<=iachar("Z") ) then
            strOut(i:i) = achar(iachar(strIn(i:i))+32)
        else
            strOut(i:i) = strIn(i:i)
        end if
    end do
end function f_lcase

program main
    use data_common
    implicit none

    integer, parameter	:: maxlen=256, len_display=50, max_constit=30, max_nr=10, max_conv=100
    real(kind=8)	:: pk
	real(kind=8)	:: epsilon=1.D-8

    integer	        :: npasimp, npasfich, npasecran, alloc_res
    real(kind=8)	:: temp, tinit, dens, ph
    real(kind=8)	:: sodium, potassium, lithium, calcium, magnesium
    real(kind=8)	:: chlorine, sulfate, nitrate, boron, silicon, alk
    real(kind=8)	:: pco2, dil, stdmax, pkmol, pkeq
    character (len=maxlen)	:: syst_S, out_units, p_S, f_S, x_S
    character (len=maxlen)	:: add_min, rem_min

    real(kind=8), allocatable, dimension (:)	:: tot

    allocate (tot(0:12))

    f_S = "./" // "example.txt"

    open (10,file=f_S,action='READ')
        read(10,*) x_S
        x_S=f_lcase(x_S)

        read(10,*) temp
        read(10,*) dens
        read(10,*) ph

        read(10,*) sodium
        read(10,*) potassium
        read(10,*) lithium
        read(10,*) calcium
        read(10,*) magnesium
        read(10,*) chlorine
        read(10,*) sulfate
        read(10,*) nitrate
        read(10,*) boron
        read(10,*) silicon

        read(10,*) alk
        read(10,*) pco2

        read(10,*) syst_S
        read(10,*) out_units
        read(10,*) dil
        read(10,*) add_min
        read(10,*) rem_min
        read(10,*) stdmax
        read(10,*) pkmol
        read(10,*) pkeq

        read(10,*) p_S
        if (p_S == "") then
            npasecran = 1
        else 
            read(p_S,*) npasecran
        end if

        read(10,*) p_S
        if (p_S == "") then
            npasimp = 1
        else 
            read(p_S,*) npasimp
        end if
        npasfich = npasimp
        
        tot(1) = sodium / 1000.
        tot(2) = potassium / 1000.
        tot(3) = lithium / 1000.
        tot(4) = calcium / 1000.
        tot(5) = magnesium / 1000.
        tot(6) = chlorine / 1000.
        tot(7) = sulfate / 1000.
        tot(8) = nitrate / 1000.
        tot(9) = boron / 1000.
        tot(10) = silicon / 1000.
        tot(12) = alk / 1000.
    close(10)

    tinit = temp

    print *, tot(1)

end program main