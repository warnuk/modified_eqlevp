program test
    implicit none
    character(len=256) :: label, system, units
    character(len=256) :: add_minerals, rem_minerals
    integer :: print_step, output_step
    real(kind=8) :: temp, dens, ph, pco2, pkmol
    real(kind=8) :: pkeq, increment, dilute, max_tds
    real(kind=8) :: na, k, li, ca, mg, cl, so4, no3, b, si, alk
    integer :: verbose, output
    
    open (20, file="example.txt", action='READ')
        read(20,*) label; read(20,*) temp; read(20,*) dens
        read(20,*) ph; read(20,*) na; read(20,*) k
        read(20,*) li; read(20,*) ca; read(20,*) mg
        read(20,*) cl; read(20,*) so4; read(20,*) no3
        read(20,*) b; read(20,*) si; read(20,*) alk
        read(20,*) pco2; read(20,*) system; read(20,*) units
        read(20,*) dilute; read(20,*) add_minerals; read(20,*) rem_minerals
        read(20,*) max_tds; read(20,*) pkmol; read(20,*) pkeq
        read(20,*) print_step; read(20,*) output_step
        read(20,*) verbose; read(20,*) output
    close(20)

    write(*,*) label; write(*,*) temp; write(*,*) dens
    write(*,*) ph; write(*,*) na; write(*,*) k
    write(*,*) li; write(*,*) ca; write(*,*) mg
    write(*,*) cl; write(*,*) so4; write(*,*) no3
    write(*,*) b; write(*,*) si; write(*,*) alk
    write(*,*) pco2; write(*,*) system; write(*,*) units
    write(*,*) dilute; write(*,*) add_minerals; write(*,*) rem_minerals
    write(*,*) max_tds; write(*,*) pkmol; write(*,*) pkeq
    write(*,*) print_step; write(*,*) output_step
    write(*,*) verbose; write(*,*) output

end program test