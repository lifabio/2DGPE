!  $2DGPE_v1.f90 
!
!  FUNCTIONS:
!  $2DGPE_v1 - 
!

!****************************************************************************
!
!  PROGRAM: $2DGPE_v1
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************
    Module dati
        implicit none
        save
	    !Double/Single precision definition
        integer, parameter :: single = selected_real_kind(6,37)
        integer, parameter :: double = selected_real_kind(15,307)

        !Global const.
        real(kind=double), parameter    :: pi=3.1415926535897932384626433832795d0
        complex(kind=double), parameter :: iu=(0.d0,1.d0)   ! imaginary Unit

	    !Global Variables

    endmodule dati
    
    program $2DGPE_v1
    use dati
    implicit none

    ! Variables

    ! Body of $2DGPE_v1
    print *, 'Hello World'

    end program $2DGPE_v1

