program in_out
	use diffur
    use odu
    use autonomus
    implicit none
    	
!	-----------------------------------------------------
!	|	Program integrates differential equation   		|
!	|	Available methods are:              			|
!	|			Runge-Kutta - rk key					|
!	|			Adams interpolation - ai key			|
!	|			Adams extrapolation - ae key			|
!	-----------------------------------------------------

    real(8), dimension(1:int(interval/h + 1), size(x_0)) :: x
    real(8) :: t
    character*6 :: method
    character*10 :: file_name
    
    call getarg(1, method)
    
    select case(trim(method))
        case('rosen')
            call Rosenbrock(x)
            file_name = 'rosen.dat'
        case('precor')
            call Predictor_corrector(x)
            file_name = 'precor.dat'
        case default
            write(*,*) 'Inkorrect method! Please enter rosen (Rosenbrock), precor (Predictor_corrector)'
            read(*,*) method
    end select
    
    t = 0
    
    open(1, file = trim(file_name))
    do i = 1, int(interval/h + 1)
        
        write(1,*) t, x(i, :)
        t = t + h
    
    enddo    
    close(1)

end program in_out
