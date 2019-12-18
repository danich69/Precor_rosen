module diffur

    implicit none
    
    integer, parameter :: dims = 4
    integer, parameter :: ord = 4
    integer :: i
    
    real(8), parameter :: interval = 2
    real(8), dimension(dims), parameter :: x_0 = (/1, -1, 1, 1/)
    real(8), parameter :: h = 1e-3

    contains
    
        function f(x)
        
            real(8), dimension(:) :: x
            real(8), dimension(size(x)) :: f
            
             f(1) = x(3)
             f(2) = x(4)
             f(3) = cos(x(3)) - x(4)**2
             f(4) = sin(x(4)) - x(1)
        
        end function f

end module diffur
