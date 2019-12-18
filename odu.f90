module odu
	use diffur
    use non_linear_system
    implicit none
    	
    contains
    
        subroutine Runge_Kutta(x)
        
            real(8), dimension(:, :) :: x
            real(8), dimension(dims) :: k1, k2, k3, k4
            
            x(1, :) = x_0
            
            do i = 1, size(x)/dims - 1
                
                k1 = h * f(x(i, :))
                k2 = h * f(x(i, :) + k1/2.0)
                k3 = h * f(x(i, :) + k2/2.0)
                k4 = h * f(x(i, :) + k3)

                x(i+1, :) = x(i, :) + (k1 + 2 * k2 + 2 * k3 + k4)/6.0
            
            enddo
            
        end subroutine Runge_Kutta
        
        function Adams_extrapolation(x, A) result(solution)
        
            real(8), dimension(:, :) :: x
            real(8), dimension(size(x,2)) :: solution
            real(8), dimension(0 : ord-1) :: A
            real(8), dimension(dims) :: s
            integer :: i, j
            
            i = ord
            
            s = 0
            
            do j = 0, ord-1
                
                s = A(j)*f(x(i - j, :)) + s
                
            enddo
                
            solution = h * s + x(i, :)
            
        
        end function Adams_extrapolation
        
        function A_coefficents(n)
        
            real(8), dimension(0 : n-1) :: A_coefficents
            real(8), dimension(0 : n-1) :: powers, new
            real(8) :: a
            integer :: n, j, i
            
            do j = 0, n-1
                
                powers = 0.0
                powers(0) = 1.0
                
                do i = 0, n-1
                
                    if (i /= j) then
                    
                        new = powers * i
                        powers = eoshift(powers, -1)
                        powers = powers + new
                        
                    endif
                    
                enddo
                
                do i = 0, n-1
                    powers(i) = powers(i)/real(i+1)
                enddo
                
                A_coefficents(j) = (-1.0_8)**j/(fact(j)*fact(n - 1 - j) )* sum(powers)
                
            enddo
        
        end function A_coefficents
        
        function Adams_interpolation(x, prediction, B) result(solution)
        
            real(8), dimension(:, :) :: x
            real(8), dimension(size(x,2)) :: solution, prediction
            real(8), dimension(-1 : ord-2) :: B
            real(8), dimension(dims) :: s
            integer :: i, j
            
            i = ord
            s = 0
            
            do j = 0, ord-2
                
                s = B(j)*f(x(i - j, :)) + s
                
            enddo
                
            call nonlinear_system_solver(prediction, funct, solution)
            
            contains
                    
                function funct(y)
        
                    real(8), dimension(:) :: y
                    real(8), dimension(size(y)) :: funct
            
                    funct = x(i, :) + h * s + h * B(-1) * f(y) - y
            
                end function funct
        
        end function Adams_interpolation
        
        function B_coefficents(n)
        
            real(8), dimension(-1 : n-2) :: B_coefficents
            real(8), dimension(-1 : n-2) :: powers, shifter
            real(8) :: a
            integer :: n, j, i
            
            do j = -1, n-2
                
                powers = 0
                powers(-1) = 1
                
                do i = -1, n-2
                
                    if (i /= j) then
                    
                        shifter = powers * i
                        powers = eoshift(powers, -1)
                        powers = powers + shifter
                        
                    endif
                    
                enddo
                
                do i = -1, n-2
                    powers(i) = powers(i)/(i+2)
                enddo
                
                B_coefficents(j) = (-1.0)**(j+1)/(fact(j+1)*fact(n - 2 - j) )* sum(powers)
                
            enddo
        
        end function B_coefficents
        
        function fact(n)
            
            integer, intent(in) :: n
            integer :: i, fact
            
            fact = 1.0
            do i = 2, n
                fact = fact * i
            enddo

        end function fact

end module odu
