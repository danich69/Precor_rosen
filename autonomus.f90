module autonomus
    use odu
    use diffur
    use linear_system
    use non_linear_system
    implicit none
    
    contains
    
        subroutine Rosenbrock(x)
        
            real(8), dimension(:, :) :: x
            real(8), dimension(dims, dims) :: Jacobi_matrix, E, S
            real(8), dimension(dims) :: solution_vector
            
            real(8) :: alpha, beta, gamma_const
            integer :: i
            
            x(1, :) = x_0
            E = 0
            forall(i = 1:dims) E(i,i) = 1.0_16
            alpha = 1.077_16
            beta = -0.372_16
            gamma_const = -0.577_16
            
            do i = 1, size(x)/dims - 1
                
                Jacobi_matrix = matrix_derivative(f, x(i, :))
                S = (E - alpha * h * Jacobi_matrix - beta * h**2 * matmul(Jacobi_matrix, Jacobi_matrix))
                solution_vector = f(x(i,:) + gamma_const * h * f(x(i,:))) * h + matmul(S, x(i, :))
                call Lin_sys_solver(S, solution_vector, x(i+1, :), dims, 'ModGauss')
            
            enddo

        end subroutine Rosenbrock
        
        subroutine Predictor_corrector(x)
        
            real(8), dimension(:, :) :: x
            real(8), dimension(size(x, 2)) :: prediction
            real(8), dimension(ord) :: A, B
            real(8) :: t1, t2, t3, t4, textr, tint
            
            A = A_coefficents(ord)
            B = B_coefficents(ord)
            
            x(1, :) = x_0
            call Runge_Kutta(x(1 : ord, :))
            textr = 0
            tint = 0
            do i = 0, size(x)/dims - 1 - ord
            
                prediction = Adams_extrapolation(x(i+1:ord + i, :), A)
                x(ord + i+1, :) = Adams_interpolation(x(i+1:ord + i, :), prediction, B)
            
            enddo
        
        end subroutine Predictor_corrector

end module autonomus
