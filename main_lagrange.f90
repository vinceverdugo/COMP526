PROGRAM lagrange_program
    IMPLICIT NONE

    INTERFACE 
        SUBROUTINE lagrange_interpolation(x_data, y_data, x, y)
            REAL, dimension(:), intent(in) :: x_data
            REAL, dimension(:), intent(in) :: y_data
            REAL, dimension(:), intent(in) :: x
            REAL, dimension(:), intent(out) :: y
        END SUBROUTINE lagrange_interpolation
    END INTERFACE
    
    !!PART 1
    !Define parameters and arrays
    REAL :: x_data(5) !array to store abscissae values from input_data.txt
    REAL :: y_data(5) !array to store ordinate values from input_data.txt
    REAL :: x(50) !array to store domain points from domain_points.txt
    REAL :: y(50) !array that holds return values 
    REAL :: j = 0.0 !Iteration variable
    INTEGER :: i = 1

    !!PART 2
    !Open input_data.txt
    OPEN(1, file='input_data.txt')
    !do loop to read values from each line of input_data.txt' into x_data and y_data
    do i = 1, 5
        READ(1, *) x_data(i), y_data(i)
    end do
    !Close input_data.txt 
    CLOSE(1)

    !!PART 3
    WRITE(*,*) x_data
    WRITE(*,*) y_data

    !!PART 4
    !Read the domain points from domain_points.txt and store them into x array
    OPEN(2, file='domain_points.txt')
    READ(2,*) x

    !!PART 5
    !Call lagrange_interpolation subroutine 
    CALL lagrange_interpolation(x_data, y_data, x, y)

    !!PART 6
    !Create new file output.txt, and write the output from y
    OPEN(3, file='output.txt', status='new')
    !Use loop to write each element on new line - will be useful for EC
    do i = 1, SIZE(y)   
        WRITE(3, *) y(i)
    end do
    CLOSE(3)
    
    100 FORMAT(F3.1, F6.3) !Format specifier for two floating point values: one 3 columns wide, with 1 decimal place, and the other one 6 columns wide (counting from the end of the previous one), with 3 decimal places
    101 FORMAT(F19.17) !Floating point format specifier 19 columns wide, with 17 decimal places
    102 FORMAT(F0.17) !the 0 here means that processor selects the smallest positive field width necessary

END PROGRAM lagrange_program

SUBROUTINE lagrange_interpolation(x_data, y_data, x, p)
    IMPLICIT NONE
    REAL, dimension(:), intent(in) :: x_data
    REAL, dimension(:), intent(in) :: y_data
    REAL, dimension(:), intent(in) :: x
    REAL, dimension(:), intent(out) :: p
    INTEGER :: i, j !Iteration variables
    !Define w
    REAL :: w = 0.0 !w represents interpolated value of x, ie. p(x) = w

    !Iterate over each value of x
    do i = 1, SIZE(x)
        !Calculate p(x(i)) = w
        do j = 1, SIZE(y_data)
            w = w + y_data(j) * ((PRODUCT(x(i) - x_data(1:j-1)))*(PRODUCT(x(i) - x_data(j+1:))) / (PRODUCT(x_data(j) - x_data(1:j-1))*PRODUCT(x_data(j) - x_data(j+1:))))
        end do
        !Update p
        p(i) = w
        !Reinitialize w
        w = 0.0
    end do
END SUBROUTINE lagrange_interpolation