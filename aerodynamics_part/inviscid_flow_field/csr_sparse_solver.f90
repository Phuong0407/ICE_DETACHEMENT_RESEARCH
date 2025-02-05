module csr_sparse_solver
    use, intrinsic :: iso_c_binding
    implicit none
    private
    real(8), save :: zero_compared_tol, matrix_norm
    integer, save :: size_matrix, non_zero_element
    real(8), allocatable, save :: val(:)
    integer, allocatable, save :: col_idx(:)
    integer, allocatable, save :: row_ptr(:)
    real(8), save :: conv_tol = 1.0E-6
    integer, save :: max_iter = 100000

    public :: init_sparse_solver, solve_sparse_system, mult_sparse_matrix

contains

    subroutine init_sparse_solver(sparse_matrix, n) bind(C, name="init_sparse_solver")
        implicit none
        integer(c_int), intent(in) :: n
        real(8), intent(in), dimension(n, n) :: sparse_matrix
        integer :: i, j, index

        non_zero_element = 0
        size_matrix = n

        matrix_norm = 0.0
        do i = 1, size_matrix
            matrix_norm = max(matrix_norm, sum(abs(sparse_matrix(i, :))))
        end do
        zero_compared_tol = epsilon(1.0) * matrix_norm

        do i = 1, size_matrix
            do j = 1, size_matrix
                if (abs(sparse_matrix(i, j)) > zero_compared_tol) then
                    non_zero_element = non_zero_element + 1
                end if
            end do
        end do

        allocate(val(non_zero_element))
        allocate(col_idx(non_zero_element))
        allocate(row_ptr(size_matrix + 1))

        row_ptr = 0
        index = 1
        row_ptr(1) = 1
        do i = 1, size_matrix
            do j = 1, size_matrix
                if (abs(sparse_matrix(i, j)) > zero_compared_tol) then
                    val(index) = sparse_matrix(i, j)
                    col_idx(index) = j
                    index = index + 1
                end if
            end do
            row_ptr(i + 1) = index
        end do
    end subroutine init_sparse_solver

    subroutine mult_sparse_matrix(x, y, n) bind(C, name="mult_sparse_matrix")
        implicit none
        integer(c_int), intent(in) :: n
        real(8), intent(in) :: x(n)
        real(8), intent(out) :: y(n)
        integer :: i, j

        y = 0.0d0
        do i = 1, n
            do j = row_ptr(i), row_ptr(i + 1) - 1
                y(i) = y(i) + val(j) * x(col_idx(j))
            end do
        end do
    end subroutine mult_sparse_matrix

    subroutine solve_sparse_system(rhs_vect, solution, n) bind(C, name="solve_sparse_system")
        implicit none
        integer(c_int), intent(in) :: n
        real(8), intent(in) :: rhs_vect(n)
        real(8), intent(out) :: solution(n)

        real(8), allocatable :: r(:), p(:), Ap(:)
        real(8) :: alpha, beta, rnorm, bnorm, tol
        integer :: iter

        allocate(r(n), p(n), Ap(n))
        tol = conv_tol

        solution = 0.0d0
        bnorm = sqrt(dot_product(rhs_vect, rhs_vect))

        call mult_sparse_matrix(solution, Ap, n)
        r = rhs_vect - Ap
        p = r
        rnorm = dot_product(r, r)

        iter = 0
        do while (sqrt(rnorm) / bnorm > tol .and. iter < max_iter)
            iter = iter + 1
            call mult_sparse_matrix(p, Ap, n)
            alpha = rnorm / dot_product(p, Ap)
            if (alpha == 0.0d0) exit

            solution = solution + alpha * p
            r = r - alpha * Ap

            if (sqrt(dot_product(r, r)) / bnorm < tol) exit
            beta = dot_product(r, r) / rnorm
            if (beta == 0.0d0) exit

            p = r + beta * p
            rnorm = dot_product(r, r)
        end do

        if (iter == max_iter) then
            print *, "FAILED TO CONVERGE IN ", max_iter, " ITERATIONS."
        end if
        deallocate(r, p, Ap)
    end subroutine solve_sparse_system

end module csr_sparse_solver