module csr_sparse_solver
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr, c_loc
    implicit none
    private
    real(8), save :: zero_compared_tol, matrix_norm
    integer, save :: size_matrix, non_zero_element
    real(8), allocatable, save :: val(:)
    integer, allocatable, save :: col_idx(:)
    integer, allocatable, save :: row_ptr(:)
    real(8), save :: conv_tol = 1.0E-6
    integer, save :: max_iter = 100000

    public :: solve_sparse_system

contains

    subroutine init_sparse_solver(sparse_matrix, n)
        implicit none
        integer(c_int), intent(in) :: n
        real(8), intent(in), dimension(n * n) :: sparse_matrix
        integer :: i, j, index

        non_zero_element = 0
        size_matrix = n

        matrix_norm = 0.0
        do i = 1, size_matrix
            matrix_norm = max(matrix_norm, sum(abs(sparse_matrix((i-1) * n + 1 : i * n))))
        end do
        zero_compared_tol = epsilon(1.0) * matrix_norm

        do i = 1, size_matrix
            do j = 1, size_matrix
                if (abs(sparse_matrix((i-1) * n + j)) > zero_compared_tol) then
                    non_zero_element = non_zero_element + 1
                end if
            end do
        end do

        if (allocated(val)) deallocate(val)
        if (allocated(col_idx)) deallocate(col_idx)
        if (allocated(row_ptr)) deallocate(row_ptr)

        allocate(val(non_zero_element))
        allocate(col_idx(non_zero_element))
        allocate(row_ptr(size_matrix + 1))

        row_ptr = 0
        index = 1
        row_ptr(1) = 1
        do i = 1, size_matrix
            do j = 1, size_matrix
                if (abs(sparse_matrix((i - 1) * n + j)) > zero_compared_tol) then
                    val(index) = sparse_matrix((i-1) * n + j)
                    col_idx(index) = j
                    index = index + 1
                end if
            end do
            row_ptr(i + 1) = index
        end do
    end subroutine init_sparse_solver

    subroutine mult_sparse_matrix(x, y, n)
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

    subroutine solve_sparse_system_from_rhs(rhs_vect, n, solution)
        implicit none
        integer, intent(in) :: n
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
            print *, "Failed to converge in ", max_iter, " iterations in FORTRAN MODULE CSR_SPARSE_SOLVER!"
        end if
        deallocate(r, p, Ap)
    end subroutine solve_sparse_system_from_rhs

    subroutine solve_sparse_system(sparse_matrix, rhs, solution_ptr, n) bind(C, name="solve_sparse_system")
        use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr, c_loc
        implicit none
        integer(c_int), intent(in) :: n
        real(c_double), intent(in), dimension(n * n) :: sparse_matrix
        real(c_double), intent(in), dimension(n) :: rhs
        type(c_ptr), intent(out) :: solution_ptr

        real(c_double), allocatable, save :: solution(:)

        if (allocated(solution)) then
            if (size(solution) /= n) then
                deallocate(solution)
            end if
        end if

        if (.not. allocated(solution)) then
            allocate(solution(n))
        end if

        call init_sparse_solver(sparse_matrix, n)
        call solve_sparse_system_from_rhs(rhs, n, solution)

        solution_ptr = c_loc(solution)

    end subroutine solve_sparse_system

end module csr_sparse_solver