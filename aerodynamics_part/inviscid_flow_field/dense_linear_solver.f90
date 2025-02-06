module dense_matrix_solver
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr, c_loc
    implicit none
    private
    public :: solve_dense_system

contains

    subroutine decompose_LU(A, idx, parity, is_singular)
        implicit none
        real(8), parameter :: tiny = 1.5D-16
        real(8), intent(inout), dimension(:,:) :: A
        integer, intent(out) :: parity, is_singular
        integer, intent(out), dimension(:) :: idx
        integer :: mat_size, i, j, k, i_max
        real(8) :: sum, row_max_val, pivot, pivot_comparator, row_swap, pivot_divisor
        real(8), allocatable :: implicit_scaling_factor(:)

        mat_size = size(A,1)
        allocate(implicit_scaling_factor(mat_size))

        parity = 1
        is_singular = 0

        do i = 1, mat_size
            row_max_val = 0.d0
            do j = 1, mat_size
                if (abs(A(i,j)) > row_max_val) row_max_val = abs(A(i,j))
            end do
            if (row_max_val < tiny) then
                is_singular = 1
                deallocate(implicit_scaling_factor)
                return
            end if
            implicit_scaling_factor(i) = 1.d0 / row_max_val
        end do

        do j = 1, mat_size
            do i = 1, j - 1
                sum = A(i,j)
                do k = 1, i - 1
                    sum = sum - A(i,k) * A(k,j)
                end do
                A(i,j) = sum
            end do

            pivot = 0.d0
            do i = j, mat_size
                sum = A(i,j)
                do k = 1, j - 1
                    sum = sum - A(i,k) * A(k,j)
                end do
                A(i,j) = sum
                pivot_comparator = implicit_scaling_factor(i) * abs(sum)
                if (pivot_comparator >= pivot) then
                    i_max = i
                    pivot = pivot_comparator
                end if
            end do

            if (j /= i_max) then
                do k = 1, mat_size
                    row_swap = A(i_max, k)
                    A(i_max, k) = A(j, k)
                    A(j, k) = row_swap
                end do
                parity = -parity
                implicit_scaling_factor(i_max) = implicit_scaling_factor(j)
            end if

            idx(j) = i_max
            if (abs(A(j, j)) < tiny) A(j, j) = tiny

            if (j /= mat_size) then
                pivot_divisor = 1.d0 / A(j, j)
                do i = j + 1, mat_size
                    A(i,j) = A(i,j) * pivot_divisor
                end do
            end if
        end do

        deallocate(implicit_scaling_factor)
    end subroutine decompose_LU

    subroutine substitute_backward_LU(A, idx, B)
        implicit none
        real(8), intent(in), dimension(:,:) :: A
        integer, intent(in), dimension(:) :: idx
        real(8), intent(inout), dimension(:) :: B
        integer :: mat_size, first_non_null_idx, row_perm_idx, i, j
        real(8) :: sum

        mat_size = size(A,1)
        first_non_null_idx = 0

        do i = 1, mat_size
            row_perm_idx = idx(i)
            sum = B(row_perm_idx)
            B(row_perm_idx) = B(i)
            if (first_non_null_idx /= 0) then
                do j = first_non_null_idx, i - 1
                    sum = sum - A(i,j) * B(j)
                end do
            else if (sum /= 0.d0) then
                first_non_null_idx = i
            end if
            B(i) = sum
        end do

        do i = mat_size, 1, -1
            sum = B(i)
            if (i < mat_size) then
                do j = i + 1, mat_size
                    sum = sum - A(i,j) * B(j)
                end do
            end if
            B(i) = sum / A(i,i)
        end do
    end subroutine substitute_backward_LU

    subroutine solve_dense_system(A, B, solution_ptr, n) bind(C, name="solve_dense_system")
        implicit none
        integer(c_int), intent(in) :: n
        real(c_double), intent(inout), dimension(n, n) :: A
        real(c_double), intent(in), dimension(n) :: B
        type(c_ptr), intent(out) :: solution_ptr

        real(c_double), allocatable, save :: solution(:)
        integer(c_int), allocatable :: idx(:)
        integer(c_int) :: is_singular, parity

        if (allocated(solution)) then
            if (size(solution) /= n) then
                deallocate(solution)
            end if
        end if
        if (.not. allocated(solution)) then
            allocate(solution(n))
        end if

        allocate(idx(n))
        solution = B

        call decompose_LU(A, idx, parity, is_singular)
        if (is_singular == 0) then
            call substitute_backward_LU(A, idx, solution)
        end if

        solution_ptr = c_loc(solution)
        deallocate(idx)
    end subroutine solve_dense_system

end module dense_matrix_solver