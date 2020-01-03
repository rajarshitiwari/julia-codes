#!/usr/bin/julia

# using Plots
using LinearAlgebra

"""
function that return the tight binding hamiltonian
"""
function generate_tb_hamiltonian_orb_sqr!(matrix::Array{Float64, 2}, hoppings::Vector{Float64}, twolsize::Int64)
    #real(8), dimension(:,:), intent(inout) :: matrix #matrix to be diagonalized
    #real(8), dimension(:), intent(in) :: hoppings    #hopping amplitudes
    #integer, intent(in) :: twolsize
    # integer :: m1,m2,co(2),xi(2),ii,jj
    # integer :: hsize

    matrix .= zero(Float64)
    nsize = twolsize * twolsize
    hsize = 2 * nsize # same as hsize = *(2, twolsize, twolsize)
    co = Vector{Int64}(undef, 2)
    indices = reshape(Vector(1:nsize), (twolsize,twolsize))
    
    #for i1 = 1:twolsize, i2 = 1:twolsize
    #    indices[i1, i2] = twolsize * (i1 - 1) + i2
    #end
    oneton = Vector(1:twolsize)
    minusplus = hcat(circshift(oneton, -1), circshift(oneton, 1))
    
    #----- < define all tight binding entries > -----#
    for i1 = 1:twolsize, i2 = 1:twolsize
        xi = [i1,i2]
        ii = indices[i1, i2]
        # first four neighbours
        for m1 = 1:2, m2 = 1:2
            co .= xi
            co[m1] = minusplus[xi[m1], m2]
            jj = twolsize * (co[1] - 1) + co[2]
            matrix[2*ii-1, 2*jj-1] = -hoppings[1]
            matrix[2*ii  , 2*jj  ] = -hoppings[1]
            #
            matrix[2*ii-1, 2*jj  ] = -hoppings[2]
            matrix[2*ii  , 2*jj-1] = -hoppings[2]
        end
    end
    # return matrix
end # subroutine generate_tb_hamiltonian_orb_sqr


function eps0_k(k)
    eps0_k = -2.0 * sum(cos, k)
    return eps0_k
end # function eps0_k

# PROGRAM TB_SPECTRUM_ORB_SQR
# REAL(8), PARAMETER :: PI = 4.0D0 * ATAN(1.0D0)
# INTEGER :: LSIZE,TWOLSIZE,HSIZE,NUM_THREAD
# REAL(8), ALLOCATABLE, DIMENSION(:) :: EIGENVALUES
# REAL(8), ALLOCATABLE, DIMENSION(:,:) :: TB_MATRIX
# REAL(8), DIMENSION(2) :: HOPPING
# REAL(8) :: KEY(2),KPLUSQ(2),G1(2),G2(2),COEFF,E0,E1,E2
# REAL(8) :: ALPHA1,ALPHA2,BETA,BETASQ,T1,T2
# INTEGER :: I1,I2,I3,I,INFO

hopping = [0.0, 1.0]       ## DIRECT, INTER-ORBITAL HOPPING
lsize = 24
twolsize = 2 * lsize
hsize = 2 * twolsize^2

eigenvalues = zeros((hsize,))

function generate_tb_spectrum(hopping::Vector{Float64}, eigenvalues::Vector{Float64}, twolsize::Int64)
    hsize = 2 * twolsize^2
    tb_matrix = zeros((hsize,hsize))

    # println(tb_matrix)
    # Plots.PyPlotBackend()
    # plotly(tb_matrix)
    println("hsize= $hsize")
    time1 = time()
    generate_tb_hamiltonian_orb_sqr!(tb_matrix, hopping, twolsize)
    eig = eigvals!(tb_matrix)
    println(size(eig))

    time2 = time()
    println("time = ", time2 - time1)
    # eigenvalues = 0.0
    # generate_tb_hamiltonian_orb_sqr(tb_matrix,hopping,twolsize)
    
    # # call syevd(tb_matrix,eigenvalues,'n','u' ,info)
    # if INFO /= 0
    #     PRINT*,'PROBLEM IN DIAGONALIZATION, INFO = ', INFO
    # end
    
    # # CALL CPU_TIME(T2)
    
    # print'(f15.10)',t2-t1
    #   open(11,file='eigenvalues_r1.dat')
    #   do i = 1,hsize
    #      write(11,'(f15.10)')eigenvalues(i)
    #   end do
    #   close(11)
    
    #   call cpu_time(t1)
    #   eigenvalues = 0.0d0
    #   g1 = 2.0d0 * [pi,0.0d0]
    #   g2 = 2.0d0 * [0.0d0,pi]
    #   i = 1
    #   do i1 = 0, twolsize - 1
    #      do i2 = 0, twolsize - 1
    #         key = i1 * g1/dble(twolsize) + i2 * g2/dble(twolsize)
    #         e0 = eps0_k(key)
    #         eigenvalues(i:i+1) = e0 * [hopping(1) + hopping(2), hopping(1) - hopping(2)]
    #         i = i + 2
    #      end do
    #   end do
    #   call cpu_time(t2)
    #   print'(f15.10)',t2-t1
    #   call dlasrt('i',hsize,eigenvalues,info)
    #   open(11,file='eigenvalues_k1.dat')
    #   do i = 1,hsize
    #      write(11,'(f15.10)')eigenvalues(i)
    #   end do
    #   close(11)
    
end # END PROGRAM TB_SPECTRUM_ORB_SQR

generate_tb_spectrum(hopping, eigenvalues, twolsize)
