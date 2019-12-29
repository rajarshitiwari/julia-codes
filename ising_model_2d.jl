#!/usr/bin/env julia

using LinearAlgebra
const PI = 4.0*atan(1.0)
# Int :: indices(:,:)
# Int :: minusplus(:,:)


function metro_polis(de, temp)::Bool
    flag = false
    if de < 0.0
        flag = true
    else
        metro_police = rand()
        compare = exp(-de/(temp + 1.0e-10))
        if metro_police < compare
            flag = true
        else
            flag = false
        end
    end
    return flag
end

function get_energy_diff(spins_lattice::Array{Int}, tmp_spin::Int,
                         jey1::Float64, jey2::Float64, magnetic_field::Float64,
                         minusplus::Array{Int}, xi::Array{Int})
    energy_difference = 0.0
    # xi = [i1,i2]
    i1, i2 = xi
    # nearest neighbor part
    for m1 = 1:2, m2 = 1:2
        co = xi
        co[m1] = minusplus[xi[m1], m2]
        energy_difference += 2.0 * jey1 * tmp_spin * spins_lattice[co[1], co[2]]
    end
    # next nearest neighbour
    for m1 = 1:2, m2 = 1:2
        co = [minusplus[i1, m1], minusplus[i2, m2]]
        energy_difference += 2.0 * jey2 * tmp_spin * spins_lattice[co[1], co[2]]
    end
    
    previous_s = spins_lattice[i1, i2]
    current_s = tmp_spin
    energy_difference += -magnetic_field * current_s + magnetic_field * previous_s
    #
    return energy_difference
end

function total_energy_ising_2d(spins_lattice::Array{Int64}, jey1::Float64, jey2::Float64, magnetic_field::Float64)::Float64
    # INTEGER, DIMENSION(:,:), INTENT(IN) :: SPINS_LATTICE
    # REAL(8), INTENT(IN) :: JEY1,JEY2
    # REAL(8), INTENT(IN) :: MAGNETIC_FIELD

    co = zeros(Int, 2)
    xi = zeros(Int, 2)
    
    lsize = size(spins_lattice)[1]
    total_energy_ising_2d = 0.0
    energy_nn = 0
    energy_nnn = 0
    for i1 = 1:lsize, i2 = 1:lsize

        xi .= [i1, i2]    #NEAREST NEIGHBOR PART
        for m1 = 1:2, m2 = 1:1
            co .= xi
            co[m1] = minusplus[xi[m1], m2]
            energy_nn = energy_nn + spins_lattice[i1, i2] * spins_lattice[co[1], co[2]]
        end

        #NEXT NEAREST NEIGHBOR PART
        for m1 = 1:2, m2 = 1:1
            co[1] = minusplus[i1, m1]
            co[2] = minusplus[i2, m2]
            energy_nnn = energy_nnn + spins_lattice[i1, i2] * spins_lattice[co[1], co[2]]
        end
    end

    total_energy_ising_2d = jey1 * energy_nn + jey2 * energy_nnn
    total_s = sum(spins_lattice)
    total_energy_ising_2d = total_energy_ising_2d - magnetic_field * total_s
    return total_energy_ising_2d
end # function total_energy_ising_2d

function calculate_average_sssq(spins::Array{Int})
    # INTEGER, DIMENSION(:,:), INTENT(IN) :: SPINS
    # REAL(8), INTENT(OUT) :: TOTAL_S,TOTAL_SSQ

    total_s::Float64 = sum(spins[:, :])
    total_ssq::Float64 = total_s * total_s
    return [total_s, total_ssq]
end # end subroutine calculate_average_sssq

function create_indices(lsize::Int)
    # INTEGER, INTENT(IN) :: LSIZE

    indices = zeros(Int, (lsize, lsize))
    for i1 = 1:lsize, i2 = 1:lsize
        indices[i1, i2] = lsize * (i1 - 1) + i2
    end
    return indices
end # END SUBROUTINE CREATE_INDICES


function generate_random_spin()
    # INTEGER, INTENT(OUT) :: SPINCONF
    return rand() >= 0.50 ? 1 : -1
end # SUBROUTINE GENERATE_RANDOM_SPIN

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ###
### TAKE THE SPIN CONFIGURATION AND CALCULATE SI.SJ      ###
function calculate_si_sj(spins::Array{Int})::Array{Float64}
    # INTEGER, DIMENSION(:,:), INTENT(IN) :: SPINS
    # INTEGER, DIMENSION(:,:), INTENT(OUT) :: SI_SJ

    lsize::Int = size(spins, 1)
    nsize::Int = lsize * lsize
    si_sj = zeros((nsize, nsize))

    #ALL THE SITES HAVE NON ZERO SPINS
    for i1 = 1:lsize, i2 = 1:lsize
        ii::Int = indices[i1, i2]
        for j1 = 1:lsize, j2 = 1:lsize
            jj::Int = indices[j1, j2]
            #if ii < jj
            #    continue
            #end
            si_sj[ii, jj] = spins[i1, i2] * spins[j1, j2]
            #si_sj[jj, ii] = si_sj[ii, jj]
        end
    end
    return si_sj
end # SUBROUTINE CALCULATE_SI_SJ

function selected_sq_from_si_sj(si_sj::Array{Float64}, que::Vector{Int}, lsize::Int)::Float64
    # INTEGER, DIMENSION(:,:), INTENT(IN) :: SI_SJ
    # INTEGER, DIMENSION(:), INTENT(IN) :: QUE
    # INTEGER, INTENT(IN) :: LSIZE
    # REAL(8), INTENT(OUT) :: STRUCTURE_FACTOR

    nsize::Int = lsize * lsize
    coeff::Float64 = 2.0 * PI/ lsize
    normalize = nsize^2

    tmpdp = 0.0
    for i = 1:nsize
        tmpdp += si_sj[i, i]
    end
    # tmpdp = sum(diag(si_sj, 0))
        
    structure_factor = tmpdp

    rij = zeros(Int, 2)
    q = que
    for i1 = 1:lsize, j1 = 1:lsize
        rij[1] = i1 - j1
        for i2 = 1:lsize, j2 = 1:lsize
            rij[2] = i2 - j2
            ii, jj = indices[i1, i2], indices[j1, j2]
            if jj >= ii
                continue
            end
            itmp = dot(q, rij)
            tmpdp = coeff * itmp
            tmpdp = si_sj[ii, jj] * 2.0 * cos(tmpdp)
            structure_factor += tmpdp
        end
    end
    structure_factor = structure_factor / normalize
    return structure_factor
end # SUBROUTINE SELECTED_SQ_FROM_SI_SJ

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ###
### GIVEN SI.SJ FOR THE SPIN CONFIGURATION, CALCULATE    ###
### THE STRUCTURE FACTOR S(Q) FROM THIS                  ###
function structure_factor_from_si_sj(si_sj::Array{Float64}, lsize::Int)::Array{Float64}
    # INTEGER, DIMENSION(:,:), INTENT(IN) :: SI_SJ
    # REAL(8), DIMENSION(:,:), INTENT(OUT) :: STRUCTURE_FACTOR

    structure_factor = zeros(lsize, lsize)
    que = zeros(Int, 2)
    
    # LOOP Q1
    for q1 = 1:lsize
        que[1] = q1 - 1
        # LOOP Q2
        for q2 = 1:lsize
            que[2] = q2 - 1
            #
            structure_factor[q1, q2] = selected_sq_from_si_sj(si_sj, que, lsize)
            #
        end # DO LOOP_Q2
    end # DO LOOP_Q1

    return structure_factor

end # SUBROUTINE STRUCTURE_FACTOR_FROM_SI_SJ

# PROGRAM COOL_ISING_SPINS_3D

#   USE MYTYPE

#   IMPLICIT NONE

const lsize = 20
const nsize = lsize * lsize
#   UNIT_ES   = 100
#   UNIT_MX   = 101
#   UNIT_SP   = 102
#   UNIT_SQ   = 103
#   UNIT_SISJ = 104
#   UNIT_EPS  = 105
#   UNIT_ACC  = 106
#   UNIT_TMP  = 111

#   DEFINE THE TIME OF SWEEP AND NUMBER OF STEPS TO EQUILIBRIATE.
global maxmcsweep = 1000 #0
global nequil = 500 #0
global ndel = 1
global number_of_temperature = 60

#   #======READ THE INPUT VARIABLES FIRSTLY FROM ARGUMENT
#   #IF THAT NOT GIVEN THEN FROM STANDARD INPUT=======#
global jey1 = -1.0
global jey2 =  0.0
#   PRINT*,"ENTER JEY2,MAX TEMPERATURE"
#   READ(*,*)JEY2,MAX_TEMPERATURE #JEY2 = 0.5D0

global magnetic_field = 0.0
global max_temperature = 3.0 * abs(jey1)
global delta_t = max_temperature / (number_of_temperature - 1)

si_sj = zeros(nsize, nsize)
global sf_ising = zeros(lsize, lsize) #   ALLOCATE(SF_ISING(0:LSIZE-1,0:LSIZE-1))

global spins_lattice = zeros(Int, (lsize, lsize))


#   ALLOCATE(MINUSPLUS(LSIZE,2),ONETON(LSIZE))

global indices = create_indices(lsize)
global minusplus = [circshift(Vector(1:lsize), 1) circshift(Vector(1:lsize), -1)]

#   OPEN(UNIT_TMP, FILE = 'ising_2d.inp')
#   WRITE(UNIT_TMP,'("JEY1            = ",F15.5)')JEY1
#   WRITE(UNIT_TMP,'("JEY2            = ",F15.5)')JEY2
#   WRITE(UNIT_TMP,'("LSIZE           = ",I15)')LSIZE
#   WRITE(UNIT_TMP,'("NSIZE           = ",I15)')NSIZE
#   WRITE(UNIT_TMP,'("SWEEP LENGTH    = ",I15)')MAXMCSWEEP
#   WRITE(UNIT_TMP,'("EQUIL LENGTH    = ",I15)')NEQUIL
#   WRITE(UNIT_TMP,'("NO. OF T        = ",I15)')NUMBER_OF_TEMPERATURE
#   WRITE(UNIT_TMP,'("MAXIMUM T       = ",F15.5)')MAX_TEMPERATURE
#   WRITE(UNIT_TMP,'("DELTA T         = ",F15.5)')DELTA_T
#   WRITE(UNIT_TMP,'("MAGNETIC FIELD  = ",3(F15.5,1X))')MAGNETIC_FIELD

#   CLOSE(UNIT_TMP)

#   WRITE(*,'("JEY1            = ",F15.5)')JEY1
#   WRITE(*,'("JEY2            = ",F15.5)')JEY2
#   WRITE(*,'("LSIZE           = ",I15)')LSIZE
#   WRITE(*,'("NSIZE           = ",I15)')NSIZE
#   WRITE(*,'("SWEEP LENGTH    = ",I15)')MAXMCSWEEP
#   WRITE(*,'("EQUIL LENGTH    = ",I15)')NEQUIL
#   WRITE(*,'("NO. OF T        = ",I15)')NUMBER_OF_TEMPERATURE
#   WRITE(*,'("MAXIMUM T       = ",F15.5)')MAX_TEMPERATURE
#   WRITE(*,'("DELTA T         = ",F15.5)')DELTA_T
#   WRITE(*,'("MAGNETIC FIELD  = ",3(F15.5,1X))')MAGNETIC_FIELD

#GENERATE THE SNAPSHOT OF SPINS (RANDOM CONFIGURATION)#

for i1 = 1:lsize, i2 = 1:lsize
    #SPINS_LATTICE(I1,I2) = 1
    spins_lattice[i1, i2] = generate_random_spin()
end
unit_tmp = open("starting_spins.dat", "w")
for i1 = 1:lsize
    for i2 = 1:lsize
        print(unit_tmp, string(spins_lattice[i1, i2], " "))
    end
    print(unit_tmp,"\n")
end
close(unit_tmp)

#   #=====================================#
#   #=======OPEN THE FILES TO WRITE=======#
#   OPEN(UNIT_ES, FILE = 'energy_susceptibility.dat')
#   OPEN(UNIT_MX, FILE = 'magnetization_chai.dat')
#   OPEN(UNIT_SP, FILE = 'spin_configurations.dat')
#   OPEN(UNIT_SQ, FILE = 'structure_factors.dat')
#   OPEN(UNIT_SISJ, FILE = 'avgsisj.dat')
#   OPEN(UNIT_ACC, FILE = 'acceptence.dat')

#   WRITE(UNIT_ES,'(3(A15,1X))')'TEMPERATURE','AVERAGE_ENERGY','SPECIFIC_HEAT'
#   WRITE(UNIT_MX,'(6(A15,1X))')'TEMPERATURE','SQ_FERRO','SQ_G_TYPE_AF',&
#        &'CHAI','S','SSQ'
#   #======================================#

total_energy = total_energy_ising_2d(spins_lattice, jey1, jey2, magnetic_field)
println("TOTAL ENERGY IN STARTING = ", total_energy / nsize)

# println(spins_lattice)
println("STARTING THE T LOOP")

# total_s = sum(spins_lattice)

previous_energy = total_energy

temperature = max_temperature
println(temperature)
count_temp = 0

while temperature::Float64 > 0.0         # temperature_loop
    #
    global mcsweep = 1
    count = 0
    oldcount = 0
    outcount = 0
    global avg_sisj = zeros(nsize, nsize)
    # global avg_sisj .= 0.0
    global avge = 0.0
    global avgesq = 0.0
    magnetization = 0.0
    sqaf = 0.0
    global average_s = 0.0
    global average_ssq = 0.0
    chai = 0.0
    global imc = 0
    global ifail = 0
    global ipass = 0
    # CALL CPU_TIME(TIME1) #TIME1 = X05BAF()
    # MC_LOOP
    @time while mcsweep <= maxmcsweep
        #
        # SWEEP_I1, SWEEP_I2
        for i1 = 1:lsize, i2 = 1:lsize
            tmp_spin = - spins_lattice[i1, i2] # FLIP THE TEMP SPIN
            ########################################
            # CALCULATE THE ENERGY DIFFERENCE NOW  #
            ########################################
            energy_difference = get_energy_diff(spins_lattice, tmp_spin, jey1, jey2, magnetic_field, minusplus, Vector{Int}([i1, i2]))
            #
            flag = metro_polis(energy_difference, temperature)
            if flag
                ipass += 1
                spins_lattice[i1, i2] = tmp_spin # update the temp spin
                final_energy = previous_energy + energy_difference
            else
                ifail += 1
                final_energy = previous_energy
            end
            imc += 1
        end                     # sweep_i1, sweep_i2
        
        #
        if mcsweep > nequil
            count += 1
            # taking_avg
            if count - oldcount == ndel
                outcount += 1

                si_sj .= calculate_si_sj(spins_lattice)
                total_s, s_square = calculate_average_sssq(spins_lattice)
                # total_s = sum(spins_lattice(:,:))
                # s_square = total_s**2
                #
                total_energy = total_energy_ising_2d(spins_lattice, jey1, jey2, magnetic_field)
                #
                average_s += total_s
                average_ssq += s_square
                avg_sisj .+= si_sj
                avgesq += total_energy^2
                avge += total_energy
                #
                oldcount = count
            end # taking_avg
        end
        mcsweep += 1
    end                         # end mc_loop while

    divide = outcount * 1.0
    avg_sisj .= avg_sisj / divide
    avge = avge / divide
    avgesq = avgesq / divide
    average_s = average_s / divide
    average_ssq = average_ssq / divide
    # write(unit_acc,'(3(f15.10,1x))')temperature,dble(ipass)/dble(imc),dble(ifail)/dble(imc)
    
    # ###########################################################
    # ###########################################################
    # call cpu_time(time2) #time2 = x05baf()#
    # write(*,"(a,2x,f15.10)",advance='no')'running program at t =',temperature
    # write(*,'(3x,"time taken =",f12.2," seconds")',advance='no')(time2-time1)
    
    # call cpu_time(tsf1)        #tsf1 = x05baf()
    @time sf_ising = structure_factor_from_si_sj(avg_sisj, lsize)
    # call cpu_time(tsf2)        #tsf2 = x05baf()
    max_sq = maximum(sf_ising)
    # write(*,'(3x,"time for s(q) =",1x,f12.2,1x,"seconds. maximum s(q) =",1x,f20.10)')tsf2-tsf1,max_sq
    println("count = ", count_temp, " Temperature = ", temperature)
    c_v = (avgesq - (avge)^2) / nsize
    c_v = c_v / (temperature^2)
    avge = avge / nsize
    
    chai = average_ssq / nsize - average_s^2 / nsize
    chai = chai / temperature
    sqaf = sf_ising[Int(lsize/2), Int(lsize/2)]
    magnetization = sf_ising[1, 1]
    #
    # write(unit_es,'(2(f20.10,1x),e20.10)')temperature,avge,c_v
    # write(unit_mx,'(6(f20.10,1x))')temperature,magnetization,sqaf,chai,average_s/dble(nsize),average_ssq/dble(nsize)
    # write(unit_sp,*)'t=',temperature
    for i1 = 1:lsize, i2 = 1:lsize
        # write(unit_sp,'(3(i2,1x))')i1,i2,spins_lattice(i1,i2)
        qstr = string(i1 - 1, "_", i2 - 1)
        # write(unit_sisj,'(3(i2,1x))')i1,i2,avg_sisj(1,indices(i1,i2))
        # write(unit_sq,'(a9,1x,1(f20.10,1x),e20.10)')qstr,temperature,sf_ising(i1-1,i2-1)
    end
    # write(unit_sp,*)
    # write(unit_sq,*)
    global temperature += -delta_t
    global count_temp += 1
end # do temperature_loop




# CLOSE(UNIT_ES)
# CLOSE(UNIT_MX)
# CLOSE(UNIT_SP)
# CLOSE(UNIT_SQ)
# CLOSE(UNIT_SISJ)
# CLOSE(UNIT_EPS)
# CLOSE(UNIT_TMP)

#   OPEN(101,FILE = 'spin_gs.dat')
#   OPEN(102,FILE = 'sisj_gs.dat')
#   DO I1 = 1, LSIZE
#      DO I2 = 1, LSIZE
#         WRITE(101,'(3(I2,1X))')I1,I2,SPINS_LATTICE(I1,I2)
#         WRITE(102,'(3(I2,1X))')I1,I2,SI_SJ(1,INDICES(I1,I2))
#      END DO;WRITE(101,*);WRITE(102,*)
#   END DO
#   CLOSE(101)
#   CLOSE(102)
#   PRINT*,'GROUND STATE SPIN SAVED IN FILE=','spin_gs.dat'
#   PRINT*,'LOWEST E',TOTAL_ENERGY_ISING_2D(SPINS_LATTICE,JEY1,JEY2,MAGNETIC_FIELD)/DBLE(NSIZE)
#   PRINT*,'NUMBER OF Ts =',COUNT_TEMP

# CONTAINS
# END PROGRAM COOL_ISING_SPINS_3D


# if energy_difference < 0.0
#     spins_lattice[i1, i2] = tmp_spin # update the temp spin
#     final_energy = previous_energy + energy_difference
#     ipass = ipass + 1
# else
#     metro_police = rand()
#     compare = energy_difference
#     compare = compare / (temperature + 1.0e-10)
#     compare = exp(-compare)
#     if metro_police < compare
#         spins_lattice[i1, i2] = tmp_spin # update the temp spin
#         final_energy = previous_energy + energy_difference
#         ipass += 1
#     else
#         final_energy = previous_energy
#         ifail += 1
#     end
# end


# energy_difference = 0.0
# xi = Vector{Int}([i1, i2]) #nearest neighbor part
# for m1 = 1:2, m2 = 1:2
#     co = xi
#     co[m1] = minusplus[xi[m1], m2]
#     energy_difference += 2.0 * jey1 * tmp_spin * spins_lattice[co[1], co[2]]
# end
# # xi = [i1,i2]
# for m1 = 1:2, m2 = 1:2
#     co = [minusplus[i1, m1], minusplus[i2, m2]]
#     energy_difference += 2.0 * jey2 * tmp_spin * spins_lattice[co[1], co[2]]
# end
# previous_s = spins_lattice[i1, i2]
# current_s = tmp_spin
# energy_difference += -magnetic_field * current_s + magnetic_field * previous_s
