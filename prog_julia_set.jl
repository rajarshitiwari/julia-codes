#!/usr/bin/env julia
const NUM_X = 500
const NUM_Y = 500
const MAX_COUNT = 10000
const NUM_COLORS = 255

function julia_set()
    
    ϕ = (sqrt(5.0) + 1.0) * 0.5
    ai = 1.0im
    # CONST = (-0.429D0,0.571D0)
    # CONST = (-0.54D0,0.55D0)
    # CONST = (-0.74543D0,0.11301D0)!0.518D0 * (-1.0D0,1.0D0)
    CC = 0.8*ai # (ϕ - 2.0D0) + AI * (ϕ - 1.0D0)
    radius = 2.0
    # OPEN(1,FILE='1.dat')
    # OPEN(2,FILE='2.dat')

    for ix = -NUM_X:NUM_X
        x = 1.45 * ix / NUM_X
        for iy = -NUM_Y:NUM_Y
        flag = 0
            y = 1.05 * iy / NUM_Y
            zed = x + ai * y
            zedsq = zed * zed
            count = 0
            while abs(zed) <= radius
                zedsq = zed * zed
                zed = zedsq + CC
                count += 1
                if count >= MAX_COUNT
                    flag = 1
                    break
                end
            end
            # WRITE(1,'(2(I4,1X),I6)')IX,IY,COUNT!MOD(COUNT,NUM_COLORS)
            if flag == 1
                # WRITE(2,'(2(F15.10,1X))')X,Y!MOD(COUNT,NUM_COLORS)
            end
        end
        # WRITE(1,*)
    end
end # PROGRAM JULIA_SET

julia_set()
