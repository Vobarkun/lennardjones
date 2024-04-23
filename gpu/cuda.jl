function spreadBits(x)
    x &= 0x0000ffff                   # x = ---- ---- ---- ---- fedc ba98 7654 3210
    x = (x | (x << 8)) & 0x00ff00ff  # x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x | (x << 4)) & 0x0f0f0f0f  # x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x | (x << 2)) & 0x33333333  # x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x | (x << 1)) & 0x55555555  # x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    x + one(x)
end

function morton(x, y)
    (spreadBits(y - one(y)) << 1) + spreadBits(x - one(x)) - 0x2
end
morton(xy) = morton(xy...)

function sortedCellList(positions; maxind = 256)
    N = length(positions)
    cellinds = ((x, i) -> (mod1(morton(ceil.(Int32, x)), Int32(maxind)), i)).(positions, 0x1:Int32(N))
    sort!(cellinds, by = first)
    firsts = CUDA.zeros(Int32, maxind)
    ((o, (i, _), (j, _), n) -> (i != j && (o[i] = n); nothing)).(Ref(firsts), view(cellinds, 2:N), view(cellinds, 1:N-1), 2:N)
    CUDA.@allowscalar firsts[cellinds[1][1]] = 1
    cellinds, firsts
end

function sortbycells(positions; maxind = 256)
    N = length(positions)
    combined = map(positions, Int32(1):Int32(length(positions))) do p, i 
        cellind = mod1(morton(ceil.(Int32, p)), Int32(maxind))
        (cellind, i)
    end
    sort!(combined, by = first)
    cellinds = first.(combined)
    perm = last.(combined)

    firsts = CUDA.zeros(Int32, maxind)
    findfirsts!(firsts, cellinds)
    cellinds, perm, firsts
end


function kernel_first!(firsts, cellinds)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = gridDim().x * blockDim().x
    cellinds = CUDA.Const(cellinds)
    @inbounds for i = index:stride:length(cellinds)
        ci = cellinds[i]
        if i == 1 || cellinds[i-1] != ci
            firsts[ci] = i
        end
    end
    nothing
end

function findfirsts!(firsts, cellinds)
    @cuda threads=512 blocks=ceil(Int, length(cellinds)/512) kernel_first!(firsts, cellinds)
end

function ljf(σ, ε, r; b = 1.0)
    rr = r ⋅ r
    σr = (σ^2 / rr) ^ 3
    (σr > 1/250) * min(1e5, 24ε / rr * (2σr^2 - b * σr)) * r
end

function kernel_force!(fx, fy, positions, cellinds, firsts)
    index = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    stride = gridDim().y * blockDim().y

    o = SA[threadIdx().x % 3 - 1, (threadIdx().x - 1) ÷ 3 - 1]
    
    for i in index:stride:length(cellinds)
        p = positions[i]

        cind = mod1(morton(ceil.(Int32, p) + o), Int32(256))
        j = firsts[cind]
        j == 0 && continue

        f = SA[0.0f0,0.0f0]
        while j <= length(positions)
            q = positions[j]
            if cind == cellinds[j]
                if norm(p - q) < 1
                    f += -ljf(0.1, 1, q - p)
                end
            else
                break
            end
            j += 1
        end
        
        CUDA.@atomic fx[i] += f[1]
        CUDA.@atomic fy[i] += f[2]
    end
    nothing
end

function force!(fx, fy, positions, cellinds, firsts)
    fx .= 0; fy .= 0
    numblocks = (1, ceil(Int, length(positions)/64))
    @cuda threads=(9, 64) blocks=numblocks kernel_force!(fx, fy, positions, cellinds, firsts)
    fx, fy
end




function kernel_count!(counts, positions, cellinds, firsts)
    index = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    stride = gridDim().y * blockDim().y

    o = SA[threadIdx().x % 3 - 1, (threadIdx().x - 1) ÷ 3 - 1]
    
    for i in index:stride:length(cellinds)
        p = positions[i]

        cind = mod1(morton(ceil.(Int32, p) + o), Int32(256))
        j = firsts[cind]
        j == 0 && continue

        c = 0
        while j <= length(positions)
            q = positions[j]
            if cind == cellinds[j]
                if norm(p - q) < 1
                    c += 1
                end
            else
                break
            end
            j += 1
        end
        
        CUDA.@atomic counts[i] += c
    end
    nothing
end

function count!(counts, positions, cellinds, firsts)
    numblocks = (1, ceil(Int, length(positions)/64))
    @cuda threads=(9, 64) blocks=numblocks kernel_count!(counts, positions, cellinds, firsts)
    counts
end