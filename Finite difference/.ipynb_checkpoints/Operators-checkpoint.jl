function diff!(du, u, axis, idx, idx_m, idx_p, Δu)
    fill!(du, 0);
    if size(u, axis)>=2
        for i in 2:size(u, axis)-1
            idx[axis]= i;
            idx_m[axis]= i-1;
            idx_p[axis]= i+1;
            du[idx...] .= (u[idx_p...].- u[idx_m...])./2Δu; # (CE) replace later by EF, EB, CEF
        end
    else
        if axis!=ESCAPE_AXIS @warn "Can't differentiate, size of array too small (<2)" end
    end
end

function diff2!(d2u, u, axis, idx, idx_p, idx_m, Δu)
    fill!(d2u, 0);
    if size(u, axis)>=2
        for i in 2:size(u, axis)-1
            idx[axis]= i;
            idx_p[axis]= i+1;
            idx_m[axis]= i-1;
            d2u[idx...] .= (u[idx_p...].- 2 .*u[idx...] .+ u[idx_m...])./(Δu^2); #  Staggered grid CE      replace later by EF, EB, CEF
        end
    else
        if axis!=ESCAPE_AXIS @warn "Can't differentiate, size of array too small (<2)" end
    end
end

# `d_x!` should take in `nx ny nz` and return `nx ny nz`

function d_z!(du, u, Δz)
    diff!(du, u, 1, [1,:,:], [1,:,:], [1,:,:], Δz);
end

function d_y!(du, u, Δy)
    diff!(du, u, 2, [:,1,:], [:,1,:], [:,1,:], Δy);
end

function d_x!(du, u, Δx)
    diff!(du, u, 3, [:,:,1], [:,:,1], [:,:,1], Δx);
end

# we might not even need t so stop worrying about it for the moment
function d_t!(du, u, component, Δt)
    diff!(du, u, 4, [:,:,:,1, component], [:,:,:,1, component], Δt);
end


# later replace the following by a combination of d_x, d_y, d_z, d_t
function d2_z!(d2u, u, Δz)
    diff2!(d2u, u, 1, [1,:,:], [1,:,:], [1,:,:], Δz);
end

function d2_y!(d2u, u, Δy)
    diff2!(d2u, u, 2, [:,1,:], [:,1,:], [:,1,:], Δy);
end

function d2_x!(d2u, u, Δx)
    diff2!(d2u, u, 3, [:,:,1], [:,:,1], [:,:,1], Δx);
end

# no worrying about time now
function d2_t!(d2u, u, component, Δt)
    diff2!(d2u, u, 4, [:,:,:,1, component], [:,:,:,1, component], [:,:,:,1, component], Δt);
end


function div!(du, du_temp, u, Δx, Δy, Δz)
    fill!(du, 0.);
    d_x!(du_temp, u[:,:,:,1], Δx);
    broadcast!(+, du, du, du_temp);
    d_y!(d2u_temp, u[:,:,:,2], Δy);
    broadcast!(+, du, du, du_temp);
    d_z!(d2u_temp, u[:,:,:,3], Δz);
    broadcast!(+, du, du, du_temp);
end

# u= nx ny nz 3
# du= nx ny nz 3

function curl!(du, du_p, du_m, u, Δx, Δy, Δz)
    # @show size.([du, du_p, du_m, u])
    
    fill!(du, 0.);
    fill!(du_p, 0.);
    fill!(du_m, 0.);
    
    d_y!(du_p, u[:,:,:,3], Δy); # ∂zy
    d_z!(du_m, u[:,:,:,2], Δz); # ∂yz
    # broadcast!(-, du[:,:,:,1], du_p, du_m);
    du[:,:,:,1].= du_p .- du_m;
    # @show norm(du_p), norm(du_m), norm(du[:,:,:,1])
    
    d_z!(du_p, u[:,:,:,1], Δz); # ∂xz
    d_x!(du_m, u[:,:,:,3], Δx); # ∂zx
    # broadcast!(-, du[:,:,:,2], du_p, du_m);
    du[:,:,:,2].= du_p .- du_m;
    # @show norm(du_p), norm(du_m), norm(du[:,:,:,2])
    
    d_x!(du_p, u[:,:,:,2], Δx); # ∂yx
    d_y!(du_m, u[:,:,:,1], Δy); # ∂xy
    # broadcast!(-, du[:,:,:,3], du_p, du_m);
    du[:,:,:,3].= du_p .- du_m;
    # @show norm(du_p), norm(du_m), norm(du[:,:,:,3])
    
end

# u= nx ny nz
# d2u= nx ny nz

function laplacian!(d2u, d2u_temp, u, Δx, Δy, Δz)
    fill!(d2u, 0.);
    fill!(d2u_temp, 0.);
    d2_x!(d2u_temp, u, Δx);
    broadcast!(+, d2u, d2u, d2u_temp);
    d2_y!(d2u_temp, u, Δy);
    broadcast!(+, d2u, d2u, d2u_temp);
    d2_z!(d2u_temp, u, Δz);
    broadcast!(+, d2u, d2u, d2u_temp);
end