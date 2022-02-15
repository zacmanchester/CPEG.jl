
struct Scaling
    dscale::Float64
    tscale::Float64
    uscale::Float64
    function Scaling()
        new(1e6,1e3,1e3)
    end
end

@inline function scale_rv(p::Scaling,r::SVector{3,T},v::SVector{3,T}) where T
    r = r/p.dscale
    v = v/(p.dscale/p.tscale)
    return r,v
end

@inline function unscale_rv(p::Scaling,r::SVector{3,T},v::SVector{3,T}) where T
    r = r*p.dscale
    v = v*(p.dscale/p.tscale)
    return r,v
end

@inline function scale_va(p::Scaling,v::SVector{3,T},a::SVector{3,T}) where T
    v = v/(p.dscale/p.tscale)
    a = a/(p.dscale/p.tscale^2)
    return v,a
end

@inline function unscale_va(p::Scaling,v::SVector{3,T},a::SVector{3,T}) where T
    v = v*(p.dscale/p.tscale)
    a = a*(p.dscale/p.tscale^2)
    return v,a
end


function unscale_X(p::Scaling,X::Vector{SVector{7,Float64}})
    N = length(X)
    Xu = [@SVector zeros(7) for i = 1:N]
    for i = 1:N
        r,v = unscale_rv(p,X[i][SA[1,2,3]],X[i][SA[4,5,6]])
        σ = X[i][7]
        Xu[i] = SVector{7}([r;v;σ])
    end
    return Xu
end
