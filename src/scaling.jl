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
    a = v/(p.dscale/p.tscale^2)
    return v,a
end

@inline function unscale_va(p::Scaling,v::SVector{3,T},a::SVector{3,T}) where T
    v = v*(p.dscale/p.tscale)
    a = v*(p.dscale/p.tscale^2)
    return v,a
end
