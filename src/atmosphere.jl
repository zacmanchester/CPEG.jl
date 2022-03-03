
struct CPEGDensityParameters
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
    f::Float64
    g::Float64
    h::Float64
    i::Float64
    function CPEGDensityParameters()
        # default values from a marsgram sample
        p = [ -4.001833776317166
              -0.015696377977412
              -0.129412788692219
               0.000199820058253
               0.010377518080309
               0.000043652882189
              -0.000539682362508
              -0.000000205086106
               0.000000874980179]
        new(p...)
    end
end

function density(p::CPEGDensityParameters, h::T) where T
   h = h / 1000

   # clamp in a forward diff friendly way
   if h > 125.0
       h = one(h)*125.0
   elseif h < 0.2
       h = one(h)*0.2
   end
   num = @evalpoly(h, p.a, p.c, p.e, p.g, p.i)
   den = @evalpoly(h, 1.0, p.b, p.d, p.f, p.h)
   exp(num/den)
end
