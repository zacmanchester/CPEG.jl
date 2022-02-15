
struct GravityParameters
    μ::Float64
    J2::Float64
    R::Float64
    function GravityParameters()
        GM_MARS = 4.2828375816e13
        J2_MARS = 1960.45e-6
        R_MARS = 3386200.0
        new(GM_MARS,J2_MARS,R_MARS)
    end
end

@inline function gravity(g::GravityParameters,r::SVector{3, T}) where T

    nr = norm(r)

    x,y,z = r

    # precompute repeated stuff
    Re_r_sqr = 1.5*g.J2*(g.R/nr)^2
    five_z_sqr = 5*z^2/nr^2

    return  (-g.μ/nr^3)*SA[x*(1 - Re_r_sqr*(five_z_sqr - 1)),
                           y*(1 - Re_r_sqr*(five_z_sqr - 1)),
                           z*(1 - Re_r_sqr*(five_z_sqr - 3))]
end

# let
#
#     g = GravityParameters()
#
#     r = 1.0*3386200.0*normalize(randn(3))
#
#     r = SVector{3}(r)
#
#     @show norm(gravity(g,r))
#     # @btime gravity($g,$r)
#
# end
