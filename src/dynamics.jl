
struct Parameters
    density::DensityParameters
    gravity::GravityParameters
    aero::AeroParameters
    function Parameters()
        new(DensityParameters(),GravityParameters(),AeroParameters())
    end
end

struct Scaling
    dscale::Float64
    tscale::Float64
    uscale::Float64
    function Scaling()
        new(1e6,1e3,1e3)
    end
end

struct Planet
    Ω::Float64
    ω::SVector{3,Float64}
    function Planet()
        Ω = 7.088218127774194e-5 # rad/s
        ω = SA[0,0,Ω]
        return new(Ω,ω)
    end
end

struct EntryVehicle
    params::Parameters
    scale::Scaling
    planet::Planet
    function EntryVehicle()
        new(Parameters(),Scaling(),Planet())
    end
end

function dynamics(ev::EntryVehicle, x::SVector{7,T}, u::SVector{1,T})::SVector{7,T} where T
#
    # dscale, tscale, uscale = ev.scale.dscale, ev.scale,tscale, ev.scale,uscale
    dscale = ev.scale.dscale
    tscale = ev.scale.tscale
    uscale = ev.scale.uscale

    # r = x[1:3]
    r = x[SA[1,2,3]]
    v = x[SA[4,5,6]]
    σ = x[7]

    # ρ = density(ev,r)
    L,D = LD_mags(ev,r,v)

    e1,e2 = e_frame(r,v)

    D_a = -(D/norm(v))*v
    L_a = L*sin(σ)*e1 + L*cos(σ)*e2

    g = gravity(ev.params.gravity,r)

    ω = ev.planet.ω
    v̇ = D_a + L_a + g - 2*cross(ω,v) - cross(ω,cross(ω,r))

    # rescale units
    v = v/(dscale/tscale)
    v̇ = v̇/(dscale/tscale^2)

#
#     r = SA[]

    # return SA[1,2,3,4,5,6,7.0]
    return SA[v[1],v[2],v[3],v̇[1],v̇[2],v̇[3],u[1]*uscale]
end

let

    ev = EntryVehicle()

    # @show ev.scale.uscale

    x = @SVector randn(7)
    u = 1.45

    @btime dynamics($ev,$x,$u)

    # A = ForwardDiff.jacobian(_x -> dynamics(ev,_x,u),x)
    @btime ForwardDiff.jacobian(_x -> dynamics($ev,_x,$u),$x)

end
