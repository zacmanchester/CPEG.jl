
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

@inline function scale_rv(ev::EntryVehicle,r::SVector{3,T},v::SVector{3,T})::Tuple{SVector{3,T},SVector{3,T}} where T
    r = r/ev.scale.dscale
    v = v/(ev.scale.dscale/ev.scale.tscale)
    return r,v
end

@inline function unscale_rv(ev::EntryVehicle,r::SVector{3,T},v::SVector{3,T})::Tuple{SVector{3,T},SVector{3,T}} where T
    r = r*ev.scale.dscale
    v = v*(ev.scale.dscale/ev.scale.tscale)
    return r,v
end

@inline function scale_va(ev::EntryVehicle,v::SVector{3,T},a::SVector{3,T})::Tuple{SVector{3,T},SVector{3,T}} where T
    v = v/(ev.scale.dscale/ev.scale.tscale)
    a = v/(ev.scale.dscale/ev.scale.tscale^2)
    return v,a
end

@inline function unscale_va(ev::EntryVehicle,v::SVector{3,T},a::SVector{3,T})::Tuple{SVector{3,T},SVector{3,T}} where T
    v = v*(ev.scale.dscale/ev.scale.tscale)
    a = v*(ev.scale.dscale/ev.scale.tscale^2)
    return v,a
end

function dynamics(ev::EntryVehicle, x::SVector{7,T}, u::SVector{1,W})::SVector{7,T} where {T,W}

    # scaling numbers
    dscale = ev.scale.dscale
    tscale = ev.scale.tscale
    uscale = ev.scale.uscale

    # scaled variables
    r = x[SA[1,2,3]]
    v = x[SA[4,5,6]]
    σ = x[7]

    # unscale
    r,v=unscale(ev,r,v)

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

    return SA[v[1],v[2],v[3],v̇[1],v̇[2],v̇[3],u[1]*uscale]
end

function rk4(ev::EntryVehicle,x_n::SVector{7,T},u::SVector{1,W},dt::Float64)::SVector{7,T} where {T,W}
    k1 = dt*dynamics(ev,x_n,u)
    k2 = dt*dynamics(ev,x_n+k1/2,u)
    k3 = dt*dynamics(ev,x_n+k2/2,u)
    k4 = dt*dynamics(ev,x_n+k3,u)
    return (x_n + (1/6)*(k1 + 2*k2 + 2*k3 + k4))
end

function rollout(ev::EntryVehicle,x0::SVector{7,T},U_in::Vector{SVector{1,T}},dt::Float64)::Tuple{Vector{SVector{7,T}},Vector{SVector{1,T}}} where T
    N = 1000
    X = [@SVector zeros(length(x0)) for i = 1:N]
    U = [@SVector zeros(length(U_in[1])) for i = 1:N]
    X[1] = x0
    end_idx = NaN
    for i = 1:N-1
        U[i] = (i>length(U_in)) ? U_in[end] : U_in[i]

        X[i+1] = rk4(ev,X[i],U[i],dt)

        @show (norm(X[i+1][1:3])*ev.scale.dscale - ev.params.gravity.R)/1e3

        if (norm(X[i+1][1:3])*ev.scale.dscale - ev.params.gravity.R) < 10e3
            @info "hit alt"
            @show i
            end_idx = i+1
            break
        end
    end

    # trim relavent
    if isnan(end_idx)
        error("didn't hit the altitude during the rollout")
    end

    X = X[1:end_idx]
    U = U[1:(end_idx-1)]
    return X, U
end


let

    ev = EntryVehicle()

    # # @show ev.scale.uscale
    #
    # x = @SVector randn(7)
    # u = SA[1.45]
    # #
    # @btime dynamics($ev,$x,$u)
    #
    # @btime rk4($ev,$x,$u,0.1)
    # #
    # A = ForwardDiff.jacobian(_x -> dynamics(ev,_x,u),x)
    # # @btime ForwardDiff.jacobian(_x -> dynamics($ev,_x,$u),$x)
    # @btime ForwardDiff.jacobian(_x -> rk4($ev,_x,$u,0.1),$x)
    # ForwardDiff.jacobian(_x -> rk4(ev,_x,u,0.1),x)
    Rm = ev.params.gravity.R
    r0 = SA[Rm+125e3, 0.0, 0.0] #Atmospheric interface at 125 km altitude
    V0 = 5.845e3 #Mars-relative velocity at interface of 5.845 km/sec
    γ0 = -15.474*(pi/180.0) #Flight path angle at interface
    v0 = V0*SA[sin(γ0), cos(γ0), 0.0]
    σ0 = deg2rad(90)

    r0,v0 = scale_rv(ev,r0,v0)

    dt = 1.0/ev.scale.tscale

    x0 = SA[r0[1],r0[2],r0[3],v0[1],v0[2],v0[3],σ0]

    # @show norm(r0)
    #
    # @show norm(r0)*ev.scale.dscale - ev.params.gravity.R
    N = 100
    U = [SA[0.0] for i = 1:N-1]
    #
    X,U = rollout(ev, x0, U, dt)
    # N = 100
    # X = [@SVector zeros(7) for i = 1:N]

    # @show typeof([r0;v0;σ0])
    # X[1] = SA[r0[1],r0[2],r0[3],v0[1],v0[2],v0[3],σ0]
    # @btime $X[1] = SA[$r0[1],$r0[2],$r0[3],$v0[1],$v0[2],$v0[3],$σ0]

    # @show X[1]
    # @show typeof(X[1])



end
