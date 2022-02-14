
function dynamics(ev::EntryVehicle, x::SVector{7,T}, u::SVector{1,W}) where {T,W}

    # scaled variables
    r_scaled = x[SA[1,2,3]]
    v_scaled = x[SA[4,5,6]]
    σ = x[7]

    # unscale
    r, v = unscale(ev.scale,r_scaled,v_scaled)

    # density
    ρ = density(ev.params.density, r)


    L, D = LD_mags(ev,r,v)

    e1, e2 = e_frame(r,v)

    D_a = -(D/norm(v))*v
    L_a = L*sin(σ)*e1 + L*cos(σ)*e2

    g = gravity(ev.params.gravity,r)

    ω = ev.planet.ω
    a = D_a + L_a + g - 2*cross(ω,v) - cross(ω,cross(ω,r))

    # rescale units
    v,a = scale_va(ev.scale,v,a)

    return SA[v[1],v[2],v[3],a[1],a[2],a[3],u[1]*ev.scale.uscale]
end

# function f(x::SVector{7,T},u::SVector{1,W})::SVector{7,T} where {T,W}


function rk4(
    ev::EntryVehicle,
    x_n::SVector{7,T},
    u::SVector{1,W},
    dt::Float64) where {T,W}

    k1 = dt*dynamics(ev,x_n,u)
    k2 = dt*dynamics(ev,x_n+k1/2,u)
    k3 = dt*dynamics(ev,x_n+k2/2,u)
    k4 = dt*dynamics(ev,x_n+k3,u)
    return (x_n + (1/6)*(k1 + 2*k2 + 2*k3 + k4))
end

function rollout(ev::EntryVehicle,x0::SVector{7,T},U_in::Vector{SVector{1,T}},dt::Float64) where T
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

function get_jacobians(
    ev::EntryVehicle,
    X::Vector{SVector{7,T}},
    U::Vector{SVector{1,T}},
    dt::Float64
    ) where T

    N = length(X)
    A = [ForwardDiff.jacobian(_x->rk4(ev,_x,U[i],dt),X[i]) for i = 1:N-1]
    B = [ForwardDiff.jacobian(_u->rk4(ev,X[i],_u,dt),U[i]) for i = 1:N-1]
    return A,B
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

    i = 7
    @show ForwardDiff.jacobian(_x->rk4(ev,_x,U[i],dt),X[i])
    @show ForwardDiff.jacobian(_u->rk4(ev,X[i],_u,dt),U[i])
    # A,B= get_jacobians(ev,X,U,dt)



end
