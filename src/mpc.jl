

using LinearAlgebra
using StaticArrays
using ForwardDiff
using SparseArrays
using SuiteSparse
using MATLAB

function eg_mpc_quad(
    ev::EntryVehicle,
    A::Vector{SMatrix{7, 7, Float64, 49}},
    B::Vector{SMatrix{7, 1, Float64, 7}},
    X::Vector{SVector{7,Float64}},
    U::Vector{SVector{1,Float64}},
    rf_s::SVector{3,Float64})

    # sizes for state and control
    nx = 7
    nu = 1
    N = length(X)

    # indicees for state and control
    idx_x = [(i-1)*(nx+nu) .+ (1:nx) for i = 1:length(X)]
    idx_u = [(i-1)*(nx+nu) .+ nx .+ (1:nu) for i = 1:(length(X)-1)]

    # constraint jacobian (in this form L≦Az≤U)
    nz = (N*nx) + (N-1)*nu
    nc = (N-1)*nx + nx #+ (N)*1
    A_eq = spzeros(nc,nz)
    A_ineq = spzeros(N,nz) # state constraints

    # dynamics constraint (equality)
    idx_c = [(i-1)*(nx) .+ (1:nx) for i = 1:(N-1)]
    for i = 1:(N-1)
        A_eq[idx_c[i],idx_x[i]]   = A[i]
        A_eq[idx_c[i],idx_u[i]]   = B[i]
        A_eq[idx_c[i],idx_x[i+1]] = -I(nx)
    end
    A_eq[(N-1)*nx .+ (1:nx), idx_x[1]] = I(nx)

    # state constraints on δσ (inequality)
    for i = 1:N
        A_ineq[i,idx_x[i][7]] = 1
    end

    # constraint bounds
    low_dyneq = zeros((N-1)*nx)
    up_dyneq = copy(low_dyneq)
    low_x0 = zeros(nx)
    up_x0 = zeros(nx)
    up_tr = deg2rad(10)*ones(N)
    low_tr = -deg2rad(10)*ones(N)

    low_eq = [low_dyneq;low_x0]
    up_eq = [up_dyneq;up_x0]
    # stack everything up

    A = [A_eq;A_ineq]
    Lo = [low_eq;low_tr]
    Up = [up_eq;up_tr]

    # cost function terms
    R = 1
    P = spzeros(nz,nz)
    q = zeros(nz)
    for i = 1:(N-1)
        P[idx_u[i],idx_u[i]] = [R]
        q[idx_u[i]] = [R*U[i][1]]
    end
    rr = normalize(rf_s)
    Qn = 1000*(I - rr*rr')
    P[idx_x[N][1:3],idx_x[N][1:3]] = Qn'*Qn
    q[idx_x[N][1:3]] = -(Qn'*Qn)'*(rf_s - X[N][1:3])

    Q = copy(P) #+ 1e-6*I
    q = copy(q)
    A = copy(A_eq)
    b = copy(up_eq)
    G = [A_ineq;-A_ineq]
    h = [up_tr;-low_tr]

    z = quadprog(Q,q,A,b,G,h)
    δx = [z[idx_x[i]] for i = 1:(N)]
    δu = [z[idx_u[i]] for i = 1:(N-1)]

    # osqp = OSQP.Model()
    # OSQP.setup!(osqp; P=P, q=q, A=A, l=(Lo), u=(Up), eps_abs = 1e-8, eps_rel = 1e-8, max_iter = 10000,polish = 1)
    # results = OSQP.solve!(osqp)

    # # recover solution
    # δx = [results.x[idx_x[i]] for i = 1:(N)]
    # δu = [results.x[idx_u[i]] for i = 1:(N-1)]

    # apply the δ's
    cX = X + δx
    cU = U + δu

    # return cU, norm(δu)

    return nothing

end

function tt()


        ev = EntryVehicle()
        Rm = ev.params.gravity.R
        r0 = SA[Rm+125e3, 0.0, 0.0] #Atmospheric interface at 125 km altitude
        V0 = 5.845e3 #Mars-relative velocity at interface of 5.845 km/sec
        γ0 = -15.474*(pi/180.0) #Flight path angle at interface
        v0 = V0*SA[sin(γ0), cos(γ0), 0.0]
        σ0 = deg2rad(3)

        r0sc,v0sc = scale_rv(ev.scale,r0,v0)

        dt = 2.0/ev.scale.tscale

        x0 = SA[r0sc[1],r0sc[2],r0sc[3],v0sc[1],v0sc[2],v0sc[3],σ0]


        N = 100
        U = [SA[0.0] for i = 1:N-1]
        #
        X,U = rollout(ev, x0, U, dt)
        @show length(X)

        A,B= get_jacobians(ev,X,U,dt)

        @show cond(A[1])
        @show A[2]

        @show B[2]

        Rf = Rm+10.0e3 #Parachute opens at 10 km altitude
        rf = Rf*SA[cos(7.869e3/Rf)*cos(631.979e3/Rf); cos(7.869e3/Rf)*sin(631.979e3/Rf); sin(7.869e3/Rf)]
        rf_s = rf/ev.scale.dscale


        eg_mpc_quad(ev,A,B,X,U,rf_s)
        # X = unscale_X(ev.scale,X)
        # # @show typeof(X[1][SA[1,2,3]])
        #
        # alt, dr, cr = postprocess(ev,X,[r0;v0])
        #
        #
        # mat"
        # figure
        # hold on
        # plot($alt/1e3)
        # hold off "
        #
        # mat"
        # figure
        # hold on
        # plot($dr/1e3,$cr/1e3)
        # hold off"

        return nothing
end



tt()
