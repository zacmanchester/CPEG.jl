using LinearAlgebra, MATLAB, ForwardDiff, SparseArrays, Printf
using Convex, Hypatia
const FD = ForwardDiff
include("/Users/kevintracy/.julia/dev/CPEG/src/qp_solver.jl")


function dynamics(x,u)
    r = x[1:3]
    v = x[4:6]
    return [v;u[1:3] + [0;0;-9.81]]
end
function rk2(x,u)
    Δt = u[4]
    k1 = Δt*dynamics(x,u)
    k2 = Δt*dynamics(x + k1/2, u)
    return x + k2
end
function rollout(x0,U)
    N = length(U)+1
    X = [zeros(length(x0)) for i = 1:N]
    X[1] = x0
    for i = 1:(N-1)
        X[i+1] = rk2(X[i],U[i])
    end
    return X
end
function get_jacobians(X,U)
    N = length(X)
    A = [FD.jacobian(_x -> rk2(_x,U[i]),X[i]) for i = 1:N-1]
    B = [FD.jacobian(_u -> rk2(X[i],_u),U[i]) for i = 1:N-1]
    return A,B
end

function mpc(X,U)
    A,B = get_jacobians(X,U)
    nx = length(X[1])
    nu = length(U[1])
    N = length(X)
    alt_target = 10.0
    dt_target = 0.5

    δX = Variable(nx,N)
    δU = Variable(nu,N-1)

    J = 0

    γ = 1e3
    for i = 1:N-1
        J += γ*quadform(U[i][1:3] + δU[1:3,i],Matrix(I(3)))
        J += square(U[i][4] + δU[4,i] - dt_target)
    end

    p= minimize(J)

    for i = 1:N-1
        p.constraints += δU[4,i] >= 0.3
        p.constraints += δU[4,i] <= 0.7
    end
    p.constraints += X[N][3] + δX[3,N] == alt_target
    p.constraints += δX[:,1] == zeros(6)
    for i = 1:N-1
        # @info "here?"
        p.constraints += δX[:,i+1] == A[i]*δX[:,i] + B[i]*δU[:,i]
        # @info "nope"
    end

    solve!(p, Hypatia.Optimizer();silent_solver = true)

    δU = evaluate(δU)
    # @show norm(δU)
    U = [U[i] + vec(δU[:,i]) for i = 1:length(U)]

    return U, norm(δU)
end

function mat_from_vec(X)
    Xm = zeros(length(X[1]),length(X))
    for i = 1:length(X)
        Xm[:,i] = X[i]
    end
    return Xm
end
function tt()



    x0 = [0;0;100;0;0;-10.0]
    # u0 = [0;0;9.7]
    u1 = [0;0;9.81;0.5]

    N = 20

    U = [copy(u1) for i = 1:N-1]

    X = rollout(x0,U)

    xm = mat_from_vec(X)


    for i = 1:5

        U, ndu = mpc(X,U)

        if ndu < 1e-1
            break
        end
        X = rollout(x0,U)
        @show X[end][3]

    end



    mat"
    figure
    hold on
    %plot3($xm(1,:),$xm(2,:),$xm(3,:))
    plot($xm(1:3,:)')
    hold off
    "


    return nothing
end

tt()
