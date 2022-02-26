# cd("/Users/kevintracy/.julia/dev/CPEG")
# Pkg.activate(".")
# cd("/Users/kevintracy/.julia/dev/CPEG/src")
# using LinearAlgebra
# using StaticArrays
# using ForwardDiff
# using SparseArrays
# using SuiteSparse
# using Printf
#
# include("qp_solver.jl")
# include("atmosphere.jl")
# include("scaling.jl")
# include("gravity.jl")
# include("aero_forces.jl")
# include("vehicle.jl")
# include("dynamics.jl")
# include("post_process.jl")

using LinearAlgebra
using StaticArrays
using ForwardDiff
using SparseArrays
using SuiteSparse
using MATLAB

function eg_mpc_quad(
    ev::CPEGWorkspace,
    A::Vector{SMatrix{7, 7, Float64, 49}},
    B::Vector{SMatrix{7, 1, Float64, 7}},
    X::Vector{SVector{7,Float64}},
    U::Vector{SVector{1,Float64}})

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
    dyn_eq = zeros((N)*nx) # N-1 dynamics constraints + N IC constraint
    # low_x0 = zeros(nx)
    # up_x0 = zeros(nx)
    up_tr = ev.cost.σ_tr*ones(N)
    low_tr = -ev.cost.σ_tr*ones(N)

    # cost function terms
    P = spzeros(nz,nz)
    q = zeros(nz)
    for i = 1:(N-1)
        P[idx_u[i],idx_u[i]] = [1.0]
        q[idx_u[i]] = [U[i][1]]
    end
    rr = normalize(ev.cost.rf/ev.scale.dscale)
    Qn = ev.cost.γ*(I - rr*rr')
    P[idx_x[N][1:3],idx_x[N][1:3]] = Qn'*Qn
    q[idx_x[N][1:3]] = -(Qn'*Qn)'*(ev.cost.rf/ev.scale.dscale - X[N][1:3])

    # Q = copy(P) #+ 1e-6*I
    # q = copy(q)
    # A = copy(A_eq)
    # b = copy(up_eq)
    G = [A_ineq;-A_ineq]
    h = [up_tr;-low_tr]

    z = quadprog(P,q,A_eq,dyn_eq,G,h; verbose = ev.solver_opts.verbose,
                                         atol = ev.solver_opts.atol,
                                    max_iters = ev.solver_opts.max_iters)
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

    return cU, norm(δu)

    # return cU

end

function main_cpeg(ev,x0_s,U)

    # vectors for storing trajectory information
    T = 20
    althist = [zeros(2) for i = 1:T]
    drhist = [zeros(2) for i = 1:T]
    crhist = [zeros(2) for i = 1:T]
    dunorm = zeros(T)
    for i = 1:T
        X,U = rollout(ev, x0_s, U)

        althist[i], drhist[i], crhist[i] = postprocess_scaled(ev,X,x0_s)
        A,B = get_jacobians(ev,X,U)
        U, ndu = eg_mpc_quad(ev,A,B,X,U)
        @show ndu
        if ndu<1e-2
            @info "success"
            println("miss distance: ",norm(X[end][1:3] - ev.cost.rf/ev.scale.dscale)*ev.scale.dscale/1e3," km")
            return U, althist[1:i], drhist[1:i], crhist[1:i]
        end
    end
    @error "CPEG failed"
end
# function main_cpeg(ev,x0_s,U,dt)
#
#     # vectors for storing trajectory information
#     T = 20
#     althist = [zeros(2) for i = 1:T]
#     drhist = [zeros(2) for i = 1:T]
#     crhist = [zeros(2) for i = 1:T]
#     dunorm = zeros(T)
#     for i = 1:T
#         X,U = rollout(ev, x0_s, U, dt)
#
#         althist[i], drhist[i], crhist[i] = postprocess_scaled(ev,X,x0_s)
#         A,B = get_jacobians(ev,X,U,dt)
#         U, ndu = eg_mpc_quad(ev,A,B,X,U)
#         @show ndu
#         if ndu<1e-2
#             @info "success"
#             println("miss distance: ",norm(X[end][1:3] - ev.cost.rf/ev.scale.dscale)*ev.scale.dscale/1e3," km")
#             return U, althist[1:i], drhist[1:i], crhist[1:i]
#         end
#     end
#     @error "CPEG failed"
# end
function tt()


        ev = CPEGWorkspace()
        ev.solver_opts.verbose = false

        Rm = ev.params.gravity.R
        r0 = SA[Rm+125e3, 0.0, 0.0] #Atmospheric interface at 125 km altitude
        V0 = 5.845e3 #Mars-relative velocity at interface of 5.845 km/sec
        γ0 = -15.474*(pi/180.0) #Flight path angle at interface
        v0 = V0*SA[sin(γ0), cos(γ0), 0.0]
        σ0 = deg2rad(3)

        r0sc,v0sc = scale_rv(ev.scale,r0,v0)

        # dt = 2.0/ev.scale.tscale

        x0_s = SA[r0sc[1],r0sc[2],r0sc[3],v0sc[1],v0sc[2],v0sc[3],σ0]


        N = 100
        U = [SA[0.0] for i = 1:N-1]
        #
        # X,U = rollout(ev, x0_s, U, dt)
        # @show length(X)

        # A,B= get_jacobians(ev,X,U,dt)

        # @show cond(A[1])
        # @show A[2]

        # @show B[2]

        Rf = Rm+10.0e3 #Parachute opens at 10 km altitude
        rf = Rf*SA[cos(7.869e3/Rf)*cos(631.979e3/Rf); cos(7.869e3/Rf)*sin(631.979e3/Rf); sin(7.869e3/Rf)]
        # @show Rf
        # @show norm(rf)
        # rf_s = rf/ev.scale.dscale


        # eg_mpc_quad(ev,A,B,X,U,rf_s)


        @info "cPEG TIME"

        # qp solver settings
        ev.solver_opts.verbose = false
        ev.solver_opts.atol = 1e-5
        ev.solver_opts.max_iters = 50

        # MPC stuff
        ev.cost.σ_tr = deg2rad(20)
        ev.cost.rf = rf
        ev.cost.γ = 1e3

        # sim stuff
        ev.dt = 2.0

        Ucpeg, althist, drhist, crhist = main_cpeg(ev,x0_s,U)

        xf_dr, xf_cr = rangedistances(ev,rf,SVector{6}([r0;v0]))


        num2plot = float(length(althist))
        # plot_groundtracks(drhist/1e3,crhist/1e3,althist/1e3,xf_dr/1e3,xf_cr/1e3,num2plot,"quad")
        # @show xf_dr/1e3
        # @show xf_cr/1e3
        # X,U = rollout(ev,x0_s,Ucpeg,dt)
        # X = unscale_X(ev.scale,X)
        # # @show typeof(X[1][SA[1,2,3]])
        #
        # alt, dr, cr = postprocess(ev,X,[r0;v0])
        #
        # σ = [X[i][7] for i = 1:length(X)]
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
        #
        # mat"
        # figure
        # hold on
        # title('Bank Angle')
        # plot(rad2deg($σ))
        # hold off
        # "

        return nothing
end

function plot_groundtracks(drhist,crhist,althist,xf_dr,xf_cr,num2plot,id)

    mat"
    figure
    hold on
    rgb1 = [29 38 113]/255;
    rgb2 = 1.3*[195 55 100]/255;
    drgb = rgb2-rgb1;
    cmap = [];
    for i = 1:round($num2plot)
        px = $drhist{i};
        py = $crhist{i};
        plot(px,py,'Color',rgb1 + drgb*(i)/($num2plot),'linewidth',3)
        cmap = [cmap;rgb1 + drgb*(i)/($num2plot)];
        %plot(px(1),py(1),'r.','markersize',20)
    end
    plot($xf_dr,$xf_cr,'g.','markersize',20)
    colormap(cmap);
    pp = colorbar;
    pp.Ticks = 0:(1/$num2plot):1;
    pp.Location = 'northoutside';
    pp.TickLabels = 0:round($num2plot);
    xlabel('downrange (km)')
    ylabel('crossrange (km)')
    %xlim([250,700])
    %ylim([0,20])
    hold off
    fleg = legend('figure()');
    set(fleg,'visible','off')
    %addpath('/Users/kevintracy/devel/WiggleSat/matlab2tikz-master/src')
    %matlab2tikz(strcat('cpeg_examples/bank_angle/tikz/',$id,'_crdr.tex'))
    %close all
    "

    # this one is for plotting
    mat"
    figure
    hold on
    rgb1 = [29 38 113]/255;
    rgb2 = 1.3*[195 55 100]/255;
    drgb = rgb2-rgb1;
    cmap = [];
    for i = 1:round($num2plot)
        px = $drhist{i};
        alt = $althist{i};
        plot(px,alt,'Color',rgb1 + drgb*(i)/($num2plot),'linewidth',3)
        cmap = [cmap;rgb1 + drgb*(i)/($num2plot)];
    end
    %plot($xf_dr,$xf_cr,'g.','markersize',20)
    colormap(cmap);
    pp = colorbar;
    pp.Ticks = 0:(1/$num2plot):1;
    pp.Location = 'northoutside';
    pp.TickLabels = 0:round($num2plot);
    plot([0,800],ones( 2,1)*10,'r' )
    plot($xf_dr,10,'g.','markersize',20)
    %xlim([400,700])
    %ylim([8,30])
    xlabel('downrange (km)')
    ylabel('altitude (km)')
    hold off
    %saveas(gcf,'alt.png')
    fleg = legend('figure()');
    set(fleg,'visible','off')
    %addpath('/Users/kevintracy/devel/WiggleSat/matlab2tikz-master/src')
    %matlab2tikz(strcat('cpeg_examples/bank_angle/tikz/',$id,'_altdr.tex'))
    %matlab2tikz('bankaoa_alt.tex')
    %close all
    "
    return nothing
end

tt()
