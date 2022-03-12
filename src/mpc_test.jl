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
    X::Vector{SVector{7,Float64}})

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
    up_tr = ev.cost.σ_tr*ones(N)
    low_tr = -ev.cost.σ_tr*ones(N)

    # cost function terms
    P = spzeros(nz,nz)
    q = zeros(nz)
    R = 0.01
    for i = 1:(N-1)
        P[idx_u[i],idx_u[i]] = [R]
        q[idx_u[i]] = [R*ev.U[i][1]]
    end

    # projection matrix onto the landing plane
    Qn = ev.cost.γ*landing_plane_proj(ev.cost.rf)

    # cost terms for terminal cost
    P[idx_x[N][1:3],idx_x[N][1:3]] = Qn'*Qn
    q[idx_x[N][1:3]] = -(Qn'*Qn)*(ev.cost.rf/ev.scale.dscale - X[N][1:3])

    # setup standard form QP
    P += 1e-6*I
    G = [A_ineq;-A_ineq]
    h = [up_tr;-low_tr]

    # solve the equality only constrained QP
    sol = lu([P A_eq';A_eq spzeros(N*nx,N*nx)])\[-q;dyn_eq]
    z = sol[1:length(q)]
    qp_iters = 1

    # if this violates the inequality constraints, then we send it to quadprog
    if sum(G*z .> h) != 0
        # solve with QP solver (src/qp_solver.jl)
        z, qp_iters = quadprog(P,q,A_eq,dyn_eq,G,h; verbose = ev.qp_solver_opts.verbose,
                                             atol = ev.qp_solver_opts.atol,
                                        max_iters = ev.qp_solver_opts.max_iters)
    end

    # pull out δx and δu from z
    δx = [z[idx_x[i]] for i = 1:(N)]
    δu = [z[idx_u[i]] for i = 1:(N-1)]

    # update U = U + δu
    ev.U += δu

    return norm(δu), qp_iters
end

function landing_plane_proj(rf)
    rr = normalize(rf)
    I - rr*rr'
end
function update_miss_distance!(ev::CPEGWorkspace,X::Vector{SVector{7,Float64}})
    ev.miss_distance = miss_distance(ev,X)
    return nothing
end
function miss_distance(ev::CPEGWorkspace,X::Vector{SVector{7,Float64}})
    Qn = landing_plane_proj(ev.cost.rf)
    rf_sim = X[end][1:3]*ev.scale.dscale
    norm(Qn*(rf_sim - ev.cost.rf))
end

function update_σ!(ev::CPEGWorkspace,X::Vector{SVector{7,Float64}})
    ev.σ = [X[i][7] for i = 1:length(X)]
end

function main_cpeg(ev, r, v, σ)

    rs, vs = scale_rv(ev.scale,r,v)
    x0_s = SVector{7}([rs;vs;σ])

    old_md = NaN
    new_md = NaN
    ndu = NaN
    qp_iters = NaN

    althist = [zeros(2) for i = 1:ev.max_cpeg_iter]
    drhist = [zeros(2) for i = 1:ev.max_cpeg_iter]
    crhist = [zeros(2) for i = 1:ev.max_cpeg_iter]
    σhist = [zeros(2) for i = 1:ev.max_cpeg_iter]

    if ev.verbose
        @printf "iter    miss (km)    Δmiss (km)     |ΔU|        N   QP iters\n"
        @printf "----------------------------------------------------------\n"
    end

    for i = 1:ev.max_cpeg_iter

        # rollout and check miss distance
        X = rollout(ev, x0_s)
        althist[i], drhist[i], crhist[i] = postprocess_scaled(ev,X,x0_s)
        σhist[i]=[X[i][7] for i = 1:length(X)]
        old_md = new_md
        new_md = miss_distance(ev,X)

        if ev.verbose && i>1
            @printf("%3d     %9.3e    %+9.3e    %9.3e   %3d  %3d\n",
                    i, new_md/1e3,(new_md - old_md)/1e3, ndu, length(X), qp_iters)
        end

        # check the flipping
        # flipping = (max(old_md,new_md)<5e3) & (i>5) & (qp_iters == 1)

        # check termination criteria
        if (abs(new_md - old_md) < ev.miss_distance_tol) | (ndu < ev.ndu_tol) | (i == ev.max_cpeg_iter)#| flipping
            update_σ!(ev,X)
            return althist[1:i], drhist[1:i], crhist[1:i], σhist[1:i]
        end

        # linearize and solve QP
        A,B = get_jacobians(ev,X)
        ndu, qp_iters = eg_mpc_quad(ev,A,B,X)
    end
    error("CPEG failed: reached max iters")
end

function scale_state(ev,r,v,σ)
    rs,vs = scale_rv(ev.scale,r,v)
    SVector{7}([rs;vs;σ])
end
function unscale_state(ev,x)
    rs = x[SA[1,2,3]]
    vs = x[SA[4,5,6]]
    r,v = unscale_rv(ev.scale,rs,vs)
    return r,v,x[7]
end
function tt()


        ev = CPEGWorkspace()

        Rm = ev.params.gravity.R
        r0 = SA[Rm+125e3, 0.0, 0.0] #Atmospheric interface at 125 km altitude
        V0 = 5.845e3 #Mars-relative velocity at interface of 5.845 km/sec
        γ0 = -15.474*(pi/180.0) #Flight path angle at interface
        v0 = V0*SA[sin(γ0), cos(γ0), 0.0]
        σ0 = deg2rad(3)

        Rf = Rm+10.0e3 #Parachute opens at 10 km altitude
        rf = Rf*SA[cos(7.869e3/Rf)*cos(631.979e3/Rf); cos(7.869e3/Rf)*sin(631.979e3/Rf); sin(7.869e3/Rf)]


        # vehicle parameters
        ev.params.aero.Cl = 0.29740410453983374
        ev.params.aero.Cd = 1.5284942035954776
        ev.params.aero.A = 15.904312808798327    # m²
        ev.params.aero.m = 2400.0                # kg

        # qp solver settings
        ev.qp_solver_opts.verbose = false
        ev.qp_solver_opts.atol = 1e-10
        ev.qp_solver_opts.max_iters = 50

        # MPC stuff
        ev.cost.σ_tr = deg2rad(10) # radians
        ev.cost.rf = rf # goal position (meters, MCMF)
        ev.cost.γ = 1e3

        # sim stuff
        ev.dt = 2.0 # seconds

        # CPEG settings
        ev.verbose = true
        ev.ndu_tol = 1e-3
        ev.max_cpeg_iter = 30
        ev.miss_distance_tol = 1e1  # m

        ev.scale.uscale = 1e1

        N = 100
        X = [@SVector zeros(7) for i = 1:N]
        U = [@SVector zeros(1) for i = 1:N-1]
        X[1] = scale_state(ev,r0,v0,σ0)
        # main_cpeg(ev,r0,v0,σ0)
        # main_cpeg(ev,unscale_state(ev,X[1])...)
        # @show unscale_state(ev,X[1])
        # @show r0
        # @show v0
        # @show σ0
        # for i = 1:N-1
            x0s = [3.394493027569473, 0.399591678638693, 0.001009766098441132, -1.456578872683738, 4.90997118202939, 0.11856550572023335, 1.002017922478615]
            ev.U = SVector{1}.([[-1.0446874855627493], [-1.5480207770930021], [-2.0792718023435577], [-2.6165103698995678], [-3.133332980329963], [-3.602272358062922], [-3.998882259324137], [-4.305331038330731], [-4.512508731782757], [-4.6202527616273725], [-4.635962501176615], [-4.5722798946790375], [-4.444568523571758], [-4.268731622216839], [-4.059636636270206], [-3.8301827417699306], [-3.590907860816463], [-3.349978236653403], [-3.1134079271625916], [-2.8853876478909783],
             [-2.6686409267757703], [-2.4647591748467277], [-2.274492048866316], [-2.0979855378265144], [-1.9349693968414852], [-1.7849000190375819], [-1.6470663824705791], [-1.5206666365829429], [-1.4048620459126708], [-1.2988138882311988], [-1.2017077791152193], [-1.1127688912168858], [-1.0312706985091609], [-0.9565392051243061], [-0.8879540970721507], [-0.8249478583856876], [-0.7670035961772956], [-0.7136520995286785], [-0.6644684966304356], [-0.6190687583989716], [-0.5771062134840089],
              [-0.5382681804201349], [-0.5022727811224021], [-0.4688659711138508], [-0.4378188022176291], [-0.40892492032561795], [-0.38199829236131694], [-0.356871151306608], [-0.3333921451660487], [-0.3114246742864395], [-0.29084540102664813], [-0.2715429160242747], [-0.2534165459799516], [-0.23637528880207825], [-0.22033686300468008], [-0.2052268593482769], [-0.1909779838052427], [-0.1775293819838559], [-0.16482603613939062], [-0.1528182268259207], [-0.141461052095098],
              [-0.13071399792837857], [-0.12054055430002916], [-0.11090787191443939], [-0.10178645524833921], [-0.09314988806239005], [-0.08497458803317286], [-0.07723958760150706], [-0.0699263385415113], [-0.06301853813160513], [-0.05650197515788373], [-0.050364394305557204], [-0.0445953777984172], [-0.0391862434320541], [-0.034129958415797114], [-0.029421068692721798], [-0.025055643647920903], [-0.021031236343760105], [-0.017346859638204308], [-0.014002978749707781],
               [-0.011001521030997203], [-0.008345903905917261], [-0.006041082110155127], [-0.004093615560086733], [-0.002511759356214488], [-0.0013055776105475364], [-0.0004870829722087122], [-7.040391306118486e-5]])
            # call cpeg
            # @show X[i]
            # @show ev.U
            # main_cpeg(ev, unscale_state(ev,X[i])...)
            # for i = 1:length(ev.U)
            #     ev.U[i] = SA[0.0]
            # end
            ev.U = [ev.U[i+1] for i = 1:length(ev.U)-1]
            r0, v0, σ0 = unscale_state(ev,x0s)
            althist, drhist, crhist, σhist = main_cpeg(ev, r0,v0,σ0)
            xf_dr, xf_cr = rangedistances(ev,rf,SVector{6}([r0;v0]))


            num2plot = float(length(althist))
            plot_groundtracks(drhist/1e3,crhist/1e3,althist/1e3,σhist,xf_dr/1e3,xf_cr/1e3,num2plot,"quad")
            # U[i] = ev.U[1]

            # dr = 100/ev.scale.dscale * @SVector randn(3)
            # X[i]
            # dv = @SVector randn(3)
            # @show X[i]
            # X[i+1] = rk4(ev,X[i],ev.U[1],ev.dt/ev.scale.tscale) #+ SA[dr[1];dr[2];dr[3];0;0;0;0]

        # end

        # alt, dr, cr = postprocess(ev,X,[r0;v0])

        # # althist, drhist, crhist = main_cpeg(ev,x0_s)
        # # main_cpeg(ev,x0_s)
        # main_cpeg(ev,r0,v0,σ0)
        # r02 = r0 - SA[100,600,4.2e3]
        # main_cpeg(ev,r02,v0,σ0)
        #
        # pp2 = ev.σ
        # pp = [ev.U[i][1] for i = 1:length(ev.U)]
        # #
        # mat"
        # figure
        # hold on
        # title('Control')
        # %plot(rad2deg($pp))
        # plot($pp)
        # hold off"
        #
        # mat"
        # figure
        # title('Bank angle deg')
        # hold on
        # plot(rad2deg($pp2))
        # hold off"


        return nothing
end

function plot_groundtracks(drhist,crhist,althist,σhist,xf_dr,xf_cr,num2plot,id)

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

    mat"
    figure
    hold on
    rgb1 = [29 38 113]/255;
    rgb2 = 1.3*[195 55 100]/255;
    drgb = rgb2-rgb1;
    cmap = [];
    for i = 1:round($num2plot)
        bank = $σhist{i};
        plot(rad2deg(bank),'Color',rgb1 + drgb*(i)/($num2plot),'linewidth',3)
        cmap = [cmap;rgb1 + drgb*(i)/($num2plot)];
    end
    colormap(cmap);
    pp = colorbar;
    pp.Ticks = 0:(1/$num2plot):1;
    pp.Location = 'northoutside';
    pp.TickLabels = 0:round($num2plot);
    %xlim([400,700])
    %ylim([8,30])
    xlabel('step')
    ylabel('bank angle (degrees)')
    hold off
    %saveas(gcf,'alt.png')
    fleg = legend('figure()');
    set(fleg,'visible','off')
    "

    return nothing
end

tt()


# # incorrect
# typeof(ev.params.aero) = AeroParameters
# typeof(ρ) = Float64
# typeof(r) = SVector{3, ForwardDiff.Dual{ForwardDiff.Tag{var"#475#479"{Int64, CPEGWorkspace, Float64}, Float64}, Float64, 7}}
# typeof(v) = SVector{3, ForwardDiff.Dual{ForwardDiff.Tag{var"#475#479"{Int64, CPEGWorkspace, Float64}, Float64}, Float64, 7}}
#
# # correct
# typeof(ev.params.aero) = AeroParameters
# typeof(ρ) = ForwardDiff.Dual{ForwardDiff.Tag{var"#477#481"{Int64, CPEGWorkspace, Vector{SVector{7, Float64}}, Float64}, Float64}, Float64, 1}
# typeof(r) = SVector{3, ForwardDiff.Dual{ForwardDiff.Tag{var"#477#481"{Int64, CPEGWorkspace, Vector{SVector{7, Float64}}, Float64}, Float64}, Float64, 1}}
# typeof(v) = SVector{3, ForwardDiff.Dual{ForwardDiff.Tag{var"#477#481"{Int64, CPEGWorkspace, Vector{SVector{7, Float64}}, Float64}, Float64}, Float64, 1}}



# # succeeding
# typeof(r) = SVector{3, ForwardDiff.Dual{ForwardDiff.Tag{var"#1006#1007"{CPEGWorkspace}, Float64}, Float64, 7}}
# typeof(h) = ForwardDiff.Dual{ForwardDiff.Tag{var"#1006#1007"{CPEGWorkspace}, Float64}, Float64, 7}
# typeof(ρ) = ForwardDiff.Dual{ForwardDiff.Tag{var"#1006#1007"{CPEGWorkspace}, Float64}, Float64, 7}
# typeof(ev.params.aero) = AeroParameters
# typeof(r) = SVector{3, ForwardDiff.Dual{ForwardDiff.Tag{var"#1006#1007"{CPEGWorkspace}, Float64}, Float64, 7}}
# typeof(v) = SVector{3, ForwardDiff.Dual{ForwardDiff.Tag{var"#1006#1007"{CPEGWorkspace}, Float64}, Float64, 7}}
#
# # failing
# typeof(r) = SVector{3, ForwardDiff.Dual{ForwardDiff.Tag{var"#1006#1007"{CPEGWorkspace}, Float64}, Float64, 7}}
# typeof(h) = ForwardDiff.Dual{ForwardDiff.Tag{var"#1006#1007"{CPEGWorkspace}, Float64}, Float64, 7}
# typeof(ρ) = Float64
# typeof(ev.params.aero) = AeroParameters
# typeof(r) = SVector{3, ForwardDiff.Dual{ForwardDiff.Tag{var"#1006#1007"{CPEGWorkspace}, Float64}, Float64, 7}}
# typeof(v) = SVector{3, ForwardDiff.Dual{ForwardDiff.Tag{var"#1006#1007"{CPEGWorkspace}, Float64}, Float64, 7}}
