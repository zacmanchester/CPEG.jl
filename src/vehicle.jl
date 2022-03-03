
struct Parameters
    density::CPEGDensityParameters
    gravity::GravityParameters
    aero::AeroParameters
    function Parameters()
        new(CPEGDensityParameters(),GravityParameters(),AeroParameters())
    end
end

struct Planet
    Ω::Float64
    ω::SVector{3,Float64}
    function Planet()
        # mars
        Ω = 7.088218127774194e-5 # rad/s
        ω = SA[0,0,Ω]
        return new(Ω,ω)
    end
end

mutable struct COST
    rf::SVector{3,Float64}
    γ::Float64
    σ_tr::Float64
    function COST()
        new(SA[1.,2.,3.], 1e3 , deg2rad(20))
    end
end

mutable struct CPEGWorkspace
    params::Parameters
    scale::Scaling
    planet::Planet
    qp_solver_opts::SolverSettings
    cost::COST
    U::Vector{SVector{1,Float64}}
    miss_distance::Float64
    dt::Float64
    σ::Vector{Float64}
    max_cpeg_iter::Int64
    miss_distance_tol::Float64
    ndu_tol::Float64
    verbose::Bool
    function CPEGWorkspace()
        U = [SA[0.0] for i = 1:1000]
        new(Parameters(), Scaling(), Planet(), SolverSettings(), COST(), U,0.0,0.0,zeros(2),20,1e3,1e-2,true)
    end
end

# function initialize_CPEG(p::)



# # density methods
# @inline function altitude(ev::CPEGWorkspace,r::SVector{3,T}) where T
#     norm(r) - ev.params.gravity.R
# end
#
# @inline function density(ev::CPEGWorkspace,r::SVector{3,T}) where T
#     return density(ev.params.density,altitude(ev,r))
# end
