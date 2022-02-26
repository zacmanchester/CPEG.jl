
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
    solver_opts::SolverSettings
    cost::COST
    U::Vector{SVector{1,Float64}}
    miss_distance::Float64
    function CPEGWorkspace()
        U = [SA[1.0], SA[1.0]]
        new(Parameters(), Scaling(), Planet(), SolverSettings(), COST(), U,0.0)
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
