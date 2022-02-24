
struct Parameters
    density::DensityParameters
    gravity::GravityParameters
    aero::AeroParameters
    function Parameters()
        new(DensityParameters(),GravityParameters(),AeroParameters())
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

struct CPEGWorkspace
    params::Parameters
    scale::Scaling
    planet::Planet
    solver_opts::SolverSettings
    function CPEGWorkspace()
        new(Parameters(),Scaling(),Planet(),SolverSettings())
    end
end




# # density methods
# @inline function altitude(ev::CPEGWorkspace,r::SVector{3,T}) where T
#     norm(r) - ev.params.gravity.R
# end
#
# @inline function density(ev::CPEGWorkspace,r::SVector{3,T}) where T
#     return density(ev.params.density,altitude(ev,r))
# end
