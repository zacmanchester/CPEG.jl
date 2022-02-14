

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




# density methods
@inline function altitude(ev::EntryVehicle,r::SVector{3,T}) where T
    norm(r) - ev.params.gravity.R
end

@inline function density(ev::EntryVehicle,r::SVector{3,T}) where T
    return density(ev.params.density,altitude(ev,r))
end
