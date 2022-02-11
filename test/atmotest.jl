


# import CPEG

p = CPEG.DensityParameters()


# this data came from a marsgramm output
hs = [5*1e3, 30e3, 50e3, 75e3, 100e3]
ρs = [8.277E-03, 7.555E-04, 7.081E-05, 5.867E-06, 3.612E-07]

for i = 1:length(hs)

    @test log(ρs[i]) ≈ log(CPEG.density(p,hs[i])) rtol = 1e-2

end
