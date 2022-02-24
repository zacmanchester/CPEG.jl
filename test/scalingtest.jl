

ev = CPEG.CPEGWorkspace()

ev.scale.dscale = 124.92
ev.scale.tscale = 0.324
ev.scale.uscale = 1.45

# scaled r,v,a
rs = SA[23.45,434.35,-53.23]
vs = SA[-.043,34.2,-34.22]
as = SA[-2.3,4.344,5.32]

# unscale rs vs
r,v = CPEG.unscale_rv(ev.scale,rs,vs)

# rescale and compare
rs2,vs2 = CPEG.scale_rv(ev.scale,r,v)
@test rs ≈ rs2 atol = 1e-12
@test vs ≈ vs2 atol = 1e-12

# unscale vs as
v2,a2 = CPEG.unscale_va(ev.scale,vs,as)

# compare to the first one
@test v2 ≈ v rtol = 1e-6

v2_2s, a2_2s = CPEG.scale_va(ev.scale,v2,a2)

@test v2_2s ≈ vs atol = 1e-12
@test a2_2s ≈ as atol = 1e-12
