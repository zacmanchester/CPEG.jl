
# ev = CPEG.CPEGWorkspace()
#
# Rm = ev.params.gravity.R
# r0 = SA[Rm+125e3, 0.0, 0.0] #Atmospheric interface at 125 km altitude
# V0 = 5.845e3 #Mars-relative velocity at interface of 5.845 km/sec
# γ0 = -15.474*(pi/180.0) #Flight path angle at interface
# v0 = V0*SA[sin(γ0), cos(γ0), 0.0]
# σ0 = deg2rad(90)
#
# r0sc,v0sc = CPEG.scale_rv(ev.scale,r0,v0)
#
# dt = 1.0/ev.scale.tscale
#
# x0 = SA[r0sc[1],r0sc[2],r0sc[3],v0sc[1],v0sc[2],v0sc[3],σ0]
#
# N = 100
# U = [SA[0.0] for i = 1:N-1]
# #
# X1,U1 = CPEG.rollout(ev, x0, U, dt)
#
# @test length(X1) == 113
# @test (norm((X1[end])[1:3]) * ev.scale.dscale - ev.params.gravity.R) / 1000.0 < 10.0
# @test (norm((X1[end-1])[1:3]) * ev.scale.dscale - ev.params.gravity.R) / 1000.0 > 10.0
#
#
# # now update the scaling
#
# ev.scale.dscale = 1.0
# ev.scale.tscale = 1.0
# ev.scale.uscale = 1.0
#
# x0 = SA[r0[1],r0[2],r0[3],v0[1],v0[2],v0[3],σ0]
# dt = 1.0/ev.scale.tscale
# X,U = CPEG.rollout(ev, x0, U, dt)
#
# @test length(X) == 113
# @test (norm((X[end])[1:3]) * ev.scale.dscale - ev.params.gravity.R) / 1000.0 < 10.0
# @test (norm((X[end-1])[1:3]) * ev.scale.dscale - ev.params.gravity.R) / 1000.0 > 10.0


# # test scaling
# for i = 1:length(X)
#
#     r1 = SVector{3}(X1[i][1:3])


function test_rollout(ev,r0,v0,σ0)

    r0sc,v0sc = CPEG.scale_rv(ev.scale,r0,v0)

    # dt = 1.0/ev.scale.tscale
    ev.dt = 1.0

    x0 = SA[r0sc[1],r0sc[2],r0sc[3],v0sc[1],v0sc[2],v0sc[3],σ0]

    N = 100
    U = [SA[0.0] for i = 1:N-1]
    #
    X,U = CPEG.rollout(ev, x0, U)

    # @show length(X)

    @test length(X) == 241
    @test (norm((X[end])[1:3]) * ev.scale.dscale - ev.params.gravity.R)  < 10.0e3
    @test (norm((X[end-1])[1:3]) * ev.scale.dscale - ev.params.gravity.R) > 10.0e3

    return X
end




ev1 = CPEG.CPEGWorkspace()

Rm = ev1.params.gravity.R
r0 = SA[Rm+125e3, 0.0, 0.0] #Atmospheric interface at 125 km altitude
V0 = 5.845e3 #Mars-relative velocity at interface of 5.845 km/sec
γ0 = -15.474*(pi/180.0) #Flight path angle at interface
v0 = V0*SA[sin(γ0), cos(γ0), 0.0]
σ0 = deg2rad(42)

X1 = test_rollout(ev1,r0,v0,σ0)

ev2 = CPEG.CPEGWorkspace()

ev2.scale.dscale = 1.0
ev2.scale.tscale = 1.0
ev2.scale.uscale = 1.0

X2 = test_rollout(ev2,r0,v0,σ0)

@test length(X1) == length(X2)
for i = 1:length(X1)

    r1_sc = SVector{3}(X1[i][1:3])
    v1_sc = SVector{3}(X1[i][4:6])

    r1,v1 = CPEG.unscale_rv(ev1.scale,r1_sc,v1_sc)

    r2_sc = SVector{3}(X2[i][1:3])
    v2_sc = SVector{3}(X2[i][4:6])

    r2,v2 = CPEG.unscale_rv(ev2.scale,r2_sc,v2_sc)

    @test r1 ≈ r2 rtol = 1e-6
    @test v1 ≈ v2 rtol = 1e-6
end
