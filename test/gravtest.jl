
cp = CPEG.CPEGWorkspace()
g = cp.params.gravity
# check gravity at surface
r = g.R*normalize(SA[1, 2, 4.0])


# check gravity at surface is around 3.7 m/s^2
a = CPEG.gravity(g, r)
@test norm(a) < 3.8
@test norm(a) > 3.6

# check J2 calculation by doing both parts seperately

# position vector
r_1 = g.R*1.1*normalize(SA[1, -2, 1.8])

# accel from the gravity function
a1 = CPEG.gravity(g, r_1)

# gravity due to central body
a1_central = -r_1*g.μ/(norm(r_1)^3)

# gravity due to J2
nr = norm(r_1)
x,y,z = r_1
a1_j2 = -1.5*g.J2*(g.μ/nr^2)*(g.R/nr)^2*SA[(1-5*(z/nr)^2)*x/nr,
                                           (1-5*(z/nr)^2)*y/nr,
                                           (3-5*(z/nr)^2)*z/nr]

# test the sum of these is equal to the output of the gravity function
@test a1 ≈ (a1_central + a1_j2) rtol = 1e-12
