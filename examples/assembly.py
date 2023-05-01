from firedrake import *
from firedrake_adjoint import *

m = UnitSquareMesh(2, 2)
fs = FunctionSpace(m, "CG", 1)

f = Function(fs)
g = Function(fs)
x, y = SpatialCoordinate(m)
g.assign(1.0)
f.interpolate(sin(2*pi*x) * cos(2*pi*y) * g)
g.assign(2.0)

J = assemble(inner(f-g, f-g)*dx) + assemble(inner(f, f)*dx)

Jhat = ReducedFunctional(J, Control(f))
print(Jhat.derivative().dat.data[:])

get_working_tape().visualise("assembly.pdf")
