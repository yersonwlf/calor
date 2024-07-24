#  geometría
from fenics import *
mesh = Mesh("copa.xml")
plot(mesh)

deltat = 0.1 
total_time = 25.0
k = 250  # coeficiente de difusividad térmica

#  espacio funcional
V = FunctionSpace(mesh, "P", 1)

# Leer el archivo para las superficies
caras = MeshFunction('size_t', mesh, "copa_facet_region.xml")

# condiciones de borde
u_entrada = 9
u_salida = 18

bc1 = DirichletBC(V, Constant(u_entrada), caras, 69)
bc2 = DirichletBC(V, Constant(u_salida), caras, 70)
bc = [bc1, bc2]

#  formulación variacional
u = TrialFunction(V)
v = TestFunction(V)
un = Function(V)

Copa = u*v*dx +  deltat*k*dot(grad(u), grad(v))*dx  - un*v*dx

a = lhs(Copa)
l = rhs(Copa)

vtkfile = File('solucion_copa/copa.pvd')
u = Function(V)
t = 0

num_iter = int(total_time/deltat)
for n in range(num_iter):
    solve(a == l, u, bc)
    un.assign(u)
    vtkfile << (u, t)
    t += deltat

