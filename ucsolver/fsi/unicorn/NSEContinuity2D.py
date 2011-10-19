from ffc import *

# Reserved variables for forms
(a, L, M) = (None, None, None)

# Reserved variable for element
element = None

# Copyright (c) 2005 Johan Hoffman 
# Licensed under the GNU GPL Version 2
#
# Modified by Anders Logg 2006
#
# First added:  2005
# Last changed: 2006-03-28
#
# The contiuity equation for the incompressible 
# Navier-Stokes equations using cG(1)cG(1)
#
# Compile this form with FFC: ffc NSEContinuity2D.form.

cell = "triangle"

scalar = FiniteElement("Lagrange", cell, 1)
# Dimension of domain
d = scalar.cell_dimension()
vector = VectorElement("Lagrange", cell, 1)
constant_scalar = FiniteElement("Discontinuous Lagrange", cell, 0)

q  = TestFunction(scalar)  # test basis function
P  = TrialFunction(scalar) # trial basis function
P0  = Function(scalar) # trial basis function
uc = Function(vector)      # linearized velocity

delta1 = Function(constant_scalar) # stabilization parameter
delta3 = Function(constant_scalar) # stabilization parameter
k = Function(constant_scalar)

rho = Function(scalar)
f = Function(vector)
fc = Function(vector)
#n = Function(vector)
n = FacetNormal(cell)
k = Function(constant_scalar)

i0 = Index() # index for tensor notation
i1 = Index() # index for tensor notation

def ugradu(u, v):
    return [dot(u, grad(v[i])) for i in range(d)]

eps = 0.0
alpha = 2*k

a = (eps*P*q + delta1*dot(grad(q), grad(P)) + alpha*dot(grad(q), grad(P)))*dx
L = (-q*div(uc))*dx + alpha*dot(grad(q), grad(P0))*dx - mult(delta1, dot( ugradu(uc, uc), grad(q)))*dx + \
    0.0*delta1*dot(mult(rho, f), grad(q))*dx

compile([a, L, M, element], "NSEContinuity2D", {'language': 'dolfin', 'blas': False, 'form_postfix': True, 'precision': '15', 'cpp optimize': False, 'split_implementation': True, 'quadrature_points': False, 'output_dir': '.', 'representation': 'quadrature', 'cache_dir': None, 'optimize': False})
