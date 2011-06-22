from ffc import *

# Reserved variables for forms
(a, L, M) = (None, None, None)

# Reserved variable for element
element = None

# Copyright (c) 2005 Johan Jansson (johanjan@math.chalmers.se)
# Licensed under the GNU GPL Version 2
#
# First added:  2005
# Last changed: 2006-03-28
#
# The bilinear form for the incompressible Navier-Stokes equations
# Compile this form with FFC: ffc Elasticity.form.

cell = "triangle"

K1 = VectorElement("Lagrange", cell, 1)
# Dimension of domain
d = K1.cell_dimension()
K2 = FiniteElement("Lagrange", cell, 1)
K3 = FiniteElement("Discontinuous Lagrange", cell, 0)
K4 = VectorElement("Discontinuous Lagrange", cell, 0)
K5 = VectorElement("Discontinuous Lagrange", cell, 0, d * d)

K = K1

v = Function(K)
dtU = Function(K)
U = Function(K)

Um = Function(K4)
P = Function(K2)
nu = Function(K3)
d1 = Function(K3)
d2 = Function(K3)
f = Function(K)

def tomatrix(q):
    return [ [q[d * i + j] for i in range(d)] for j in range(d) ]

def ugradu(u, v):
    return [dot(u, grad(v[i])) for i in range(d)]

SD = mult(d1, dot(ugradu(Um, U), ugradu(Um, v))) + \
    + mult(d1, dot(grad(P), ugradu(Um, v))) + \
    mult(d2, dot(div(U), div(v)))

S = mult(P, Identity(d)) - mult(nu, grad(U))

M = (   -dot(dtU, v) +
	-dot(ugradu(Um, U), v) + \
	dot(S, grad(v)) + \
	dot(f, v) +
	SD
	 ) * dx

compile([a, L, M, element], "Drag2D", {'language': 'dolfin', 'blas': False, 'form_postfix': True, 'precision': '15', 'cpp optimize': False, 'split_implementation': True, 'quadrature_points': False, 'output_dir': '.', 'representation': 'quadrature', 'cache_dir': None, 'optimize': False})
