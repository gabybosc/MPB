    import numpy as np
from scipy.spatial import KDTree  # kd-tree for quick nearest-neighbor lookup

# This class provides an index into a set of k-dimensional points which can be
# used to rapidly look up the nearest neighbors of any point.

n = 3
k = 32

# fill the cube with random points
data = np.random.rand(10000,n)
kdtree = KDTree(data)

# pick a point (at the center of the cube)
point = 0.5 * np.ones((1,n))

# Coords of k-Nearest Neighbors
dists, idxs = kdtree.query(point, k)
idxs = idxs[0]
X = data[idxs,:]

# Calculate coefficients
C = (np.dot(point.T, np.ones((1,k)))).T  # central node
dX = X - C  # diffs from central node
G = np.dot(dX, dX.T)
F = np.multiply(G, G)
v = np.diag(G)
N = np.eye(k) - np.dot(np.linalg.pinv(G), G)  # pinv = Moore-Penrose inverse
a = np.dot(np.dot(N, np.linalg.pinv(np.dot(F,N))), v)  # these are the coeffs

#  Temperature distribution is  T = 25.4 * r^2
X2 = np.multiply(X, X)  # x^2
C2 = np.multiply(C, C)  # C²
T = 25.4 * n * np.mean(X2, 1).T  # T es 25.4 * x² en 3D
Tc = 25.4 * n * np.mean(C2, 1).T  # central node
dT = T - Tc  # diffs from central node

# Analytical gradient ==>  gradT = 2*25.4* x
g = 2 * 25.4 * point
print("g[]: %s" % (g))

# Estimated gradient
y = np.dot(dX.T, np.multiply(a, dT))
print("y[]: %s,   Relative Error = %.8f" % (y, np.linalg.norm(g-y)/np.linalg.norm(g)))

"""
Given a set of M vectors in n-dimensions (call them b_k), find a set of
coeffs (call them a_k) which yields the best estimate of the identity
matrix and the zero vector

                                 M
 (1) min ||E - I||,  where  E = sum  a_k b_k b_k
     a_k                        k=1

                                 M
 (2) min ||z - 0||,  where  z = sum  a_k b_k
     a_k                        k=1


Note that the basis vectors {b_k} are not required
to be normalized, orthogonal, or even linearly independent.

First, define the following quantities:

  B             ==> matrix whose columns are the b_k
  G = B'.B      ==> transpose of B times B
  F = G @ G     ==> @ represents the hadamard product (multiplica coef a coef)
  v = diag(G)   ==> vector composed of diag elements of G

The above minimizations are equivalent to this linearly constrained problem

  Solve  F.a = v
  s.t.   G.a = 0

Let {X} denote the Moore-Penrose inverse of X.
Moore-Penrose inverse = (A'A)⁻¹A'
Then the solution of the linear problem can be written:

  N = I - {G}.G       ==> projector into nullspace of G
  a = N . {F.N} . v

The utility of these coeffs is that they allow you to write
very simple expressions for the derivatives of a tensor field.

Let D be the nabla operator and d be the difference operator with respect to
the central (aka 0th) node so that, for any scalar/vector/tensor quantity Y,
we have:
  dY = Y - Y_0

Let x_k be the position vector of the kth node.
And for our basis vectors, take
  b_k = dx_k  =  x_k - x_0.

Assume that each node has a field value associated with it
(e.g. temperature), and assume a quadratic model [about x = x_0]
for the field [g=gradient, H=hessian, ":" is the double-dot product]
Hessian is a square matrix of second-order partial derivatives of a
scalar-valued function, or scalar field
A double dot product is a double contraction over the last two indices of the
first tensor and the first two indices of the second tensor.

    Y = Y_0 + (x-x_0).g + (x-x_0)(x-x_0):H/2
    dY = dx.g + dxdx:H/2
    D2Y = I:H            ==> Laplacian of Y

Evaluate the model at the kth node

    dY_k = dx_k.g  +  dx_k dx_k:H/2

Multiply by a_k and sum

     M               M                  M
    sum a_k dY_k =  sum a_k dx_k.g  +  sum a_k dx_k dx_k:H/2
    k=1             k=1                k=1

                 =  0.g   +  I:H/2
                 =  D2Y / 2

Thus, we have a second order estimate of the Laplacian

                M
   Lap(Y_0) =  sum  2 a_k dY_k
               k=1


Now play the same game with a linear model
    dY_k = dx_k.g

But this time multiply by (a_k dx_k) and sum

     M                    M
    sum a_k dx_k dY_k =  sum a_k dx_k dx_k.g
    k=1                  k=1

                      =  I.g
                      =  g


In general, the derivatives at the central node can be estimated as

           M
    D#Y = sum  a_k dx_k#dY_k
          k=1

           M
    D2Y = sum  2 a_k dY_k
          k=1

where
   # stands for the {dot, cross, or tensor} product
       yielding the {div, curl,  or grad} of Y
 and
   D2Y stands for the Laplacian of Y
   D2Y = D.DY = Lap(Y)
"""
