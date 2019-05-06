import numpy as np
import matplotlib.pyplot as plt
alpha = 5
beta = 3
N = 500
DIM = 2

np.random.seed(2)

# Generate random points on the unit circle by sampling uniform angles
theta = np.random.uniform(0, 2*np.pi, (N,1))
eps_noise = 0.2 * np.random.normal(size=[N,1])
circle = np.hstack([np.cos(theta), np.sin(theta)])

# Stretch and rotate circle to an ellipse with random linear tranformation
B = np.random.randint(-3, 3, (DIM, DIM))
noisy_ellipse = circle.dot(B) + eps_noise

mag = np.loadtxt('../../datos/MAG_1s/mvn_mag_l2_2016092ss1s_20160401_v01_r01.sts', skiprows=148) #en mi compu
posicion = np.zeros((len(mag[:,0]), 3))
for j in range(11,14):
    posicion[:,j-11] = mag[:, j]

orbita = posicion / 3390 #radios marcianos

# Extract x coords and y coords of the ellipse as column vectors
X = noisy_ellipse[:,0:1]
Y = noisy_ellipse[:,1:]

# Formulate and solve the least squares problem ||Ax - b ||^2
A = np.hstack([X**2, X * Y, Y**2, X, Y])
b = np.ones_like(X)
x = np.linalg.lstsq(A, b)[0].squeeze()

# Print the equation of the ellipse in standard form
print('The ellipse is given by {0:.3}x^2 + {1:.3}xy+{2:.3}y^2+{3:.3}x+{4:.3}y = 1'.format(x[0], x[1],x[2],x[3],x[4]))

# Plot the noisy data
plt.scatter(X, Y, label='Data Points')
