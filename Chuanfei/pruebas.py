"""
No tiene nada útil pero está bueno para probar cosas en la terminal sin correr
todo de nuevo
"""


# MAVEN
x0 = 0.78
e = 0.9

R = [1.082, -0.064, 0.515]
theta = np.linspace(0, np.pi * 2, 100)

r0 = R - np.array([x0, 0, 0])
theta0 = np.arccos(r0[0] / np.linalg.norm(r0))

L0 = np.linalg.norm(r0) * (1 + e * np.cos(theta0))
r1 = L0 / (1 + e * np.cos(theta))

X1_M = x0 + r1 * np.cos(theta)
Y1_M = r1 * np.sin(theta)


def beta(x0, e):

  R = [1.082, -0.064, 0.515]
  theta = np.linspace(0, np.pi * 2, 100)

  r0 = R - np.array([x0, 0, 0])
  theta0 = np.arccos(r0[0] / np.linalg.norm(r0))

  L0 = np.linalg.norm(r0) * (1 + e * np.cos(theta0))
  r1 = L0 / (1 + e * np.cos(theta))

  X1_M = x0 + r1 * np.cos(theta)
  Y1_M = r1 * np.sin(theta)
  return(X1_M, Y1_M)


# betas
fig, ax = plt.subplots(1, 1)
xy = np.column_stack([x.flat, y.flat])  # Create a (N, 2) array of (x, y) pairs.
grid_x, grid_y = np.mgrid[x.min() : x.max() : 1000j, y.min() : y.max() : 1000j]
grid_0 = scipy.interpolate.griddata(
    xy, np.log(beta_str_z0), (grid_x, grid_y), method="linear"
)

for x0 in [0.5, 0.55, 0.6]:
  for ep in [0.95]:
    ax.plot(beta(x0, ep)[0], beta(x0, ep)[1], linestyle="--", label="beta=1")
    ax.plot(X1_M, Y1_M, c="k", linestyle="-", label="MAVEN")

ax.set_xlim([1.0, 2])
ax.set_ylim([-0.5, 0.5])

ax.pcolormesh(
    grid_x, grid_y, ma.masked_invalid(grid_0), cmap="coolwarm", vmin=-4, vmax=4
)
ax.set_title(r"Z=0 log($\beta*$)")
ax.set_aspect("equal", "box")
plt.show()

fig = plt.figure(1, (1., 2.))
grid = AxesGrid(fig, 111,  # similar to subplot(142)
                nrows_ncols=(1, 2),
                axes_pad=0.5,
                share_all=True,
                label_mode="L",
                cbar_location="right",
                cbar_mode="single",
                )

ax0 = grid[0]
ax1 = grid[1]
ax0.pcolormesh(
    grid_x, grid_y, ma.masked_invalid(grid_0), cmap="coolwarm", vmin=-4, vmax=4
)
ax0.set_title(r"Z=0 log($\beta*$)")
ax1.pcolormesh(
    grid_x, grid_y, ma.masked_invalid(grid_1), cmap="coolwarm", vmin=-4, vmax=4
)
ax1.set_title(r"Y=0 log($\beta*$)")
ax0.set_aspect("equal", "box")
ax1.set_aspect("equal", "box")
# cbar_ax = fig.add_axes([0.9, 0.1, 0.04, 0.7])  # [left, bottom, width, height]
ax0.set_ylabel(r"y (R$_M$)")
ax0.set_xlabel(r"x (R$_M$)")
ax1.set_xlabel(r"x (R$_M$)")
ax0.set_xlim([1.0, 2])
ax0.set_ylim([-0.5, 0.5])

ax1.set_xlim([1.0, 2])
ax1.set_ylim([-0.5, 0.5])
# grid.cbar_axes[0].colorbar(im)
fig.colorbar(cm.ScalarMappable(norm=Normalize(-4, 4),
             cmap="coolwarm"), cax=grid.cbar_axes[0])

figure.set_size_inches(9, 6)
# plt.savefig("../../../Dropbox/Paper2/beta_2d.png", dpi=600)
plt.show()


xy = np.column_stack([x.flat, y.flat])  # Create a (N, 2) array of (x, y) pairs.
# Interpolate and generate heatmap:
grid_x, grid_y = np.mgrid[x.min() : x.max() : 1000j, y.min() : y.max() : 1000j]

# interpolation method can be linear or cubic


def data(i, xy, z, grid_x, grid_y, zmin, zmax, metodo="linear", colormap="coolwarm"):
    grid_z = scipy.interpolate.griddata(
      xy,  B_y0[:, 0], (grid_x, grid_y), method=metodo)
    grid[i].pcolormesh(
        grid_x, grid_y, ma.masked_invalid(grid_z), cmap=colormap, vmin=zmin, vmax=zmax
    )


fig = plt.figure(1, (2., 3.))
grid = AxesGrid(fig, 111,  # similar to subplot(142)
                nrows_ncols=(2, 3),
                axes_pad=0.1,
                share_all=True,
                label_mode="L",
                cbar_location="right",
                cbar_mode="single",
                )

data(0, xy, , grid_x, grid_y, -50, 50)
data(1, xy, z, grid_x, grid_y, -50, 50)
data(2, xy, z, grid_x, grid_y, -50, 50)
data(3, xy, z, grid_x, grid_y, -50, 50)
data(4, xy, z, grid_x, grid_y, -50, 50)
data(5, xy, z, grid_x, grid_y, -50, 50)
for i in range(6):
    grid[i].plot(X1, Y1, c="k", linestyle="--", label="beta=1")
    grid[i].plot(X1_M, Y1_M, c="k", linestyle="-", label="MAVEN")
    grid[i].set_xlim([1.0, 2])
    grid[i].set_ylim([-0.5, 0.5])
    grid[i].set_title("titulo")
plt.show()

subplot_2d(x, y, B_z0[:, 2], -50, 50, r"Z=0 B$_z$", "coolwarm")
for i in [0, 1, 2]:
    plt.setp(ax[0, i].get_xticklabels(), visible=False)
    ax[1, i].set_xlabel(r"x (R$_M$)")
    ax[0, i].set_aspect("equal", "box")
    ax[1, i].set_aspect("equal", "box")
for i in [0, 1]:
    plt.setp(ax[i, 1].get_yticklabels(), visible=False)
    plt.setp(ax[i, 2].get_yticklabels(), visible=False)
ax[0, 0].set_ylabel(r"z (R$_M$)")
ax[1, 0].set_ylabel(r"y (R$_M$)")

cbar_ax = fig.add_axes([0.9, 0.1, 0.04, 0.85])  # [left, bottom, width, height]
cb = fig.colorbar(
    cm.ScalarMappable(norm=Normalize(-50, 50), cmap="coolwarm"), cax=cbar_ax
)
cb.ax.set_title("B (nT)")
figure = plt.gcf()  # get current figure
figure.set_size_inches(9, 6)
# when saving, specify the DPI
# plt.savefig("../../../Dropbox/Paper2/B_2d.png", dpi=600)
plt.show()
