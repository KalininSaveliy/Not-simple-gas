import xarray as xr
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from config.cfg import Data


def read(data: Data, i:int)->xr.DataArray:
        ar = data.get_xy_data(i)
        data = xr.DataArray(ar, dims=("y", "x"), coords=data.get_coords("x", "y"))
        return data


cfg_name = "main"  # TODO: make it not hard const
data = Data("config/" + cfg_name + ".json")

# find min and max for cbar
data_min = +10e10
data_max = -10e10
for i in range(len(data)):
    d = data.get_xy_data(i).reshape(-1)
    data_min = min(data_min, d.min())
    data_max = max(data_max, d.max())


plate = [[0, data.cfg["real_plate_len"]], [0, 0]]

cbar = None
fig, ax = plt.subplots()
ax.plot(*plate, color="black", linewidth=2)

def update(frame_ind):
    print("frame ", frame_ind)
    x = read(data, frame_ind)
    mesh = x.plot.pcolormesh(
        ax=ax,
        x='x', y='y',
        shading='auto',
        add_colorbar=False,
        cmap='viridis',
        vmin=data_min,
        vmax=data_max
    )
    # color bar (should be one for all pictures)
    if frame_ind == 0:
        cbar = fig.colorbar(mesh, ax=ax)
        cbar.set_label("n / n0")
    it = frame_ind * data.cfg["s_t"]
    time = round(it * data.cfg["real_plate_len"] * data.cfg["Knudsen"] / 290 * 1000, 1)  # miliSeconds
    temperature = "T1/T0=" +  str(data.cfg["T1"]) + ", T2/T0=" + str(data.cfg["T2"])
    ax.set_title(str(it) + " (" + str(time) + "ms, " + temperature + ")")


ani = FuncAnimation(fig=fig, func=update, frames=len(data), interval=1000, repeat=False, init_func=lambda : None)
ani.save("graphs/eval.gif", writer='pillow', fps=1)
