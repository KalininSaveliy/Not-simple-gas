import pandas as pd
import xarray as xr
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from config.cfg import Config, path_list


class Data:
    def __init__(self, cfg:Config):
        self.sep = ','
        self.header = None
        self.data_dir = "data\\" + cfg.cfg["save_folder"] + "\\"
        self.search_data()
        self.read_coords()
        self.find_min_max_value()

    def search_data(self)->None:
        self.data_list = [self.data_dir + str(n) + ".csv" for n in range(0, cfg.cfg["n_t"], cfg.cfg["s_t"])]
        self.length = len(self.data_list)
        print(self.data_list)
        self.grid_file = self.data_dir + "grid.csv"

    def read_coords(self)->None:
        self.coords = dict()
        df = pd.read_csv(self.grid_file)
        for dim in ["x", "y"]:
            dim_coords = df[dim]
            dim_coords = dim_coords[dim_coords.notna()]
            self.coords[dim] = dim_coords
    
    def read_table(self, i:int)->pd.DataFrame:
        return pd.read_csv(self.data_list[i], sep = self.sep, header=self.header)

    def find_min_max_value(self)->None:
        self.data_min = +10e10  # TODO: rename variables
        self.data_max = -10e10
        for i in range(len(self.data_list)):
            data = self.read_table(i).to_numpy().reshape(-1)
            self.data_min = min(self.data_min, data.min())
            self.data_max = max(self.data_max, data.max())

    def read(self, i:int)->xr.DataArray:
        df = self.read_table(i)
        data = xr.DataArray(df, dims=("y", "x"), coords=self.coords)
        return data


cfg_name = "main"  # TODO: make it not hard const

cfg = Config("config\\" + cfg_name + ".json")
plate = [[0, cfg.cfg["real_plate_len"]], [0, 0]]
data = Data(cfg)

cbar = None
fig, ax = plt.subplots()
ax.plot(*plate, color="black", linewidth=2)

def update(frame_ind):
    print("frame ", frame_ind)
    x = data.read(frame_ind)
    mesh = x.plot.pcolormesh(
        ax=ax,
        x='x', y='y',
        shading='auto',
        add_colorbar=False,
        cmap='viridis',
        vmin=data.data_min,
        vmax=data.data_max
    )
    # color bar (should be one for all pictures)
    if frame_ind == 0:
        cbar = fig.colorbar(mesh, ax=ax)
        cbar.set_label("n / n0")
    it = frame_ind * cfg.cfg["s_t"]
    time = round(it * cfg.cfg["real_plate_len"] * cfg.cfg["Knudsen"] / 290 * 1000, 1)  # miliSeconds
    temperature = "T1/T0=" +  str(cfg.cfg["T1"]) + ", T2/T0=" + str(cfg.cfg["T2"])
    ax.set_title(str(it) + " (" + str(time) + "ms, " + temperature + ")")

ani = FuncAnimation(fig=fig, func=update, frames=data.length, interval=1000, repeat=False, init_func=lambda : None)
ani.save("graphs/eval.gif", writer='pillow', fps=1)

# plt.show()
