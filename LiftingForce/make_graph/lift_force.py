import numpy as np
from matplotlib import pyplot as plt
from config.cfg import Data


data = Data()

def get_plate_indexes() -> tuple[int, int, int]:
    x_coords = data.get_coords("x")["x"]
    x_beg = np.where(x_coords > 0.0)[0][0]
    x_end = np.where(x_coords > data.cfg.real_plate_len)[0][0]
    y_coords = data.get_coords("y")["y"]
    y_end = np.where(y_coords > 0.0)[0][0]
    return x_beg, x_end, y_end - 1, y_end + 1

x_beg, x_end, y_beg, y_end = get_plate_indexes()
print(x_beg, x_end, y_beg, y_end)
v_coords = data.get_coords("v")["v"]

def calc_liftig_force(ind):
    distribution = np.sum(data.get_data(ind)[:, :, y_beg:y_end, x_beg:x_end], axis=data.get_axis("vx"))
    dp_up   = 0.0
    dp_down = 0.0
    print(distribution.shape)
    for i in range(distribution.shape[2]):
        for j, v in enumerate(v_coords):
            print(j, v, i)
            if (v > 0):
                dp_up += v * distribution[j, 1, i]
            else:
                dp_down += v * distribution[j, 0, i]
    return -2 * (dp_down + dp_up)

iter = np.arange(len(data)) * data.cfg.s_t
liftintg_force = np.zeros(len(data))
for i in range(len(data)):
    liftintg_force[i] = calc_liftig_force(i)

plt.plot(iter, liftintg_force)
plt.show()
