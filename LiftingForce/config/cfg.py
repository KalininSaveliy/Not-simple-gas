import json
from collections.abc import Mapping
from types import MappingProxyType
import os
from glob import glob


def path_list(folder:str, pattern:str)->list:
    return glob(os.path.join(folder, pattern))

class Config(Mapping):
    def __init__(self, filename: str):
        self._filename = filename
        with open(filename, "r") as file:
            self._data = MappingProxyType(json.load(file))

    def __getitem__(self, key):
        return self._data[key]

    def __getattr__(self, name):
            return self._data[name]

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    # # ---- запрет изменения ----
    # def __setattr__(self, name, value):
    #     if name.startswith("_"):
    #         super().__setattr__(name, value)
    #     else:
    #         raise TypeError("Config is read-only")

    # def __setitem__(self, key, value):
    #     raise TypeError("Config is read-only")

    # def __delattr__(self, name):
    #     raise TypeError("Config is read-only")

    # def __delitem__(self, key):
    #     raise TypeError("Config is read-only")

    def make_cpp_config(self):
        data = [
            self.n_t,
            self.s_t,
            self.n_x,
            self.n_y,
            self.n_v,
            self.v_cut,
            self.alpha,
            self.beta,
            self.Knudsen,
            self.Mach,
            self.T1,
            self.T2,
            self.real_plate_len,
            int(self.debug),
            self.save_folder
        ]
        c_cfg_filename = self._filename.split(".")[0] + "_Cpp.txt"
        with open(c_cfg_filename, "w") as file:
            for el in data:
                file.write(str(el) + '\n')


if __name__ == "__main__":
    cfg_list = path_list(os.path.dirname(__file__), "*.json")
    for cfg_name in cfg_list:
        cfg = Config(cfg_name)
        cfg.make_cpp_config()
