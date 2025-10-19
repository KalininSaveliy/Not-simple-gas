import json
import os
from glob import glob


def path_list(folder:str, pattern:str)->list:
    return glob(os.path.join(folder, pattern))

class Config:
    def __init__(self, filename: str):
        self.filename = filename
        with open(filename, "r") as file:
            self.cfg = json.load(file)

    def make_cpp_config(self):
        data = [
            self.cfg["n_t"],
            self.cfg["s_t"],
            self.cfg["n_x"],
            self.cfg["n_y"],
            self.cfg["n_v"],
            self.cfg["v_cut"],
            self.cfg["alpha"],
            self.cfg["beta"],
            self.cfg["Knudsen"],
            self.cfg["Mach"],
            self.cfg["T1"],
            self.cfg["T2"],
            self.cfg["real_plate_len"],
            int(self.cfg["debug"]),
            self.cfg["save_folder"]
        ]
        c_cfg_name = self.filename.split(".")[0] + "_Cpp.txt"
        with open(c_cfg_name, "w") as file:
            for el in data:
                file.write(str(el) + '\n')


if __name__ == "__main__":
    cfg_list = path_list(os.path.dirname(__file__), "*.json")
    for cfg_name in cfg_list:
        cfg = Config(cfg_name)
        cfg.make_cpp_config()
