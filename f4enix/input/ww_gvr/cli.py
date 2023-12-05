from copy import deepcopy
from enum import Enum
import os
from pathlib import Path
from typing import Dict, Optional, Tuple

from f4enix.input.ww_gvr import WW


MAIN_MENU = """
 ***********************************************
        Weight window manipulator and GVR
 ***********************************************

 * Open weight window file   (open)   
 * Display ww information    (info)   
 * Write ww                  (write)
 * Export as VTK             (vtk)
 * Plot                      (plot)    
 * Weight window operation   (operate)
 * GVR generation            (gvr)
 * Exit                      (end)    
"""

OPERATE_MENU = """             
 * Softening and normalize   (soft)
 * Add                       (add)
 * Remove                    (rem)
 * Mitigate long histories   (mit)
 * Exit                      (end)
"""


class Command(Enum):
    OPEN = "open"
    INFO = "info"
    WRITE = "write"
    VTK = "vtk"
    PLOT = "plot"
    OPERATE = "operate"
    GVR = "gvr"
    END = "end"
    EMPTY = "empty"


class CommandOperate(Enum):
    SOFT = "soft"
    ADD = "add"
    REMOVE = "remove"
    MITIGATE = "mitigate"
    END = "end"
    EMPTY = "empty"


class Menu:
    def __init__(self) -> None:
        self.weight_windows: Dict[str, WW] = dict()
        self.extra_text: str = ""
        self.command_map = {
            Command.OPEN: self.handle_open,
            Command.INFO: self.handle_info,
            Command.WRITE: self.handle_write,
            Command.VTK: self.handle_vtk,
            Command.PLOT: self.handle_plot,
            Command.OPERATE: self.handle_operate,
            Command.GVR: self.handle_gvr,
            Command.END: exit,
        }
        self.command_operate_map = {
            CommandOperate.SOFT: self.handle_soft,
            CommandOperate.ADD: self.handle_add,
            CommandOperate.REMOVE: self.handle_remove,
            CommandOperate.MITIGATE: self.handle_mitigate,
            CommandOperate.END: self.display_menu,
        }
        self.main_menu_loop()

    def main_menu_loop(self):
        while True:
            self.display_menu()
            self.get_command()

    def display_menu(self):
        self.clear_screen()
        print(MAIN_MENU)
        print(self.extra_text)

    def display_operate_menu(self):
        self.clear_screen()
        print(OPERATE_MENU)
        print(self.extra_text)

    def get_command(self):
        while True:
            attempt = input(" Enter command: ").lower()
            try:
                function = self.command_map[Command(attempt)]
                function()
                break
            except ValueError:
                print(" Invalid command...")

    def handle_open(self):
        file_path = Path(input(" Enter file path: "))

        if file_path.stem in self.weight_windows:
            print(" Weight window already loaded...")
            self.handle_open()

        try:
            weight_window = WW.load_from_ww_file(file_path)
            self.weight_windows[file_path.stem] = weight_window
            self.extra_text = f" {file_path.stem} loaded..."
        except FileNotFoundError:
            self.extra_text = " File not found..."

    def handle_info(self):
        ww_key = self.select_ww_key()
        if ww_key is None:
            return
        ww = self.weight_windows[ww_key]
        self.extra_text = ww.info

    def handle_write(self):
        ww_key = self.select_ww_key()
        if ww_key is None:
            return
        ww = self.weight_windows[ww_key]
        ww.write_to_ww_file()
        self.extra_text = f" {ww_key}_written saved!"

    def handle_vtk(self):
        ww_key = self.select_ww_key()
        if ww_key is None:
            return
        ww = self.weight_windows[ww_key]
        ww.export_as_vtk()
        self.extra_text = f" {ww_key}.vtk saved!"

    def handle_plot(self):
        ww_key = self.select_ww_key()
        if ww_key is None:
            return
        ww = self.weight_windows[ww_key]
        ww.geometry.plot()

    def handle_gvr(self):
        file_path = Path(input(" Enter Meshtally file path: "))

        if file_path.stem in self.weight_windows:
            print(" Weight window already loaded with that file name...")
            self.handle_gvr()

        maximum_splitting, softening = self._ask_gvr_parameters()

        try:
            weight_window = WW.create_gvr_from_meshtally_file(
                file_path, maximum_splitting
            )
            weight_window.soften(softening)
            self.weight_windows[file_path.stem] = weight_window
            self.extra_text = f" {file_path.stem} loaded..."
        except FileNotFoundError:
            self.extra_text = " File not found..."

    def select_ww_key(self) -> Optional[str]:
        if len(self.weight_windows) == 0:
            self.extra_text = " No weight windows loaded..."
            return None

        if len(self.weight_windows) == 1:
            return list(self.weight_windows.keys())[0]

        for i, ww_key in enumerate(self.weight_windows.keys()):
            print(f" {i+1}. {ww_key}")
        while True:
            ww_key = input(" Select weight window by number: ")
            if ww_key.isdigit():
                ww_key = int(ww_key)
                if ww_key in range(1, len(self.weight_windows) + 1):
                    return list(self.weight_windows.keys())[ww_key - 1]
            print(" Invalid weight window number...")

    def _ask_gvr_parameters(self) -> Tuple[float, float]:
        try:
            maximum_splitting = float(
                input(" Enter maximum splitting ratio (default 5): ")
            )
        except ValueError:
            maximum_splitting = 5
        try:
            softening = float(input(" Enter softening factor (default 1): "))
        except ValueError:
            softening = 1.0

        return maximum_splitting, softening

    def handle_operate(self):
        self.extra_text = ""
        self.display_operate_menu()
        self.get_operate_command()

    def get_operate_command(self):
        while True:
            attempt = input(" Enter command: ").lower()
            try:
                function = self.command_operate_map[CommandOperate(attempt)]
                function()
                break
            except ValueError:
                print(" Invalid command...")

    def handle_soft(self):
        ww_key = self.select_ww_key()
        if ww_key is None:
            return
        ww = deepcopy(self.weight_windows[ww_key])
        softening = float(input(" Enter softening factor: "))
        ww.soften(softening)
        ww.file_path = ww.file_path.with_name(f"{ww_key}_softened")
        self.weight_windows[ww.file_path.stem] = ww
        self.extra_text = f" {ww.file_path.stem} created..."

    def handle_add(self):
        ww_key = self.select_ww_key()
        if ww_key is None:
            return
        ww = deepcopy(self.weight_windows[ww_key])
        normalization = float(input(" Enter normalization factor: "))
        softening = float(input(" Enter softening factor: "))
        ww.add_particle(norm=normalization, soft=softening)
        ww.file_path = ww.file_path.with_name(f"{ww_key}_2_particles")
        self.weight_windows[ww.file_path.stem] = ww
        self.extra_text = f" {ww.file_path.stem} created..."

    def handle_remove(self):
        ww_key = self.select_ww_key()
        if ww_key is None:
            return
        ww = deepcopy(self.weight_windows[ww_key])
        ww.remove_particle()
        ww.file_path = ww.file_path.with_name(f"{ww_key}_1_particle")
        self.weight_windows[ww.file_path.stem] = ww
        self.extra_text = f" {ww.file_path.stem} created..."

    def handle_mitigate(self):
        ww_key = self.select_ww_key()
        if ww_key is None:
            return
        ww = deepcopy(self.weight_windows[ww_key])
        max_ratio = float(input(" Enter the maximum ratio allowed: "))
        ww.mitigate_long_histories(max_ratio)
        ww.file_path = ww.file_path.with_name(f"{ww_key}_mitigated")
        self.weight_windows[ww.file_path.stem] = ww
        self.extra_text = f" {ww.file_path.stem} created..."

    @staticmethod
    def clear_screen():
        if os.name == "nt":
            os.system("cls")
        else:
            os.system("clear")
