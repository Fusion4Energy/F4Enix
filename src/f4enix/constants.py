"""
Constants and patterns used throughout F4Enix

Attributes
----------
ITER_Z_LEVELS: list[list[str, float, float, float, float, float, float]]
    list of parameters for the general slicing at Z Iter levels (e.g. B2L)
    format is [name, origin[x,y,x], norm[x, y, z]]
ITER_TOROIDAL_SLICES: list[list[str, float, float, float, float, float, float]]
    list of parameters for the general toroidal slicing (e.g. Port #1/#10)
    format is [name, origin[x,y,x], norm[x, y, z]]

"""

"""
Copyright 2019 F4E | European Joint Undertaking for ITER and the Development of
Fusion Energy (‘Fusion for Energy’). Licensed under the EUPL, Version 1.2 or - 
as soon they will be approved by the European Commission - subsequent versions
of the EUPL (the “Licence”). You may not use this work except in compliance
with the Licence. You may obtain a copy of the Licence at: 
    https://eupl.eu/1.2/en/  
Unless required by applicable law or agreed to in writing, software distributed
under the Licence is distributed on an “AS IS” basis, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the Licence permissions
and limitations under the Licence.
"""

import os
import re
from pathlib import Path

# --- Typing ---
PathLike = Path | str | os.PathLike

# --- PATTERNS ---
PAT_BREAK = re.compile(r"[\s\t]*\n")
PAT_BLANK = re.compile(r"[\s\tCc]*\n")
PAT_SPACE = re.compile(r"[\s\t]+")
PAT_DIGIT = re.compile(r"\d+")
PAT_COMMENT = re.compile(r"[Cc][\s\t\n]+")
PAT_COMMENT_TEXT = re.compile(r"[Cc]\s.*\n")
PAT_CARD_KEY = re.compile(r"[A-Za-z]+\d+")
PAT_F_TR_CARD_KEY = re.compile(r"(FMESH|TR|F)\d+")
PAT_DOLLAR_COMMENT = re.compile(r"\$.*\n")
PAT_MAT = re.compile(r"[\s\t]*[mM]\d+")
PAT_MX = re.compile(r"[\s\t]*mx\d+", re.IGNORECASE)
SCIENTIFIC_PAT = re.compile(r"[-+]*\d.\d+E[+-]\d+")
PAT_FMESH_KEY = re.compile(r"^FMESH\d+")
PAT_NP = re.compile(r"(?<=:)[nN,pP]+")

# --- Plotter ---
# coordinates are in meters
# ITER_Z_LEVELS = [
#     ["B2_L (pz=-11.500 m)", 0, 0, -11.500, 0, 0, 1],
#     ["B2_I (pz=-8.425 m)", 0, 0, -8.425, 0, 0, 1],
#     ["B2_U (pz=-4.250 m)", 0, 0, -4.250, 0, 0, 1],
#     ["B1_L (pz=-5.150 m)", 0, 0, -5.150, 0, 0, 1],
#     ["B1_I (pz=-2.625 m)", 0, 0, -2.625, 0, 0, 1],
#     ["B1_U (pz=-1.000 m)", 0, 0, -1.000, 0, 0, 1],
#     ["L1_L (pz=-0.100 m)", 0, 0, 0.100, 0, 0, 1],
#     ["L1_I (pz=2.690 m)", 0, 0, 2.690, 0, 0, 1],
#     ["L1_U (pz=4.380 m)", 0, 0, 4.380, 0, 0, 1],
#     ["L2_L (pz=5.480 m)", 0, 0, 5.480, 0, 0, 1],
#     ["L2_I (pz=7.970 m)", 0, 0, 7.970, 0, 0, 1],
#     ["L2_U (pz=9.560 m)", 0, 0, 9.560, 0, 0, 1],
#     ["L3_L (pz=10.660 m)", 0, 0, 10.660, 0, 0, 1],
#     ["L3_I (pz=14.8375 m)", 0, 0, 14.8375, 0, 0, 1],
#     ["L3_U (pz=16.200 m)", 0, 0, 16.200, 0, 0, 1],
#     ["L4_L (pz=19.215 m)", 0, 0, 19.215, 0, 0, 1],
#     ["L4_I (pz=22.8375 m)", 0, 0, 22.8375, 0, 0, 1],
#     ["L4_U (pz=25.560 m)", 0, 0, 25.560, 0, 0, 1],
# ]
# updated values from ITER_D_3FM52L (mid values)
# see discussion in issue #71
ITER_Z_LEVELS = [
    ["B2_Floor (pz=-12.62 m)", 0, 0, -12.62, 0, 0, 1],
    ["B2_Mid (pz=-10.77 m)", 0, 0, -10.77, 0, 0, 1],
    ["B2_Ceiling (pz=-7.99 m)", 0, 0, -7.99, 0, 0, 1],
    ["B1_Floor (pz=-6.30 m)", 0, 0, -6.30, 0, 0, 1],
    ["B1_Mid (pz=-4.56 m)", 0, 0, -4.56, 0, 0, 1],
    ["B1_Ceiling (pz=-2.82 m)", 0, 0, -2.82, 0, 0, 1],
    ["L1_Floor (pz=-1.02 m)", 0, 0, -1.02, 0, 0, 1],
    ["L1_Mid (pz=0.81 m)", 0, 0, 0.81, 0, 0, 1],
    ["L1_Ceiling (pz=2.64 m)", 0, 0, 2.64, 0, 0, 1],
    ["L2_Floor (pz=4.34 m)", 0, 0, 4.34, 0, 0, 1],
    ["L2_Mid (pz=6.09 m)", 0, 0, 6.09, 0, 0, 1],
    ["L2_Ceiling (pz=7.84 m)", 0, 0, 7.84, 0, 0, 1],
    ["L3_Floor (pz=9.58 m)", 0, 0, 9.58, 0, 0, 1],
    ["L3_Mid (pz=12.58 m)", 0, 0, 12.58, 0, 0, 1],
    ["L3_Ceiling (pz=15.58 m)", 0, 0, 15.58, 0, 0, 1],
    ["L4_Crane_Hall (pz=19.01 m)", 0, 0, 19.01, 0, 0, 1],
    ["L4_Mid (pz=20.99 m)", 0, 0, 20.99, 0, 0, 1],
    ["L4_Ceiling (pz=22.98 m)", 0, 0, 22.98, 0, 0, 1],
]
ITER_TOROIDAL_SLICES = [
    ["Ports #1/#10", 0, 0, 0, -0.2727659196338417, 1.5469324010307122, 0.0],
    ["Ports #2/#11", 0, 0, 0, -0.785398163397448, 1.3603495231756635, 0.0],
    ["Ports #3/#12", 0, 0, 0, -1.2032997974129325, 1.0096884162048878, 0.0],
    ["Ports #4/#13", 0, 0, 0, -1.4760657170467744, 0.5372439848258247, 0.0],
    ["Ports #5/#14", 0, 0, 0, -1.5707963267948966, 1.9236706937217898e-16, 0.0],
    ["Ports #6/#15", 0, 0, 0, -1.4760657170467746, -0.5372439848258244, 0.0],
    ["Ports #7/#16", 0, 0, 0, -1.2032997974129327, -1.0096884162048876, 0.0],
    ["Ports #8/#17", 0, 0, 0, -0.785398163397449, -1.360349523175663, 0.0],
    ["Ports #9/#18", 0, 0, 0, -0.2727659196338418, -1.5469324010307122, 0.0],
]
# Color Maps
TID_CATEGORIES = {
    "colors": ["green", "orange", "red"],
    "values": [1, 10],  # Gy
    "categories": ["< 1 Gy", "1-10 Gy", "> 10 Gy"],
}

TNF_CATEGORIES = {
    "colors": ["green", "orange", "red"],
    "values": [0.01, 100],  # n/cm^2/s
    "categories": ["< 0.01 n/cm^2/s", "0.01-100 n/cm^2/s", "> 100 n/cm^2/s"],
}

SDDR_CATEGORIES = {
    "colors": ["white", "blue", "green", "yellow", "orange", "red"],
    "values": [1.32e-3, 7.5e-3, 25e-3, 2, 100],  # mSv/h
    "categories": [
        "< 1.32 uSv/h",
        "1.32-7.5 uSv/h",
        "7.5-25 uSv/h",
        "0.025-2 mSv/s",
        "2-100 mSv/h",
        "> 100 mSv/h",
    ],
}

# --- Constants ---
ITER_PLASMA_360_500MW = 1.7757e20


# --- Variables for E-Lite ---
SECTOR_BOUNDARIES = {
    1: [427024, 427016],
    "2 & 3": [-437544, 437543],
    4: [451024, 451016],
    5: [459024, 459016],
    6: [467024, 467016],
    7: [475024, 475016],
    8: [483024, 483016],
    9: [491024, 491016],
}

SECTOR_BOUNDARIES_ANGLES = {
    1: [10, 50],
    "2 & 3": [50, 130],
    4: [130, 170],
    5: [170, -150],
    6: [-150, -110],
    7: [-110, -70],
    8: [-70, -30],
    9: [-30, 10],
}

PLASMA_CELLS = {
    1: 427001,
    "2 & 3": 435001,
    4: 451001,
    5: 459001,
    6: 467001,
    7: 475001,
    8: 483001,
    9: 491001,
}

SECTOR_NAMES = [1, "2 & 3", 4, 5, 6, 7, 8, 9]

# TALLY_CARD_TYPES = ['F', '+F', '*F', 'FC', 'E', 'T', 'C', 'FQ',
#                     'FM', 'DE', 'DF', 'EM', 'TM', 'CM', 'CF',
#                     'SF', 'FS', 'SD', 'FU', 'FT', 'TF', 'NOTRN']

UNION_INTERSECT_SYMBOLS = {"union": ":", "intersect": ""}

CONV = {
    "Result": "Value",
    "Rel": "Error",
    "R": "Cor A",
    "Z": "Cor B",
    "Theta": "Cor C",
    "Energy": "Energy",
}

# isotopes with negative heating in FENDL 3.1x
FAULTY_ISOTOPES = [
    "31069",
    "31071",
    "35079",
    "35081",
    "42092",
    "42094",
    "42095",
    "42096",
    "42097",
    "42098",
    "42100",
    "45103",
    "48108",
    "48110",
    "48112",
    "48113",
    "48114",
    "48116",
    "51121",
    "51123",
    "56130",
    "56132",
    "56134",
    "56135",
    "56136",
    "56137",
    "58142",
    "72174",
    "72176",
    "72177",
    "72178",
    "72179",
    "72180",
]

AVOGADRO_NUMBER = 6.0220434469282e23
