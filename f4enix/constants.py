import re

# --- PATTERNS ---
PAT_BREAK = re.compile(r'[\s\t]*\n')
PAT_DIGIT = re.compile(r'\d+')
PAT_COMMENT = re.compile(r'[Cc][\s\t]+')
PAT_COMMENT_TEXT = re.compile(r'[Cc]\s.*\n')
PAT_CARD_KEY = re.compile(r'[A-Za-z]+\d+')
PAT_DOLLAR_COMMENT = re.compile(r'\$.*\n')
PAT_MAT = re.compile(r'[\s\t]*[mM]\d+')
PAT_MX = re.compile(r'[\s\t]*mx\d+', re.IGNORECASE)
SCIENTIFIC_PAT = re.compile(r'-*\d.\d+E[+|-]\d+')

# --- Plotter ---
# coordinates are in meters
ITER_Z_LEVELS = [
                ['B2_L', 0, 0, -11.500, 0, 0, 1],
                ['B2_I', 0, 0, -8.425, 0, 0, 1],
                ['B2_U', 0, 0, -4.250, 0, 0, 1],
                ['B1_L', 0, 0, -5.150, 0, 0, 1],
                ['B1_I', 0, 0, -2.625, 0, 0, 1],
                ['B1_U', 0, 0, -1.000, 0, 0, 1],
                ['L1_L', 0, 0, 0.100, 0, 0, 1],
                ['L1_I', 0, 0, 2.690, 0, 0, 1],
                ['L1_U', 0, 0, 4.380, 0, 0, 1],
                ['L2_L', 0, 0, 5.480, 0, 0, 1],
                ['L2_I', 0, 0, 7.970, 0, 0, 1],
                ['L2_U', 0, 0, 9.560, 0, 0, 1],
                ['L3_L', 0, 0, 10.660, 0, 0, 1],
                ['L3_I', 0, 0, 14.8375, 0, 0, 1],
                ['L3_U', 0, 0, 16.200, 0, 0, 1],
                ['L4_L', 0, 0, 19.215, 0, 0, 1],
                ['L4_I', 0, 0, 22.8375, 0, 0, 1],
                ['L4_U', 0, 0, 25.560, 0, 0, 1]
                ]
# Color Maps
TID_COLORS = ['green', 'orange', 'red']
TID_LIMITS = [1, 10]  # Gy
TID_LABELS = ['< 1 Gy', '1-10 Gy', '> 10 Gy']

TNF_COLORS = ['green', 'orange', 'red']
TNF_LIMITS = [0.01, 100]  # n/cm^2/s
TNF_LABELS = ['< 0.01 n/cm^2/s', '0.01-100 n/cm^2/s', '> 100 n/cm^2/s']

SDDR_COLORS = ['white', 'blue', 'green', 'yellow', 'orange', 'red']
SDDR_LIMITS = [1.32e-3, 7.5e-3, 25e-3, 2, 100]  # mSv/h
SDDR_LABELS = ['< 1.32 uSv/h', '1.32-7.5 uSv/h', '7.5-25 uSv/h',
               '0.025-2 mSv/s', '2-100 mSv/h', '> 100 mSv/h']
