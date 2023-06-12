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