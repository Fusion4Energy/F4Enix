from f4enix.output.outputAPI import Output
import logging

# you can use loggining to have a more detailed output
logging.basicConfig(level=logging.INFO)

file = r'D:\DATA\laghida\Documents\02_tasks\WAVS\debug\wavs_ref.o'

# Initialize the output object
outp = Output(file)

# print the lost particle debug
outpath = r'D:\DATA\laghida\Documents\02_tasks\WAVS\debug'
outp.print_lp_debug(outpath)
