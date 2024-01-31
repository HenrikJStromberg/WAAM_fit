import timeit

# Useful Constants for benchmarking
INPUT_PATH = "examples/Testobjekt_aehnlich_paper.stp"
OUTPUT_PATH  = "./output/Testobjekt"
STEP_SIZE = 0.5

# Arguments for __evalRadius()__
setup = """
import WAAMEvaluator as WAAM_py
import WAAMEvaluator_cy as WAAM_cy
import numpy as np

PATH_TO_N = "examples/benchmark/N"
PATH_TO_cnts = "examples/benchmark/cnts"

index = 0
r = 1.0
N = np.loadtxt(PATH_TO_N)
cnts = np.loadtxt(PATH_TO_cnts)

python = lambda: WAAM_py.__evalRadius__(index, cnts, N, r)
cython = lambda: WAAM_cy.__evalRadius__(index, cnts, N, r)
"""

evalRadius_python = min(timeit.Timer("python()", setup=setup).repeat(100, 100))
evalRadius_cython = min(timeit.Timer("cython()", setup=setup).repeat(100, 100))

print("Best Python run took: " + str(evalRadius_python) + "s")
print("Best Cython run took: " + str(evalRadius_cython) + "s")