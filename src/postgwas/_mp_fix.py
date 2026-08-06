import multiprocessing as mp

try:
    mp.set_start_method("fork", force=True)
except RuntimeError:
    pass
