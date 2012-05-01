import time

def get_time_remainaing(elapsed, ncomplete, ntotal):
    time_per_n = elapsed / ncomplete
    remaining = (ntotal - ncomplete) * time_per_n
    return "%.3f"%remaining