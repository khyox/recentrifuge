import multiprocessing
import os

METHODS = ['spawn', 'fork', 'forkserver']

if __name__ == '__main__':
    for method in METHODS:
        context = multiprocessing.get_context(method)
        with context.Pool(processes=4) as pool:
            results = [pool.apply_async(os.getpid, ()) for i in range(4)]
            print([res.get(timeout=1) for res in results], method)

