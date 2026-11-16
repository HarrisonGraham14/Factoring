# test code for algorithm, results and code for generation of figures is in results.py
import time
import random
import math
from primality import primality
import multiprocessing
import matplotlib.pyplot as plt

def randomPrime(min, max):
    p = 1
    while primality.isprime(p) == False:
        p = random.randint(min, max)
    return p

def timeAlgorithm(algorithm, args, queue):
    t_start = time.time()
    algorithm(*args)
    t_end = time.time()
    queue.put(t_end - t_start)

def plotAsymptotics(algorithm, targetTime = 3600, printLogs = True):
    t_start = time.time()

    # axes for plot
    N_values = []
    time_values = []

    # phase 1: factoring products of d-digit primes
    digits = 0
    if printLogs:
        print("\nPhase 1 started")
    while True:
        digits += 0.5
        if printLogs:
            print(f"\nFactoring product of primes with {digits} digits")

        pMin = math.floor(10**(digits - 0.5))
        pMax = math.ceil(10**(digits))

        time_allowed = (targetTime - (time.time() - t_start)) / 100

        p1 = randomPrime(pMin, pMax)
        p2 = randomPrime(pMin, pMax)
        N = p1 * p2
        if printLogs:
            print(f"{p1} * {p2}... ", end="")

        queue = multiprocessing.Queue()
        process = multiprocessing.Process(target = timeAlgorithm, args = (algorithm, (N,), queue))
        process.start()
        process.join(time_allowed + 5)
        if not process.is_alive():
            sampleTime = queue.get()
            N_values.append(N)
            time_values.append(sampleTime)
            if printLogs:
                print(f"finished in {sampleTime}s")
            if sampleTime >= time_allowed:
                break
        else:
            if printLogs:
                print(f"not finished")
            break

    # phase 2: fill in the plot with 100 more evenly spaced samples
    max_value = max(N_values)
    bin_width = max_value / 101
    N_values = [] # discard, since multithreading gives different timing, making noise on the graph
    time_values = []
    BSGS_count_values = []
    warnings = []
    if printLogs:
        print("\n\nPhase 2 started")
    for i in range(100):
        if printLogs:
            print(f"\nFactoring a number in bin {i + 1}/100")

        pMin = min(2, math.floor(math.sqrt(i * bin_width) - 100))
        pMax = math.floor(math.sqrt((i + 1) * bin_width)) + 100

        p1 = randomPrime(pMin, pMax)
        p2 = randomPrime(pMin, pMax)
        N = p1 * p2
        if printLogs:
            print(f"{p1} * {p2}... ", end="")

        t_start = time.time()
        x, BSGS_count = algorithm(N)
        sampleTime = time.time() - t_start

        N_values.append(N)
        time_values.append(sampleTime)
        BSGS_count_values.append(BSGS_count)

        if N % x != 0:
            warnings.append(f"extracted non-factor {x} from {N}")
        if x == 1 or x == N:
            warnings.append(f"extracted trivial factor {x} from {N}")

        if printLogs:
            print(f"finished in {round(sampleTime)}s")


    print(f"\nN_values = {N_values}")
    print(f"time_values = {time_values}")
    print(f"BSGS_count_values = {BSGS_count_values}")

    plt.plot(N_values, time_values, "ko")
    plt.show()

    for warning in warnings:
        print("WARNING: " + warning)

        

from factor import factor
if __name__ == "__main__":
    plotAsymptotics(factor, 3*60*60)