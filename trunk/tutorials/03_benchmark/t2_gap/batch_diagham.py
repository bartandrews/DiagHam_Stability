import sys
import csv
import subprocess as sp
from time import perf_counter


if __name__ == "__main__":

    batchfile = sys.argv[1]

    executable = "/home/bart/DiagHam_Stability/trunk/build/FTI/src/Programs/FCI/FCIHofstadterModel"

    with open(batchfile, 'r') as csvfile:
        file = csv.reader(csvfile, delimiter=',')
        for j, row in enumerate(file):
            flags = f" --boson --flat-band --lanczos-precision 1e-10 -n 2 -m 30000 --use-lapack" \
                    f" -p {row[0]} -x {row[1]} -y {row[2]} -X {row[3]} -Y {row[4]}" \
                    f" --t2 {row[5]} --gamma-y {row[6]}"
            command = executable+flags  # define the command
            print(f"Running {command}")
            t0 = perf_counter()  # start the timer
            sp.run(command, stdout=sp.DEVNULL, stderr=sp.DEVNULL, shell=True)  # execute the command
            print(f"Run {j} complete. Time taken: {perf_counter() - t0:.1f} s.")  # end the timer
