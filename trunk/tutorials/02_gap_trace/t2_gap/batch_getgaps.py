import sys
import csv
import os

if __name__ == "__main__":

    batchfile = sys.argv[1]
    outfile = sys.argv[2]

    data = [["Np", "UCarea", "t2", "aniso", "twist", "gap"]]

    with open(batchfile, 'r') as csvfile:  # open batchfile
        file = csv.reader(csvfile, delimiter=',')
        for i, row in enumerate(file):
            evals = []
            n, x, y, X, Y, t2, gy = int(row[0]), int(row[1]), int(row[2]), int(row[3]), int(row[4]), \
                                    float(row[5]), float(row[6])
            filename = (f"bosons_hofstadter_X_{X}_Y_{Y}_q_1_"
                        f"n_{n}_x_{x}_y_{y}_t2_{t2:g}_t3_0_alpha_1_u_1_gx_0_gy_{gy:g}.dat")
            if os.path.isfile(filename):
                print(f"Opening {filename}")
                with open(filename, 'r') as csvfile2:  # open dat file
                    file2 = csv.reader(csvfile2, delimiter=' ')
                    for row2 in file2:
                        if row2[1].isnumeric():
                            evals.append(float(row2[2]))
                    s = sorted(evals)
                    delta = float(s[1] - s[0])
                    data.append([n, X*Y, t2, (y*Y)/(x*X), gy, X*Y*delta])
            else:
                print(f"File {filename} not found.")

    with open(outfile, 'w') as csvfile3:  # open outfile
        file3 = csv.writer(csvfile3, delimiter=',')
        for i in data:
            file3.writerow(i)
