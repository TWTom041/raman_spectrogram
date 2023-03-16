import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter
from scipy.integrate import trapz
import os


def peak_stripping(spectrum, window):
    stripping = 0
    y = savgol_filter(spectrum, window, 0)
    n = len(spectrum)
    baseline = np.zeros(n)
    for i in range(n):
        if spectrum[i] > y[i]:
            stripping = 1
            baseline[i] = y[i]
        else:
            baseline[i] = spectrum[i]
    return baseline, stripping


def baseline(spectrum):
    length_of_spectrum = len(spectrum)
    lp = int(np.ceil(0.5 * length_of_spectrum))
    initial_spectrum = np.concatenate((np.ones(lp) * spectrum[0], spectrum, np.ones(lp) * spectrum[-1]))
    l2 = len(initial_spectrum)
    S = initial_spectrum
    n = 1
    flag1 = 0
    A = np.zeros(int(np.floor((l2 - 1) / 2)) + 1)
    stripped_spectrum = {}
    while flag1 == 0:
        n = n + 2
        i = int((n - 1) / 2)
        baseline, stripping = peak_stripping(S, n)
        A[i - 1] = trapz(S - baseline)
        stripped_spectrum[i] = baseline
        S = baseline
        if i > 3:
            if A[i - 2] > A[i - 3] and A[i - 2] > A[i - 1]:
                i_min = i - 2
                flag1 = 1
    base = stripped_spectrum[i_min]
    corrected_spectrum = initial_spectrum - base
    corrected_spectrum = corrected_spectrum[lp:lp + length_of_spectrum]
    base = base[lp:lp + length_of_spectrum]
    return base, corrected_spectrum


def plot(filename):
    current_folder = r"C:\Users\Tom\Downloads\Compressed\222\21\Exp_MatLab\Exp1\PDMS-CW"
    os.makedirs(current_folder, exist_ok=True)
    modes = os.path.splitext(filename)[0].split("_")
    pre, ex = modes[2].split("10")
    pre = pre if pre else "1"

    scale = 1.0

    x = []
    y = []

    with open(os.path.join(current_folder, filename), "r") as f:
        for lines in f.read().split("\n"):
            if lines != "":
                xr, yr = lines.split()
                x.append(float(xr))
                y.append(float(yr))

    x = np.array(x) * scale
    y = np.array(y) * scale

    base, corrected = baseline(y)

    plt.figure()
    plt.plot(x, corrected, label=f"{modes[1]}, scale={scale}")
    plt.plot(x, base, label="baseline")
    plt.plot(x, y, label="original")

    plt.xlabel("Raman Shift ($cm^{-1}$)")
    plt.ylabel("Intensity (arbitrary)")
    plt.title(f"{modes[0]}:{modes[1]} with {modes[3]}=${pre}\\times10^{ex}M$")

    plt.legend()
    plt.show(block=False)
    return


def multi_plot(filenames_path):
    with open(filenames_path, "r") as f:
        list_of_filenames = f.read().split("\n")
        for i in list_of_filenames:
            plot(i.strip() + ".txt")
    return


def main():
    multi_plot(r"C:\Users\Tom\Downloads\Compressed\222\21\Exp_MatLab\Exp1\PDMS-CW_files.txt")
    plt.show()


if __name__ == "__main__":
    main()
