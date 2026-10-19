import os, sys, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

ordner = r"D:\Chrome\FP_08_26\FP_08_26"
os.makedirs("plots", exist_ok=True)
ergebnisse = open("ergebnisse.txt", "w", encoding="utf-8")


parula = LinearSegmentedColormap.from_list("parula", ["#352A87", "#0F5CDD", "#1481D6", "#06A4CA",
                                                      "#2EB7A4", "#87BF77", "#D1BB59", "#FEC832", "#F9FB0E"])
okabe_ito = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#000000"]

druecke = [10, 20, 30, 40]
hoehen = [0.5, 1,2,3,4] #mm


ref_size = 5
left_edge = 42
right_edge = 160
mmpp = ref_size / (right_edge - left_edge)

def write(text):
    print(text)
    ergebnisse.write(text + "\n")

def save_fig(fig, title, name, unit, mm_achsen = False):
    plt.figure()
    if mm_achsen:
        rows, columns = fig.shape
        width = columns * mmpp
        hoehe = rows * mmpp
        plt.imshow(fig, cmap=parula, aspect="auto", extent=[-width / 2, width / 2, 0, hoehe])
        plt.xlabel("x [mm]")
        plt.ylabel("z [mm]")
    else:
        plt.imshow(fig, cmap=parula, aspect="auto")
        plt.xlabel("Spalte [px]")
        plt.ylabel("Zeile [px]")
    plt.colorbar(label=unit)
    plt.title(title)
    plt.savefig("plots/" + name + ".pdf", bbox_inches="tight")
    plt.close()


ref_fig = np.loadtxt(glob.glob(ordner + "/Cam/Results/INT*.txt")[0])
plt.figure()
plt.imshow(ref_fig, cmap=parula, aspect="auto")
plt.colorbar(label="Intensitaet [counts]")
plt.axvline(left_edge, color="red")
plt.axvline(right_edge, color="red")
plt.title("Pixel zu mm: 5 mm = 118 px")
plt.xlabel("Spalte [px]")
plt.ylabel("Zeile [px]")
plt.savefig("plots/A3_pixel_zu_mm.pdf", bbox_inches="tight")
plt.close()

write("Aufgabe 3: Pixelgroesse = " + str(round(mmpp, 4)) + " mm/px")

pfad = ordner + "/Processed_data/10barGasDruck/Results/PHA SID4 13h42m03s"
image_gas = np.loadtxt(pfad + "265ms.txt")
image_nogas = np.loadtxt(pfad + "225ms.txt")
image_nogas_after = np.loadtxt(pfad + "305ms.txt")

save_fig(image_gas, "Gasbild ohne Optimierung (10 bar)", "A5_ohne_optimierung", "[lambda]")
save_fig(image_gas - image_nogas, "Gas minus Bild davor (10 bar)", "A5_minus_davor", "[lambda]")
save_fig(image_gas - image_nogas_after, "Gas minus Bild danach (10 bar)", "A5_minus_danach", "[lambda]")

fig, axes = plt.subplots(2, 2, figsize=(9, 6))
for i in range(4):
    druck = druecke[i]
    bild = np.loadtxt(ordner + "/Abbildungen/" + str(druck) + "bar/substracted 1.txt", delimiter=",")
    ax = axes[i // 2][i % 2]
    im = ax.imshow(bild, cmap=parula, aspect="auto", vmin=-0.02, vmax=0.24)
    ax.set_title(str(druck) + " bar")
fig.colorbar(im, ax=axes, label="[lambda]")
plt.savefig("plots/A5_optimiert_alle_druecke.pdf", bbox_inches="tight")
plt.close()


def kosinus(r,n, R):
    if n == 0:
        return np.ones_like(r)
    return 1 - (-1)**n * np.cos(n * np.pi * r / R)

def abel_projektion(y , n, R):
    if abs(y) >= R:
        return 0.0
    t = np.linspace(0, np.sqrt(R ** 2 - y ** 2), 400)
    r = np.sqrt(y ** 2 + t **2)
    return 2 * np.trapezoid(kosinus(r,n,R), t)

def density_calc(druck, number_of_cosines):
    data = glob.glob(ordner + "/Abbildungen/" + str(druck) + "bar/substracted *.txt")
    sum = 0
    for datei in data:
        sum = sum + np.loadtxt(datei, delimiter=",")
    phasedata_average = sum / len(data)

    phasedata_average = phasedata_average[5:123, 28:179]
    rows, cols = phasedata_average.shape

    new_values_phase = np.zeros((rows, cols))
    x = np.arange(cols)
    for i in range(rows):
        gerade = np.polyfit(x, phasedata_average[i], 1)
        new_values_phase[i] = phasedata_average[i] - np.polyval(gerade, x)
    new_values_phase = new_values_phase * 2 * np.pi
    new_values_phase = new_values_phase - new_values_phase.max()

    if number_of_cosines == 10:
        save_fig(phasedata_average, "Phase data cut, " + str(druck) + " bar",
                       "A6_" + str(druck) + "bar_phase_cut", "[lambda]")
        save_fig(new_values_phase, "Phase data normed, " + str(druck) + " bar",
                       "A6_" + str(druck) + "bar_phase_normed", "[rad]", mm_achsen=True)

    mid = 75
    for i in range(rows):
        row = new_values_phase[i]
        beste_spalte = mid
        kleinste_differenz = 1e9
        for c in range(1, cols - 1):
            differenz = abs(row[:c].mean() - row[c:].mean())
            if differenz < kleinste_differenz:
                kleinste_differenz = differenz
                beste_spalte = c
        new_values_phase[i] = np.roll(row, mid - beste_spalte)

    R = cols / 2
    y = np.arange(cols) - cols // 2
    A = np.zeros((cols, number_of_cosines + 1)) 
    B = np.zeros((cols, number_of_cosines + 1))

    for j in range(cols):
        for n in range(number_of_cosines + 1):
            A[j, n] = abel_projektion(y[j], n, R)
            B[j, n] = kosinus(np.array(abs(y[j]), dtype=float), n, R)

    ablinv_whole = np.zeros((rows, cols))
    for i in range(rows):
        koeff = np.linalg.lstsq(A, new_values_phase[i], rcond=None)[0]
        ablinv_whole[i] = B @ koeff


    lam = 0.635            # um
    xi = mmpp * 1000       # um pro Pixel
    density = -ablinv_whole * 1.66541e22 * lam / xi

    name = "A6_" + str(druck) + "bar_N" + str(number_of_cosines)

    save_fig(density, "Gas density, " + str(druck) + " bar, N = " + str(number_of_cosines),
                   name + "_dichte", "Gas density [#/cm^3]", mm_achsen=True)
    write(str(druck) + " bar, N = " + str(number_of_cosines)
             + ": maximale Dichte = " + "%.2e" % density.max() + " cm^-3")


    x_achse = (np.arange(cols) - cols // 2) * mmpp
    z_achse = (rows - 1 - np.arange(rows)) * mmpp

    plt.figure()
    for k in range(len(hoehen)):
        h = hoehen[k]
        zeile = np.argmin(abs(z_achse - h))
        profil = density[zeile - 2:zeile + 3].mean(axis=0)
        ueber_halb = np.where(profil > profil.max() / 2)[0]
        fwhm = (ueber_halb[-1] - ueber_halb[0] + 1) * mmpp
        write("   z = " + str(h) + " mm: n_max = " + "%.2e" % profil.max()
                 + " cm^-3, FWHM = " + "%.2f" % fwhm + " mm")
        plt.plot(x_achse, profil, color=okabe_ito[k], label=str(h) + " mm")


    plt.grid()
    plt.legend()
    plt.xlabel("radius r [mm]")
    plt.ylabel("Gas density [#/cm^3]")
    plt.title("Gas density profiles " + str(druck) + " bar, N = " + str(number_of_cosines))
    plt.savefig("plots/" + name + "_profile.pdf", bbox_inches="tight")
    plt.close()

for druck in druecke:
    density_calc(druck, 10)

for n in [5, 20, 40]:
    density_calc(40, n)

relativer_fehler = (471 - 464) / 467.5



for z in hoehen:
    error = np.sqrt((z * relativer_fehler) ** 2 + mmpp ** 2 / 12)
    write("Hoehenfehler bei z = " + str(z) + " mm: +- " + "%.3f" % error + " mm")

ergebnisse.close()