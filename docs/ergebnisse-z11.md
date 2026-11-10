# Ergebnisse Stufe Z11: Kern auf CPU und GPU

> Server-Auftrag #19 (2026-11-10), Server herkules3: 32 Kerne, RTX 6000 Ada (48 GB), Treiber
> 590.48, nvcc 12.9, CMake 3.31, GCC 14.3, Python 3.12. Stand PR #18.

## Bau und Tests

Bau mit CUDA ohne Fehler. T39 und T40 grün (12 Tests), gesamte Suite 191 Tests grün. Alle Backends
fangen bei beiden Größen dieselben Zellen. Mit einem conda-Python, dessen libstdc++ älter als der
Compiler ist, scheitert der Import (`GLIBCXX_3.4.31 not found`); dafür gibt es die Option
`BINDOPTIX_STATIC_LIBSTDCXX` (wie BWPM).

## Leistung (720 Schritte, CPU mit 16 Threads, unter Fremdlast)

| Zellen | Referenz (Python) | CPU-Kern | CUDA-Kern |
| --- | --- | --- | --- |
| 20 000 | 24,0 s (0,6 M/s) | 0,73 s (19,6 M/s) | 0,59 s (24,3 M/s) |
| 200 000 | 240 s (0,6 M/s) | 8,1 s (17,8 M/s) | 5,7 s (25,3 M/s) |

(M/s = Millionen Zellschritte je Sekunde.)

## Bewertung

Der CPU-Kern ist 30- bis 33-mal schneller als die Referenz. Die GPU gewinnt nur Faktor 1,2–1,4
gegen 16 CPU-Threads. Hauptgrund ist die doppelte Genauigkeit: Die RTX 6000 Ada rechnet FP64 mit
1/64 der FP32-Rate, und je Zellschritt fallen rund 20 Exponential-, Sinus-, Kosinus- und
Logarithmusauswertungen in FP64 an. Dazu kam Fremdlast auf GPU und CPU. Ein Modus mit einfacher
Genauigkeit auf der GPU (Zufallszahlen und Zähler weiter exakt, Physik in FP32, statistisch statt
bitgenau geprüft) verspricht ein bis zwei Größenordnungen; er ist der nächste Schritt, falls die
Rechenzeit drängt.
