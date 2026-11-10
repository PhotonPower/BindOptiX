# Abhängigkeiten

> Grundsatz des Maintainers: so wenig externe Abhängigkeiten wie möglich; bevorzugt die Programme
> im eigenen GitHub-Konto (PhotonPower), notfalls werden sie erweitert.

## Erlaubt ohne Rückfrage

| Paket | Zweck |
| --- | --- |
| Python ≥ 3.10, NumPy, SciPy | Referenzimplementierung, Numerik |
| pytest, ruff, mypy | Entwicklung (`requirements-dev.txt`) |
| scikit-build-core, pybind11, CMake ≥ 3.24, C++20-Compiler, OpenMP (optional) | Kern `bindoptix._core` (Z11), Aufbau wie BWPM |
| CUDA ≥ 12.6 | CUDA-Backend des Kerns (Z11, nur auf dem großen Rechner, `BINDOPTIX_ENABLE_CUDA=ON`) |

## Programme des Maintainers (lesend; Erweiterungen als Vorschlag)

| Programm | Nutzung in BindOptiX |
| --- | --- |
| BioOptoMD | Dipolensembles Q je Rezeptorzustand, Brechzahlprofil N (über Ausgabeformate, wie der RTT-Fusion-Adapter) |
| hp-FEM | Feldkarte F, Absorption und Erwärmung (M9), Dipol im Schichtstapel (M17, M18) für die Antwortkarte; Vorschlag: Stokes für die Strömung (D7) |
| CliffordBem | Feldkarte F (Nahfeld), Zerfallsraten und Fluoreszenzverstärkung, optische Kräfte auf Zellen (Z8) |
| RTT-Fusion, Raytatouille, Fieldtracing | Kette Fernfeld → SIL → Relay-Optik → Detektor (Pilot K5), Port-Typen, Konventionen |
| Nurbforge | Geometrie des Chips (Kanal, Spots, Kavitäten), Designparameter für die Optimierung |
| Voxelspaetzle | Form und Brechzahl der gedruckten SIL (V6, V10) |
| Scattuccino | Untergrund aus trüben Proben und Autofluoreszenz (später) |
| SpeckleBem | Streulicht rauer Flächen (später) |

## Bewusst nicht verwendet

- Fremde CFD- und FEM-Pakete (OpenFOAM, FEniCS): stattdessen analytische Strömung, Vorschlag an
  hp-FEM oder kleiner eigener Löser (D7).
- Laufzeitkompilierung (NVRTC) für Reaktionsmodelle: Modelle als Template-Parameter, wie die
  Wandformen in BioOptoMD (ADR 0007).
