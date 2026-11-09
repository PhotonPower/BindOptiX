# Erste Ergebnisse Stufe Z6: Photonenbudget einer gefangenen Zelle

> Stand 2026-11-09. Skript: `experiments/z6_photonenbudget.py`. Modell: `bindoptix.signal`
> (T17–T20). Optische Karten sind Platzhalter (`placeholder_response_table`,
> `uniform_field_map`), bis hp-FEM und CliffordBem exportieren (Z2); alle Werte **P**.

## Annahmen

Cy5-ähnlicher Farbstoff (ε = 250 000 /(M·cm), q₀ = 0,27, Anregung 640 nm, Emission 670 nm,
Bleichausbeute 10⁻⁶), 10 mW gleichmäßig auf das Gesichtsfeld Ø 5 mm (509 W/m²), Sammeleffizienz
von SIL und Relay-Optik 30 %, Kamera: 2048 Pixel über 5 mm (2,44 µm je Pixel im Objekt),
Quanteneffizienz 0,8, Lesefehler 1,5 e⁻, Punktbild σ = 1,2 µm. Aptamerdichte 10⁴ /µm². Weitere
Untergründe (Photolumineszenz des Golds, Lack der SIL, Glas) zusammen 0,2 Photonen/(µm²·s).

## Ergebnis

Ein Farbstoff wird 155-mal je Sekunde angeregt, davon kommen 12,6 Photonen je Sekunde an. Bei
dauernder Beleuchtung bleicht er im Mittel nach 1,8 h.

Kürzeste Belichtung, bei der eine Zelle mit 95 % Wahrscheinlichkeit erkannt wird, bei höchstens
0,05 Fehlalarmen je mm² durch Photonen- und Lesefehler:

| eingeschaltete Farbstoffe je Zelle | Restleuchten aus: 0,1 % | 1 % | 5 % |
| --- | --- | --- | --- |
| 20 | 10 s | 100 s | 300 s |
| 200 | 0,1 s | 1 s | 3 s |
| 2 500 | 0,01 s | 0,01 s | 0,03 s |

(2 500 entspricht einer Kontaktfläche von 0,5 µm², von deren Aptameren die Hälfte schaltet.)

## Befunde

1. **Photonenrauschen ist bei diesen Annahmen nicht das Problem.** Schaltet eine Zelle Hunderte
   bis Tausende Aptamere, genügen Belichtungen unter einer Sekunde.
2. **Die beiden entscheidenden Unbekannten** sind die Zahl der Aptamere, die eine Zelle
   tatsächlich einschaltet, und das Restleuchten der ausgeschalteten Aptamere (Löschung durch das
   Gold). Fallen beide ungünstig aus, steigt die nötige Belichtung um vier Größenordnungen. Mit
   echten Antwortkarten folgt das Restleuchten aus den Höhen der BioOptoMD-Ensembles über dem Gold
   (Z2); messen lässt es sich direkt (Messplan M7).
3. **Bleichen** begrenzt die Dauerbeleuchtung auf rund eine Stunde; Bilder in Abständen (z. B.
   1 s alle 5 min) kosten praktisch nichts.
4. **Fehlalarme kommen deshalb vor allem aus Strukturen**, nicht aus dem Rauschen:
   autofluoreszente Partikel, Aggregate, ungleichmäßiges Restleuchten über den Spot. Gegenmittel
   sind Differenzbilder gegen eine Aufnahme zu Beginn und ein zeitliches Kriterium (ein Signal
   erscheint und bleibt). Das ist Aufgabe der Auslese in Z7.
