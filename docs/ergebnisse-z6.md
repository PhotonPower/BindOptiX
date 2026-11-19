# Erste Ergebnisse Stufe Z6: Photonenbudget einer gefangenen Zelle

> Stand 2026-11-09. Skript: `experiments/z6_photonenbudget.py`. Modell: `bindoptix.signal`
> (T17–T20). Optische Karten sind Platzhalter (`placeholder_response_table`,
> `uniform_field_map`), bis hp-FEM und CliffordBem exportieren (Z2); alle Werte **P**. Mit Karten
> aus hp-FEM: Abschnitt am Ende.

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

## Mit Karten aus hp-FEM (2026-11-19)

Skript: `experiments/z6_echte_karten.py`. Anordnung nach D4: Anregung von der Lösungsseite
senkrecht (640 nm, x-polarisiert), Sammlung durch das Substrat. Stapel Wasser / 50 nm Gold (P,
bis D4) / Glas 1,52, Gold nach johnson1972 aus der Materialbibliothek von hp-FEM. Feldkarte und
Antworttabelle (24 Höhen 1–300 nm) aus hp-FEM Stufe 2 (#23); die Feldkarte stimmt mit einer
eigenen Airy-Rechnung überein (|E(0)|/E₀ = 0,708, R = 0,872). Farbstoff, Leistung, Kamera und
Höhen weiter **P**; Farbstoffe isotrop gemittelt; SIL und Relay-Optik ohne Verluste außer der
Winkelauswahl (vorher pauschal 30 %).

| Höhe über dem Gold | Anregung gegen ohne Gold | Quantenausbeute (q₀ = 0,27) | Anteil von A_rad ins Substrat, NA < 1,333 / ≥ 1,333 | erkannte Photonen je s und Farbstoff |
| --- | --- | --- | --- | --- |
| 2 nm | 0,57 | 0,002 | 0,025 / 0,56 | 0,10 |
| 5 nm | 0,68 | 0,031 | 0,024 / 0,55 | 1,5 |
| 10 nm | 0,88 | 0,117 | 0,023 / 0,52 | 7,2 |
| 20 nm | 1,32 | 0,190 | 0,022 / 0,47 | 15,9 |

(Platzhalterkarten vorher: 12,6 je s mit 30 % Sammlung.)

Kürzeste Belichtung wie oben, Farbstoffe 5 nm über dem Gold:

| eingeschaltete Farbstoffe je Zelle | Restleuchten aus: 0,1 % | 1 % | 5 % |
| --- | --- | --- | --- |
| 20 | 100 s | – | – |
| 200 | 1 s | 10 s | 30 s |
| 2 500 | 0,03 s | 0,1 s | 0,3 s |

(– : nicht erreicht bis 300 s.)

### Befunde

1. **Das Gold löscht unterhalb von 10 nm stark.** Die Quantenausbeute fällt von 0,27 auf 0,12
   bei 10 nm, 0,03 bei 5 nm und 0,002 bei 2 nm. Der Abstand des Farbstoffs vom Gold (Länge von
   Anker und Aptamer) wird damit zum wichtigsten Parameter des Signals; er kommt aus den
   BioOptoMD-Ensembles.
2. **Stehwelle bei Anregung von oben.** Das Gold spiegelt 87 %; an der Oberfläche ist die
   Intensität nur halb so groß wie ohne Gold, bei 20 nm 1,3-mal so groß.
3. **Gesammelt wird fast nur SPCE.** Von der Leistung im Substrat liegen über 95 % bei NA ≥ 1,333
   (Kopplung über das Oberflächenplasmon), unter NA 1,333 nur 2–3 % von A_rad. Eine SIL, die nur
   bis NA 1,33 aufnimmt, verlöre also rund das Zwanzigfache. Die NA der SIL gehört deshalb in D4.
4. **Restleuchten aus dem Gold allein.** Säße der ausgeschaltete Farbstoff 2 nm und der
   eingeschaltete 10 nm über dem Gold (P), wäre das Verhältnis 0,10/7,2 ≈ 1,4 %, also in der
   Spalte „1 %“ der Tabellen. Ob ein Aptamer beim Schalten den Abstand so ändert, ist offen
   (Messplan M7).

## Reihe über Golddicke und Haftschicht (2026-11-19)

Exporte aus hp-FEM für 15 Stapel (#34): Wasser / Au d / [Ti 2 nm oder Cr 2 nm] / Glas,
d = 5–50 nm; Ti und Cr nach johnson1974. Annahmen wie oben; `python experiments/z6_echte_karten.py
DIR`. Erkannte Photonen je s und eingeschaltetem Farbstoff, SIL nimmt das ganze Substrat auf
(„alle“) bzw. nur NA < 1,333:

| Gold | 5 nm Höhe: alle / < 1,333 | 10 nm | 20 nm |
| --- | --- | --- | --- |
| 5 nm | 0,48 / 0,29 | 1,08 / 0,67 | 3,8 / 2,4 |
| 10 nm | 1,28 / 0,67 | 3,07 / 1,63 | 6,9 / 3,8 |
| 20 nm | 1,90 / 0,55 | 6,78 / 2,02 | 14,0 / 4,4 |
| 30 nm | **2,03** / 0,28 | **8,69** / 1,26 | **18,5** / 2,9 |
| 50 nm | 1,52 / 0,07 | 7,23 / 0,33 | 15,9 / 0,8 |
| 30 nm + Ti 2 nm | 1,73 / 0,22 | 7,48 / 0,97 | 16,4 / 2,3 |
| 30 nm + Cr 2 nm | 1,73 / 0,22 | 7,41 / 1,00 | 16,2 / 2,3 |

Bei 2 nm Höhe liegen alle Stapel bei 0,08–0,17 je s (Löschung). Mit Haftschicht für die
übrigen Dicken siehe Skriptausgabe; Ti und Cr unterscheiden sich um höchstens 8 %.

### Befunde

1. **30 nm Gold geben das stärkste Signal**, wenn die SIL das ganze Substrat aufnimmt; 20 nm
   liegen 6–24 % darunter, 50 nm 14–25 %. Dünnes Gold (5–10 nm) ist schlechter, obwohl es weniger
   spiegelt: Es löscht stärker (A_tot bei 20 nm Höhe 27 für 5 nm Gold gegen 3,2 für 30 nm).
2. **Die NA der SIL entscheidet.** Nimmt sie nur bis NA 1,333 auf, verliert 30 nm Gold einen
   Faktor 6–7, und dünneres Gold (10–20 nm) wird besser, bleibt aber 3–4-mal unter dem vollen
   Wert. Eine gedruckte SIL auf dem Substrat sollte deshalb bis NA ≈ 1,45 sammeln.
3. **Eine Haftschicht von 2 nm kostet bei 30 nm Gold etwa 11–15 %**, bei dünnerem Gold bis knapp 40 %;
   Ti und Cr sind gleichwertig. Ein molekularer Haftvermittler vermiede diesen Verlust.
4. Die Volumendaten gelten für geschlossene Filme; die Werte für 5 nm (und eingeschränkt 10 nm)
   sind deshalb nur grob.
