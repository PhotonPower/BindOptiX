# Erste Ergebnisse Stufe Z10: Klassifikation der Stämme

> Stand 2026-11-09. Skript: `experiments/z10_klassifikation.py`. Modell: Trajektorienbibliothek
> (`bindoptix.classify`, geprüft durch T38 gegen direkte Läufe), Bindungsmatrix, Klassifikator
> nach der multinomialen Likelihood (T25). Synthetische Bindungsmatrix: spezifische Fangrate
> 0,01 /s, Kreuzreaktivität als Anteil davon, lognormal gestreut (σ = 1). Alle Werte **P**.

## Aufbau

Schräges 5 × 5-Grid, ein Kanal je Gridreihe, Ablauf wie Z7 (0,2 µL Konzentrat je Kanal laden,
60 min Inkubation ohne Strömung, 10 min Waschen), Bibliothek aus 600 Zellen je Kanal.
„Genauigkeit“ ist der Anteil richtig erkannter Proben (ohne Entscheidung zählt als falsch);
„Zellen“ ist die Zahl der Zellen im Konzentrat (über alle Kanäle).

## 25 Stämme, je ein Aptamer pro Spot

| Fluss | Kreuzreaktivität | gefangen je Zelle | Genauigkeit 50 Zellen | 200 Zellen |
| --- | --- | --- | --- | --- |
| 5 nL/min | 1 % | 3,3 % | 19 % | 50 % |
| 5 nL/min | 10 % | 9,1 % | 11 % | 27 % |
| 5 nL/min | 30 % | 11,5 % | 8 % | 15 % |

1,5 und 15 nL/min geben ähnliche Werte (±3 Prozentpunkte).

## 5 Klassen, alle Aptamere in jeder Reihe (lateinisches Quadrat)

| Kreuzreaktivität | gefangen je Zelle | Genauigkeit 10 Zellen | 50 Zellen | 200 Zellen |
| --- | --- | --- | --- | --- |
| 1 % | 5,3 % | 22 % | 60 % | 91 % |
| 10 % | 10,3 % | 18 % | 33 % | 48 % |
| 30 % | 11,9 % | 17 % | 27 % | 38 % |

## Befunde

1. **Der Engpass ist die Zahl gefangener Zellen.** Selbst wenn jeder Spot jede Zelle fängt,
   erreichen in diesem Ablauf nur 13–14 % der Zellen überhaupt einen Spot (auch mit Strömung statt
   Inkubation). Am Maßnahmenwert (100 KBE/100 mL, 50 % Rückhalt → etwa 50 Zellen) sind das rund
   sieben gefangene Zellen. Für den Nachweis „Legionellen ja/nein“ mit einem Aptamer auf allen
   Spots reicht das (P(mindestens eine) ≈ 99,9 %), für 25 Subtypen nicht.
2. **Im Layout „ein Kanal je Reihe“ sieht jede Zelle nur 5 der 25 Aptamere.** Das kostet beim
   spezifischen Fang den Faktor 5. Aptamere in jeder Reihe zu wiederholen (lateinisches Quadrat)
   hebt die Genauigkeit für 5 Klassen auf 91 % bei 200 Zellen.
3. **Kreuzreaktivität ist der zweite Hebel.** Schon 10 % halbieren die Genauigkeit. Lange
   Verweilzeiten auf einem Spot (Gleiten, ≈ 46 min je Spot) machen auch schwache Bindung
   wahrscheinlich und verschärfen das; kürzere Verweilzeiten würden spezifisch und
   kreuzreaktiv besser trennen, kosten aber Fänge.
4. **Folgerung für das Vorgehen:** gestufte Klassifikation (Legionellen ja/nein → L. pneumophila
   → Serogruppe 1) mit wiederholten Aptameren je Reihe; die Unterscheidung von 25 Subtypen braucht
   deutlich mehr Zellen (Anreicherung) oder eine Kanalführung, bei der jede Zelle alle Aptamere
   passiert (D15). Die Bindungsmatrix aus M5 ersetzt die synthetische.
