# Erste Ergebnisse Stufe Z7: Ablauf und digitale Auslese

> Stand 2026-11-09. Skript: `experiments/z7_ende_zu_ende.py`. Kette: Protokoll
> (`bindoptix.protocol`) → Partikel mit Einfang (Z3) → Bilder (`bindoptix.imaging`, Photonen nach
> Z6) → Differenzbild, Schwelle, Persistenz (`bindoptix.readout`) → Zählung je Spot. Geprüft durch
> T21, T22, T37. Alle Parameter **P**.

## Ablauf

Mittlerer Kanal des schrägen Grids. Laden 40 min bei 5 nL/min mit 100 Zellen, gleichmäßig
eingespült; Inkubation 60 min ohne Strömung; Waschen 10 min bei 20 nL/min. Aufnahmen (1 s) alle
10 min und eine Bestätigungsaufnahme 1 min nach dem Ende. Signal und Untergrund wie Z6 mit
200 eingeschalteten Aptameren je Zelle und 1 % Restleuchten; Ungleichmäßigkeit der Spots 30 %.

## Ergebnis

| Zeit [min] | gefangen (wahr) | gezählt |
| --- | --- | --- |
| 20 | 1 | 1 |
| 50 | 5 | 5 |
| 80 | 8 | 8 |
| 100 | 11 | 11 |
| 111 | 15 | 15 |

Je Spot (stromauf → stromab): wahr 6, 3, 4, 2, 0 – gezählt 6, 3, 4, 2, 0.

## Befunde

1. **Die Auslese zählt richtig, wenn drei Dinge stimmen:** Differenzbild gegen eine Aufnahme vor
   dem Laden (sonst erzeugt die Ungleichmäßigkeit der Spots über 100 Fehlalarme in einem kleinen
   Ausschnitt, T21), Persistenz über zwei Aufnahmen (sonst zählen kurz sichtbare Partikel mit) und
   eine Bestätigungsaufnahme nach dem Ende.
2. **Viele Zellen werden an der Vorderkante eines Spots gefangen**, weil sie gleitend ankommen
   und sofort binden. Ihre gemessene Lage liegt dann oft knapp außerhalb des Spots. Ohne
   Toleranz beim Zuordnen fehlten im ersten Lauf zwei von drei Zellen eines Spots; mit 5 µm Rand
   stimmt jede Zählung. Für die Optik heißt das: Die Spotkanten müssen scharf abgebildet und gut
   bekannt sein (Registrierung), und Kanten sind Orte hoher Signaldichte.
3. **15 von 100 Zellen** werden in diesem Ablauf gefangen, im Rahmen der Z3-Ergebnisse
   (5 nL/min, Fang sofort, kein Gleiten während der Inkubation ohne Strömung).

## Mit dem Photonenbudget nach D4 (2026-11-19)

Signal je eingeschaltetem Farbstoff 7,1 Photonen/s statt 12,6 (Z6: 30 nm Gold, SIL bis NA 1,45,
Farbstoff 10 nm über dem Gold, `docs/ergebnisse-z6.md`); 200 eingeschaltete Aptamere je Zelle,
1 % Restleuchten und die Höhe bleiben **P**. Ablauf und Auslese wie oben.

| Belichtung | gezählt am Ende (wahr 15) | je Spot (wahr 6, 3, 4, 2, 0) |
| --- | --- | --- |
| 1 s | 13 | 5, 3, 4, 1, 0 |
| 2 s | 15 | 6, 3, 4, 2, 0 |
| 3 s, 10 s | 15 | 6, 3, 4, 2, 0 |

Das Skript belichtet jetzt 2 s; damit stimmt jede Zählung zu jedem Zeitpunkt wie im ersten Lauf.
Korrektur dabei: Die Photonen der Zellen wurden nicht mit der Belichtungszeit multipliziert (der
Untergrund schon); bei 1 s ohne Folgen, die Ergebnisse oben bleiben gültig.

**Befund:** Mit dem echten Budget reicht 1 s nicht mehr sicher; 2 s genügen bei 200
eingeschalteten Aptameren. Sitzen die Farbstoffe nur 5 nm über dem Gold (1,7 Photonen/s), wären
es nach Z6 rund 10 s. Die Höhe aus BioOptoMD bestimmt also direkt die Belichtung und damit, wie
oft während des Ablaufs aufgenommen werden kann (Bleichen, Z6).
