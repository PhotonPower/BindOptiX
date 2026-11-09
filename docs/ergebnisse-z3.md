# Erste Ergebnisse Stufe Z3: Einfang im Kanal

> Stand 2026-11-09. Skript: `experiments/z3_einfang_flussrate.py` (nicht Teil der Tests).
> Modell: `bindoptix.particles` (Langevin, Sedimentation, Strömung im Rechteckkanal, Gleiten am
> Boden, Einfang an Spots). Geprüft durch T7–T11, T27, T28, T31. Parameter mit **P** sind
> Platzhalter.

## Aufbau

Ein Kanal 200 × 100 µm, 10 mm lang, Reihe aus 5 Spots (200 × 200 µm über die volle Breite,
Rastermaß 800 µm, erster Spot bei x = 3,4 mm). 1 500 nominale Legionellen (0,5 × 2,5 µm,
Δρ = 80 kg/m³ **P**, nicht motil) werden zur Zeit 0 flussgewichtet am Einlass gestartet und 6 h
verfolgt. Fangrate auf den Spots unendlich (ideales Aptamer). Zwei Grenzfälle am Boden:
„bleibt liegen“ (unspezifische Haftung sofort) und „gleitet“ (keine Haftung außerhalb der Spots,
Gleitgeschwindigkeit u_wall = s · a · ∂u/∂z mit s = 0,5 **P**).

## Ergebnis

| Q je Kanal | mittlere Geschwindigkeit | u_wall | Anteil auf Spots, liegen | Anteil auf Spots, gleiten |
| --- | --- | --- | --- | --- |
| 0,5 nL/min | 0,42 µm/s | 0,007 µm/s | 0,9 % | 1,3 % |
| 1,5 nL/min | 1,25 µm/s | 0,022 µm/s | 7,5 % | 24 % |
| 5 nL/min | 4,2 µm/s | 0,072 µm/s | 6,4 % | 28 % |
| 15 nL/min | 12,5 µm/s | 0,22 µm/s | 1,9 % | 18 % |

Statistischer Fehler bei 1 500 Zellen: ±0,5 bis ±1,2 Prozentpunkte.

## Befunde

1. **Gleiten ist langsam.** Eine abgesunkene Zelle gleitet mit etwa 6 s a/h ≈ 1,7 % der mittleren
   Strömungsgeschwindigkeit (T31). Bei 5 nL/min braucht sie für 1 mm fast 4 h. Gleiten hilft
   deshalb nur Zellen, die kurz vor den Spots landen; die ganze Kanallänge wirkt nicht als
   Absetzstrecke. Die Abschätzung (Faktor 10) war zu optimistisch: Hier sind es Faktor 3 bis 9.
2. **Es gibt eine beste Flussrate.** Zu langsam: Die Zellen landen vor dem Gesichtsfeld. Zu
   schnell: Sie verlassen den Kanal, bevor sie absinken. Für diesen Kanal liegt das Optimum bei
   etwa 2 bis 5 nL/min; dort erreichen höchstens rund 30 % der Zellen einen Spot.
3. **Die Lage der Spots im Kanal ist ein Designparameter.** Die Landeorte verteilen sich über
   Millimeter; Spots weiter stromab oder ein Kanal, der vor dem Feld flacher wird, verschieben das
   Optimum. Das ist ein Fall für die Optimierung (Z10).
4. **Unbekannte Wandparameter** (Gleitfaktor, unspezifische Haftung, Fangrate je Kontakt) ändern
   das Ergebnis um mehr als einen Faktor 3. Sie gehören auf die Liste der ersten Messungen (D10,
   D17).

## Reihenfolge der Spots (D15)

Skript `experiments/z3_spotreihenfolge.py`: ganzer Chip über `bindoptix.chipsim.run_chip`,
schräges Grid (tan θ = 1/5) mit einem Kanal je Gridreihe, 5 nL/min je Kanal, 1 500 Zellen je
Kanal, 6 h, alle Spots fangen sofort (jedes Aptamer passt).

| Boden | eingefangen gesamt | Spot 1 (stromauf) | Spot 2 | Spot 3 | Spot 4 | Spot 5 |
| --- | --- | --- | --- | --- | --- | --- |
| Zellen bleiben liegen | 6,2 % | 1,6 % | 1,1 % | 1,3 % | 1,1 % | 1,1 % |
| Zellen gleiten | 29,9 % | 10,5 % | 5,2 % | 5,0 % | 4,5 % | 4,7 % |

(Anteile an allen gestarteten Zellen, Summe über die fünf Kanäle.)

Gleiten die Zellen, fängt der erste Spot jeder Reihe etwa doppelt so viel wie die folgenden:
Er sammelt alles, was stromauf landet und herangleitet, und schattet die übrigen ab. Für
spezifische Aptamere bedeutet das: Ein kreuzreaktives Aptamer weiter stromauf nimmt Zellen weg,
die ein passendes Aptamer weiter stromab gefunden hätten. Die Reihenfolge der Aptamere in einer
Reihe ist damit ein Designparameter, und die Klassifikation (Z10) muss die Lage berücksichtigen.

## Noch nicht im Modell

Behinderte Diffusion und Reibung nahe der Wand (Brenner, Goldman–Cox–Brenner), Ausrichtung der
Stäbchen in der Scherströmung, Motilität am Boden über das Ablösen hinaus, Pulsform und Dauer der
Probenzugabe, Kontinuum für Proteine (Z4).
