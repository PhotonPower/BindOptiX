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

## Noch nicht im Modell

Behinderte Diffusion und Reibung nahe der Wand (Brenner, Goldman–Cox–Brenner), Ausrichtung der
Stäbchen in der Scherströmung, Motilität am Boden über das Ablösen hinaus, Pulsform und Dauer der
Probenzugabe, Kontinuum für Proteine (Z4).
