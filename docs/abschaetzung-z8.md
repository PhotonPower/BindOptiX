# Abschätzung Stufe Z8: Erwärmung und optische Kräfte

> Stand 2026-11-09. Einstiegskriterium der Roadmap: Z8 wird gebaut, wenn Erwärmung oder Kräfte
> eine Zielgröße um mehr als 5 % ändern. Annahmen wie Z6 (10 mW auf Ø 5 mm, 509 W/m²); alle
> Werte **P**. Ergebnis: **Z8 entfällt** für diese Anregung.

## Erwärmung

Absorbiert der dünne Goldfilm 30–50 % der Anregung, fließen q = 150–250 W/m² gleichmäßig über
die Kreisfläche (Radius R = 2,5 mm) in Glas und Wasser. Für eine gleichmäßig beheizte Kreisfläche
auf einem Halbraum ist die Erwärmung in der Mitte ΔT = q R / k (carslaw1959); mit Glas und Wasser
auf beiden Seiten k ≈ 1,0 + 0,6 W/(m·K):

| Absorption | q | ΔT in der Mitte | Änderung k_on (E_a = 50 kJ/mol) | Änderung Viskosität |
| --- | --- | --- | --- | --- |
| 30 % | 153 W/m² | 0,24 K | 1,7 % | 0,6 % |
| 50 % | 255 W/m² | 0,40 K | 2,8 % | 1,0 % |

Die SIL unter dem Glas leitet zusätzlich Wärme ab; die Werte sind eine obere Schranke. Lokale
Überhöhungen in plasmonischen Kavitäten bleiben bei 509 W/m² im Millikelvin-Bereich, weil die
absorbierte Leistung je Kavität winzig ist. 5 % Änderung der Kinetik verlangen ΔT ≈ 0,7 K, also
etwa **18 mW** Anregung auf dem Feld (bei 50 % Absorption).

## Optische Kräfte

Bei gleichmäßiger Beleuchtung gibt es keine Gradientenkraft. Der Strahlungsdruck auf eine Zelle,
F = n I C_pr / c mit C_pr ≈ 0,1 × geometrischer Querschnitt (schwach streuende Zelle), beträgt
etwa 3·10⁻¹⁹ N, das Gewicht im Wasser 2,6·10⁻¹⁶ N: Verhältnis 10⁻³. Feldüberhöhungen an Kavitäten
wirken nur auf Nanometer-Volumina; BioOptoMD nennt für Moleküle erst ab etwa 10³-facher Überhöhung
merkliche Kräfte. Thermophorese fehlt mangels Temperaturgradient.

## Folgerung

Z8 entfällt bei der jetzigen Anregung. Wiedereinstieg, wenn die Anregung über etwa 15 mW auf dem
Feld steigt, fokussiert oder gescannt wird, oder wenn die Kavitäten stark absorbieren (dann mit
der Absorptionsdichte q aus der Feldkarte F, hp-FEM M9).
