# Messplan: Parameter für Transport und Einfang

> Stand 2026-11-09. Ziel: die Platzhalter aus Z3 (D9, D10, D17) mit wenigen Versuchen festlegen.
> Reihenfolge nach Einfluss auf das Ergebnis (`docs/ergebnisse-z3.md`): Wandparameter ändern den
> Einfang um mehr als Faktor 3, die Zellgröße um Faktor 3 bis 4, Wandhydrodynamik kaum.
> *Legionella pneumophila* ist Risikogruppe 2: Arbeiten im S2-Labor.

## Gemeinsamer Aufbau

Ein Kanal 200 × 100 µm wie im Chip (oder breiter für bessere Statistik), Boden aus demselben
Schichtstapel (Gold mit Passivierung, ohne bzw. mit Aptamer-Spots), Spritzenpumpe 0,5–20 nL/min,
Mikroskop mit 10×–40× Objektiv von unten, Zellen angefärbt (Nukleinsäure-Farbstoff) oder ein
fluoreszierender Stamm. Bildfolgen 1–20 Bilder/s, Verfolgung der Zellen in 2D, Fokuslage für z.
Die Wandscherrate folgt aus Q und dem Kanal (`bindoptix.flow.DuctFlow.floor_shear_rate_s`).

## Versuche

| Nr. | Größe im Modell | Versuch | Auswertung | Umfang |
| --- | --- | --- | --- | --- |
| M1 | Zellmaße, Δρ (`CellModel`) | ruhende Kammer, Zellen sinken aus homogener Verteilung | Länge und Breite aus Bildern; Ankunftsrate am Boden über der Zeit gegen das Modell (Sinken + Diffusion, T8/T10) → v_s, daraus Δρ | 200 Zellen je Stamm |
| M2 | Motilität: Anteil, v, τ, α (D9) | ruhende Kammer, Bahnverfolgung im Volumen nahe dem Boden | Anteil schwimmender Zellen, Laufgeschwindigkeit, Laufzeit, Wendewinkel; je Wachstumsphase und nach der Vorkonzentrierung | 100 Bahnen |
| M3 | Gleitfaktor s (D17) | abgesunkene Zellen auf passiviertem Boden, Wandscherung 0,05–1 s⁻¹ | Geschwindigkeit am Boden gegen a · ∂u/∂z → s; prüfen, ob Zellen gleichmäßig gleiten oder ruckweise | 50 Zellen je Scherrate |
| M4 | unspezifische Haftung k_ns (D17) | wie M3, über 30–60 min | Anteil der Zellen, die liegen bleiben, über der Zeit → Rate; Vergleich verschiedener Passivierungen | 200 Zellen je Passivierung |
| M5 | Fangrate k_cap je Spot, Kreuzreaktivität (D3, D13, D17) | gleitende Zellen überqueren einen Aptamer-Spot; passender und nicht passender Stamm | Fangwahrscheinlichkeit P = 1 − exp(−k_cap · T), T = Spotlänge / Gleitgeschwindigkeit; T über die Scherrate variieren → k_cap und seine Scherabhängigkeit | 100 Überquerungen je Paar |
| M6 | Rückhalt η der Vorkonzentrierung (D12) | Membranfiltration dotierter Proben, Resuspension auf ≈ 1 µL | gefundene / eingesetzte Zellen (Kultur, Mikroskopie oder qPCR) | 5 Proben je Konzentration |

Größenordnungen aus dem Modell: Bei 5 nL/min ist die Wandscherung in Kanalmitte etwa 0,34 s⁻¹
und eine gleitende Zelle mit s = 0,5 bewegt sich mit rund 0,07 µm/s; in 10 min sind das 40 µm,
gut messbar. Ein Spot von 200 µm wird in etwa 46 min überquert; damit ist M5 für k_cap zwischen
etwa 10⁻⁴ und 10⁻¹ s⁻¹ empfindlich.

## Rückfluss ins Modell

Jede gemessene Größe ersetzt einen Platzhalter mit Kennzeichnung „Kalibrierung“ (AGENTS.md
Regel 9) und Unsicherheit; danach werden `experiments/z3_einfang_flussrate.py` und
`experiments/z3_spotreihenfolge.py` neu gerechnet. M5 liefert die ersten Einträge der
Bindungsmatrix für die Klassifikation (Z10).
