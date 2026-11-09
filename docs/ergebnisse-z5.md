# Erste Ergebnisse Stufe Z5: Fouling und Einfang

> Stand 2026-11-09. Modell: `chipsim.run_channel_coupled` (Kontinuum Z4 → freie Plätze θ(t) je
> Spot → Fangrate k_cap = k_bond · θ · A_contact der Partikel Z3), geprüft durch T36. Alle
> Parameter sind Platzhalter (**P**); die Zahlen zeigen das Verhalten, nicht den Sensor.

## Aufbau

Mittlerer Kanal des schrägen 5 × 5-Grids, 5 nL/min, Γ_max = 10⁴ /µm², k_bond so, dass ein
freier Spot eine aufliegende Zelle im Mittel nach 1 s fängt. Störprotein wie Albumin
(D = 60 µm²/s) mit 10 /µm³ (≈ 17 nM) im Zulauf, bindet irreversibel und unspezifisch mit
k_on = 1 µm³/s (**P**).

## Ergebnis

Freie Plätze je Spot (Anteil von Γ_max), Spot 1 stromauf:

| Zeit | Spot 1 | Spot 2 | Spot 3 | Spot 4 | Spot 5 |
| --- | --- | --- | --- | --- | --- |
| 5 min | 97 % | 99 % | 100 % | 100 % | 100 % |
| 15 min | 72 % | 88 % | 95 % | 98 % | 99 % |
| 25 min | 25 % | 63 % | 84 % | 94 % | 98 % |
| 3 h | 0 % | 0 % | 0 % | 0 % | 0 % |

Fänge in 3 h aus 300 Zellen: 68 ohne, 25 mit Störprotein.

## Befunde

1. **Die Füllfront läuft stromab.** Der erste Spot einer Reihe verliert seine Plätze zuerst.
   Zusammen mit der Abschattung aus Z3 (der erste Spot fängt die meisten Zellen) ist er doppelt
   empfindlich; ein Referenz- oder Opferspot stromauf könnte die übrigen schützen.
2. **Der Wettlauf ist eine Frage der Zeitskalen.** Zellen brauchen rund eine Stunde, um zu
   sinken und die Spots zu erreichen; ein Störprotein mit nanomolarer Konzentration belegt die
   Spots transportbegrenzt in einer halben Stunde. In Phase 1 (gereinigtes Konzentrat im
   Aptamer-Puffer) ist das unkritisch, für Leitungswasser und Proben mit Proteinen (Phase 2)
   entscheidet die Passivierung (Messplan M4).
3. Die Messung der unspezifischen Bindung von Proteinen an die Spots gehört deshalb zu M4.
