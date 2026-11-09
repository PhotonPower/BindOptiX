# Abschätzung: Probenvolumen und Einfangzeit für Legionellen

> Stand 2026-11-09, Stufe Z0/Z3. Überschlag zur Orientierung, keine Abnahme. Formeln mit Quelle
> (`docs/quellen.md`), Parameter gekennzeichnet: **L** Literaturwert, **P** Platzhalter (Quelle
> offen, AGENTS.md Regel 9). Die Zahlen werden in Z3 durch die Tests T27 und T28 abgesichert.

## Anforderung

Trinkwasserinstallation (TrinkwV, technischer Maßnahmenwert): 100 KBE/100 mL. Verdunstungskühlanlagen
(42. BImSchV): Prüfwerte 100 und 1 000, Maßnahmenwert 10 000 KBE/100 mL. Die Nachweisgrenze im
Ausgangswasser muss also bei **c₀ ≈ 1 KBE/mL** liegen (= 10⁻¹² µm⁻³).

KBE zählt nur kultivierbare Zellen; Aptamere binden auch tote und VBNC-Zellen. Sensorwerte und
KBE sind daher nicht eins zu eins vergleichbar (Validierung, Z9).

## Parameter

| Größe | Wert | Status |
| --- | --- | --- |
| Zellform | gestrecktes Rotationsellipsoid, Breite 0,3–0,9 µm, Länge 2–5 µm (Filamente länger) | P |
| Dichteüberschuss Δρ | 60–100 kg/m³ | P |
| Wasser bei 20 °C | η = 1,002 mPa·s | L |
| Reibung des Ellipsoids | Perrin bzw. Happel–Brenner, Orientierung gemittelt: 1/f = (1/f∥ + 2/f⊥)/3 | L (happel1983) |

Sinkgeschwindigkeit v_s = Δρ V g / f, Diffusion D = k_B T / f (Kugelgrenzfall gegen Stokes und
Stokes–Einstein geprüft, Abweichung < 10⁻⁴).

| Fall | v_s | D | Sinkzeit über 100 µm |
| --- | --- | --- | --- |
| klein, 0,3 × 2 µm, Δρ = 60 | 0,008 µm/s | 0,56 µm²/s | 3,6 h |
| nominal, 0,5 × 2,5 µm, Δρ = 80 | 0,026 µm/s | 0,40 µm²/s | 65 min |
| groß, 0,9 × 5 µm, Δρ = 100 | 0,11 µm/s | 0,21 µm²/s | 16 min |

Sedimentations-Péclet-Zahl über 100 µm: v_s h / D ≈ 1,4 (klein) bis 50 (groß). Sedimentation
dominiert, Diffusion ist bei kleinen Zellen nicht vernachlässigbar.

**Motilität** (P): Mit Geißel und Run-and-Tumble (v = 10–30 µm/s, τ ≈ 1 s, α ≈ 0,3) wäre
D_eff = v²τ/(3(1 − α)) ≈ 50–400 µm²/s, also 100- bis 1000-mal größer als die Brownsche Diffusion.
Ob Zellen in der Probe motil sind, hängt von der Wachstumsphase ab (D9). Unten ist Motilität nicht
berücksichtigt; die Zahlen sind die ungünstige Seite.

## Szenarien (nominale Zelle)

**A. Ohne Vorkonzentrierung, Überströmung.** Sedimentationsfluss auf den Sensor: c₀ v_s ≈
**10⁻⁴ Zellen pro mm² und Stunde**. Nicht nachweisbar.

**B. Vorkonzentrat, Überströmung.** Bei laminarer Strömung über eine Fläche A ist der eingefangene
Anteil min(1, v_s A / Q), unabhängig vom Strömungsprofil und von der Kanalhöhe (Absetzbecken nach
Hazen). Vollständiger Einfang verlangt Q ≤ v_s A:

| Sensorfläche A | Q_max | Zeit für 10 µL | Zeit für 100 µL |
| --- | --- | --- | --- |
| 0,25 mm² | 0,4 nL/min | 440 h | 4 400 h |
| 1 mm² | 1,5 nL/min | 110 h | 1 100 h |
| 10 mm² | 15 nL/min | 11 h | 110 h |

**C. Vorkonzentrat, Stillstand über dem Sensor.** Das ganze Konzentrat muss über der Sensorfläche
stehen: A = V_c / h, Zeit h / v_s.

| Volumen V_c | Höhe h | nötige Fläche A | Zeit |
| --- | --- | --- | --- |
| 10 µL | 50 µm | 200 mm² | 33 min |
| 10 µL | 100 µm | 100 mm² | 65 min |
| 100 µL | 100 µm | 1 000 mm² | 65 min |

Zeiten sind brauchbar, aber die Fläche liegt bei 1–10 cm². Das Gesichtsfeld von SIL und
Relay-Optik ist vermutlich um Größenordnungen kleiner (D14).

**D. Durchströmung des Sensors (perforierter Goldfilm).** Jede Zelle, die durch die Fläche fließt,
wird zurückgehalten; eingefangen werden η · c₀ · V_Probe Zellen (η: Rückhalt, P). Bei 100 mL und
1 KBE/mL sind das bis zu 100 Zellen. Größenordnung des Durchflusses: Lochströmung durch eine dünne
Platte nach Sampson, q = r³ Δp / (3η) je Loch; r = 100 nm, Rastermaß 500 nm, Δp = 1 bar ergibt
rund 0,1 mL/s je mm². Das ist eine Obergrenze ohne Filterkuchen, Stützstruktur und Verstopfung.
Folgen für den Aufbau: Die Probe muss unter dem Film abfließen; das verträgt sich nicht mit einer
massiven Glasplatte mit aufgedruckter SIL und verlangt eine Membran (z. B. auf SiN) und einen
anderen Ort für die Sammeloptik.

**Nachweis:** Bei Poisson-verteilten Zellen braucht man im Mittel λ ≈ 3 eingefangene Zellen für
95 % Wahrscheinlichkeit, mindestens eine zu sehen.

## Folgerungen

1. **Phase 1 (vorkonzentriert, D12):** Das Problem ist die Fläche, nicht die Zeit. Entweder das
   Konzentrat auf wenige µL bringen und eine große Fläche abbilden (Mosaik, Scan), oder die Zellen
   aktiv auf das Gesichtsfeld konzentrieren (Hydrodynamik, Dielektrophorese, Magnetpartikel).
   Damit wird Stufe Z8 (Kräfte) wahrscheinlich nötig.
2. **Phase 2 (Leitungswasser direkt, D12):** Überströmung scheidet aus (Szenario A, B). Möglich ist
   nur Einfang beim Durchströmen (D) oder eine vorgeschaltete Konzentrierung im Gerät.
3. **Motilität** kann den Einfang um Größenordnungen verbessern, wenn die Zellen motil sind und
   beim Kontakt binden. Das muss gemessen werden (D9, D10).
4. Die Optik entscheidet über das Gesichtsfeld und damit über die Einfangfläche. Diese Kopplung
   ist die erste Designfrage für die gemeinsame Optimierung (Z10).
