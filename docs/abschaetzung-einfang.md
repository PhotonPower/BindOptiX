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

## Gesichtsfeld 5 mm (D14)

Gesichtsfeld von SIL und Relay-Optik: Durchmesser 5 mm, also A ≈ 19,6 mm². Stillstand über dem
Sensor (Szenario C), das ganze Konzentrat steht über dem Gesichtsfeld (V_c = A · h):

| Kammerhöhe h | Volumen V_c | Sinkzeit nominal | klein | groß |
| --- | --- | --- | --- | --- |
| 100 µm | 2,0 µL | 1,1 h | 3,6 h | 0,3 h |
| 250 µm | 4,9 µL | 2,7 h | 9,0 h | 0,6 h |
| 500 µm | 9,8 µL | 5,4 h | 18 h | 1,3 h |
| 1 mm | 19,6 µL | 11 h | 36 h | 2,6 h |

Diffusion über die Kammerhöhe (h²/2D) dauert 3 h bei 100 µm und 87 h bei 500 µm; die Sedimentation
dominiert. Überströmung (Szenario B) erlaubt höchstens 30 nL/min (nominal), also 5,5 h für 10 µL.

**Erwartete Zellzahl:** Aus 100 mL bei 1 KBE/mL und Rückhalt η der Vorkonzentrierung (P) kommen
100 · η Zellen ins Gesichtsfeld: 30 bis 70 bei η = 0,3 bis 0,7. Mit λ = 3 für 95 % Nachweis liegt die
Nachweisgrenze rechnerisch bei 0,04 bis 0,1 KBE/mL, also unter der Anforderung.

**Anforderung an die Fehlalarme:** Bei 30 bis 70 echten Zellen auf 19,6 mm² müssen falsch erkannte
Spots (unspezifische Bindung, Partikel, Untergrund) deutlich unter 1 je Messung bleiben, also unter
etwa 0,05 je mm². Das ist eine harte Vorgabe für das Signalmodell (Z6) und für die Optik.

**Sammeleffizienz über das Feld:** Eine SIL ist nur in der Mitte aplanatisch. Über 5 mm ändert sich
die Sammeleffizienz mit dem Ort; die Antwortkarte (Z2) muss das ganze Feld abdecken, und der
Feldverlauf wird ein Designparameter (Raytatouille).

## Chip mit 5 × 5 Spots und Kanälen ≤ 200 µm

Vorgabe des Maintainers (2026-11-09): 5 × 5 Sensorflächen im Gesichtsfeld, Kanäle höchstens
200 µm breit. Geometrie in `bindoptix.chip` (Layout `rows`, D15): je Spotreihe ein gerader Kanal,
Spots 200 × 200 µm, Rastermaß 800 µm, Kanalhöhe 100 µm, Kanallänge 10 mm (alle P). Eine 5 × 5-Matrix
passt nur mit Rastermaß ≤ 833 µm in den Kreis von 5 mm (äußerste Ecke bei √2 · (2p + s/2)); bei
1 mm Rastermaß läge sie außerhalb.

| Größe | Wert |
| --- | --- |
| Spotfläche gesamt | 1,0 mm² |
| Kanalboden im Gesichtsfeld | 4,4 mm² |
| Flüssigkeit über dem Kanalboden im Gesichtsfeld | 0,44 µL |
| Volumen aller Kanäle (10 mm lang) | 1,0 µL |

**Das ändert die Folgerung aus dem Abschnitt „Gesichtsfeld 5 mm“:** Über dem Gesichtsfeld stehen
jetzt nur 0,44 µL statt bis zu 20 µL, und eingefangen wird nur auf 1 mm² Spotfläche. Damit gibt es
zwei Grenzfälle:

| Annahme | wirksame Fläche | Q_max (nominal) | Zeit für 1 µL | für 10 µL |
| --- | --- | --- | --- | --- |
| Zellen bleiben, wo sie absinken | Spots, 1 mm² | 1,5 nL/min | 11 h | 109 h |
| abgesunkene Zellen gleiten am passivierten Boden weiter bis zu einem Spot | Kanalboden, 10 mm² | 15 nL/min | 1,1 h | 11 h |

(klein bis groß: Faktor 3 langsamer bis 4 schneller.)

Ob Legionellen nach dem Absinken am Boden mit der Strömung weitergleiten und an einem passenden
Spot binden (Wandbewegung einer Kugel in Scherströmung: goldman1967), entscheidet über einen
Faktor 10 und ist die wichtigste offene Frage für Z3 und für die Messungen (D10). Gleiten die Zellen,
wirkt die ganze Kanallänge stromauf als Absetzstrecke; das Konzentrat sollte dann bei rund 1 µL
liegen (Faktor 10⁵ aus 100 mL).

**Klassifikation:** Im Layout `rows` trifft eine Zelle höchstens die fünf Spots ihrer Reihe. Für
25 Subtypen müsste jede Zelle alle Aptamere passieren können: ein Mäander durch alle 25 Spots,
oder dieselben Aptamere in jeder Reihe in anderer Reihenfolge (D15).

> **Nachtrag aus Z3 (`docs/ergebnisse-z3.md`):** Gleiten ist mit etwa 1,7 % der mittleren
> Strömungsgeschwindigkeit langsam. Statt Faktor 10 bringt es im Modell Faktor 3 bis 9, und pro Kanal
> erreichen bei der besten Flussrate (2–5 nL/min) höchstens rund 30 % der Zellen einen Spot.

## Folgerungen

1. **Phase 1 (vorkonzentriert, D12):** Mit 5 × 5 Spots in Kanälen ≤ 200 µm braucht es ein
   Konzentrat von rund 1 µL (Faktor 10⁵ aus 100 mL) und Durchfluss von etwa 15 nL/min, wenn
   abgesunkene Zellen am Boden bis zu den Spots gleiten; sonst zehnmal langsamer. Ob aktive Kräfte
   (Z8) nötig sind, hängt an dieser Frage. (Ältere Fassung ohne Kanäle: Konzentrat ≤ 10 µL.)
2. **Phase 2 (Leitungswasser direkt, D12):** Überströmung scheidet aus (Szenario A, B). Möglich ist
   nur Einfang beim Durchströmen (D) oder eine vorgeschaltete Konzentrierung im Gerät.
3. **Motilität** kann den Einfang um Größenordnungen verbessern, wenn die Zellen motil sind und
   beim Kontakt binden. Das muss gemessen werden (D9, D10).
4. Die Optik entscheidet über das Gesichtsfeld und damit über die Einfangfläche. Diese Kopplung
   ist die erste Designfrage für die gemeinsame Optimierung (Z10).
