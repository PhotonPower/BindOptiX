# Plan und Stufen

> Erstfassung 2026-11-09. Z0 abgenommen am 2026-11-09 (PR #1). Entscheidungen am Ende: D1–D5 und
> D12–D14 entschieden, D15, D16 neu offen, die übrigen offen. Abschätzung zum Einfang: `docs/abschaetzung-einfang.md`.

## Rolle im Framework

BindOptiX ist die **Zeit-, Transport- und Biologieschicht** über der Fusion (RTT-Fusion,
Raytatouille). Es rechnet kein Licht. Die Optik kommt als vorab gerechnete Karten aus den
vorhandenen Lösern, BindOptiX verbindet sie mit dem zeitlichen Ablauf eines Assays:

| Aufgabe | Wer rechnet | Was BindOptiX daraus macht |
| --- | --- | --- |
| Rezeptorschicht, Fluorophor-Dipole je Bindungszustand | BioOptoMD (Q, N) | gewichtet die Zustände mit der Belegung (Port B) |
| Anregungsfeld, Absorption, Erwärmung | hp-FEM, CliffordBem (Port F) | Anregungsrate, Bleichen, T-abhängige Kinetik |
| Emission, Quenching, Sammlung durch SIL und Relay-Optik | K5-Kette: CliffordBem/hp-FEM → `rtt-field` → Raytatouille (Antwortkarte, Darstellung M) | Photonenrate je Pixel und Spot |
| Strömung, Transport, Einfang, Bindung, Assay-Ablauf | **BindOptiX** | Belegung B(t), Ereignisliste |
| Kamera, Zählung, Klassifikation, Zielfunktion | **BindOptiX** | Nachweisgrenze, Zeit bis zum Ergebnis, Konfusionsmatrix |

Grundlage ist die Trennung der Zeitskalen: Optik zeitharmonisch (fs), Fluoreszenzzyklus (ns),
Transport und Bindung (s bis h). Die Optik wird je Design einmal gerechnet, nicht je Zeitschritt.
Die Kopplung ist einseitig und damit ein azyklischer Ausführungsgraph
(F → Kinetik → B → Q-Gewichte → Signal); Ausnahmen in ADR-Entwurf Datenports, Abschnitt „Grenzen“.

## Zielgrößen

Anwendung: Legionellen in Wasser (D2). Anforderung aus dem Regelwerk: technischer Maßnahmenwert
der TrinkwV 100 KBE/100 mL; 42. BImSchV Prüfwerte 100 und 1 000, Maßnahmenwert 10 000 KBE/100 mL.
Ziel der Nachweisgrenze im Ausgangswasser daher **1 KBE/mL**.

| Größe | Definition | Stufe |
| --- | --- | --- |
| Nachweisgrenze | kleinste Konzentration (KBE/mL), bei der ein Stamm mit vorgegebener Sicherheit erkannt wird | Z8, Z10 |
| Zeit bis zum Ergebnis | Zeit vom Einspülen bis zur Entscheidung | Z7, Z8 |
| Stammerkennung | Konfusionsmatrix über die Stämme aus D2, inklusive „kein Bakterium“ | Z10 |
| SNR je gebundenem Bakterium | Photonen des Spots gegen Untergrund und Rauschen in einem Kamerabild | Z6 |
| Fehlalarmrate | falsch erkannte Spots je mm² und Stunde (unspezifische Bindung, Untergrund) | Z6, Z8 |

## Stufen

| Stufe | Inhalt | Abnahme |
| --- | --- | --- |
| **Z0 Spezifikation** ✅ | Zielgrößen, Konventionen (`docs/konventionen.md`), Datenformate, Entscheidungen D1–D11; Inventareintrag im Fusionsplan; ADR-Entwurf Datenports (F, B, Antwortkarte); Volltext der Kernquellen | Spezifikation und ADR-Entwurf freigegeben |
| **Z1 Sensorgeometrie und Bezugssysteme** ✅ (Kavitäten folgen mit D4) | Chip über Nurbforge (Kanal, Spots, Kavitäten, SIL als Verweis), Sensorsystem (Oberfläche z = 0, Lösung z > 0), Pose zum RTT-System mit Spiegelregel, Lage der Schwerkraft (Einbaulage) | T1, T2, T29, T30 |
| **Z2 Optische Antwortkarten** ◐ (Datenmodell, Interpolation, Format, T6 ✅; Export aus hp-FEM/CliffordBem vorgeschlagen) | Import der Karten aus dem Ausführungsgraphen der Fusion (Pilot K5): Anregungsfeld F, Zerfalls- und Detektionstensoren je Emitterort; Interpolation zwischen Stützstellen; nichts davon wird hier gerechnet | T3–T6 |
| **Z3 Partikel (Bakterien)** ✅ (Zellmodell, Kanalströmung, Langevin mit Brücke, Gleiten, Spots, Chiplayout, Wandhydrodynamik; Ausrichtung der Stäbchen orientierungsgemittelt, Verfeinerung bei Bedarf) | Stäbchen als Ellipsoid (Reibung nach Perrin), Langevin-Dynamik mit Drift aus Strömung und Sedimentation, Motilität (Run-and-Tumble), sterische Wände, Einfang mit Wahrscheinlichkeit je Wandkontakt; Strömung zunächst analytisch (D7) | T7–T11, T27, T28, T31–T34 |
| **Z4 Kontinuum (Proteine, Störstoffe)** ✅ (2D-Schnitt und 3D-Kanalvolumen, mehrere Spezies, kompetitive Langmuir-Kinetik, Spots aus dem Chiplayout; Kopplung an die Partikel in Z5) | Konvektion–Diffusion im Kanal, mehrere Spezies, Oberflächenreaktion | T12–T14, T35 |
| **Z5 Bindung und Rezeptorzustände** | kompetitive Langmuir-Kinetik, irreversibler Einfang mit Avidität für Bakterien, unspezifische Bindung; Rezeptorzustände (frei, Ziel gebunden, Störstoff gebunden) als Port B; Modelle als Template-Parameter (kein NVRTC) | T15, T16 |
| **Z6 Signalmodell** | Anregungsrate aus F, Quantenausbeute und Sammlung aus der Antwortkarte, je Stichprobe des Q-Ensembles (nicht über gemittelte Tensoren, siehe ADR-Entwurf); Bleichen; Untergrund (Gold-Photolumineszenz, Autofluoreszenz von 2PP-Lack, Glas und Bakterien) als Parameter oder aus Scattuccino; Kameramodell (Schrot-, Lese-, Dunkelrauschen) | T17–T20 |
| **Z7 Assay-Protokoll und digitale Auslese** | Abläufe (Einspülen, Inkubation, Waschen), synthetische Kamerabilder, Spot-Erkennung, Zählung je Spot und Zeit | T21, T22 |
| **Z8 Thermik und Kräfte (nur nach Abschätzung)** | ΔT aus F (hp-FEM M9) → k(T), D(T), Thermophorese; Kraftkarten für Bakterien aus CliffordBem (Spannungstensor) als Drift. Einstieg nur, wenn die Abschätzung in Z6 einen Einfluss > 5 % auf eine Zielgröße zeigt | T23, T24 |
| **Z9 Kalibrierung** | Einfangwahrscheinlichkeit, Kinetikparameter, Untergrund aus Messdaten (D10) mit Unsicherheiten | Vorhersage an nicht zur Kalibrierung genutzten Messungen |
| **Z10 Klassifikation und Optimierung** | Kreuzreaktivitätsmatrix des Panels, synthetische Trainingsdaten, Konfusionsmatrix; Designoptimierung über Parameter der Nurbforge-Geometrie und des Schichtstapels (BoTorch wie hp-FEM, D11) | T25, T26 |
| **Z11 C++/CUDA-Kern** | Partikel und Kontinuum auf der GPU (RTX 6000 Ada, RTX 3090; FP32 mit sorgfältiger Akkumulation wie BioOptoMD ADR 0002) | Sollwertdateien aus der Python-Referenz, Laufzeit |
| **Z12 Fusion und Pilot** | Adapter für RTT-Fusion; Pilot K8 „Biosensor über der Zeit“: BioOptoMD + hp-FEM/CliffordBem + Raytatouille + BindOptiX, mit gedruckter SIL aus Voxelspaetzle (V10) | Konformitätssatz der Fusion; Vergleich mit Messung |

Reihenfolge: Z2 vor Z3, weil die optischen Karten zeigen, welche Bereiche der Oberfläche
überhaupt Signal liefern; das begrenzt Transport und Bindung auf die relevanten Flächen. Z3 vor
Z4, weil Bakterien bei realistischen Konzentrationen einzelne Teilchen sind (10³ KBE/mL ≈ 10 Zellen
in 10 µL), kein Konzentrationsfeld.

## Entscheidungen

| Nr. | Frage | Entscheidung bzw. Vorschlag | Stand |
| --- | --- | --- | --- |
| D1 | Name | BindOptiX überall; GitHub-Beschreibung angepasst | entschieden 2026-11-09 |
| D2 | Zielorganismen, Probe, Konzentration | Legionellen; zuerst vorkonzentrierte Proben im Aptamer-Puffer, danach Leitungswasser direkt (D12); Zielkonzentration 1 KBE/mL im Ausgangswasser (Regelwerk) | entschieden 2026-11-09 |
| D3 | Panel und Signal | Aptamere aus der Literatur, Turn-on; 25 Subtypen gestuft (D13). Literatur bisher nur für L. pneumophila bzw. Serogruppe 1, Kreuzreaktivität zwischen Subtypen unbekannt | entschieden 2026-11-09 |
| D4 | Schichtstapel und Anregung | zuerst Gold (Thiol-Anker), danach Graphit; Transmissionsmessung; Wellenlänge folgt aus dem Fluorophor, sichtbar (Vorschlag: rot, ≥ 600 nm). Offen: Golddicke oder perforierter Film; Graphit nur als dünne Lage transmissionstauglich | teilweise entschieden |
| D5 | Fluidik | Durchfluss, Kanal < 1 mm, horizontal. Offen: Sensorfläche als Kanalboden (Annahme) | teilweise entschieden |
| D6 | Auslese | digitales Zählen je Spot | Vorschlag |
| D7 | Strömung | analytisch (Poiseuille, Hele-Shaw) bis Z7; danach Erweiterungsvorschlag an hp-FEM (Stokes, Taylor-Hood) unter `docs/vorschlaege/` oder kleiner eigener Löser; keine fremden CFD-Pakete (`docs/abhaengigkeiten.md`) | Vorschlag |
| D8 | Sprache und Kern | Python-Referenz (Z1–Z10), C++/CUDA-Kern danach (Z11) | Vorschlag |
| D9 | Motilität | Run-and-Tumble als Parameter, lebend und tot getrennt; Motilität von Legionellen in der Probe messen | Vorschlag |
| D10 | Messdaten für die Kalibrierung | offen; bis dahin Literaturwerte mit Kennzeichnung | offen |
| D11 | Optimierer | BoTorch (wie `hpfem.opt`), Mehrziel über die Zielgrößen | Vorschlag |
| D12 | Probenvorbereitung | Phase 1: Vorkonzentrierung im Labor (z. B. Membranfiltration, dabei Wechsel in den Aptamer-Puffer); Phase 2: Leitungswasser direkt. Phase 2 verlangt Einfang beim Durchströmen oder Konzentrierung im Gerät (Abschätzung, Szenario D) | entschieden 2026-11-09 |
| D13 | Klassen | gestuft: Legionella spp. gegen Rest, L. pneumophila, Serogruppe 1, dann weitere Subtypen bis 25; die Simulation liefert die Anforderungen an das Panel | Vorschlag |
| D14 | Gesichtsfeld von SIL und Relay-Optik gegen nötige Einfangfläche | Gesichtsfeld Ø 5 mm (19,6 mm²) mit 5 × 5 Spots, Kanäle ≤ 200 µm (2026-11-09). Rastermaß ≤ 833 µm, damit die Matrix ins Feld passt. Einfang siehe `docs/abschaetzung-einfang.md`: Konzentrat ≈ 1 µL; Faktor 10 hängt daran, ob abgesunkene Zellen bis zu den Spots gleiten | entschieden 2026-11-09 |
| D15 | Kanalführung | Maintainer 2026-11-09: schräg stehendes 5 × 5-Grid, endgültiges Layout offen. Die Geometrie erlaubt jeden Kippwinkel und frei liegende gerade Kanäle und meldet je Kanal die überdeckten Spots in Strömungsrichtung. Bei tan θ = 1/5 liegt jeder Spot quer zur Strömung an eigener Stelle (Abstand p/√26 ≈ 157 µm bei p = 800 µm, also kleiner als ein 200-µm-Spot) | teilweise entschieden |
| D16 | Maße | Spotgröße, Rastermaß, Kanalhöhe und -länge, Substratdicke; bis dahin Platzhalter 200 µm, 800 µm, 100 µm, 10 mm | offen |
| D17 | Wandparameter: Gleitfaktor am Boden, unspezifische Haftung, Fangrate je Kontakt, Kontaktabstand | Platzhalter (0,5; 0 bzw. ∞; ∞; 20 nm). Sie ändern den Einfang um mehr als Faktor 3 (`docs/ergebnisse-z3.md`); Messplan `docs/messplan-wandparameter.md` (M3–M5) | offen |
