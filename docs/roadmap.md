# Plan und Stufen

> Erstfassung 2026-11-09. Entscheidungen des Maintainers (offen: D1–D11) stehen am Ende.

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
| **Z0 Spezifikation** | Zielgrößen, Konventionen (`docs/konventionen.md`), Datenformate, Entscheidungen D1–D11; Inventareintrag im Fusionsplan; ADR-Entwurf Datenports (F, B, Antwortkarte); Volltext der Kernquellen | Spezifikation und ADR-Entwurf freigegeben |
| **Z1 Sensorgeometrie und Bezugssysteme** | Chip über Nurbforge (Kanal, Spots, Kavitäten, SIL als Verweis), Sensorsystem (Oberfläche z = 0, Lösung z > 0), Pose zum RTT-System mit Spiegelregel, Lage der Schwerkraft (Einbaulage) | T1, T2 |
| **Z2 Optische Antwortkarten** | Import der Karten aus dem Ausführungsgraphen der Fusion (Pilot K5): Anregungsfeld F, Zerfalls- und Detektionstensoren je Emitterort; Interpolation zwischen Stützstellen; nichts davon wird hier gerechnet | T3–T6 |
| **Z3 Partikel (Bakterien)** | Langevin-Dynamik mit Drift aus Strömung und Sedimentation, Motilität (Run-and-Tumble), sterische Wände, Einfang mit Wahrscheinlichkeit je Wandkontakt; Strömung zunächst analytisch (D7) | T7–T11 |
| **Z4 Kontinuum (Proteine, Störstoffe)** | Konvektion–Diffusion im Kanal, mehrere Spezies, Oberflächenreaktion | T12–T14 |
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

## Entscheidungen (offen, Issues mit Label `frage`)

| Nr. | Frage | Vorschlag |
| --- | --- | --- |
| D1 | Name: Repository „BindOptiX“, GitHub-Beschreibung „PhotoKinX“ | BindOptiX überall; Beschreibung anpassen |
| D2 | Zielorganismen und Stämme, Probenmatrix, Konzentrationsbereich | offen – bestimmt Größe, Motilität, Störstoffe, Klassen |
| D3 | Aptamer-Panel, Multiplex (räumlich getrennte Spots oder Farben), Fluorophore, Signalprinzip (Turn-on durch Konformation oder Ablösung, Turn-off) | offen; bis dahin ein Aptamer, ein Farbstoff, Turn-on |
| D4 | Schichtstapel und Anregung: Golddicke, Graphen ja/nein, Kavitäten, Wellenlänge, Seite, Winkel, Leistung | offen; Referenz: Glas / Au 50 nm / Wasser wie K4 |
| D5 | Fluidik: Durchfluss oder stehender Tropfen, Kanalmaße, Einbaulage | offen – bestimmt Sedimentation und Einfang |
| D6 | Auslese: Kamera, Belichtung, digitales Zählen oder integriertes Signal | digitales Zählen je Spot |
| D7 | Strömung | analytisch (Poiseuille, Hele-Shaw) bis Z7; danach Erweiterungsvorschlag an hp-FEM (Stokes, Taylor-Hood) unter `docs/vorschlaege/` oder kleiner eigener Löser; keine fremden CFD-Pakete (`docs/abhaengigkeiten.md`) |
| D8 | Sprache und Kern | Python-Referenz (Z1–Z10), C++/CUDA-Kern danach (Z11), wie BWPM und Voxelspaetzle; README entsprechend anpassen |
| D9 | Motilität: lebende und tote Zellen getrennt modellieren | ja, Run-and-Tumble als Parameter je Stamm |
| D10 | Messdaten für die Kalibrierung (Bindungskurven, Einfang je Kontakt, Untergrund) | offen; bis dahin Literaturwerte mit Kennzeichnung |
| D11 | Optimierer | BoTorch (wie `hpfem.opt`), Mehrziel über die Zielgrößen |
