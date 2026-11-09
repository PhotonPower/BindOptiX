# Testfälle

Status: ✅ umgesetzt und grün, ◐ teilweise, ○ geplant. Sollwerte nie aus dem geprüften Code.

| Nr. | Prüfung | Referenz | Stufe | Status |
| --- | --- | --- | --- | --- |
| T1 | Pose Sensorsystem ↔ RTT-System, Rundreise (Punkte, polare und axiale Vektoren, Tensoren) | Identität bis 1e−12; Drehungen nach Rodrigues im Test | Z1 | ✅ `tests/test_frames.py` |
| T2 | Spiegelregel z → −z (Beleuchtung von der Lösungsseite): polare und axiale Vektoren, Tensoren | Befund K4 (RTT-Fusion); Kreuzprodukt ist axial, T(a × b) = T(a) × T(b); Invarianz von μᵀAμ | Z1 | ✅ `tests/test_frames.py` |
| T3 | Anregungsfeld an einem Goldfilm (Kretschmann, p-Pol.) aus F | Transfer-Matrix, Resonanzwinkel wie K4 | Z2 | ○ |
| T4 | Zerfallsratentensor eines Dipols über Glas / Au / Wasser über dem Abstand | Chance–Prock–Silbey bzw. Lukosz–Kunz (unabhängige Python-Rechnung im Test) | Z2 | ○ |
| T5 | Detektionstensor: Abstrahlung in das Substrat, überkritischer Anteil, SPCE-Winkel = SPR-Winkel | Lukosz–Kunz; Lakowicz 2004 | Z2 | ○ |
| T6 | Eigenschaften der Tensoren: symmetrisch, positiv semidefinit, A_det ≤ A_rad ≤ A_tot; Mittel über Orientierungen = Spur/3; Interpolation erhält beides | analytisch; Ikosaeder als sphärisches 5-Design | Z2 | ✅ `tests/test_ports.py` |
| T7 | Freie Diffusion: mittlere quadratische Verschiebung 6Dt | Einstein, Standardfehler 4σ | Z3 | ○ |
| T8 | Sedimentation: Driftgeschwindigkeit und Gleichgewicht c(z) ∝ exp(−z·v_s/D) | Stokes; Mason–Weaver | Z3 | ○ |
| T9 | Run-and-Tumble: effektive Diffusion v²τ/(3(1−α)) | Lovely–Dahlquist 1975 | Z3 | ○ |
| T10 | Einfang an absorbierender Wand aus homogener Anfangsverteilung: Γ(t) = 2c₀√(Dt/π) | Smoluchowski | Z3 | ○ |
| T11 | Determinismus: gleicher Seed, gleiche Ergebnisse, unabhängig von der Threadzahl | Rundreise | Z3 | ○ |
| T12 | Transportbegrenzter Einfang im Kanal | Lévêque-Lösung (Graetz-Problem) | Z4 | ○ |
| T13 | Übergang reaktions- zu transportbegrenzt über die Damköhler-Zahl | Squires, Messinger, Manalis 2008 | Z4 | ○ |
| T14 | Massenerhaltung (Volumen + Oberfläche + Auslass) | Bilanz bis 1e−10 relativ | Z4 | ○ |
| T15 | Langmuir und kompetitive Langmuir-Kinetik, Gleichgewicht und Zeitverlauf | analytische Lösungen | Z5 | ○ |
| T16 | Rezeptorzustände: Anteile summieren sich zu 1; ohne Analyt = Zustand „frei“ | Grenzfall | Z5 | ○ |
| T17 | Einzelner Emitter im homogenen Medium: Rate = Anregung · q₀ · Sammeleffizienz | analytisch | Z6 | ○ |
| T18 | Modifizierte Quantenausbeute aus den Tensoren | Formel mit γ_nr⁰ = γ_r⁰(1 − q₀)/q₀, unabhängige Rechnung | Z6 | ○ |
| T19 | Bleichen: Monoexponentiell bei konstanter Anregung | analytisch | Z6 | ○ |
| T20 | Kamera: Mittel und Varianz der Zählrate | Poisson plus Gauß-Leserauschen | Z6 | ○ |
| T21 | Spot-Erkennung auf synthetischen Bildern mit bekannter Wahrheit | Trefferquote, Fehlalarme gegen Vorgabe | Z7 | ○ |
| T22 | Protokoll: Waschschritt entfernt reversibel Gebundenes nach k_off | analytisch | Z7 | ○ |
| T23 | Erwärmung: Bilanz ∫q dV = absorbierte Leistung aus F | Energieerhaltung | Z8 | ○ |
| T24 | Optische Kraft auf eine Kugel aus der Kraftkarte | Mie (CliffordBem-Sollwerte) | Z8 | ○ |
| T25 | Klassifikation: Spielproblem mit bekannter Konfusionsmatrix | analytisch | Z10 | ○ |
| T26 | Optimierer findet das bekannte Optimum eines Spielproblems | analytisch | Z10 | ○ |
| T27 | Sinkgeschwindigkeit und Diffusion eines Stäbchens (Rotationsellipsoid, Orientierung gemittelt) | Perrin, Happel–Brenner; Kugelgrenzfall Stokes und Stokes–Einstein | Z3 | ○ |
| T28 | Einfang im laminaren Durchfluss über einer Fläche: Anteil min(1, v_s A/Q), unabhängig vom Profil | Hazen, analytisch | Z3 | ○ |
| T29 | Chipgeometrie: Volumen der Kanäle (Nurbforge) gegen l·b·h, Spotflächen, Kanalboden im Gesichtsfeld gegen numerische Quadratur der Sehnenlänge, Spots in ihrem Kanal, STEP-Rundreise, ungültige Layouts | analytisch, `scipy.integrate.quad` | Z1 | ✅ `tests/test_chip.py` |
| T30 | Schräges Grid und freie Kanäle: Querpositionen bei tan θ = 1/5 alle verschieden mit Abstand p/√26, Passung drehinvariant, Überdeckung Spot/Kanal gegen Monte Carlo, gedrehte Kanalkörper, sich kreuzende Kanäle abgelehnt | analytisch, Monte Carlo mit Seed | Z1 | ✅ `tests/test_chip.py` |

T3–T5 prüfen exportierte Karten aus hp-FEM bzw. CliffordBem und warten auf deren Export
(`docs/vorschlaege/export-antwortkarte.md`). Datenmodell, Interpolation, Ausbeute und Format sind
in `tests/test_ports.py` geprüft.
