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
| T7 | Freie Diffusion: mittlere quadratische Verschiebung 6Dt | Einstein, Standardfehler 4σ | Z3 | ✅ `tests/test_transport.py` |
| T8 | Sedimentation: Driftgeschwindigkeit und Gleichgewicht c(z) ∝ exp(−z·v_s/D) | Stokes; Mason–Weaver | Z3 | ✅ `tests/test_transport.py` |
| T9 | Run-and-Tumble: effektive Diffusion v²τ/(3(1−α)) | Lovely–Dahlquist 1975 | Z3 | ✅ `tests/test_transport.py` |
| T10 | Einfang an absorbierender Wand aus homogener Anfangsverteilung: Γ(t) = 2c₀√(Dt/π); Wandkontakt auch zwischen den Schritten (Brownsche Brücke) | Smoluchowski | Z3 | ✅ `tests/test_transport.py` |
| T11 | Determinismus: gleicher Seed, gleiche Ergebnisse, unabhängig von der Aufteilung in Pakete (Zufallszahlen bitgleich, Gleitkomma bis auf Rundung wegen SIMD-Pfaden von NumPy, Zustände gleich) | Rundreise | Z3 | ✅ `tests/test_transport.py` |
| T12 | Transportbegrenzter Einfang im Kanal: lineare Scherung gegen Lévêque (1 %); Parabelprofil etwa 2 % darunter | Lévêque-Lösung, Gesamtfluss 0,8075 c₀ (D²γ̇L²)^(1/3) | Z4 | ✅ `tests/test_continuum.py` |
| T13 | Übergang reaktions- zu transportbegrenzt über die Damköhler-Zahl: Grenzfälle (2 %, 4 %), Da = 1 gegen Reihenschaltung J_R/(1 + Da) (15 %); reaktionsbegrenzte Langmuir-Kinetik | squires2008; analytische Langmuir-Lösung | Z4 | ✅ `tests/test_continuum.py` |
| T14 | Massenerhaltung (Volumen + Oberfläche + Auslass), zwei konkurrierende Spezies, Belegung ≤ Γ_max auch bei sehr langen Schritten (Teilschritte) | Bilanz bis 1e−10 relativ | Z4 | ✅ `tests/test_continuum.py` |
| T15 | Langmuir und kompetitive Langmuir-Kinetik mit Gewichten α: Zeitverlauf einer Spezies analytisch, Gleichgewicht in geschlossener Form, Zeitverlauf mehrerer Spezies gegen ODE-Löser, irreversible Bindung | analytisch, `scipy.integrate.solve_ivp` | Z5 | ✅ `tests/test_binding.py` |
| T16 | Rezeptorzustände (Port B): Anteile je Spot summieren sich zu 1, liegen in [0, 1], erster Zustand „free“; Dateiformat-Rundreise mit Zellereignissen | Grenzfall, Rundreise | Z5 | ✅ `tests/test_binding.py` |
| T17 | Einzelner Emitter im homogenen Medium: Rate = σ I/(hν) · q₀ · η; Orientierungsfaktor 3 für Dipol parallel zur Polarisation | analytisch, σ = 3,82·10⁻²¹ cm² · ε | Z6 | ✅ `tests/test_signal.py` |
| T18 | Modifizierte Quantenausbeute aus den Tensoren | Formel mit γ_nr⁰ = γ_r⁰(1 − q₀)/q₀, unabhängige Rechnung | Z6 | ✅ `tests/test_ports.py`, `tests/test_signal.py` |
| T19 | Bleichen: monoexponentiell bei konstanter Anregung, Photonen im Zeitfenster | analytisch, Quadratur | Z6 | ✅ `tests/test_signal.py` |
| T20 | Kamera: Mittel und Varianz; Überschreitungswahrscheinlichkeit Poisson + Gauß gegen Monte Carlo; Schwelle für eine Fehlalarmdichte | Poisson plus Gauß-Leserauschen, Monte Carlo | Z6 | ✅ `tests/test_signal.py` |
| T21 | Spot-Erkennung auf synthetischen Bildern mit bekannter Wahrheit: Differenzbild + Persistenz findet ≥ 24/25 Zellen ohne Fehlalarm; ohne Referenzbild > 10 Fehlalarme durch die Ungleichmäßigkeit, ohne Persistenz zählen flüchtige Partikel; Zählung je Spot mit Randtoleranz; Photonenerhaltung beim Abbilden | bekannte Wahrheit | Z7 | ✅ `tests/test_readout.py` |
| T22 | Protokoll: Waschschritt entfernt reversibel Gebundenes nach exp(−k_off t) (2 %) | analytisch | Z7 | ✅ `tests/test_readout.py` |
| T23 | Erwärmung: Bilanz ∫q dV = absorbierte Leistung aus F | Energieerhaltung | Z8 | ○ |
| T24 | Optische Kraft auf eine Kugel aus der Kraftkarte | Mie (CliffordBem-Sollwerte) | Z8 | ○ |
| T25 | Klassifikation: Konfusion bei einer gefangenen Zelle in geschlossener Form; exakte Aufzählung gegen Monte Carlo; Unabhängigkeit von der Konzentration | analytisch, Aufzählung | Z10 | ✅ `tests/test_classify.py` |
| T26 | Optimierer findet das bekannte Optimum eines Spielproblems | analytisch | Z10 | ○ |
| T27 | Sinkgeschwindigkeit und Diffusion eines Stäbchens (Rotationsellipsoid, Orientierung gemittelt) | geschlossene Form nach Perrin (unabhängig von den Achsenreibungen nach Happel–Brenner); Kugelgrenzfall Stokes und Stokes–Einstein | Z3 | ✅ `tests/test_transport.py` |
| T28 | Einfang im laminaren Durchfluss: Anteil ∫ min(q(y), L v_s) dy / Q mit dem Säulenfluss q(y) des Rechteckkanals (Hazen je Säule) | Hazen, Reihenlösung termweise integriert | Z3 | ✅ `tests/test_transport.py` |
| T29 | Chipgeometrie: Volumen der Kanäle (Nurbforge) gegen l·b·h, Spotflächen, Kanalboden im Gesichtsfeld gegen numerische Quadratur der Sehnenlänge, Spots in ihrem Kanal, STEP-Rundreise, ungültige Layouts | analytisch, `scipy.integrate.quad` | Z1 | ✅ `tests/test_chip.py` |
| T30 | Schräges Grid und freie Kanäle: Querpositionen bei tan θ = 1/5 alle verschieden mit Abstand p/√26, Passung drehinvariant, Überdeckung Spot/Kanal gegen Monte Carlo, gedrehte Kanalkörper, sich kreuzende Kanäle abgelehnt | analytisch, Monte Carlo mit Seed | Z1 | ✅ `tests/test_chip.py` |
| T31 | Gleiten am Boden mit u_wall = s · a · ∂u/∂z und Einfang beim Erreichen eines Spots | Wandscherung als Differenzenquotient einer langen Reihe; Ankunftszeit analytisch | Z3 | ✅ `tests/test_transport.py` |
| T32 | Strömung im Rechteckkanal: Durchfluss aus der Reihe, Grenzfall breiter Kanal (ebene Poiseuille-Strömung), Wandscherung, Haftbedingung | bruus2008, Quadratur | Z3 | ✅ `tests/test_transport.py` |
| T33 | Chip und Partikel: Sensor- ↔ Kanalsystem (Rundreise), Spots der Gridreihen im Kanalsystem, beschnittene Spotflächen = Überdeckung, Schwerkraft gedreht, Lauf über das Layout = direkter Lauf, Fänge je Gridspot | analytisch; Gleichheit mit dem direkten Lauf | Z3 | ✅ `tests/test_chipsim.py` |
| T34 | Wandhydrodynamik: λ⊥ gegen Brenners exakte Reihe (0,05 ≤ h/a ≤ 10, 2 %), Grenzwerte und Ableitung, λ∥ nach Faxén, Boltzmann-Gleichgewicht mit ortsabhängigem D (Itô-Drift), gebremste Annäherung an den Boden | brenner1961 (Reihe im Test), bevan2000, Quadratur | Z3 | ✅ `tests/test_wall.py` |
| T35 | Kontinuum im Kanalvolumen (3D): quer gleichförmiger Fall = 2D-Lösung × Breite (1e−10), absorbierender Spot im Rechteckkanal gegen lokale Lévêque-Lösung mit der Wandscherung γ̇(y) (4 %), Massenbilanz und Symmetrie mit zwei Spezies, Spots des Chiplayouts als Reaktionsflächen | T12, Lévêque, Bilanz | Z4 | ✅ `tests/test_continuum3d.py` |
| T36 | Kopplung Fouling → Einfang: ohne Fouling identisch zum Lauf mit konstanter Fangrate; mit Fouling weniger Fänge, keine Fänge nach Füllung der Spots, Füllfront wandert stromab, Port B konsistent; Affinität je Spot (Zeile der Bindungsmatrix) | Gleichheit mit dem ungekoppelten Lauf, Grenzfälle | Z5 | ✅ `tests/test_binding.py` |
| T37 | Protokollmechanik: Fortsetzen eines Laufs bitgleich zum Durchlauf, Einspülplan, wartende Zellen bewegen sich nicht, Inkubation ohne Strömung ohne Drift | Gleichheit, analytisch | Z7 | ✅ `tests/test_readout.py` |
| T38 | Trajektorienbibliothek: vorhergesagte Fänge je Spot erklären direkte Läufe mit endlichen Fangraten und unspezifischer Haftung (gleiche Zufallsströme, Poisson-Binomial-Streuung) | direkter Partikellauf | Z10 | ✅ `tests/test_classify.py` |

T3–T5 prüfen exportierte Karten aus hp-FEM bzw. CliffordBem und warten auf deren Export
(`docs/vorschlaege/export-antwortkarte.md`). Datenmodell, Interpolation, Ausbeute und Format sind
in `tests/test_ports.py` geprüft.
