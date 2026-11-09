# Testfälle

Status: ✅ umgesetzt und grün, ◐ teilweise, ○ geplant. Sollwerte nie aus dem geprüften Code.

| Nr. | Prüfung | Referenz | Stufe | Status |
| --- | --- | --- | --- | --- |
| T1 | Pose Sensorsystem ↔ RTT-System, Rundreise | Identität bis 1e−12 | Z1 | ○ |
| T2 | Spiegelregel z → −z (Beleuchtung von der Lösungsseite): polare Vektoren (Dipole) und Felder | Befund K4 (RTT-Fusion), Vorzeichen analytisch | Z1 | ○ |
| T3 | Anregungsfeld an einem Goldfilm (Kretschmann, p-Pol.) aus F | Transfer-Matrix, Resonanzwinkel wie K4 | Z2 | ○ |
| T4 | Zerfallsratentensor eines Dipols über Glas / Au / Wasser über dem Abstand | Chance–Prock–Silbey bzw. Lukosz–Kunz (unabhängige Python-Rechnung im Test) | Z2 | ○ |
| T5 | Detektionstensor: Abstrahlung in das Substrat, überkritischer Anteil, SPCE-Winkel = SPR-Winkel | Lukosz–Kunz; Lakowicz 2004 | Z2 | ○ |
| T6 | Eigenschaften der Tensoren: symmetrisch, positiv semidefinit; Mittel über Orientierungen = Spur/3 | analytisch | Z2 | ○ |
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
