# Changelog

## Unreleased

- Z6, Teil 1: `bindoptix.signal` (Farbstoff, Ensembles je Rezeptorzustand, Anregungsrate aus der
  Feldkarte mit Orientierungsfaktor, Ausbeute und Sammlung aus der Antwortkarte je Stichprobe,
  Bleichen, Kamera mit Poisson- und Lesefehler, Schwelle für eine Fehlalarmdichte),
  Platzhalterkarten ohne eigene Optik, Konstanten h, c, ε₀; Testfälle T17–T20; Photonenbudget
  `docs/ergebnisse-z6.md`; Entscheidung D18, Messung M7.

- Z5: `bindoptix.binding` (kompetitive Langmuir-Kinetik exakt über die Matrixexponentielle,
  Gleichgewicht in geschlossener Form, Einfang ganzer Zellen mit Avidität, Port B `Occupancy`),
  zeitabhängige Fangraten in `ChannelModel` (`spot_rates_fn`), Kopplung Fouling → Einfang
  `chipsim.run_channel_coupled` mit Affinität je Spot; Testfälle T15, T16, T36.
- Korrektur Z4: Lange Zeitschritte konnten mehr binden als freie Plätze vorhanden waren (freie
  Plätze aus dem Vorschritt); `step` teilt solche Schritte jetzt konservativ, bis die Kapazität
  eingehalten ist (T14 erweitert).

- Z4, Teil 2: `bindoptix.continuum3d.ChannelVolume` (Konvektion–Diffusion im ganzen Kanalvolumen mit
  u(y, z), Spots als Polygone, BiCGSTAB mit LU-Vorkonditionierer), `chipsim.channel_volume` für
  einen Kanal des Chiplayouts; Testfall T35. Z4 abgeschlossen.

- Z4, Teil 1: `bindoptix.continuum` (Konvektion–Diffusion im Längsschnitt des Kanals, mehrere
  Spezies, kompetitive Langmuir-Kinetik auf Spotbereichen, implizites Euler-Verfahren, exakte
  Massenbilanz, stationäre Lösung für lineare Oberflächensenken), `DuctFlow.width_averaged_velocity_um_s`;
  Testfälle T12–T14.

- Z3, Teil 3: `bindoptix.wall` (behinderte Diffusion und Sedimentation nahe dem Boden, Itô-Drift),
  Option `wall_hydrodynamics` mit Kontaktabstand in `ChannelModel`; Testfall T34; Messplan für die
  Wand- und Zellparameter (`docs/messplan-wandparameter.md`). Z3 abgeschlossen.

- Z3, Teil 2: `bindoptix.chipsim` verbindet `ChipLayout` und Partikellauf (Kanalsysteme, beschnittene
  Spots, Schwerkraft je Kanal, `run_chip` mit Fängen je Gridspot); Testfall T33; Ergebnis zur
  Reihenfolge der Spots in `docs/ergebnisse-z3.md`.

- Z3, Teil 1: `bindoptix.cells` (Stäbchen als Ellipsoid, Sedimentation, Diffusion, Run-and-Tumble),
  `bindoptix.flow` (Poiseuille-Strömung im Rechteckkanal, Wandscherung, Säulenfluss),
  `bindoptix.rng` (zählerbasierte Zufallszahlen), `bindoptix.particles` (Langevin mit Brownscher
  Brücke am Boden, Gleiten, Spots, unspezifische Haftung); Testfälle T7–T11, T27, T28, T31, T32;
  erste Ergebnisse `docs/ergebnisse-z3.md`; Entscheidung D17.

- Z2, Teil 1: `bindoptix.ports` mit `EmitterResponseTable` (Tensoren A_tot, A_rad, A_det,d,
  Prüfung auf Symmetrie, Semidefinitheit und A_det ≤ A_rad ≤ A_tot, Quantenausbeute,
  Sammelanteil) und `FieldMap`/`FieldMapStack` (Port F, geschachtelte Gitter, Leistungsskalierung,
  Schutz nichtlinearer Karten), trilineare Interpolation; `bindoptix.results` (Hülle
  `raytatouille-result`); Testfall T6; Erweiterungsvorschlag an hp-FEM und CliffordBem.

- Chipgeometrie für ein schräg stehendes Grid (`grid_angle_rad`) und frei liegende gerade Kanäle
  (`Channel`), Überdeckung Spot/Kanal (`coverage`, `uncovered_spots`), Querpositionen der Spots;
  Testfall T30.

- Z1, Teil 2: `bindoptix.chip` mit 5 × 5 Spots im Gesichtsfeld Ø 5 mm, Kanäle ≤ 200 µm, Körper
  über Nurbforge, STEP-Export; Testfall T29; Abschätzung zum Einfang für dieses Layout; D15, D16.

- Z1, Teil 1: Python-Paket `bindoptix` mit CI (ruff, mypy strict, pytest), `bindoptix.units`
  (Teilchenzahlen, KBE, k_on), `bindoptix.frames` (Pose Sensorsystem → RTT-System mit
  Spiegelregel, polare und axiale Vektoren, Tensoren; Einbaulage mit Schwerkraftrichtung);
  Testfälle T1, T2.

- Entscheidungen D1–D5, D12, D13 (Legionellen, Gold zuerst, Vorkonzentrierung, gestufte Klassen);
  Abschätzung zu Probenvolumen und Einfangzeit; Testfälle T27, T28; Z0 abgenommen.

- Spezifikation Stufe Z0 (Entwurf): Roadmap, Konventionen, Testfälle, Abhängigkeiten, Quellen,
  ADR-Entwurf Datenports (F, B, Antwortkarte als M), Inventareintrag für den Fusionsplan,
  Regeln für Programmieragenten.
