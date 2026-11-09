# Changelog

## Unreleased

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
