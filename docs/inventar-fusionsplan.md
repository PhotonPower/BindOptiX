# Inventareintrag für den Fusionsplan

> Vorlage: Fieldtracing, Branch `docs-fusion-plan`, `docs/rtt-fusion/plan.md`, Anhang A und
> Abschnitte 2.1, 8, 9. Einzufügen per PR dort. Stand BindOptiX: Commit vom 2026-11-03
> (README, LICENSE).

## Anhang A: BindOptiX

| Feld | Inhalt |
| --- | --- |
| Name, Repository, Version | BindOptiX, `PhotonPower/BindOptiX`, 0.0 (nur README und LICENSE; Name in der GitHub-Beschreibung noch „PhotoKinX“, D1) |
| Verfahren, Gültigkeitsbereich | Kein Lichtlöser. Transport (Langevin-Partikel für Zellen, Konvektion–Diffusion für Moleküle), Oberflächenkinetik, Assay-Abläufe, Signal- und Kameramodell aus optischen Karten. Quasistatisch: Optik einmal je Design, Biologie über s bis h |
| Sprache, Python-Anbindung | Python-Referenz (NumPy, SciPy), danach C++/CUDA-Kern mit Python-Anbindung (D8) |
| Lizenz | MIT |
| Eingaben | Geometrie (Nurbforge), Q je Rezeptorzustand (BioOptoMD), F (hp-FEM, CliffordBem), Antwortkarte als M (K5-Kette), Strömung (analytisch, später extern), Assay-Protokoll |
| Ausgaben | B (Belegung über der Zeit), Photonenraten je Spot und Pixel, synthetische Kamerabilder, Kennzahlen (Nachweisgrenze, Zeit, Konfusionsmatrix) |
| Konventionen | keine Wellenoptik intern; F mit e^(−iωt), je mW; Sensorsystem z = 0, Lösung z > 0 (wie BioOptoMD); µm, s, K intern, mm im Austausch; Teilchenzahlen statt mol |
| Linear / reziprok / deterministisch | nichtlinear (Kinetik), stochastisch mit Seed; Ergebnisse mit Standardfehler |
| Ressourcen | Referenz auf der CPU; GPU (RTX 6000 Ada, RTX 3090) für > 10⁶ Partikel |
| Validierung | geplant, `docs/testfaelle.md` (T1–T26), Konformität DP1–DP8 im ADR-Entwurf |

## Zeile für Abschnitt 2.1 (Inventar der Module)

| Modul | Verfahren | Stand | Sprache, Anbindung | GPU |
| --- | --- | --- | --- | --- |
| **BindOptiX** (geplant) | Zeit-, Transport- und Biologieschicht: Partikel und Kontinuum, Kinetik, Signal aus F und Antwortkarten | Z0 (Spezifikation) | Python, später C++/CUDA | vorgesehen (Z11) |

## Zeile für die Konventionstabelle in 2.1

| Modul | Zeitfaktor | Einheiten | Geometrie und Polarisation | Anpassung im Adapter |
| --- | --- | --- | --- | --- |
| BindOptiX | – (nur Daten); F mit e^(−iωt) | µm, s, K, mW; Teilchen/µm³ | Sensoroberfläche z = 0, Lösung z > 0; Schwerkraft als Vektor | Pose und Spiegelregel zum RTT-System |

## Zeile für Abschnitt 8 (Pilotfälle)

| Nr. | Inhalt | Module | Pfad | Prüfung | Phase |
| --- | --- | --- | --- | --- | --- |
| K8 | Biosensor über der Zeit: Einfang von Bakterien, Rezeptorzustände, Fluoreszenz durch gedruckte SIL mit Spiegelrand und Relay-Optik auf die Kamera | BioOptoMD + hp-FEM/CliffordBem + `rtt-field` + Raytatouille + Voxelspaetzle + BindOptiX | Q je Zustand, F, Antwortkarte (M) → BindOptiX → B → Signal | Teilprüfungen DP1–DP8, T3–T10; ohne Analyt = Untergrund allein | P5, nach K5 |

## Zeile für Abschnitt 9 (Entscheidungen)

| Nr. | Frage | Vorschlag |
| --- | --- | --- |
| D13 | Datenports F und B, Antwortkarte als M | ADR-Entwurf aus BindOptiX (`docs/adr-entwurf-datenports.md`) |
