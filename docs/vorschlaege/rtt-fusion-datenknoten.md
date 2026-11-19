# Vorschlag: BindOptiX als Datenknoten und Zeitschleife in RTT-Fusion

> Erweiterungsvorschlag nach AGENTS.md Regel 4 (BindOptiX liest RTT-Fusion nur). Bezug: Stufe Z12,
> ADR-Entwurf Datenports (`docs/adr-entwurf-datenports.md`), in RTT-Fusion `docs/fusionsplan.md`
> Abschnitt 6.2 (D13), 8 (K5, K8) und 9, `docs/konzept-node-workflow.md` Abschnitt 7 und L6,
> Stufe R20. Stand gelesen: RTT-Fusion `main` 2f4db25 (2026-11-19). Eingereicht als
> RTT-Fusion #106; Antworten am Ende.

## Ausgangslage

- **D13 ist im Fusionsplan noch Vorschlag** (Abschnitt 6.2 „Vorgeschlagen (D13)“, Tabelle 9). Die
  Namen der Daten-Sockets stehen aber schon fest (L6) und im Code (`graph/model.py`,
  `DATA_TYPES`): `field_map` (F), `occupancy` (B), `response_map` (M).
- **Erzeuger gibt es inzwischen.** hp-FEM exportiert F und M (hp-FEM #214, PRs #220–#224);
  BindOptiX hat sie abgenommen (Testfälle T3–T5 mit Karten aus hp-FEM, #23, #34, #37). BindOptiX
  liest Q aus dem `DipoleAndDistanceTracker` von BioOptoMD und erzeugt B (T43,
  `bindoptix.signal`, `bindoptix.binding.Occupancy`). CliffordBem exportiert noch nicht
  (CliffordBem #1).
- **BindOptiX ist kein Lichtlöser.** Der Adaptervertrag `describe`/`prepare`/`apply` gilt für
  Licht an Ports; BindOptiX ist ein Datenknoten und, für den Ablauf über der Zeit, ein Meta-Block
  (Konzept Abschnitt 7: „Zeitschleife … BindOptiX: Optik einmal je Design, Biologie über s bis h“).

## Vorschlag

### 1. D13 übernehmen

F, B und die Antwortkarte als Datenports wie N und Q, ohne neue Lichtdarstellung. Dateiformat ist
der Umschlag `raytatouille-result` (Schema 0.1.x) mit den Typen `FieldMap`, `Occupancy` und
`EmitterResponseTable`, wie ihn hp-FEM bereits schreibt; Inhalt, Einheiten, Normierung und
Vorzeichen nach dem ADR-Entwurf und `docs/konventionen.md` von BindOptiX. Die Zeilen „Vorgeschlagen
(D13)“ in Abschnitt 6.2 würden Zeilen der Datenport-Tabelle mit den Socket-Namen aus L6.

### 2. Datenknoten `bindoptix.signal`

Nach dem Muster von `biooptomd.layer` (Kategorie DATA bzw. ANALYSIS, Adapter
`rtt_fusion/adapters/bindoptix.py`):

| Socket | Richtung | Datentyp | Inhalt |
| --- | --- | --- | --- |
| `excitation` | ein | `field_map` | F bei der Anregungswellenlänge |
| `response` | ein | `response_map` | Antwortkarte bei der Emissionswellenlänge |
| `q_<zustand>` | ein, je Zustand | `Q` | Ensemble je Rezeptorzustand (dynamische Sockets) |
| `occupancy` | ein | `occupancy` | Belegung B je Spot über der Zeit |
| `rates` | aus | `result` | erkannte Rate je Rezeptor, Spot, Zeit und Detektorbereich |

Parameter: Farbstoff (ε, q₀, Anregung, Emission, Bleichausbeute), Leistung, Bestrahlungsstärke je
mW mit Status L/K/P, Detektorbereiche der Optik. Die Rechnung ist `bindoptix.signal.spot_rates_s`;
der Adapter importiert BindOptiX verzögert wie `adapters/hpfem.py` hp-FEM (BindOptiX braucht ohne
C++-Kern nur NumPy und SciPy, `pip install … -C cmake.define.BINDOPTIX_BUILD_CORE=OFF`).

**Q ohne doppelte Umrechnung.** `biooptomd.layer` gibt Q im RTT-System aus (mm, gespiegelt je nach
`light_from`). BindOptiX rechnet im Sensorsystem (µm, Lösung bei z > 0) und liest die
Tracker-Ausgabe direkt (`EmitterEnsemble.from_tracker`). Vorschlag: Der Datentyp `Q` trägt die
Tracker-Reihe (n_replicas, n_frames, 7) in nm mit der Konvention von BioOptoMD; die Umrechnung
ins RTT-System macht der Knoten, der Licht damit rechnet. Die Pose Sensor ↔ RTT mit Spiegelregel
(DP8) bleibt an einer Stelle.

### 3. Meta-Block `bindoptix.time_loop` (mit R20)

Umschließt den optischen Untergraphen (F, M; einmal je Design gerechnet und gecacht) und führt
Protokoll, Transport, Einfang und Bindung von BindOptiX über der Zeit aus; Ausgaben `occupancy`,
Bilder über der Zeit und Kennzahlen (`result`). Neu gerechnet wird die Optik nur bei
label-freien Sensoren (ADR-Entwurf, Grenzen). Er setzt die Meta-Blöcke aus R20 voraus und ist
deshalb der zweite Schritt.

### 4. Konformität

Für Datenknoten gibt es keinen Konformitätssatz wie für Lichtlöser (Plan 6.5). Für BindOptiX
gelten die Prüfungen DP1–DP8 des ADR-Entwurfs (Konzept L3); welche Testfälle von BindOptiX sie
abdecken, steht dort. In RTT-Fusion genügte `tests/test_adapter_bindoptix.py` mit:

- Lesen der Formate (Umschlag, Version, Typ) und Ablehnung falscher Typen;
- DP7 (Mischung der Q-Ensembles nach B, ohne Analyt nur Q_frei) über den Knoten;
- DP8 (Spiegelung von Q zwischen Sensor- und RTT-System) gegen `biooptomd.dipoles_from_tracker`;
- einer Kette aus `biooptomd.layer` → `bindoptix.signal` mit einer Antwortkarte aus hp-FEM.

## Reihenfolge

1. D13 im Fusionsplan entscheiden (Abschnitt 6.2, Tabelle 9).
2. Datenknoten `bindoptix.signal` mit Adapter und Test.
3. K5 (BioOptoMD → Löser → Fernfeld → Raytatouille), weil K8 laut Plan danach kommt.
4. Meta-Block Zeitschleife mit R20, dann Pilot K8.

## Offene Fragen an RTT-Fusion

- Soll `Q` die Tracker-Reihe in der Konvention von BioOptoMD tragen (Vorschlag oben) oder die
  Dipole im RTT-System, wie `biooptomd.layer` sie heute ausgibt?
- Liegt der Adapter in `rtt_fusion.adapters` (wie alle bisherigen) oder in einem eigenen Paket
  (Plan: „Adapter in getrennten Paketen“)?
- Wie werden die Detektorbereiche einer Antwortkarte (z. B. `substrate_na_lt_1.45`) an die
  Detektoren des optischen Graphen gebunden: über den Namen oder über einen eigenen Parameter?

## Antworten aus RTT-Fusion #106 (Maintainer, 2026-11-19)

- **D13 wird übernommen**, die Socket-Namen bleiben.
- **Q** trägt die Tracker-Reihe (n_replicas, n_frames, 7) = (R_F [nm], μ̂, d [nm]) im System von
  BioOptoMD (Oberfläche z = 0, Lösung z > 0) mit Metadaten des Laufs (Quelle, Zeitschritt). Ins
  RTT-System rechnet nur der Lichtknoten um (`dipoles_from_tracker`); DP8 bleibt an einer Stelle.
  `EmitterEnsemble.from_tracker` liest diese Reihe bereits.
- **Adapter** in `rtt_fusion/adapters/bindoptix.py` mit verzögertem Import (BindOptiX ist MIT).
- **Detektorbereiche** über einen eigenen Parameter von `bindoptix.signal`: Bereichsname der
  Antwortkarte → Detektor bzw. Ausgang des Graphen; fehlende oder doppelte Namen ergeben eine
  Diagnose; ohne Parameter nur bei genau einem Bereich.
- **Reihenfolge in RTT-Fusion:** (1) D13 und Vertrag für Q im Fusionsplan, (2) Evaluator für
  `biooptomd.layer`, (3) `bindoptix.signal` mit Adapter und `tests/test_adapter_bindoptix.py`,
  (4) K5, (5) `bindoptix.time_loop` (R20) und K8; Beginn nach der Stufe R32.
- Hinweis: `biooptomd.layer` hat noch keinen Evaluator; Q wird im Graphen bisher weder erzeugt
  noch verbraucht.
