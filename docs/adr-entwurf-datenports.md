# ADR 00xx (Entwurf): Datenports F (Feldkarte) und B (Belegung), Antwortkarte als Tabelle M

> **Entwurf aus BindOptiX (Stufe Z0) für den Fusionsplan** (Fieldtracing, Branch
> `docs-fusion-plan`, `docs/rtt-fusion/plan.md`, Abschnitte 6.2 und 6.3). Abschnitt 6.2 verlangt
> für neue Darstellungen ein ADR. Bezug: Pilot K4 (RTT-Fusion, abgenommen), Pilot K5 (geplant, P5).

**Status:** freigegeben in BindOptiX am 2026-11-09 (PR #1); im Fusionsplan Entscheidung D13, dort offen

## Kontext

Die Fusion koppelt Licht über die Darstellungen R, H, G, S, M, P und liefert Strukturen und
Quellen über die Datenports N und Q **einseitig** an Lichtgebiete (BioOptoMD → Löser). Ein
Biosensor braucht zusätzlich die Gegenrichtung und die Zeit:

1. **Licht → Biologie.** Anregungsrate und Bleichen der Fluorophore, Erwärmung (hp-FEM M9) und
   optische Kräfte hängen vom Feld im Volumen ab. G ist dafür ungeeignet: G ist ein Querfeld in
   einer Ebene mit Richtung und Leistung, kein Volumenfeld.
2. **Zeit.** Bindung und Transport laufen über Sekunden bis Stunden, die Optik ist zeitharmonisch.
   Je Zeitschritt neu zu rechnen ist weder nötig noch bezahlbar. Was sich mit der Zeit ändert, ist
   nur, **welche** Rezeptoren in welchem Zustand sind.
3. **Q je Zustand.** BioOptoMD liefert Dipolensembles für einen Rezeptorzustand (frei, Ziel
   gebunden, …). Das Signal entsteht aus der zeitabhängigen Mischung dieser Ensembles über der
   Oberfläche.
4. **Sammlung.** Emitter strahlen untereinander inkohärent (Plan 7.3). Die Kette Emitter →
   Nahfeld → SIL → Relay → Detektor ist für einen gegebenen Aufbau linear in der Leistung jedes
   Emitters und hängt nur von Ort und Orientierung ab. Sie lässt sich einmal tabellieren.

## Entscheidung

### 1. Neuer Datenport F (Feldkarte)

- **Inhalt:** komplexes elektrisches Feld **E** ∈ ℂ³ auf einem regulären kartesischen Gitter im
  Sensorsystem, bei einer Vakuumwellenlänge; optional abgeleitete Größen: Absorptionsdichte q
  (W/m³ je mW) und Temperaturerhöhung ΔT (K je mW).
- **Normierung:** je 1 mW einfallender Leistung am Eingangsport der Fusion; Zeitfaktor e^(−iωt).
  E skaliert mit √P, q mit P, ΔT mit P, solange keine Temperaturrückkopplung gerechnet wurde.
  Mit Rückkopplung (hp-FEM, ε(T)) trägt die Karte die Leistung, für die sie gilt, und das Flag
  `nonlinear`; dann ist sie nur bei dieser Leistung gültig.
- **Gitter:** mehrere geschachtelte Gitter je Karte (fein an Kavitäten und Hotspots, grob im
  Kanal). Ausgewertet wird das feinste Gitter, das den Punkt enthält. Begründung: Nanometer an
  den Kavitäten gegen Mikrometer im Kanal; ein gemeinsames Netz für Optik und Fluidik ist
  ausdrücklich nicht vorgesehen.
- **Metadaten:** Pose zum RTT-System, Medium je Gitterpunkt (Materialkennung), Herkunft (Löser,
  Version, Inhalts-Hash, Konvergenz), wie die Beugungstabelle.
- **Semantik:** Datenport ohne Richtung und ohne Leistungsfluss; im Ausführungsgraphen eine
  Abhängigkeit, kein Pfadabschnitt (wie N und Q).
- **Erzeuger:** hp-FEM (Nédélec-Feld, `absorbed_power_density`, `physics::Thermal`), CliffordBem
  (`nearfield`), SpeckleBem (Nahfeld, im Adapter konjugiert).

### 2. Neuer Datenport B (Belegung)

- **Inhalt:** je Spot und Zeitpunkt die Anteile der Rezeptorzustände f_s(x, y, t) auf einem
  Flächengitter (Kontinuum) **und** eine Ereignisliste für Zellen (Ort, Zeit, Spezies, Zustand),
  weil gebundene Bakterien diskrete Objekte sind.
- **Kopplung an Q:** BioOptoMD liefert je Zustand s ein Ensemble Q_s (Stichproben von Ort relativ
  zum Rezeptor, Orientierung, Abstand zur Oberfläche). Das wirksame Ensemble an einem Ort ist die
  Mischung Σ_s f_s(x, y, t) · Q_s. Gewichtet werden die Stichproben, nicht die Mittelwerte.
- **Kopplung an N:** nur für molekulare Schichten; lateral gemittelt wie in K4. Für Zellen gilt N
  nicht (Mikrometerobjekte, keine Schicht). Zellen als Streukörper gehen später über Geometrie an
  CliffordBem oder SpeckleBem, nicht über diesen ADR.
- **Erzeuger:** BindOptiX. Umgesetzt (Z5, `bindoptix.binding.Occupancy`) zunächst als Mittel je Spot
  über der Zeit plus Ereignisliste; ein Flächengitter folgt, wenn Z6 die Verteilung innerhalb der Spots
  braucht.

### 3. Antwortkarte als Tabelle der Darstellung M (keine neue Darstellung)

Für einen Dipol der Orientierung μ̂ am Ort r ist jede Leistung, die er erzeugt, eine
quadratische Form in μ̂. Für reelles μ̂ genügt je Größe eine reelle, symmetrische 3×3-Matrix
(sechs Zahlen je Ort und Wellenlänge):

| Tensor | Bedeutung | aus |
| --- | --- | --- |
| A_tot(r, λ_em) | gesamte abgegebene Leistung (strahlend + absorbiert), normiert auf den freien Dipol im Medium am Ort r | Im G(r, r) bzw. Poynting-Fluss um die Quelle |
| A_rad(r, λ_em) | ins Fernfeld abgestrahlte Leistung, normiert wie oben | Fernfeld über die ganze Kugel |
| A_det,d(r, λ_em) | Leistung, die Detektorbereich d erreicht (Pixel oder Spot), normiert wie oben | Fernfeld → S → `rtt-field` → Raytatouille (K5) |

Daraus folgen modifizierte Quantenausbeute und Sammeleffizienz mit der intrinsischen
Quantenausbeute q₀ des Farbstoffs:

q = a_rad / (a_tot + (1 − q₀)/q₀),  η_d = a_det,d / a_rad,  mit a = μ̂ᵀ A μ̂.

Die Anregungsrate folgt aus F: R_exc ∝ |μ̂_abs · E(r)|². Das Signal je Emitter ist damit
R_exc · q · η_d.

**Wichtig:** Das Produkt ist nicht quadratisch in μ̂. Gemittelt wird deshalb über die
Stichproben des Q-Ensembles, nicht mit einem vorher gemittelten Tensor ⟨μ̂μ̂ᵀ⟩. Die Tensoren
selbst sind exakt; nur die Mittelung muss spät erfolgen.

Format: Tabelle in `raytatouille-result` (ADR 0023), Typ `EmitterResponseTable`, Achsen
(Ort, λ_em, Detektorbereich), Interpolation nur über glatte Bereiche (analog Fieldtracing ADR 0003);
an Kavitäten adaptive Stützstellen. Herkunft wie bei der Beugungstabelle.

## Konformitätstests

| Nr. | Prüfung | Referenz |
| --- | --- | --- |
| DP1 | F an einer Platte und am Goldfilm (Kretschmann) | Airy, Transfer-Matrix, K4-Resonanz |
| DP2 | ∫ q dV = absorbierte Leistung aus der Leistungsbilanz des Lösers | Energieerhaltung |
| DP3 | F aus hp-FEM gegen F aus CliffordBem für eine Goldkugel | Querprobe, Mie |
| DP4 | A_tot, A_rad über Glas / Au / Wasser über dem Abstand | Chance–Prock–Silbey, Lukosz–Kunz |
| DP5 | A_det für ein ideales Objektiv: Anteil in den Raumwinkel der NA, überkritischer Anteil, SPCE-Winkel | Lukosz–Kunz, Lakowicz 2004 |
| DP6 | Tensoren symmetrisch, positiv semidefinit; Mittel über die Kugel = Spur/3 | analytisch |
| DP7 | B → Q: Gewichte summieren zu 1; ohne Analyt nur Q_frei | Grenzfall |
| DP8 | Spiegelung z → −z an F, Q und den Tensoren | Vorzeichen analytisch (Befund K4) |

## Grenzen

- **Quasistatik:** Die Optik ändert sich durch die Bindung nicht. Für Fluoreszenz gilt das, weil
  Δn der Rezeptorschicht klein ist (K4: ~0,5° SPR-Verschiebung je 0,36 nm optischer Dicke).
  Für label-free Sensoren ist die Rückwirkung das Signal selbst; dort wird nicht F aktualisiert,
  sondern N aus B neu gebildet und die Lichtrechnung je Zeitpunkt wiederholt.
- **Orientierung während des Zyklus:** Die Tensoren setzen voraus, dass sich der Fluorophor
  zwischen Absorption und Emission nicht dreht, oder dass μ̂_abs und μ̂_em getrennt aus Q kommen.
  Rotationsdiffusion im ns-Bereich wird über die Ensemble-Stichproben abgebildet, nicht hier.
- **Sättigung:** linear in der Anregung; Sättigung und Triplett-Blinken sind Modellteile von
  BindOptiX, nicht des Ports.
- **Zellen als Streuer:** nicht Teil dieses ADR.

## Verworfene Alternativen

- **Optik in BindOptiX rechnen** (de-Feijter-, SERS-, Fluoreszenz-Kernel im eigenen Code):
  verdoppelt hp-FEM, CliffordBem und die K5-Kette.
- **G für Volumenfelder verwenden:** G ist ein Lichtport in einer Ebene mit Richtung und Leistung;
  die Bedeutung würde verwischt.
- **Feldkarten je Zeitschritt:** unnötig wegen der Trennung der Zeitskalen.
- **Gemittelte Orientierungstensoren:** falsch, weil das Signal nicht quadratisch in μ̂ ist.

## Folgen

- Fusionsplan 6.2: zwei neue Zeilen (F, B), Antwortkarte als Inhalt von M.
- Fusionsplan 6.3: Umwandlungen „F → Anregungsrate“, „B + Q → gewichtetes Q-Ensemble“,
  „Fernfeld eines Dipols → S → Detektor → A_det“.
- Fusionsplan 8: neuer Pilot K8 (siehe `docs/inventar-fusionsplan.md`).
- Adapter: hp-FEM und CliffordBem brauchen einen Export von F und der Tensoren; beide haben die
  Rechnungen schon (Nahfeld, Dipolquelle, Zerfallsraten). Erweiterungen werden dort als Issues
  vorgeschlagen, wie in Entscheidung D12 des Fusionsplans.
