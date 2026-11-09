# Vorschlag: Export von Antwortkarte und Feldkarte aus hp-FEM und CliffordBem

> Erweiterungsvorschlag nach AGENTS.md Regel 4 (BindOptiX liest dort nur). Bezug: ADR-Entwurf
> Datenports (F, Antwortkarte als M), Fusionsplan D12 und D13, Stufe Z2, Testfälle T3–T5 und
> DP1–DP8. Dateiformat und Prüfungen: `bindoptix.ports` (`EmitterResponseTable`, `FieldMap`).

## Was gebraucht wird

| Größe | Inhalt | Erzeuger |
| --- | --- | --- |
| A_tot(r) | abgegebene Leistung eines Dipols am Ort r, Tensor 3×3, normiert auf den freien Dipol im Medium bei r | hp-FEM (M17, M18), CliffordBem (Zerfallsraten) |
| A_rad(r) | ins Fernfeld abgestrahlte Leistung, über alle Richtungen beider Halbräume | wie oben |
| A_det,d(r) | Leistung, die Detektorbereich d erreicht | Fernfeld des Lösers, Gewichte J_d(k̂) aus der K5-Kette (`rtt-field`, Raytatouille) |
| E(r) | Anregungsfeld je mW auf einem Gitter um die Spots | hp-FEM, CliffordBem (`nearfield`) |

## Drei Rechnungen je Ort genügen

Die Löser sind linear. Mit den Feldern E⁽ʲ⁾, H⁽ʲ⁾ der Einheitsdipole entlang x, y, z am Ort r
ergibt sich jede Orientierung als Überlagerung; alle Leistungen sind quadratische Formen:

- A_tot,ij aus dem Poynting-Fluss durch eine geschlossene Fläche um die Quelle,
  ½ Re ∮ ½ (E⁽ⁱ⁾ × H⁽ʲ⁾* + E⁽ʲ⁾ × H⁽ⁱ⁾*) · n dS (besser konditioniert als das Eigenfeld bei r, wie in
  hp-FEM M17 beschrieben), geteilt durch die Leistung des freien Dipols;
- A_rad,ij = Re ∫ E∞⁽ⁱ⁾* · E∞⁽ʲ⁾ dΩ über das Fernfeld, normiert wie oben;
- A_det,d,ij = Re ∫ (J_d E∞⁽ⁱ⁾)* · (J_d E∞⁽ʲ⁾) dΩ, mit J_d(k̂) der Jones-Matrix vom Fernfeld in die
  Detektorebene des Bereichs d (null außerhalb der Richtungen, die d erreichen).

Die Symmetrie A_ij = A_ji folgt aus der Reziprozität; die Asymmetrie der berechneten Matrix ist
ein Konvergenzmaß und gehört in die Herkunftsdaten. BindOptiX lehnt Tabellen ab, die nicht
symmetrisch und positiv semidefinit sind oder A_det ≤ A_rad ≤ A_tot verletzen.

**Gültigkeit von A_det:** Richtungen werden inkohärent summiert. Das gilt, wenn der Detektorbereich
groß gegen die Punktbildfunktion ist (Spot 200 µm, Bild auf der Kamera viele Pixel). Für einzelne
Pixel ist das kohärente Feld in der Bildebene nötig; dann liefert die K5-Kette dieses Feld je
Einheitsdipol, und A_det wird über die Pixelfläche integriert.

## Gitter

Positionen auf einem achsparallelen Gitter im Sensorsystem (µm, Oberfläche z = 0, Lösung z > 0).
Für ebene Schichtstapel genügt eine Achse z (x und y invariant, Achsenlänge 1). Über Kavitäten
feine Gitter, darüber grobe; BindOptiX interpoliert trilinear (konvexe Kombination, erhält
positive Semidefinitheit). Für das Feld auf 5 mm: Gitter je Spot statt eines Gitters über das ganze
Feld; der Verlauf der Sammeleffizienz über das Feld kommt aus J_d(k̂) je Spotlage.

## Format

Hülle `raytatouille-result` (Raytatouille ADR 0023), Typen `EmitterResponseTable` und `FieldMap`
mit den Schlüsseln aus `bindoptix.ports`. Herkunft (`source`): Löser, Version, Commit,
Inhalts-Hash der Struktur, Konvergenz (Netzfeinheit, Asymmetrie, Energiebilanz).

## Abnahme im Löser

- Dipol über Glas/Au/Wasser gegen Chance–Prock–Silbey bzw. Lukosz–Kunz (T4, DP4).
- Dipol vor idealem Spiegel gegen die geschlossene Spiegeldipol-Lösung:
  senkrecht 1 + 3 (sin x/x³ − cos x/x²), parallel 1 − (3/2)(sin x/x + cos x/x² − sin x/x³),
  x = 2kd (Grenzwerte 2 und 0 für d → 0).
- Ohne Schicht: A_tot = A_rad = I.
