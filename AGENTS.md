# Regeln für Programmieragenten

Diese Datei liest jeder Agent zu Beginn jeder Sitzung. Sie gilt zusammen mit `docs/roadmap.md`,
`docs/konventionen.md`, `docs/abhaengigkeiten.md`, `docs/quellen.md` und `docs/testfaelle.md`.

## Projekt in einem Satz

BindOptiX ist die Zeit-, Transport- und Biologieschicht des PhotonPower-Frameworks für optische
Biosensoren: Einfang von Bakterien und Molekülen, Bindung und Rezeptorzustände, Assay-Abläufe und
das Signal auf der Kamera, gekoppelt an die Optik über vorab gerechnete Karten aus der Fusion
(BioOptoMD, hp-FEM, CliffordBem, Raytatouille). BindOptiX rechnet kein Licht.

## Bauen und testen

Ab Stufe Z1 (Paketgerüst folgt):

```bash
pip install -r requirements-dev.txt && pip install -e .
pytest
ruff check . && ruff format --check . && mypy
```

## Arbeitsregeln

1. **Ein Auftrag = ein Ziel einer Stufe** aus `docs/roadmap.md`. Stufen nicht vorziehen.
2. **Erst Test, dann Code.** Sollwerte aus analytischen Lösungen, unabhängigen Rechnungen oder
   veröffentlichten Messungen, nie aus dem geprüften Code. Tests und Sollwerte nur mit Freigabe
   des Maintainers ändern.
3. **Konventionen sind Gesetz** (`docs/konventionen.md`): an den Grenzen die von Raytatouille und
   der Fusion.
4. **Keine Optik im eigenen Code.** Felder, Zerfallsraten und Sammeleffizienz kommen als Datenports
   aus den Lösern (ADR-Entwurf `docs/adr-entwurf-datenports.md`). Fehlt dort etwas, wird ein
   Erweiterungsvorschlag unter `docs/vorschlaege/` geschrieben.
5. **Abhängigkeiten minimal** (`docs/abhaengigkeiten.md`): NumPy und SciPy; sonst bevorzugt die
   Programme des Maintainers (PhotonPower, dort nur lesend). Neue externe Abhängigkeiten nur mit
   Freigabe.
6. **Fremder Code:** GPL- oder unklar lizenzierter Code wird nicht übernommen, nur als Vergleich
   über veröffentlichte Ergebnisse genutzt.
7. **Formeln mit Quelle** (`docs/quellen.md`); Gleichungsnummern nur aus dem Volltext. Status
   jeder Quelle festhalten (Volltext gelesen / nur zitiert).
8. **Determinismus:** Zufall nur mit explizitem Seed; Monte-Carlo-Ergebnisse mit Standardfehler.
9. **Biologische Parameter kennzeichnen:** Literaturwert, Kalibrierung oder Platzhalter; Platzhalter
   nie stillschweigend in Ergebnisse übernehmen.
10. **Keine Geheimnisse** in Dateien, Commits, Remote-URLs.
11. **Unklarheit = Frage**: Issue mit Label `frage`; Abweichung von einer Quelle: Label `paper`.

## Sprache und Stil

- Code, Kommentare, Docstrings, Commit-Messages englisch; Projektdokumentation deutsch.
- ruff (Zeilenlänge 100), mypy strict; Commits im Imperativ, ein Thema pro Commit.

## Definition of Done (PR-Checkliste)

- [ ] CI grün
- [ ] Referenztests mit unabhängigem Sollwert
- [ ] `docs/testfaelle.md` aktualisiert (Status, Messwerte)
- [ ] neue Quellen in `docs/quellen.md`
- [ ] Docstrings mit Einheiten, Vorzeichen und Verweis
- [ ] `CHANGELOG.md` unter „Unreleased“
