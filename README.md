# BindOptiX

Digitaler Zwilling optischer Biosensoren über der Zeit: Einfang von Bakterien und Molekülen,
Bindung an Rezeptoren, Assay-Abläufe und das Signal auf der Kamera.

**Rolle im Framework:** BindOptiX ist die Zeit-, Transport- und Biologieschicht über der Fusion
(`PhotonPower/RTT-Fusion`, Raytatouille). Licht rechnen die vorhandenen Löser: BioOptoMD liefert
die Rezeptorschicht und die Fluorophor-Dipole je Bindungszustand, hp-FEM und CliffordBem das
Anregungsfeld, Zerfallsraten und Quenching, Raytatouille mit Field Tracing die Abbildung durch SIL
und Relay-Optik. BindOptiX verbindet diese Karten mit Strömung, Transport (Langevin-Partikel für
Zellen, Kontinuum für Moleküle), Bindungskinetik und Kameramodell und liefert die Zielgrößen:
Nachweisgrenze, Zeit bis zum Ergebnis, Erkennung der Stämme.

**Erster Anwendungsfall:** Fluoreszenz-Biosensor für Legionellen in Wasser, gestuft bis zu
25 Subtypen (Aptamere mit Fluorophor auf Gold, später Graphit, optional plasmonische Kavitäten, Glassubstrat,
2PP-gedruckte SIL mit Spiegelrand, Relay-Optik, Kamera).

**Stand:** Spezifikation abgenommen (Z0); Abschätzung zum Einfang in
[docs/abschaetzung-einfang.md](docs/abschaetzung-einfang.md). Siehe [docs/roadmap.md](docs/roadmap.md).

## Dokumentation

- [Plan und Stufen](docs/roadmap.md)
- [Konventionen](docs/konventionen.md)
- [Testfälle](docs/testfaelle.md)
- [Abhängigkeiten](docs/abhaengigkeiten.md)
- [Quellen](docs/quellen.md)
- [ADR-Entwurf: Datenports F, B und Antwortkarte](docs/adr-entwurf-datenports.md)
- [Inventareintrag für den Fusionsplan](docs/inventar-fusionsplan.md)
- [Regeln für Programmieragenten](AGENTS.md)

## Lizenz

MIT, siehe [LICENSE](LICENSE).
