# Konventionen

Verbindlich; an allen Grenzen die von Raytatouille und der Fusion (RTT-Fusion
`docs/konventionen.md`). Intern gelten die folgenden Festlegungen.

- **Einheiten intern:** Längen in µm, Zeit in s, Temperatur in K, Leistung in mW. Im Austausch mit
  der Fusion Längen in mm; interne Größen im Namen kenntlich (`_um`, `_s`).
- **Stoffmengen als Teilchenzahlen:** Konzentrationen in 1/µm³, Belegungen in 1/µm². So gelten
  dieselben Größen für Moleküle (1 nM ≈ 0,602 /µm³) und für Zellen (1 KBE/mL = 10⁻¹² /µm³).
  Umrechnungen nur über `bindoptix.units`; ein Test prüft die Konstanten.
- **Kinetik:** k_on in µm³/s je Bindungsplatz, k_off in 1/s; Literaturwerte in 1/(M·s) werden
  beim Einlesen umgerechnet.
- **Sensorsystem:** rechtshändig, Sensoroberfläche bei z = 0, Lösung bei z > 0, Normale +z, wie
  BioOptoMD. Die Pose ins RTT-System ist explizit. Beleuchtung von der Lösungsseite ist eine
  **Spiegelung z → −z** (Determinante −1), keine Drehung (Befund K4); Dipole transformieren als
  polare Vektoren.
- **Schwerkraft:** als Vektor im Sensorsystem angegeben (Einbaulage), nie implizit −z.
- **Felder** (Port F): ∝ exp(i(k·r − ωt)), Zeitfaktor e^(−iωt), normiert auf 1 mW einfallende
  Leistung am Eingangsport der Fusion. Löser mit anderem Zeitfaktor (SpeckleBem e^(+jωt)) werden
  im Adapter konjugiert. Die Löser dürfen stattdessen auf die Bestrahlungsstärke I₀ einer
  einfallenden ebenen Welle normieren (`irradiance_w_m2`, W/m², senkrecht zum Strahl im
  Einfallsmedium, I₀ = n_inc |E₀|²/(2 Z₀)); BindOptiX rechnet solche Karten mit der lokalen
  Bestrahlungsstärke am Spot je mW am Eingangsport um (`FieldMap.per_port_power`, mit Status
  L/K/P). Dazu je Karte die Einfallsrichtung `k_hat` und die Polarisation `polarisation`
  im Sensorsystem; Anregung im Aufbau von der Lösungsseite (k̂ zeigt nach −z). Die Phase ist
  auf (0, 0, 0) bezogen (`phase_reference`). Medium je Gitterpunkt als `medium`/`media`;
  Knoten auf einer Grenzfläche tragen den Grenzwert der Seite, deren Medium sie tragen, und
  über Medien hinweg wird nicht interpoliert (#23).
- **Zufall:** nur mit explizitem Seed; Monte-Carlo-Ergebnisse immer mit Standardfehler.
- **Ergebnisse:** Format `raytatouille-result` (ADR 0023) für Kennzahlen und Tabellen; Felder und
  Zeitreihen als reguläre Gitter (`.npy`, HDF5) mit Metadaten nach demselben Schema.
