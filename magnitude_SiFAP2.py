#!/usr/bin/env python3
"""SiFAP2: airmass, magnitudine calibrata e flusso ottico da riga di comando.

Per ogni osservazione calcola:
  1. airmass iniziale e medio (da data di inizio e durata)
  2. magnitudine osservata dalla retta di calibrazione, dal rate
  3. magnitudine corretta per estinzione atmosferica (k_ext * airmass)
  4. magnitudine corretta anche per estinzione interstellare (A_V)
  5. flusso ottico (erg s^-1 cm^-2)
  6. (opzionale) confronto con una magnitudine di riferimento: Delta e sigma

Retta di calibrazione (rate R in kcps):
    -2.512 * log10(R_kcps) = m * mag_obs + q
    =>  mag_obs = (-2.512 * log10(R_kcps) - q) / m

Tre modi d'uso (le opzioni da riga di comando hanno priorita' sul file di config):

  # a) singola osservazione da terminale
  python magnitude_SiFAP2.py --config config.json --date "2025-08-10 22:24:20" \
      --duration 5406.7 --rate 1234.5 --rate-err 12.3 --k-ext 0.12

  # b) piu' osservazioni da CSV (colonne: id,date,duration,rate,rate_err[,k_ext,k_ext_err])
  python magnitude_SiFAP2.py --config config.json --obs-file osservazioni.csv --output risultati.csv

  # c) interattivo: il programma chiede i dati uno alla volta
  python magnitude_SiFAP2.py --config config.json --interactive

Le date sono in UTC; il rate e' in counts/s.
"""
import argparse
import csv
import json
import sys
from dataclasses import dataclass
from typing import Optional

import numpy as np
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time, TimeDelta

# ============================================================================
# AIRMASS
# ============================================================================

# Sito di default: Telescopio Nazionale Galileo (La Palma), come nello script originale.
DEFAULT_SITE = {"lat": 28.754, "lon": -17.8890555555556, "height": 2387.2}


def make_location(lat_deg, lon_deg, height_m):
    """EarthLocation da latitudine/longitudine in gradi e quota in metri."""
    return EarthLocation(lat=lat_deg * u.deg, lon=lon_deg * u.deg, height=height_m * u.m)


def make_target(ra, dec):
    """SkyCoord da RA/DEC.

    Accetta stringhe sessagesimali ('18h08m27.647s', '-36d58m43.90s')
    oppure gradi decimali (numeri o stringhe numeriche).
    """
    try:
        return SkyCoord(float(ra) * u.deg, float(dec) * u.deg)
    except (TypeError, ValueError):
        return SkyCoord(str(ra), str(dec))


def airmass_from_altitude(alt_deg, model="secz"):
    """Airmass da altezza in gradi.

    model = 'secz'  : sec(z), come nello script originale (buona per z < ~60 gradi)
    model = 'kasten': Kasten & Young (1989), piu' accurata a bassa altezza
    """
    alt = np.asarray(alt_deg, dtype=float)
    if model == "secz":
        return 1.0 / np.sin(np.radians(alt))
    if model == "kasten":
        return 1.0 / (np.sin(np.radians(alt)) + 0.50572 * (alt + 6.07995) ** -1.6364)
    raise ValueError(f"Modello di airmass sconosciuto: {model!r} (usa 'secz' o 'kasten')")


def compute_airmass(coord, location, start_utc, duration_s, step_s=30.0, model="secz"):
    """Airmass all'inizio e medio su un'osservazione.

    Parametri
    ---------
    coord       : SkyCoord del target
    location    : EarthLocation del sito
    start_utc   : data/ora di inizio in UTC (stringa 'YYYY-MM-DD HH:MM:SS' o astropy Time)
    duration_s  : durata dell'osservazione in secondi
    step_s      : passo di campionamento per la media (default 30 s)
    model       : 'secz' (default) o 'kasten'

    Ritorna un dict con airmass_start, airmass_mean, alt_start, alt_mean (gradi).
    """
    t0 = start_utc if isinstance(start_utc, Time) else Time(start_utc, scale="utc")
    n_steps = max(int(duration_s // step_s), 1)
    times = t0 + TimeDelta(step_s * u.s) * np.arange(n_steps + 1)

    altaz = coord.transform_to(AltAz(location=location, obstime=times))
    alt = altaz.alt.deg
    X = airmass_from_altitude(alt, model=model)

    return {
        "airmass_start": float(X[0]),
        "airmass_mean": float(np.mean(X)),
        "alt_start": float(alt[0]),
        "alt_mean": float(np.mean(alt)),
    }

# ============================================================================
# FOTOMETRIA
# ============================================================================

POGSON = 2.512
C_ANG = 2.99792458e18  # velocita' della luce in Angstrom/s

# Parametri di banda (default: Johnson-Cousins/Bessell). Verifica che siano
# quelli che vuoi usare per il tuo lavoro prima di pubblicare risultati.
#   F_nu0        : flusso di punto zero in erg s^-1 cm^-2 Hz^-1 (3640 Jy = 3.640e-20)
#   lambda_eff   : lunghezza d'onda effettiva in Angstrom
#   delta_lambda : larghezza di banda in Angstrom
BANDS = {
    "V": {"F_nu0": 3.640e-20, "lambda_eff": 5500.0, "delta_lambda": 890.0},
    "R": {"F_nu0": 3.080e-20, "lambda_eff": 6400.0, "delta_lambda": 1580.0},
    "B": {"F_nu0": 4.260e-20, "lambda_eff": 4400.0, "delta_lambda": 890.0},
}


@dataclass
class Calibration:
    """Retta di calibrazione e correzioni."""
    m: float
    q: float
    err_m: float = 0.0
    err_q: float = 0.0
    A_V: float = 0.0       # estinzione interstellare (mag)
    err_AV: float = 0.0
    # Riferimento opzionale (es. fotometria LCO) per il confronto
    mag_ref: Optional[float] = None
    err_mag_ref: float = 0.0


def compute_mag(R, eR, airmass, k_ext, cal, err_k_ext=0.0, label=""):
    """Magnitudine corretta per estinzione atmosferica e interstellare.

    R, eR     : rate e relativo errore in counts/s (cps); internamente convertiti in kcps
    airmass   : airmass (iniziale o medio, a scelta dell'utente)
    k_ext     : coefficiente di estinzione atmosferica (mag/airmass)
    cal       : Calibration

    Ritorna un dict con mag_obs, mag_atm (corretta solo per atmosfera),
    mag_corr (corretta anche per A_V), err_mag e, se presente il riferimento,
    delta / err_tot / sigma. Se R <= 0 ritorna None.
    """
    if R <= 0:
        print(f"  [!] [{label}] R = {R:.2f} <= 0, magnitudine non definita.")
        return None

    R_kcps = R / 1e3
    eR_kcps = eR / 1e3
    log_R = np.log10(R_kcps)

    mag_obs = (-POGSON * log_R - cal.q) / cal.m
    mag_atm = mag_obs - k_ext * airmass          # corretta per estinzione atmosferica
    mag_corr = mag_atm - cal.A_V                 # + correzione interstellare

    dfdR = -POGSON / (np.log(10) * R_kcps * cal.m)
    dfdm = -(-POGSON * log_R - cal.q) / cal.m**2
    dfdq = -1.0 / cal.m
    err_mag = np.sqrt(
        (dfdR * eR_kcps) ** 2
        + (dfdm * cal.err_m) ** 2
        + (dfdq * cal.err_q) ** 2
        + (airmass * err_k_ext) ** 2
        + cal.err_AV ** 2
    )

    out = {
        "mag_obs": float(mag_obs),
        "mag_atm": float(mag_atm),
        "mag_corr": float(mag_corr),
        "err_mag": float(err_mag),
    }
    if cal.mag_ref is not None:
        delta = mag_corr - cal.mag_ref
        err_tot = np.sqrt(err_mag**2 + cal.err_mag_ref**2)
        out.update(delta=float(delta), err_tot=float(err_tot), sigma=float(abs(delta) / err_tot))
    return out


def optical_flux(mag, smag, band="V"):
    """Flusso ottico (erg s^-1 cm^-2) e relativo errore dalla magnitudine.

    F = F_nu0 * 10^(-0.4 mag) * c / lambda^2 * delta_lambda
    (flusso in banda; formula identica a quella del codice originale)
    """
    if band not in BANDS:
        raise ValueError(f"Banda {band!r} non definita. Disponibili: {', '.join(BANDS)}")
    b = BANDS[band]
    F = b["F_nu0"] * 10 ** (-0.4 * mag) * C_ANG / b["lambda_eff"] ** 2 * b["delta_lambda"]
    sF = 0.4 * np.log(10) * F * smag
    return F, sF

# ============================================================================
# CONFIG
# ============================================================================

def load_config(path):
    if not path:
        return {}
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def get(args_val, cfg, section, key, default=None):
    """Priorita': argomento da terminale > config > default."""
    if args_val is not None:
        return args_val
    return cfg.get(section, {}).get(key, default)


def build_setup(args, cfg):
    """Raccoglie target, sito, calibrazione e opzioni dai flag e dal config."""
    setup = {
        "ra": get(args.ra, cfg, "target", "ra"),
        "dec": get(args.dec, cfg, "target", "dec"),
        "lat": get(args.lat, cfg, "site", "lat", DEFAULT_SITE["lat"]),
        "lon": get(args.lon, cfg, "site", "lon", DEFAULT_SITE["lon"]),
        "height": get(args.height, cfg, "site", "height", DEFAULT_SITE["height"]),
        "m": get(args.m, cfg, "calibration", "m"),
        "q": get(args.q, cfg, "calibration", "q"),
        "err_m": get(args.err_m, cfg, "calibration", "err_m", 0.0),
        "err_q": get(args.err_q, cfg, "calibration", "err_q", 0.0),
        "A_V": get(args.av, cfg, "calibration", "A_V", 0.0),
        "err_AV": get(args.err_av, cfg, "calibration", "err_AV", 0.0),
        "mag_ref": get(args.mag_ref, cfg, "reference", "mag"),
        "err_mag_ref": get(args.err_mag_ref, cfg, "reference", "err", 0.0),
        "k_ext": get(args.k_ext, cfg, "extinction", "k_ext"),
        "err_k_ext": get(args.err_k_ext, cfg, "extinction", "err_k_ext", 0.0),
        "band": get(args.band, cfg, "options", "band", "V"),
        "model": get(args.airmass_model, cfg, "options", "airmass_model", "secz"),
        "use": get(args.airmass_use, cfg, "options", "airmass_use", "mean"),
        "step": get(args.step, cfg, "options", "step", 30.0),
    }
    missing = [k for k in ("ra", "dec", "m", "q") if setup[k] is None]
    if missing:
        sys.exit(
            "Mancano parametri obbligatori: " + ", ".join(missing)
            + "\nIndicali nel file di config (--config) o come opzioni da terminale."
        )
    return setup

# ============================================================================
# INPUT DELLE OSSERVAZIONI
# ============================================================================

def _num(s):
    return float(s) if s not in (None, "") else None


def observations_from_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        for i, row in enumerate(csv.DictReader(f), start=1):
            yield {
                "id": row.get("id") or f"obs{i}",
                "date": row["date"],
                "duration": float(row["duration"]),
                "rate": float(row["rate"]),
                "rate_err": float(row["rate_err"]),
                "k_ext": _num(row.get("k_ext")),
                "err_k_ext": _num(row.get("k_ext_err")),
            }


def _ask(prompt, cast=str, default=None, required=True):
    suffix = f" [{default}]" if default is not None else ""
    while True:
        raw = input(f"{prompt}{suffix}: ").strip()
        if raw == "" and default is not None:
            return default
        if raw == "" and not required:
            return None
        if raw == "":
            continue
        try:
            return cast(raw)
        except ValueError:
            print("  Valore non valido, riprova.")


def observations_interactive(setup):
    print("Inserimento interattivo (invio vuoto sull'ID per terminare). Date in UTC.\n")
    n = 0
    while True:
        n += 1
        oid = input(f"ID osservazione [obs{n}] (vuoto = fine): ").strip()
        if oid == "" and n > 1:
            return
        oid = oid or f"obs{n}"
        yield {
            "id": oid,
            "date": _ask("  Data inizio UTC (YYYY-MM-DD HH:MM:SS)"),
            "duration": _ask("  Durata (s)", float),
            "rate": _ask("  Rate (cps)", float),
            "rate_err": _ask("  Errore sul rate (cps)", float),
            "k_ext": _ask("  Coefficiente di estinzione k (mag/airmass)", float,
                          default=setup["k_ext"]),
            "err_k_ext": setup["err_k_ext"],
        }

# ============================================================================
# CALCOLO E OUTPUT
# ============================================================================

def process(obs, setup, coord, location, cal):
    am = compute_airmass(coord, location, obs["date"], obs["duration"],
                         step_s=setup["step"], model=setup["model"])
    X = am["airmass_mean"] if setup["use"] == "mean" else am["airmass_start"]

    k = obs["k_ext"] if obs.get("k_ext") is not None else setup["k_ext"]
    if k is None:
        sys.exit(f"[{obs['id']}] Manca il coefficiente di estinzione k_ext "
                 "(--k-ext, config, o colonna k_ext nel CSV).")
    ek = obs["err_k_ext"] if obs.get("err_k_ext") is not None else setup["err_k_ext"]

    res = compute_mag(obs["rate"], obs["rate_err"], X, k, cal, err_k_ext=ek, label=obs["id"])
    row = {
        "id": obs["id"], "date": obs["date"], "duration_s": obs["duration"],
        "rate_cps": obs["rate"], "rate_err_cps": obs["rate_err"],
        "k_ext": k, "airmass_start": am["airmass_start"], "airmass_mean": am["airmass_mean"],
        "alt_start_deg": am["alt_start"], "airmass_used": X,
    }
    if res is None:
        return row
    row.update({k_: res[k_] for k_ in ("mag_obs", "mag_atm", "mag_corr", "err_mag")})
    F, sF = optical_flux(res["mag_corr"], res["err_mag"], band=setup["band"])
    row.update(flux=F, flux_err=sF)
    for key in ("delta", "err_tot", "sigma"):
        if key in res:
            row[key] = res[key]
    return row


def print_table(rows, with_ref):
    head = (f"{'ID':<6}| {'Inizio (UTC)':<19} | {'Dur(s)':>8} | {'X ini':>6} | {'X med':>6} "
            f"| {'m_obs':>7} | {'m_atm':>7} | {'m_corr':>7} | {'err':>6} | {'Flusso':>10} | {'err F':>10}")
    if with_ref:
        head += f" | {'Delta':>7} | {'sigma':>5}"
    print(head)
    print("-" * len(head))
    for r in rows:
        if "mag_corr" not in r:
            print(f"{r['id']:<6}| {r['date']:<19} | {r['duration_s']:>8.1f} | {r['airmass_start']:>6.3f} "
                  f"| {r['airmass_mean']:>6.3f} | (magnitudine non definita)")
            continue
        line = (f"{r['id']:<6}| {r['date']:<19} | {r['duration_s']:>8.1f} | {r['airmass_start']:>6.3f} "
                f"| {r['airmass_mean']:>6.3f} | {r['mag_obs']:>7.3f} | {r['mag_atm']:>7.3f} "
                f"| {r['mag_corr']:>7.3f} | {r['err_mag']:>6.3f} | {r['flux']:>10.3e} | {r['flux_err']:>10.3e}")
        if with_ref and "delta" in r:
            line += f" | {r['delta']:>7.3f} | {r['sigma']:>5.2f}"
        print(line)


def write_csv(rows, path):
    fields = []
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"\nRisultati salvati in {path}")

# ============================================================================
# CLI
# ============================================================================

def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", help="file JSON con target, sito, calibrazione (vedi config_example.json)")
    p.add_argument("--obs-file", help="CSV con le osservazioni")
    p.add_argument("--interactive", action="store_true", help="chiede i dati da terminale")
    p.add_argument("--output", help="salva i risultati in CSV")

    g = p.add_argument_group("singola osservazione")
    g.add_argument("--id", default="obs1")
    g.add_argument("--date", help="inizio osservazione, UTC, 'YYYY-MM-DD HH:MM:SS'")
    g.add_argument("--duration", type=float, help="durata (s)")
    g.add_argument("--rate", type=float, help="rate (cps)")
    g.add_argument("--rate-err", type=float, help="errore sul rate (cps)")

    g = p.add_argument_group("estinzione")
    g.add_argument("--k-ext", type=float, help="coefficiente di estinzione atmosferica (mag/airmass)")
    g.add_argument("--err-k-ext", type=float, help="errore su k_ext")
    g.add_argument("--av", type=float, help="estinzione interstellare A_V (mag)")
    g.add_argument("--err-av", type=float, help="errore su A_V")

    g = p.add_argument_group("calibrazione (-2.512*log10(R_kcps) = m*mag + q)")
    g.add_argument("--m", type=float)
    g.add_argument("--q", type=float)
    g.add_argument("--err-m", type=float)
    g.add_argument("--err-q", type=float)

    g = p.add_argument_group("target e sito")
    g.add_argument("--ra", help="es. '18h08m27.647s' oppure gradi")
    g.add_argument("--dec", help="es. '-36d58m43.90s' oppure gradi")
    g.add_argument("--lat", type=float, help="latitudine (gradi)")
    g.add_argument("--lon", type=float, help="longitudine (gradi, Est positiva)")
    g.add_argument("--height", type=float, help="quota (m)")

    g = p.add_argument_group("riferimento (opzionale, es. LCO)")
    g.add_argument("--mag-ref", type=float)
    g.add_argument("--err-mag-ref", type=float)

    g = p.add_argument_group("opzioni")
    g.add_argument("--band", choices=list(BANDS), help="banda per il flusso (default V)")
    g.add_argument("--airmass-model", choices=["secz", "kasten"], help="default secz")
    g.add_argument("--airmass-use", choices=["mean", "start"],
                   help="airmass usato per la correzione: medio (default) o iniziale")
    g.add_argument("--step", type=float, help="passo di campionamento airmass in secondi (default 30)")
    return p.parse_args()


def main():
    args = parse_args()
    cfg = load_config(args.config)
    setup = build_setup(args, cfg)

    coord = make_target(setup["ra"], setup["dec"])
    location = make_location(setup["lat"], setup["lon"], setup["height"])
    cal = Calibration(m=setup["m"], q=setup["q"], err_m=setup["err_m"], err_q=setup["err_q"],
                      A_V=setup["A_V"], err_AV=setup["err_AV"],
                      mag_ref=setup["mag_ref"], err_mag_ref=setup["err_mag_ref"])

    if args.obs_file:
        observations = list(observations_from_csv(args.obs_file))
    elif args.interactive:
        observations = list(observations_interactive(setup))
    elif args.date:
        if None in (args.duration, args.rate, args.rate_err):
            sys.exit("Per una singola osservazione servono --date, --duration, --rate, --rate-err.")
        observations = [{"id": args.id, "date": args.date, "duration": args.duration,
                         "rate": args.rate, "rate_err": args.rate_err,
                         "k_ext": None, "err_k_ext": None}]
    else:
        sys.exit("Specifica --date (singola), --obs-file (CSV) oppure --interactive. Usa -h per l'aiuto.")

    rows = [process(o, setup, coord, location, cal) for o in observations]
    print()
    print_table(rows, with_ref=setup["mag_ref"] is not None)
    print(f"\nBanda: {setup['band']} | modello airmass: {setup['model']} | "
          f"airmass usato per la correzione: {setup['use']}")
    if args.output:
        write_csv(rows, args.output)


if __name__ == "__main__":
    main()
