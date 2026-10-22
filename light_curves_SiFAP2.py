#!/usr/bin/env python3
"""SiFAP2: curva di luce netta con sky shift ADDITIVO via nodding.

COME SI USA
  1. Modifica i parametri nel blocco "CONFIGURAZIONE" qui sotto (e solo quello).
  2. Lancia:   python light_curves_SiFAP2.py

COSA FA (per una coppia di file .qdp target + sky, 4 colonne: t, dt, rate, err)
  1. trova gli episodi di nodding (automatico, con soglia sulla mediana, oppure manuale)
  2. esclude dai calcoli i burst (restano visibili nel plot)
  3. calcola il fattore additivo k
         k = <target | nodding> - <sky | riferimento>
     dove il riferimento e' una finestra prima/dopo il nodding (K_MODE = "window")
     oppure il cielo durante il nodding stesso (K_MODE = "simultaneous")
  4. costruisce il cielo traslato (sky + k) e la curva netta (target - sky traslato)
  5. calcola i rate medi fuori nodding (target, sky + k, netto)
  6. salva il plot in PDF e, a richiesta, la curva netta e un riepilogo in CSV

Per la magnitudine usa il rate netto con magnitude_SiFAP2.py.
"""
import csv
import os
import sys
import warnings

import numpy as np

# =============================================================================
#                              CONFIGURAZIONE
#        Cambia solo i valori in questo blocco. I tempi sono in secondi.
# =============================================================================

# ── 1. File ──────────────────────────────────────────────────────────────────
OBS_ID      = None          # nome dell'osservazione; None -> dal nome del file (es. "N3")
TARGET_FILE = "target.qdp"  # curva del rivelatore target
SKY_FILE    = None          # curva dello sky; None -> "<TARGET_FILE>_sky.qdp"
OUTPUT_DIR  = None          # dove salvare i risultati; None -> cartella del file target

# ── 2. Nodding ───────────────────────────────────────────────────────────────
NODDING = None              # None -> rilevamento automatico
                            # oppure intervalli a mano [(t_inizio, t_fine), ...]
                            #   es. [(1063, 1093), (6420, 6456)]
NODDING_THRESHOLD = 0.95    # auto: rate < soglia * mediana -> candidato nodding
NODDING_MIN_BINS  = 5       # auto: bin contigui minimi per un episodio valido

# ── 3. Burst da escludere dai calcoli (restano visibili nel plot) ────────────
BURST = []                  # es. [(4880, 5000)]; [] se non ce ne sono

# ── 4. Calcolo del fattore k ─────────────────────────────────────────────────
K_MODE = "window"           # "window"       -> sky preso in una finestra pre/post-nodding
                            # "simultaneous" -> sky preso durante il nodding stesso
                            #                   (consigliato se la curva ha un trend)
K_WINDOWS = None            # solo "window": finestra per ogni episodio, nell'ordine
                            #   es. [("post", 250), ("pre", 250)]
                            #   ("pre" = secondi prima del nodding, "post" = dopo)
K_WINDOW_DEFAULT = ("pre", 250.0)   # per gli episodi senza voce in K_WINDOWS
                                    # None -> usa tutta la baseline fuori nodding

# ── 5. Output ────────────────────────────────────────────────────────────────
START_UTC    = None         # opzionale: inizio osservazione (es. "2025-08-10 22:24:20"),
                            # scritto nella colonna "date" del riepilogo
SAVE_PLOT    = True         # salva il plot in PDF
SHOW_PLOT    = False        # mostra il plot a schermo (serve un'interfaccia grafica)
SAVE_NET_CSV = False        # salva la curva netta in CSV
SUMMARY_CSV  = None         # es. "riepilogo.csv": una riga per osservazione
PLOT_YLIM    = None         # es. (23000, 37000) per il pannello superiore; None -> automatico
PLOT_LABEL_DY     = 5.0     # altezza delle etichette nel pannello superiore (in sigma)
PLOT_LABEL_DY_NET = 3.5     # altezza dell'etichetta nel pannello della curva netta

# ── 6. Piu' osservazioni in una sola esecuzione ──────────────────────────────
# Lascia [{}] per analizzare una sola osservazione con i parametri qui sopra.
# Per piu' osservazioni elenca, per ciascuna, solo i parametri che cambiano:
#   OBSERVATIONS = [
#       dict(TARGET_FILE="oss_A.qdp", NODDING=[(1063, 1093), (6420, 6456)]),
#       dict(TARGET_FILE="oss_B.qdp", NODDING=None, BURST=[(4880, 5000)], K_MODE="simultaneous"),
#   ]
OBSERVATIONS = [{}]

# =============================================================================
#                         FINE DELLA CONFIGURAZIONE
# =============================================================================

# Nomi dei parametri di configurazione (servono per le OBSERVATIONS)
_CONFIG_KEYS = [k for k in list(globals())
                if k.isupper() and not k.startswith("_") and k != "OBSERVATIONS"]

# ─────────────────────────────────────────────────────────────────────────────
# LETTURA FILE
# ─────────────────────────────────────────────────────────────────────────────

def read_qdp(filepath):
    """Legge un .qdp a 4 colonne (t, dt, rate, err). Tiene solo bin finiti con errore > 0."""
    times, time_errs, rates, rate_errs = [], [], [], []
    with open(filepath, "r") as fh:
        for line in fh:
            line = line.strip()
            if (not line or line.startswith("!") or line.startswith("@")
                    or line.upper().startswith("READ")
                    or line.upper().startswith("NO ")):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                times.append(float(parts[0]))
                time_errs.append(float(parts[1]))
                rates.append(float(parts[2]))
                rate_errs.append(float(parts[3]))
            except ValueError:
                continue

    t, dt = np.array(times), np.array(time_errs)
    r, dr = np.array(rates), np.array(rate_errs)
    good = np.isfinite(r) & np.isfinite(dr) & (dr > 0)
    t, dt, r, dr = t[good], dt[good], r[good], dr[good]
    order = np.argsort(t, kind="stable")
    return t[order], dt[order], r[order], dr[order]


def obs_label(filepath):
    """Nome dell'osservazione dal nome del file (es. 'J1808_N3.qdp' -> 'N3')."""
    base = os.path.splitext(os.path.basename(filepath))[0]
    for p in base.split("_"):
        if p.upper().startswith("N") and p[1:].isdigit():
            return p.upper()
    return base

# ─────────────────────────────────────────────────────────────────────────────
# NODDING, MASCHERE, MEDIE
# ─────────────────────────────────────────────────────────────────────────────

def detect_nodding_auto(times, rates, threshold_fraction=0.95, min_bins=5):
    """Episodi in cui il rate resta sotto soglia * mediana per almeno min_bins bin."""
    med = np.median(rates)
    low = rates < threshold_fraction * med
    intervals, mask = [], np.zeros(len(times), dtype=bool)
    i = 0
    while i < len(low):
        if low[i]:
            j = i
            while j < len(low) and low[j]:
                j += 1
            if (j - i) >= min_bins:
                intervals.append((times[i], times[j - 1]))
                mask[i:j] = True
            i = j
        else:
            i += 1
    return intervals, mask


def build_mask(times, intervals):
    """Maschera bool da una lista di (t_start, t_end)."""
    mask = np.zeros(len(times), dtype=bool)
    for (ts, te) in intervals:
        mask |= (times >= ts) & (times <= te)
    return mask


def weighted_mean(rates, errs):
    """Media pesata (peso 1/sigma^2) e relativo errore."""
    if len(rates) == 0:
        raise ValueError("nessun bin selezionato per la media pesata")
    w = 1.0 / errs**2
    mean = np.sum(w * rates) / np.sum(w)
    err = 1.0 / np.sqrt(np.sum(w))
    return mean, err


def window_mask(t, ts, te, side, window):
    """Bin nella finestra prima ('pre') o dopo ('post') un episodio di nodding (ts, te)."""
    if side == "pre":
        return (t >= ts - window) & (t < ts)
    return (t > te) & (t <= te + window)


def check_pair(value):
    """('pre', 250) -> ('pre', 250.0); None -> None."""
    if value is None:
        return None
    side, sec = value
    side = str(side).strip().lower()
    if side not in ("pre", "post"):
        raise ValueError(f"finestra k: '{side}' non valido (usa 'pre' o 'post')")
    if float(sec) <= 0:
        raise ValueError("finestra k: i secondi devono essere > 0")
    return (side, float(sec))

# ─────────────────────────────────────────────────────────────────────────────
# ANALISI
# ─────────────────────────────────────────────────────────────────────────────

def analyze(P, obs):
    """Nodding, burst, k, curva netta e rate medi. Ritorna un dict di risultati."""
    target_path = os.path.expanduser(P["TARGET_FILE"])
    sky_path = os.path.expanduser(P["SKY_FILE"] or os.path.splitext(P["TARGET_FILE"])[0] + "_sky.qdp")
    for fp in (target_path, sky_path):
        if not os.path.isfile(fp):
            raise FileNotFoundError(f"File non trovato: {fp}")

    t_tgt, _, r_tgt, dr_tgt = read_qdp(target_path)
    t_sky, _, r_sky, dr_sky = read_qdp(sky_path)

    print(f"\nTarget : {len(t_tgt)} bin  |  t = [{t_tgt.min():.1f}, {t_tgt.max():.1f}] s")
    print(f"Sky    : {len(t_sky)} bin  |  t = [{t_sky.min():.1f}, {t_sky.max():.1f}] s")

    # ── 1. Nodding: manuale o automatico ─────────────────────────────────────
    if P["NODDING"] is not None:
        nodding_intervals = [(float(a), float(b)) for a, b in P["NODDING"]]
        mask_nod_tgt = build_mask(t_tgt, nodding_intervals)
        print(f"\nNodding MANUALE: {len(nodding_intervals)} episodi")
    else:
        nodding_intervals, mask_nod_tgt = detect_nodding_auto(
            t_tgt, r_tgt, P["NODDING_THRESHOLD"], int(P["NODDING_MIN_BINS"]))
        print(f"\nNodding AUTO: {len(nodding_intervals)} episodi  "
              f"(soglia {P['NODDING_THRESHOLD'] * 100:.0f}% mediana)")

    mask_nod_sky = build_mask(t_sky, nodding_intervals)

    for i, (ts, te) in enumerate(nodding_intervals):
        n_ep = ((t_tgt >= ts) & (t_tgt <= te)).sum()
        print(f"  [{i + 1}]  t = {ts:.1f} → {te:.1f} s  (Δt = {te - ts:.1f} s,  {n_ep} bin target)")

    if not mask_nod_tgt.any():
        raise ValueError("nessun bin di nodding trovato: indica NODDING a mano o abbassa "
                         "NODDING_THRESHOLD")

    ok_tgt = ~mask_nod_tgt
    ok_sky = ~mask_nod_sky

    # ── 2. Burst: esclusi dai calcoli, non dal plot ──────────────────────────
    bursts = [(float(a), float(b)) for a, b in (P["BURST"] or [])]
    mask_burst_tgt = build_mask(t_tgt, bursts)
    mask_burst_sky = build_mask(t_sky, bursts)
    ok_tgt = ok_tgt & ~mask_burst_tgt
    ok_sky = ok_sky & ~mask_burst_sky

    if mask_burst_tgt.any():
        print("\nBurst esclusi dai calcoli:")
        for i, (ts, te) in enumerate(bursts):
            n = ((t_tgt >= ts) & (t_tgt <= te)).sum()
            print(f"  [{i + 1}]  t = {ts:.1f} → {te:.1f} s  ({n} bin)")

    # Bin di nodding usati per k (senza eventuali burst sovrapposti)
    use_nod_tgt = mask_nod_tgt & ~mask_burst_tgt
    use_nod_sky = mask_nod_sky & ~mask_burst_sky

    # ── 3. Finestre per k (solo K_MODE = "window") ───────────────────────────
    k_mode = P["K_MODE"]
    if k_mode not in ("window", "simultaneous"):
        raise ValueError(f"K_MODE '{k_mode}' non valido (usa 'window' o 'simultaneous')")

    windows = []   # (ts, te, side, secondi) per ogni episodio con finestra
    mask_pre_sky = np.zeros(len(t_sky), dtype=bool)
    win_label = "simultaneo"

    if k_mode == "window":
        per_ep = [check_pair(v) for v in (P["K_WINDOWS"] or [])]
        default = check_pair(P["K_WINDOW_DEFAULT"])
        for i, (ts, te) in enumerate(nodding_intervals):
            pair = per_ep[i] if i < len(per_ep) else default
            if pair is None:
                continue
            side, window = pair
            ep_mask = window_mask(t_sky, ts, te, side, window) & ok_sky
            mask_pre_sky |= ep_mask
            windows.append((ts, te, side, window))
            print(f"  Finestra k ep.[{i + 1}] ({side}-nod {window:.0f} s): {ep_mask.sum()} bin sky")

        if not windows:
            win_label = "intera baseline"
        elif len({(w[2], w[3]) for w in windows}) == 1:
            win_label = f"{windows[0][2]}-nod {windows[0][3]:.0f} s"
        else:
            win_label = "config mista"

        if not windows or mask_pre_sky.sum() < 3:
            if windows:
                print("  ⚠  Troppo pochi bin nelle finestre; uso intera baseline.")
            mask_pre_sky = ok_sky.copy()
            win_label = "intera baseline"
            windows = []

    # ── 4. Offset k ──────────────────────────────────────────────────────────
    mean_tgt_nod, err_tgt_nod = weighted_mean(r_tgt[use_nod_tgt], dr_tgt[use_nod_tgt])
    if k_mode == "simultaneous":
        mean_sky_ref, err_sky_ref = weighted_mean(r_sky[use_nod_sky], dr_sky[use_nod_sky])
        k_ref_label = "sky | nodding (simultaneo)"
    else:
        mean_sky_ref, err_sky_ref = weighted_mean(r_sky[mask_pre_sky], dr_sky[mask_pre_sky])
        k_ref_label = f"sky | {win_label}"

    k = mean_tgt_nod - mean_sky_ref
    err_k = np.sqrt(err_tgt_nod**2 + err_sky_ref**2)

    r_sky_sh = r_sky + k
    dr_sky_sh = dr_sky

    print("\n── Offset additivo k ────────────────────────────────────────")
    print(f"  <Target | nodding>   = {mean_tgt_nod:.2f} ± {err_tgt_nod:.2f}  ct/s")
    print(f"  <{k_ref_label}> = {mean_sky_ref:.2f} ± {err_sky_ref:.2f}  ct/s")
    print(f"  k                    = {k:.2f} ± {err_k:.2f}  ct/s")

    # ── 5. Curva netta e rate medi fuori nodding ─────────────────────────────
    mean_tgt_out, err_tgt_out = weighted_mean(r_tgt[ok_tgt], dr_tgt[ok_tgt])
    mean_sky_sh_out, err_sky_sh_out = weighted_mean(r_sky_sh[ok_sky], dr_sky_sh[ok_sky])

    if len(t_tgt) == len(t_sky) and np.allclose(t_tgt, t_sky, atol=1e-3):
        r_sky_on_tgt, dr_sky_on_tgt = r_sky_sh, dr_sky_sh
    else:   # griglie temporali diverse: interpola lo sky traslato sui tempi del target
        r_sky_on_tgt = np.interp(t_tgt, t_sky, r_sky_sh)
        dr_sky_on_tgt = np.interp(t_tgt, t_sky, dr_sky_sh)
    r_net = r_tgt - r_sky_on_tgt
    dr_net = np.sqrt(dr_tgt**2 + dr_sky_on_tgt**2)

    mean_net, err_net = weighted_mean(r_net[ok_tgt], dr_net[ok_tgt])
    t_exp = float(t_tgt[-1] - t_tgt[0])

    print("\n── Rate medi fuori nodding ──────────────────────────────────")
    print(f"  R_target          = {mean_tgt_out:.2f} ± {err_tgt_out:.2f}  ct/s")
    print(f"  R_sky + k         = {mean_sky_sh_out:.2f} ± {err_sky_sh_out:.2f}  ct/s")
    print(f"  R_net (tgt−sky−k) = {mean_net:.2f} ± {err_net:.2f}  ct/s")

    return {
        # le prime colonne sono quelle richieste da magnitude_SiFAP2.py --obs-file
        "id": obs, "date": P["START_UTC"] or "", "duration": t_exp,
        "rate": float(mean_net), "rate_err": float(err_net),
        "n_nodding": len(nodding_intervals), "k_mode": k_mode,
        "k": float(k), "err_k": float(err_k),
        "R_target": float(mean_tgt_out), "err_R_target": float(err_tgt_out),
        "R_sky_shifted": float(mean_sky_sh_out), "err_R_sky_shifted": float(err_sky_sh_out),
        "_k_ref_label": k_ref_label, "_target_path": target_path,
        "_nodding": nodding_intervals, "_bursts": bursts, "_windows": windows,
        "_arr": dict(t_tgt=t_tgt, r_tgt=r_tgt, dr_tgt=dr_tgt, t_sky=t_sky, r_sky=r_sky,
                     dr_sky=dr_sky, r_sky_sh=r_sky_sh, dr_sky_sh=dr_sky_sh,
                     r_sky_on_tgt=r_sky_on_tgt, r_net=r_net, dr_net=dr_net,
                     ok_tgt=ok_tgt, ok_sky=ok_sky, mask_nod_tgt=mask_nod_tgt,
                     mask_burst_tgt=mask_burst_tgt),
    }

# ─────────────────────────────────────────────────────────────────────────────
# PLOT E SALVATAGGI
# ─────────────────────────────────────────────────────────────────────────────

def curve_label(ax, t, r, mask_ok, text, color, frac=0.2, dy_sigma=4.0):
    t_ok, r_ok = t[mask_ok], r[mask_ok]
    if len(t_ok) == 0:
        return
    idx = max(0, min(int(frac * len(t_ok)), len(t_ok) - 1))
    i0, i1 = max(0, idx - 25), min(len(r_ok), idx + 25)
    y0 = np.median(r_ok[i0:i1])
    sig = np.std(r_ok[i0:i1])
    ax.text(t_ok[idx], y0 + dy_sigma * sig, text, color=color, fontsize=8.5, fontweight="bold",
            ha="center", va="bottom", zorder=10,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.6))


def make_plot(res, P, outdir):
    import matplotlib
    if not P["SHOW_PLOT"]:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    a = res["_arr"]
    obs, k = res["id"], res["k"]
    dy = float(P["PLOT_LABEL_DY"])

    fig, (ax_raw, ax_net) = plt.subplots(
        2, 1, figsize=(12, 7), sharex=True,
        gridspec_kw={"height_ratios": [2, 1], "hspace": 0.05})

    C_TGT, C_SKY, C_SH, C_NET, C_BURST = "#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e"
    ek = dict(capsize=0, elinewidth=0.5)

    # Finestre usate per k (giallo, sotto al nodding in grigio)
    for (ts, te, side, window) in res["_windows"]:
        if side == "pre":
            x0, x1 = max(a["t_tgt"].min(), ts - window), ts
        else:
            x0, x1 = te, min(a["t_tgt"].max(), te + window)
        ax_raw.axvspan(x0, x1, color="lightyellow", alpha=0.7, zorder=0)

    for (ts, te) in res["_nodding"]:
        ax_raw.axvspan(ts, te, color="lightgray", alpha=0.55, zorder=1)
        ax_net.axvspan(ts, te, color="lightgray", alpha=0.55, zorder=1)

    for (ts, te) in res["_bursts"]:
        for ax in (ax_raw, ax_net):
            ax.axvspan(ts, te, color=C_BURST, alpha=0.25, zorder=1)
            ax.axvline(ts, color=C_BURST, lw=0.8, ls="--", zorder=2)
            ax.axvline(te, color=C_BURST, lw=0.8, ls="--", zorder=2)

    ax_raw.errorbar(a["t_tgt"], a["r_tgt"], yerr=a["dr_tgt"], fmt=".", ms=2, lw=0,
                    color=C_TGT, ecolor=C_TGT, alpha=0.8, zorder=3, **ek)
    ax_raw.errorbar(a["t_sky"], a["r_sky"], yerr=a["dr_sky"], fmt=".", ms=2, lw=0,
                    color=C_SKY, ecolor=C_SKY, alpha=0.7, zorder=2, **ek)
    ax_raw.errorbar(a["t_sky"], a["r_sky_sh"], yerr=a["dr_sky_sh"], fmt=".", ms=2, lw=0,
                    color=C_SH, ecolor=C_SH, alpha=0.75, zorder=2, **ek)

    curve_label(ax_raw, a["t_tgt"], a["r_tgt"], a["ok_tgt"], "Target detector", C_TGT, dy_sigma=dy)
    curve_label(ax_raw, a["t_sky"], a["r_sky"], a["ok_sky"], "Sky detector", C_SKY, dy_sigma=dy)
    curve_label(ax_raw, a["t_sky"], a["r_sky_sh"], a["ok_sky"],
                f"Shifted sky  (k = {k:.0f} ct/s)", C_SH, dy_sigma=dy)

    ax_raw.set_ylabel(r"Count rate (ct s$^{-1}$)", fontsize=10)
    if P["PLOT_YLIM"]:
        ax_raw.set_ylim(*P["PLOT_YLIM"])
    ax_raw.set_title(f"{obs} light curve", fontsize=12, fontweight="normal")
    ax_raw.tick_params(which="both", direction="in", top=True, right=True)

    ax_net.errorbar(a["t_tgt"], a["r_net"], yerr=a["dr_net"], fmt=".", ms=2, lw=0,
                    color=C_NET, ecolor=C_NET, alpha=0.8, zorder=3, **ek)
    curve_label(ax_net, a["t_tgt"], a["r_net"], a["ok_tgt"], "Net  (target - shifted sky)",
                C_NET, dy_sigma=float(P["PLOT_LABEL_DY_NET"]))
    ax_net.set_xlabel("Time (s)", fontsize=10)
    ax_net.set_ylabel(r"Net count rate (ct s$^{-1}$)", fontsize=10)
    ax_net.tick_params(which="both", direction="in", top=True, right=True)

    fig.align_ylabels([ax_raw, ax_net])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.tight_layout()

    if P["SAVE_PLOT"]:
        out_pdf = os.path.join(outdir, f"LC_SiFAP2_{obs}_scaled.pdf")
        fig.savefig(out_pdf, dpi=150, bbox_inches="tight")
        print(f"\nPlot salvato: {out_pdf}")
    if P["SHOW_PLOT"]:
        plt.show()
    plt.close(fig)


def save_net_csv(res, outdir):
    a = res["_arr"]
    path = os.path.join(outdir, f"LC_SiFAP2_{res['id']}_net.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["time_s", "rate_target", "err_target", "rate_sky_shifted",
                    "rate_net", "err_net", "in_nodding", "in_burst"])
        for i in range(len(a["t_tgt"])):
            w.writerow([f"{a['t_tgt'][i]:.4f}", f"{a['r_tgt'][i]:.4f}", f"{a['dr_tgt'][i]:.4f}",
                        f"{a['r_sky_on_tgt'][i]:.4f}", f"{a['r_net'][i]:.4f}", f"{a['dr_net'][i]:.4f}",
                        int(a["mask_nod_tgt"][i]), int(a["mask_burst_tgt"][i])])
    print(f"Curva netta salvata: {path}")


def print_summary(res):
    print("\n" + "=" * 65)
    print(f"RIEPILOGO  [{res['id']}]")
    print("=" * 65)
    print(f"  k  ({res['_k_ref_label']})  = {res['k']:.2f} ± {res['err_k']:.2f}  ct/s")
    print(f"  R_target  (fuori nod) = {res['R_target']:.2f} ± {res['err_R_target']:.2f}  ct/s")
    print(f"  R_sky+k   (fuori nod) = {res['R_sky_shifted']:.2f} ± {res['err_R_sky_shifted']:.2f}  ct/s")
    print(f"  R_net     (fuori nod) = {res['rate']:.2f} ± {res['rate_err']:.2f}  ct/s")
    print(f"  T_exp = {res['duration']:.1f} s")
    print("-" * 65)
    print("  Per la magnitudine (magnitude_SiFAP2.py):")
    print(f"    --rate {res['rate']:.2f} --rate-err {res['rate_err']:.2f} --duration {res['duration']:.1f}")
    print("=" * 65)


def write_summary(results, path):
    rows = [{k: v for k, v in r.items() if not k.startswith("_")} for r in results]
    fields = []
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    with open(path, "w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=fields)
        wr.writeheader()
        wr.writerows(rows)
    print(f"\nRiepilogo salvato in {path}")

# ─────────────────────────────────────────────────────────────────────────────
# ESECUZIONE
# ─────────────────────────────────────────────────────────────────────────────

def run_observation(P):
    obs = P["OBS_ID"] or obs_label(P["TARGET_FILE"])
    print("=" * 65)
    print(f"SiFAP2 Light Curve Analyzer — {obs}  (sky shift additivo)")
    print("=" * 65)

    res = analyze(P, obs)

    outdir = P["OUTPUT_DIR"] or os.path.dirname(os.path.abspath(res["_target_path"]))
    os.makedirs(outdir, exist_ok=True)
    if P["SAVE_PLOT"] or P["SHOW_PLOT"]:
        make_plot(res, P, outdir)
    if P["SAVE_NET_CSV"]:
        save_net_csv(res, outdir)

    print_summary(res)
    return res


def main():
    defaults = {k: globals()[k] for k in _CONFIG_KEYS}
    observations = OBSERVATIONS or [{}]

    results, failed = [], []
    for override in observations:
        unknown = [k for k in override if k not in defaults]
        if unknown:
            sys.exit(f"OBSERVATIONS: parametri sconosciuti: {', '.join(unknown)}\n"
                     f"Parametri validi: {', '.join(sorted(defaults))}")
        P = {**defaults, **override}
        label = P["OBS_ID"] or obs_label(P["TARGET_FILE"])
        try:
            results.append(run_observation(P))
        except (FileNotFoundError, ValueError) as exc:
            print(f"\n  ✗ [{label}] {exc}\n")
            failed.append(label)

    if defaults["SUMMARY_CSV"] and results:
        write_summary(results, defaults["SUMMARY_CSV"])
    if failed:
        print(f"\nOsservazioni non completate: {', '.join(failed)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
