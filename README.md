# SiFAP2 tools

Raccolta di strumenti per l'analisi e la pianificazione delle osservazioni con SiFAP2.

| Strumento | A cosa serve |
|---|---|
| [`magnitude_SiFAP2.py`](#1-magnitude_sifap2py--airmass-magnitudine-e-flusso) | da data, durata e rate: airmass, magnitudine calibrata e corretta, flusso ottico |
| [`tempo_obs_SiFAP2_sens.py`](#2-tempo_obs_sifap2_senspy--tempo-di-osservazione-necessario) | da magnitudine e ampiezza attesa: tempo di osservazione necessario per rivelare un segnale pulsato |
| [`light_curves_SiFAP2.py`](#3-light_curves_sifap2py--nodding-burst-k-e-curva-netta) | da file `.qdp` target + sky: nodding, burst, fattore k, curva netta e rate medi |

## Requisiti

```
pip install astropy numpy scipy matplotlib
```

Ogni script è autonomo e ha bisogno solo di alcune librerie: `astropy` serve a
`magnitude_SiFAP2.py`, `scipy` a `tempo_obs_SiFAP2_sens.py`, `matplotlib` ai plot di
`light_curves_SiFAP2.py` (che per il resto usa solo `numpy`).

## File

| File | Contenuto |
|---|---|
| `magnitude_SiFAP2.py` | airmass, magnitudine calibrata e flusso ottico, in un unico script |
| `config_example.json` | target, sito, calibrazione: **copialo e compilalo** |
| `osservazioni_esempio.csv` | esempio di input multiplo |
| `tempo_obs_SiFAP2_sens.py` | stima interattiva del tempo di osservazione |
| `light_curves_SiFAP2.py` | analisi delle curve di luce target + sky: nodding, burst, k, curva netta, rate |

## Convenzione di calibrazione (comune ai due strumenti)

Entrambi usano la retta di calibrazione, con il rate R espresso in **kcps**:

```
-2.512 * log10(R_kcps) = m * mag + q
```

- `magnitude_SiFAP2.py` la inverte per ottenere la magnitudine dal rate (`m` e `q` vengono dal tuo
  file di config).
- `tempo_obs_SiFAP2_sens.py` la usa nel verso opposto, dalla magnitudine al rate (con `m` e `q`
  scritti nello script, vedi sotto).

---

## 1. `magnitude_SiFAP2.py` — airmass, magnitudine e flusso

Per ogni osservazione, a partire da data, durata e rate, calcola:

1. airmass iniziale e medio
2. magnitudine osservata dalla retta di calibrazione
3. magnitudine corretta per estinzione atmosferica (`k_ext * airmass`)
4. magnitudine corretta anche per estinzione interstellare (`A_V`)
5. flusso ottico (erg s⁻¹ cm⁻²)
6. (opzionale) confronto con una magnitudine di riferimento (es. LCO): Δ e σ

### Configurazione

Copia `config_example.json`, inserisci target, sito e retta di calibrazione
(`m`, `q` e relativi errori, `A_V`). `m` e `q` sono obbligatori e sono lasciati a `null`
di proposito: i valori dipendono dal tuo strumento e dalla tua calibrazione.

Le opzioni da terminale hanno priorità sul file di config.

### Uso

Singola osservazione:

```
python magnitude_SiFAP2.py --config config.json \
    --date "2025-08-10 22:24:20" --duration 5406.7 \
    --rate 1234.5 --rate-err 12.3 --k-ext 0.12
```

Più osservazioni da CSV (`id,date,duration,rate,rate_err[,k_ext,k_ext_err]`):

```
python magnitude_SiFAP2.py --config config.json --obs-file osservazioni.csv --output risultati.csv
```

Interattivo:

```
python magnitude_SiFAP2.py --config config.json --interactive
```

Tutte le date sono in **UTC**; il rate è in **counts/s**.

### Opzioni utili

- `--airmass-use mean|start`: airmass usato per la correzione (default: medio)
- `--airmass-model secz|kasten`: `secz` riproduce lo script originale; a bassa altezza
  (come per bersagli australi da La Palma) Kasten & Young è più accurato
- `--band V|R|B`: banda per il flusso (default V)
- `--mag-ref`, `--err-mag-ref`: magnitudine di riferimento per Δ e σ
- `--err-k-ext`: errore su k, propagato nell'errore sulla magnitudine

### Note

- I parametri di banda in `magnitude_SiFAP2.py` (dizionario `BANDS`) sono valori standard
  Johnson/Bessell: verifica che siano quelli giusti per il tuo lavoro.
- Lo script è diviso in sezioni (airmass, fotometria, config, input, CLI). Le funzioni si possono
  anche importare in altri codici:
  `from magnitude_SiFAP2 import compute_airmass, compute_mag, optical_flux, Calibration`.
- Il sito di default è il TNG (La Palma); cambialo con `--lat/--lon/--height` o nel config.

---

## 2. `tempo_obs_SiFAP2_sens.py` — tempo di osservazione necessario

Stima per quanto tempo bisogna osservare una sorgente di data magnitudine V per rivelare una
modulazione periodica di data ampiezza, con una probabilità di falso allarme e di falsa
esclusione fissate. Si usa **prima** di proporre o eseguire un'osservazione.

```
python tempo_obs_SiFAP2_sens.py
```

Il programma è interattivo e chiede, in ordine:

| Domanda | Significato | Default |
|---|---|---|
| `Enter magnitude (V band)` | magnitudine V della sorgente | — |
| `Enter target bkg-subtracted amplitude (in %)` | ampiezza della modulazione che si vuole rivelare, **dopo sottrazione del fondo**, in % | — |
| `Assume default bkg rate of 15 kcps?` | rate di fondo; rispondendo `n` si inserisce il valore in kcps | 15 kcps |
| `Is this a Single Trial search?` | ricerca a singolo trial, oppure numero di trial indipendenti (es. frequenze o periodi indipendenti cercati) | single trial |

Le risposte affermative accettate sono invio, `y`, `yes`, `s`, `si`. Un valore non valido per
magnitudine o ampiezza interrompe il programma; per fondo e trial si torna al default.

### Cosa calcola

1. **Rate della sorgente** dalla magnitudine, con la retta di calibrazione:
   `R_net = 10^(-(m·mag + q)/2.512)` kcps, con `m = 0.942` e `q = -18.67`
   (valori fissati in `main()`; modificali se cambia la calibrazione).
   Il rate totale è `R_tot = R_net + R_bkg`.
2. **Soglia di potenza** `P_thr` per la probabilità di falso allarme `pfa = 0.0027`
   (Vaughan et al. 1994). Con `N_trials > 1` la soglia viene alzata di `2.24·ln(N_trials)`, in modo
   che la probabilità di falso allarme valga `pfa` sull'insieme dei trial.
3. **Sensibilità in potenza** `P`: potenza del segnale che supera la soglia con probabilità
   `1 − pfd`, con `pfd = 0.05` (probabilità di falsa esclusione al 5%, cioè rivelazione al 95%).
   È il parametro di non-centralità `chndtrinc(P_thr, 2·N_ps, pfd)`, con `N_ps = 1` spettro di
   potenza.
4. **Ampiezza diluita dal fondo**: `A = A_sub · R_net / R_tot`.
5. **Tempo di osservazione**: `T_obs = 1.61² · P / (R_tot · A²)`.
6. Conteggi integrati (`N_net`, `N_bkg`, `N_tot`) e un controllo di coerenza: l'ampiezza
   sottratta del fondo ricalcolata da `N_tot` deve tornare uguale a quella inserita.

### Esempio

V = 15, ampiezza 5%, fondo di default, single trial (`15`, `5`, invio, invio):

```
Net Source Rate (R_net)     : 64168.94 c/s
Background Rate (R_bkg)     : 15000.00 c/s
Total Count Rate (R_tot)    : 79168.94 c/s
Threshold Power (P_thr)     : 13.25
Search power sensitivity (P): 26.72
REQUIRED OBS. TIME (Tobs)   : 0.53 seconds
Amplitude A (including bkg) : 4.0527 %
(Double check) Bkg-subtracted amplitude: 5.0000 %
```

Lo stesso caso con 1000 trial indipendenti alza `P_thr` a 28.72 e il tempo richiesto a 0.96 s (invece di 0.53 s).
Il tempo è stampato in secondi, minuti o ore a seconda dell'ordine di grandezza.

### Parametri che si possono modificare nello script

Nel codice sono commentate le estensioni previste:

- `pfa`, `pfd`: probabilità di falso allarme e di falsa esclusione (ora `0.0027` e `0.05`)
- `N_ps`: numero di spettri di potenza sommati (ora 1). Per più spettri c'è già la riga
  commentata con `chdtri`
- `m`, `q`: retta di calibrazione
- `bkg_default`: rate di fondo di default (`1.5e4` c/s)

### Note

- Le costanti `1.12` (compare come `2.24 = 2 × 1.12` nella soglia) e `1.61` (nel tempo di
  osservazione) sono scritte nello script senza riferimento: citare qui la fonte aiuta chi usa
  lo strumento.
- Se `sp.special` non viene trovato con versioni vecchie di SciPy, sostituire
  `import scipy as sp` con `from scipy import special` e usare `special.chndtrinc`.

---

## 3. `light_curves_SiFAP2.py` — nodding, burst, k e curva netta

Analizza le curve di luce SiFAP2 con il metodo del **sky shift additivo via nodding**. È uno
script unico: i parametri si cambiano nel blocco **CONFIGURAZIONE** all'inizio del file, poi si
lancia

```
python light_curves_SiFAP2.py
```

Per una coppia di file `.qdp` (target e sky, 4 colonne: tempo, errore sul tempo, rate, errore sul
rate) lo script:

1. trova gli episodi di **nodding**, in modo automatico o dagli intervalli indicati
2. esclude i **burst** dai calcoli (restano visibili nei plot)
3. calcola il **fattore additivo k**
4. costruisce il cielo traslato e la **curva netta**
5. calcola i **rate medi** fuori nodding (target, sky + k, netto)
6. salva il plot in PDF e, a richiesta, la curva netta e un riepilogo in CSV

### Metodo

```
k            = <target | nodding> − <sky | riferimento>
sky traslato = sky + k
curva netta  = target − sky traslato
```

I rate medi sono medie pesate (peso 1/σ²), calcolate fuori dai nodding e fuori dai burst. Se
target e sky hanno griglie temporali diverse, lo sky traslato viene interpolato sui tempi del target.

Il riferimento per lo sky si sceglie con `K_MODE`:

| `K_MODE` | Riferimento per lo sky | Quando usarlo |
|---|---|---|
| `"window"` (default) | bin dello sky in una finestra prima (`"pre"`) o dopo (`"post"`) ogni episodio di nodding | curva senza trend marcati |
| `"simultaneous"` | bin dello sky durante il nodding stesso | curva con trend: target e sky sono simultanei, quindi il trend si cancella |

In modalità `"window"` la finestra si sceglie episodio per episodio (`K_WINDOWS`). Se le finestre
contengono meno di 3 bin dello sky si usa l'intera baseline.

### Parametri

Tutti nel blocco CONFIGURAZIONE, con commenti e valori di esempio. I tempi sono in secondi.

| Parametro | Significato | Default |
|---|---|---|
| `OBS_ID` | nome dell'osservazione (titolo del plot, nomi dei file di output) | dal nome del file (es. `N3`) |
| `TARGET_FILE` | file `.qdp` del rivelatore target | — |
| `SKY_FILE` | file `.qdp` dello sky | `<TARGET_FILE>_sky.qdp` |
| `OUTPUT_DIR` | cartella dove salvare i risultati | cartella del file target |
| `NODDING` | `None` = rilevamento automatico, oppure `[(t1, t2), ...]` | `None` |
| `NODDING_THRESHOLD` | auto: rate < soglia × mediana → candidato nodding | `0.95` |
| `NODDING_MIN_BINS` | auto: bin contigui minimi per un episodio valido | `5` |
| `BURST` | `[(t1, t2), ...]` esclusi dai calcoli | `[]` |
| `K_MODE` | `"window"` o `"simultaneous"` | `"window"` |
| `K_WINDOWS` | solo `"window"`: finestra per episodio, nell'ordine, es. `[("post", 250), ("pre", 250)]` | `None` |
| `K_WINDOW_DEFAULT` | finestra per gli episodi senza voce in `K_WINDOWS`; `None` = intera baseline | `("pre", 250)` |
| `START_UTC` | opzionale: inizio osservazione, scritto nella colonna `date` del riepilogo | `None` |
| `SAVE_PLOT` / `SHOW_PLOT` | salva il PDF / mostra il plot a schermo | `True` / `False` |
| `SAVE_NET_CSV` | salva la curva netta in CSV | `False` |
| `SUMMARY_CSV` | file CSV con una riga per osservazione | `None` |
| `PLOT_YLIM` | limiti dell'asse y del pannello superiore, es. `(23000, 37000)` | automatico |
| `PLOT_LABEL_DY`, `PLOT_LABEL_DY_NET` | altezza delle etichette nel plot (in sigma) | `5.0`, `3.5` |

Esempi:

```python
# nodding manuale (due episodi) e un burst, con il cielo preso durante il nodding
NODDING = [(1063, 1093), (6420, 6456)]
BURST   = [(4880, 5000)]
K_MODE  = "simultaneous"

# nodding automatico, con il cielo preso 250 s dopo il primo episodio e 250 s prima del secondo
NODDING   = None
K_MODE    = "window"
K_WINDOWS = [("post", 250), ("pre", 250)]
```

### Più osservazioni

Per analizzarne più di una in un'esecuzione, elenca in `OBSERVATIONS` solo i parametri che
cambiano rispetto al blocco principale:

```python
OBSERVATIONS = [
    dict(TARGET_FILE="oss_A.qdp", NODDING=[(1063, 1093), (6420, 6456)]),
    dict(TARGET_FILE="oss_B.qdp", NODDING=None, BURST=[(4880, 5000)], K_MODE="simultaneous"),
]
```

Con il valore di default `[{}]` viene analizzata una sola osservazione.

### Output

Nella cartella del file target, oppure in `OUTPUT_DIR`:

- `LC_SiFAP2_<id>_scaled.pdf`: curve target, sky e sky traslato in alto, curva netta in basso;
  nodding in grigio, finestre per k in giallo, burst in arancione
- `LC_SiFAP2_<id>_net.csv` (con `SAVE_NET_CSV`): tempo, rate del target, sky traslato, netto e flag
  `in_nodding` / `in_burst` per ogni bin
- il file `SUMMARY_CSV`: una riga per osservazione con k e rate medi

### Dal rate alla magnitudine

Al termine di ogni osservazione lo script stampa il rate netto con la durata, pronti per
`magnitude_SiFAP2.py`:

```
--rate 5953.52 --rate-err 2.85 --duration 5406.7
```

Le prime colonne del riepilogo CSV (`id, date, duration, rate, rate_err`) sono quelle richieste
da `magnitude_SiFAP2.py --obs-file`: imposta `START_UTC` (o compila la colonna `date`) e passa il
file direttamente a quello script. La durata è quella tra il primo e l'ultimo bin della curva.

### Note

- Se un burst si sovrappone a un nodding, i suoi bin sono esclusi anche dal calcolo di k.
- Il plot usa un backend non interattivo: senza `SHOW_PLOT` non serve un'interfaccia grafica.
- Se un'osservazione fallisce (file mancante, nessun nodding trovato), le altre vengono
  analizzate comunque e l'errore compare nel riepilogo finale.
