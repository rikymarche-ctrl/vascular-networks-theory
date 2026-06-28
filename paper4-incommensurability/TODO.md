# TODO — Lista azioni CONSOLIDATA (file UNICO)

> **Questo è l'unico file TODO.** Unisce il referee report (Blocchi A–F, tutti risolti) e
> l'audit interno "Opus Pre-Giacomin" (Blocco G, 6 residui aperti — vedi in fondo).
> Supera ed elimina logicamente `TODO.md` e `TODO 2.md` in `Progetti IA\branching papers\`.

---

# ⭐ 0. PEER REVIEW UFFICIALE — *Transport Phenomena* (TP-1083) — ricevuta 2026-06-21 — **PRIORITÀ ASSOLUTA**

> Questa è la **vera** peer review (decisione: REVISIONE). Sostituisce in priorità l'audit
> interno (Blocchi A–G sotto, che restano come lavoro autoriale già fatto / **G1–G6 ancora aperti**).
> **Stato verifica highlight 2026-06-21:** 81 span colorati nel sorgente (tutti `\chA` magenta = autore),
> 4 occorrenze "note that"/"notable"…, 2 "etc."/"and so on".
>
> **Formato richiesto dall'editor:** copia evidenziata con **un colore diverso per ciascun reviewer**
> + risposte punto-per-punto separate ("Response to Reviews"), **una per referee**.
> **Mappa colori** (toggle `\ifhighlight` già nei due preamboli): 🔵 `\chR`=Reviewer 1 · 🔴 `\chRB`=Reviewer 2 ·
> 🟢 `\chE`=Editoriale · 🟣 `\chA`=autore (audit interno A–G, già applicato). Le nuove modifiche da R1/R2/editor
> vanno avvolte nel COLORE del rispettivo richiedente.

## 0.E — Richieste EDITORIALI (🟢 `\chE`)
- [x] **E1 — Copia evidenziata multi-colore + copia pulita** dal sorgente unico (toggle `\ifhighlight`). R1=blu, R2=rosso, editor=verde, autore=magenta.
  - **FATTO 2026-06-21:** Toggle reso controllabile da CLI in entrambi i preamboli: `\newif\ifhighlight \ifdefined\CLEANCOPY\highlightfalse\else\highlighttrue\fi`. **Default invariato** = copia evidenziata (4 colori), prodotta da `build.ps1` (ZERO regressione). Nuovo `build-clean.ps1` produce la **copia pulita** passando `\def\CLEANCOPY{}` a pdflatex (output con suffisso "(Clean Copy)"). Testato: clean main compila exit 0, nessun errore, toggle risolto correttamente. `build.ps1` lasciato intatto.
- [x] **E2 — "Response to Reviews" separate**, una per R1 e una per R2, punto-per-punto (documento a parte, non nel manoscritto).
  - **FATTO 2026-06-21:** Create due lettere LaTeX standalone in nuova cartella `response/`: `Response-to-Reviewer-1.tex` (1 pag.) e `Response-to-Reviewer-2.tex` (3 pag.). Ciascuna: intestazione TP-1083/journal, intro di ringraziamento, **point-by-point** con il commento del referee (parafrasato fedelmente dal report, NON quote verbatim inventate) seguito da "Response." + posizione/colore esatti della modifica. R1: solo R1.1 (riformulazione gittata, Wo∝M^{1/4}, additiva, blu) + nota di consistenza §5–6. R2: R2.1 (sottosez. assunzioni §7), R2.2 (test discriminante cross-regno §6+S6), **R2.3 segnalato come punto centrale** (sottosez. S7 + tabella 6 funzionali × 4 assiomi), R2.4 (passerelle §2+§5); chiusura che segnala anche i refinements autore-iniziati (magenta: citazione Fåhræus–Lindqvist corretta, disambiguazione η_d/η e β/β_M). Colori coerenti col manoscritto (R1→blu, R2→rosso). **Compilano puliti** (pdflatex ×2): R1 1 pag, R2 3 pag, **0 errori `!`**, 0 undefined. Ausiliari ripuliti; `response/` = solo 2 `.tex` + 2 `.pdf`. → **E2 COMPLETO. Pacchetto di revisione TP-1083 completo.**
- [x] **E3 — Completare TUTTE le liste che finiscono in "etc."/"and so on"/"and the like"** dove la continuazione non è ovvia (2 occorrenze: supplemental, Section-03). 🟢
  - **FATTO 2026-06-21 (🟢 `\chE`):** Section-03:70 `(fixed total volume, hierarchical structure, etc.)` → `(fixed total volume and hierarchical connectivity)` (lista chiusa). supplemental:784 `H_1..H_3, etc.` → `..., with successive harmonics (n≥4) carrying progressively smaller weights`. Re-grep `etc.|and so on` su tutti i .tex = **0**.
- [x] **E4 — Eliminare ovunque** "notable/notably/noteworthy/note that/(it should be) noted that/it is noteworthy" e simili (4 occorrenze: supplemental, Section-03×2, Section-05). 🟢
  - **FATTO 2026-06-21 (🟢 `\chE`):** Section-05:55 `Note that the` → `The`; Section-03:46 `It is crucial to note that these` → `These`; Section-03:176 `though we note that alternative` → `though alternative`; supplemental:352 `Note that the combined shift` → `The combined shift`. Re-grep parole bandite su tutti i .tex = **0**.
- [x] **E5 — Eliminare parole vuote** ("process", "system", e simili) dove possibile senza perdere senso. 🟢
  - **FATTO 2026-06-21 (🟢 `\chE`):** Revisione di tutte le ~42 occorrenze. "process": unica occorrenza in tutto il progetto era Section-07:218 `angiogenic and remodeling processes` → `angiogenesis and remodeling` (riempitivo eliminato); re-grep `process` = **0**. "system": le restanti ~37 sono usi TECNICI legittimi (physical/biological/composite/sub/vascular/arterial/physiological/organ systems, "composite system" nelle dim. di separabilità) — tagliarle perderebbe senso, quindi conservate per scelta motivata. (Section-06:384 "systems and regimes" e Section-03:11 "any system" valutati e tenuti: non ridondanti.)
- [x] **E6 — Figure ad alta risoluzione:** 300 dpi (foto), 600 dpi (linework/grafici) alle dimensioni finali, testo nitido. fig_phase_diagram + fig_transition sono TikZ→vettoriali (OK, fornire PDF/EPS vettoriale); **womersley_verification.png da rigenerare a ≥600 dpi** (`generate_womersley_figure.py`, parametro dpi).
  - **FATTO 2026-06-21:** `generate_womersley_figure.py` ora salva **PDF vettoriale** (`womersley_verification.pdf`, qualità press, risoluzione infinita) + PNG a **600 dpi** (fallback). Riferimento in `supplemental.tex:750` aggiornato a `.pdf`. Ricompilato: supplemento exit 0, figura inclusa come `Graphic file (type pdf)`, nessun warning file mancante. fig_phase_diagram + fig_transition restano TikZ vettoriali (già OK).

## 0.R1 — Reviewer #1 (🔵 `\chR`) — *semplice*
- [x] **R1.1 — Apertura topo-vs-balena (600 vs 20 bpm):** confrontare il **VOLUME di pompaggio** (gittata/portata) invece della sola frequenza cardiaca. Riformulare l'esempio iniziale (Section-01) in termini di portata; verificare che il punto sul numero di Womersley regga (balena: bassa frequenza ma gittata enorme; topo: alta frequenza, gittata minima).
  - **FATTO 2026-06-21 (🔵 `\chR`, Section-01:6-20):** Riformulato l'incipit. La frequenza cardiaca da sola è FUORVIANTE: poiché $\mathrm{Wo}=R\sqrt{\omega/\nu}$ *cresce* con la frequenza, i 600 bpm del topo lo renderebbero più inerziale della balena — l'opposto della tesi. Correzione rigorosa: gittata cardiaca $\dot Q\propto M^{3/4}$, gittata sistolica $\propto M$ (una contrazione di balena espelle ~$10^7\times$ il volume del topo); il raggio aortico $R\propto M^{3/8}$ domina su $f\propto M^{-1/4}$, quindi $\mathrm{Wo}\propto M^{3/8}\sqrt{M^{-1/4}}=M^{1/4}$ → **cresce** con la massa nonostante il ritmo più lento. Risolve il paradosso del reviewer (topo veloce ma viscoso).
  - **Coerenza interna verificata:** identica derivazione già provata in §6 (Section-06:120-124); aggiunto cross-ref `(\S\ref{sec:womersley_minimax})`. **Check numerico pulito** (mouse→human, entrambi seguono l'allometria standard): $M^{1/4}$ predice $\mathrm{Wo}_\text{umano}/\mathrm{Wo}_\text{topo}=(70000/20)^{1/4}\approx 7.69$ vs tabella §6 ($17.9/2.5\approx 7.16$), entro ~7%. (NB: il check umano/colibrì NON è valido — il colibrì viola di proposito $f\propto M^{-1/4}$: è l'override di frequenza, Wo=1.76 da 1000 bpm espliciti.)
  - **Revisione resa ADDITIVA (richiesta utente):** il testo originale è preservato verbatim (nero, intatto); l'unica modifica è l'inserimento blu → tutte le modifiche evidenziate, nessuna cancellazione non tracciabile.
  - **Chiarimento di coerenza con §6 (🟣 `\chA`, Section-01:21-25):** la frase originale "the mouse's is governed by the viscous linearity of Murray's law" è leggibile come "architettura topo = Murray α=3.0", ma §6 (tab:hummingbird) predice topo α*=2.60 (transition), "NOT α≈3.0". Aggiunta una frase 🟣 che inquadra i due regimi come i **limiti di singolo-meccanismo IN ISOLAMENTO** della Lagrangiana unificata — attrattore d'onda $\alpha_w=\VarAlphaW$ (balena inerziale) e ottimo viscoso $\alpha_t\approx\VarAlphaT$ (analogo network di Murray, topo) — così il "Yet, α* nearly identical ($\approx\VarAlphaStar$)" successivo legge come il paradosso, e Murray NON viene letto come l'esponente realizzato del topo. Coerente con §5 (corollario limiti) e §6 (valori realizzati). Compila pulito (EXIT 0, 0 errori).

## 0.R2 — Reviewer #2 (🔴 `\chRB`) — *molte modifiche di testo, da rendere perfetto*
- [x] **R2.1 — Discussione estesa delle ASSUNZIONI e delle ALTERNATIVE.** Le derivazioni poggiano su assunzioni (ottimizzazione evolutiva, vincoli informazionali, fitness landscape metabolici) non univocamente garantite: aggiungere discussione esplicita di queste assunzioni e delle possibili alternative (anche se il test empirico non è sempre fattibile). (Discussion §7.)
  - **FATTO 2026-06-21 (🔴 `\chRB`, Section-07 sottosezione finale `sec:assumptions` "Underlying Assumptions and Alternative Hypotheses"):** Nuova sottosezione di chiusura della Discussione che enumera le **5 assunzioni portanti** con, per ciascuna, l'alternativa principale e l'osservazione discriminante: **(A1) Ottimalità** (alt: canalizzazione genetica / auto-organizzazione fisica / evoluzione neutrale; (ii) NON è in conflitto — i modelli dinamici Hu-Cai/Ronellenfitsch convergono allo stesso ottimo, "meccanismo vs target"; (i) blueprint genetico è il bersaglio di falsificazione del test ontogenetico §8); **(A2) Minimax vs average-case** (alt: Bayesian expected-cost / single-objective / satisficing; firma robustezza esclude le ultime due; **non-unicità onesta**: un prior largo simmetrico converge a ~stesso α*, quindi minimax adottato come scelta conservativa assumption-light); **(A3) Bound informazionale** (alt: coordinamento long-range / pre-specifica genetica; vincolato dal phase-lag blind spot §2, argomento sulla *disponibilità informativa* robusto alla modalità di segnalazione); **(A4) Currency metabolica** (alt: altre currency selettive — robustezza sviluppo, fatica, termoregolazione, riparo, riproduzione; mitigato da `thm:rigidity` (α* indipendente dalle costanti metaboliche) + un *terzo* canale incommensurabile esce dal minimax a due canali, cfr. rene); **(A5) Incommensurabilità stessa** (alt: budget energetico unico → basterebbe Murray/WBE; è la distinguibilità-sotto-scaling di `rem:incommensurability`, ciò che `thm:gauge` prova non rimovibile; firma empirica = struttura dual-attractor). Paragrafo finale **"Epistemic status"**: quali assunzioni sono falsificabili (A1,A3), quali supportate (A5), quali delimitano il dominio (A4), quale meno separabile (A2). Aggiunta `\label{rem:incommensurability}` (invisibile, no-wrap) al remark di distinguibilità per citarlo. Titolo sottosezione lasciato in nero (convenzione del manoscritto: 0 titoli colorati, 0 `\texorpdfstring` → evita warning hyperref bookmark); corpo **interamente `\chRB`**. `build.ps1` exit 0, entrambi i PDF prodotti; check 2-passate: tutte le label R2.1 (`prop:minimax_saddle`, `rem:gauge_parsimony`, `rem:incommensurability`, `thm:nogo`, `tab:ontogenetic`, `thm:rigidity`, `thm:gauge`) **risolte**.
- [x] **R2.2 — Validazione empirica più ampia.** Discutere/rafforzare i test del Principio di Incommensurabilità e dell'attrattore minimax su una diversità più ampia di sistemi biologici (espandere oltre le 27 reti; collega a C6).
  - **FATTO 2026-06-21 (🔴 `\chRB`, Section-06 §"Independent Datasets" + intro S6 supplemento):** Invece di gonfiare la tabella con dati deboli/inventati (scientificamente improprio), **riformulato le 27 reti come TEST DISCRIMINANTE cross-regno** del principio. Aggiunto paragrafo `\chRB` in §6: il principio NON predice un esponente universale ma che l'esponente realizzato segua *quali canali di costo incommensurabili sono fisicamente attivi* → regola di smistamento *parameter-free*: (i) reti **senza canale d'onda pulsatile** (xilema liana/conifera/angiosperma, vie bronchiali/tracheali, venatura fogliare, slime mold + miceli, ritorno venoso quasi-statico → **10 sistemi**, tutti in 2.80–3.15) all'attrattore viscoso α_t→3.0; (ii) arterie mammifere pulsatili (entrambi i canali) al minimax α*≈\VarAlphaStar; (iii) reti wave-dominated/planari (retina d=2, carotide/condotti) verso α_w≈\VarAlphaW. Lo smistamento regge su **3 regni** (animale/vegetale/fungino-protista), **entrambe le dimensioni** d, e stati **patologici** (ipertensione/tumore → limite wave-skewed) — span irriproducibile da una teoria a singolo attrattore (Murray α=3 o area-preserving α=2). **Falsificabilità bidirezionale**: una rete non-pulsatile a α≈\VarAlphaStar, o un condotto pulsatile ad alto Wo a α≈3.0, contraddirebbero la regola. Chiusura: direzioni di estensione *oltre* le 27 (duty cycle divergenti — circolazioni aviarie/rettiliane; sistemi circolatori aperti invertebrati; organ-on-chip a pulsatilità programmabile). Cornice `\chRB` gemella nell'intro S6 (catalogo assemblato *come test discriminante*, non lista di accordi). `build.ps1` exit 0, entrambi i PDF. Solo macro esistenti, nessuna nuova label.
- [x] **R2.3 — Alternative al "linear functional excess" (PUNTO CENTRALE di R2).** Il claim che l'eccesso funzionale lineare sia l'**unico** funzionale di costo ammissibile va supportato meglio: considerare ESPLICITAMENTE formulazioni di costo alternative, valutarne i limiti, e giustificare l'esclusione a favore del lineare. (Section-03 Gauge + supplement S7 Jensen/Onsager; collega a G1/G2.)
  - **FATTO 2026-06-21 (🔴 `\chRB`):** Il materiale era sparso (esclusione log §3:498, Jensen non-lineare S7:1235, non-sostituibilità §3:582, terzo costo §3:599). Consolidato in nuova sottosezione S7 `sec:S7_alternatives` "Systematic Exclusion of Alternative Cost Functionals": 4 assiomi di ammissibilità (scale-invariance, consistenza compositiva/Jensen [G2], Onsager linear-response [G1], non-sostituibilità) + **tabella `tab:cost_alternatives`** che testa lineare/quadratico/log/power-law/absolute-threshold/moltiplicativo → solo il lineare passa tutti e 4. Frase conclusiva: linearità *forced* dai 4 requisiti, rilassabile solo rigettandone uno (target di falsificazione). Puntatore in §3 (remark non-sostituibilità) → "Supplemental Material (Section S7)". Supplement compila exit 0, `\ref`/`\eqref` risolvono.
- [x] **R2.4 — Densità espositiva.** Espandere i passaggi intermedi e aggiungere interpretazione biologica più profonda dei risultati formali (pass di accessibilità Sezioni 2–6, costrutti introdotti troppo rapidamente).
  - **FATTO 2026-06-21 — prima versione (2 passerelle §2+§5), poi AMPLIATO 2026-06-21 su richiesta esplicita utente** ("testo meno denso con più step concettuali su TUTTE le §2–6, parti dall'inizio del file, 1 sezione alla volta, fatto bene; abbiamo più testo, non rinunciare a nulla"). La prima versione era sotto-dimensionata (avevo razionalizzato col "non gonfiare" — scusa, non vera soluzione). **Pass completo (🔴 `\chRB`), una sezione alla volta dall'inizio, +10 passerelle concettuali (totale `\chRB` 21→31):**
    - **§1 Introduction (×2):** (a) lettura *in parole* del Kinematic Matching Theorem tra enunciato e dimostrazione (cos'è Wo_c come "interruttore"; sotto soglia = resistori viscosi/Murray, sopra = guide d'onda pulsatili; due soglie = incommensurabilità); (b) glossa intuitiva del *modo evanescente* nel Livello 3.
    - **§2 Scaling-Conflict (×1 nuova, +1 preesistente):** analogia del muratore per *cosa afferma* la Proposizione no-go (serve sapere $g$ = "altezza nel muro" → serve $r_0$, info non-locale) che prepara lo Step 3.
    - **§3 Gauge/ATP (×3):** (a) riconciliazione "costo ∝ eccesso *assoluto*" vs "eccesso *frazionario*" (within- vs across-organism); (b) lettura in parole del Teorema di Onsager (motore del funzionale lineare: penalità lineare nell'eccesso frazionario, quadratica nella deviazione geometrica → tollera eterogeneità); (c) glossa di *cosa chiede* l'equazione di Jensen pesata (score combinato = media degli score, indipendente dal confine).
    - **§4 Architectural-Invariance (×2):** (a) **gioco contro la Natura** — la rete sceglie α, l'avversario (stato fisiologico) sceglie η peggiore, il saddle è la geometria *immune* al regime (cuore concettuale del minimax); (b) glossa del teorema dell'inviluppo agganciata all'intuizione del gioco.
    - **§5 Single-Mechanism-Limits (×1, fatto nella prima versione):** interpretazione biologica dei due limiti (statico→Murray vasi piccoli/non-pulsatili; wave→impedance matching grandi condotti; minimax→letti medio-calibro campionati).
    - **§6 Architectural-Transition (×1 nuova, +1 preesistente R2.2):** lettura in parole del **Phase Decoupling** (fluido e onda "si accendono" a taglie diverse perché $Y_c\propto\sqrt{Y_L}$ → disallineamento permanente = incommensurabilità).
    - **§9 Retinal-Paradox (×1):** enunciato esplicito del *paradosso* (diametri 2D α≈2.0 vs angoli lontani da 90° planare) e della risoluzione (canali *decoupled*: diametri→fluido/onda dimension-sensitive, angoli→equilibrio meccanico) in apertura.
    - **§7/§8 NON modificate**: sono sintesi/discussione (già prosa interpretativa ricca; §7 ha R2.1 + molti remark). Non sono "costrutti introdotti troppo in fretta" → aggiungervi testo sarebbe padding. (Da confermare con utente se vuole comunque un pass di fluidità.)
    - Nessuna nuova `\ref`/`\label` introdotta (solo prosa + macro/simboli esistenti). `build.ps1` exit 0 ai checkpoint (§1–3 e §4/6/9), entrambi i PDF. Lettera R2 (Comment 4) aggiornata di conseguenza. → **Blocco R2 (R2.1–R2.4) COMPLETO.**

## 0.CV — VERIFICA del *VALORE* CITATO (🟣 `\chA`) — **richiesta esplicita utente 2026-06-21**
> NON i metadati (volume/anno) ma **il valore/affermazione attribuito alla fonte**: ogni numero o claim
> con `\cite` va verificato contro il **CONTENUTO** del PDF in `branching papers\…\ARTICOLI CITATI`.
> **Le 2 citazioni AGGIUNTE da me sono ora verificate leggendo i PDF per intero (nessun dato contraddetto;
> esiti sotto). Tutte le ALTRE citazioni pre-esistenti: già controllate dall'utente.**
- [x] **CV1 — c_wave = 5–8 m/s → `huokassab2006`: CONFERMATO PER DERIVAZIONE.** Il PDF non stampa "5–8 m/s" esplicitamente, MA l'appendice (H1085–86) dà `c₀=√(Eh/ρR)` con il valore misurato **E·h/R = 4.0×10⁵ dyn/cm²**; con ρ=1.06 g/cm³ → **c₀ ≈ 6.1 m/s** (e `c=√(1−F10)·c₀ ≤ c₀`), che cade dentro 5–8 m/s. Citazione VALIDA. **Fix applicato:** tabella Sec-07:715 tipo "Measured"→`\chA{Literature}` (la velocità d'onda è modellata da E,h,R misurati, non misurata direttamente).
- [x] **CV2 — Fåhræus-Lindqvist "~30% a r<10 µm" → `fahraeuslindqvist1931`: ERRORE DI SCALA, CORRETTO.** Letto il PDF intero (p.565–567): dati solo IN VITRO in tubi di vetro fino a **0.04 mm (40 µm)** minimo; viscosità relativa ~4.6 a 0.3 mm → ~2.8 a 40 µm (~39%), e gli autori dichiarano esplicitamente **"~50% a 0.03 mm"**. Quindi: (a) **r<10 µm è FUORI dai loro dati**; (b) il loro numero è ~40–50%, **non 30%**; (c) è in vitro (in vivo è attenuato). Il "~30% a r<10 µm" attribuito a FL1931 **non era sostenibile**. **Fix applicato (Sec-09:86-93):** riscritto citando il loro reale risultato (deviazione da Poiseuille sotto ~0.3 mm; ~½ a ~30 µm in vitro), retina = "still narrower microvessels", e α_eff≈2.84 ribattezzato "conservative in-vivo estimate" (in vivo < in vitro). Tutto `\chA`.
  - ⚠️ **Nota rigore (separata):** `\VarAlphaFahraeusEff = 2.84` è **hardcoded** in compute.py:861 (stima d'autore, non derivata). Valutare se derivarla o dichiararla stima (simile a C3). NON bloccante per CV2.
- [x] **CV3 / CV4 — tutte le altre citazioni pre-esistenti**: l'utente conferma di averle **già verificate** personalmente (`kassab1993`, `jiang1994`, `taylor2024`, `Guo2003`, ecc.). Fuori dal mio scope: ho verificato SOLO le 2 che avevo inserito io.

## 0.AUDIT — Confermare che il lavoro A–G dichiarato è REALMENTE nel sorgente
- [x] **AUDIT1** — L'utente percepisce poche modifiche evidenziate. Dato: 81 span `\chA`. Incrociare ogni fix dichiarato (A–F, C, D) con la sua reale presenza nel `.tex` e col colore giusto; nessun fix "solo dichiarato nel TODO" ma assente nel testo.
  - **FATTO 2026-06-21 — AUDIT SUPERATO** (conteggio `\chRB` aggiornato dopo l'ampliamento R2.4: ora **31**, era 21). Censimento span su tutti i `.tex`: **`\chA`=112** (autore A–G; cresciuto da 81 perché include ora G3–G6), **`\chR`=1** (Reviewer 1 = R1.1, unica inserzione additiva §1, reale e sostanziale), **`\chRB`=31** (Reviewer 2: 21 R2.1–R2.3 + 10 passerelle R2.4), **`\chE`=7** (editoriale). **Riconciliazione esatta**: `\chRB`=21 = R2.3(4: §3:596 + S7 supp 1254/1282/1301) + R2.1(13: §7 framing + 5 paragrafi ×(titolo+corpo) + epistemic status) + R2.2(2: §6:409 + S6 supp:956) + R2.4(2: §2:54 + §5:41); `\chE`=7 = E3(2)+E4(4)+E5(1). **Macro**: tutti e 4 definiti IDENTICI in entrambi i preamboli (main:36-39, supp:33-36) → blu/rosso/verde/magenta. **Toggle** default = copia evidenziata (main:35, supp:32: `\ifdefined\CLEANCOPY...`). **Check "TUTTE"** via git diff: le righe aggiunte apparentemente non-avvolte sono **artefatti di reflow** (es. §4: l'unica modifica è `duty cycle`→`\chA{marginal-balance weight}`, la prosa adiacente "The numerator derives…" è preesistente, ricomparsa nel diff solo per ri-spezzettamento riga). **Nessun fix dichiarato-ma-assente; nessuna prosa nuova non evidenziata.**
  - ⚠️ **Nota cosmetica (opzionale, non bloccante):** `\chE` usa `\color{green}` puro (basso contrasto su bianco). Il piano suggeriva `green!55!black` per leggibilità. Lasciato com'è (scelta sessioni precedenti, identico nei due preamboli, l'utente non l'ha segnalato). Cambiarlo è 1 riga ×2 preamboli se desiderato.

---


> Manoscritto: *The Incommensurability Principle in Biological Transport*
> (arXiv:2605.03219), sottomesso a **Transport Phenomena** (De Gruyter).
> Review indipendente eseguita su `manuscript/` + `supplements/` + `scripts/` + `figures/`
> (lettura del README fatta **dopo** l'analisi, come richiesto: conferma B1 e l'angolo retinico stantio).
>
> **Verdetto referee: MAJOR REVISION.**
> Lavoro ambizioso, matematicamente ricco, in gran parte algebricamente coerente
> (verificati a mano: Q⁻¹=6/Wo², radice 1.740, sensibilità App. A, Wo₀∝M¹ᐟ⁴, colibrì 1.76, 1490 m⁻¹, r_crit 1.16 mm…).
> MA: un probabile **errore matematico nel cuore della tesi** (A1) + diverse
> **contraddizioni numeriche interne** che un referee trova subito.
>
> Legenda severità: 🔴 CRITICO · 🟠 MAJOR · 🟡 MINOR
> Stato: [ ] da fare · [~] in corso · [x] fatto

---

## A. PROBLEMI CRITICI (sostanza scientifica)

### ✅ A1 — RITIRATO: il PAPER HA RAGIONE (verifica numerica Bessel — esito riassunto qui sotto; nota di lavoro `A1_VERIFICATION.md` rimossa il 2026-06-21)
> **ESITO 2026-06-10:** la critica A1 era SBAGLIATA. L'admittanza caratteristica corretta è
> `Y_c ∝ √(1−F10)` (Womersley/McDonald), la cui radice Bessel esatta è **2.144** (≈ 3/√2 = 2.121,
> stesso ~1% di accordo del threshold fluido). Il mio errore: lettura letterale di "Y_c ∝ √Y_L" (→ 2.927),
> trascurando il fattore `√(iωρ)=e^{iπ/4}` (il π/4 mancante). Soglia 2D→0 anch'essa CONFERMATA.
> Resta solo un suggerimento minore/costruttivo: scrivere `Y_c ∝ √(1−F10)` (non `√Y_L`) e aggiungere
> la verifica numerica della soglia d'onda (NOI l'abbiamo fatta: conferma). **Testo originale superato sotto.**

**[SUPERATO] È il punto più importante: ci poggia tutto il "doppio-soglia" e il Paradosso Retinico.**
File: `supplements/supplemental.tex:77-102` (S1); `manuscript/sections/Section-06-Architectural-Transition.tex:184-238`; `Section-09-Retinal-Paradox.tex:27-40`.

Dettaglio matematico:
- A basso Wo ho ricalcolato esplicitamente: `Y_L ∝ Wo²/8 − i·Wo⁴/48`.
- Argomento (fase standard) di Y_L: `arg(Y_L) = arctan(−Wo²/6)`, in modulo **arctan(Wo²/6)** (piccolo → Y_L quasi reale/resistivo).
- `Y_c = √Y_L` dimezza QUESTO argomento ⇒ `Q⁻¹_{Y_c} = cot(½·arctan(Wo²/6))`.
- Nel paper invece: `Q⁻¹_{Y_c} = cot(½·arctan(6/Wo²))`. Ma `arctan(6/Wo²) = π/2 − arctan(Wo²/6)`:
  è stato dimezzato l'**angolo di perdita** (il complemento), NON l'argomento complesso che la √ dimezza.

Conseguenze (non cosmetiche):
- Scelta del paper: Q⁻¹_{Y_c} **cresce** con Wo; `Wo_c^wave(d=2)=0` → "onda sempre sopra soglia in 2D" (Paradosso Retinico).
- Dimezzamento standard: Q⁻¹_{Y_c} **decresce** con Wo (coerente con Y_L); `Wo_c^wave(3D)≈2.83` (non 2.121);
  **`Wo_c^wave(d=2)→∞`** → la conclusione retinica si **RIBALTA**.

Aggravante: la soglia d'onda 2.121 **non è MAI verificata numericamente** contro Bessel.
In `scripts/main/compute.py:22` è hardcoded `Wo_c = math.sqrt(3)`; `scripts/verification/generate_womersley_figure.py`
verifica solo la soglia **fluida** (Q⁻¹=2 → 1.740). La soglia d'onda regge su quel singolo passaggio di mezzo-angolo.

- [x] **FATTO 2026-06-22:** Implementato `scripts/verification/verify_wave_threshold.py` (Bessel esatto via mpmath). Esito definitivo: **lettura A** (paper, `Y_c ∝ √(1−F10)`, Womersley-McDonald) → **Wo_c^wave(3D)=2.144** (vs 3/√2=2.121, ~1%, stesso accordo della soglia fluida 1.740 vs √3); **d=2 → nessuna radice finita, Q⁻¹_{Y_c}≥1 ∀Wo>0 e CRESCE con Wo → Wo_c^wave(d=2)=0** (NON ∞ → il Paradosso Retinico **regge**, non si ribalta). **Lettura B** ingenua (`√Y_L`, omette il fattore `e^{iπ/4}` di `√(iωρ)`) → 2.927: **è l'errore, esclusa**. Conferma piena dell'esito A1 del 2026-06-10.
- [x] **FATTO 2026-06-22:** Non serve correggere nulla a valle (2.121 confermato, collasso 2D confermato). Aggiunta nota di verifica numerica `\chA` nel supplemento (sottosezione "Numerical Verification of the Critical Threshold", accanto a quella fluida): Wo_c^wave≈2.144 (~1%), collasso 2D→0, lettura √Y_L (2.93) esclusa. Il disclaimer `\chA` di §6 (estensione del criterio a Y_c = *postulato*) resta corretto: verificata la **conseguenza numerica** del postulato, non il postulato stesso. Riga "2.121 mai verificata numericamente" (sopra) ora **superata**.

### 🟡 A2 — [DECLASSATO a chiarezza dopo verifica A1] La classificazione dei regimi in 2D è confondente ma NON errata
> Dopo la verifica A1 (esito riassunto nel blocco A1 sopra): il risultato d=2 (soglia d'onda→0) è numericamente CONFERMATO. Resta solo un
> punto di chiarezza: Re/|Im| ha significato opposto per Y_L (alto=overdamped) e Y_c (alto=più propagante).
> Riusare le stesse etichette è confondente, ma non è un errore. Testo originale sotto.
File: `Section-06-Architectural-Transition.tex:232-238`; def. regimi in `Section-01-Introduction.tex:145-158` e `supplemental.tex` S4.
- Per costruzione: `Q⁻¹ > d−1` = **sovrasmorzato → NESSUNA propagazione d'onda**.
- In 2D scrivi "Q⁻¹_{Y_c} ≥ 1 per ogni Wo>0" e lo chiami **wave-dominated**. Ma `Q⁻¹ ≥ d−1` è esattamente *sovrasmorzato* (niente onda) secondo la tua stessa definizione. Logica regime↔soglia invertita (indipendente da A1).
- [x] **FATTO 2026-06-20:** Chiarita la convenzione di segno in Section-06:232. Il criterio $\mathcal{Q}^{-1}_{Y_c}\ge d-1$ è applicato all'**ammettenza caratteristica (onda trasmessa) $Y_c$**, per cui un valore alto = canale d'onda permanentemente impedance-matched/propagante — **senso opposto** alla convenzione di smorzamento del modo longitudinale (dove alto = overdamped). Aggiunta nota esplicita che i due sensi non vanno confusi → risolve l'apparente inversione regime↔soglia. `\chA`.

### 🔴 A3 — "α=2.39 sta tra 2.72 e 3.0": aritmeticamente FALSO
File: `Section-07-Discussion.tex:147-193`.
- Spieghi il pooled di Taylor (2.39) come "regime mixing" tra grossi vasi (α≈2.72, Wo>2) e arteriole (α≈2.85–3.0), dicendo che 2.39 "sta tra le due". Una media pesata di valori **tutti ≥2.72 non può dare 2.39**.
- Perché il mixing dia 2.4 serve una frazione consistente a **α≈2.0** (limite d'onda), che contraddice "grossi vasi a 2.72" nella stessa sottosezione (e collega ad A4).
- [x] **FATTO 2026-06-20:** Corretto l'errore aritmetico in Section-07:161-201. Il mixing ora copre **[α_w≈2.0, 3.0]** (NON [2.72,3.0]): onda α_w≈2.0 nei conduit più grandi (Wo≫√3), minimax α\*≈2.72 in zona di transizione (Wo~√3), viscoso →3.0 (arteriole). Aggiunto esplicitamente che un pooled di 2.39 < 2.72 **richiede** un contributo sostanziale dai conduit prossimali a basso α (onda) e NON può venire da un range [2.72,3.0]. Declassato "strong confirmation"/"quantitatively confirmed" → "consistent with / compatible with", con caveat che serve dato Wo-stratificato e che l'esponente diameter-weighted pooled non è direttamente confrontabile con l'attrattore area-preserving. `\chA`. Coerente con A4.

### 🔴 A4 — Kleiber 3/4 richiede α≈2.0 dominante, in tensione con α*≈2.72 universale
File: `Section-07-Discussion.tex:589-627`.
- Con `b = dα/(2d+α)` (Paper III): in 3D `b=3/4` SOLO per α=2.0; per α=2.72 → `b≈0.94`.
- La derivazione di Kleiber richiede che i conduit volumetricamente dominanti stiano a α→2.0; ma è la stessa popolazione (grandi coronarie) che altrove validi a ≈2.7 vs Kassab. Doppia identità "grandi vasi = 2.0 (per Kleiber) / = 2.72 (per validazione)" non riconciliata.
- [x] **FATTO 2026-06-20:** Risolta la doppia-identità in Section-07:621-626. Distinti esplicitamente: i vasi volume-dominanti per Kleiber = **aorta + grandi arterie (Wo≫√3, α→2.0)**, mentre i letti coronarici/polmonari misurati morfometricamente = **medio-calibro (Wo~√3, α≈2.7)**. Non è lo stesso segmento: sono regimi di Wo diversi. La α\*(M) per-organismo (2.629) è il valore effettivo della zona di transizione; entro l'albero α(g) varia 2.0(prossimale)→3.0(distale), e il volume-weighted (dominato dai più grandi a α≈2.0) dà b=3/4. `\chA`. Coerente con A3.

### 🔴 A5 — La massa di transizione M* perde nitidezza (0.84 g vs 30–60 g)
File: `Section-06-Architectural-Transition.tex:308-312` vs headline M*≈0.84 g (abstract, Conclusion, App. C).
- Ammetti che il modello morfometrico completo "sposta la transizione nei 30–60 g": spostamento ~50× che indebolisce le predizioni colibrì (4 g) e neonatale (~1 g), entrambe dipendenti dalla posizione di M*.
- [x] **FATTO 2026-06-20:** Chiarito in Section-06:308. Le predizioni-bandiera (colibrì 4g, neonato ~1g) sono **ancorate al numero di Womersley** (Wo₀ vs Wo_c=√3), NON alla massa di transizione assoluta. M\* è un landmark **model-dependent**: la stima single-scale 0.84g sale a decine di grammi quando taper+asimmetria riscalano la mappa Wo↔massa (Paper II). I test colibrì/neonato dipendono solo dal crossing di Wo₀ rispetto a √3 → **robusti** alla rifinitura di M\*. Il testo già ammetteva il gap 0.84g vs 30-60g; aggiunta la frase che àncora le predizioni al Wo. `\chA`.

---

## B. CONTRADDIZIONI NUMERICHE INTERNE (un referee le trova in mezz'ora)

> Causa comune probabile: valori stantii non rigenerati da un unico `compute.py`.
> Azione trasversale: [x] **FATTO 2026-06-20** — rigenerato `dynamic_variables.tex` da un'unica run di `compute.py` e propagato ovunque (manoscritto, README, figure, tabelle). **TUTTO il blocco B (B1–B12) è ora [x].** Infrastruttura revisione 4-colori (`\chR`/`\chRB`/`\chE`/`\chA` = R1/R2/Editorial/Author) + toggle `\ifhighlight` (highlighttrue=colorato / highlightfalse=pulito) in entrambi i preamboli. Tutti gli edit del blocco B avvolti in `\chA` (magenta = author-initiated). **Build verificato:** `main.pdf` (503 KB) e `supplemental.pdf` (417 KB) compilano con latexmk, **0 "Undefined control sequence"**, nessun errore LaTeX.

### 🟠 B1 — Wo del COLIBRÌ: 7.0 vs 1.76 (+ α 2.72/2.60/2.80; massa 4 g vs 2 g)
- 7.0: `Section-01-Introduction.tex:59-62` + README riga 87.
- 1.76: `Section-06-Architectural-Transition.tex:326-341`, `Section-07-Discussion.tex:577`.
- α predetto: 2.72 (intro) / 2.60 (tab. colibrì) / 2.80 (Disc. "Falsifiable Predictions" riga 730); massa 4 g vs 2 g.
- Nota: con mouse(20g,600)→Wo=2.5, lo scaling dà 1.76 per il colibrì (verificato). **7.0 è incompatibile** col 2.5 del mouse → è semplicemente sbagliato.
- [x] Uniformare a 1.76 ovunque (intro, README, tabelle) + un solo α predetto + una sola massa. **FATTO 2026-06-20:** Wo=1.76, 4 g, α≈2.6 (transition) in README:87, Section-01:60, Section-07:578+730; re-grep a zero.

### ✅ 🟠 B2 — Duty cycle η*: 0.777 vs 0.979 — [x] RISOLTO 2026-06-20 (test shielding η-free, canonico)
- 0.777: `manuscript/figures/fig_phase_diagram.tex` (`\VarEtaStar`) + README riga 51.
- 0.979: `supplements/supplemental.tex:825` ("network-level Lagrangian at η*=0.979").
- **DIAGNOSI 2026-06-20 (non è un semplice value-swap):**
  - `compute.py:687` calcola `alpha_h, eta_h = find_minimax(M0)` ⇒ canonico `\VarEtaStar = eta_h = 0.777` (saddle minimax, via marginal-balance `1/(1+|Δw/Δv|)`).
  - `verify_shielding.py:142` ha `eta=0.979` **hardcoded di default** (origine dello 0.979 nel SM).
  - **PROBLEMA PROFONDO:** i due script NON concordano sul modello. A η*=0.777, `verify_shielding.minimax_alpha` dà **α*=3.0** (degenerato, fissato al bound) per ogni Wo₀; `compute.find_minimax` dà **α*=2.627** allo stesso η*=0.777. Convenzioni di costo diverse.
  - **tab:shielding (SM:843-845) è STANTIA rispetto al suo stesso script:** tabella dice α*=2.7789/2.6645 (Δ=0.0093/0.0088); lo script attuale a η=0.979 dà 2.9312/2.7700 (Δ=0.019/0.014). Inoltre il testo "|Δα*|<0.01" è FALSO (reale max 0.019; soglia script 0.02).
  - **FATTO 2026-06-20:** Reso il test di shielding **η-free e canonico**. Aggiunte a `compute.py`: `lame_factor`, `C_wave_lame` (C_wave multi-armonico con correzione thick-wall Lamé opzionale, profilo h/r 0.08→0.42), `find_minimax_lame` (saddle `C_wave_lame = C_visc`, **nessun η fisso**). Lo 0.979 era l'`eta=0.979` hardcoded di `minimax_alpha` (minimizzazione a η fisso) — **eliminato**. `verify_shielding.py` riscritto come thin layer sopra `compute.py` (DRY, stesse funzioni → stessi numeri della tabella). tab:shielding ora **macro-driven** (`\VarShieldThin/Lame/Delta{Two,Four,Eight}`, `\VarShieldDeltaMax`) → mai più stantia. Ancorato al modello canonico: la colonna "thin" riproduce l'attrattore allometrico (use_lame=False ≡ C_wave): Wo₀=2(9g)→3.000, Wo₀=4(148g)→2.705, Wo₀=8(2369g)→2.663; **|Δα\*|_max=0.0032 < 0.01** (claim ora VERO). Testo SM:832 riscritto: rimosso "η\*=0.979", descritto come saddle η-free, bound = `\VarShieldDeltaMax`. **η\* canonico = 0.814** (`\VarEtaStar`, B4) ovunque; lo 0.979 non esiste più. Rimosso `\VarShieldingErr` (dead). Re-grep 0.979: zero in .tex/script.

### ✅ 🟠 B3 — Attrattore d'onda elastico: due formule diverse
- `(5−p)/2 = 2.115`: `Section-05-Single-Mechanism-Limits.tex:47`.
- `(2+p)/(1+p) = 1.565`: `supplements/supplemental.tex:789` (e dice "derived in the main text", ma il main usa l'altra).
- [x] Decidere la formula corretta dell'attrattore d'onda elastico e renderla unica main↔SM. **FATTO 2026-06-20:** **(5−p)/2 = 2.115 è corretta** — derivata da Moens-Korteweg `c=√(Eh/2ρr)`, `h∝r^p` ⇒ `c∝r^{(p-1)/2}` ⇒ `Z_c∝r^{(p-5)/2}`; matching `Z_p=Z_d/N` con `β=N^{-1/α}` ⇒ `α_w=(5−p)/2`. Coincide con main + tutte le figure (`\VarAlphaWTwo=2.115`). La `(2+p)/(1+p)` del SM (S5:791) era ERRATA (incoerente con la SUA stessa velocità Moens-Korteweg) e l'attribuzione "derived in the main text" falsa. Corretto supplemental:791 → `(5−p)/2`.

### ✅ 🟠 B4 — α* simmetrico: 2.626 / 2.627 / 2.629 / 2.63 / 2.65 — [x] RISOLTO 2026-06-20
- 2.626 (Conclusion table, Disc. param table), 2.627 (`\VarAlphaStarModel`), 2.629 (`\VarAlphaHeteroBase` in S3), 2.63 (tab. ext. SM riga 970), 2.65 (README riga 46/69 "geometric model").
- Problema: la decomposizione delle eterogeneità riporta shift a 0.001 (−0.003, −0.006, −0.009) ma il baseline è ambiguo a livello 0.003.
- **DIAGNOSI 2026-06-20:** la causa è in `compute.py` — `\VarAlphaStarModel=2.627` viene da `find_minimax(M0)` (riga 687) mentre `\VarAlphaHeteroBase=2.629` viene da `find_minimax_hetero(70, A=1, taper=1)` (riga 731). **Con A=1, taper=1 NON c'è eterogeneità ⇒ dovrebbero dare lo STESSO valore**, ma differiscono di 0.002: discrepanza di code-path tra `C_wave` e `C_wave_hetero`. Gli shift (−0.003 ecc.) sono calcolati dal baseline 2.629, non da 2.627.
- [x] **FATTO 2026-06-20:** Causa radice = il global `A_ratio=0.82` guidava `find_minimax`→`\VarAlphaStarModel` (etichettato "symmetric" ma calcolato ASIMMETRICO), mentre `\VarAlphaHeteroBase` usava A=1. **Fix in compute.py:** `A_ratio=1.0` (baseline GENUINAMENTE simmetrico) + nuovo `A_emp=0.85` (unica asimmetria empirica misurata, Kassab1993/Kaimovitz2008). `beta_asymmetric`, `\VarAsymmetryA`, le 3 chiamate hetero, la curva full-model della figura e la sensibilità α_w ora usano `A_emp`. **Risultato:** `\VarAlphaStarModel`≡`\VarAlphaHeteroBase`=**2.629** (un solo valore simmetrico); `\VarAlphaStarElastic`≡`\VarAlphaHeteroBaseElastic`=**2.655**. Shift hetero invariati (−0.003/−0.006/−0.009). **CAMBIO MATERIALE (riportato): η\* 0.777→0.814** (il peso marginale del saddle dipende dall'asimmetria del baseline; con baseline simmetrico è 0.814 — rafforza B12: η\* NON è un duty cycle [0.30,0.50]). Manoscritto: hardcoded `2.626`→`\VarAlphaStarModel` (Sec-07:699, Sec-08:65, decomposizione Sec-07:708), `0.82`/Horsfield→`\VarHeteroAsymmetryA`/Kassab1993 (Sec-07:691), SM:980/1123 `2.63`→`\VarAlphaStarModel`, SM:1101 `0.06`→`0.07`. README: `2.626`→`2.629`, rimosso "2.65 geometric" (assente dal testo), η\* 0.777→0.814 + "duty cycle"→"saddle weight". Tutti gli span hardcoded modificati avvolti in `\chA`. Narrativa elastica Sec-07:505-546 è 100% macro-driven → ora auto-coerente. Re-grep: residui solo legittimi (`\VarAlphaHeteroAsym`=2.626 asimmetrico, `\VarBetaShrewObs`=0.82 beta osservato).

### 🟠 B5 — Offset angolo retinico mal etichettato
File: `Section-09-Retinal-Paradox.tex:76-78`.
- Dici "offset dal limite [90°] di solo 0.20σ", ma 3.6°/18° = 0.20σ è l'offset dal **75°** (Murray 3D), non da 90°. Da 90° l'offset è 11.4° = 0.63σ. (I dati 78.6° in realtà supportano la tua tesi "angoli restano 3D Murray", ma la frase dice il contrario.)
- [x] Correggere: il piccolo offset 0.20σ è da 75°, non da 90°. **FATTO 2026-06-20:** Section-09 ora "departing from the 3D Murray prediction by 0.20σ (gap to 90° is 0.63σ)".

### 🟠 B6 — "Wo_c² = d" presentata come legge di embedding ma vale solo a d=3
File: `supplements/supplemental.tex:196,926`.
- `√(6/(d−1)) = √d` solo per d=3. A d=2 darebbe √2, ma usi √6. La presentazione come "isotropic power flow embedding condition in d dimensions" è fuorviante (coincidenza a d=3).
- [x] Riformulare come coincidenza a d=3, non come legge generale. **FATTO 2026-06-20:** supplemental.tex:194 e :928 ora `Wo_c²=6/(d−1)`, forma `=d` dichiarata coincidenza a d=3, non generalizza.

### 🟠 B7 — La soglia fluida da energia dà 2.35, non √3
File: `supplements/supplemental.tex:188-208` (S2).
- L'equipartizione cinetico-viscosa dà `Wo_c^eq ≈ 2.35`, poi scartata e sostituita dalla "geometria" √3. Presentarla tra le "conferme indipendenti" di √3 è fuorviante (dà 2.35, non 1.73).
- [x] Chiarire che l'equipartizione NON corrobora √3 (le conferme reali sono 1.740 esatta vs 1.732 asintotica). **FATTO 2026-06-20:** supplemental.tex:190 riformulato — 2.35 è "crude estimate" che NON coincide con √3 e non lo giustifica.

### 🟠 B8 — Cross-ref Onsager errato
File: `Section-03-Gauge-Invariance.tex:89` ("see Supplemental Material, Section S5 for full derivation").
- S5 è "Limits of the Thin-Wall Approximation". La prova completa di linearità/Onsager è in **S7** (Mathematical Foundations / segmentation).
- [x] Correggere il rimando S5 → S7. **FATTO 2026-06-20:** verificato che S7 ("Mathematical Foundations of Optimization and Invariance") è la sede reale della prova di unicità/linearità; Section-03:89 ora punta a S7.

### ✅ 🟡 B9 — Conteggio reti: 27 vs 28 vs 30
- "30 distinct networks" `Section-06:373`; tabella SM "27 Cases"; testo SM "3 out of 28".
- Anche "Table S2" citato ma label reale `tab:extended_valid`.
- [x] Uniformare il numero e i riferimenti tabella. **FATTO 2026-06-20:** ground truth = **27** (la tabella `tab:extended_valid` ha fisicamente righe 1→27). Corretti: Section-06:373 "30 distinct"→"27 distinct"; Section-06:410 "30+ systems in tab:multivalid" (la tab:multivalid del main ha solo 4 righe!) → "27 systems in the extended validation set (SM, Section S6)"; supplemental:1129 "3 out of 28"→"3 out of 27". La fragile "Table S2" (la tabella è in realtà la 4ª del SM, e l'autore numera "Table SN" a mano) sostituita con il riferimento di sezione robusto **"Section S6"** (la 27-case è una subsection di `sec:S6`, già citata altrove come "Sections S6.1–S6.2"). Re-grep globale: zero residui.

### ✅ 🟡 B10 — Attrattore statico α_t: 2.89 / 2.90 / 2.92 / 2.95 / 3.0 a intercambio
- Bound `(5+p)/2≈2.89` (`Section-04:101`); `\VarAlphaT=2.920`; "α_t≈3.0" (`Section-06:77`); E10 a 2.95 (tab. ontogenetica). L'uso di "3.0" (Murray) mina il tuo stesso raffinamento a 2.90.
- [x] Distinguere chiaramente "Murray 3.0" dal "static optimum 2.90" e usarli in modo consistente. **FATTO 2026-06-20:** principio — il simbolo **α_t = attrattore statico raffinato ≈ 2.92** (`\VarAlphaT`, banda 2.90–2.94, bound (5+p)/2≈2.89<α_t<3); **3.0 = Murray idealizzato**, usato solo per sistemi statici non-pulsatili (xylem/airway/fungi nella tab. SM) ed etichettato come tale. Corretti i 4 usi che chiamavano "α_t≈3.0": Section-06:77 e :91 (teorema+prova, contraddicevano la PROPRIA figura a 2.92)→`\VarAlphaT`; Section-06:355 (caption)→`\VarAlphaT`; Section-07:411 ("transport-**metabolic** attractor")→`\VarAlphaT`; Section-09:89 "α_t=3.0 (Murray)"→rietichettato "the classical Murray exponent (α=3.0)" (riferimento Murray genuino, ma non è α_t). E10 ontogenetico 2.95 (>banda 2.94) → 2.90. Re-grep globale (incl. SM): nessun residuo erroneo; i 3.00 SM sono dati di sistemi statici (corretti).

### ✅ 🟡 B11 — "Wo₀" di riferimento incoerente
- Aorta umana Wo₀≈17.9 (`supplemental.tex:225`) vs 6.8 (tab. colibrì `Section-06:327`). Il vaso/generazione di riferimento per "Wo₀" non è mai definito → predizioni falsificabili ambigue.
- [x] Definire esplicitamente il vaso di riferimento per Wo₀ e rendere coerenti i valori delle tabelle. **FATTO 2026-06-20:** canonico **Wo₀ = numero di Womersley alla radice aortica (gen. 0)** — confermato in 6 punti (Section-06:117/119, Section-07:571, compute.py:170, verify_shielding.py:90, supplemental:227). Il **6.8** della tab. colibrì era ERRATO (≈ valore delle grandi arterie r₃, Wo₃≈6.0): tre calcoli indipendenti danno l'aortica umana ≈18 (diretto r₀=12mm×1490=17.9; scaling dal mouse 2.5×3500^{3/8}×(70/600)^{1/2}≈18.2). Corretto 6.8→**17.9** (coerente col supplemento) e aggiunta definizione esplicita nella caption della tab. colibrì. Re-grep: nessun 6.8 residuo.

### ✅ 🟠 B12 — Collisione del simbolo η ("duty cycle") — [x] RISOLTO 2026-06-20
File: `Section-07-Discussion.tex:62-74` (cardiaco 0.30–0.50) vs η* minimax (0.78/0.98).
- Stesso simbolo e nome per due grandezze diverse: la frazione sistolica cardiaca e il peso dell'avversario del minimax. L'"empirical signature" (η varia 40%, α fisso) usa quello **cardiaco** → sembra un gioco di parole.
- **DIAGNOSI 2026-06-20:** lo stesso simbolo η è usato come peso del Lagrangiano con DUE range diversi: `[0.30,0.50]` (duty cycle cardiaco, Section-07 proposizione + uncertainty set :80-83) e `[0,1]` (peso costo-mixing, fig_phase_diagram: η→1 onda, η→0 trasporto, saddle η*≈0.777). η*=0.777 è **fuori** [0.30,0.50] ⇒ NON può essere un valore del duty cycle cardiaco. Quindi non è solo notazione: è un nodo di **interpretazione del modello** (il peso d'avversario è il duty cycle fisiologico oppure un peso libero?), accoppiato a B2.
- [x] **FATTO 2026-06-20:** Tre grandezze distinte, ora separate: **(1)** `η` = peso costo-mixing astratto ∈[0,1] (Lagrangiano, asse fig_phase_diagram — invariato, già corretto); **(2)** `η_d` = duty cycle cardiaco = t_sys/T ∈[0.30,0.50] (fisiologico, il range dell'avversario); **(3)** `η*` = **peso marginale d'equilibrio** emergente (Eq. eta_star = (1+|C'_w/C'_t|)⁻¹) = **0.814** (B4) — diagnostico DERIVATO, NON un duty cycle, NON vincolato a [0.30,0.50]. Il modello era già corretto (α* dove C_wave=C_transport è robusto al duty cycle; η* è il peso emergente): era solo un'etichetta sbagliata ("duty cycle η*"). Corretto: Section-04:32/73/88 "emergent duty cycle η\*"→"marginal-balance weight η\*"; Section-05:8 "duty cycle"→"cost-mixing weight", :128/129 →"marginal-balance weight"/"cost-mixing weight"; Section-07:39 intro η (cost-mixing weight, motivato da η_d), :44/52/230 η\*→"marginal-balance weight", :62-71/90/146/206 duty cycle cardiaco→`η_d`, :421 "duty-cycle ratio"→"marginal-balance ratio"; **aggiunta frase chiarificatrice** Section-07:54-59 (η_d fisiologico [0.30,0.50] motiva il peso, α\* robusto alla sua variazione; η\*=0.814 è diagnostico derivato non vincolato al range). Section-03:499/564 "metabolic duty cycle"→"metabolic scale"/"normalization" (terzo uso improprio, gauge). fig_phase_diagram già corretto (η∈[0,1], η\* al saddle). SM:781 "duty cycle 1/3" = forma d'onda del flusso (concetto diverso, lasciato). Tutti gli span avvolti in `\chA`. Re-grep "duty cycle $\eta^*$": zero.

---

## C. GAP DI RIGORE / CLAIM DA RIDIMENSIONARE

### 🟠 C1 — Corollario informazionale conta male i bit
File: `Section-02-Scaling-Conflict.tex:166-188`.
- Per sapere a quale **generazione** sei servono ~`log₂G` bit (G livelli), non `g·log₂N` (= indirizzo completo del nodo). La cellula, per regolare μ(g), ha bisogno solo dell'indice di generazione. Il "bound di Shannon" è retorico, non un argomento di canale.
- [x] **FATTO 2026-06-20:** Riscritto `cor:information` (Section-02:166-188) → "Information-Accounting Argument". Corretto il conteggio: l'indice di generazione costa **~log₂G bit** (uno fra G livelli), NON `g·log₂N` (indirizzo completo del nodo). L'ostacolo non è la *magnitudine* ma la **non-località**: l'indice di generazione non è ricostruibile dai soli rapporti emodinamici locali (invarianti per riscalamento globale dell'albero), e lo shielding topologico ne blocca la trasmissione dalla radice. Rimosso il "Shannon communication bound" (era retorico). `\chA`. README riallineato (C1↔D1).

### 🟡 C2 — "Mathematical Proof of Global Convexity" è locale + numerica
File: `supplements/supplemental.tex:1233-1250` (S7 convexity, "Lemma 1").
- L'argomento è locale (vicino al minimo) + tabella numerica Hessiana; non è una prova globale.
- [x] **FATTO 2026-06-20:** Declassato in supplemental.tex (S7 convexity). Titolo "Mathematical Proof of Global Convexity" → "**Local Convexity Argument and Numerical Verification**"; "We first present the mathematical proof" → "argument valid in a neighbourhood of the minimum... we do not claim a closed-form global proof"; "Lemma 1" → "**Local argument (near the minimizer)**" con esplicito "This argument is local; global strict convexity over [2,3] is established by the numerical Hessian scan". `\chA`.

### 🟠 C3 — Appendice A: relazione `a = 1/(1−α²) − 2` asserita senza derivazione
File: `Section-Appendix-A-Fragility.tex:16-18`.
- Genera l'amplificazione 27.7× e δa<0.0018 (headline "fragilità"), ma l'origine fisica della "vessel-wall metabolic optimality condition" non è mostrata. (L'aritmetica a valle torna: a=−2.156, σ=2.77, δa=0.0018 verificati.)
- [x] **FATTO 2026-06-20:** Etichettata come **ansatz fenomenologico** in Appendix-A:16-18, non risultato di prima principi. Aggiunta frase: è la mappa a singolo parametro più semplice con il fixed point richiesto e un polo vicino all'esponente fisiologico; usata solo per esporre il meccanismo di amplificazione della sensibilità. La conclusione (fragilità) è qualitativamente robusta alla forma funzionale precisa (ogni mappa liscia con near-pole in quel range dà σ grande). `\chA`. L'aritmetica a valle (a=−2.156, σ=2.77, δa=0.0018) resta valida.

### 🔴/🟠 C4 — `Q⁻¹ = d−1` applicato a Y_c è postulato, non derivato
- L'argomento cinematico/evanescente (Livello 3) vale solo per il modo longitudinale. Riapplicarlo identico a Y_c è asserito (collegato ad A1).
- [x] **FATTO 2026-06-20:** Dichiarato esplicitamente **postulato/ipotesi** in Section-06:188-192. Il criterio isotropo `Q⁻¹=d−1` è derivato in S1 per il modo longitudinale (bulk fluid); la sua applicazione all'ammettenza caratteristica `Y_c` è estesa **per analogia** (base fisica: phase-halving `Y_c∝√Y_L`), adottata come working hypothesis non come risultato derivato. Aggiunto che la conclusione qualitativa (due soglie incommensurabili) segue dal solo phase-halving ed è robusta al valore preciso della soglia d'onda. `\chA`. Coerente con A1 (verifica Bessel interna conferma ≈3/√2).

### 🟠 C5 — La soglia d'onda non entra nei numeri (è decorativa)
- I calcoli di α*(M) usano SOLO `Wo_c=√3` (hardcoded). La soglia 3/√2 non entra mai nel solver → il "doppio-soglia incommensurabile" è decorativo per le predizioni quantitative.
- [x] **FATTO 2026-06-20:** Dichiarato esplicitamente in Section-06 (dopo il paragrafo del gap [√3,3/√2]). Le predizioni **quantitative** (M\*, curva α\*(M), esponenti per-specie) entrano nel solver SOLO tramite la soglia fluida `Wo_c=√3` (crossover viscoso→inerziale). La soglia d'onda `3/√2` fissa il **bordo superiore del gap di incommensurabilità** e la struttura qualitativa dei regimi (incl. il collasso retinico d=2), ma NON compare nell'ottimizzazione numerica. Statement netto: "the headline numbers depend only on √3". `\chA`.

### 🟠 C6 — Dipendenza pesante da Paper I/II/III non pubblicati
- Il valore-bandiera α*=2.72 dipende dal solver coerente di **Paper II**; QUESTO paper da solo dà 2.626 (combacia col ratto, non con 2.72). L'abstract a tratti promette più di quanto questo lavoro dimostri standalone.
- [x] **FATTO 2026-06-20:** Abstract (main.tex) reso standalone-onesto: "The symmetric model derived **standalone** here yields α\*_model≈\VarAlphaStarModel (2.629)... heterogeneities **(resolved by the coherent impedance solver of the companion work)** shift large-mammal values toward \VarAlphaStar (2.72)". Esplicitato che il 2.72 dipende dal solver completo di Paper II, mentre questo lavoro da solo dà 2.629 (in accordo col ratto). `\chA`. Coerente con il remark Section-05 (baseline α_w=2.0) già presente.

### 🟠 C7 — p (esponente parete) = massima sensibilità ED è circolare
- `S_p≈1.11` (la più alta della tabella) e p è ricavato dagli stessi dati Kassab usati per validare. Ammesso ("Circular Validation Transparency"), ma la magnitudine merita più peso di una nota.
- [x] **FATTO 2026-06-20:** Rafforzato il paragrafo "Circular Validation Transparency" (Section-07). Aggiunto: p è il **parametro strutturale più sensibile** della Table~\ref{tab:sensitivity} (|S_p|≈\VarSP = **1.17**, la voce più alta), quindi la circolarità NON è benigna: il valore symmetric-model di α\* è in parte **calibrato** dal dataset porcino, non indipendente da esso. L'accordo va letto come consistency check interno su quel tessuto, non come validazione pienamente indipendente. `\chA`. (Le predizioni falsificabili indipendenti restano non intaccate.)

---

## D. MINORI / PRESENTAZIONE

### 🟡 D1 — README disallineato col manoscritto
File: `README.md`.
- Cita `Section-02-No-Go-Theorem.tex` (reale: `Section-02-Scaling-Conflict.tex`).
- Chiama "No-Go **Theorem**" ciò che il paper ha ammorbidito in "**Proposition** (regole simmetriche)".
- Angolo retinico **71.3°** (paper: 75° calc / 78.6° oss.) → conferma intro "71°" stantio.
- Colibrì **Wo=7.0** (corpo: 1.76) → conferma B1.
- "geometric model α*≈2.65" assente dal testo.
- Stato "Peer review completed, all revisions applied / audit 2026-05-24" da aggiornare.
- [x] **FATTO 2026-06-20:** README riallineato. Filename `Section-02-No-Go-Theorem.tex`→`Section-02-Scaling-Conflict.tex`; "No-Go **Theorem**"→"No-Go **Proposition**" (sottotitolo, lista teoremi → "one proposition, two theorems", blocco struttura "Proposition 1"); status "Peer review completed... audit 2026-05-24"→"Major revision in progress (Transport Phenomena), 2026-06-20"; bit count `G log₂N / Shannon`→`~log₂G, non-local` (C1↔D1, riga 26/48/227); "minimax **duty cycle** η\*"→"minimax **saddle weight** η\*" (B12); human Wo 6.8→**17.9** (B11); α_t 2.90→**2.92** (B10). Già OK da B1/B4/B5: colibrì 1.76, α\*=2.629, η\*=0.814, retina 75°, rimosso "2.65 geometric". (README è markdown → no `\chA`.)

### 🟡 D2 — Figure
- `manuscript/figures/fig_transition.tex:13`: i `ytick` saltano **2.6** (dove cadono 2.627 e 2.72).
- `fig_phase_diagram` usa α_w=2.115; `fig_transition` usa α_w=2.0 → uniformare.
- [x] **FATTO 2026-06-20:** fig_transition.tex: aggiunto **2.6** alla `ytick` principale (rimosso il duplicato 2.60 dagli extra ticks). Per α_w: i due valori NON sono un errore ma una distinzione **già spiegata** in Section-05:22-26 (rigido α_w=\VarAlphaW=2.0 vs elastico α_w=\VarAlphaWFig=2.115; il phase diagram usa l'elastico). Risolto rendendo la distinzione **esplicita** anche nella caption di fig_transition: aggiunto che la linea d'onda mostrata è il baseline rigido (2.0) usato per calcolare la curva symmetric minimax, distinto dall'elastico 2.115 del phase diagram. `\chA`. (Forzare lo stesso numero sarebbe scorretto: la curva è calcolata a α_w=2.0.)

### 🟡 D3 — Riferimenti citati solo come testo, ASSENTI dalla bibliografia
- Hughes 2000 (`Section-07:738`), Nichols 2011, Kassab 1997, Fåhræus 1929, Horsfield 1971 (nel main).
- (Nota build positiva: tutte le chiavi `\cite` usate risolvono nei `.bib`; nessun riferimento mancante a livello di compilazione. Questi sono testo semplice.)
- [x] **FATTO 2026-06-20 (parziale, onesto):** Le citazioni testuali del main erano tutte nella tabella parametri (Section-07:706-715) + Hughes (riga 762). Aggiunte a references.bib **2 voci verificabili**: `nichols2011` (McDonald's Blood Flow in Arteries, 6th ed) e `fahraeuslindqvist1931` (Am J Physiol 96:562-568) — citate con `\cite`. Tabella: "Kassab 1993"→`\cite{kassab1993}` (già in bib), "Kassab 1997" (taper)→`\cite{kassab1993}` (la morfometria reale contiene il taper; evita di fabbricare una voce 1997), "Nichols 2011"→`\cite{nichols2011}`. **Corretto un errore di contenuto:** la riga "$\alpha_w$ (Fåhræus) ... Fåhræus 1929" era mislabellata (α_w è l'attrattore d'onda, NON un effetto Fåhræus) → rietichettata "(wave attractor)", source `\cite{paperII}`, tipo "Theory". Citato `fahraeuslindqvist1931` al primo uso sostanziale dell'effetto (Section-09:87). **Horsfield 1971** NON più presente nel main (rimosso in B4). **Hughes 2000** (riga 762): NON ho fabbricato i metadati (non verificabili offline) — la voce completa va fornita dall'autore; lasciato come attribuzione testuale. Tutti gli span `\chA`.
  - **AGGIORNAMENTO 2026-06-21 (verifica citazioni contro ARTICOLI CITATI):** `fahraeuslindqvist1931` **VERIFICATO** sul PDF reale (Robin Fåhræus & Torsten Lindqvist, Am J Physiol, recv. 6 Dic 1930→1931, p.562) → nomi completi nel `.bib`. **`nichols2011` (McDonald's) RIMOSSO** — non presente in ARTICOLI CITATI, non verificabile dall'autore. Sostituito (scelta utente) con **`huokassab2006`** = "Pulsatile blood flow in the entire coronary arterial tree: theory and experiment", Huo & Kassab, Am J Physiol Heart Circ Physiol **291:H1074–H1087, 2006**, doi:10.1152/ajpheart.00200.2006 — **VERIFICATO** sul PDF (fit migliore: la sua glossary definisce `c=√(1−F10)·c₀`, la stessa wave speed del paper). Tabella Sec-07 `c_wave`: `\cite{nichols2011}`→`\cite{huokassab2006}`. Build ri-verificato: 0 undefined.

### 🟡 D4 — Abstract: indipendenza da "ATP stoichiometry"/"cardiac output" non riflessa in tabella
- La tabella di sensibilità elenca solo b, Q₀, ℓ₀ (no parametro ATP; cardiac output solo via Q₀).
- [x] **FATTO 2026-06-20:** Allineato via mapping esplicito nella caption di `tab:sensitivity` (Section-05). Aggiunto: i tre parametri metabolici della tabella (b, Q₀, ℓ₀) sono **esattamente** le grandezze citate nell'abstract — Q₀ rappresenta la **cardiac output**, b codifica la **stechiometria ATP/ossigeno** — quindi la tabella sostanzia direttamente l'indipendenza dichiarata. `\chA`.

### 🟡 D5 — Intro "Critical Falsifiable Predictions" con valori stantii
- Colibrì Wo≈7.0→α 2.72; retina θ≈71°.
- [x] **FATTO 2026-06-20 (già via B1/B5):** Verificato che l'intro "Critical Falsifiable Predictions" (Section-01:53-77) è già allineata: colibrì Wo≈**1.76**, α\*≈2.6 (transition), angolo retinico via macro `\VarAngleRetinalCalcDeg`. Nessun valore stantio residuo. `\chA` già presenti.

---

## E. NOTE POSITIVE (da preservare)
- Tutte le chiavi `\cite` usate risolvono nei rispettivi `.bib` (nessun missing citation a compilazione).
- PDF in `output/` esistono → compila.
- Script di verifica presenti (`verify_heterogeneity.py`, `verify_shielding.py`, `generate_womersley_figure.py`).
- Algebra verificata e corretta in molti punti: Q⁻¹=6/Wo²; radice 1.740; Wo₀∝M¹ᐟ⁴; ratio M*(d=2)/M*(d=3)=4; colibrì 1.76; coeff. 1490 m⁻¹; r_crit 1.16 mm; App. A (a=−2.156, σ=2.77, δa=0.0018); algebra di `Wo_c^wave=√(3d(d-2)/(d-1))` (dati i loro assunti).

---

## F. PRIORITÀ D'AZIONE (ordine consigliato)
1. **A1** — risolvere la derivazione della soglia d'onda con check numerico Bessel di `Q⁻¹_{Y_c}` (esistenziale). Prima lo script, poi il testo.
2. **A2 / A3 / A4** — logica regime-2D; aritmetica del 2.39; tensione α=2.0 ↔ 2.72 (Kleiber vs universalità).
3. **B1–B12** — eliminare le contraddizioni numeriche rigenerando i `dynamic_variables` da un'unica `compute.py` e propagando (i più imbarazzanti: colibrì 7.0/1.76; η* 0.777/0.979).
4. **C1–C7** — ridimensionare i claim "rigorosi" che sono ansatz/argomenti locali; dichiarare che la soglia d'onda non entra nei numeri.
5. **D1–D5** — README, figure, riferimenti mancanti.

### Suggerimento di partenza a basso rischio / alto impatto
- Blocco **B (allineamento numerico)** + **D3 (riferimenti mancanti)**: meccanico, sicuro.
- Blocco **A1**: delicato → prima implementare lo script di verifica numerica di `Q⁻¹_{Y_c}` per decidere se il valore corretto è 2.121 o ~2.83 (e se d=2 dà 0 o ∞), POI toccare il testo.

---

## G. RESIDUI dall'audit interno "Opus Pre-Giacomin" (ex-`TODO.md`) — DA FARE

> **CONSOLIDAMENTO 2026-06-21:** questo file è ora l'**unica** lista azioni. Unisce il
> referee report (ex-`TODO 2.md`, Blocchi A–F: **tutti [x]**, build verificato) e l'audit
> interno "Opus Pre-Giacomin" (ex-`TODO.md`). I punti di quest'ultimo già coperti dai
> Blocchi A–D **non** sono ripetuti; restano aperti SOLO i 6 seguenti (verificati con grep
> contro il manoscritto corrente). Le copie in `Progetti IA\branching papers\` sono superate.

- [x] **G1 — Onsager (Thm 6): "lineare al prim'ordine" è in realtà quadratico.** `(Φ−Φ*)/Φ*` è O(‖δx‖²); chiamarla dipendenza lineare al prim'ordine è errato. Correggere la terminologia (Section-03 ~263/306/317; supplemental S7). *Piccolo.*
  - **FATTO 2026-06-21 (🟣 `\chA`, Section-03):** Chiarito che la penalità è **lineare nell'eccesso di dissipazione ΔΦ**, non in δx. Riga 52-53 riformulata; eq.~Onsager: resto `O(‖δx‖²)`→`O(‖δx‖³)` (il termine principale ΔΦ/Φ* è già quadratico in δx); riga 81 "linearity holds to first order in ‖δx‖" → "penalty linear in ΔΦ; near optimum ΔΦ is O(‖δx‖²) since ∇Φ(x*)=0". Coerente con Step 1 (Hessiana) e Step 2 ("quadratic in deviations"). Main compila exit 0.
- [x] **G2 — Sion / quasiconvessità non si somma.** "convessa + quasiconvessa ⇒ quasiconvessa" è FALSO. Riformulare: Hessiana numericamente PD ⇒ **convessa** ⇒ Sion banale. (Il C2 ha solo declassato il titolo della prova di convessità, non questa affermazione specifica.) *Piccolo.*
  - **FATTO 2026-06-21 (🟣 `\chA`, Section-04):** Eliminata l'inferenza errata. Proof (38-44): `L_net` ora dichiarata **convex** in α perché ∂²C_wave/∂α²>0 sul range fisiologico (C_wave convessa lì) + C_transport strict convex ⇒ combinazione non-negativa convessa; Sion soddisfatto *a fortiori*. Remark `rem:quasiconvex` (117-119, 143-147): rimossa la frase contraddittoria "While not strictly convex"; ora "stronger property: convex over the physiological range ⇒ quasiconvex". Main compila exit 0.
- [x] **G3 — Analogie di fisica teorica (gauge / Anderson / RG).** Aggiungere disclaimer espliciti "euristica, non teorema". *Piccolo.*
  - **FATTO 2026-06-21 (🟣 `\chA`):** (1) Anderson localization — "exactly as in Anderson localization" → "in heuristic analogy with the onset of Anderson localization … (a qualitative parallel … not a formal correspondence)" sia in `Section-01:176` sia in `supplemental:520`. (2) Gauge — `Section-07:280` "The analogy with gauge theories … is precise" → "is heuristic, not formal—a global scaling symmetry rather than the local gauge symmetry of field theory (cf.\ Remark~\ref{rem:gauge_terminology})". (3) RG — `Section-07:294` aggiunto disclaimer in testa alla sottosezione: vocabolario RG usato euristicamente, nessuna trasformazione RG formale costruita; chiarito che la correzione finite-size α*(G) è invece una *vera* conseguenza della curvatura del saddle. Main+supp compilano exit 0.
- [x] **G4 — Cross-ref equazione di Jensen.** Numerazione incoerente (Eq 38 vs 45 vs 18 tra main e supplement); inoltre è la generalizzazione convessa pesata, non la "Jensen functional equation" standard. Correggere i `\ref`. *Piccolo.*
  - **FATTO 2026-06-21 (🟣 `\chA`):** (a) Rinomina: `Section-03:267` "This is Jensen's functional equation." → "the weighted Jensen (affinity) equation: … equality for all weights …, not the classical midpoint Jensen equation F((x+y)/2)=(F(x)+F(y))/2"; `supplemental:1206` "generalized Jensen Functional Equation" → "generalized weighted Jensen (affinity) equation (or the Cauchy weighted-average equation)". (b) Cross-ref: "Equation (38)" hard-coded → `\eqref{eq:jensen_weighted_supp}` con nuovo `\label` sull'eq. di consistenza compositiva. Grep globale: ZERO riferimenti `Eq./Equation N` hard-coded residui in tutto il progetto. `\eqref` risolve (nessun Reference undefined).
- [x] **G5 — Overloading di β.** β = rapporto di ramificazione / esponente di Kleiber (3/4) / % di taper in punti diversi (Table 2 illeggibile). Rinominare le variabili (es. `b` per Kleiber, `r_ratio` per branching). *Medio.*
  - **FATTO 2026-06-21 (🟣 `\chA`):** β mantenuto come **unico** significato = *branching ratio* (uso dominante e standard nella letteratura vascolare Huang/Kassab; tabella tab:beta_comparison, tabella supplement, `compute.py beta_symmetric/asymmetric` invariati → nessuna rigenerazione macro). Disambiguati i due usi in conflitto: (1) esponente metabolico/Kleiber → `\beta_M` in `Section-01:60` ("$\beta_M=3/4$") e `Section-03:540` ("$\beta_M(\alpha,d)=d\alpha/(2d+\alpha)$ (… not to be confused with the branching ratio $\beta$)"); (2) taper% in `Section-07:546` "$\beta\approx\VarHeteroTaperPercent\%$" → "a per-generation radius reduction of $\approx\VarHeteroTaperPercent\%$" (β rimosso). `\tau,\kappa,\nu,b` erano tutti già occupati → `\beta_M` scelta non-collidente. Ora β denota esattamente una grandezza.
- [x] **G6 — Tabella ontogenetica (Table 7).** Lo scaling reale dei dati è ~M^0.43, non M^{1/4}; la transizione cade ~2 g, non 0.84 g. Chiarire che embrioni ≠ adulti oppure ricalcolare i dati teorici (può toccare `compute.py`/`.dat`). *Medio.*
  - **FATTO 2026-06-21 (🟣 `\chA`, `Section-08` `tab:ontogenetic`):** Verificate entrambe le discrepanze: i dati (M,Wo₀) della tabella implicano pendenza log-log ≈0.43 (E10→E15) e la soglia Wo₀=√3 cade a ~2 g. **Causa corretta:** lo scaling *ontogenetico* (intraspecifico, singolo organismo in crescita) è legittimamente più ripido dello scaling *interspecifico* adulto (Q̇∝M^{3/4}, f∝M^{−1/4} ⟹ Wo∝M^{1/4}) usato per derivare M*≈0.84 g. La traiettoria ripida è anzi *necessaria* alla narrativa "parte viscoso → transita". Aggiunta nota dopo la tabella: chiarisce ontogenetico≠interspecifico, che la schedule è illustrativa (~M^{0.4} embrionale), che la massa di transizione (~2 g) è trajectory-specific e NON deve coincidere con M*, e che il contenuto falsificabile robusto è qualitativo (α: 3.0→α*≈2.72 quando Wo₀ supera la soglia). Nessuna rigenerazione `compute.py`/`.dat` necessaria. Main compila exit 0.
  - **→ Blocco G (G1–G6) COMPLETO.** Prossimo: 0.R2, poi 0.AUDIT, infine E2 (lettere di risposta).

> Nota fix-residuo (mouse 2.741): già risolto — grep a zero nel testo corrente.
> Nota path: gli edit vivono nel repo git `vascular-networks-theory\`; PDF/articoli in `branching papers\`. Confermare quale copia del manoscritto è canonica prima della sottomissione.
