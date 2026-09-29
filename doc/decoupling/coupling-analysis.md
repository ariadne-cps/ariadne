# Accoppiamento delle sottocartelle di source/ — Ariadne

Analisi del 29 settembre 2026, branch `master`, commit `949d7e044ae65837fc02e6387701b10e7c15ddc6`.

## Risultato

Le dieci cartelle non costituiscono attualmente dieci componenti indipendenti. Il grafo degli include diretti ha **55 archi orientati, 13 coppie con dipendenza mutua e una sola componente fortemente connessa contenente tutte le dieci cartelle**. Classificazione finale: **9 dipendenze alte, 35 medie, 11 basse**.

Questo non impedisce di usare repository diversi: si possono mantenere dipendenze forti verso pacchetti versionati. Impedisce invece di ottenere isolamento semplicemente spostando ogni directory in un repository. Occorre rimuovere i ritorni ciclici e definire dipendenze di build e contratti pubblici.

Il caso peggiore per estensione è **function che dipende da algebra**: 142 direttive include, 41 file consumatori, 31 header consumatori e 21 header di algebra distinti. Non è necessariamente il lavoro di refactoring più lungo in assoluto: è il riferimento misurato per calibrare le fasce.

## Perimetro e metodo

Sono stati letti tutti i file C++ delle dieci cartelle e i CMake principali e dei test. La misura include **374 file: 280 header e 94 sorgenti .cpp elencati nei target CMake**, con **1.043 direttive include intercartella**. Sono inclusi tutti gli header, anche non raggiunti dalla build corrente, perché l'installazione esporta indiscriminatamente i file .hpp. I due header aggregatori nella radice di source/ sono esaminati ma non trattati come undicesima cartella.

Sono esclusi dai conteggi sette .cpp non elencati nei target: algebra/dense_differential.cpp; dynamics/1D_pde.cpp e 2D_pde.cpp; geometry/polyhedron.cpp, polytope.cpp e zonotope.cpp; numeric/mpfr_array.cpp.

Il rilevamento rimuove i commenti, legge include tra virgolette o parentesi angolari e risolve percorsi relativi e percorsi rispetto a source/. Conta direttive testuali, incluse ripetizioni; non conta come archi separati le dipendenze transitive. Non valuta le condizioni del preprocessore: gli include condizionali sono dipendenze potenziali. Le dipendenze dai test sono valutate a livello di configurazione della build, senza aggiungerle al grafo di produzione.

Le dipendenze esterne utility/logging/threading e configuration, GMP/MPFR e backend grafici non sono nodi del grafo richiesto. Restano necessarie per estrarre pacchetti utilizzabili.

È un'analisi statica degli include, dei tipi e della build, non una prova di compilazione o di link separato. Non sono stati compilati i componenti né eseguiti i test; gli include candidati alla rimozione richiedono una verifica successiva. Include transitivi, forward declaration, istanziazioni template, simboli linkati e macro possono introdurre dipendenze semantiche non rappresentate da un include diretto. Un arco assente significa «nessun include diretto rilevato», non indipendenza dimostrata.

Sono inoltre presenti due include interni non risolti nell'albero: function/calculus_base.hpp → calculus_interface.hpp e numeric/flt64.hpp → numeric/is_number.hpp. Non generano archi inventati. config.hpp è generato dalla build e resta fuori dal grafo intercartella.

## Soglie

Per ciascuna relazione **consumatore C dipende da fornitore P**, misuro:

- I: numero di direttive include verso P nei file di C.
- F: numero di file distinti di C coinvolti.
- H: numero di header distinti di C coinvolti (compresi .tpl.hpp, .inl.hpp e .decl.hpp).
- T: numero di header distinti di P referenziati.
- W = I + 2H: indice euristico che attribuisce un costo aggiuntivo alla propagazione attraverso gli header.
- R = W / 204: indice relativo al massimo osservato, function ← algebra, con W = 142 + 2×31 = 204.

Il coefficiente e le soglie sono una scelta di questa valutazione architetturale, non una metrica standard né una stima di giorni di lavoro. F, T e penetrazione percentuale sono riportati per distinguere ampiezza e dimensioni del componente.

| Fascia iniziale | Soglia rispetto al caso massimo | Interpretazione operativa |
|---|---|---|
| Basso (B) | R < 5%, ossia W ≤ 10 | Confine ristretto; tipicamente include superfluo, dichiarazioni da estrarre, adattatore localizzato o piccolo nucleo condiviso. |
| Medio (M) | 5% ≤ R < 25%, ossia W tra 11 e 50 | Tipi concreti, algoritmi o responsabilità mescolate; per eliminare la dipendenza servono interventi strutturali. |
| Alto (A) | R ≥ 25%, ossia W ≥ 51 | Dipendenza forte e diffusa nel modello matematico e/o nelle interfacce; conviene mantenerla esplicita o inizialmente co-localizzata. |

La classificazione finale corregge verso M sei relazioni che il volume da solo classificherebbe B: **algebra ← function; foundation ← numeric; foundation ← symbolic; geometry ← solvers; io ← symbolic; symbolic ← geometry**. Le motivazioni sono riportate nella tabella completa. Non promuovo automaticamente tutti gli archi di un ciclo: basta anche un arco basso a chiudere un grande ciclo.

«Disaccoppiare» qui significa rendere il confine gestibile ed eliminare i ritorni verso implementazioni superiori, anche spostando contratti o adattatori. Non significa poter riscrivere senza sforzo la semantica matematica del fornitore. Il grado basso non rende automaticamente l'intera cartella pronta per un repository autonomo.

## Matrice completa

**Riga = cartella dipendente; colonna = cartella fornitrice.** La freccia nel grafo va quindi dalla colonna alla riga. A/M/B sono le fasce finali; — indica assenza di include diretto. * indica promozione motivata da B a M.


| Dipendente ↓ / fornitore → | foundation | numeric | algebra | function | symbolic | geometry | solvers | io | dynamics | hybrid |
|---|---|---|---|---|---|---|---|---|---|---|
| foundation | · | M* | — | — | M* | — | — | — | — | — |
| numeric | A | · | — | — | B | — | — | — | — | — |
| algebra | B | A | · | M* | — | M | — | B | — | — |
| function | B | A | A | · | M | M | — | — | — | — |
| symbolic | B | M | M | M | · | M* | — | — | — | — |
| geometry | M | M | M | M | M | · | M* | M | — | B |
| solvers | B | M | A | A | B | M | · | — | — | — |
| io | — | M | B | M | M* | M | — | · | — | — |
| dynamics | — | M | M | A | M | A | M | M | · | B |
| hybrid | B | M | M | M | M | A | M | M | M | · |

## Confini che meritano attenzione

1. **foundation / numeric / symbolic.** I paradigmi e le logiche sono fondamentali per numeric, ma logical.cpp dipende a sua volta da Sequence, Natural e Symbolic. `symbolic/templates.hpp` include operatori, interi e sequenze numeriche. Spostare tutto questo header in foundation trasferirebbe il ciclo: bisogna separare il nucleo generico/logico dalle specializzazioni numeriche e dalle strutture simboliche superiori. La coppia foundation/numeric è un buon candidato a un primo pacchetto comune, mantenendo sottotarget distinguibili.

2. **algebra / function.** Il legame dominante è algebra → function nel verso del grafo. Il ritorno non è solo cosmetico: graded.hpp implementa compute_procedure e algebra_operations.tpl.hpp conosce TaylorSeries e AnalyticFunction. Valutare di spostare queste estensioni nel livello function; non basta una dichiarazione anticipata.

3. **geometry non è un unico livello.** Interval e Box sono primitivi richiesti da algebra/function, mentre function_set, affine_set e paver richiedono funzioni e risolutori. Anche box.cpp contiene integrazioni con function. Prima di estrarre geometry intera conviene distinguere primitivi geometrici, insiemi funzionali e algoritmi di paving/ottimizzazione. Il fatto che ariadne-core compili direttamente geometry/interval.cpp è un riscontro concreto di questo confine già attraversato.

4. **io / geometry.** La grafica è una dipendenza di ritorno delle strutture matematiche: Drawable e draw sono presenti in box, grid paving e altri tipi. In direzione opposta, figure e drawer conoscono geometria e funzioni concrete. Separare interfacce minime e adattatori di disegno dai backend Cairo/Gnuplot richiede una ristrutturazione distribuita.

5. **dynamics / hybrid.** Il legame hybrid che dipende da dynamics è sostanziale: HybridEnclosure contiene un LabelledEnclosure. Il ritorno dynamics che dipende da hybrid è un solo include in enclosure.cpp, senza uso esplicito di DiscreteEvent rilevato. Anche geometry/list_set.hpp include discrete_location.hpp pur non usandone concretamente il tipo. Sono i primi due candidati per una pulizia verificabile con compilazione e test.

6. **function / symbolic.** Formula, Procedure ed Expression si attraversano; esistono sia template generici riutilizzati in numeric/foundation sia ponti tra espressioni, funzioni e insiemi. Servono un nucleo comune e adattatori separati; il solo numero di include sottostima la propagazione.

Evidenze principali:

- [source/foundation/logical.cpp:32](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/foundation/logical.cpp#L32)
- [source/symbolic/templates.hpp:33](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/symbolic/templates.hpp#L33)
- [source/algebra/graded.hpp:538](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/algebra/graded.hpp#L538)
- [source/algebra/algebra_operations.tpl.hpp:29](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/algebra/algebra_operations.tpl.hpp#L29)
- [source/function/domain.hpp:33](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/function/domain.hpp#L33)
- [source/geometry/function_set.hpp:39](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/geometry/function_set.hpp#L39)
- [source/io/graphics_interface.hpp:32](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/io/graphics_interface.hpp#L32)
- [source/geometry/grid_paving.hpp:108](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/geometry/grid_paving.hpp#L108)
- [source/dynamics/enclosure.cpp:68](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/dynamics/enclosure.cpp#L68)
- [source/geometry/list_set.hpp:38](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/geometry/list_set.hpp#L38)
- [source/hybrid/hybrid_enclosure.hpp:125](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/hybrid/hybrid_enclosure.hpp#L125)
- [CMakeLists.txt:99](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/CMakeLists.txt#L99)

## Build e test: cosa manca per repository separati

- I dieci target di cartella sono OBJECT library; source/CMakeLists.txt non dichiara i loro archi di dipendenza. La radice rende visibile tutto source/ con include_directories e aggrega gli oggetti in tre librerie condivise.
- ariadne-core contiene foundation, numeric, algebra, oggetti esterni e geometry/interval.cpp. ariadne-kernel aggiunge function, geometry, solvers, io e symbolic. ariadne aggiunge dynamics e hybrid. Sono aggregati sovrapposti, non tre livelli già impacchettati con dipendenze esplicite.
- È installato/esportato il target ariadne; non ci sono dieci pacchetti CMake indipendenti. Il config.hpp generato nella directory sorgente e gli include globali sono altri contratti da trasformare in proprietà dei target.
- I test numeric/algebra linkano ariadne-core; function/symbolic/geometry/solvers/io linkano ariadne-kernel; dynamics/hybrid linkano ariadne. Perciò i test attuali non provano il link isolato delle singole cartelle.
- Esiste tests/foundation/CMakeLists.txt, ma tests/CMakeLists.txt non aggiunge foundation: test_logical non è registrato attraverso quella catena.
- Il commit analizzato introduce setup_standalone_project_tests per Ariadne nel suo insieme. È un supporto all'integrazione come progetto, non una separazione dei test delle dieci cartelle.
- utility/logging/threading provengono dall'albero di submodule; la build richiede anche coerenza del commit configuration tra dipendenze. Per ogni futuro repository occorre fissare versioni compatibili e distinguere test di componente da test di integrazione.

Evidenze di build e test:

- [source/CMakeLists.txt:1](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/CMakeLists.txt#L1)
- [source/function/CMakeLists.txt:1](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/function/CMakeLists.txt#L1)
- [CMakeLists.txt:87](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/CMakeLists.txt#L87)
- [CMakeLists.txt:99](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/CMakeLists.txt#L99)
- [CMakeLists.txt:156](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/CMakeLists.txt#L156)
- [CMakeLists.txt:191](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/CMakeLists.txt#L191)
- [tests/CMakeLists.txt:1](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/tests/CMakeLists.txt#L1)
- [tests/numeric/CMakeLists.txt:21](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/tests/numeric/CMakeLists.txt#L21)
- [tests/function/CMakeLists.txt:15](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/tests/function/CMakeLists.txt#L15)
- [tests/dynamics/CMakeLists.txt:19](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/tests/dynamics/CMakeLists.txt#L19)
- [.gitmodules:1](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/.gitmodules#L1)

## Ordine suggerito

1. **Costruire i confini nel monorepo.** Header pubblici controllati, target con dipendenze PUBLIC/PRIVATE/INTERFACE, configurazione generata nel build tree, test linkati al componente e alle sole dipendenze dichiarate. Controlli di compilazione degli header e consumatori esterni di prova. Nessuna migrazione di repository necessaria per iniziare.
2. **Tagliare i ritorni piccoli verificabili.** Prima i due include verso hybrid, poi i contratti dichiarativi e il nucleo dei template simbolici, evitando di creare nuovi cicli. Spostare grafica di Tensor ed estensioni function presenti in algebra dove opportuno.
3. **Definire il nucleo matematico.** Iniziare con foundation+numeric come unità di rilascio, stabilizzare algebra, individuare un modulo di primitivi geometrici. Non assumere che spostare soltanto Interval/Box risolva tutti gli archi: controllare implementazioni e istanziazioni.
4. **Separare responsabilità intermedie.** Insiemi funzionali, solver e rendering/adattatori grafici; separare template simbolici generici dalle conversioni Expression/Formula. È la parte più strutturale del lavoro.
5. **Estrarre repository dopo la verifica dei confini.** dynamics e hybrid possono essere pacchetti applicativi distinti che consumano versioni del nucleo; hybrid manterrà una dipendenza esplicita da dynamics. Come prima prova operativa, hybrid è un candidato dopo la rimozione dei due ritorni dal basso, ma dipenderà ancora da gran parte del nucleo.

Per uno sviluppo parallelo rapido, è ragionevole mantenere inizialmente un repository del nucleo (con sottotarget e test separati) e separare i consumatori applicativi. Creare dieci repository subito moltiplicherebbe i vincoli di versione senza rimuovere l'accoppiamento.

## Tutte le relazioni e relative motivazioni

Il verso scritto è sempre **fornitore → dipendente**, come richiesto. I/H/F sono conteggi statici, non chiamate a runtime. Ogni riga include un esempio verificabile; include-evidence.csv contiene tutte le occorrenze, inclusi gli include interni alla stessa cartella.

| Fornitore → dipendente | Fascia | I / F / H | R | Motivazione ed esempio |
|---|---|---|---|---|

| numeric → foundation | M* | 2 / 1 / 0 | 1.0% | logical.cpp implementa congiunzioni/disgiunzioni su Sequence e usa Natural; ciclo sostanziale anche se limitato a un file. [source/foundation/logical.cpp:32](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/foundation/logical.cpp#L32) |
| symbolic → foundation | M* | 1 / 1 / 0 | 0.5% | LogicalExpression è costruita sui template Symbolic e sugli operatori numerici: estrarre il nucleo logico dei template insieme al riordino di numeric. [source/foundation/logical.cpp:33](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/foundation/logical.cpp#L33) |
| foundation → numeric | A | 53 / 40 / 33 | 58.3% | Logical e paradigm entrano in 40 file e 33 header: forte dipendenza dei contratti numerici. [source/numeric/approximate_real.hpp:39](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/numeric/approximate_real.hpp#L39) |
| symbolic → numeric | B | 2 / 2 / 1 | 2.0% | Solo templates.hpp, da number_wrapper.hpp e real.cpp. Estrarre il sottoinsieme di expression template numerici; evitare di trasferire il ciclo in foundation. [source/numeric/number_wrapper.hpp:55](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/numeric/number_wrapper.hpp#L55) |
| foundation → algebra | B | 3 / 2 / 2 | 3.4% | Paradigmi e dichiarazioni logiche concentrati nei contratti di base; separabili come piccolo contratto condiviso. [source/algebra/declarations.hpp:16](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/algebra/declarations.hpp#L16) |
| numeric → algebra | A | 41 / 24 / 18 | 37.8% | Tipi scalari, operatori, arrotondamento e concetti numerici attraversano algebra e template. [source/algebra/algebra.cpp:25](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/algebra/algebra.cpp#L25) |
| function → algebra | M* | 2 / 2 / 2 | 2.9% | Graded calcola Procedure; le operazioni di algebra compongono TaylorSeries. Servono separazione delle estensioni e revisione dei template. [source/algebra/algebra_operations.tpl.hpp:29](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/algebra/algebra_operations.tpl.hpp#L29) |
| geometry → algebra | M | 6 / 6 / 3 | 5.9% | Interval entra in interfacce, differential, expansion e matrix; occorre estrarre una base per intervalli, non solo cambiare include. [source/algebra/algebra_interface.hpp:39](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/algebra/algebra_interface.hpp#L39) |
| io → algebra | B | 2 / 1 / 1 | 2.0% | Grafica concentrata in Tensor: estrarre interfacce minime/adattatore grafico, verificando ereditarietà e compatibilità. [source/algebra/tensor.hpp:34](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/algebra/tensor.hpp#L34) |
| foundation → function | B | 4 / 3 / 3 | 4.9% | Dichiarazioni, paradigmi e representation concentrati nei contratti; isolabili conservando un pacchetto di base. [source/function/calculus_base.hpp:32](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/function/calculus_base.hpp#L32) |
| numeric → function | A | 57 / 39 / 29 | 56.4% | Scalari, precisioni e aritmetica sono incorporati in API e template delle funzioni. [source/function/affine.cpp:25](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/function/affine.cpp#L25) |
| algebra → function | A | 142 / 41 / 31 | 100.0% | Caso massimo: vettori, matrici, differenziali, espansioni, interfacce e mixin attraversano 41 file e 31 header. [source/function/affine.cpp:29](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/function/affine.cpp#L29) |
| symbolic → function | M | 6 / 5 / 3 | 5.9% | Formula e Procedure condividono template simbolici; function.cpp usa Expression. Occorre separare albero generico e conversioni. [source/function/formula.cpp:27](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/function/formula.cpp#L27) |
| geometry → function | M | 21 / 8 / 6 | 16.2% | Domain usa Interval e Box; multifunction e funzioni misurabili usano set e function_set. Separare domini primitivi da insiemi funzionali. [source/function/domain.hpp:32](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/function/domain.hpp#L32) |
| foundation → symbolic | B | 4 / 3 / 3 | 4.9% | Logical e Tribool in tre header; preservare un contratto comune minimo. [source/symbolic/expression.hpp:45](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/symbolic/expression.hpp#L45) |
| numeric → symbolic | M | 13 / 7 / 7 | 13.2% | Tipi di numero, operatori e sequenze compaiono nei template e nella valutazione simbolica. [source/symbolic/assignment.hpp:43](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/symbolic/assignment.hpp#L43) |
| algebra → symbolic | M | 9 / 5 / 4 | 8.3% | Valutazione delle espressioni su Algebra e predicate su vettori, matrici e differenziali. [source/symbolic/expression.cpp:29](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/symbolic/expression.cpp#L29) |
| function → symbolic | M | 13 / 7 / 5 | 11.3% | Conversioni Expression/Formula, vincoli e function_expression legano il livello simbolico ai modelli funzionali. [source/symbolic/expression.cpp:32](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/symbolic/expression.cpp#L32) |
| geometry → symbolic | M* | 5 / 3 / 2 | 4.4% | Expression_set espone Box e costruisce insiemi funzionali: la piccola superficie nasconde un ponte semantico da separare. [source/symbolic/expression.hpp:40](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/symbolic/expression.hpp#L40) |
| foundation → geometry | M | 12 / 11 / 11 | 16.7% | Predicati logici e paradigmi fanno parte delle interfacce di insiemi e paving. [source/geometry/box.hpp:34](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/geometry/box.hpp#L34) |
| numeric → geometry | M | 24 / 14 / 10 | 21.6% | Intervalli, box e set dipendono da tipi numerici e arrotondamento. [source/geometry/affine_set.cpp:28](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/geometry/affine_set.cpp#L28) |
| algebra → geometry | M | 19 / 12 / 7 | 16.2% | Vettori, matrici e algebra sostengono punti, box, insiemi affini e vincoli. [source/geometry/affine_set.cpp:29](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/geometry/affine_set.cpp#L29) |
| function → geometry | M | 35 / 14 / 6 | 23.0% | Insiemi funzionali/affini e curve memorizzano funzioni e modelli; anche box.cpp usa modelli Taylor. [source/geometry/affine_set.cpp:25](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/geometry/affine_set.cpp#L25) |
| symbolic → geometry | M | 5 / 4 / 3 | 5.4% | Punti etichettati, curve e insiemi funzionali usano identificatori, spazi, assegnamenti e template simbolici. [source/geometry/curve.cpp:41](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/geometry/curve.cpp#L41) |
| solvers → geometry | M* | 5 / 3 / 0 | 2.5% | Paver e insiemi affini/funzionali invocano programmazione lineare/nonlineare e constraint solver: serve un confine algoritmico, pur con pochi include. [source/geometry/affine_set.cpp:31](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/geometry/affine_set.cpp#L31) |
| io → geometry | M | 17 / 14 / 10 | 18.1% | Disegno ed ereditarietà da Drawable sono distribuiti in numerosi tipi geometrici; separare interfacce e renderer. [source/geometry/affine_set.cpp:40](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/geometry/affine_set.cpp#L40) |
| hybrid → geometry | B | 1 / 1 / 1 | 1.5% | Solo list_set.hpp include discrete_location.hpp; nel file resta una dichiarazione anticipata ma non un impiego concreto. Candidato a pulizia locale. [source/geometry/list_set.hpp:38](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/geometry/list_set.hpp#L38) |
| foundation → solvers | B | 3 / 3 / 1 | 2.5% | Tribool concentrato nei solver di vincoli/nonlineari; piccolo contratto logico condivisibile. [source/solvers/constraint_solver.cpp:30](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/solvers/constraint_solver.cpp#L30) |
| numeric → solvers | M | 19 / 15 / 9 | 18.1% | Scalari e precisioni sono presenti in 15/21 file e 9/12 header; alta penetrazione relativa, volume inferiore al caso massimo. [source/solvers/constraint_solver.cpp:31](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/solvers/constraint_solver.cpp#L31) |
| algebra → solvers | A | 46 / 17 / 9 | 31.4% | Algoritmi e API usano algebra, differenziali, vettori e matrici, anche template di implementazione. [source/solvers/bounder.hpp:34](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/solvers/bounder.hpp#L34) |
| function → solvers | A | 62 / 18 / 9 | 39.2% | Integratori, bounder e solver operano su modelli, formule, procedure e patch concrete. [source/solvers/bounder.cpp:26](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/solvers/bounder.cpp#L26) |
| symbolic → solvers | B | 1 / 1 / 1 | 1.5% | Un include di expression_set.hpp in inclusion_integrator.hpp: isolare il ponte per le inclusioni simboliche. [source/solvers/inclusion_integrator.hpp:40](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/solvers/inclusion_integrator.hpp#L40) |
| geometry → solvers | M | 10 / 9 / 6 | 10.8% | Box, intervalli e grid_paving oltre alle dichiarazioni; separare primitivi e funzioni di paving. [source/solvers/constraint_solver.cpp:34](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/solvers/constraint_solver.cpp#L34) |
| numeric → io | M | 8 / 6 / 3 | 6.9% | GraphicsBoundingBoxType e API grafiche espongono tipi numerici; serve una rappresentazione grafica più neutra. [source/io/cairo.cpp:30](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/io/cairo.cpp#L30) |
| algebra → io | B | 3 / 3 / 3 | 4.4% | Solo algebra/declarations.hpp in tre header: estrarre dichiarazioni minime, senza dipendenza dalle implementazioni algebriche. [source/io/drawer.hpp:34](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/io/drawer.hpp#L34) |
| function → io | M | 7 / 7 / 3 | 6.4% | Backend e figure conoscono funzioni concrete; separare disegno delle funzioni dal backend grafico. [source/io/cairo.cpp:31](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/io/cairo.cpp#L31) |
| symbolic → io | M* | 8 / 4 / 1 | 4.9% | Figure etichettate e backend usano Variable, Space ed expression_set; separare la grafica etichettata. [source/io/cairo.cpp:34](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/io/cairo.cpp#L34) |
| geometry → io | M | 12 / 7 / 3 | 8.8% | Drawer e backend manipolano box, punti, paving e insiemi funzionali; separare adattatori. [source/io/cairo.cpp:32](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/io/cairo.cpp#L32) |
| numeric → dynamics | M | 17 / 12 / 9 | 17.2% | Scalari e bounds concreti attraversano header e implementazioni. [source/dynamics/1D_pde.hpp:1](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/dynamics/1D_pde.hpp#L1) |
| algebra → dynamics | M | 25 / 15 / 6 | 18.1% | Enclosure e algoritmi PDE usano vettori, matrici, differenziali e algebra concreta. [source/dynamics/1D_pde.hpp:3](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/dynamics/1D_pde.hpp#L3) |
| function → dynamics | A | 50 / 20 / 9 | 33.3% | Funzioni, modelli Taylor, patch e vincoli sono parte della rappresentazione e degli evolutori. [source/dynamics/differential_inclusion.cpp:27](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/dynamics/differential_inclusion.cpp#L27) |
| symbolic → dynamics | M | 29 / 14 / 7 | 21.1% | Variabili, spazi, assegnamenti ed espressioni definiscono i sistemi continui e le enclosure. [source/dynamics/differential_inclusion.cpp:26](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/dynamics/differential_inclusion.cpp#L26) |
| geometry → dynamics | A | 41 / 14 / 11 | 30.9% | Enclosure, orbite, griglie e raggiungibilità espongono insiemi e paving. [source/dynamics/differential_inclusion.hpp:33](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/dynamics/differential_inclusion.hpp#L33) |
| solvers → dynamics | M | 21 / 12 / 6 | 16.2% | Evolutori legati a integratori, solver di vincoli e configurazioni; occorrono contratti e factory più stretti. [source/dynamics/differential_inclusion_evolver.cpp:30](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/dynamics/differential_inclusion_evolver.cpp#L30) |
| io → dynamics | M | 11 / 7 / 4 | 9.3% | Draw e interfacce grafiche sono presenti negli header di enclosure, orbite e storage; estrarre adattatori. [source/dynamics/enclosure.cpp:64](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/dynamics/enclosure.cpp#L64) |
| hybrid → dynamics | B | 1 / 1 / 0 | 0.5% | Un solo include di discrete_event.hpp in enclosure.cpp; nessun uso esplicito di DiscreteEvent rilevato nello stesso file. Candidato a rimozione, da compilare. [source/dynamics/enclosure.cpp:68](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/dynamics/enclosure.cpp#L68) |
| foundation → hybrid | B | 1 / 1 / 1 | 1.5% | Una dipendenza diretta da Tribool; contratto minimo condivisibile. [source/hybrid/hybrid_set_interface.hpp:37](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/hybrid/hybrid_set_interface.hpp#L37) |
| numeric → hybrid | M | 17 / 13 / 7 | 15.2% | Tempi, enclosure e set usano scalari e bounds concreti. [source/hybrid/hybrid_enclosure.cpp:28](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/hybrid/hybrid_enclosure.cpp#L28) |
| algebra → hybrid | M | 11 / 9 / 5 | 10.3% | Contenitori e strutture algebriche nei tipi ibridi; dipendenza diffusa ma più ristretta di function/algebra. [source/hybrid/hybrid_enclosure.cpp:29](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/hybrid/hybrid_enclosure.cpp#L29) |
| function → hybrid | M | 30 / 17 / 8 | 22.6% | Modelli, funzioni e vincoli entrano in automi, evoluzione ed enclosure. [source/hybrid/hybrid_automaton-composite.cpp:25](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/hybrid/hybrid_automaton-composite.cpp#L25) |
| symbolic → hybrid | M | 28 / 17 / 11 | 24.5% | Spazi, assegnamenti ed espressioni definiscono automi, variabili e insiemi ibridi. [source/hybrid/discrete_location.hpp:34](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/hybrid/discrete_location.hpp#L34) |
| geometry → hybrid | A | 43 / 14 / 8 | 28.9% | Set, paving, griglie e box sono incorporati nella rappresentazione ibrida. [source/hybrid/hybrid_enclosure.cpp:37](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/hybrid/hybrid_enclosure.cpp#L37) |
| solvers → hybrid | M | 12 / 7 / 3 | 8.8% | Evoluzione e simulazione usano integratori e risolutori, con configurazioni esposte. [source/hybrid/hybrid_enclosure.cpp:46](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/hybrid/hybrid_enclosure.cpp#L46) |
| io → hybrid | M | 14 / 8 / 3 | 9.8% | Grafica in moduli dedicati ma anche in set, paving, orbite ed enclosure; separare renderer e contratti. [source/hybrid/hybrid_enclosure.cpp:51](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/hybrid/hybrid_enclosure.cpp#L51) |
| dynamics → hybrid | M | 9 / 8 / 6 | 10.3% | HybridEnclosure contiene LabelledEnclosure; evolutori e reachability riusano infrastrutture continue. Non basta una forward declaration. [source/hybrid/hybrid_enclosure.hpp:51](https://github.com/ariadne-cps/ariadne/blob/949d7e044ae65837fc02e6387701b10e7c15ddc6/source/hybrid/hybrid_enclosure.hpp#L51) |