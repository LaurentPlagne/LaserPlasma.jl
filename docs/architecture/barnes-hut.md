# Barnes-Hut — arbre d'interaction (fiche de sous-système)

Fiche ouverte au démarrage du chantier Barnes-Hut (phase 1, 2026-09-13). Elle fige la
lecture du code C de 2000 et le contrat numérique à reproduire. **Lire ce document avant
de toucher au portage de l'arbre.**

Référence C : `~/these_postdoc/postdoc/eelinux04/C++/new/` (mai 2000, version des slides)
et `eelinux05/...` (mars 2000). Les fichiers de l'arbre (`load.c`, `grav_s.c`,
`multi_tree.c`, `multistep.c`) sont identiques dans les deux copies ; seul `ppbs.c`
diverge (`do_scale_veloc=1` en e04, `0` en e05) — voir §8 pour ce que cela implique.

## 1. Rôle et chemins vivants

- L'arbre accélère la force électrons-électrons (et électrons-ions, les ions étant des
  corps de l'arbre quand `bckgrd_ion==2`). `use_tree=1` + `stepper=3`.
- Chemin vivant pour les runs des figures :
  `multistep_rec` (multistep.c:1556) → `euler1_rec` (:1086) → `build_force_i` (:168) →
  `build_force_ele_tree_i` (multi_tree.c:283) → `hackgrav_i_plasma` (grav_s.c:256) →
  `gravsub_i_plasma` (grav_s.c:327).
- **Morts ou non liés** (ne pas les porter) :
  - `cc.c` n'est **pas** dans le makefile ; il définit un `build_force_ele_tree_i` et un
    `mkbody_tab_rec` fantômes (électrons seuls, sans Ewald).
  - `multistep_tree` (multi_tree.c:147) et `mk_tree_rec` (:62) n'ont aucun appelant.
  - `force_manybody_tree` (ppbs.c:1219) n'est utilisé que par `dodeint`/`verlet`
    (`stepper` 0/1/2) ; les runs des figures passent par `multistep_rec`.
  - `treescan_s` appelle `gravsub_s` en commentaire : la marche ne calcule rien, elle
    **remplit la liste**.

## 2. Structures (defs.h)

- `node` : index, type (`BODY=01`, `CELL=02`), `mass`, `charge`, `pos`, `next`.
  `BODY` ajoute `vel`, `acc`, `phi` ; `CELL` ajoute `rcrit2`, `more`, `subp[8]`, `dip`,
  `quad`.
- L'ordre d'insertion des corps est **électrons 0..N-1 puis ions** (si
  `bckgrd_ion==2`), masse et charge `q_ion` pour les ions, `-1` pour les électrons.
- `nbody` compte électrons + ions ; les listes d'interaction n'existent que pour les
  électrons (`interaction_list[l]`, `l=0..nb_elec-1`).

## 3. Construction (`maketree`, load.c:28)

1. `newtree()` recycle les cellules (liste chaînée `freecell`) ; `cellused` remis à 0.
2. Racine au centre 0 (`CLRV`), `rsize` initialisé à `lenght_cell = L` par
   `init_plasma` (ewald.c:2569 — le `4·rjel` de `init_part` (ppbs.c:342) est le chemin
   cluster) et **doublé** jusqu'à contenir la boîte (`expandbox`, load.c:105) ; `rsize`
   est global et persistant (il ne rétrécit jamais). En plasma `2·xyzmax < L = rsize`,
   donc `expandbox` ne double jamais.
3. `loadbody` descend par `subindex` (octant : bit k si `pos[k] >= centre[k]`), taille de
   cellule divisée par 2 à chaque niveau.
4. `hackcofm(root, rsize)` calcule masse, charge, centre de masse récursivement ;
   **vérifie** que le cm reste dans la cellule (erreur sinon) ; appelle `setrcrit` avec
   le centre géométrique encore en place ; **écrase ensuite `Pos(p)` par le cm**.
5. `threadtree(root, NULL)` installe `Next`/`More` (parcours linéaire : `More` = premier
   enfant, `Next` = frère suivant ou nœud suivant du parent).
6. `usequad>0` → `hackdip(root)` ; `usequad>1` → `hackquad(root)`. Les moments se
   recalculent par somme des enfants et décalage au cm.

### Critère d'ouverture (`setrcrit`, load.c:173)

`ioptions` est lu dans l'input et transformé en chaîne `options` (`mkstring.c`) :
`ioptions==0` → `"bh86"`, sinon `"sw93"`. **Tous les runs archivés ont `ioptions=1`
(sw93) et `theta=0.5`.** Le critère vivant est donc :

```
bmax2 = Σ_k max(dmin_k, psize - dmin_k)²,   dmin_k = cm_k - (centre_k - psize/2)
rc    = sqrt(bmax2) / theta
Rcrit2(cellule) = rc²
```

(Le cm et le centre géométrique entrent tous deux dans `dmin`.) Le `psize` est le côté
de la cellule ; `theta=0` ouvrirait tout (force exacte).

### Ouverture effective (`subdivp_s`, grav_s.c:165)

`dr = Pos(q) - pos0`, puis **image minimale** `BICV` si `do_plasma==1` (composante par
composante, `while` sur ±L/2 ; voir §6 pour le signe), `drsq = |dr|²` **sans
adoucissement**, ouvert si `drsq < Rcrit2` (strict).

### Marche (`treescan_s`, grav_s.c:101)

Itérative sur `Next`/`More` : cellule trop proche → `More` ; sinon on **accepte** le nœud
dans `interaction_list_i` (auto-interaction exclue, `skipself`) et on compte
`n2bterm`/`nbcterm`. `hackgrav_s` ne fait que remplir et compter (le calcul de force
`gravsub_s` est en commentaire).

## 4. Listes d'interaction (`multi_tree.c`)

- `make_interaction_list()` (:72) : pour chaque électron, `hackgrav_s` puis copie de
  `interaction_list_i[0..size-1]` dans un tableau par électron.
- `make_interaction_list_stable(stable)` (:97) : idem, **sauf** pour les électrons
  stables (`stable[l+1]==1`) où `size_list[l]=1` et la copie recopie un pointeur
  périmé (`interaction_list_i[0]` du dernier appel). **Ce stub n'est jamais lu** : les
  forces des stables ne sont pas évaluées (cf. §5).
- `free_interaction_list()` (:127) libère les `nb_elec` listes (stubs compris).

## 5. Articulation arbre ↔ multistep (multistep.c) — sémantique à respecter

- `multistep_rec` (:1616-1624), **une fois par macro-pas** :
  `mkbody_tab_rec(posdeb)` → `maketree(bodytab,nbody)` → `do_tree=1`, `rebuild_tree=0`
  → `make_interaction_list()` (tous les électrons).
- `step_rec` (:1428-1432), **à chaque niveau de récursion** (donc aussi à l'entrée,
  en doublon du `make_interaction_list` ci-dessus) :
  `free_interaction_list()` → `maketree(bodytab,nbody)` → `make_interaction_list_stable(stable)`.
  Les listes sont **gelées** pendant les deux `rk4_rec` du niveau, puis reconstruites au
  niveau suivant avec le masque `stable` courant.
- `euler1_rec` (:1095-1107), **à chaque étage RK4** : `mkbody_tab_rec(posdeb)` remet les
  positions des corps à celles de l'étage, puis, `rebuild_tree` valant toujours 0,
  `hackcofm_nl` / `hackdip_nl` / `hackquad_nl` **rafraîchissent les moments seuls** : la
  topologie (qui est dans quelle cellule) ne bouge pas dans l'étage. `hackcofm_nl` ne
  refait pas la vérification de bornes.
- `build_force_i` (:181-210) : en mode arbre, la force e-e vient de la liste et **la
  somme séparée sur les ions est désactivée** (`fion=0`) : les ions sont dans l'arbre.
- `do_tree` ne redescend à 0 que si `population<0` (:1449) — jamais vrai
  (`population` est un compteur ≥0) : branche morte.
- Le premier `rk4_rec` de `multistep_rec` (avant `step_rec`) consomme les listes du
  `make_interaction_list` initial ; `step_rec` les libère immédiatement et reconstruit.

## 6. Évaluation de la force (`gravsub_i_plasma`, grav_s.c:327)

Pour un nœud accepté q et un point `pos0` :

```
dr   = pos0 - Pos(q)  ; BICV(dr) si plasma ; drsq = |dr|² + eps_tree²
drab = 1/sqrt(drsq)
phii = Charge(q)/drab ; ai = dr·(phii/drsq)      # monopôle adouci
phi0 += phii ; acc0 += ai
pot, fcmono, fcdipx/y/z = calc_ewald_sum(dr)     # usequad==2 ; sinon calc_ewald
acc0 += Charge(q)·fcmono ; phi0 += Charge(q)·pot # correction périodique
si CELL :
  # dipôle (usequad>0) : terme adouci + phi0 += Dip·fcmono (Hessienne Ewald)
  # quadrupôle + Hessienne (usequad>1) : fcorrdip[k] = Dip·fcdip_k ; terme quadrupôle
```

- `usequad=2` pour les deux runs des figures : Hessiennes requises (§phase 2).
- Le rayon non adouci `eps_tree` du monopôle vaut `r_pot` dans les runs archivés
  (0.05 pour mesure79, 0.144 pour mesure157) ; **en mode arbre c'est `eps_tree` qui
  adoucit aussi les ions**, alors que le direct utilise `r_pot` — les deux sont égaux
  dans les runs archivés, ne pas le présumer ailleurs.
- L'auto-interaction est exclue au moment de la **construction de la liste**, pas ici.

## 7. Contrat d'identité numérique (pour la validation)

- Même intégrateur, mêmes tirages, mêmes forces → comparer à la précision machine.
- L'ordre d'accumulation de `gravsub_i_plasma` est la référence (somme séquentielle des
  listes dans l'ordre de la marche) ; le portage doit le conserver pour les tests
  pas-à-pas.
- Les stubs des stables, le `newtree`/recyclage des cellules, `rsize` persistant,
  l'ordre `Next`/`More` et la sémantique « listes gelées / moments rafraîchis » doivent
  être reproduits **à l'identique** : ce sont des choix qui changent le pas-à-pas.
- `energy.dat` sur un run arbre n'est **pas** calculé par l'arbre :
  `make_epot_bkg_tree`/`make_epot_ele_tree` (multi_tree.c:305/:345) forcent `theta=0` et
  appellent les routines directes (la branche `else` est morte). L'énergie est donc
  O(N²) même en mode arbre — c'est une sortie, pas un test de l'arbre.

## 8. Pièges constatés (mesurés)

- **Thermostat des archives** : l'archive `mesure79` (et les scans N=500 vérifiés)
  colle au binaire `do_scale_veloc=1` (e04), pas à e05. Premier pas : ekin
  21379,742 → 20851,91, identique à notre recompilation e04 à 2e-7 près, alors que le
  binaire e05 donne 21383,38. La classification de
  `docs/reprise/2026-09-13.md` §7 (« mesure79 sans rescaling ») était fondée sur les
  valeurs d'ekin mi-cycle, qui sont l'oscillation de quiver (jusqu'à ~7e5) et non la
  température. Aux frontières de cycle l'ekin reste ~2,1e4 et ΔU n'est pas cumulatif :
  **les deux runs des figures ont le thermostat actif**.
- **Divergence archive/recompilation** : dès t=0 le potentiel diffère de 5,9e-3
  (1e-7 relatif) entre l'archive et notre recompilation — graine d'une divergence
  chaotique (ΔE ~1e-2 à t=1,6 ; ~1e2 à t=8,4). Les sous-runs archivés 25.02.00
  (eps=1e-6) et 26.02.00 (eps=1e-7) donnent ΔU cycle 1 de 627 et 484 contre 507 pour
  le run publié (eps=1,748e-3) : le ΔU par cycle a une **dispersion ~25 % due à la
  tolérance d'intégration seule**. La comparaison pas-à-pas doit se faire contre
  **notre oracle recompilé**, pas contre l'archive ; l'archive sert aux figures
  (observables statistiques).
- **`magic number` d'ouverture** : le pas de temps et `eps` changent la population, donc
  les listes et le résultat ; comparer à pas identique (le multistep adaptatif ne donne
  pas le même découpage si on part d'un état légèrement différent).
- **Coût des sorties** : `energy.dat` (direct, O(N²)) coûte ~0,13 s/pas à N=2000 en
  mode arbre ; `multi_distribution`+`make_field_ion_dire` ~0,17 s/pas de plus. Un
  baseline « intégrateur seul » et un baseline « production » diffèrent d'un facteur
  ~1,4.
- **`rsize` ne rétrécit pas** entre deux `maketree` : une boîte qui grandit puis
  rétrécit garde la grande racine.
- ⚠️ **`rsqrt` est `sqrt`** (`real.h:80`) : `drab = rsqrt(drsq)` vaut `r`, pas
  `1/r`. Écrire `1/sqrt` inverse le monopôle (et fausse `dr3inv`/`dr5inv`) — piège
  qui ne se voit pas sur les nœuds BODY seuls si l'on valide contre le direct.
- **Signes opposés entre marche et force** : `subdivp_s` prend
  `dr = Pos(q) − pos0` (seul le carré compte), `gravsub_i_plasma` prend
  `dr = pos0 − Pos(q)` ; les deux passent par `BICV`.
- **Arbre vs direct : effet réel, mais dominé par le chaos.** Mêmes états initiaux
  (nos binaires, N=2000, config mesure79), ΔU cycle 1/2 :
  arbre 686,7/1491,7 vs direct 658,1/1110,7 (thermostat OFF) ;
  arbre 677,6/754,2 vs direct 547,0/715,4 (thermostat ON).
  L'écart cycle 1 va de +4 % (OFF) à +24 % (ON) : le lissage θ=0,5 change le
  pas-à-pas, et la divergence s'amplifie. Une comparaison arbre/direct n'a de sens
  que pas-à-pas contre le même code, ou en moyenne d'ensemble.

## 9. Baseline mesurée (phase 1, darwin arm64, clang -O2)

Configuration mesure79 : N=2000, r_ws=1,442, kT=6,935, F0=79, ω=3, nb_step=32,
eps=1,748e-3, θ=0,5, sw93, `usequad=2`, `eps_tree`=0,05, dumps retirés.

| binaire (sorties) | pas | arbre | direct | gain |
|---|---|---|---|---|
| `ppbs_bench` (aucune) | 30 | 24,00 s (0,80 s/pas) | 79,93 s (2,66 s/pas) | ×3,3 |
| `ppbs_min` (energy.dat) | 30 | 27,94 s (0,93 s/pas) | 83,66 s (2,79 s/pas) | ×3,0 |
| production (energy+temperature+power+deltaU, e04) | 397 | 439,5 s (1,11 s/pas) | 1233,7 s (3,11 s/pas) | ×2,8 |

Les sorties pèsent lourd : `energy.dat` seul +0,13 s/pas (arbre, N=2000) et
`multi_distribution`+`make_field_ion_dire` +0,18 s/pas de plus. Le « ×3 » porte
donc sur l'intégrateur, pas sur le run complet (×2,8).

N=500 (même config) : `ppbs_bench` 30 pas → arbre 5,20 s (0,173), direct 5,50 s
(0,183) → ×1,06 ; `ppbs_min` 91 pas → arbre 15,23 s (0,167), direct 16,19 s (0,178).

**Conséquence** : l'arbre ne sert qu'à N≳2000. Vérifié dans les archives : les scans
N=500 (a00/a01/a10) et N=1000 sont en **DIRE** ; l'arbre n'apparaît qu'à N=2000
(p. ex. `r14i2000e2000q001a10/F0000W0100N020/TREE_RL_ELE`). La campagne des scans
N=500 n'a pas besoin du tree pour la vitesse.

Repro : `ppbs_bench` = `ppbs_min` avec l'appel `energy()` retiré ; compilation
`make CC=clang CFLAGS="-O2 -DDOUBLEPREC -std=gnu89 -Wno-implicit-function-declaration
-Wno-format-invalid-specifier" LDFLAGS=""` ; entrée `ppbs.inp` = input archive mesure79
converti, `Final time` ajusté (`maxnbstep=(int) xfinal/delta_x` tronque **d'abord**
`xfinal`).

## 10. Suite du chantier (phases 2-5 du plan)

2. ✅ **Table d'Ewald + Hessiennes — fait (2026-09-13).**
   `build_ewald_hessian_table(L)` remplit un `FieldJacobian` via
   `ewald_real_hessian` / `ewald_recip_hessian` / `ewald_direct_hessian`
   (= `Sum_Ewald_1/2/direct`, ewald.c:1818/:1925/:2051) ;
   `ewald_force_jacobian` lit les 13 composantes (= `calc_ewald_sum`, ewald.c:2389).
   Validé contre un dump de l'oracle (`calc_ewald_sum` après `build_tab_ewald(16)`) :
   7 points × 2 boîtes (L=8 et L=8.5635), max |écart| ≤ 6e-16 sur les 13 composantes ;
   test dédié dans `runtests.jl` (+ garde explicite sur la table sans Hessienne).
   Coûts mesurés : table 0,34 s contre 0,24 s (monopôle) ; lecture 287 ns contre
   221 ns par appel. ⚠️ `calc_ewald` et `calc_ewald_sum` somment les 8 nœuds dans un
   ordre différent : écart à l'ulp, pas bit-à-bit. Les coordonnées doivent être dans
   `[-L/2, L/2[` (hors borne haute, lecture hors table — jamais atteint en pratique,
   `BICV`/`min_image` replient avant).
3. ✅ **Octree — fait (2026-09-13).** `src/tree.jl` : `Octree` (corps 1..nbody puis
   cellules ajoutées, 0 = aucun), `set_bodies!`/`maketree!` (= `newtree`/`makecell`/
   `expandbox`/`loadbody`/`subindex`/`hackcofm`/`setrcrit` sw93/`threadtree`/
   `hackdip`/`hackquad`). Le recyclage des cellules du C devient une troncature du
   vecteur ; `rsize` persiste et ne fait que croître. Validé contre un dump oracle
   de `maketree` (10 électrons + 10 ions, `test/data/tree_ref.txt`) : ordre DFS,
   profondeurs, type, indice, liens `more`/`next` **identiques**, masses/charges/
   positions/`Rcrit2`/dipôle/quadrupôle à ≤2e-15. 16 `@test` + `@inferred`
   (`runtests.jl`).
4. ✅ **Marche + interaction — fait (2026-09-13).** `build_interaction_list!`
   (= `hackgrav_s`/`treescan_s`/`subdivp_s`, `BICV` inclus), `gravsub_plasma!`
   (= `gravsub_i_plasma`) et `tree_force_i` (= `hackgrav_i_plasma`). Validés contre
   deux dumps de l'oracle (`init_plasma(sample=2)` → `maketree` → listes → `Acc`/`Phi`) :
   - θ=0,5 (listes surtout de corps) : listes **identiques** (taille et ordre des
     indices pré-ordre), forces ≤ 8,9e-16 ;
   - θ=2,0 (cellules à dipôle non nul acceptées, tous les termes exercés) : listes
     identiques, acc ≤ 8,9e-16, φ ≤ 2,2e-16.
   L'état initial `init_plasma!` est reproduit **bit-à-bit**. Références
   `test/data/tree_interl{,_coarse}_ref.txt` ; 54 `@test` + `@inferred`, plus un
   contrôle de proximité arbre/direct (`|acc + F_direct| < 0,05`).
   Reste la phase 5 (intégration `multistep` par dispatch).
5. ✅ **Intégration `multistep` — faite (2026-09-13).** `ForceScheme` +
   `DirectScheme`/`TreeScheme` (listes par électron, masses/charges) et hooks
   `prepare_step!`/`prepare_stage!`/`prepare_level!` dans `multistep.jl` ;
   `hackcofm_nl!`/`refresh_moments!` dans `tree.jl` ; `run_multistep(...;
   use_tree = true)`. La sémantique du §5 est conservée (arbre + listes
   reconstruits par niveau, moments rafraîchis par étage, ions dans l'arbre).
   - Validé contre `ppbs_min` en mode arbre (N=50, ε=0.144, F0=1.57, eps=0.01,
     7 pas) : `energy.dat` ≤ 3,4e-13, `pos_0001` ≤ 2,7e-15, `mome_0001` ≤ 3,1e-15.
   - Vitesse (mesure79, intégrateur seul, **30 pas** depuis t=0, machine au
     repos ; C = `ppbs_bench` clang `-O2`, mono-thread) :

     | N | C (s/pas) | Julia 1 thread | Julia 4 threads |
     |---|---|---|---|
     | 2000 | 0,813 | **0,64** | **0,218** |
     | 4000 | 1,945 | **1,665** | **0,537** |

     Julia est devant en mono-thread (×0,79 et ×0,86 du C) ; le multithread
     (validé **bit-à-bit** 1 vs 4 threads sur 30 pas) donne ×3,0–3,6 vs C. Le
     direct Julia à N=2000 est 2,807 s/pas (C 2,664), donc l'arbre mono-thread
     garde un gain ~×4. ⚠️ Mesurer sur 30 pas, pas sur les 4 premiers (le coût
     adaptatif croît avec les rencontres proches) et vérifier qu'aucun eval ne
     tourne en fond.
   - **SIMD (fait)** : le reader d'Ewald est une mise à jour rang-1 sur les
     13 composantes (poids des 8 nœuds identiques). Une copie **node-major**
     `packed` (16 lignes par nœud : pot, champ, 9 dérivées, paddée) rend chaque
     nœud contigu, et `ewald_force_jacobian_packed` combine `@ntuple`/`@ntuple`
     à 16 lanes — LLVM vectorise en `<2 x double>` (NEON f64), **2×** sur le
     reader, **bit-à-bit identique**. Gain bout-en-bout ~28 % à N=2000 (1,097 →
     0,78). ⚠️ Piège mesuré : la même formulation avec des closures (`ntuple(q ->
     …)`) est 50× plus lente (boxing) ; il faut la forme littérale `@ntuple`.
     Le reader de référence reste `ewald_force_jacobian` (grilles séparées, utilisé
     par les tests et le direct).
   - Coût mémoire : table Hessienne ~0,63 Mo pour L≈8,6 (`packed` en plus des
     13 grilles). Allocations du pas : ~1 Ko (arène plate), GC négligeable ;
     `@inferred` OK sur tout le chemin chaud.
   - **Optimisations structurelles (faites)** : (1) nœuds en **SoA** (`ntype`,
     `mass`, `charge`, `pos::3×n`, `next`, `more`, `subp::8×n`, `rcrit2`,
     `dip::3×n`, `quad::9×n`) — les boucles chaudes lisent des colonnes
     contiguës au lieu de références `Vector{TreeNode}` ; (2) walk à **curseur**
     (pas de `push!`, arène dimensionnée par doublement, `resize!` final) ;
     (3) `usequad` en **paramètre de type** (`Octree{Q}`), branches du noyau
     supprimées à la compilation ; (4) `hackquad!` avec `@ntuple` (9 composantes) ;
     (5) **multithreading par électron** : `fill_lists!` en une arène par thread
     et `euler1_rec!` sur des plages d'électrons (chaque électron n'écrit que ses
     colonnes ; aucun réduction). Validé bit-à-bit 1 vs 4 threads (30 pas).
   - Gain mono-thread de ce lot : N=2000 0,798 → 0,64 ; N=4000 2,073 → 1,665
     (le C reprenait à N=4000 à cause de la localité des nœuds et du walk ;
     c'est corrigé). Avec 4 threads : 0,218 / 0,537 s/pas.
   - Pistes restantes : invariants de type pour `Q` dans `TreeScheme`, seuil
     `PAR_THRESHOLD` à affiner, et — si un jour utile — l'évaluation par lots
     (SIMD inter-électrons) plutôt que par interaction.
   - Allocations : **~1 Ko/pas en régime établi** (pratiquement nulles). Deux
     étapes : (1) réutilisation des objets cellule (`newtree!`/`makecell!`, 9,5 →
     4,6 Mo/pas) ; (2) arène plate — les listes sont empilées dans
     `TreeScheme.arena::Vector{Int}` avec `starts[i]`/`lens[i]` par électron
     (`fill_lists!`), au lieu d'un `Vector{Int}` par électron (4,6 Mo → ~1 Ko).
     `@inferred` OK sur `bicv`, `fill_interaction_list!`, `gravsub_plasma!`,
     `prepare_step!`, `prepare_stage!` (pas d'instabilité détectée). Le reliquat
     ×1,44 vs C est du calcul (reader Hessien, rafraîchissements de moments,
     marche), plus des allocations.
