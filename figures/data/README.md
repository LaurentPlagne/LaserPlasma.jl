# Historical figure data

These small compressed files came from the 2000 C-program archive. They are included so the Julia scripts in `figures/` work without the original `~/these_postdoc` tree. Their numerical contents are unchanged. The scripts redraw the curves with CairoMakie; large simulation states and the C binaries are not distributed.

| Bundled file | Source below `~/these_postdoc/postdoc/` |
|---|---|
| `plasma.tex.gz` | `eelinux04/tex/plasma/plasma.tex.gz`; slide source and complete figure list |
| `expew3.xmgr.gz`, `expew5.xmgr.gz` | `eeux03/C++/new/decker/xmgr/` |
| `fig5b.xmgr.gz` | `eeux03/C++/new/fig5b.xmgr.gz` |
| `fleur.xmgr.gz` | `eeux03/C++/old/new/ddata/plasmaj/r10i1000e1000q001a00/F0100W0400N010/DIRE_RL_ELE/` |
| `mesure157.xmgr.gz`, `traj157.xmgr.gz` | `eeux03/C++/new/ddata/plasmaj/r14i2000e2000q001a00/F0157W0300N010/TREE_RL_ELE/` |
| `mesure79_power.dat.gz`, `mesure79_deltaU.gz`, `traj79.xmgr.gz` | `eelinux05/C++/new/ddata/plasmaj/r14i2000e2000q001a00/F7900W0300N010/TREE_RL_ELE/` (`power.dat.gz`, `deltaU.gz`, `traj79.xmgr.gz`) |
| `indiv.agr.gz`, `indiv2.agr.gz` | `eelinux04/tex/plasma/` |

**A provenance detail:** the Grace project directory for `mesure157` and `traj157` is marked `a00`, while the identified C input for the `mesure157` run in `eelinux05` is marked `a01` (`r_pot = 0.144`). The Grace project is the source of the plotted values; its directory name alone does not establish the run parameters. See the [figure gallery](../../docs/src/gallery.md).
