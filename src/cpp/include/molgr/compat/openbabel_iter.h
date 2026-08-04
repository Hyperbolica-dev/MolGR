#pragma once

#include <openbabel/obiter.h>

// OpenBabel 3.2 switched these convenience macros to range-based wrappers that
// expand through an unqualified `impl::...` namespace. That works inside
// OpenBabel's own namespace, but breaks when the macros are expanded from our
// nested namespaces during normal project builds. Rebind them to the classic
// fully-qualified iterator types so both 3.1 and 3.2 behave the same way.
#ifdef FOR_ATOMS_OF_MOL
#undef FOR_ATOMS_OF_MOL
#endif
#define FOR_ATOMS_OF_MOL(a, m) for (OpenBabel::OBMolAtomIter a(m); a; ++a)

#ifdef FOR_BONDS_OF_MOL
#undef FOR_BONDS_OF_MOL
#endif
#define FOR_BONDS_OF_MOL(b, m) for (OpenBabel::OBMolBondIter b(m); b; ++b)

#ifdef FOR_NBORS_OF_ATOM
#undef FOR_NBORS_OF_ATOM
#endif
#define FOR_NBORS_OF_ATOM(a, p) for (OpenBabel::OBAtomAtomIter a(p); a; ++a)

#ifdef FOR_BONDS_OF_ATOM
#undef FOR_BONDS_OF_ATOM
#endif
#define FOR_BONDS_OF_ATOM(b, p) for (OpenBabel::OBAtomBondIter b(p); b; ++b)

#ifdef FOR_RESIDUES_OF_MOL
#undef FOR_RESIDUES_OF_MOL
#endif
#define FOR_RESIDUES_OF_MOL(r, m) for (OpenBabel::OBResidueIter r(m); r; ++r)

#ifdef FOR_ATOMS_OF_RESIDUE
#undef FOR_ATOMS_OF_RESIDUE
#endif
#define FOR_ATOMS_OF_RESIDUE(a, r) for (OpenBabel::OBResidueAtomIter a(r); a; ++a)

#ifdef FOR_DFS_OF_MOL
#undef FOR_DFS_OF_MOL
#endif
#define FOR_DFS_OF_MOL(a, m) for (OpenBabel::OBMolAtomDFSIter a(m); a; ++a)

#ifdef FOR_BFS_OF_MOL
#undef FOR_BFS_OF_MOL
#endif
#define FOR_BFS_OF_MOL(a, m) for (OpenBabel::OBMolAtomBFSIter a(m); a; ++a)

#ifdef FOR_BONDBFS_OF_MOL
#undef FOR_BONDBFS_OF_MOL
#endif
#define FOR_BONDBFS_OF_MOL(b, m) for (OpenBabel::OBMolBondBFSIter b(m); b; ++b)

#ifdef FOR_RINGS_OF_MOL
#undef FOR_RINGS_OF_MOL
#endif
#define FOR_RINGS_OF_MOL(r, m) for (OpenBabel::OBMolRingIter r(m); r; ++r)

#ifdef FOR_ANGLES_OF_MOL
#undef FOR_ANGLES_OF_MOL
#endif
#define FOR_ANGLES_OF_MOL(a, m) for (OpenBabel::OBMolAngleIter a(m); a; ++a)

#ifdef FOR_TORSIONS_OF_MOL
#undef FOR_TORSIONS_OF_MOL
#endif
#define FOR_TORSIONS_OF_MOL(t, m) for (OpenBabel::OBMolTorsionIter t(m); t; ++t)

#ifdef FOR_PAIRS_OF_MOL
#undef FOR_PAIRS_OF_MOL
#endif
#define FOR_PAIRS_OF_MOL(p, m) for (OpenBabel::OBMolPairIter p(m); p; ++p)
