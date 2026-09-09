"""Run every check. Exit status 0 iff all pass."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import test_basic, test_aromas, test_partitioned, test_weil, test_contractions, test_cotangent, test_exotic, test_vdegree, test_n0, test_round4, test_affine, test_spectrum, test_presentation, test_relations, test_aut, test_jettransport, test_order2

MODULES = [("basic (T) checks", test_basic), ("aroma scaling", test_aromas),
           ("partitioned RK", test_partitioned), ("Weil algebra R[e]/(e^3)", test_weil),
           ("contraction trichotomy", test_contractions), ("cotangent / reverse mode", test_cotangent),
           ("exotic class / definability", test_exotic),
           ("v-degree of the defect", test_vdegree),
           ("N_0 stress test", test_n0),
           ("round 4: N_0 law, defect characterisation, coherence", test_round4),
           ("affine rigidity and the collision trichotomy", test_affine),
           ("the closure spectrum", test_spectrum),
           ("presentations and the gluing problem", test_presentation),
           ("relations between algebras", test_relations),
           ("what Aut-equivariance imposes", test_aut),
           ("jet transport: prerequisites and the shape", test_jettransport),
           ("order-2 classification: detector and simultaneous retraction", test_order2)]

def main():
    allok = True
    for title, mod in MODULES:
        print(f"\n=== {title} ===")
        r = mod.main(); allok &= bool(r)
    print("\n" + ("ALL CHECKS PASSED" if allok else "SOME CHECKS FAILED"))
    return allok

if __name__ == "__main__": raise SystemExit(0 if main() else 1)
