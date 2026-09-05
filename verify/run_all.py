"""Run every check. Exit status 0 iff all pass."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import test_basic, test_aromas, test_partitioned, test_weil, test_contractions, test_cotangent, test_exotic

MODULES = [("basic (T) checks", test_basic), ("aroma scaling", test_aromas),
           ("partitioned RK", test_partitioned), ("Weil algebra R[e]/(e^3)", test_weil),
           ("contraction trichotomy", test_contractions), ("cotangent / reverse mode", test_cotangent),
           ("exotic class / definability", test_exotic)]

def main():
    allok = True
    for title, mod in MODULES:
        print(f"\n=== {title} ===")
        r = mod.main(); allok &= bool(r)
    print("\n" + ("ALL CHECKS PASSED" if allok else "SOME CHECKS FAILED"))
    return allok

if __name__ == "__main__": raise SystemExit(0 if main() else 1)
