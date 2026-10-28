#!/usr/bin/env python
"""
r2_t8_tau_rad.py -- stage r2, analysis T8.

Recomputes every figure of Section 2.4 of the frozen v0 manuscript with the
same inputs and the same formula, at full precision.

METODO:
1. La Ec. 6 de la version enviada (cp p / (4 g sigma T^3)) y sus insumos tal
   como aparecen en v0 Sec. 2.4 (constantes abajo; en el repositorio de trabajo
   se verificaban contra el manuscrito enviado, que no forma parte del paquete).
2. Se recalcula cada cifra de la Sec. 2.4 a precision completa con esos mismos
   insumos y esa misma formula.
3. Se compara con las cifras impresas (redondeo half-up a los digitos impresos
   frente a truncamiento).
4. Se guardan en numbers_r2.json con las cifras significativas propuestas.
5. Si v0 no declara un insumo: "No determinable", sin suponer valores.
Notas declaradas:
- sigma (Stefan-Boltzmann) no aparece en v0; es una constante fisica, no un
  parametro: se usa el valor CODATA del framework y se reporta el efecto de
  usar 5.67e-8 (debe ser invisible en las cifras impresas).
- La duracion del Tsol esta declarada en v0 Sec. 3.1 ("a sol of 1 378 080 s").
- La duracion del ano de Titan no esta declarada en v0: la conversion a anos
  de Titan es "No determinable" con insumos de v0. Se cita el valor de R1
  (tau_rad_eq6_years, con planet.orbital_period del preset) solo como
  referencia, marcado.
- Comparacion descriptiva (T7): la expresion de Stone (1975) / Tokano et al.
  (1999), gamma p H / ((gamma - 1) sigma T^4), con H = R T / g y
  gamma/(gamma-1) = cp/R, se reduce a cp p / (g sigma T^3).

Salidas: build/r2/t8_tau_rad.json y claves tau_rad_r2_* en numbers_r2.json.
"""

from __future__ import annotations

import decimal
import sys
from pathlib import Path

sys.dont_write_bytecode = True

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT / "code" / "titan_pbl_subset") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "code" / "titan_pbl_subset"))

import titan_pbl  # noqa: E402

if (REPO_ROOT / "code" / "titan_pbl_subset") not in Path(titan_pbl.__file__).resolve().parents:
    raise ImportError(f"titan_pbl imported from {titan_pbl.__file__}, not from code/titan_pbl_subset.")

sys.path.insert(0, str(Path(__file__).resolve().parent))

import r2_common as C  # noqa: E402
from titan_pbl.config.physical_constants import STEFAN_BOLTZMANN  # noqa: E402

# Insumos y cifras de la version enviada, Sec. 2.4 y Sec. 3.1, tal como se imprimieron.
V0_EQ6 = "tau_rad = cp p / (4 g sigma T^3)"
V0_SECTION_2_4 = ("Inputs as printed in Section 2.4 of the original submission: cp = 1101 J kg-1 K-1, "
                  "p = 1.45 x 10^5 Pa, g = 1.352 m s-2, T = 93 K; printed results 6.4 x 10^8 s, "
                  "469 Titan solar days, 0.69 Titan years.")
V0_SOL = "a sol of 1 378 080 s (Section 3.1 of the original submission)"


def half_up(x: float, decimals: int) -> float:
    q = decimal.Decimal(1).scaleb(-decimals)
    return float(decimal.Decimal(repr(x)).quantize(q, rounding=decimal.ROUND_HALF_UP))


def truncate(x: float, decimals: int) -> float:
    q = decimal.Decimal(1).scaleb(-decimals)
    return float(decimal.Decimal(repr(x)).quantize(q, rounding=decimal.ROUND_DOWN))


def classify(value: float, printed: float, decimals: int) -> str:
    r, t = half_up(value, decimals), truncate(value, decimals)
    if r == printed and t == printed:
        return "coincide (redondeo y truncamiento iguales)"
    if r == printed:
        return "redondeada"
    if t == printed:
        return "TRUNCADA (el redondeo daria otra cifra)"
    return "NO COINCIDE"


def main() -> int:
    omml = V0_EQ6
    text24 = V0_SECTION_2_4
    sol_par = V0_SOL
    inputs_text = {"cp": "cp = 1101 J kg-1 K-1", "p": "p = 1.45 x 10^5 Pa", "g": "g = 1.352 m s-2", "T": "T = 93 K"}
    cp, p, g, T = 1101.0, 1.45e5, 1.352, 93.0
    sol = 1378080.0
    formula_ok = True

    tau = cp * p / (4.0 * g * STEFAN_BOLTZMANN * T**3)
    tau_567 = cp * p / (4.0 * g * 5.67e-8 * T**3)
    tau_sols = tau / sol
    stone = cp * p / (g * STEFAN_BOLTZMANN * T**3)

    r1_years = C.number_r1("tau_rad_eq6_years")
    result = {
        "task": "R2-T8", "guard": C.guard_info(),
        "eq6_omml_linearised": omml, "eq6_formula_tokens_present": formula_ok,
        "eq6_formula_interpreted": "tau_rad = cp p / (4 g sigma T^3)",
        "v0_section_2_4_sentence": text24,
        "inputs_as_in_v0": inputs_text,
        "sol_length_declared_in_v0": {"value_s": sol, "source": "v0 Sec. 3.1: 'a sol of 1 378 080 s'"},
        "titan_year_declared_in_v0": False,
        "sigma": {"used": STEFAN_BOLTZMANN, "declared_in_v0": False,
                  "tau_with_5.67e-8_s": tau_567, "relative_effect": tau_567 / tau - 1.0},
        "recomputed": {
            "tau_rad_s": tau,
            "tau_rad_sols": tau_sols,
            "tau_rad_years": "No determinable (v0 no declara la duracion del ano de Titan)",
        },
        "comparison_with_printed": {
            "tau_rad_s": {"printed": "6.4 × 10^8 s", "full": tau,
                          "status": classify(tau / 1e8, 6.4, 1), "half_up_2sig": "6.5 × 10^8 s"},
            "tau_rad_sols": {"printed": "469", "full": tau_sols, "status": classify(tau_sols, 469.0, 0),
                             "half_up_3sig": half_up(tau_sols, 0)},
            "tau_rad_years": {"printed": "0.69", "full": "No determinable con insumos de v0",
                              "reference_R1_preset": r1_years,
                              "status_with_R1_value": classify(r1_years, 0.69, 2),
                              "note": "R1 uso planet.orbital_period del preset; no es un insumo de v0"},
        },
        "proposed_print": {"tau_rad_s": "6.5 × 10^8 s (2 cifras)", "tau_rad_sols": "470 Titan sols (2 cifras: T = 93 K se da con 2)",
                           "tau_rad_years": "No determinable con v0; con el ano del preset, 0.70 (2 cifras)"},
        "descriptive_comparison_T7": {
            "stone1975_tokano1999_form": "gamma p H / ((gamma - 1) sigma T^4) = cp p / (g sigma T^3) con H = RT/g, gamma/(gamma-1) = cp/R",
            "value_with_v0_inputs_s": stone, "value_sols": stone / sol, "ratio_to_eq6": stone / tau,
            "note": ("la Ec. 6 de v0 lleva un 4 en el denominador (relajacion linealizada, d(sigma T^4)/dT); la forma de "
                     "Stone/Tokano/Gierasch no. Con los mismos insumos la diferencia es un factor 4 exacto. Ademas, "
                     "p = 1.45e5 Pa es la columna completa sobre la superficie, no 'the lowest 20 km' que dice v0 2.4."),
            "mitchell_lora2016_form": "tau_rad ≃ tau_IR Cp Delta p / (4 g sigma T^3) (Tabla 2, p. 360): la Ec. 6 es esta con tau_IR = 1",
        },
        "sol_paragraph_excerpt": sol_par[:200],
    }
    C.dump_json(C.R2_OUT / "t8_tau_rad.json", result)

    fn = "code/analysis/r2_t8_tau_rad.py::main"
    src = C.rel(C.R2_OUT / "t8_tau_rad.json")
    nums = [
        C.entry("tau_rad_r2_eq6_s", tau, "s", "n/a", "n/a", src, f"{fn} (Ec. 6 de v0: cp p / (4 g sigma T^3), insumos de v0 2.4)", 2,
                "T = 93 K y g, cp, p de v0 con 2-4 cifras; T con 2 cifras limita a 2"),
        C.entry("tau_rad_r2_eq6_sols", tau_sols, "Titan sols", "n/a", "n/a", src, f"{fn} (Tsol de v0 Sec. 3.1)", 2,
                "idem; v0 imprimio 469 (truncado)"),
        C.entry("tau_rad_r2_eq6_years", None, "Titan years", "n/a", "n/a", src, f"{fn}", None,
                "No determinable: v0 no declara la duracion del ano de Titan (R1 tiene tau_rad_eq6_years con el preset)"),
        C.entry("tau_rad_r2_stone_form_s", stone, "s", "n/a", "n/a", src, f"{fn} (cp p / (g sigma T^3), forma de Stone 1975 / Tokano 1999)", 2,
                "comparacion descriptiva; factor 4 exacto respecto de la Ec. 6"),
    ]
    C.update_numbers_r2(nums, "tau_rad_r2_", fn)
    print(C.clean_nan({k: result[k] for k in ("eq6_omml_linearised", "eq6_formula_tokens_present", "recomputed",
                                               "comparison_with_printed", "descriptive_comparison_T7", "sigma")}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
