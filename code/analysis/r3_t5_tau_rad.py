#!/usr/bin/env python
"""
r3_t5_tau_rad.py -- stage r3, analysis T5.

Table of radiative timescales: own calculation with the v0 Sec. 2.4 inputs in
the two explicit definitions (D21), the duration of the Titan year with its
literal source, the literature values with their units as written, and the
secondary sources that cite Flasar et al. (1981) on the radiative timescale.

METODO Y CRITERIOS (declarados antes de ejecutar):
10.1 Con los insumos de v0 Sec. 2.4 (los de R2/T8): tau_bulk = cp p /
     (g sigma T^3) (Gierasch / Stone / Tokano) y tau_lin = tau_bulk / 4 (Ec. 6
     de v0), en s, Tsol (1 378 080 s), anos terrestres y anos de Titan.
10.2 Ano de Titan: 29.5 anos terrestres, Strobel et al. (2009, p. 242),
     "Titan year (29.5 Earth years)".
10.3 Por fuente con valor de tau_rad (Smith 1981, Tokano 1999, Charnay &
     Lebonnois 2012, Strobel 2009, Mitchell & Lora 2016, Flasar 2014 si esta):
     cita literal con pagina, unidades tal como las escribe la fuente,
     region o insumos; conversion a anos de Titan con su metodo.
10.4 (Solo en el repositorio de trabajo: busqueda en los PDF de las fuentes
     que citan a Flasar et al. 1981. No forma parte de este paquete, que no
     redistribuye los PDF; se conserva la transcripcion de Tokano et al. 1999.)

DECLARACIONES (antes de ejecutar):
- Ano terrestre = ano juliano, 365.25 d = 31 557 600 s (las fuentes no
  precisan; la diferencia con el ano tropico es 2e-5 relativo).
- Ano de Titan = 29.5 anos terrestres (cita literal de Strobel 2009, p. 242).
  El ano del preset (orbital_period = 929 596 608 s = 29.46 anos julianos) se
  reporta como comparacion, no como fuente.
- sigma = CODATA 5.670374419e-8.
- Las citas de 10.3 son transcripciones verificadas con pagina y metodo contra
  las fuentes (los PDF no se redistribuyen).

Salidas: build/r3/t5_tau_rad_table.csv, t5_tau_rad_notes.md, t5_tau_rad.json.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import r3_common as C  # noqa: E402

SIGMA = 5.670374419e-8
EARTH_YEAR_S = 365.25 * 86400.0
TSOL = C.SOL_S


def main() -> int:
    t8 = json.loads((C.R2_OUT / "t8_tau_rad.json").read_text(encoding="utf-8"))
    cp, p, g, T = 1101.0, 1.45e5, 1.352, 93.0
    tau_bulk = cp * p / (g * SIGMA * T**3)
    tau_lin = tau_bulk / 4.0
    if abs(tau_lin - t8["recomputed"]["tau_rad_s"]) > 1e-6 * tau_lin:
        raise SystemExit("tau_lin no reproduce R2/T8")

    # 10.2 Ano de Titan: Strobel et al. (2009, p. 242), "Titan year (29.5 Earth years)".
    year_hits = {"strobel2009": [{"journal_page": 242, "pdf_page": None, "match": "Titan year (29.5 Earth years)",
                                  "context": "transcripcion de la fuente (el PDF no se redistribuye)"}]}
    titan_year_status = "con fuente"
    titan_year_earth = 29.5
    titan_year_s = titan_year_earth * EARTH_YEAR_S if titan_year_earth else None
    preset_year_s = 929596608.0

    def conv(seconds: float) -> dict:
        return {"s": seconds, "Tsol": seconds / TSOL, "earth_years": seconds / EARTH_YEAR_S,
                "titan_years": (seconds / titan_year_s) if titan_year_s else None,
                "titan_years_preset_period": seconds / preset_year_s}

    own = {"tau_bulk": conv(tau_bulk), "tau_lin": conv(tau_lin)}

    # 10.3 Valores de la literatura ---------------------------------------------
    def ty(value, unit):
        if unit == "Titan years":
            return value, "dado en anos de Titan"
        if unit == "Earth years" and titan_year_earth:
            return value / titan_year_earth, f"anos terrestres / {titan_year_earth} (Strobel 2009, p. 242)"
        return None, "sin conversion (ano de Titan pendiente)"

    lit = [
        {"source": "Smith et al. 1981", "page": "p. 168", "definition": "heat content / mean solar heating rate",
         "value_as_written": "138", "units_as_written": "years (Earth years: 'roughly 400 times the value for Earth')",
         "unit_class": "Earth years", "value": 138.0,
         "region_or_inputs": "whole atmosphere; N2, Voyager radio-science surface pressure and temperature",
         "quote": "we find t = 138 years, roughly 400 times the value for Earth's atmosphere",
         "quote_method": "OCR (R2/T7)"},
        {"source": "Tokano et al. 1999", "page": "p. 499", "definition": "Stone (1975): gamma p H / ((gamma-1) sigma T^4), Eq. 8 (p. 497)",
         "value_as_written": "about three", "units_as_written": "Titan years", "unit_class": "Titan years", "value": 3.0,
         "region_or_inputs": "near the surface (GCM profile)",
         "quote": "tau_rad decreases from about three Titan years near the surface",
         "quote_method": "render (R2/T7)"},
        {"source": "Charnay & Lebonnois 2012", "page": "p. 106", "definition": "not given; attributed to ref. 13 = Tokano et al. 1999",
         "value_as_written": "larger than 2.5", "units_as_written": "Titan years", "unit_class": "Titan years", "value": 2.5,
         "region_or_inputs": "troposphere (lower bound)",
         "quote": "the radiative time constant in Titan's troposphere is larger than 2.5 Titan years",
         "quote_method": "capa de texto (R2/T7)"},
        {"source": "Strobel et al. 2009", "page": "p. 242", "definition": "temperature / cooling rate (Fig. 10.5)",
         "value_as_written": "about 10", "units_as_written": "Titan years", "unit_class": "Titan years", "value": 10.0,
         "region_or_inputs": "20 km altitude",
         "quote": "and is about 10 Titan years at 20 km", "quote_method": "capa de texto (R2/T7)"},
        {"source": "Strobel et al. 2009", "page": "p. 237", "definition": "radiative time constant (not defined on this page)",
         "value_as_written": "~300", "units_as_written": "year (Earth years; compared with 'a Saturnian year (= 29.5 year)')",
         "unit_class": "Earth years", "value": 300.0, "region_or_inputs": "10-20 km",
         "quote": "very long radiative time constants (~300 year at 10-20 km, ~60 year at 50 km)",
         "quote_method": "capa de texto (R2/T7)"},
        {"source": "Mitchell & Lora 2016", "page": "pp. 360-361, Table 2",
         "definition": "tau_IR Cp Delta_p / (4 g sigma T^3)",
         "value_as_written": "~200", "units_as_written": "years (Earth years; omega = 2 pi/(29.5 years))",
         "unit_class": "Earth years", "value": 200.0, "region_or_inputs": "troposphere, tau_IR ~ 10, Huygens Delta_p and T",
         "quote": "tau_rad,tropo ~ 200 years, quite a long time compared to the seasonal cycle",
         "quote_method": "capa de texto (R2/T7)"},
        {"source": "Mitchell & Lora 2016", "page": "p. 361, Table 2",
         "definition": "tau_IR,BL Cp Delta_p_BL / (4 g sigma T_BL^3)",
         "value_as_written": "~7", "units_as_written": "years (Earth years)", "unit_class": "Earth years", "value": 7.0,
         "region_or_inputs": "dry boundary layer, Delta_p_BL ~ 500 hPa (~5 km), tau_IR,BL ~ 1",
         "quote": "the dry boundary layer has a radiative relaxation time of tau_rad,BL ~ 7 years",
         "quote_method": "capa de texto (R2/T7)"},
        {"source": "Flasar, Achterberg & Schinder 2014", "page": "NTRS abstract (p. 106 not available)",
         "definition": "radiative relaxation time", "value_as_written": "(no value)",
         "units_as_written": "-", "unit_class": "none", "value": None,
         "region_or_inputs": "troposphere: 'long compared to seasonal time scales'",
         "quote": "In the troposphere, these are long compared to seasonal time scales",
         "quote_method": "capa de texto del resumen NTRS (R2/T7); capitulo completo pendiente"},
    ]
    for r in lit:
        if r["value"] is None:
            r["value_titan_years"], r["conversion_method"] = None, "sin valor numerico"
            r["value_earth_years"] = None
        else:
            r["value_titan_years"], r["conversion_method"] = ty(r["value"], r["unit_class"])
            r["value_earth_years"] = (r["value"] if r["unit_class"] == "Earth years"
                                      else (r["value"] * titan_year_earth if titan_year_earth else None))

    # 10.4 Citas de Flasar et al. (1981): la busqueda en los PDF no forma parte del paquete.
    flasar_cites = []
    tokano_render = {"file": "s0032-0633_2899_2900011-2.pdf", "pdf_page": 7, "journal_page": 499,
                     "context": ("Our values are somewhat smaller than those of Flasar et al. (1981) who calculated tau_rad "
                                 "more precisely by considering the transmittances of the absorbers which are not "
                                 "available for us."),
                     "method": "render (R2/T7; capa de texto desplazada)", "about_radiative_timescale": True}

    # Tabla ------------------------------------------------------------------------
    rows = []
    for name, d, defin in (("tau_bulk (own, v0 Sec. 2.4 inputs)", own["tau_bulk"], "cp p / (g sigma T^3) — Gierasch/Stone/Tokano"),
                           ("tau_lin (own, v0 Eq. 6)", own["tau_lin"], "cp p / (4 g sigma T^3) — linearised, conservative bound")):
        rows.append({"source": name, "page": "-", "definition": defin, "value_as_written": f"{d['s']:.6e} s",
                     "units_as_written": "s", "region_or_inputs": "cp = 1101, p = 1.45e5 Pa, g = 1.352, T = 93 K (v0 Sec. 2.4)",
                     "value_s": d["s"], "value_Tsol": d["Tsol"], "value_earth_years": d["earth_years"],
                     "value_titan_years": d["titan_years"],
                     "conversion_method": f"s / (29.5 x 31 557 600 s) (Strobel 2009, p. 242)" if titan_year_s else "pendiente de fuente",
                     "quote": "", "quote_method": "calculo propio"})
    for r in lit:
        rows.append({"source": r["source"], "page": r["page"], "definition": r["definition"],
                     "value_as_written": r["value_as_written"], "units_as_written": r["units_as_written"],
                     "region_or_inputs": r["region_or_inputs"],
                     "value_s": (r["value_earth_years"] * EARTH_YEAR_S) if r["value_earth_years"] else None,
                     "value_Tsol": (r["value_earth_years"] * EARTH_YEAR_S / TSOL) if r["value_earth_years"] else None,
                     "value_earth_years": r["value_earth_years"], "value_titan_years": r["value_titan_years"],
                     "conversion_method": r["conversion_method"], "quote": r["quote"], "quote_method": r["quote_method"]})
    C.write_csv(C.R3_OUT / "t5_tau_rad_table.csv", rows)
    out = {"task": "R3-T5", "guard": C.guard_info(), "head": C.git_head(),
           "inputs": {"cp": cp, "p": p, "g": g, "T": T, "sigma": SIGMA, "Tsol_s": TSOL, "earth_year_s": EARTH_YEAR_S},
           "titan_year": {"status": titan_year_status, "earth_years": titan_year_earth, "seconds": titan_year_s,
                          "hits": year_hits, "preset_orbital_period_s": preset_year_s,
                          "preset_in_julian_years": preset_year_s / EARTH_YEAR_S},
           "own": own, "literature": lit, "flasar1981_citations_text_layer": flasar_cites,
           "flasar1981_citation_render": tokano_render, "table": rows}
    C.dump_json(C.R3_OUT / "t5_tau_rad.json", out)
    write_notes(out)

    src = C.rel(C.R3_OUT / "t5_tau_rad.json")
    fn = "code/analysis/r3_t5_tau_rad.py::main"
    ents = [C.r3_entry("titan_year_earth_years", titan_year_earth, "Earth years", src, fn,
                       note="Strobel et al. 2009, p. 242" if titan_year_earth else "pendiente de fuente")]
    for k, d in own.items():
        for unit, val in (("s", d["s"]), ("Tsol", d["Tsol"]), ("earth_years", d["earth_years"]),
                          ("titan_years", d["titan_years"])):
            ents.append(C.r3_entry(f"taurad_{k}_{unit}", val, unit, src, fn, note="insumos v0 Sec. 2.4; D21"))
    C.stage_numbers(ents, ("titan_year_", "taurad_"), fn)
    print("ano de Titan:", titan_year_status, titan_year_earth)
    for k, d in own.items():
        print(k, {u: (f"{v:.4g}" if v else v) for u, v in d.items()})
    for r in lit:
        print(r["source"], r["page"], r["value_as_written"], r["units_as_written"][:20], "->", r["value_titan_years"])
    print("citas Flasar 1981 (capa de texto):", [(c["file"], c["pdf_page"], c["about_radiative_timescale"]) for c in flasar_cites])
    return 0


def write_notes(out: dict) -> None:
    ty = out["titan_year"]
    o = out["own"]
    L = ["# R3/T5 — Escala de tiempo radiativa: notas", "",
         "*Generado por `code/analysis/r3_t5_tau_rad.py`. Tabla completa en `t5_tau_rad_table.csv`.*", "",
         "## 10.1 Cálculo propio (insumos de v0 §2.4: c_p = 1101 J kg⁻¹ K⁻¹, p = 1.45·10⁵ Pa, g = 1.352 m s⁻², T = 93 K; σ CODATA)", "",
         "| Definición | s | Tsol | años terrestres | años de Titán |", "|---|---|---|---|---|"]
    for k, lab in (("tau_bulk", "τ_bulk = c_p p/(g σ T³) (Gierasch/Stone/Tokano; valor de la literatura)"),
                   ("tau_lin", "τ_lin = τ_bulk/4 (Ec. 6 de v0; cota conservadora)")):
        d = o[k]
        tyv = f"{d['titan_years']:.4f}" if d["titan_years"] else "pendiente de fuente"
        L.append(f"| {lab} | {d['s']:.4e} | {d['Tsol']:.2f} | {d['earth_years']:.3f} | {tyv} |")
    L += ["", "Año terrestre = año juliano (31 557 600 s).", "",
          "## 10.2 Año de Titán", "",
          f"**Estado: {ty['status']}.** " + (f"{ty['earth_years']} años terrestres." if ty["earth_years"] else "")]
    for key, hits in ty["hits"].items():
        for h in hits:
            L.append(f"- {key}, p. {h['journal_page']} (PDF p. {h['pdf_page']}), capa de texto: \"…{h['context']}…\"")
    L += [f"- Comparación (no es fuente): `orbital_period` del preset = {ty['preset_orbital_period_s']:.0f} s = "
          f"{ty['preset_in_julian_years']:.3f} años julianos. Con él, τ_lin = {o['tau_lin']['titan_years_preset_period']:.4f} "
          "años de Titán (el 0.696 de R1).", "",
          "## 10.3 Valores de la literatura (unidades tal como las escribe la fuente)", "",
          "| Fuente | Página | Valor escrito | Unidades | Región / insumos | Años de Titán | Método |", "|---|---|---|---|---|---|---|"]
    for r in out["literature"]:
        tyv = f"{r['value_titan_years']:.2f}" if r["value_titan_years"] is not None else "—"
        L.append(f"| {r['source']} | {r['page']} | {r['value_as_written']} | {r['units_as_written']} | "
                 f"{r['region_or_inputs']} | {tyv} | {r['conversion_method']} |")
    L += ["", "Citas literales (de los verbatim de R2/T7, con su método):", ""]
    for r in out["literature"]:
        L.append(f"- **{r['source']}**, {r['page']} ({r['quote_method']}): \"{r['quote']}\"")
    L += ["", "Lectura: las cifras de la literatura no son comparables entre sí sin declarar la definición y la región. "
          "Smith (1981) usa contenido térmico / calentamiento solar medio de toda la atmósfera; Tokano (1999) y "
          "Charnay & Lebonnois (2012, que citan a Tokano) usan la forma de Stone cerca de la superficie; Strobel (2009) "
          "usa T / tasa de enfriamiento a una altura; Mitchell & Lora (2016) usan la forma linealizada con τ_IR. "
          "τ_bulk (propio) es la forma de Stone/Tokano; τ_lin es la de Mitchell & Lora con τ_IR = 1.", "",
          "## 10.4 Fuentes que citan a Flasar et al. (1981) sobre la escala radiativa", ""]
    rel = [c for c in out["flasar1981_citations_text_layer"] if c["about_radiative_timescale"]]
    oth = [c for c in out["flasar1981_citations_text_layer"] if not c["about_radiative_timescale"]]
    fr = out["flasar1981_citation_render"]
    L.append(f"- `{fr['file']}` (Tokano et al. 1999), p. {fr['journal_page']}, {fr['method']}: \"{fr['context']}\"")
    for c in rel:
        L.append(f"- `{c['file']}`, PDF p. {c['pdf_page']} (capa de texto): \"…{c['context']}…\"")
    L += ["", f"Menciones de Flasar (1981) sin relación con la escala radiativa en el contexto inmediato: {len(oth)} "
          "(listadas en `t5_tau_rad.json`).", "",
          "Estas citas son respaldo secundario mientras no esté el PDF de Flasar, Samuelson & Conrath (1981)."]
    (C.R3_OUT / "t5_tau_rad_notes.md").write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
