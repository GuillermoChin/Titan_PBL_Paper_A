# Datos externos

Este directorio mezcla dos cosas con estatus muy distinto.

## Versionado en el repositorio

| Archivo | Origen | Uso |
|---|---|---|
| `huygens/HASI_L4_ATMO_PROFILE_DESCEN_TAB_v2.txt` | PDS, HASI Level 4 | Perfil de inicializacion del modelo |
| `Soundings/input_sounding.*` | Archivo de Hartwick et al. (2026) | Sondeos de entrada de mtWRF |
| `input_sounding.Titan_HASI-sfc_1466hPa-0u-gcmsq_corrected` | idem | Copia del sondeo por defecto |
| `make_paper_A_figures.py` | Este proyecto | Figuras y tablas del Paper A |

El perfil HASI esta duplicado a proposito en `titan_pbl/data/huygens/`: esa es
la copia que viaja con el paquete instalado y la que hace que la suite de tests
pase sin depender del directorio de trabajo. Las dos copias deben ser identicas
byte a byte, y `tests/unit/test_paths.py` lo verifica.

> **Nota sobre el nombre.** El producto del PDS se distribuye como
> `HASI_L4_ATMO_PROFILE_DESCEN.TAB`. Aqui se guarda como
> `..._DESCEN_TAB_v2.txt` porque el punto intermedio hacia que la extension
> aparente fuera `.TAB_v2`, lo que confunde a las herramientas que despachan por
> extension. El contenido no se toca.

## NO versionado: el archivo de simulaciones

Los cuatro directorios `Dry Simulations/`, `SPS Simulations/`,
`NSMT Simulations/` y `MRPS Simulations/` (~1.5 GB de NetCDF) estan excluidos
via `.gitignore`. Son dato externo, inmutable y ya citable por DOI; versionarlos
duplicaria un archivo publico y haria el repositorio inmanejable.

**Origen:** Hartwick, Rafkin & Soto (2026), archivo completo de simulaciones,
Zenodo [doi:10.5281/zenodo.17246644](https://doi.org/10.5281/zenodo.17246644).

**Como reponerlo:** descargar y descomprimir de modo que los cuatro directorios
queden como hijos directos de `data/`. El script de figuras los busca ahi por
defecto; `TITAN_DATA_DIR` permite apuntar a otro volumen.

### Advertencia sobre nombres de archivo en Windows

Los wrfout se publican como `wrfout_d01_2007-00001_09:30:00`. Los dos puntos son
ilegales en NTFS, de modo que al descomprimir en Windows quedan como
`wrfout_d01_2007-00001_09_30_00`. `resolve_archive()` en
`make_paper_A_figures.py` tolera ambas formas; no hace falta renombrar nada.

### Que corrida es la baseline seca

El archivo trae varias corridas secas y ninguna se identifica de forma
inequivoca como "la baseline del paper". El discriminante es `TMN`, la
temperatura del subsuelo profundo, impuesta fija en la capa mas honda. Ver
**H-23** en `docs/HALLAZGOS.md` para el inventario completo y la justificacion
del default. Se selecciona con `TITAN_BASELINE_RUN=TMN91|baseline5sol`.

## Uso

```
python data/make_paper_A_figures.py
```

Variables de entorno reconocidas:

| Variable | Efecto |
|---|---|
| `TITAN_DATA_DIR` | Raiz del archivo de Zenodo. |
| `PAPER_A_OUT` | Directorio de salida. Por defecto `./paper_A` |
| `TITAN_BASELINE_RUN` | `TMN91` (default) o `baseline5sol` |
| `SKIP_ARCHIVE` | `1` omite la Fuente B y produce solo Figuras 1-4 y Tablas 1-2 |
