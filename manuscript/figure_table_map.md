# Mapa de ubicación de figuras y tablas — Paper A

Cada fila dice dónde va el elemento en el manuscrito, qué marcador buscar en el
`.docx`, y qué archivo exacto corresponde. Los marcadores están en línea propia,
en rojo y negrita, para que sean imposibles de pasar por alto.

**Fuente A** = framework propio (`titan_pbl`), diagnósticos sin integración
temporal. **Fuente B** = archivo publicado de Hartwick et al. (2026),
doi:10.5281/zenodo.17246644. La separación es deliberada y debe conservarse.

| Elemento | Sección | Marcador en el texto | Archivo | Fuente |
|---|---|---|---|---|
| Figure 1 | §2.1 | `[[INSERTAR FIGURA 1 AQUÍ — archivo: figures/fig01_schematic.png]]` | `figures/fig01_schematic.png` | A |
| Table 1 | §3.1 | `[[INSERTAR TABLA 1 AQUÍ — archivo: tables/table01_configuration.csv]]` | `tables/table01_configuration.csv` | A |
| Figure 2 | §3.2 | `[[INSERTAR FIGURA 2 AQUÍ — archivo: figures/fig02_heating_profiles.png]]` | `figures/fig02_heating_profiles.png` | A |
| Figure 3 | §3.3 | `[[INSERTAR FIGURA 3 AQUÍ — archivo: figures/fig03_destabilisation.png]]` | `figures/fig03_destabilisation.png` | A |
| Table 2 | §3.4 | `[[INSERTAR TABLA 2 AQUÍ — archivo: tables/table02_parameter_sweep.csv]]` | `tables/table02_parameter_sweep.csv` | A |
| Figure 4 | §3.5 | `[[INSERTAR FIGURA 4 AQUÍ — archivo: figures/fig04_mixing_time.png]]` | `figures/fig04_mixing_time.png` | A |
| Figure 5 | §3.6 | `[[INSERTAR FIGURA 5 AQUÍ — archivo: figures/fig05_validation.png]]` | `figures/fig05_validation.png` | B |
| Table 3 | §3.6 | `[[INSERTAR TABLA 3 AQUÍ — archivo: tables/table03_validation.csv]]` | `tables/table03_validation.csv` | B |
| Figure 6 | §4.4 | `[[INSERTAR FIGURA 6 AQUÍ — archivo: figures/fig06_timeseries.png]]` | `figures/fig06_timeseries.png` | B |
| Table 4 | §4.4 | `[[INSERTAR TABLA 4 AQUÍ — archivo: tables/table04_contrasts.csv]]` | `tables/table04_contrasts.csv` | B |

## Formatos disponibles

Cada figura existe en tres formatos con el mismo nombre base:

| Extensión | Uso |
|---|---|
| `.png` | 600 dpi, para insertar en Word |
| `.pdf` | vectorial, para la versión de producción |
| `.eps` | vectorial, si la revista lo pide |

Las tablas están en `.csv` (dato crudo) y en Markdown dentro de `captions.md`,
listas para pegar y formatear en Word.

## Tablas auxiliares

`tables/table_aux_thresholds.csv` y `tables/table_aux_mixing_times.csv` respaldan
los números citados en §3.3 y §3.5. En el manuscrito aparecen como tablas dentro
del texto, sin número de figura; no requieren inserción aparte.
