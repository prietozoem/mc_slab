# Monte Carlo transport: slab shielding (example 1) and neutron point source in a water sphere (example 2)

A standalone web app (plain HTML/CSS/JS, no build step, no dependencies) that teaches the basics of Monte Carlo radiation transport.
A parallel beam (or a line source) of monoenergetic photons **or neutrons** hits an infinite slab and the transmitted flux is tallied. Choose the particle with the *Photons / Neutrons* buttons at the top of the control panel.
Interface in English and Spanish (toggle at the top right).

## Modes

- **Automatic**: runs N histories, draws the trajectories of a sample of them and accumulates the tallies.
- **Step by step**: a floating (draggable, resizable) window shows each stage of a history:
  1. total cross section at the photon energy (point marked on the XCOM curve),
  2. cumulative distribution of the flight path with the random number and the sampled distance,
  3. transport to the collision site (or escape through a face),
  4. partial cross sections (coherent, incoherent, photoelectric, pair) with the values marked at the photon energy,
  5. the cumulative vector of partial cross sections with the random number that selects the reaction,
  6. the reaction itself (polar plot of the scattering-angle distribution with the sampled direction).

  Buttons: **Next step** (Space), **Next particle** (P: finishes the current history and starts the next), **Run continuously** (with an adjustable delay).

## Physics and approximations

- Geometry: slab 0 ≤ z ≤ T, infinite in x and y (in the code the thickness axis is called x, and the displayed x, y are code y, z). Source on the entrance face, monodirectional along +z: a **pencil beam** (point) or a **line source**: a segment of length ΔX on the entrance face (z = 0), emission points uniform on it. Views: 2D z–x, 2D z–y (wheel: zoom about the cursor, drag: pan, double-click: reset) and an orbitable 3D view of the full 3D transport.
- Cross sections: NIST XCOM, 1 keV to 100 MeV, log–log interpolation, absorption edges kept (the table rows with repeated energy).
  Materials: water, aluminum, iron, ordinary concrete (NIST composition), tungsten, lead.
- Incoherent scattering: Klein–Nishina on a free electron at rest (Kahn's rejection method).
- Coherent scattering: Thomson × a simple Gaussian form factor (qualitative, forward peaked).
- Photoelectric absorption and pair production absorb the photon locally (no electrons, fluorescence or annihilation photons).
- Tallies at the exit face: transmitted fraction (uncollided and scattered, with relative errors), the uncollided part compared with exp(−Σt·T), the energy spectrum, and the scalar-flux estimator Σ 1/|μ|.

## Run locally

Open `index.html` in a browser, or serve the folder: `python3 -m http.server` and visit http://localhost:8000.

## Publish on GitHub Pages

1. Create a repository and put the contents of this folder in the repository root (so `index.html` is at the top level).
2. On GitHub: **Settings → Pages → Build and deployment → Deploy from a branch**, choose the `main` branch and the `/ (root)` folder.
3. After a minute the app is at `https://<user>.github.io/<repository>/`.

## Neutron mode (ENDF/B-VII.1)

Cross sections come from the evaluated library **ENDF/B-VII.1** (pointwise, 0 K), retrieved from NNDC Sigma and stored in `js/data_n.js` (loaded only when *Neutrons* is selected; about 1.2 MB). The file must sit next to `data.js`; if it is missing the app shows a message and stays in photon mode.

Per nuclide the app uses elastic (MT 2), inelastic (MT 4), radiative capture (MT 102) and "other absorption" (MT 1 − MT 2 − MT 4 − MT 102, i.e. (n,2n), (n,p), (n,α)…, smoothed and treated as absorption). At each collision the nuclide is chosen with probability Σᵢ/Σt and then the reaction with probability σᵢ,ᵣ/σᵢ (an extra *Choosing the nuclide* step appears in step mode).

Simplifications (this is a teaching code): targets at rest, no thermal scattering law (neutron energy floored at 0.0253 eV); elastic scattering isotropic in the centre of mass; inelastic scattering isotropic in the lab with an evaporation spectrum above the first threshold (no discrete levels); absorption reactions simply remove the neutron (no secondary particles, no fission); minor isotopes lumped into the main one; concrete as NIST "Concrete, Ordinary"; carbon as natural C. Energy range 10⁻¹¹–20 MeV; the spectrum is shown per unit lethargy on a log energy axis.

### Thermal scattering law S(α,β) for hydrogen in water

For **water**, a checkbox switches on the ENDF/B-VII.1 thermal scattering law of H bound in H₂O (`tsl_0001_H(H2O)`, MF7/MT4, 293.6 K, symmetric S(α,β), LAT = 1). Below 4 eV the elastic cross section of hydrogen is replaced by the thermal (inelastic incoherent) cross section computed from the table, and the outgoing energy and angle are sampled from it (up- and down-scattering, so neutrons thermalise to a Maxwellian near kT instead of being stopped at 0.0253 eV). The table is lazy-loaded from `js/data_tsl.js` (about 36 kB: 138 β × 94 α points, thinned 2× from the evaluation, ln S quantised).

Method: dσ/dE′ = σ_b·A₀/(4E) · e^(−β/2) ∫ S(α,β) dα with σ_b = σ_free((A₀+1)/A₀)²; the table is integrated over the allowed α range for a grid of incident energies; the energy transfer is sampled from the tabulated β grid (with jitter inside each bin) and α is then sampled from S(α,β) on [α₋, α₊] to obtain the scattering cosine. Not included: coherent elastic (none for liquid water), the oxygen of water is still a free nucleus (no O-in-H₂O law in VII.1), no temperature other than 293.6 K, ENDF-style continuity at 4 eV is not enforced (small step in Σ).

Nuclides: H-1, C (nat), O-16, Na-23, Mg-24, Al-27, Si-28/29/30, K-39/41, Ca-40, Fe-54/56/57/58, W-182/183/184/186, Pb-204/206/207/208. Data thinned in (ln E, ln σ) to about 7 %.

## Example 2: neutron point source in a water sphere

The buttons at the top of the page switch between **Example 1 (slab)** and **Example 2 (point source in a sphere)**; `index.html#sphere` opens example 2 directly, which is convenient for a link to students. Example 2 is neutrons only (ENDF/B-VII.1, same data and physics as the neutron slab, including the optional S(α,β) law for H in water).

- **Geometry**: homogeneous sphere of radius R (default water, 30 cm) in vacuum, isotropic point source of monoenergetic neutrons at the centre (cos θ = 2ξ₁ − 1, φ = 2πξ₂). A neutron that leaves the sphere never returns. The 2D views are projections of the 3D paths on a plane; there is also the 3D view with orbit camera. The step-by-step window and the automatic mode work as in example 1.
- **Bottom plot, flux spectrum in the sphere**: volume-averaged scalar flux per unit lethargy, E·φ(E) in 1/cm² per source neutron, from the **track-length estimator**: every flight adds its length (clipped at the surface) to the energy bin of the neutron during that flight, and the sum over histories is divided by (V · N · Δln E). The ±1σ bars come from the history-to-history scatter of the bin scores (`fl` and `fl2` in the tally). Linear or log y axis. Dashed grey line: kT = 0.0253 eV; red line: source energy. In the plot per unit lethargy the slowing-down region is flat (1/E flux) and a thermal neutron population is a Maxwellian peak at about 2kT (E·φ ∝ E² e^(−E/kT)); with the S(α,β) law switched off the neutrons pile up at the 0.0253 eV floor instead.
- **Right-hand plot**: energy spectrum of the neutrons that leak out of the sphere (collided ones, per unit lethargy), with the uncollided fraction marked; the table gives the leaked/absorbed fractions, the analytic uncollided leakage e^(−Σt·R), the mean number of collisions per history and the mean flux in the sphere.
- **Checks** (`node test/test_sphere.js`): distance to the surface; in a nearly empty sphere every neutron flies a straight line of length R, so the mean flux must equal R/V = 3/(4πR²) and sit in the bin of E₀ (exact); leaked + absorbed = N; uncollided leakage = e^(−Σt·R); the track-length estimator reproduces the collision count; isotropy of the source (⟨μ⟩ = 0, ⟨μ²⟩ = 1/3); generator = fast loop.

## Tests

`node test/test_neutron.js [path/to/data_n.js]` checks the neutron engine (uncollided = exp(−Σt·T), elastic kinematics, thermal floor, generator = fast loop, ...). Without an argument it uses the small subset `test/data_n_mini.js`; pass `js/data_n.js` to also test every material.

`node test/test_sphere.js [path/to/data_n.js]` checks the sphere example (see above).

`node test/test_tsl.js [js/data_tsl.js]` checks the thermal scattering code (cross section against the analytic free-gas formula with a synthetic table, detailed balance of the scattering chain, kinematics, conservation, generator = fast loop); pass the real file to run the checks that apply to it.

`node test/test_engine.js` checks the engine without a browser: table interpolation, uncollided transmission against exp(−Σt·T) for four material/energy cases, the Klein–Nishina sampler against numerical integration, and that the step-by-step generator and the fast loop give identical tallies for the same seed.

## Files

```
index.html          page
css/style.css       styles
js/data.js          XCOM tables (partial cross sections, cm²/g) and densities
js/data_n.js        ENDF/B-VII.1 neutron cross sections (lazy-loaded; see Neutron mode)
js/data_tsl.js      ENDF/B-VII.1 S(alpha,beta) of H in H2O (lazy-loaded when water + neutrons)
js/engine.js        Monte Carlo engine (no DOM; works in the browser and in Node)
js/plot.js          small canvas plotting helper
js/i18n.js          English / Spanish text
js/app.js           user interface
test/test_engine.js engine tests (Node)
test/test_neutron.js neutron tests (Node); test/data_n_mini.js tiny test data
test/test_tsl.js    S(alpha,beta) tests (Node)
test/test_sphere.js sphere example tests (Node)
```

## Data source

M. J. Berger, J. H. Hubbell, S. M. Seltzer, J. Chang, J. S. Coursey, R. Sukumar, D. S. Zucker and K. Olsen,
*XCOM: Photon Cross Sections Database*, NIST Standard Reference Database 8 (XGAM), https://physics.nist.gov/xcom (tables retrieved 2026-11-15).

M. B. Chadwick et al., *ENDF/B-VII.1 Nuclear Data for Science and Technology*, Nucl. Data Sheets 112 (2011) 2887.
Retrieved from NNDC Sigma, https://www.nndc.bnl.gov/sigma (2026-11-15).

Thermal scattering law: ENDF/B-VII.1 `tsl_0001_H(H2O)` (MacFarlane, LANL evaluation), taken from the IAEA-NDS mirror https://www-nds.iaea.org/public/download-endf/ENDF-B-VII.1/tsl/ (2026-11-15).
