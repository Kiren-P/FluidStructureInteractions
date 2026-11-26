---
title: "FSI Group Assignment, Part I: Question 7"
subtitle: "A numerical model we would use, and one we would not use, to calculate the wave forces on the wall"
---

# What the model has to capture

## Wall and wave conditions

The wall is a vertical, impermeable wall standing on the bed in $h = 20$ m of water, with a width $W = 60$ m. The long-crested waves approach it at normal incidence. We used the sea state of Questions 3 to 6: JONSWAP with $H_{m0} = 6$ m, $T_p = 8.8$ s and $\gamma = 3.3$. At the peak period $k_p h = 1.23$ and $L_p = 102$ m (intermediate water), $H_{m0}/h = 0.30$ and the Ursell number is 7.8 (Question 4). In the 3 h record of Question 4 the largest wave is $H_{max} = 11.2$ m (linear) to 11.8 m (with second-order bound waves), with a highest crest of 8.1 m. For the wall load in Part II the current is neglected, so we leave it out here as well.

Table 1 lists the physical processes that set the force on this wall, and what each one demands from a numerical model.

Table 1. Processes that determine the wave force on the wall and the consequence for the model.

| Process | Evidence for this wall | What the model must do |
|---|---|---|
| Reflection, standing wave | A vertical impermeable wall reflects almost all energy ($K_r \approx 1$); at the wall $\eta$ is up to twice the incident value (Lecture 7, slide 3). | Resolve the phase of incident and reflected waves; absorb the reflected waves offshore. |
| Nonlinearity | Stokes range for $H_{m0}$, $H_{max}$ near the breaking limit (Le Méhauté); bound waves raise the highest crest from 5.7 to 8.1 m (Question 4b). | Nonlinear free surface, valid at $kh \approx 1$ and for the higher harmonics. |
| Breaking, impulsive loads | Highest waves at 92 to 96 % of Miche's limit; impacts on the deck at +10 m (Part II, Question 4). | Overturning surface, trapped air, a reliable force impulse. |
| End effects | $W/L_p = 0.59$: diffraction around the wall ends. | Three dimensions, or a correction for the finite width. |

**Nonlinearity.** In the Le Méhauté diagram ($h/gT_p^2 = 0.026$) $H_{m0}$ ($H/gT^2 = 0.008$) lies in the Stokes range and $H_{max}$ ($H/gT^2 = 0.015$) close to the breaking limit (Lecture 6, slides 13 and 23). With second-order bound waves the mean crest is 40 % higher than the mean trough (Question 4b). In front of the wall the standing wave doubles the surface elevation, so the waves at the wall are more nonlinear than the incident waves.

**Breaking and impacts.** Miche's limit at $L_p$ is $0.142 \tanh(k_p h) L_p = 12.2$ m (Lecture 6, slide 22), so the highest waves of the record (11.2 to 11.8 m) reach 92 to 96 % of the limit before they meet the wall. In the parameter map of PROVERBS (Oumeraci et al., 2001) a vertical wall on the bed with $H_s/h = 0.30 < 0.35$ falls in the "small wave" class, where quasi-standing, pulsating loads dominate. The highest individual waves ($H/h \approx 0.56$) exceed 0.35, however, and can break against the wall, giving short (milliseconds) impact pressures far above the pulsating ones. The deck of Part II, Question 4 adds confined impacts.

The force on the wall follows from the pressure integrated over the wetted height,
$$F_x(t) = W \int_{-h}^{\eta_w(t)} p(z,t)\,\mathrm{d}z,$$
with $\eta_w$ the surface elevation at the wall (or, in 3D, the integral over the wall face). The model must therefore deliver $p(z,t)$ at the wall face, not only wave heights.

# Model we would use: a two-phase RANS–VOF numerical wave flume

## The model

We would use a two-phase Reynolds-averaged Navier–Stokes model with a volume-of-fluid interface (RANS–VOF), for example OpenFOAM with the waves2Foam (Jacobsen et al., 2012) or olaFlow/IHFOAM (Higuera et al., 2013) wave toolboxes. Lecture 8 (slide 4) gives exactly this class of model for "breaking-wave impact (slamming) on a pile or wall": two-phase CFD (RANS) coupled to a potential-flow far field, together with lab tests. It solves the incompressible, Reynolds-averaged Navier–Stokes equations for water and air as one mixture, and tracks the surface with the volume fraction $\alpha$ (1 in water, 0 in air) (Lecture 8, slides 1 and 2).

**What it resolves.** The full nonlinear free surface, including overturning and jets, so the standing wave, the higher harmonics, the run-up against the wall, wave breaking at or in front of the wall, overtopping (Part II, Question 3) and the impact on a deck (Part II, Question 4). The pressure field follows from the momentum and continuity equations, so $p(z,t)$ on the wall, and with it $F_x(t)$, the overturning moment and the force impulse $\int F_x\,\mathrm{d}t$, come directly out of the model. Because the air is a phase in the model, the air pocket that a plunging breaker traps against the wall is represented.

**What it models.** Turbulence is not resolved but modelled with a closure such as $k$–$\omega$ SST. Standard closures over-produce turbulence under (non-breaking) waves and damp the waves artificially, so we would use a stabilised closure (Larsen and Fuhrman, 2018; Lecture 8, slide 1). The air is normally treated as incompressible.

## Set-up for this wall

- **2DV numerical flume** in a vertical slice normal to the wall through its middle: a generation and absorption zone, about two to three wavelengths of propagation (200 to 300 m), and the wall at the end. The air layer above the still water level must be high enough for the run-up of the standing wave (more than about 15 m).
- **Input from the Part I model.** At the inflow boundary we impose the incident surface elevation and velocities of Question 4 (linear plus second-order bound waves, same seed), or couple to a fully nonlinear potential-flow model in the far field (Lecture 8). The force record then belongs to the same waves as the rest of the assignment, and the linear standing-wave load of Part II, Question 2 is a direct check.
- **Active absorption at the generation boundary.** This is essential for a reflecting wall: without it the reflected waves are reflected again at the wavemaker and the energy in the flume builds up.
- **Resolution.** Near the free surface about 10 to 20 cells per wave height and at least 100 cells per wavelength, so roughly $\Delta z \approx 0.3$ m and $\Delta x \approx 0.5$ to 1 m, with a Courant number below about 0.5 ($\Delta t \approx 0.02$ s). This gives about $10^5$ to $3 \cdot 10^5$ cells and $5 \cdot 10^5$ time steps for the 3 h record: days on a few tens of cores. Impact pressures need local cells of a few centimetres at the wall and $\Delta t \sim 10^{-4}$ s.
- **Two-step strategy.** Run the whole 3 h record at standard resolution for the statistics of the quasi-standing force. Then select the largest wave groups from the Question 4 time series (or design them as focused groups) and rerun only these windows of a few minutes on the fine grid for the impulsive loads.
- **Verification and validation** (Lecture 8, slide 4). Grid and time-step convergence of the force for a regular wave; comparison with linear and second-order standing-wave theory (Sainflou, Goda) for non-breaking waves; validation of breaking cases with flume tests, for example in the facility of Part II, Question 8.

## Limits

- **Cost.** RANS–VOF is roughly 4000 times as expensive as a spectral model per simulated hour over the same coastline and is only practical on domains of 10 to 100 m (Lecture 8, slide 3). A 2DV flume is affordable; a 3D model of the 60 m wall with its ends needs about $10^7$ cells and is only feasible for a few selected events.
- **2DV ignores the wall ends.** With $W/L_p = 0.59$ diffraction around the ends lowers the standing wave over a large part of the wall, so a 2DV result is conservative. We would correct it with a linear diffraction (BEM) calculation of the reduction along the wall, or with 3D runs of the design events.
- **Impact pressures are uncertain.** The peak pressure depends on the grid, the time step, the shape of the breaker and the compressibility of the trapped air (Bredmose et al., 2009), and also scatters strongly in laboratory tests. The force impulse is much more robust than the peak, so we would report the quasi-static force and the impulse of the impact separately.
- **Numerical diffusion** reduces the wave height over long propagation distances; the domain must be kept short and the incident waves checked against the input.

For the non-breaking part of the load a cheaper wave-resolving model is also adequate: the non-hydrostatic model SWASH (Zijlema et al., 2011) with a few vertical layers handles reflection, dispersion at $kh \approx 1.2$ and the nonlinear harmonics, but it cannot represent an overturning surface or air, so it misses the impacts. Following Lecture 8, we would combine models rather than stretch one: potential flow or SWASH for the propagation and screening of many waves, RANS–VOF for the design events at the wall.

# Model we would not use: a spectral (phase-averaged) wave model

We would not use a spectral wave model such as SWAN (Booij et al., 1999) or WAVEWATCH III to calculate the forces on the wall. These models solve the wave action balance for the variance density $E(f,\theta)$ and are the right tool for the wave climate over a domain of 10 to 1000 km (Lecture 8, slides 2 to 4). For a wall force they fail for the following reasons:

- **No phase, so no force in time.** The output is $H_{m0}$, $T_p$, the mean direction and the spectrum. There is no $\eta(t)$ and no $p(z,t)$, so no force time series, no peak force and no impulse. A force can only follow afterwards from an empirical formula (Goda) applied to $H_{m0}$, and then the force comes from that formula, not from the model.
- **No standing wave.** Reflection is only a coefficient that sends energy back. Without phases the model cannot form the nodes and antinodes in front of the wall, so it does not give the doubled surface elevation and pressure at the wall face. Diffraction around the wall ends is only approximated.
- **No nonlinearity and breaking at wave level.** Triad interactions and depth-induced breaking are parameterised as source terms for the spectrum. The crest shape, the bound harmonics and the breaking of individual waves against the wall, which set the extreme and impulsive loads, are absent (Lecture 8, slide 3).
- **Wrong scale.** The grid spacing of a spectral model is tens to hundreds of metres; the whole wall fits in one or two cells.

A spectral model is still useful one step earlier: it can transform an offshore wave climate to the site and so provide $H_{m0}$, $T_p$ and the spectrum that we prescribed in Question 3. It should not be used for the force itself.

Table 2 places these two models among the other models of Lecture 8. Two more models are clearly unsuitable for this wall. Direct numerical simulation is out of reach: the Reynolds number of the orbital flow is about $u_{rms} h/\nu \approx 1.5 \cdot 20/10^{-6} = 3 \cdot 10^7$, and the cost grows as $Re^3$. The hydrostatic shallow-water equations have no frequency dispersion; at $k_p h = 1.23$ they give a celerity $\sqrt{gh} = 14.0$ m/s instead of 11.6 m/s (21 % too fast) and steepen the waves into bores, with a hydrostatic pressure on the wall.

Table 2. Suitability of the models of Lecture 8 for the force on this wall.

| Model | Standing wave | Nonlinear, breaking waves | Impact, air | Cost | Use for the wall force |
|---|---|---|---|---|---|
| Spectral (SWAN, WW3) | no (energy only) | parameterised | no | very low | **no**; only the sea state at the site |
| Shallow-water equations | yes | no dispersion at $kh = 1.2$ | no | low | no |
| Linear diffraction (BEM) | yes, in 3D | linear only | no | low | quasi-standing load, end effects |
| Non-hydrostatic (SWASH) | yes | yes, single-valued surface | no | moderate | quasi-standing load, screening |
| Two-phase RANS–VOF | yes | yes, overturning | yes (incompressible air) | high | **yes**, recommended |
| DNS | yes | yes | yes | impossible at $Re \sim 10^7$ | no |

# Conclusion

To calculate the forces on the wall we would use a two-phase RANS–VOF numerical wave flume, driven at its offshore boundary by the incident wave field of Question 4 and with active absorption of the reflected waves. It is the only model in the lecture overview that resolves the standing wave, the nonlinear crests, breaking at the wall and the trapped air together, and it delivers the pressure on the wall face directly. Its price is a high computational cost, a sensitivity of the impact peaks to resolution and air compressibility, and the need for validation with lab tests; we would therefore run the full record in 2DV, refine only the design events, and correct for the finite width of the wall. A spectral wave model such as SWAN we would not use for the force, because it carries no phase information and so cannot represent the standing wave, the pressure in time or the impacts; it belongs one step earlier in the chain, to define the sea state at the wall.

# References

- Course slides CIEM4000-26, Fluid-Structure Interaction: Lecture 6 (Nonlinear waves), Lecture 7 (Waves in varying water depth), Lecture 8 (Computational fluid dynamics).
- Booij, N., Ris, R.C. and Holthuijsen, L.H. (1999). A third-generation wave model for coastal regions: 1. Model description and validation. Journal of Geophysical Research 104(C4), 7649–7666.
- Bredmose, H., Peregrine, D.H. and Bullock, G.N. (2009). Violent breaking wave impacts. Part 2: modelling the effect of air. Journal of Fluid Mechanics 641, 389–430.
- Higuera, P., Lara, J.L. and Losada, I.J. (2013). Realistic wave generation and active wave absorption for Navier–Stokes models: application to OpenFOAM®. Coastal Engineering 71, 102–118.
- Jacobsen, N.G., Fuhrman, D.R. and Fredsøe, J. (2012). A wave generation toolbox for the open-source CFD library: OpenFoam®. International Journal for Numerical Methods in Fluids 70(9), 1073–1088.
- Larsen, B.E. and Fuhrman, D.R. (2018). On the over-production of turbulence beneath surface waves in Reynolds-averaged Navier–Stokes models. Journal of Fluid Mechanics 853, 419–460.
- Oumeraci, H., Kortenhaus, A., Allsop, W., de Groot, M., Crouch, R., Vrijling, H. and Voortman, H. (2001). Probabilistic Design Tools for Vertical Breakwaters (PROVERBS). Balkema, Lisse.
- Zijlema, M., Stelling, G. and Smit, P. (2011). SWASH: an operational public domain code for simulating wave fields and rapidly varied flows in coastal waters. Coastal Engineering 58(10), 992–1012.
