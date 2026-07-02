# Multi-rate (isoconversional) TGA in one sitting

A short reference on multi-rate thermogravimetry — what it is, why it
works, and where to read more. Companion to `Methods.tex`, kept
separate because it is background on a technique the `tga-analyze`
script does not itself implement.

## The basic idea

Run the same material (ideally from the same homogenized batch) under
identical atmosphere and sample geometry, but at several different
heating rates $\beta$ — typically a fourfold spread is enough, e.g.,
$\beta \in \{2.5, 5, 10, 20, 40\}\,\mathrm{K\,min^{-1}}$. Each ramp
produces a different $W(T)$ curve. Then for any chosen fractional
conversion $\alpha$ (say $0.1, 0.2, \ldots, 0.9$), read off the
temperature $T_\alpha(\beta)$ at which that conversion occurs in each
ramp. The function $T_\alpha(\beta)$ is the kinetic fingerprint of the
event at extent $\alpha$.

## Why it works

A solid-state reaction obeys

$$
\frac{\mathrm{d}\alpha}{\mathrm{d}t} = A\,e^{-E_\alpha / RT}\, f(\alpha),
$$

with $A$ a pre-exponential, $E_\alpha$ an (in general
$\alpha$-dependent) activation energy, and $f(\alpha)$ a reaction-model
function. On a linear ramp $T = T_0 + \beta t$, so

$$
\frac{\mathrm{d}\alpha}{\mathrm{d}T} = \frac{A}{\beta}\,e^{-E_\alpha/RT}\, f(\alpha).
$$

The crucial consequence: at higher $\beta$, the reaction has less time
per Kelvin to advance, so any given $\alpha$ is reached at a higher
$T$. The temperature shift with $\beta$, at fixed $\alpha$, depends
only on the activation energy of the reaction at that $\alpha$ — not
on the (often poorly known) functional form $f(\alpha)$. That is the
"model-free" promise of isoconversional kinetics.

## Why it separates overlapping events

Two chemical processes with different activation energies shift along
the temperature axis by different amounts when you change $\beta$. If
you watch two peaks that overlap at $10\,\mathrm{K\,min^{-1}}$, then
ramp at $2.5$ and $40\,\mathrm{K\,min^{-1}}$ as well, peaks with
different $E_a$ slide apart at different speeds. Even if they always
overlap at any single rate, their $T_\alpha(\beta)$ traces fan out
distinguishably across rates. That fan-out is information you cannot
extract from a single ramp at any quality.

## The standard analyses

Three families, all model-free, all of the form "linear-regress
$T_\alpha$ on $\beta$ in some clever way":

- **Friedman (1964)** — differential, takes logs of
  $\mathrm{d}\alpha/\mathrm{d}t$ directly. Cleanest but noise-sensitive
  because of the derivative.
- **Ozawa–Flynn–Wall (OFW, 1965/66)** — integral, uses Doyle's
  approximation. Simple, but the approximation systematically biases
  $E_\alpha$ at extreme conversions.
- **Vyazovkin advanced (1997, 2001)** — nonlinear integral, minimizes
  a model-free objective. This is the form the ICTAC committee
  recommends.

## The reference to start with

Vyazovkin, S., Burnham, A. K., Criado, J. M., Pérez-Maqueda, L. A.,
Popescu, C., Sbirrazzuoli, N. (2011). "ICTAC Kinetics Committee
recommendations for performing kinetic computations on thermal
analysis data." *Thermochimica Acta* **520**, 1–19.

This is the document of record — it specifies what to do, what not to
do, what to report, and which methods are obsolete. Forty pages,
designed to be readable. The follow-up Vyazovkin et al. (2014)
*Thermochim. Acta* **590**, 1–23, covers multi-step kinetics
specifically and is where the "isoconversional separation of
overlapping events" case is made explicitly. The 2020 update,
*Thermochim. Acta* **689**, 178597, is the latest revision.

## Book-length treatment

Vyazovkin, S., *Isoconversional Kinetics of Thermally Stimulated
Processes* (Springer, 2015). Single-author, sits in the middle ground
between a textbook and a research monograph. The first three chapters
are a clean introduction to the framework; chapters 5–6 are about
overlapping reactions specifically.

## Historical / textbook background

Brown, M. E., *Introduction to Thermal Analysis: Techniques and
Applications*, 2nd ed. (Springer/Kluwer, 2001). Chapter 10 ("Reaction
Kinetics from Thermal Analysis") is the gentlest correct introduction
to all of this. Brown also has a series of polemical papers in JTAC in
the 1990s/2000s about the misuse of single-rate Arrhenius fits in
solid-state kinetics, worth reading if you ever feel the urge to fit
kinetics to a single ramp.

## For cement specifically

Most cement-hydration TGA work *does not* do isoconversional analysis
— it reports masses and characteristic temperatures and stops there.
There is a small literature on isoconversional analysis of portlandite
dehydration and of carbonate decomposition in carbonated cement, but
the field has not really adopted it as standard practice. If you
wanted a methodological angle distinct from a "method-agreement as
diagnostic" paper, "isoconversional kinetics on PC + zeolite ternary
carbonation" would be unusual and citable — but it requires the
multi-rate experimental campaign, which is several days of instrument
time per sample. Worth it only if the kinetic separation is the
question, not just the mass balance.

## One pitfall to know about

Multi-rate work is extremely sensitive to *thermal lag* and
*self-cooling* — the temperature the sample experiences is not exactly
the programmed temperature, and the offset grows with $\beta$ and with
sample mass. The ICTAC recommendations spend several pages on this.
Practical rule: keep sample mass small ($\lesssim 5\,\mathrm{mg}$ for
$\beta \le 20\,\mathrm{K\,min^{-1}}$) and run a temperature-calibration
sample (indium melt, etc.) at every $\beta$ you use, not just at one.
