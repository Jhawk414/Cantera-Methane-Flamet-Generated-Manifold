# Cantera-Methane-Flamet-Generated-Manifold
Purpose: Creates a Flamelet-Generated Manifold (FGM) using a progress-variable (PV) approach for a CH₄-air diffusion flame in Cantera/Python. Code is available for both the GRI-Mech 1.2 and 3.0 mechanisms. The code allows user inputs of mean mixture fraction 𝑍, normalized progress variable 𝐶_norm, and mean mixture fraction variance 𝑍′′ to calculated beta-PDF filtered combustion properties (T, ρ, Yₖ, for k species).

Status: Mostly-complete since May 2025. Current/future work is not being performed.

Recommendations:
1) Vary flame strain-rate instead of ambient pressures. This allows for simulation of near-extinguished flamelets to perfectly-combusted flamelets, which expands the condition of the manifold.
2) Filter manifold completely. Only the manifold values "looked-up" are β-PDF filtered based on the user-input of mixture fraction variance 𝑍′′. 

This was for my graduate-level combustion class (AE 774) completed in Spring 2025.

==== NOTES ====
1) A full *filtered* manifold (lookup table) of properties is not performed. Rather, an "unfiltered" manifold is built, from which properties are interpolated for from user-inputs.
2) The progress variable 𝐶 is normalized (𝐶_norm = (𝐶 − 𝐶_min)/(𝐶_max − 𝐶_min) to allow for the manifold to be built on a regular grid. This is done since each flamelet simulated has a varying maximum and minimum 𝐶.
