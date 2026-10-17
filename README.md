# GRS.jl_paper

This repository contains scripts and simulation schedule specifications that can be used to reproduce the results and figures presented in:

> "A generative model to simulate the dynamics of gene regulatory systems at cellular scale"

Note that the individual artifacts were generated with different software environments, as recorded in their respective subdirectories. For convenience, the following table lists the concrete versions of [Julia](https://github.com/julialang/julia) and [GeneRegulatorySystems.jl](https://github.com/drostlab/GeneRegulatorySystems.jl) that will reproduce the figures.

| Figure | Julia | GeneRegulatorySystems.jl |
|------------------------|----------------|--------------------------------|
| 02 Fano | `v1.10.4` | `v0.3.1-27-g274ef70` |
| 03 bulk / ecoli | `v1.12.6` | `v0.5.0` |
| 03 bulk / human | `v1.12.6` | `v0.5.0` |
| 04 scheduling | `v1.10.4` | `v0.3.1-28-gfba2261` |
| 05 AC-DC | `v1.10.6` | `v0.3.2-5-g46b20fe` |
| 06 differentiation / proportion | `v1.11.8` | `v0.3.2-6-gdaf682b` |
| 06 differentiation / timing | `v1.11.5` | `v0.3.2-6-gdaf682b` |
| 06 differentiation / flag | `v1.11.6` | `v0.3.2-6-gdaf682b` |
| 06 differentiation / grid | `v1.11.8` | `v0.3.2-6-gdaf682b` |
| 07 random differentiation | `v1.11.7` | `v0.3.2-6-gdaf682b` |
| 08 benchmark | `v1.12.4` | [nifets/GeneRegulatorySystems.jl](https://github.com/nifets/GeneRegulatorySystems.jl)\#`d8a5a364` |

(The benchmark was based on a development branch containing optimizations that have not been merged into the upstream GeneRegulatorySystems.jl repository yet.)
