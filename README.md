# Shilling

Bucket full of solvers for nth order non-linear ODEs all implemented in MATLAB.

## List of Available Solvers

### Fixed Step Solvers

Name       | Description| Properties | Order | Stages   | Registers | CFL  | Reference
---        | ---        | ---        | ---   | ---      | ---       | ---  | ---
`euler` | Euler |  | 1 | 1 | 1 | 1.0 | [Euler (1768)](https://archive.org/details/institutionescal020326mbp)
`midpoint` | Midpoint |  | 2 | 2 | 2 |  | ?
`heun` | Heun |  | 2 | 2 | 2 |  | ?
`rkssp22` | 2-stage, 2nd order TVD Runge-Kutta Shu-Osher | Strong stability preserving | 2 | 2 | 1 | 1.0 | [Shu & Oscher (1988)](https://ntrs.nasa.gov/api/citations/19880014833/downloads/19880014833.pdf)
`rk3` | 3th order Runge-Kutta |  | 3 | 3 | 3 |  | ?
`rkssp33` | 3-stage, 3rd order TVD Runge-Kutta Shu-Osher | Strong stability preserving | 3 | 3 | 1 | 1.0 | [Shu & Oscher (1988)](https://ntrs.nasa.gov/api/citations/19880014833/downloads/19880014833.pdf)
`rkssp53` | 5-stage, 3rd order SSP Runge-Kutta Spiteri-Ruuth | Strong stability preserving | 3 | 5 | 2 | 2.65 | [Ruuth (2006)](https://www.ams.org/journals/mcom/2006-75-253/S0025-5718-05-01772-2/S0025-5718-05-01772-2.pdf)
`rk4` | Classic 4th order Runge-Kutta |  | 4 | 4 | 4 |  | [Kutta (1901)](https://archive.org/stream/zeitschriftfrma12runggoog#page/n449/mode/2up)
`rks4` | 4th order Runge-Kutta Shanks |  | 4 | 4 | 4 |  | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)
`rkr4` | 4th order Runge-Kutta Ralston |  | 4 | 4 | 4 |  | [Ralston (1962)](https://doi.org/10.1090%2FS0025-5718-1962-0150954-0)
`rkls44` | 4-stage, 4th order low storage non-TVD Runge-Kutta Jiang-Shu | Low storage | 4 | 4 | 2 |  | [Jiang and Shu (1988)](https://ntrs.nasa.gov/api/citations/19960007052/downloads/19960007052.pdf)
`rkls54` | 5-stage, 4th order low storage Runge-Kutta Carpenter-Kennedy | Low storage | 4 | 5 | 2 | 0.32 | [Carpenter & Kennedy (1994)](https://ntrs.nasa.gov/api/citations/19940028444/downloads/19940028444.pdf)
`rkssp54` | 5-stage, 4th order SSP Runge-Kutta Spiteri-Ruuth | Strong stability preserving | 4 | 5 | 4 | 1.51 | [Ruuth (2006)](https://www.ams.org/journals/mcom/2006-75-253/S0025-5718-05-01772-2/S0025-5718-05-01772-2.pdf)
`rks5` | 5th order Runge-Kutta Shanks |  | 5 | 5 | 5 |  | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)
`rk5` | 5th order Runge-Kutta |  | 5 | 6 | 6 |  | ?
`rkc5` | 5th order Runge-Kutta Cassity |  | 5 | 6 | 6 |  | [Cassity (1966)](https://epubs.siam.org/doi/10.1137/0703052)
`rkl5` | 5th order Runge-Kutta Lawson |  | 5 | 6 | 6 |  | [Lawson (1966)](https://epubs.siam.org/doi/abs/10.1137/0703051)
`rklk5a` | 5th order Runge-Kutta Luther-Konen 1 |  | 5 | 6 | 6 |  | [Luther & Konen (1965)](https://epubs.siam.org/doi/abs/10.1137/1007112)
`rklk5b` | 5th order Runge-Kutta Luther-Konen 2 |  | 5 | 6 | 6 |  | [Luther & Konen (1965)](https://epubs.siam.org/doi/abs/10.1137/1007112)
`rkb6` | 6th order Runge-Kutta Butcher |  | 6 | 7 | 7 |  | [Butcher (1963)](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/40DFE501CAB781C9AAE1439B6B8F481A/S1446788700023387a.pdf/div-class-title-on-runge-kutta-processes-of-high-order-div.pdf)
`rk7` | 7th order Runge-Kutta Shanks |  | 7 | 9 | 9 |  | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)
`rk8_10` | 10-stage, 8th order Runge-Kutta Shanks |  | 8 | 10 | 10 |  | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)
`rkcv8` | 11-stage, 8th order Runge-Kutta Cooper-Verner |  | 8 | 11 | 11 |  | [Cooper & Verner (1972)](https://epubs.siam.org/doi/abs/10.1137/0709037)
`rk8_12` | 12-stage, 8th order Runge-Kutta Shanks |  | 8 | 12 | 12 |  | [Shanks (1965)](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19650022581.pdf)
`rkz10` | 10th order Runge-Kutta Zhang |  | 10 | 16 | 16 |  | [Zhang (2019)](https://arxiv.org/abs/1911.00318)
`rko10` | 10th order Runge-Kutta Ono |  | 10 | 17 | 17 |  | [Ono (2003)](http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK10/RKcoeff10f_1.pdf)
`rkh10` | 10th order Runge-Kutta Hairer |  | 10 | 17 | 17 |  | [Hairer (1978)](https://www.researchgate.net/publication/31221486_A_Runge-Kutta_Method_of_Order_10)

