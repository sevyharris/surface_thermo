#!/usr/bin/env python
# encoding: utf-8

name = "Pt_fcc100_gemnet"
shortDesc = ""
longDesc = """

"""
entry(
    index = 0,
    label = "X",
    molecule = 
"""
1 X u0 p0 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0,0,0,0,0,0,0], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[0,0,0,0,0,0,0], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
)

entry(
    index = 1,
    label = "NH_ads",
    molecule = 
"""
1 N u0 p1 c0 {2,S} {3,D}
2 H u0 p0 c0 {1,S}
3 X u0 p0 c0 {1,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.522095,0.0230275,-4.1015e-05,3.45378e-08,-1.10063e-11,1882.64,0.702895], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.00965,-0.00193793,3.39605e-06,-1.75081e-09,3.0497e-13,839.832,-25.4807], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc100",
)

entry(
    index = 2,
    label = "O_ads",
    molecule = 
"""
1 O u0 p2 c0 {2,D}
2 X u0 p0 c0 {1,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.706915,0.0100544,-1.82511e-05,1.53113e-08,-4.88378e-12,-17082.3,-3.65367], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.93245,-0.000233882,4.44579e-07,-2.53166e-10,4.76949e-14,-17490.8,-14.1315], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc100",
)

entry(
    index = 3,
    label = "H_ads",
    molecule = 
"""
1 H u0 p0 c0 {2,S}
2 X u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.968124,0.00972652,-9.98146e-06,4.51738e-09,-6.53052e-13,-4605.67,3.68566], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.63591,-0.00139077,2.56746e-06,-1.43976e-09,2.68426e-13,-5529.04,-14.612], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc100",
)

entry(
    index = 4,
    label = "H2_ads",
    molecule = 
"""
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
3 X u0 p0 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.46517,0.00265521,-5.3514e-06,4.83897e-09,-1.53528e-12,-1672.76,-17.3649], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.07305,-0.000633302,9.14036e-07,-3.14515e-10,3.34021e-14,-1750.23,-20.0826], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc100",
)

entry(
    index = 5,
    label = "OH_ads",
    molecule = 
"""
1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 X u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.469717,0.0190287,-3.40462e-05,2.86609e-08,-9.09792e-12,-26339.2,-3.36602], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.00598,-0.00164051,2.81911e-06,-1.40845e-09,2.39008e-13,-27179.4,-24.776], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc100",
)

entry(
    index = 6,
    label = "NH3_ads",
    molecule = 
"""
1 N u0 p1 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
5 X u0 p0 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.69902,0.0163997,-1.80086e-05,1.16626e-08,-3.06679e-12,-14295,-7.00396], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.67066,-0.00661833,1.15383e-05,-5.93809e-09,1.03349e-12,-16026.4,-42.0759], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc100",
)

entry(
    index = 7,
    label = "C_ads",
    molecule = 
"""
1 C u0 p0 c0 {2,Q}
2 X u0 p0 c0 {1,Q}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.190483,0.0133617,-2.34968e-05,1.92688e-08,-6.04361e-12,2173.45,-0.0312176], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.89295,-0.000386881,7.30147e-07,-4.14715e-10,7.80127e-14,1585.18,-14.665], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc100",
)

entry(
    index = 8,
    label = "N_ads",
    molecule = 
"""
1 N u0 p1 c0 {2,T}
2 X u0 p0 c0 {1,T}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.0447309,0.0122877,-2.15056e-05,1.75782e-08,-5.50052e-12,9336.07,-0.814305], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.8985,-0.000368407,6.94634e-07,-3.94381e-10,7.41681e-14,8788.47,-14.3743], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc100",
)

entry(
    index = 9,
    label = "N2_ads",
    molecule = 
"""
1 N u0 p1 c0 {2,D} {3,S}
2 N u0 p1 c0 {1,D} {4,S}
3 X u0 p0 c0 {1,S}
4 X u0 p0 c0 {2,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[3.09239,0.00335031,-5.45516e-06,5.62256e-09,-2.18272e-12,-1510.27,-7.12614], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.41874,-0.00158778,2.85897e-06,-1.54044e-09,2.77884e-13,-1840.32,-13.7347], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc100",
)

entry(
    index = 10,
    label = "NH2_ads",
    molecule = 
"""
1 N u0 p1 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 X u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.71998,0.0311438,-5.08684e-05,4.10223e-08,-1.27021e-11,-10278.7,5.25427], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.85773,-0.00423287,7.40259e-06,-3.81725e-09,6.65435e-13,-12055.4,-36.1526], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc100",
)

entry(
    index = 11,
    label = "H2O_ads",
    molecule = 
"""
1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 X u0 p0 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.67218,0.000469858,2.35789e-06,-2.41112e-09,8.91355e-13,-30327,-11.2457], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.97814,-0.00300636,5.0503e-06,-2.44566e-09,4.03741e-13,-30686.1,-18.0035], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc100",
)

