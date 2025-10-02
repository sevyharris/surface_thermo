#!/usr/bin/env python
# encoding: utf-8

name = "Fe_bcc100_gemnet"
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
            NASAPolynomial(coeffs=[2.02633,0.0168053,-2.0116e-05,1.40864e-08,-3.97032e-12,-27815.6,-9.21021], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.81632,-0.00639004,1.11462e-05,-5.73993e-09,9.99219e-13,-29466,-43.1831], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "bcc100",
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
            NASAPolynomial(coeffs=[0.969165,0.00854348,-1.50779e-05,1.23998e-08,-3.89789e-12,-44284.5,-2.93426], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.9324,-0.00024322,4.59286e-07,-2.60907e-10,4.9083e-14,-44657.7,-12.2442], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "bcc100",
)

entry(
    index = 3,
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
            NASAPolynomial(coeffs=[3.04494,0.00338008,-5.08056e-06,5.07474e-09,-1.96649e-12,-19150,-7.06769], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.44395,-0.00157361,2.8433e-06,-1.53988e-09,2.78962e-13,-19508.9,-14.0964], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "bcc100",
)

entry(
    index = 4,
    label = "NH_ads",
    molecule = 
"""
1 N u0 p1 c0 {2,S} {3,D}
2 H u0 p0 c0 {1,S}
3 X u0 p0 c0 {1,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.19088,0.0221757,-4.01666e-05,3.42131e-08,-1.09888e-11,-25296.4,-1.05806], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.03047,-0.00183407,3.20273e-06,-1.64088e-09,2.84277e-13,-26260.4,-25.6677], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "bcc100",
)

entry(
    index = 5,
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
            NASAPolynomial(coeffs=[-0.0929436,0.0261281,-4.40101e-05,3.65199e-08,-1.15449e-11,-28315.2,-1.69488], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.98405,-0.00385498,6.71083e-06,-3.43445e-09,5.94635e-13,-29746.9,-35.6752], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "bcc100",
)

entry(
    index = 6,
    label = "N_ads",
    molecule = 
"""
1 N u0 p1 c0 {2,T}
2 X u0 p0 c0 {1,T}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.953106,0.00916984,-1.68866e-05,1.43101e-08,-4.59806e-12,-26051.5,-5.12958], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.94349,-0.000190406,3.63564e-07,-2.07345e-10,3.9096e-14,-26409.9,-14.4637], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "bcc100",
)

entry(
    index = 7,
    label = "H_ads",
    molecule = 
"""
1 H u0 p0 c0 {2,S}
2 X u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.508452,0.00788838,-6.78363e-06,1.88813e-09,1.77443e-13,-9436.54,1.883], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.65759,-0.00131485,2.4254e-06,-1.35972e-09,2.5346e-13,-10273.9,-14.3249], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "bcc100",
)

entry(
    index = 8,
    label = "C_ads",
    molecule = 
"""
1 C u0 p0 c0 {2,Q}
2 X u0 p0 c0 {1,Q}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.656339,0.0104766,-1.92657e-05,1.63103e-08,-5.2371e-12,-16296.1,-4.11367], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.93484,-0.00022018,4.20212e-07,-2.39613e-10,4.51759e-14,-16707.1,-14.8031], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "bcc100",
)

entry(
    index = 9,
    label = "H2_ads",
    molecule = 
"""
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
3 X u0 p0 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.61478,0.0019549,-4.02158e-06,3.68574e-09,-1.15869e-12,-10709.4,-16.3311], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.07593,-0.000616349,8.81771e-07,-2.9619e-10,2.99846e-14,-10761.6,-18.3669], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "bcc100",
)

entry(
    index = 10,
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
            NASAPolynomial(coeffs=[3.94233,0.00840669,-1.22498e-05,9.94258e-09,-3.06725e-12,-45886.2,-15.9978], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.96842,-0.00318047,5.39156e-06,-2.64754e-09,4.42562e-13,-46553.6,-30.8156], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "bcc100",
)

entry(
    index = 11,
    label = "OH_ads",
    molecule = 
"""
1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 X u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.51373,0.0154215,-2.86085e-05,2.47021e-08,-7.96829e-12,-52098.2,-6.9479], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.05307,-0.00140532,2.37873e-06,-1.15764e-09,1.91796e-13,-52720.3,-23.4836], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "bcc100",
)

