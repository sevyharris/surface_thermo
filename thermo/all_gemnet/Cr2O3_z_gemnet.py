#!/usr/bin/env python
# encoding: utf-8

name = "Cr2O3_z_gemnet"
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
    label = "N_ads",
    molecule = 
"""
1 N u0 p1 c0 {2,T}
2 X u0 p0 c0 {1,T}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.757413,0.00614765,-7.32242e-06,4.15175e-09,-9.06723e-13,27062.6,-3.60006], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.8247,-0.000677818,1.25379e-06,-7.04592e-10,1.31582e-13,26555.7,-13.9838], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Cr2O3",
    facet = "z",
)

entry(
    index = 2,
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
            NASAPolynomial(coeffs=[2.80481,0.0137511,-2.17107e-05,1.76443e-08,-5.44997e-12,-7835.02,-10.387], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.0315,-0.00331682,5.67404e-06,-2.82778e-09,4.78661e-13,-8723.68,-30.8665], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Cr2O3",
    facet = "z",
)

entry(
    index = 3,
    label = "O_ads",
    molecule = 
"""
1 O u0 p2 c0 {2,D}
2 X u0 p0 c0 {1,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.814215,0.00889449,-1.53337e-05,1.2398e-08,-3.84852e-12,-20885.6,-3.60159], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.92024,-0.000293262,5.51599e-07,-3.12857e-10,5.87997e-14,-21296.6,-13.6442], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Cr2O3",
    facet = "z",
)

entry(
    index = 4,
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
            NASAPolynomial(coeffs=[4.05626,0.00335133,-2.73546e-06,1.63373e-09,-3.29277e-13,-31991.5,-9.92775], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.97392,-0.00293873,4.92655e-06,-2.37618e-09,3.91015e-13,-32459.7,-19.5589], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Cr2O3",
    facet = "z",
)

entry(
    index = 5,
    label = "N2_ads",
    molecule = 
"""
1 N u0 p1 c0 {2,T}
2 N u0 p1 c0 {1,T}
3 X u0 p0 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.1149,-0.000917745,1.70469e-06,-1.32223e-10,-3.71218e-13,-7334.14,-8.6787], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.39272,-0.00150971,2.68499e-06,-1.42268e-09,2.53094e-13,-7444.99,-10.2219], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Cr2O3",
    facet = "z",
)

entry(
    index = 6,
    label = "H2_ads",
    molecule = 
"""
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
3 X u0 p0 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.33761,0.00321492,-6.3559e-06,5.6654e-09,-1.79559e-12,-574.053,-15.982], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.06654,-0.000623425,8.96613e-07,-3.05097e-10,3.18001e-14,-674.174,-19.2719], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Cr2O3",
    facet = "z",
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
            NASAPolynomial(coeffs=[1.86449,0.00521339,-9.75608e-06,8.35988e-09,-2.70779e-12,55782.3,-8.41988], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.97099,-9.40357e-05,1.807e-07,-1.03269e-10,1.94942e-14,55587.6,-13.5851], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Cr2O3",
    facet = "z",
)

entry(
    index = 8,
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
            NASAPolynomial(coeffs=[2.00178,0.0120935,-1.13292e-05,6.49511e-09,-1.4864e-12,-10858.4,-3.34964], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.76254,-0.00622785,1.08114e-05,-5.5281e-09,9.56875e-13,-12321.5,-32.4955], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Cr2O3",
    facet = "z",
)

entry(
    index = 9,
    label = "OH_ads",
    molecule = 
"""
1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 X u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.81671,0.0138275,-2.53095e-05,2.15987e-08,-6.89118e-12,-33551.3,-7.28957], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.03622,-0.00129132,2.15564e-06,-1.02456e-09,1.66259e-13,-34122.1,-22.3662], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Cr2O3",
    facet = "z",
)

entry(
    index = 10,
    label = "H_ads",
    molecule = 
"""
1 H u0 p0 c0 {2,S}
2 X u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.958096,0.0143488,-2.27605e-05,1.7317e-08,-5.13751e-12,9055.89,3.07181], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.80291,-0.000746877,1.39264e-06,-7.863e-10,1.47326e-13,8260.53,-15.176], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Cr2O3",
    facet = "z",
)

entry(
    index = 11,
    label = "NH_ads",
    molecule = 
"""
1 N u0 p1 c0 {2,S} {3,D}
2 H u0 p0 c0 {1,S}
3 X u0 p0 c0 {1,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.71386,0.0144137,-2.66659e-05,2.30582e-08,-7.44676e-12,12679.5,-7.64077], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.0658,-0.00144569,2.45717e-06,-1.20523e-09,2.01026e-13,12085.8,-23.3214], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Cr2O3",
    facet = "z",
)

