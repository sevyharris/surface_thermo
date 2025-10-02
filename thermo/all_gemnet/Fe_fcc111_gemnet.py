#!/usr/bin/env python
# encoding: utf-8

name = "Fe_fcc111_gemnet"
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
            NASAPolynomial(coeffs=[4.42485,0.00623854,-8.29297e-06,6.59305e-09,-1.99125e-12,-34148.1,-17.8522], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.96733,-0.00310051,5.23623e-06,-2.5558e-09,4.24984e-13,-34726.8,-30.395], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "fcc111",
)

entry(
    index = 2,
    label = "H2_ads",
    molecule = 
"""
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
3 X u0 p0 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.76907,0.001243,-2.69662e-06,2.56408e-09,-7.97976e-13,-2075.51,-16.0013], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.08271,-0.000627099,9.00634e-07,-3.06441e-10,3.17467e-14,-2101.5,-17.3469], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "fcc111",
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
            NASAPolynomial(coeffs=[0.229424,0.012321,-2.25778e-05,1.90667e-08,-6.11093e-12,-42699.8,-2.41892], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.92182,-0.000266065,5.07211e-07,-2.89116e-10,5.44983e-14,-43187.8,-15.0622], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "fcc111",
)

entry(
    index = 4,
    label = "C_ads",
    molecule = 
"""
1 C u0 p0 c0 {2,Q}
2 X u0 p0 c0 {1,Q}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[0.87269,0.00943744,-1.72694e-05,1.4571e-08,-4.66733e-12,5222.58,-4.8425], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.93937,-0.000206896,3.94198e-07,-2.24647e-10,4.23399e-14,4847.17,-14.5516], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "fcc111",
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
            NASAPolynomial(coeffs=[1.22968,0.0171027,-3.21563e-05,2.79263e-08,-9.04986e-12,-38186.8,-7.16133], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.04415,-0.0013448,2.26514e-06,-1.09094e-09,1.79083e-13,-38842.8,-24.9087], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "fcc111",
)

entry(
    index = 6,
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
            NASAPolynomial(coeffs=[2.85494,0.0100167,-1.43782e-05,1.01449e-08,-2.84436e-12,-8764.23,-13.066], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.78863,-0.000798994,1.48235e-06,-8.3378e-10,1.55781e-13,-9432.28,-27.5362], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "fcc111",
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
            NASAPolynomial(coeffs=[-2.26757,0.0149491,-1.82793e-05,1.05629e-08,-2.33554e-12,-7094.95,9.43879], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.62352,-0.00147558,2.73081e-06,-1.53635e-09,2.87182e-13,-8276.81,-15.0506], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "fcc111",
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
            NASAPolynomial(coeffs=[2.00857,0.0158234,-1.76524e-05,1.17485e-08,-3.18011e-12,-15542.5,-8.89485], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.73261,-0.00652757,1.13827e-05,-5.85974e-09,1.01994e-12,-17207.4,-42.6901], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "fcc111",
)

entry(
    index = 9,
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
            NASAPolynomial(coeffs=[1.1018,0.0206524,-3.43167e-05,2.88469e-08,-9.25021e-12,-15064.4,-6.235], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[7.01888,-0.00393395,6.86713e-06,-3.53173e-09,6.13797e-13,-16294.4,-34.7875], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "fcc111",
)

entry(
    index = 10,
    label = "N_ads",
    molecule = 
"""
1 N u0 p1 c0 {2,T}
2 X u0 p0 c0 {1,T}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.150037,0.0136806,-2.46776e-05,2.06125e-08,-6.55405e-12,-14123.8,-0.749779], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.90435,-0.000334459,6.34663e-07,-3.61178e-10,6.8018e-14,-14689.1,-15.1541], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "fcc111",
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
            NASAPolynomial(coeffs=[-1.49384,0.027853,-5.03846e-05,4.26798e-08,-1.36569e-11,-18640.9,4.08743], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.98286,-0.00193469,3.39003e-06,-1.74403e-09,3.03318e-13,-19834.7,-26.4337], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Fe",
    facet = "fcc111",
)

