#!/usr/bin/env python
# encoding: utf-8

name = "Pt_fcc111_gemnet"
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
    label = "O_ads",
    molecule = 
"""
1 O u0 p2 c0 {2,D}
2 X u0 p0 c0 {1,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.07907,0.00873283,-1.62372e-05,1.38517e-08,-4.47218e-12,-15524.8,-5.87362], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.94938,-0.00016687,3.19784e-07,-1.82596e-10,3.44524e-14,-15857.1,-14.6206], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc111",
)

entry(
    index = 2,
    label = "N_ads",
    molecule = 
"""
1 N u0 p1 c0 {2,T}
2 X u0 p0 c0 {1,T}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.292773,0.0144046,-2.61058e-05,2.18752e-08,-6.97133e-12,6615,-0.201509], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.90242,-0.000338785,6.43724e-07,-3.66521e-10,6.90454e-14,6027.34,-15.2507], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc111",
)

entry(
    index = 3,
    label = "NH_ads",
    molecule = 
"""
1 N u0 p1 c0 {2,S} {3,D}
2 H u0 p0 c0 {1,S}
3 X u0 p0 c0 {1,D}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-2.16659,0.0301398,-5.37406e-05,4.50264e-08,-1.42966e-11,1003.81,7.05631], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[4.94688,-0.00206307,3.62774e-06,-1.87756e-09,3.2827e-13,-330.037,-26.5847], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc111",
)

entry(
    index = 4,
    label = "OH_ads",
    molecule = 
"""
1 O u0 p2 c0 {2,S} {3,S}
2 H u0 p0 c0 {1,S}
3 X u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[1.05554,0.0172148,-3.16303e-05,2.71313e-08,-8.7222e-12,-19305.5,-5.71676], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.04054,-0.00150157,2.56144e-06,-1.2631e-09,2.11784e-13,-20017.8,-24.3924], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc111",
)

entry(
    index = 5,
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
            NASAPolynomial(coeffs=[1.33423,0.0183674,-2.1765e-05,1.48261e-08,-4.06441e-12,-11954.4,-6.1116], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[8.68239,-0.00657569,1.1463e-05,-5.89767e-09,1.02622e-12,-13744.2,-42.9066], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc111",
)

entry(
    index = 6,
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
            NASAPolynomial(coeffs=[5.59624,0.000794857,1.79031e-06,-1.96846e-09,7.59656e-13,-30363.4,-16.177], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.97053,-0.00298969,5.01762e-06,-2.42575e-09,3.99901e-13,-30734.9,-23.257], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc111",
)

entry(
    index = 7,
    label = "H2_ads",
    molecule = 
"""
1 H u0 p0 c0 {2,S}
2 H u0 p0 c0 {1,S}
3 X u0 p0 c0
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[4.59243,0.00206258,-4.2361e-06,3.88277e-09,-1.22523e-12,-896.806,-16.786], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.07715,-0.000631491,9.09861e-07,-3.11931e-10,3.28531e-14,-952.759,-18.9295], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc111",
)

entry(
    index = 8,
    label = "H_ads",
    molecule = 
"""
1 H u0 p0 c0 {2,S}
2 X u0 p0 c0 {1,S}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-1.1431,0.014667,-2.28431e-05,1.71374e-08,-5.02943e-12,-4513.96,3.90795], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.78177,-0.00082871,1.54337e-06,-8.70747e-10,1.63061e-13,-5356.27,-15.197], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc111",
)

entry(
    index = 9,
    label = "C_ads",
    molecule = 
"""
1 C u0 p0 c0 {2,Q}
2 X u0 p0 c0 {1,Q}
""",
    thermo = NASA(
        polynomials = [
            NASAPolynomial(coeffs=[-0.820753,0.0159051,-2.78556e-05,2.27778e-08,-7.12933e-12,9387.25,2.26826], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[2.86943,-0.000473747,8.93393e-07,-5.07272e-10,9.54045e-14,8679.82,-15.2628], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc111",
)

entry(
    index = 10,
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
            NASAPolynomial(coeffs=[4.18307,0.00305244,-5.09837e-06,5.42276e-09,-2.13766e-12,1971.21,-16.5074], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[5.4138,-0.001579,2.83897e-06,-1.52637e-09,2.74851e-13,1665.21,-22.633], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc111",
)

entry(
    index = 11,
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
            NASAPolynomial(coeffs=[-1.13005,0.0297961,-4.96755e-05,4.0691e-08,-1.27412e-11,-2285.82,2.50155], Tmin=(298,'K'), Tmax=(1000,'K')),
            NASAPolynomial(coeffs=[6.92291,-0.00403448,7.04111e-06,-3.61782e-09,6.28643e-13,-3925.78,-36.2269], Tmin=(1000,'K'), Tmax=(2000,'K')),
        ],
        Tmin = (298,'K'),
        Tmax = (2000,'K'),
    ),
    shortDesc = """""",
    longDesc = 
"""

""",
    metal = "Pt",
    facet = "fcc111",
)

