# script to check that NASA polynomial for H2O matches this one:

import rmgpy.thermo

thermo_H2O_katrin = rmgpy.thermo.NASA(
    polynomials = [
        rmgpy.thermo.NASAPolynomial(coeffs=[2.53777142E+00, 9.45372403E-03, -1.41325692E-05, 1.16730947E-08, -3.67657436E-12,
                                    -3.28227904E+04, -5.36547962E+00], Tmin=(300.0,'K'), Tmax=(1000.0, 'K')),
        rmgpy.thermo.NASAPolynomial(coeffs=[5.84789892E+00, -3.31527411E-03, 5.62019396E-06, -2.75865167E-09, 4.61279520E-13,
                                    -3.35523074E+04, -2.15622951E+01], Tmin=(1000.0,'K'), Tmax=(2000.0, 'K')),
    ],
    Tmin = (300.0, 'K'),
    Tmax = (2000.0, 'K'),
)


coeff1 = [2.53777142E+00, 9.45372403E-03,
                 -1.41325692E-05, 1.16730947E-08, -3.67657436E-12,
                 -3.28227904E+04, -5.36547962E+00]
coeff2 = [5.84789892E+00, -3.31527411E-03,
                 5.62019396E-06, -2.75865167E-09, 4.61279520E-13,
                 -3.35523074E+04, -2.15622951E+01]
other_thermo_H2O = rmgpy.thermo.NASA(
    polynomials = [
        rmgpy.thermo.NASAPolynomial(coeffs=coeff1, Tmin=(300.0,'K'), Tmax=(1000.0, 'K')),
        rmgpy.thermo.NASAPolynomial(coeffs=coeff2, Tmin=(1000.0,'K'), Tmax=(2000.0, 'K')),
    ],
    Tmin = (300.0, 'K'),
    Tmax = (2000.0, 'K'),
)

assert thermo_H2O_katrin.is_similar_to(other_thermo_H2O)

