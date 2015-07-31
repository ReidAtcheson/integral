# integral
Black box integration routines for C++


I often run into the need of wanting an integration routine
which simply gives me an answer, no messing around with
quadrature parameters. While the latter is more useful
in numerical analysis because the numerical properties of
a method can then be analyzed as functions of integration
parameters, sometimes I just want a number and a reasonable
confidence level of its correctness.

## Some quesitons of interest:

The error estimate
==========


The error estimate on the [wikipedia page](https://en.wikipedia.org/wiki/Gauss%E2%80%93Kronrod_quadrature_formula)
is a little odd looking. Since I don't have a good reason yet to believe it, I am resorting to a simple
approximation to relative error. Is there a good estimate to use?


Bisection strategy
==================
A "dorfler" marking style strategy may be more efficient than simply bisecting any interval which
locally fails the error check.

