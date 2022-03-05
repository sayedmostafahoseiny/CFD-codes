In this Project we solve 1D-steady heat transfer equation in a pin fin with uniform cross section by finite difference method

The code includes below parts :
1-Definition of constant parameters
2-Assembling coefficient matrix (A) in a sparse form
3-Assembling right hand side matrix (B) in a sparse form
4-solving the linear algebraic sparse system of equations by Gmres method
5-Error Analysis (calculation of first norm of error)
6-validation (comparision with Analytical solution)
7-PostProcessing ( Heat flux - fin efficiency - fin performance coefficient )


constant parameters :

R = fin radius = 1cm
L = fin length = 20cm
k = fin conductivity = 200 w/(m.k)
h= ambient convective heat transfer coefficient = 25 w/(m^2.k)
T_inf = ambient temperature = 300k
T_base = temperature of fin base = 500k
