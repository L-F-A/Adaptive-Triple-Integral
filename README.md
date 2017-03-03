A simple c++ class that implements an adaptive normalized integrator in 3d over a cube for functions that are complex. For function purely real, more efficient to save a slightly modified version where one changes "dcomplex" by "double", but everything else stay the same. 

More precisely it tries to integrate directly

<a href="https://www.codecogs.com/eqnedit.php?latex=I=\frac{1}{(2h)^3}\int_{-h}^{h}\int_{-h}^{h}\int_{-h}^{h}dx_1dx_2dx_3&space;f(x_1,x_2,x_3)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?I=\frac{1}{(2h)^3}\int_{-h}^{h}\int_{-h}^{h}\int_{-h}^{h}dx_1dx_2dx_3&space;f(x_1,x_2,x_3)" title="I=\frac{1}{(2h)^3}\int_{-h}^{h}\int_{-h}^{h}\int_{-h}^{h}dx_1dx_2dx_3 f(x_1,x_2,x_3)" /></a>

To do so, it uses the adaptive 5th order Gaussian quadrature as presented in Appendix A of Louis-François Arsenault, Patrick Sémon, and A.-M. S. Tremblay, Phys. Rev. B 86, 085133 (2012) (slightly different version available on ArXiv: https://arxiv.org/abs/1202.5814)
