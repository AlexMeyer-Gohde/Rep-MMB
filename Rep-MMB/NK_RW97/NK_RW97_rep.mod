//**************************************************************************
// A New Comparative Approach to Macroeconomic Modeling and Policy Analysis
//
// Volker Wieland, Tobias Cwik, Gernot J. Mueller, Sebastian Schmidt and
// Maik Wolters
//
// Working Paper, 2009
//**************************************************************************

// Model: NK_RW97

// Further references:
// Rotemberg, J., and M. Woodford. 1997. "An Optimization-Based Econometric Framework for the Evaluation of Monetary Policy."
// NBER Macroeconomics Annual 12, pp. 297-346.

// Last edited: 10/08/26 by S. Schmidt

// See Woodford(2003) p.246 for the model equations.


var pi y ynat rnat i x u g;



varexo u_ g_;




parameters

 beta sigma alpha theta omega kappa rhou rhog stdinflation_ stdfiscal_ phipi phix;

beta = 1/(1+0.035/4);  // 0.9913
sigma= 6.25;
alpha= 0.66;
theta= 7.66;
omega= 0.47;
kappa= (((1-alpha)*(1-alpha*beta))/alpha)*(((1/sigma)+omega)/(1+omega*theta));
rhou=0;
stdinflation_=0.154;
rhog= 0.8;
stdfiscal_=1.524;
phipi = 1.1;
phix = 1;

model(linear);



// Original Model Code:

pi   =  beta * pi(+1)+ kappa*x+ u;
u=rhou*u(-1)+u_;
x  =  x(+1) - sigma *( i - pi(+1) - rnat) ;
rnat = sigma^(-1)*((g-ynat)- (g(+1)-ynat(+1)));
ynat = sigma^(-1)*g /(sigma^(-1)+omega);
x = y-ynat;
g = rhog*g(-1) + g_;
i=phipi*pi + phix*x;
end;


shocks;
var g_= 1.524^2;
var u_=0.154^2;
end;
stoch_simul (AR=100,IRF=0, noprint,nograph);
//stoch_simul (irf = 0, ar=100, noprint);
