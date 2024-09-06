
//**************************************************************************
// A New Comparative Approach to Macroeconomic Modeling and Policy Analysis
//
// Volker Wieland, Tobias Cwik, Gernot J. Mueller, Sebastian Schmidt and
// Maik Wolters
//
// Working Paper, 2009
//**************************************************************************

// Model: NK_RW06

// Further references:
// Ravenna, Federico and Walsh, Carl, E. (2006). Optimal monetary policy with the cost channel
// Journal of Monetary Economics 53, 199-216.

//The parameter values are calibrated as described in section 4; p.212-213.

// Last edited: 30/07/11 by M. Jancokova

var x pi R;



varexo u;



parameters

sigma eta beta omega kappa phipi phix;


sigma=1.5;
eta=1;
beta=0.99;
omega=0.75;
kappa=(1-omega)*(1-omega*beta)/omega;
phipi = 1.1;
phix = 1;



model(linear);

// Policy Rule
R = phipi*pi+phix*x;

//IS curve
x=x(+1)-(1/sigma)*(R-pi(+1))+u;

//Phillips curve
pi=beta*pi(+1)+kappa*(sigma+eta)*x+kappa*R;

end;



shocks;
var u =1;
//var interest_ =1;
end;
stoch_simul (AR=100,IRF=0, noprint,nograph);