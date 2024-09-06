//**************************************************************************
// A New Comparative Approach to Macroeconomic Modeling and Policy Analysis
//
// Volker Wieland, Tobias Cwik, Gernot J. Mueller, Sebastian Schmidt and
// Maik Wolters
//
// Working Paper, 2009
//**************************************************************************

// Model: NK_LWW03

// Further references:
// Levin, A., V. Wieland, and J. Williams. 2003. "The Performance of Forecast-Based Monetary Policy Rules under Model Uncertainty."
// American Economic Review 93(3), pp. 622-645.

// Last edited: 10/09/07 by S. Schmidt


var ygap pdot rff rstar drff pdotsh;



varexo rstar_ pdotsh_;



parameters

discountt sigma phi wtrl rhorstar rhopish;


discountt = 0.990;
sigma	 = 1/(0.157);
phi	 = 0.024;
wtrl	 = 0.975;
rhorstar = 0.35;
rhopish  = 0;



model(linear);



// Original Model Code:

ygap  =  ygap(+1) - 0.25*sigma *( rff - pdot(+1) -rstar);
pdot  =  discountt*pdot(+1) + 4*phi*ygap + pdotsh;
drff  =  rff - rff(-1);
rstar =  rhorstar*rstar(-1)+ rstar_;
pdotsh = rhopish*pdotsh(-1) + pdotsh_;
rff = 1.5*pdot;


end;


shocks;
var pdotsh_=(1-rhopish^2)*2.25^2;
//var interest_=0;   //interest rate shock is added
var rstar_=(1-rhorstar^2)*3.72^2;
end;
stoch_simul (AR=100,IRF=0, noprint,nograph);
//stoch_simul (irf = 0, ar=100, noprint);
