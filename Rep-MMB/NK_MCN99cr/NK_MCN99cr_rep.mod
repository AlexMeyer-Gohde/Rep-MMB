//**************************************************************************
// A New Comparative Approach to Macroeconomic Modeling and Policy Analysis
//
// Volker Wieland, Tobias Cwik, Gernot J. Mueller, Sebastian Schmidt and
// Maik Wolters
//
// Working Paper, 2009
//**************************************************************************

// Model: NK_MCN99cr

// Further references:
// McCallum, B., and E. Nelson. 1999. "Performance of Operational Policy Rules in an Estimated Semiclassical Structural Model."
// in: John B. Taylor (ed.), Monetary Policy Rules. Chicago: University of Chicago Press for NBER, pp. 15-45.

// Last edited: 10/09/07


var pi p y R v m i eta ytilde ybar;



varexo u_ e_ ey_ ev_;




parameters

sigm CssYss YssCss gam IssYss Rss rhov rhoeta gk stigma rhoybar bet thetac1 mu1 mu2 mu3;


sigm = .203;
CssYss = .81;
YssCss = 1/CssYss;
gam = 6.579;
IssYss = .19;
Rss = .014;
rhov = .3233;
rhoeta = .9346;
gk = 0;//.0073;
stigma = 0;// .0073;
rhoybar = 1;
bet = .99;
thetac1 = .3;

// Coeffs interest rule
mu1 = 1.5;
mu2 = 0;
mu3 = 0;



model(linear);






// original interest rule
R = mu1*pi + mu2*ytilde + mu3*R(-1); // table 1.1
//R = mu1*pi(-1) + mu2*ytilde(-1) + mu3*R(-1); // table 1.2



// Original Model Code:

// AD
pi     = p - p(-1);

y =     y(+1) - sigm*CssYss*(R - (pi(+1))) + CssYss*v;
m - p = (1/(sigm*gam))*YssCss*(y - IssYss*i) - (1/gam) * (1/Rss) * R + eta;

v =     rhov*v(-1) + ev_;
eta =   rhoeta*eta(-1) + u_;

i = gk + i(-1) + e_;
ybar = stigma + rhoybar*ybar(-1)+ ey_;


// AS
ytilde = y - ybar;
pi =    bet*pi(+1) + thetac1*ytilde; // Calvo-Rotemberg price setting

end;

shocks;
var ev_ = 0.0114^2;
var u_ =  0.0225^2;
var e_ = 0.025^2;
var ey_ = 0.007^2;

end;
stoch_simul (AR=100,IRF=0, noprint,nograph);
//stoch_simul (periods = 200, irf = 60);
