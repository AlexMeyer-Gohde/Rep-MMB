//**************************************************************************
// New Keynesian model with search and matching frictions in the labor market.
// Christoffel, Kuester (JME, 2008)
//**************************************************************************

// Model: NK_CK08.mod

// Further references:
// Christoffel, Kai, and Keith Kuester (2008),
// "Resuscitating the wage channel in models with unemployment
//  fluctuations," Journal of Monetary Economics 55, p. 865� 887
//
// Note: the model original model runs at a MONTHLY frequency (available upon request)
//       The file here presents a recalibration of the model to a QUARTERLY frequency
//
//
// Last edited: 2011/05/10 by K. Kuester


var ct deltaFt deltaWt Deltastart ht Jstart
    lambdat mct mt nt Pit Piannt  qt
    Rt st ut vt wstart wt xLt yt
    // shocks
    ebt emoneyt gt zt;

varexo inno_ebt inno_zt interest_ g_;




parameters

           bet epsilon habit sig vphi
           omega price_index
           xi eta gamma wage_index vtheta alp
           rho_eb rho_emoney rho_g rho_z
           gamma_R gamma_Pi gamma_y
           sig_innoeb sig_monpol sig_innog sig_innoz
           kappa Phi Phi_y Afactor
           cbar Deltabar deltaWbar gbar hbar Jbar
           mbar mcbar mrsbar nbar qbar
           sbar ubar vbar wbar ybar b_wh Pibar;

// ***** parameter values
bet     = 0.9939; // time-discount factor
epsilon = 11;                // own price elasticity of demand for differentiated good
habit   = 0.7;               // habit persistence
sig     = 1.5;               // CRRA
vphi    = 2;                 // inverse of Frisch elasticity

omega       = 0.5;  // Calvo price stickiness
price_index = 0;    // degree of price indexation

xi          = 0.5;  // weight on unemployment in matching function
eta         = 0.5;  // bargaining power of workers
gamma       = 0.5;  // Calvo wage stickiness
wage_index  = 0;    // wage indexation
vtheta      = 0.03*3; // rate of separation
alp         = 0.99; // labor elasticity of production

gamma_R     = 0.85;	 // interest rate smoothing
gamma_Pi    = 1.5;   // weight on inflation
gamma_y     = 0.5;   // weight on output

rho_eb      = 0.9;   // persistence of preference shock
rho_emoney  = 0;     // persistence of mon pol shock
rho_g       = 0.7912; // persistence of gov spending shock
rho_z       = 0.6712; // persistence of productivity shock

sig_innoeb =   0.3640;
sig_monpol =   0.25;
sig_innog  =   0.8716;
sig_innoz  =   0.6849;


// ********* calibration targets
ybar    = 1;                    // steady state output
hbar    = 1/3;                  // steady state hours worked
Phi_y   = 0.00863;              // fraction of output lost due to fixed costs.
ubar    = 0.1;                  // steady state unemployment rate
qbar    = 0.7;                  // probability of finding a worker in a given month
gbar    = 0.347026648444032;    // share of government spending in GDP (y=c+g)
b_wh    = 0.4;                  // replacement rate
Pibar   = 1;                    // zero inflation steady state



// ********* ss for the given targets
nbar        = 1-ubar;
mbar        = vtheta*nbar;
vbar        = mbar/qbar;
sbar        = mbar/ubar;
sigmam      = mbar*(ubar^xi*vbar^(1-xi))^(-1);
thetabar    = vbar/ubar;

mcbar       = (epsilon-1)/epsilon;
xLbar       = mcbar;
Rbar        = 1/bet;

zbar        = ybar/(nbar*hbar^(alp));

wbar        = xLbar*zbar*alp*hbar^(alp-1);
Phi         = Phi_y*ybar/nbar;

PsiLbar     = xLbar*zbar*hbar^(alp) - wbar*hbar - Phi;
PsiCbar     = (1-mcbar)*ybar;

Jbar		= 1/(1-bet*(1-vtheta))*PsiLbar;
deltaFbar   = 1/(1-bet*(1-vtheta)*gamma)*wbar*hbar;

b           = b_wh*wbar*hbar;

mrsbar      = (Jbar*eta/(1-bet*(1-vtheta)*gamma)*alp/(alp-1)*hbar*wbar - (1-eta)*deltaFbar/(1-bet*(1-vtheta-sbar))*(wbar*hbar-b))/ (Jbar*hbar*eta/(1-bet*(1-vtheta)*gamma)*(-1)/(1-alp) - (1-eta)*deltaFbar*hbar/(1-bet*(1-vtheta-sbar))*1/(1+vphi)  );

deltaWbar   = 1/(1-bet*(1-vtheta)*gamma)*hbar*1/(1-alp)*(-alp*wbar - (-1)*mrsbar);
Deltabar    = 1/(1-bet*(1-vtheta-sbar))*(wbar*hbar - mrsbar*hbar/(1+vphi) - b);

kappa       = qbar*bet*Jbar;
cbar        = ybar - kappa*vbar - Phi*nbar-gbar;

lambdabar   = (cbar*(1-habit))^(-sig);
kappaL      = mrsbar*lambdabar/hbar^vphi;


// - factor of proportionality between period profits and wages
Afactor =  (1-alp)/alp*wbar*hbar/(PsiLbar);



model(linear);


// Original Model Code:

// Monetary policy rule in the paper
Rt         =     (1-gamma_R)*gamma_Pi/4*(   Pit(-1) + Pit(-2) + Pit(-3)
                                           + Pit(-4) )
                                    + (1-gamma_R)*gamma_y/4*yt
                                    + gamma_R*Rt(-1)
                                    + emoneyt;

// consumption Euler equation
lambdat    =    lambdat(+1)  + Rt + ebt - Pit(+1);

// marginal utility of consumption
lambdat    =     - sig/(1-habit)*( ct - habit*ct(-1)  );

// New Keynesian Phillips curve (with inflation indexation)
Pit        =      price_index/(1+bet*price_index)*Pit(-1)
                + bet/(1+bet*price_index)*Pit(+1)
                + 1/(1+bet*price_index)*(1-omega)*(1-omega*bet)/omega*(mct);

// marginal cost
mct        =     xLt;

// new matches
mt         =     xi*ut
              + (1-xi)*vt;

// employment
nt         =     (1-vtheta)*nt(-1)
                + mbar/nbar*mt(-1);

// unemployment
nt         =   - ubar/(1-ubar)*ut;

// job-filling rate
qt         =     mt-vt;

// job-finding rate
st         =     mt-ut;

// newly optimized wage (wage setting FOC)
Jstart + deltaWt = Deltastart + deltaFt;

// hours FOC
wt         =     xLt + zt + (alp-1)*ht ;

// evolution of aggregate wage
wt         = 	 gamma*( wt(-1) - Pit + wage_index*Pit(-1) ) + (1-gamma)*wstart;

// deltaFt (-\partial surplus of firm/\partial wage)
deltaFt    = 	    (1-bet*(1-vtheta)*gamma)*( -alp/(1-alp)*wstart + 1/(1-alp)*(xLt + zt) )
        										+ bet*(1-vtheta)*gamma*(   -alp/(1-alp)*(wstart + wage_index*Pit - wstart(+1) - Pit(+1) )
                                                                           + deltaFt(+1) + lambdat(+1) - lambdat  );

// deltaWt ( \partial surplus of worker/\partial wage)
deltaWbar*deltaWt	=         -alp/(1-alp)*wbar*hbar*( -alp/(1-alp)*wstart + 1/(1-alp)*(xLt+zt) )
        						+ 1/(1-alp)*mrsbar*hbar*( (-1)*(1+vphi)/(1-alp)*wstart - lambdat + (1+vphi)/(1-alp)*(xLt+zt) )
                                + bet*(1-vtheta)*gamma/(1-bet*(1-vtheta)*gamma)*( (alp/(1-alp))^2*wbar*hbar - (1+vphi)/(1-alp)^2*mrsbar*hbar )
                                     *( wstart + wage_index*Pit - wstart(+1) - Pit(+1) )
                                + bet*(1-vtheta)*gamma*deltaWbar*( lambdat(+1) - lambdat + deltaWt(+1)  );


// value of firm that resets wage
Jbar*Jstart =     wbar*hbar/alp*( -alp*wstart + xLt + zt )
                + bet*(1-vtheta)*gamma/(1-bet*(1-vtheta)*gamma)*wbar*hbar*( wstart(+1) + Pit(+1) - wstart - wage_index*Pit  )
        		+ bet*(1-vtheta)*Jbar*( lambdat(+1)- lambdat + Jstart(+1) );

// surplus of worker that resets wage
Deltabar*Deltastart =   wbar*hbar*1/(1-alp)*( - alp*wstart + xLt + zt )
        			  - mrsbar*hbar/(1+vphi)*( (1+vphi)/(1-alp)*( (-1)*wstart+xLt+zt) - lambdat )
        			  + bet*(1-vtheta)*gamma/(1-bet*(1-vtheta)*gamma)*(wbar*hbar*alp/(1-alp)-1/(1-alp)*mrsbar*hbar)*( wstart(+1) + Pit(+1) - wstart - wage_index*Pit )
        			  - bet*gamma*sbar/(1-bet*(1-vtheta)*gamma)*(wbar*hbar*alp/(1-alp)-1/(1-alp)*mrsbar*hbar)*( wstart(+1) + Pit(+1) - wt - wage_index*Pit )
    	              + bet*(1-vtheta-sbar)*Deltabar*(  lambdat(+1) - lambdat + Deltastart(+1)  )
                      - bet*Deltabar*sbar*st;

// vacancy posting condition
kappa/qbar*( -qt )
            =     bet*gamma/(1-bet*(1-vtheta)*gamma)*wbar*hbar*( wstart(+1) - wt - wage_index*Pit + Pit(+1) )
                + bet*Jbar*( lambdat(+1) - lambdat  + Jstart(+1) );

// resource constraint
ybar*yt    =     cbar*ct + gbar*gt + vbar*kappa*vt + Phi*nbar*nt;

// production
yt         =     nt + zt  + alp*ht;


// annual inflation
Piannt     =         Pit     + Pit(-1) + Pit(-2) + Pit(-3) ;

// ****************** shocks ************************
// shock to discount factor
ebt        =     rho_eb*ebt(-1) +  sig_innoeb*inno_ebt;

// government spending
gt         =     rho_g*gt(-1) +  g_ ;

// monetary policy shock
emoneyt    =     rho_emoney*emoneyt(-1) + sig_monpol*interest_;

// productivity shock
zt         =     rho_z*zt(-1) +  sig_innoz*inno_zt;


end;

shocks;
var inno_ebt  = 1;
var inno_zt   = 1;
var interest_ = 1;
var g_ = sig_innog^2;
end;

//check;
//stoch_simul (irf = 0, ar=100, noprint);
stoch_simul (AR=100,IRF=0, noprint,nograph);