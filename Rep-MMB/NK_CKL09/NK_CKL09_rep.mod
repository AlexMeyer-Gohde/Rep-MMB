//**************************************************************************
// New Keynesian model with search and matching frictions in the labor market.
// Calibrated to the euro area in Christoffel, Kuester, Linzert (EER, 2009)
//**************************************************************************

// Model: NK_CKL09.mod

// Further references:
// Christoffel, K., Kuester, K., and T. Linzert (2009),
// "The role of labor markets for euro area monetary policy," European Economic Review 53, 908�936
//
// This is the version of the model calibrated in Section 3 of the paper.
//
// Last edited: 2011/05/11 by K. Kuester


var ct deltaFt deltaWt Deltastart ht Jstart
    lambdat mct mt nt Pit Piannt  qt
    Rt st ut vt wstart wt xLt yt
    // *** flex price/flex wage
    cflext deltaFflext deltaWflext Deltastarflext
    hflext Jstarflext lambdaflext mflext qflext
    sflext vflext wflext xLflext yflext
    // shocks
    ebargaint ekappat esept ebt
    eCt emoneyt gt zt;

    


varexo inno_ebargaint inno_ekappat inno_esept inno_ebt
       inno_eCt  inno_zt interest_ g_;



parameters

           bet epsilon habit sig vphi
           omega price_index
           xi eta gamma wage_index vtheta alp
           rho_eb rho_emoney rho_g rho_z
           rho_ebargain rho_ekappa rho_esep rho_eC
           gamma_R gamma_Pi gamma_y gamma_dy
           sig_innoeb sig_monpol sig_innog sig_innoz
           sig_innoeC sig_innoebargain sig_innoekappa sig_innoesep
           kappa Phi PhiL  Afactor b_wh
           // steady-state values
           cbar Deltabar deltaWbar gtildebar gbar hbar Jbar
           mbar mcbar mrsbar nbar qbar
           sbar ubar vbar wbar ybar
           // steady-state value flex-price economy
           deltaWflexbar;


// ***** parameter values
gamma_R     = 0.85;	// interest rate smoothing
gamma_Pi    = 1.5;  // weight on inflation
gamma_y     = 0.5;  // weight on output
gamma_dy    = 0.0;  //

bet     = 0.992;     // time-discount factor
epsilon = 11;        // own price elasticity of demand for differentiated good
vphi    = 2;         // inverse of Frisch elasticity
sig     = 1.5;       // CRRA
habit   = 0.6;       // habit persistence

omega       = 0.75;  // Calvo price stickiness
price_index = 0.0;  // degree of price indexation

xi          = 0.6;   // weight on unemployment in matching function
eta         = 0.5;   // bargaining power of workers
gamma       = 0.83;  // Calvo wage stickiness
wage_index  = 0.0;  // wage indexation
vtheta      = 0.03;  // rate of separation
alp         = 0.66;  // labor elasticity of production

rho_eb      = 0.85;     // persistence of preference shock
rho_g       = 0.79;     // persistence of gov spending shock
rho_z       = 0.64;     // persistence of productivity shock
rho_ebargain = 0;    // persistence of bargaining power shock

rho_ekappa   = 0;    // persistence of vacancy posting costs shock
rho_esep     = 0;    // persistence of separation rate shock
rho_eC       = 0;       // persistence of cost-push shock
rho_emoney   = 0;        // persistence of mon pol shock


sig_innoeb          = 0.191; // varphi = 1: 0.218; varphi=2: 0.191; varphi=10: 0.18
sig_innoeC          = 0.00; // std dev of cost-push shock
sig_innog           = 0.47; // std dev of inno to g shock (calculated below (depends on sample))
sig_monpol          = 0.1;  // std dev of monpol shock  (calculated below from Taylor rule)
sig_innoz           = 0.39; // std dev of tech shock (calculated below (depends on alpha))
sig_innoebargain    = 0.00; // std deviation of bargaining power shock
sig_innoekappa      = 0.00; // std deviation of vacancy posting shock
sig_innoesep        = 0.00; // std deviation of separation shock


// ********* calibration targets
ybar    = 1;                    // steady state output
hbar    = 1/3;                  // steady state hours worked
Afactor = 107.9365;             // value of Afactor
ubar    = 0.0916;               // steady state unemployment rate
qbar    = 0.7;                  // probability of finding a worker in a given month
gtildebar = 0.2;             // share of government spending in GDP in data;
b_wh    = 0.65;                 // replacement rate



// ********* ss for the given targets
PhiK_xLzhalp = 1-alp-0.01;

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

PsiLbar = (1-alp)/alp*wbar*hbar/(Afactor);
PhiK    = PhiK_xLzhalp*xLbar*zbar*hbar^alp;
PhiL    = xLbar*zbar*hbar^(alp) - wbar*hbar - PhiK - PsiLbar;
Phi     = PhiL+PhiK;

PsiLbar     = xLbar*zbar*hbar^(alp) - wbar*hbar - Phi;
PsiCbar     = (1-mcbar)*ybar;

Jbar		= 1/(1-bet*(1-vtheta))*PsiLbar;
deltaFbar   = 1/(1-bet*(1-vtheta)*gamma)*wbar*hbar;

b           = b_wh*wbar*hbar;

mrsbar      = (Jbar*eta/(1-bet*(1-vtheta)*gamma)*alp/(alp-1)*hbar*wbar - (1-eta)*deltaFbar/(1-bet*(1-vtheta-sbar))*(wbar*hbar-b))/ (Jbar*hbar*eta/(1-bet*(1-vtheta)*gamma)*(-1)/(1-alp) - (1-eta)*deltaFbar*hbar/(1-bet*(1-vtheta-sbar))*1/(1+vphi)  );

deltaWbar   = 1/(1-bet*(1-vtheta)*gamma)*hbar*1/(1-alp)*(-alp*wbar - (-1)*mrsbar);
deltaWflexbar = 1/(1-bet*(1-vtheta)*0)*hbar*1/(1-alp)*(-alp*wbar - (-1)*mrsbar);
Deltabar    = 1/(1-bet*(1-vtheta-sbar))*(wbar*hbar - mrsbar*hbar/(1+vphi) - b);

kappa       = qbar*bet*Jbar;
cbar        = (ybar - kappa*vbar - PhiL*nbar)*(1-gtildebar);      if cbar<=0; error('cbar<0'); end
gbar        = gtildebar/(1-gtildebar)*cbar;

lambdabar   = (cbar*(1-habit))^(-sig);
kappaL      = mrsbar*lambdabar/hbar^vphi;




model(linear);


// Original Model Code:

// Monetary policy rule in the paper
Rt         =      (1-gamma_R)*gamma_Pi/4*( Piannt)
                                   + (1-gamma_R)*gamma_y/4*(yt-1*yflext)
                                   + gamma_R*Rt(-1)
                                   + gamma_dy*(yt-yt(-1))
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
mct        =     eCt + xLt;

// new matches
mt         =     xi*ut
              + (1-xi)*vt;

// employment
nt         =     (1-vtheta)*nt(-1)
                + mbar/nbar*mt(-1)
                - vtheta*esept;

// unemployment
nt         =   - ubar/(1-ubar)*ut;

// job-filling rate
qt         =     mt-vt;

// job-finding rate
st         =     mt-ut;

// newly optimized wage (wage setting FOC)
Jstart + deltaWt = Deltastart + deltaFt - 1/(1-eta)*ebargaint;

// hours FOC
wt         =     xLt + zt + (alp-1)*ht ;

// evolution of aggregate wage
wt         = 	 gamma*( wt(-1) - Pit + wage_index*Pit(-1) ) + (1-gamma)*wstart;

// deltaFt (-\partial surplus of firm/\partial wage)
deltaFt    = 	    (1-bet*(1-vtheta)*gamma)*( -alp/(1-alp)*wstart + 1/(1-alp)*(xLt + zt) )
        		  + bet*(1-vtheta)*gamma*(   - alp/(1-alp)*(wstart + wage_index*Pit - wstart(+1) - Pit(+1) )
                                             + deltaFt(+1) + lambdat(+1) - lambdat  - vtheta/(1-vtheta)*esept(+1) );

// deltaWt ( \partial surplus of worker/\partial wage)
deltaWbar*deltaWt	=         -alp/(1-alp)*wbar*hbar*( -alp/(1-alp)*wstart + 1/(1-alp)*(xLt+zt) )
        						+ 1/(1-alp)*mrsbar*hbar*( (-1)*(1+vphi)/(1-alp)*wstart - lambdat + (1+vphi)/(1-alp)*(xLt+zt) )
                                + bet*(1-vtheta)*gamma/(1-bet*(1-vtheta)*gamma)*( (alp/(1-alp))^2*wbar*hbar - (1+vphi)/(1-alp)^2*mrsbar*hbar )
                                     *( wstart + wage_index*Pit - wstart(+1) - Pit(+1) )
                                + bet*(1-vtheta)*gamma*deltaWbar*( lambdat(+1) - lambdat + deltaWt(+1) - vtheta/(1-vtheta)*esept(+1) );


// value of firm that resets wage
Jbar*Jstart =     wbar*hbar/alp*( -alp*wstart + xLt + zt )
                + bet*(1-vtheta)*gamma/(1-bet*(1-vtheta)*gamma)*wbar*hbar*( wstart(+1) + Pit(+1) - wstart - wage_index*Pit  )
        		+ bet*(1-vtheta)*Jbar*( lambdat(+1)- lambdat + Jstart(+1) - vtheta/(1-vtheta)*esept(+1) );

// surplus of worker that resets wage
Deltabar*Deltastart =   wbar*hbar*1/(1-alp)*( - alp*wstart + xLt + zt )
        			  - mrsbar*hbar/(1+vphi)*( (1+vphi)/(1-alp)*( (-1)*wstart+xLt+zt) - lambdat )
        			  + bet*(1-vtheta)*gamma/(1-bet*(1-vtheta)*gamma)*(wbar*hbar*alp/(1-alp)-1/(1-alp)*mrsbar*hbar)*( wstart(+1) + Pit(+1) - wstart - wage_index*Pit )
        			  - bet*gamma*sbar/(1-bet*(1-vtheta)*gamma)*(wbar*hbar*alp/(1-alp)-1/(1-alp)*mrsbar*hbar)*( wstart(+1) + Pit(+1) - wt - wage_index*Pit )
    	              + bet*(1-vtheta-sbar)*Deltabar*(  lambdat(+1) - lambdat + Deltastart(+1)  )
                      - bet*Deltabar*sbar*st
                      - bet*Deltabar*vtheta*esept(+1);

// vacancy posting condition
kappa/qbar*( ekappat-qt )
            =     bet*gamma/(1-bet*(1-vtheta)*gamma)*wbar*hbar*( wstart(+1) - wt - wage_index*Pit + Pit(+1) )
                + bet*Jbar*( lambdat(+1) - lambdat  + Jstart(+1) );

// resource constraint
ybar*yt    =     cbar*ct + gbar*gt + vbar*kappa*(vt+ekappat) + PhiL*nbar*nt;

// production
yt         =     nt + zt  + alp*ht;


// annual inflation
Piannt     =         Pit     + Pit(-1) + Pit(-2) + Pit(-3) ;


// ********** flex price, flex wage economy ****************

// Consumption Euler equation
lambdaflext    =     - sig/(1-habit)*( cflext - habit*ct(-1)  );

// xLflext
0              =      eCt + xLflext;

// matches
mflext         =     xi*ut
                   + (1-xi)*vflext;

// job-filling rate
qflext         =     mflext-vflext;

// job-finding rate
sflext         =     mflext-ut;

// wage FOC
Jstarflext + deltaWflext = Deltastarflext + deltaFflext - 1/(1-eta)*ebargaint;



// hours FOC
wflext         =     xLflext + zt + (alp-1)*hflext;

// deltaFflext
deltaFflext    =       ( -alp/(1-alp)*wflext + 1/(1-alp)*(xLflext + zt) );

// deltaWflext
deltaWflexbar*deltaWflext
               =         -alp/(1-alp)*wbar*hbar*( -alp/(1-alp)*wflext + 1/(1-alp)*(xLflext+zt) )
                                - (-1)/(1-alp)*mrsbar*hbar*( (-1)*(1+vphi)/(1-alp)*wflext - lambdaflext + (1+vphi)/(1-alp)*(xLflext+zt) );

// Jstarflext
Jbar*Jstarflext =     wbar*hbar/alp*( -alp*wflext + xLflext + zt )
                                 + bet*(1-vtheta)*Jbar*( lambdaflext(1)- lambdaflext + Jstarflext(1) - vtheta/(1-vtheta)*esept(1) );

// Deltastarflext
Deltabar*Deltastarflext = wbar*hbar*1/(1-alp)*( - alp*wflext + xLflext+zt )
                                                - mrsbar*hbar/(1+vphi)*( (1+vphi)/(1-alp)*( (-1)*wflext+xLflext+zt)  - lambdaflext )
                                                + bet*(1-vtheta-sbar)*Deltabar*
                                                        (  lambdaflext(1) - lambdaflext + Deltastarflext(1)  )
                                                - bet*Deltabar*sbar*sflext
                                                - bet*Deltabar*vtheta*esept(1);

// vacancy posting
kappa/qbar*( ekappat-qflext )
                        =    bet*Jbar*( lambdaflext(1) - lambdaflext  + Jstarflext(1) );

// resource constraint
ybar*yflext    =     cbar*cflext + gbar*gt + vbar*kappa*(vflext+ekappat) + PhiL*nbar*nt;

// production
yflext         =     nt + zt  + alp*hflext;


// ****************** shocks ************************
// shock to discount factor
ebt        =     rho_eb*ebt(-1) +  sig_innoeb*inno_ebt;

// government spending
gt         =     rho_g*gt(-1) +  g_ ;

// monetary policy shock
emoneyt    =     rho_emoney*emoneyt(-1) + sig_monpol*interest_;

// productivity shock
zt         =     rho_z*zt(-1) +  sig_innoz*inno_zt;

// shock to bargaining power
ebargaint  =     rho_ebargain*ebargaint(-1) + sig_innoebargain*inno_ebargaint ;

// shock to vacancy posting costs
ekappat    =     rho_ekappa*ekappat(-1) +  sig_innoekappa*inno_ekappat;

// shock to separation rate
esept      =     rho_esep*esept(-1) +  sig_innoesep*inno_esept;

// cost-push shock
eCt        =     rho_eC*eCt(-1) +  sig_innoeC*inno_eCt;

end;

shocks;
var inno_ebt  = 1;
var inno_zt   = 1;
var interest_ = 1;
var g_   = sig_innog^2;
var inno_ebargaint = 1;
var inno_ekappat = 1;
var inno_esept = 1;
var inno_eCt   = 1;
end;

//check;
//stoch_simul (irf = 0, ar=100, noprint);
stoch_simul (AR=100,IRF=0, noprint,nograph);