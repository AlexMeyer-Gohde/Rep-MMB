//**************************************************************************
// This file has been transformed from its MMB version to be a replication 
// .mod file in June 2022.
//
//**************************************************************************
// A New Comparative Approach to Macroeconomic Modeling and Policy Analysis
//
// Volker Wieland, Tobias Cwik, Gernot J. Mueller, Sebastian Schmidt and
// Maik Wolters
//
// Working Paper, 2011
//**************************************************************************

// Model: HK_FPP11

// Further references: Funke, Paetz and Pytlarczyk (2011)
// Stock market wealth effects in an estimated DSGE model for Hong Kong
// In Economic Modelling 28 (2011) 316-334


// Model replication coded by: Xiaobei He, xihe@wiwi.uni-frankfurt.de
//                             Tong Wang, Tong.Wang@hof.uni-frankfurt.de

// Note: There have been two typos in the original code confirmed by the authors. Therefore, the
// impulse response functions from this model file do in general not match the impulse responses
// reported in the paper. The two adjustments refer to i) the parameter value of gamma (0.13 instead of 0.03), and
// ii) the sign of the foreign demand shock in the equation for the natural
// level of output (minus sign instead of + sign in front of the cofficient).

// Last edited: Aug. 19, 2011


var x q_dach r rr pi_H a y_stern pi s e y y_n shock_eta;




varexo epsa epsy mu_p epseta;




parameters

klein_omega gamma alpha gross_theta sigma_alpha psi gross_gamma_a
gross_gamma_y gross_gamma_null gross_epsilon beta betaS lambda_x kappa_alpha
phi_pi phi_x rho_eta rho_a rho_y rho_q r_steady lambda nu mu theta omega_C phi
tau zeta eta phi2;


// all parameter are calibrated according to the original paper, see p723


rho_eta = 0.775;
klein_omega = 0.33;
phi = 0.33;
tau = 0.1234;
theta =0.6658;
gamma =0.13;
gross_epsilon = 0.015;
beta = 0.995;
rho_a =  0.8549;  //0.808;
rho_q =  0.7750;  //0.808;
rho_y =  0.8306;  //0.683;
alpha = 0.5047;   //0.6;
mu = 0.1;
nu = mu + log (1 - alpha);
r_steady = 0.01;
eta = 0.5;
zeta = 0.5;
omega = zeta + (1 - alpha) * (eta - 1);
gross_theta = omega - 1;
sigma_alpha = 1 / ((1 - alpha) + alpha * omega);
omega_C = (((1 + r_steady) * (1 + gross_epsilon) * mu) / ((r_steady + (1 + r_steady) * gross_epsilon) * (1 + mu) * (1 - klein_omega)));
psi = gamma * ((1 - beta * (1 - gamma)) / (1 - gamma)) * omega_C;
betaS = beta / (1 + psi);
gross_gamma_a = (1 + phi) / (sigma_alpha + phi);
gross_gamma_y = gross_theta * sigma_alpha / (sigma_alpha + phi);
gross_gamma_null = (1 + psi - alpha * gross_theta * sigma_alpha);
lambda_x = (1 + gross_epsilon - betaS) * (((1 + phi - mu) / mu) + alpha * sigma_alpha);
phi2 = (theta + tau  * (1 - theta * (1 - betaS)))^(-1);
lambda = (1 - tau) * (1 - theta) * (1 - betaS * theta) * phi2;
kappa_alpha = lambda * (sigma_alpha + phi);
phi_pi = 1.5;
phi_x = 0;



model(linear);




// Original Model Code:

//euler equation
x = (sigma_alpha / gross_gamma_null ) * x(+1) + (psi / gross_gamma_null) * q_dach
    - (1 / gross_gamma_null) * (r - pi_H(+1) - rr);

//natural interest rates
rr = (sigma_alpha * rho_a + psi - gross_gamma_null) * gross_gamma_a * a
    + ((sigma_alpha * rho_y + psi - gross_gamma_null) * gross_gamma_y + (gross_theta * sigma_alpha * (rho_y - 1))) * alpha * y_stern;

//stock price dynamics
q_dach = (betaS / (1 + gross_epsilon))  * q_dach(+1) - (lambda_x / (1 + gross_epsilon)) * x(+1) - (r - pi_H(+1) - rr) + shock_eta;

//NKPC
pi_H = phi2 * (theta * betaS * pi_H(+1) + tau * pi_H(-1)) + kappa_alpha * x + mu_p;

//monetary policy
// e = 0;
r = phi_pi * pi_H;

//consumer prices
//pi = p - p(-1);

//producer prices
//pi_H = p_H - p_H(-1);

//terms of trade
s = sigma_alpha * (y - y_stern);

//natural output
y_n = gross_gamma_a * a - alpha * gross_gamma_y * y_stern;

//output gap
x = y - y_n;

//CPI
pi = pi_H + alpha * (s - s(-1));

//exchange rate
e-e(-1) = s-s(-1) + pi_H;

//shocks
a = rho_a * a(-1) + epsa;
y_stern = rho_y * y_stern(-1) + epsy;
shock_eta = rho_eta * shock_eta(-1) + epseta;

end;


shocks;
var epsa; stderr 2.6684;
var epsy; stderr 0.6148;
var epseta; stderr 2.9396;
var mu_p; stderr 0.9667;
end;

%stoch_simul (irf = 0, ar=100);
stoch_simul (AR=100,IRF=0, noprint,nograph);