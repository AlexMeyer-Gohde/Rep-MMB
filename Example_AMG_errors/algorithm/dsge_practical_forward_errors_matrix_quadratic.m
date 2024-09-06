function [results] = dsge_practical_forward_errors_matrix_quadratic(inputs)
%Matrix quadratic backward errors and conditioning number follow Higham and Kim (2001)
%0=A*X^2+B*X+C
%X is the n by n solvent
%A, B, and C are n by n matrices
%Alexander Meyer-Gohde
%24/03/2022
if inputs.nstatic>0
A=[inputs.A_static;inputs.AA];%sparse(inputs.A_full);
B=[inputs.B_static inputs.B_rest;zeros(inputs.ndynamic,inputs.nstatic) inputs.BB];%sparse(inputs.B_full);
C=[inputs.C_static;inputs.CC];%sparse(inputs.C_full);
else
    A=inputs.AA;
    B=inputs.BB;
    C=inputs.CC;
end
A=[zeros(inputs.endo_nbr,inputs.npred+inputs.nstatic) A];
C=[zeros(inputs.endo_nbr,inputs.nstatic) C zeros(inputs.endo_nbr,inputs.nfwrd)];
% matrix_quadratic.nstatic=nstatic;
% matrix_quadratic.endo_nbr=n;
% matrix_quadratic.nfwrd=nfwrd;
% matrix_quadratic.nstatic=nstatic;
% matrix_quadratic.nfwrd=nfwrd;
% matrix_quadratic.npred=npred;
% matrix_quadratic.nboth=nboth;
% matrix_quadratic.nsfwrd=nsfwrd;
% matrix_quadratic.nspred=npred + nboth;
% matrix_quadratic.ndynamic=ndynamic;
% 
%  A=matrix_quadratic.A_full=[zeros(n,npred+nstatic) A];
%  B=matrix_quadratic.B_full=B;
%  C=matrix_quadratic.C_full=[zeros(n,nstatic) C zeros(n,nfwrd)];
%ny=inputs.endo_nbr;
P=inputs.X;%sparse(inputs.X);
%[ny,~]=size(P);
F=A*P+B;
results=NaN(8,5);
P_F=norm(P,'fro');



R_P=F*P+C;%A*P^2+B*P+C;
R_P_F=norm(R_P,'fro');


%H_P=kron(speye(ny),F)+kron(P',A);

%[P Q_1 Q_2 Q_3 PQ]%
%[be_lb be_ub factor psi phi sep p_psi p_phi



%sep
%P sep
%results(4,1)=svds(H_P,1,'smallest');

[A1,D1,P,Q] = qz(full(F),full(A));
[U,B1] = schur(full(-P),'complex');

C1=zeros(size(A1,1),size(B1,1));
F1=zeros(size(A1,1),size(B1,1));
E1=eye(size(B1));
separation = lapack('ZTGSYL', 'N', 3, size(A1,2), size(B1,2), A1, size(A1,1), B1, size(B1,1), C1, size(C1,1), D1, size(D1,1), E1, size(E1,1), F1, size(F1,1), 1, 1, 1, 1, ones(size(A1,2)+size(B1,2)+6,1), 0);
results(4,1)=separation{18};





%p_psi
%P p_psi
%temp=H_P\R_P(:);

%results(7,1)=normest(temp)/P_F;

C1=P*R_P*U;
temp2=lapack('ZTGSYL', 'N', 0, size(A1,2), size(B1,2), A1, size(A1,1), B1, size(B1,1), C1, size(C1,1), D1, size(D1,1), E1, size(E1,1), F1, size(F1,1), 1, 1, 1, 1, ones(size(A1,2)+size(B1,2)+6,1), 0);
temp2=Q*temp2{9}*U';
temp2=temp2(:);
results(7,1)=normest(temp2)/P_F;
%p_phi
%P p_phi
results(8,1)=(1/results(4,1))*R_P_F/P_F;


