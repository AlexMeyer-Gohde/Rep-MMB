function [results] = dsge_backward_errors_condition(inputs)
%Matrix quadratic backward errors and conditioning number follow Higham and Kim (2001)
%0=A*X^2+B*X+C
%X is the n by n solvent
%A, B, and C are n by n matrices
%Alexander Meyer-Gohde
%24/03/2022
A=sparse(inputs.A);
B=sparse(inputs.B);
C=sparse(inputs.C);
D=sparse(inputs.D);
P=sparse(inputs.P);
Q=sparse(inputs.Q);
[ny,ne]=size(Q);
F=A*P+B;
results=NaN(15,5);
alpha=norm(A,'fro');
beta=norm(B,'fro');
gamma=norm(C,'fro');
delta=norm(D,'fro');
phi=norm(F,'fro');
P_F=norm(P,'fro');
xi=norm(P,'fro');
A_F=norm(A,'fro');
A_2=normest(A);
P_2=P^2;
P_2_F=norm(P_2,'fro');
P_Q=P*Q;
P_Q_F=norm(P_Q,'fro');
Q_F=norm(Q,'fro');
m_P_Q=[P Q];
m_P_Q_F=norm(m_P_Q,'fro');
P_m_P_Q=P*[P Q];
P_m_P_Q_F=norm(P_m_P_Q,'fro');



R_P=A*P^2+B*P+C;
R_Q=F*Q+D;
R_P_F=norm(R_P,'fro');
R_Q_F=norm(R_Q,'fro');
R_P_Q=[R_P R_Q];
R_P_Q_F=norm(R_P_Q,'fro');

min_svd_P_2=svds(P_2,1,'smallest')^2;
min_svd_P=svds(P,1,'smallest')^2;
min_svd_Q=svds(Q,1,'smallest')^2;
min_svd_P_Q=svds(P_Q,1,'smallest')^2;
min_svd_F=svds(F,1,'smallest')^2;
min_svd_A=svds(A,1,'smallest')^2;

%V=kron(speye(ny),F)+kron(P',A);
%I_F=kron(speye(ne),F);
%W=kron(speye(ne+ny),F)+kron([[P';Q'],sparse(ny+ne,ne)],A);

%[P Q_1 Q_2 Q_3 PQ]%
%[R RR H+2 H+2est H+2Factor mu V/F-1 diff-1est eig_sep psi psiest phi phiest  FE1 FE2


%R
%P R
results(1,1)=R_P_F;

%Q_1 R
results(1,2)=R_Q_F;

%Q_2 R
results(1,3)=results(1,2);

%Q_3 R
results(1,4)=results(1,3);

%PQ R
results(1,5)=R_P_Q_F;

%RR
%P RR
results(2,1)=results(1,1)/(alpha*P_2_F+beta*P_F+gamma);

%Q_1 RR
results(2,2)=results(1,2)/(phi*Q_F+delta);

%Q_2 RR
results(2,3)=results(1,3)/(alpha*P_Q_F+beta*Q_F+delta);

%Q_3 RR
results(2,4)=results(2,3);

%PQ RR
results(2,5)=results(1,5)/(alpha*P_m_P_Q_F+beta*m_P_Q_F+(gamma^2+delta^2)^(1/2));%Check

%H+2svds(A,1,'smallest')
%P H+2
results(3,1)=1/(alpha^2*min_svd_P_2+beta^2*min_svd_P+gamma^2)^(1/2);

%Q_1 H+2
results(3,2)=1/(phi^2*min_svd_Q+delta^2)^(1/2);

%Q_2 H+2
results(3,3)=1/(alpha^2*min_svd_P_Q+beta^2*min_svd_Q+delta^2)^(1/2);

%Q_3 H+2


%PQ H+2
results(3,5)=1/(alpha^2*(min_svd_P_Q+min_svd_P_2)+beta^2*(min_svd_Q+min_svd_P)+gamma^2+delta^2)^(1/2);

%H+2 est
%P H+2
%H_P=[alpha*kron(P_2',speye(ny)) beta*kron(P',speye(ny)) gamma*speye(ny^2)];
%results(4,1)=1/svds(H_P,1,'smallest');

%Q_1 H+2
%H_Q1=[phi*kron(Q',speye(ny)) delta*speye(ny*ne)];
%results(4,2)=1/svds(H_Q1,1,'smallest');

%Q_2 H+2
%H_Q2=[alpha*kron(P_Q',speye(ny)) beta*kron(Q',speye(ny)) delta*speye(ny*ne)];
%results(4,3)=1/svds(H_Q2,1,'smallest');

%Q_3 H+2


%PQ H+2
%H_PQ=[alpha*kron([P_2';P_Q'],speye(ny)) beta*kron([P';Q'],speye(ny)) gamma*kron([speye(ny);sparse(ne,ny)],speye(ny)) delta*kron([sparse(ny,ne);speye(ne)],speye(ny))];
%results(4,5)=1/svds(H_PQ,1,'smallest');

%H+2Factor
%P H+2Factor
results(5,1)=results(3,1)*(alpha*P_2_F+beta*P_F+gamma);

%Q_1 H+2Factor
results(5,2)=results(3,2)*(phi*Q_F+delta);

%Q_2 H+2Factor
results(5,3)=results(3,3)*(alpha*P_Q_F+beta*Q_F+delta);

%Q_3 H+2Factor
results(5,4)=results(5,3);

%PQ H+2Factor
results(5,5)=results(3,5)*(alpha*P_m_P_Q_F+beta*m_P_Q_F+(gamma^2+delta^2)^(1/2));%Check;


%mu
%P mu
results(6,1)=(alpha*P_2_F+beta*P_F+gamma)/(alpha^2*min_svd_P_2+beta^2*min_svd_P+gamma^2)^(1/2);

%Q_1 mu
results(6,2)=(phi*Q_F+delta)/(phi^2*min_svd_Q+delta^2)^(1/2);

%Q_2 mu
results(6,3)=(alpha*P_Q_F+beta*Q_F+delta)/(alpha^2*min_svd_P_Q+beta^2*min_svd_Q+delta^2)^(1/2);

%Q_3 mu
results(6,4)=results(6,3);

%PQ mu
results(6,5)=(alpha*P_m_P_Q_F+beta*m_P_Q_F+(gamma^2+delta^2)^(1/2))/(alpha^2*(min_svd_P_2+min_svd_P_Q)+beta^2*(min_svd_P+min_svd_Q)+gamma^2+delta^2)^(1/2);%Check;


%diff-1 est
%P V/F-1
%results(7,1)=1/svds(V,1,'smallest');

%Q_1 V/F-1
%results(7,2)=1/svds(F,1,'smallest');

%Q_2 V/F-1
%results(7,3)=results(7,2);

%Q_3 V/F-1
%results(7,4)=results(7,3);

%PQ V/F-1
%results(7,5)=1/svds(W,1,'smallest');

%diff-1 est CHECK
%P diff-1 est
[temp_P_FE in difp]=SYLG_allocated(ny,ny,full(real(F)),eye(ny),real(inputs.A),real(inputs.P)',full(real(R_P)),1);
results(8,1)=difp;

%Q_1 diff-1 est
[temp_Q_FE in difq]=SYLG_allocated(ny,ne,full(real(F)),eye(ne),zeros(ny,ny),zeros(ne,ne),full(real(R_Q)),1);
%results(8,2)=results(7,2);
results(8,2)=difq;

%Q_2 diff-1 est
results(8,3)=results(8,2);

%Q_3 diff-1 est
%results(7,4)=results(7,3);

%PQ diff-1 est
[temp_PQ_FE in difpq]=SYLG_allocated(ny,ny+ne,full(real(F)),eye(ny+ne),real(inputs.A),[[real(inputs.P)';real(inputs.Q)'],zeros(ny+ne,ne)],full(real(R_P_Q)),1);

results(8,5)=difpq;

%eigensep
eig_s=eig(full(P));
%X_unstable=-F\A;
%eig_u=eig(X_unstable);
%eig_u=1./eig_u;
eig_u=eig(full(F),-full(A));%eig(A,-F); check
eigval=[eig_s;eig_u];
l_1=eigval(abs(eigval)>1);
l_2=eigval(abs(eigval)<1);
results(9,1)=min(min(abs(l_1-l_2')));



%psi
%P psi
%temp_P=V\H_P;
%results(10,1)=norm(temp_P)/P_F;
%results(11,1)=normest(temp_P)/P_F;

%Q_1 psi
%temp=I_F\H_Q1;
%results(10,2)=norm(temp)/Q_F;
%results(11,2)=normest(temp)/Q_F;

%Q_2 psi
%temp=I_F\H_Q2;
%results(10,3)=norm(temp)/Q_F;
%results(11,3)=normest(temp)/Q_F;

%Q_3 psi
%temp=kron(Q',A)*[temp_P sparse(ny*ny,ny*ne)];
%temp=temp-[alpha*kron(P_Q',speye(ny)) beta*kron(Q',speye(ny)) sparse(ny*ne,ny*ny) delta*speye(ny*ne)];
%temp=I_F\temp;
%results(10,4)=norm(temp)/Q_F;
%results(11,4)=normest(temp)/Q_F;

%PQ psi
%temp=[alpha*kron([P_2';P_Q'],speye(ny))  beta*kron([P';Q'],speye(ny)) gamma*kron([speye(ny);sparse(ne,ny)],speye(ny)) delta*kron([sparse(ny,ne);speye(ne)],speye(ny))];
%temp=W\temp;
%results(10,5)=norm(full(temp))/m_P_Q_F;
%results(11,5)=normest(temp)/m_P_Q_F;


%phi
%P phi
results(12,1)=results(7,1)*(alpha*P_2_F+beta*P_F+gamma)/P_F;
results(13,1)=(1/(min_svd_F+min_svd_P*min_svd_A)^(1/2))*(alpha*P_2_F+beta*P_F+gamma)/P_F; %CHECK

%Q_1 phi
results(12,2)=results(7,2)*(phi*Q_F+delta)/Q_F;
results(13,2)=results(8,2)*(phi*Q_F+delta)/Q_F;%CHECK

%Q_2 phi
results(12,3)=results(7,3)*(xi*Q_F*A_F+alpha*P_Q_F+beta*Q_F+delta)/Q_F;
results(13,3)=results(8,3)*(xi*Q_F*A_F+alpha*P_Q_F+beta*Q_F+delta)/Q_F;%CHECK

%Q_3 phi
results(12,4)=results(7,1)*(alpha*P_Q_F+beta*Q_F)/Q_F+results(7,1)*results(7,2)*gamma*Q_F*A_F/Q_F+results(7,2)*delta/Q_F;
results(13,4)=results(8,1)*(alpha*P_Q_F+beta*Q_F)/Q_F+results(8,1)*results(8,2)*gamma*Q_F*A_F/Q_F+results(8,2)*delta/Q_F;%CHECK

%PQ phi
results(12,5)=results(7,1)*(alpha*P_m_P_Q_F+beta*m_P_Q_F+(gamma^2+delta^2)^(1/2))/m_P_Q_F;
results(13,5)=results(8,1)*(alpha*(P_2_F+P_Q_F)+beta*(P_F+Q_F)+(gamma+delta))/m_P_Q_F;%CHECK;%results(7,5)*(alpha*(P_2_F+P_Q_F)+beta*(P_F+Q_F)+(gamma+delta))/m_P_Q_F;%CHECK


%FE
%P FE
%temp_P=V\R_P(:);
%results(14,1)=normest(temp_P)/P_F;
results(15,1)=results(8,1)*(R_P_F/P_F); %CHECK

%Q_1 FE
%temp_Q=I_F\R_Q(:);
%results(14,2)=normest(temp_Q)/Q_F;
results(15,2)=results(8,2)*(R_Q_F/Q_F); %CHECK

%Q_2 FE
%results(14,3)=results(14,2);
results(15,3)=results(15,2); %CHECK

%Q_3 FE
%temp=kron(Q',F\A)*temp_P-temp_Q;
%results(14,4)=normest(temp,1)/Q_F;
results(15,4)=results(8,2)*((R_Q_F/Q_F)+results(7,1)*A_2*R_P_F);

%PQ FE
%temp=[kron([speye(ny);sparse(ne,ny)],speye(ny)) kron([sparse(ny,ne);speye(ne)],speye(ny))]*[R_P(:);R_Q(:)];
%results(14,5)=normest(W\temp)/m_P_Q_F;
results(15,5)=results(8,1)*R_P_Q_F/m_P_Q_F;%results(7,5)*R_P_Q_F/m_P_Q_F;
