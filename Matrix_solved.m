function[A B C Indicator_Variables NY NX NK ]= Matrix(v)%
%% Legend:

%y__= y
%pi_ = \pi_t
%i__=i_t

%% Indicator
% CONTROLS ARE in t and t+1
% ENDOGENOUS STATES are in t-1 and t
% EXOGENOUS STATES are in t and t+1
Indicator_Controls=['y_','pi','i_'];
Indicator_En_States=['']; % this is foreward looking model
Indicator_Ex_States=['a_','nu'];


Indicator_Variables=[Indicator_Controls,Indicator_En_States,Indicator_Ex_States];

%% Dimensions
NY = length(Indicator_Variables)/2;
NX = length(Indicator_Ex_States)/2;
NK = NX+length(Indicator_En_States)/2;

%% Give the numbers

for i=0:length(Indicator_Variables)/2-1
eval(['ind_' Indicator_Variables(2*i+1:2*(i+1)) '=' num2str(i+1) ';']); % y is first, pi is second...
eval(['ind_' num2str(i+1) '=' num2str(i+1) ';' ]);

end

%%

%% Initialize the Matrices
A = zeros(NY,NY);
B = zeros(NY,NY);
C = zeros(NY,NX);


%% Euler Equation
% y_t+1+ ((1+varphi)/(sigma+varphi)) a_{t+1}+ 1/sigma*pi_{t+1}=y_t+ ((1+varphi)/(sigma+varphi)) a_{t}+ 1/sigma*i_{t}
A(ind_1,ind_y_)=1;
A(ind_1,ind_pi)=1/v.sigma;
A(ind_1,ind_a_)=((1+v.varphi)/(v.sigma+v.varphi));

B(ind_1,ind_y_)=1;
B(ind_1,ind_i_)=1/v.sigma;
B(ind_1,ind_a_)=((1+v.varphi)/(v.sigma+v.varphi));

%% Philips Curve
% beta*pi_{t+1}=pi_t - lambda*(sigma+varphi)y_t
A(ind_2,ind_pi)=v.beta;

B(ind_2,ind_y_)=-v.lambda*(v.sigma+v.varphi);
B(ind_2,ind_pi)=1;

%% Monetary Policy Rule
% 0 = i_t - phi_y y_t - phi_pi *pi_t (- nu)

B(ind_3,ind_i_)=1;
B(ind_3,ind_y_)=-v.phi_y;
B(ind_3,ind_pi)=-v.phi_pi;
B(ind_3,ind_nu)=-1;

%% Technological Shock
%a_{t+1}= rho_a a_t + eps_a
A(ind_4,ind_a_)=1;

B(ind_4,ind_a_)=v.rho_a;

C(ind_4,1)=1;

%% Monetary Policy Shock
%nu_{t+1}= rho_nu nu_t + eps_nu
A(ind_5,ind_nu)=1;

B(ind_5,ind_nu)=v.rho_nu;

C(ind_5,2)=0.25;

end
