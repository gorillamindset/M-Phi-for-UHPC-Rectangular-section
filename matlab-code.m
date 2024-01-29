clc
format long
storeMu = [];
storePhi = [];
%%%%-----  TO PLOT MULTIPLE GRAPHS -------%%%%%%
%%% MPhi Analysis for DR Rectangular section - Using Layer by Layer approach %%%
%%% compressive stregnth of concrete is as per UHPC Multilinear Model %%%
%%% Tensile strength of Concrete is considered as per UHPC BiLinear model %%%
%%% Tensile stregnth of UHPC after localisation is considered to be Dercreasing exp %%%

noofplots = input("enter no of plots:");
count = 0;
while count<=noofplots
    Ptc = input("enter Ptc:");
    Ptt = input("enter Ptt:");
    if Ptc==-1
        break
    end
    hold on;
    myfun(Ptc,Ptt);
end

function [list_phi,list_M] = myfun(Ptc,Ptt)
%Geometric Properties
B = 200;                    % Width of beam (mm)
D = 500;                    % Total depth of the beam (mm)
dc = 50;                    % Effective cover to the R/f (mm)
dcr = 25;                   % c/c distance between two layer of R/f (mm)
d = D - dc;                 % Effective depth of the beam (mm)
Asc = Ptc * B * D;          % Area of Compression steel reinforcement (mm^2)
Ast = Ptt * B * D;          % Area of Tension steel reinforcement (mm^2)
Ast1 = Ast/2;               % Area of Tension steel in Layer 1 (mm^2)
Ast2 = Ast/2;               % Area of Tension steel in Layer 2 (mm^2)
A = B*D;                         % Area of cross section (mm^2)
I = B*D.^3/12;                   % Moment of Inertia of the section (mm^4)
Y = D/2;                         % Max depth of N.A. (mm)
Z = I/Y;                         % Section modulus (mm^3)
t = 1;                           % Thickness of each layer (mm)

%Material Properties
Es = 200*1e3;                     % Young's modulus (MPa)
fc = 143.5;                       % Concrete compressive strength (MPa) 
fy = 500;                         % Yield strength (MPa)
ec = 0.002;                       % Ultimate concrete strain
fck = (1.25 * fc - 1.65 * 5);     % Characteristic comp strength (MPa)
Ect = 0.85 * 5000 * sqrt(fc);     % Modulus of elasticity of Concrete (MPa)
fcr = 0.7 * sqrt(fck);            % modulus of rupture
ecr = fcr / Ect;                  % rupture strain (at first tensile crack)
eccmax = 0.015;                  % Max compressive strain in concrete
ectmax = -0.025;                  % Max Tensile strain in concrete
estmax = 0.12;                    % Max strain in Tension R/f
escmax = 0.12;                    % Max strain in Compression R/f 
Vf = 2;                           % Percentage volumn of Fibers
Lf = 13;                          % Length of Fibers (mm)
Df = 0.2;                         % Dia of Fibers (mm)
S = 1;                            % Fiber Shape factor

%Other required Parameters
fcte = -0.65*sqrt(fc);
ecte = 1.1*(0.95*fcte)/Ect;
RI = Vf*Lf*S/Df/100;
fctl = 0.95*fcte*(RI)^0.25;
ectl = 20*ecte*(RI)^0.25;
fcm = fc + (9*RI);
Ecc = 4520*sqrt(fcm);        % Modulus of elasticity of Concrete (MPa)
e = 0.75 * (fcm/Ecc);
ecc1 = e;
ecc2 = 1.7*e;
ecc3 = 3.4*e;
ecc4 = 5.1*e;
ecc5 = 6.8*e;
fcc1 = 0.7*fcm;
fcc2 = fcm;
fcc3 = 0.38*fcm*RI;
fcc4 = 0.23*fcm*RI;
fcc5 = 0.14*fcm*RI;

nl = D / t;                % No of layers
ectop = 0;                 % Starting with ectop = 0
xu = d;                    % Starting iteration with xu(depth of NA) = d
list_phi = [];             % defining null matrix of phi
list_M = [];
while ectop < eccmax

    while true
        sum_Cc = 0;           % Defining zeros matrix for sum_Cc
        sum_Ct = 0;
        sum_Ts = 0;
        sum_M1 = 0;
        counter = 0;
        while counter <= nl
            counter = counter + 1;       
%      () = () + 1 - To start next loop or take next value
%       Counter - Iterative Integer varying from 1 to no. of layers           
            y_bar = xu - ((counter-1) * t + (t / 2));
            y_bar_cr = (xu - dc);
            y_bar_tr1 = (d - xu);
            y_bar_tr2 = (d - dcr - xu);
            ecl = (ectop / xu) * y_bar;
%             ecl = max(ecl,ectmax);
            es_cr = ((ectop / xu) * y_bar_cr);
            es_tr1 = ((ectop / xu) * y_bar_tr1);
            es_tr2 = ((ectop / xu) * y_bar_tr2);
            
            es_cr = min(es_cr,escmax);
            es_tr1 = min(es_tr1,estmax);
            es_tr2 = min(es_tr2,estmax);
%       y_bar = Centroidal distance of each layer from Xu
%       ybar_cR = Centroidal distance of compressive Reinforcement from xu
%       ybar_tR = Centroidal distance of Tensile Reinforcement from xu
%       ecl = Strain in each layer of conrete
%       es_cr/es_tr = strain in Comp and Tensile R/f

            if ecl >= 0 && ecl < ecc1
               fcc = (fcc1/ecc1) * ecl;
            elseif ecl < ecc2 && ecl > ecc1
                  fcc = fcc1 + (fcc2-fcc1)/(ecc2-ecc1)*(ecl-ecc1);
            elseif ecl < ecc3 && ecl > ecc2
                      fcc = fcc2 + (fcc3-fcc2)/(ecc3-ecc2)*(ecl-ecc2);
            elseif ecl < ecc4 && ecl > ecc3
                          fcc = fcc3 + (fcc4-fcc3)/(ecc4-ecc3)*(ecl-ecc3);
            elseif ecl < ecc5 && ecl > ecc4
                              fcc = fcc4 + (fcc5-fcc4)/(ecc5-ecc4)*(ecl-ecc4);
            else
               fcc = 0;
            end

%             if ecl > 0
%                 fcc = fc * (((2 * (ecl / ec)) - ((ecl / ec).^2)));
%             else
%                 fcc = 0;
%             end
            
            if ecl <= 0 && ecl > ecte
                   fct = (fcte/ecte) * ecl;
            elseif ecl >= ectl && ecl < ecte
                   fct = fcte+(fctl-fcte)/(ectl-ecte)*(ecl-ecte);
            elseif ecl >= ectmax && ecl < ectl
                   fct = fctl*(1+((ecl-ectl)/ectmax).^3)*exp(-2.2*((ecl-ectl)/ectmax));
            end
           
%      - Using Hognestad's model for calculating compressive
%        stresses in concrete i.e fcc      
%      - Assuming UHPC BiLinear model for concrete in Tension            

            if es_tr1 < (fy / Es)
                fst1 = es_tr1 * Es;
            elseif es_tr1 > (fy/Es) && es_tr1 < estmax
                fst1 = fy;
            else
                fst1 = 0;
            end

%             if es_tr1 < (fy / Es)
%                 fst1 = es_tr1 * Es;
%             else
%                 fst1 = fy;
%             end

%              if es_tr2 < (fy / Es)
%                 fst2 = es_tr2 * Es;
%               elseif es_tr2 > (fy/Es) && es_tr2 < estmax
%                 fst2 = fy;
%               else
%                  fst2 = 0;
%              end
             
            if es_tr2 < (fy / Es)
                fst2 = es_tr2 * Es;
            else
                fst2 = fy;
            end
            
            if es_cr < (fy / Es)
                fsc = es_cr * Es;
            else
                fsc = fy;
            end
%        fst, fsc = Stress in Tensile and Comp R/f.              

            Cc = fcc * B * t;
            Ct = fct * B * t;
            Cs = fsc * Asc;
            Ts1 = fst1 * Ast1;
            Ts2 = fst2 * Ast2;
            
            M1 = (Cc * y_bar) + (Ct * y_bar);
            sum_M1 = sum_M1 + M1;

            sum_Cc = sum_Cc + Cc;
            sum_Ct = sum_Ct + Ct;
%       Cc = Compressive force in each layer of Concrete
%       Ct = Tensile force in each layer of Concrete 
%       Cs / Ts = Comp and Tensile force in R/f
%       Mu1 = Moment due to force in each layer of concrete ( Cc and Ct)
%       sum_Cc = sum of all compressive forces in Concrete
%       sum_Ct = sum of all Tensile forces in Concrete
%       sum_M1 = sum of all moments due to each layer in concrete
        end
        
        M2 =  (Cs * y_bar_cr);
        M3 = (Ts1 * y_bar_tr1);
        M4 = (Ts2 * y_bar_tr2);
        M = (sum_M1 + M2 + M3 + M4) / 1e6;
%       M = Total moment due to Reinforcements and concrete.
%       Lever arm is conisdered from xu.

        if abs(sum_Cc + sum_Ct) + abs(Cs) - abs(Ts1) - abs(Ts2) < 1
       %fprintf('for ectop:%f, parameters: xu:%f, sum_Cc:%f, sum_Ct:%f, Mu:%f\n', ectop, xu, sum_Cc, sum_Ct, Mu);
%      Equilibrium condition of forces - C=T -->> Cc+Cs-Ct-Ts = 0
            break;
        end
        xu = xu - 0.0001 * d;
%       Increment in xu is in 1e-4 mm.
    end

    
    phi = (ectop / xu) * 1e6;
    fprintf('for ectop:%f, ecl:%f, es_tr1:%f, parameters: xu:%f, sum_Cc:%f, sum_Ct:%f, Mu:%f, phi:%f\n', ectop, ecl, es_tr1, xu, sum_Cc, sum_Ct, M,phi);
    list_phi = [list_phi, phi];
    list_M = [list_M, M];
    ectop = ectop + 0.00001;
    n_i = eccmax / 0.00001;
%       phi = Curvature due to each xu - from varying top strain (ectop)    
%       list_phi = values of phi for strain varying from 0 - ecmax
%       list_Mu = values of Mu for strain varying from 0 - ecmax
%       With increment of 1e-5. 
end


%       Plot = Graph of Moment vs Curvature for different strains
hold on;
plot(list_phi, list_M);
xlabel('Phi - 1/mm 10^-6');
ylabel('M - KNm');
title('M vs Phi for Rect section - UHPC - Layer by Layer');
grid on;
end



