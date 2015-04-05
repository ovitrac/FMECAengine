function P = ModifiedGrainMethod(P0,T0,T,Tb,Tm)
%MODIFIEDGRAINMETHOD calculates the liquid vapor pressure as MPBPWIN (from EPIsuite v. 4.x or above)
%   P = ModifiedGrainMethod(PT0,T0,T,Tb [,Tm])
%   T0 = reference temperature (K)
%   P0 = liquid vapor pressure at T0  (atm)
%   T  = requested temperature (K)
%   Tb = normal boiling temperature (K)
%   Tm = melting temperature (default = 0) (K)
%   P  = liquid vapor pressure at T  (atm)
%
% EXAMPLE of anisole:
%         Experimental Database Structure Match:
%           Name     :  ANISOLE
%           CAS Num  :  000100-66-3
%           Exp MP (deg C):  -37.5 
%           Exp BP (deg C):  153.7 
%           Exp VP (mm Hg):  3.54E+00 
%                  (Pa   ):  4.72E+002
%           Exp VP (deg C):  25 
%           Exp VP ref    :  AMBROSE,D ET AL. (1976) 
% 
%         SMILES : O(c(cccc1)c1)C
%         CHEM   : Benzene, methoxy-
%         MOL FOR: C7 H8 O1 
%         MOL WT : 108.14
%         ------------------------ SUMMARY MPBPWIN v1.43 --------------------
%         Boiling Point:  149.16 deg C (Adapted Stein and Brown Method)
%         Melting Point:  -55.86 deg C (Adapted Joback Method)
%         Melting Point:  -26.57 deg C (Gold and Ogle Method)
%         Mean Melt Pt :  -41.21 deg C (Joback; Gold,Ogle Methods)
%           Selected MP:  -41.21 deg C (Mean Value)
% 
%         Vapor Pressure Estimations (25 deg C):
%           (Using BP: 153.70 deg C (exp database))
%           (MP not used for liquids)
%             VP:  3.67 mm Hg (Antoine Method)
%               :  489 Pa  (Antoine Method)
%             VP:  3.08 mm Hg (Modified Grain Method)
%               :  411 Pa  (Modified Grain Method)
%             VP:  4.52 mm Hg (Mackay Method)
%               :  603 Pa  (Mackay Method)
%           Selected VP:  3.38 mm Hg (Mean of Antoine & Grain methods)
%                      :  450 Pa  (Mean of Antoine & Grain methods)
%
%{
    % estimations at -20°C
    T0 = convertunit('C','K',25,'abstemp')
    P0 = convertunit('mmHg','atm',3.08)
    Tb = convertunit('C','K',149.16,'abstemp')
    Tm = convertunit('C','K',-55.86,'abstemp')
    T  = convertunit('C','K',-20,'abstemp')
    P = ModifiedGrainMethod(P0,T0,T,Tb,Tm)
    PinPa = convertunit('atm','Pa',P)
%}
%
%  Values at -20°C from EPIsuite --- SUMMARY MPBPWIN v1.43 ---
%         Vapor Pressure Estimations (-20 deg C):
%           (Using BP: 153.70 deg C (exp database))
%           (MP not used for liquids)
%             VP:  0.103 mm Hg (Antoine Method)
%               :  13.8 Pa  (Antoine Method)
%             VP:  0.0868 mm Hg (Modified Grain Method)
%               :  11.6 Pa  (Modified Grain Method)
%             VP:  0.147 mm Hg (Mackay Method)
%               :  19.6 Pa  (Mackay Method)
%           Selected VP:  0.0951 mm Hg (Mean of Antoine & Grain methods)
%                      :  12.7 Pa  (Mean of Antoine & Grain methods)
%
%
%
% Modified Grain Method
% This method is a modification and significant improvement of the modified Watson method.
% It is applicable to solids, liquids and gases.  The modified Grain method equations are:
% Source: Lyman, W.J.  1985.   In: Environmental Exposure From Chemicals. Volume I.
%                                  Neely,W.B. and Blau,G.E. (eds), Boca Raton, FL: CRC Press, Inc., Chapter 2.

% INRA\migration v 2.0 - 04/04/2015 - Olivier Vitrac - rev.

% Revision history

% arg check
if nargin<5, Tm = 0; end
lnP0_liq_rhs = log(P0) - solidcorrection(T0);
lnP_liq = lnP0_liq_rhs .* RHSliquid(T)./RHSliquid(T0);
P = exp(lnP_liq + solidcorrection(T));


% RHS (convent within bracket) liquid vapor pressure
function rhs = RHSliquid(TK)
    Tp  = TK./Tb;
    m = 0.4133 - 0.2575 .* Tp;
    rhs = 1 - ((3-2*Tp).^m)./Tp - (2*m*(3-2*Tp).^(m-1))*log(Tp);
end

% solid vapor pressure
function lnPsolid = solidcorrection(TK)
    R = 82.057; % cm3 atm/mol K
    Tpm = TK./Tm;
    m   = 0.4133 - 0.2575 .* Tpm;
    if TK<Tm
        lnPsolid = 0.6*log(R*Tm) .* ( 1 - ((3-2*Tpm).^m)./Tpm  - (2*m*(3-2*Tpm).^(m-1))*log(Tpm) );
    else
        lnPsolid = 0;
    end
end


end
