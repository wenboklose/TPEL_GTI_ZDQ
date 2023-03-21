% Calculates the DQ impedance from the ABC impedance
% Ex:
% R = 30;
% C = 101e-6;
% Habc = R/(R*s*C+1);
% HDQ = JF_DQFromABC(Habc, wline)
%
function HDQ = JF_DQFromABC(Habc, wline)

HC = JF_CMOD(Habc, wline);
HS = JF_SMOD(Habc, wline);


HDQ = [HC, -HS; HS, HC]

%Chop the extra states out since the transfer function matrix should have
%the same denominator.
HDQ = minreal(HDQ);         %

HDQ = balreal(HDQ);         %