function hcmod = JF_SMOD(h,w)

if( nargin < 2)
    w = 500*2*pi;
end

if( nargin < 1 )
    s = tf('s');
    h = ss(1/s);
end

%Convert / force to statespace
h = ss(h);

%Shift frequency parameter
a = w;

%Element that performs frequency shift
modvector = diag(diag(h.a)*0+a);

%Anew = [ h.a - modvector, h.a*0; h.a*0, h.a+modvector];
Anew = [ h.a,   modvector;  -modvector, h.a];
Bnew = [ h.b; h.b];
%Cnew = [ -h.c/(2j), h.c/(2j)];
Cnew = [ h.c/(2), -h.c/(2)];
Dnew = h.d*0;

%Make resultant state space object
hcmod  = ss(Anew, Bnew, Cnew, Dnew);
