P = -1e3;
Q = -6;
Vse=57.5;
Vsm=Vse*sqrt(2);
Vsdq=[sqrt(3/2)*Vsm; 0];
Vsd_s=Vsdq(1);
Vsq_s=Vsdq(2)+0.01;
Id_ref=(P*Vsd_s-Q*Vsq_s)/(Vsd_s^2+Vsq_s^2);
Iq_ref=(P*Vsq_s+Q*Vsd_s)/(Vsd_s^2+Vsq_s^2);
Id = 10.04;%Id_ref
Iq = 0.06;%Iq_ref


Zdd=(P/(Id^2 + Iq^2) - (2*Id*(Id*P + Iq*Q))/(Id^2 + Iq^2)^2)
Zdq=(Q/(Id^2 + Iq^2) - (2*Iq*(Id*P + Iq*Q))/(Id^2 + Iq^2)^2)
Zqd=(- Q/(Id^2 + Iq^2) - (2*Id*(Iq*P - Id*Q))/(Id^2 + Iq^2)^2)
Zqq=(P/(Id^2 + Iq^2) - (2*Iq*(Iq*P - Id*Q))/(Id^2 + Iq^2)^2)



Zdd=20*log10(abs(P/(Id^2 + Iq^2) - (2*Id*(Id*P + Iq*Q))/(Id^2 + Iq^2)^2))
Zdq=20*log10(abs(Q/(Id^2 + Iq^2) - (2*Iq*(Id*P + Iq*Q))/(Id^2 + Iq^2)^2))
Zqd=20*log10(abs(- Q/(Id^2 + Iq^2) - (2*Id*(Iq*P - Id*Q))/(Id^2 +Iq^2)^2))
Zqq=20*log10(abs(P/(Id^2 + Iq^2) - (2*Iq*(Iq*P - Id*Q))/(Id^2 + Iq^2)^2))

