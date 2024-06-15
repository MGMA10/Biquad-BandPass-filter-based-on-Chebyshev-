clear all;
clc;
close all;
cap=10^(-11);
Max_PBF_rip=1.4;
Min_PBF_rip=20;
Fpl=10*10^6;
Fph=20*10^6;
Fsl=(20/3)*10^6;
Fsh=30*10^6;
% Passband Bandwidth @w, Center Frequency @center, and Quality Factor Qbp
W_BW = 2 * pi * ( Fph - Fpl );
W_c = 2* pi * sqrt (Fph * Fpl);
Q_bp= ( W_c/W_BW);
% Step 2: Normalized Passband Frequency ωp, & Normalized Stopband Frequency ωs, of the Prototype Lowpass Filter
ohm_p= 1;
ohm_s= (2 * pi *(Fsh - Fsl))/W_BW;
% Step 3: Filter Order n of the Prototype Lowpass Filter
E_p = sqrt (10^(Max_PBF_rip/10)-1);
E_s = sqrt (10^(Min_PBF_rip/10)-1);
n=ceil(acosh(E_s/E_p)/acosh(ohm_s/ohm_p));
% Step 4: Normalized Chebyshev Poles of the Prototype Lowpass Filter
a=(1/n)*asinh(1/E_p);
for i=1:n
    p_real(i)= sin ((2*i-1)*pi/2/n)*sinh(a);
    p_imag(i)= cos ((2*i-1)*pi/2/n)*cosh(a);
    p_k (i)=p_real(i)+j* p_imag(i);

end

A=real(p_k (2));
B=p_k (1)+p_k (3);
C=p_k (1)*p_k (3);
D=A*C;
E=1/W_BW ;
F= W_c^2 /W_BW ;

S6=1;
S5=A/E+B/E;
S4=((2*F*E+C)*E+E^2*F+B*E*A)/E^3;
S3=((2*F*E+C)*A+B*E*F+F*B*E)/E^3;
S2=((2*F*E+C)*F+E*F*F+A*F*B)/E^3;
S1=(F*F*B+A*F*F)/E^3;
S0=(F*F*F)/E^3;
G = nthroot((D/E^3),3);
p = [S6 S5 S4 S3 S2 S1 S0];
r = roots(p)

numerator = [G];
denominator1 = [1,-1*(r(1)+r(2)),r(1)*r(2)];
H1 = tf(numerator, denominator1)

denominator2 = [1,-1*(r(3)+r(4)),r(3)*r(4)];
H2 = tf(numerator, denominator2)

denominator3 = [1,-1*(r(6)+r(5)),r(5)*r(6)];
H3 = tf(numerator, denominator3)

% for Biquad filter topology
% let R1 = 1k
R1=1000;
Cap1=1/R1/sqrt(r(1)*r(2))
Cap2=1/R1/sqrt(r(3)*r(4))
Cap3=1/R1/sqrt(r(5)*r(6))

R2S1=1/(Cap1*(r(1)+r(2))*-1)
R2S2=1/(Cap2*(r(3)+r(4))*-1)
R2S3=1/(Cap3*(r(5)+r(6))*-1)


