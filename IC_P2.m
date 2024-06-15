
% Define the numerator and denominator coefficients for the transfer function H(s)
% These coefficients are specified for a high-order system.
numerator = [1.0055e+23, 0, 0, 0];   % Coefficients of the numerator polynomial of H(s)
denominator = [1,54356436.5638108,2.81252429808991e+16,9.58909261436301e+23,2.22068017524791e+32,3.38867908961806e+39,4.92231267110556e+47]; % Coefficients of the denominator polynomial of H(s)

% Create the transfer function H(s) using the 'tf' function
H = tf(numerator, denominator);

% Compute the poles of the transfer function H(s) and convert them to Hz
% Poles are the roots of the denominator polynomial of the transfer function.
poles_in_Hz = pole(H) / (2 * pi);

% Display the computed poles in Hz
disp('Poles of the transfer function in Hz:');
disp(poles_in_Hz);

% Plot the Bode plot of the transfer function H(s)
% The Bode plot shows the magnitude and phase response of the system over a range of frequencies.
figure;
bode(H), grid on;
title('Bode Plot of the Transfer Function H(s)');
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');

% Advanced Functionality:
% Plotting the Pole-Zero Map to visualize the location of poles and zeros in the complex plane.
figure;
pzmap(H);
grid on;
title('Pole-Zero Map of the Transfer Function H(s)');
xlabel('Real Part');
ylabel('Imaginary Part');

% Frequency Response Analysis: Gain and Phase Margin
% This provides insight into the stability and performance margins of the system.
[Gm, Pm, Wcg, Wcp] = margin(H);

% Display Gain Margin (Gm) and Phase Margin (Pm) along with their crossover frequencies (Wcg, Wcp)
disp('Gain Margin (dB):');
disp(20*log10(Gm)); % Convert gain margin from absolute to dB
disp('Phase Margin (degrees):');
disp(Pm);
disp('Gain Crossover Frequency (Hz):');
disp(Wcg / (2 * pi));
disp('Phase Crossover Frequency (Hz):');
disp(Wcp / (2 * pi));

% Nyquist Plot for Additional Stability Analysis
% This plot helps in determining the stability of the system by inspecting the encirclements around (-1,0).
figure;
nyquist(H);
grid on;
title('Nyquist Plot of the Transfer Function H(s)');
xlabel('Real Part');
ylabel('Imaginary Part');

% Frequency Response Analysis: Step Response
% This plot shows how the system responds to a step input, which is crucial for time-domain performance assessment.
figure;
step(H);
grid on;
title('Step Response of the Transfer Function H(s)');
xlabel('Time (seconds)');
ylabel('Amplitude');
