function [Fn,Fn2,fshift]= CalcFourierModes(f,t)
% CalcFourierModes - Calculates the Fourier series components of a uniformly sampled signal.
% This function was extracted from FourierAnalysisSVD2.m to be used globally.
% IMPORTANT: This function assumes the input signal 'f' is already sampled at a uniform time step.
%
% Inputs:
%   f - The uniformly sampled signal (e.g., flow waveform) as a row or column vector.
%   t - The corresponding time vector, assumed to be uniform.
%
% Outputs:
%   Fn      - The raw FFT output.
%   Fn2     - The shifted FFT output (zero-frequency component in the middle).
%   fshift  - The shifted frequency vector.

    a=size(t);
    N=a(1);
    if N==1
        t = t'; % Transpose if it's a row vector
        f = f';
        N = size(t,1);
        disp('Input vectors were transposed to be column vectors.')
    end

    % --- FFT Calculation ---
    % The FFT is performed directly on the input data, assuming it's uniformly sampled.
    T = t(end) - t(1);
    fs = N/T; % Corrected based on working version from user
    
    Fn=fft(f)/N;
    Fn2 = fftshift(Fn);
    
    % Create the frequency vector
    if mod(N, 2) == 0
        % Even number of points
        fshift = (-N/2:N/2-1)*(fs/N);
    else
        % Odd number of points
        fshift = (-(N-1)/2:(N-1)/2)*(fs/N);
    end
end