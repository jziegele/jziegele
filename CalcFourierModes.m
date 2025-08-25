function [Fn,Fn2,fshift]= CalcFourierModes(f,t)
% CalcFourierModes - Calculates the Fourier series components of a signal.
% This function was extracted from FourierAnalysisSVD2.m to be used globally.
%
% Inputs:
%   f - The signal (e.g., flow waveform) as a row or column vector.
%   t - The corresponding time vector.
%
% Outputs:
%   Fn      - The raw FFT output.
%   Fn2     - The shifted FFT output (zero-frequency component in the middle).
%   fshift  - The shifted frequency vector.

    a=size(t);
    N_orig=a(1);
    if N_orig==1
        t = t'; % Transpose if it's a row vector
        f = f';
        N_orig = size(t,1);
        disp('Input vectors were transposed to be column vectors.')
    end

    % --- Resampling to a uniform time grid for FFT ---
    % The FFT algorithm assumes uniformly sampled data. If the input time vector t
    % is non-uniform, we must first resample the signal f onto a uniform grid.
    
    T = t(end) - t(1); % Total duration of the signal
    N = N_orig; % Use the same number of points for the uniform grid
    fs = (N-1)/T; % Correct sampling frequency for the uniform grid
    t_uniform = linspace(t(1), t(end), N)'; % Create the uniform time vector
    
    % Interpolate the signal onto the uniform time grid
    f_uniform = interp1(t, f, t_uniform, 'linear');

    % --- FFT Calculation ---
    % Now, perform the FFT on the uniformly sampled data.
    Fn=fft(f_uniform)/N;
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