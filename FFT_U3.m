clear
% Load data

n = 1;
Regular = load('Merged.csv');

while n>0

    clf
    % Model parameters
    dt = 5e-7;             % Sampling period
    N = 2401;             % Length of signal
    Fs = 1/(dt);            % Sampling frequency
    t = (0:N-1)*dt;        % Time vector

    %FFT
    Y_Regular = fft(Regular(:,n+1));
    A2_Regular = abs(Y_Regular/N);  % Double side spectrum
    % A1 = A2(1:N/2+1); % Even signal length
    A1_Regular = A2_Regular(1:(N+1)/2); % Odd signal length
    A1_Regular(2:end-1) = 2*A1_Regular(2:end-1); % Single side spectrum

    % f = Fs*(0:(N/2))/N;   % Even
    f = Fs*(0:((N-1)/2))/N;  % Odd
    grid on;

    % Plot freq. domain
    hold on
    plot(f, A1_Regular, 'Linewidth',2)
    title('Amplitude Spectrum of X(t)')
    xlabel('f (Hz)')
    ylabel('|A1(f)|')
    legend(['60120HZ-Node',num2str(n) ])
    set(gca, 'fontsize',16)
%     xlim([0 100000])
    saveas(gcf,['60120HZ-Node',num2str(n),'.png'])

    n = input('Enter node: ');
%     pause(1)
end
