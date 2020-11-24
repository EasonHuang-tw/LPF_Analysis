%% Sensor Data
Fs = 400;                    % Sampling frequency
T = 1/Fs;                     % Sampling period
L = 10000;                     % Length of signal
t = (0:L-1)*T;                % Time vector
Fd = 40;                      %cut-off frequency

%% signal for analysis
x1 = cos(2*pi*20*t);          % First row wave
x2 = cos(2*pi*50*t);         % Second row wave
x3 = 0.3*cos(2*pi*75*t);         % Third row wave
x4 = cos(2*pi*100*t);         % Third row wave
x6=wgn(1,L,0.001);
x5 =  x1+x3;          %mix

%% ================= 1st order filter ===============
RC = 1/(2*pi*Fd);
alpha = T/(RC + T);

%% bode analysis
H = tf([1],[rc 1]);
[mag_1,phase_1,wout_1] = bode(H,{0,100*2*pi});
% wrange = [min(wout); max(wout)];                            % Units: rad/sec
% frange = wrange/(2*pi);                                     % Units: Hz

%% num analysis


y1 = zeros(1,L);
for i = 1:L
    if i > 1
        y1(i) = y1(i-1)*(1-alpha)+x5(i)*alpha;
        %disp(y1(i));
        %pause(0.1);
    end
end

%% ================= 2nd order filter ==============
Wd = 2*pi*(Fd/Fs);            %Wd = 2.pi.fd/fs
alpha = tan(Wd/2);            %alpha = tan(wd/2)
k = alpha^2/(1+sqrt(2)*alpha + alpha^2);%gain
D = 1+sqrt(2)*alpha+alpha^2; %Denominator
a1 = 2*(alpha^2-1)/D;
a2 = (1-sqrt(2)*alpha+alpha^2)/D;
b1 = 2;
b2 = 1;

%% bode analysis
z = tf('z',T);
H = k*(z^2+b1*z+b2)/(z^2+a1*z+a2);
[mag_2,phase_2,wout_2] = bode(H,{0,100*2*pi});

%% num analysis


y2 = zeros(1,L);
for i = 1:L
    if i > 2
        y2(i) = k*x5(i)+k*b1*x5(i-1)+k*b2*x5(i-2)-a1*y2(i-1)-a2*y2(i-2); %output
        %disp(y1(i));
        %pause(0.1);
    end
end

%% plot
% raw filt ideal
X = [x5;x1;y1;y2];

figure('Name', 'plot');
subplot(3,1,1)
plot(t(1:L/50),X(1,1:L/50));
title(['Raw data'])

subplot(3,1,2)
hold on;
plot(t(1:L/50),X(3,1:L/50));
plot(t(1:L/50),X(2,1:L/50));
legend('filted','ideal') ;
title('compare "1st filter" vs ideal ')

subplot(3,1,3)
hold on;
plot(t(1:L/50),X(4,1:L/50));
plot(t(1:L/50),X(2,1:L/50));
hold off;
legend( 'filted','ideal') ;
title('compare "2nd filter" vs ideal ')

%% Bode plot

figure('Name', 'bode');

subplot(4,1,1)
plot(wout_1/(2*pi), squeeze(mag_1))
grid
xlabel('Frequency (Hz)')
ylabel('magnitude (x)')
title('1st order Bode Plot')
subplot(4,1,2)
plot(wout_1/(2*pi), squeeze(phase_1))
grid
xlabel('Frequency (Hz)')
ylabel('phase (deg)')


subplot(4,1,3)
plot(wout_2/(2*pi), squeeze(mag_2))
grid
xlabel('Frequency (Hz)')
ylabel('magnitude (x)')
title('2 order Bode Plot')
subplot(4,1,4)
plot(wout_2/(2*pi), squeeze(phase_2))
grid
xlabel('Frequency (Hz)')
ylabel('phase (deg)')



%% FFT plot
figure('Name', 'FFT');
n = 2^nextpow2(L);
dim = 2;
Y = fft(X,n,dim);
P2 = abs(Y/L);
P1 = P2(:,1:n/2+1);
P1(:,2:end-1) = 2*P1(:,2:end-1);


subplot(4,1,1)
plot(0:(Fs/n):(Fs/2-Fs/n),P1(1,1:n/2))
title('Raw data in Frequency Domain')

subplot(4,1,2)
plot(0:(Fs/n):(Fs/2-Fs/n),P1(2,1:n/2))
title('Ideal data in Frequency Domain')
for i=3:4
    subplot(4,1,i)
    plot(0:(Fs/n):(Fs/2-Fs/n),P1(i,1:n/2))
    title([num2str(i-2),'order filter in the Frequency Domain'])
end
