% -------------------- Parámetros de configuración --------------------
fm = 100000; % Frecuencia de muestreo (Hz)
tm = 1/fm; % Período de muestreo (s)
ls = 200; % Largo de la señal
f_c = 1000; % Frecuencia de la señal (Hz)
f_s = 5000; % Frecuencia de muestreo teórico (Hz)
t_s = 1/f_s; % Período de muestreo teórico (s)
tau = 0.5*t_s; % Duración del pulso (s)
d = tau/t_s; % Ciclo de trabajo
% -------------------- Vectores de tiempo y señal --------------------
t = (0:ls-1)*tm;
m_t = sin(2*pi*f_c*t);
% -------------------- Auxiliares --------------------
r = floor(t_s/tm);
s = floor(tau/tm);
% -------------------- Muestreo Natural --------------------
s_nat = zeros(1,length(t));
for i = 1:r:length(t)
    if i+s <= length(t)
        s_nat(i:i+s) = 1;
    else
        s_nat(i:end) = 1;
    end
end
m_t_nat = m_t .* s_nat;
% -------------------- Muestreo Instantáneo --------------------
m_t_inst = zeros(1, length(t));
indices_muestreo = 1:r:length(t);
m_t_inst(indices_muestreo) = m_t(indices_muestreo);
% -------------------- PCM (Modulación por Pulsos Codificados) --------------------
N = 4; % Número de bits por palabra PCM (configurable)
Nniveles = 2^N; % Cantidad de niveles de cuantización
m_min = min(m_t_inst);
m_max = max(m_t_inst);
% Cuantización uniforme
delta = (m_max - m_min) / (Nniveles - 1); % Paso de cuantización
m_t_cuant = round((m_t_inst - m_min) / delta) * delta + m_min; % Señal cuantizada
% Error de cuantización
e_cuant = m_t_inst - m_t_cuant;
% -------------------- FFT --------------------
L = length(t); % Número de muestras
frecuencia = (0:L/2) * (fm / L); % Frecuencia en Hz
% FFT de cada señal
T_m_t = fft(m_t);
T_m_t_nat = fft(m_t_nat);
T_m_t_inst = fft(m_t_inst);
% Cálculo del espectro unilateral
P2_m_t = abs(T_m_t/L);
P2_m_t_nat = abs(T_m_t_nat/L);
P2_m_t_inst = abs(T_m_t_inst/L);
P1_m_t = P2_m_t(1:L/2+1);
P1_m_t_nat = P2_m_t_nat(1:L/2+1);
P1_m_t_inst = P2_m_t_inst(1:L/2+1);
% Duplicar amplitudes de los componentes no DC
P1_m_t(2:end-1) = 2*P1_m_t(2:end-1);
P1_m_t_nat(2:end-1) = 2*P1_m_t_nat(2:end-1);
P1_m_t_inst(2:end-1) = 2*P1_m_t_inst(2:end-1);
% -------------------- GRÁFICAS --------------------
% Figura 1: Señales en el tiempo
figure;
subplot(3,1,1);
plot(t, m_t, 'b');
title('Señal Original');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;
subplot(3,1,2);
plot(t, m_t_nat, '-r');
title('Muestreo Natural');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;
subplot(3,1,3);
stem(t(indices_muestreo), m_t(indices_muestreo), 'or');
title('Muestreo Instantáneo (Puntos)');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;
% Figura 2: Transformadas de Fourier
figure;
subplot(3,1,1);
plot(frecuencia, P1_m_t, 'b');
title('Espectro de la Señal Original');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
grid on;
subplot(3,1,2);
plot(frecuencia, P1_m_t_nat, '-r');
title('Espectro del Muestreo Natural');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
grid on;
subplot(3,1,3);
plot(frecuencia, P1_m_t_inst, '-g');
title('Espectro del Muestreo Instantáneo');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
grid on;
% Figura 3: PCM y cuantización
figure;
subplot(2,1,1);
plot(t, m_t, 'b', t(indices_muestreo), m_t_inst(indices_muestreo), 'or', t(indices_muestreo), m_t_cuant(indices_muestreo), 'xg');
title('Señal Original, PAM Instantáneo y Cuantizado');
legend('Original', 'Muestreo Instantáneo', 'Cuantizado');
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on;
subplot(2,1,2);
stem(t(indices_muestreo), e_cuant(indices_muestreo), 'r');
title('Error de Cuantización PCM');
xlabel('Tiempo (s)');
ylabel('Error');
grid on;