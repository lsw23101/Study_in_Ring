% Example: Inverted pendulum 의 SISO 모델


clear all;

M = 0.5;
m = 0.2;
b = 0.1;
I = 0.006;
g = 9.8;
l = 0.3;

p = I*(M+m)+M*m*l^2; %denominator for the A and B matrices

A0 = [0      1              0           0;
     0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
     0      0              0           1;
     0 -(m*l*b)/p       m*g*l*(M+m)/p  0];
B0 = [     0;
     (I+m*l^2)/p;
          0;
        m*l/p];

C = [1 0 0 0];

D = [0;
     0];


% sampling time
Ts = 0.05;


% discretize
sysC = ss(A0,B0,C,[]);
sysD = c2d(sysC, Ts);
A = sysD.A;
B = sysD.B;

% dimensions

n=4; m=1; l=1;

% controller design
Q = eye(n);
R1 = eye(m);
R2 = eye(l);
[~, K, ~] = idare(A,B,Q,R1,[],[]);
K = -K;
[~, L, ~] = idare(A.', C.', Q, R2, [], []);
L = L.';

% (F,G,H): resulting controller
F = A + B*K - L*C;
G = L;
H = K;

% plant initial state
xp0 = [1; 1; 0.1; 0.1];
% controller initial state
xc0 = [0; 0; 0; 0];

%%%% controller conversion %%%%
% observability matrix
On = obsv(F,H);

% Toeplitz matrix
Tn = zeros(m, n*l);
for i = 1:n-1
    tmp = [H*F^(i-1)*G, Tn(m*(i-1)+1:end,1:l*(n-1))];
    Tn = [Tn; tmp];
end

% (flipped) controllability matrix [F^(n-1)G, ..., FG, G]
Cn = F^(n-1)*G;
for i = 2:n
    Cn = [Cn, F^(n-i)*G];
end

% converted form: u(k)=Hu*[u(k-n);...;u(k-1)]+Hy*[y(k-n);...;y(k-1)]
Hu = H*F^n*pinv(On)
Hy = H*(Cn - F^n*pinv(On)*Tn)


%%%% find initial input-output trajectory %%%%
% yini = [y(-n);...;y(-1)], uini = [u(-n);...;u(-1)]
yini = Cn\xc0;
uini = Tn*yini;
Yini = reshape(yini,[],n);
Uini = reshape(uini,[],n);


%%%%%%%%%%%% 기존 conversion 파일 그대로 %%%%%%%%%%%% 


%% Simulation
iter = 500;

% variables for simulation with original controller
xp = xp0;
xc = xc0;
u = [];
y = [];

% quantization parameters
r = 0.0001;
s = 0.0001;

% quantization of control parameters
qHu = round(Hu/s);
qHy = round(Hy/s);

% 1x4 행렬 >> 4x1 열벡터 형태로
qHu = qHu(:);
qHy = qHy(:);


% variables for simulation with converted & quantized controller

qXp = xp0;
qXp_ = xp0;

qU = round(Uini/r);
qY = round(Yini/r);

qU_ = round(Uini/r);
qY_ = round(Yini/r);


rY = [];
rU = [];

rY_ = [];
rU_ = [];


U_cin = [];
exp_U_cin = [];


% 여기서 Hu와 Hy 행렬 늘리기

% 정수 부분
int_part_Hu = fix(Hu);
int_part_Hy = fix(Hy);

% 소수 부분 4자리 정수화
frac_part_Hu = round((Hu - int_part_Hu) / s); 
frac_part_Hy = round((Hy - int_part_Hy) / s);

% 정수 - 소수 - 정수 - 소수 붙히기 >> 8x1 열벡터
Hu_ = reshape([int_part_Hu; frac_part_Hu], [], 1);  % 8x1
Hy_ = reshape([int_part_Hy; frac_part_Hy], [], 1);  % 8x1

% 정수 정수 소수 소수 형태로 늘리기 >> 16x1 열벡터
Hu_expanded = repelem(Hu_, 2);  % 16x1
Hy_expanded = repelem(Hy_, 2);  % 16x1


% e.g.
%
% Hy =
% 
%    1.0e+03 *
% 
%     0.5652   -1.8822    2.0701   -0.7531
% 
% 
% Hy_expanded =
% 
%          565
%          565
%         2060
%         2060
%        -1882
%        -1882
%        -2499
%        -2499
%         2070
%         2070
%         1209
%         1209
%         -753
%         -753
%         -757
%         -757
% 



for i = 1:iter
    %%%%%%%%% plant + original controller %%%%%%%%%
    y = [y, C*xp(:,i)];
    u = [u, H*xc(:,i)];
    xp = [xp, A*xp(:,i) + B*u(:,i)];
    xc = [xc, F*xc(:,i) + G*y(:,i)];





    %%%%%%%%% ARX controller %%%%%%%%%
    
    % cotroller (Hadamard product)
    U_cin_temp = qHu .* reshape(qU(:, end-n+1:end), [], 1) + qHy .* reshape(qY(:, end-n+1:end), [], 1);
    U_cin = [U_cin, U_cin_temp];           % U_cin 저장 

    U_cin_scalar = sum(U_cin_temp);          % inner sum

    % actuator
    rU = [rU, r * s * U_cin_scalar];           % 마지막 값에 스케일링

    % sensor
    rY = [rY,C*qXp(:,i)];

    % state update
    qY = [qY,round(rY(:,end)/r)];
    qU = [qU,round(rU(:,end)/r)];
    qXp = [qXp,A*qXp(:,i)+B*rU(:,end)];



    %%%%%%%%% Expanded Controller %%%%%%%%%

   % 인아웃풋 데이터 처리

   % [ u(-4)의 정수 u(-4)의 소수 u(-3)의 정수 u(-3)의 소수 ...  u(-1)의 정수 u(-1)의 소수 ]
   % 16x1 사이즈 열벡터
    
    % === u 처리 ===
    % 마지막 n개의 u값 추출
    qU_latest = qU_(:, end-n+1:end);  % 4x1
    
    % 실수값 복원
    u_latest_real = r * qU_latest;
    
    % 정수 및 소수 분리
    int_part_u = fix(u_latest_real);
    frac_part_u = round((u_latest_real - int_part_u) / r);
    
    % 8x1 만들기
    expanded_u = reshape([int_part_u; frac_part_u], [], 1);
    
    % 각 항목을 두 번 반복해서 16x1 만들기
    expanded_u_16 = repelem(expanded_u, 2);  
    
    
    % === y 처리 ===
    % 마지막 n개의 y값 추출
    qY_latest = qY_(:, end-n+1:end);  % 4x1
    
    % 실수값 복원
    y_latest_real = r * qY_latest;
    
    % 정수 및 소수 분리
    int_part_y = fix(y_latest_real);
    frac_part_y = round((y_latest_real - int_part_y) / r);
    
    % 8x1 만들기
    expanded_y = reshape([int_part_y; frac_part_y], [], 1);
    
    % 각 항목을 두 번 반복해서 16x1 만들기
    expanded_y_16 = repelem(expanded_y, 2); 


    
    % cotroller (Hadamard product)
    exp_U_cin_temp = Hu_expanded .* expanded_u_16 + Hy_expanded .* expanded_y_16;
    exp_U_cin = [exp_U_cin, exp_U_cin_temp];

    scaling_factors = [r*s; s; r; 1];               % 4가지 스케일링 계수
    scaling_vector = repmat(scaling_factors, 4, 1); % 16x1 벡터로 확장
    
    exp_U_cin_scalar = sum(exp_U_cin_temp .* scaling_vector);

    % actuator
    rU_ = [rU_, r * s * U_cin_scalar];           % 마지막 값에 스케일링

    % sensor
    rY_ = [rY_, C*qXp_(:,i)];


    
    % 입출력과 상태 업데이트
    qY_ = [qY,round(rY_(:,end)/r)];
    qU_ = [qU,round(rU_(:,end)/r)];
    qXp_ = [qXp_,A*qXp_(:,i)+B*rU_(:,end)];

end


figure(1)
plot(Ts*(0:iter-1), u, 'LineWidth', 1)       
hold on
plot(Ts*(0:iter-1), rU_, '-', 'LineWidth', 1) 
plot(Ts*(0:iter-1), rU, '--', 'LineWidth', 1) 


title('Control input u')
legend('original','quantized', 'expanded')

figure(2)
plot(Ts*(0:iter-1), y, 'LineWidth', 1)        
hold on
plot(Ts*(0:iter-1), rY_, '-', 'LineWidth', 1) 
plot(Ts*(0:iter-1), rY, '--', 'LineWidth', 1) 


title('Plant output y')
legend('original 1', 'quantized 1', 'expanded')



% q 사이즈 정해야하는 u/(rs) 크기 비교 

max(U_cin(:)) 
max(exp_U_cin(:))


% Todo : lattigo 라이브러리에 적용