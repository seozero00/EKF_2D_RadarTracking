%% w07 Extended Kalman Filter Example
%% 2D Target Tracking by a radar (nonlinear measurement eq)
%% Dr Seungkeun Kim (skim78@cnu.ac.kr)
clear all; close all; clc
d2r = pi/180; r2d = 1/d2r;
%% Dynamics
T = 0.5;         % sampling time: taking a picture of a ball every 1sec
F = [            % system matrix
    1 T T^2/2
    0 1 T
    0 0 1
    ];
Gamma = [T^3/6; T^2/2; T];  % Gamma (related to system noise)
% H = [1 0];      % Measurement matrix: ONLY POSITION 
%% extension to 2D
F = [
    F zeros(size(F));
    zeros(size(F)) F
    ];
Gamma = [
    Gamma zeros(size(Gamma));
    zeros(size(Gamma)) Gamma
    ];
% H = [
%     H zeros(size(H))
%     zeros(size(H)) H
%     ];
%% Radar position
% radar 1
% xr = 1500; yr = 0; % radar position
xr = 400; yr = 4000; % radar position
%% Initial Conditions
x(:,1) = [0;15;0; 0;0;1];        % true initial state
xp(:,1) = [-100;0;0; -100;0;0];  % guess of initial posteriori estimation
y(:,1) = fn_hx(xr,yr,xp(:,1));
nx = length(xp(:,1));           % number of state
Pp = 1e3*eye(nx);               % guess of initial error covariance
sigma_Pp(:,1) = sqrt(diag(Pp));
%% Noise
sigma_w = 0.1;                       % system noise (std of acceleration)
sigma_v = [30; 5*d2r];               % measurement noise (std of position sensor)
Q = 2*sigma_w^2*eye(size(Gamma,2));  % system noise covariance matrix
% cf) size(Q): 2x2/size(Gamma): 6x2
R = diag(sigma_v.^2);                % measurement noise covariance matrix
%% KF Routine
t = 0:T:100;
for i = 1:length(t)-1
    %% True dynamics
    x(:,i+1) = F*x(:,i) + Gamma*randn(size(Gamma,2),1)*sigma_w*0;  % system dynamics
    %% ============================================================
    %% Time update
    %% ============================================================
    %% Equation 1: Prediction of state
    x_ = F*xp(:,i);
    %% Equation 2: Prediction of covariance
    P_ = F*Pp*F' + Gamma*Q*Gamma';
    %% ============================================================
    %% Measurement update
    %% ============================================================
    %% measurement generation
    y(:,i+1) = fn_hx(xr,yr,x(:,i+1)) + diag(randn(size(sigma_v,1),1))*sigma_v;
    % cf) size(sigma_v): 2x1 / randn(size(sigma_v,1),1)=randn(2,1): 2x1
    % 이므로 곱셈을 위해 2x2를 만들어줘야하기 때문에 diag을 맨 앞에 붙여준 것임
    %% IMPORTANT!! - ANGLE MANIPULATION
    %% make the absolute value of angles under 180deg
    if abs(y(2,i+1)) > pi
        y(2,i+1) = y(2,i+1) - 2*pi*sign(y(2,i+1));
    end
    %% Compute Jacobian
    H = fn_JH(xr,yr,x_);
    %% Equation 3: Innovation Covariance
    S = H*P_*H'+R;
    %% Equation 4: Residual
    nu = y(:,i+1)-fn_hx(xr,yr,x_);
    %% Equation 5: Kalman gain
    K = P_*H'*S^-1;
    %% Equation 6: State update
    xp(:,i+1) = x_ + K*nu;
    %% Equation 7: Covariance update
    Pp = (eye(size(K*H)) - K*H)*P_;
    %% ============================================================
    %% storing error covariance for plotting
    sigma_Pp(:,i+1) = sqrt(diag(Pp));
    %% storing Kalman gain for ploting
    K_store(i) = norm(K);
    %% Plot
    if rem(i,10) == 0
        hold on
        plot(x(1,i),x(4,i),'k.')
        plot(xp(1,i+1),xp(4,i+1),'bs')
        ellipse(xp(1,i+1), xp(4,i+1), 3*sigma_Pp(1,i+1), 3*sigma_Pp(4,i+1))
        pause(0.2)
    end
end

%% Plot: position estimate x_1
figure; clf; hold on;
plot(xr,yr,'o','markersize',5);
plot(x(1,:),x(4,:),'.-b');         % plot the real plant behavior
plot(xp(1,:),xp(4,:),'.-r');       % plot the Kalman filter prediction over the plant
xlabel('x, m'); ylabel('y, m')
legend('Radar','real','KF','location','best')
%% Plot: Kalman gain norm
figure;
plot(t(1:end-1),K_store);
xlabel('time, s'); ylabel('Kalman gain norm')
%% Plot: measurement
figure;
subplot(121); plot(t,y(1,:))
xlabel('time, s'); ylabel('range, m')
subplot(122); plot(t,y(2,:)*r2d)
xlabel('time, s'); ylabel('azimuth, deg')
%% Plot: estimation error
figure;
subplot(121); plot(t,xp(1,:)-x(1,:))
xlabel('time, s'); ylabel('x position error, m')
subplot(122); plot(t,xp(4,:)-x(4,:))
xlabel('time. s'); ylabel('y position error, m')