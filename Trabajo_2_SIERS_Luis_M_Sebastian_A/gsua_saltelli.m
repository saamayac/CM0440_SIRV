function S = gsua_saltelli(S)
% S = gsua_saltelli(S)
% S                         Structure array with input and output fields
% FIELDS
% S.model                   Simulation model
% S.M                       Sample matrix (NxNp) of parameters (one parameter by column and one sample by row)
% S.ynom                    Reference vectorial (1xNt) or scalar output: experimental or nominal (ynom = [] for nominal)
% scalar_characteristic     Scalar function for computing scalar sensitivity
% S.N                       Sample size (number of sets of parameters to test)
% S.Np                      Number of factors
% S.S1                      First-order global sensitivity indices with output = y (NpxNt vectorial or Npx1 scalar output)
% S.S1s                     Scalar first-order global sensitivity indices with y = scalar_fun (Npx1)
% S.ST                      Total global sensitivity indices with output = y (NpxNt vectorial or Npx1 scalar)
% S.STs                     Scalar total global sensitivity indices with y = scalar_fun (Npx1)
% S.Y                       Output for every factor sample: vectorial time response (NxNt or N*NpxNt) or scalar (N*Npx1)
% S.Ys                      Scalar outputs for y = scalar_fun (Nx1 or N*Npx1)
% S.ynom                    Nominal time response (if input S.ynom = [])
% S.t                       Time simulation vector with fixed step (1xNt)
% S.Nsim                    Number of simulations
% S.tsim                    Simulation time in datatime format
%
% Global sensitivity and uncertainty analysis using GSUA Toolbox
% https://bit.ly/Matlab_GSUA
% (c) Carlos Mario Vélez S. 2022
% Universidad EAFIT, Medellin, Antioquia, Colombia
% https://sis-control.blogspot.com/

M = S.M;
ynom = S.ynom; % It is used intrinsically in eval(S.scalar_characteristic)
t = S.t;
N = S.N;
Np = S.Np;
N = ceil(N/2); % dim(M) = 2N x Np, dim(A) = dim(B) = dim(AB) = N x Np
S.Nsim = (2+Np)*N;
Nt = length(t);
Y = zeros(2*N,Nt); % Vectorial output
Ys = zeros(2*N,1); % Scalar output
if strcmp(S.model_type,'Simulink model') == 1
    load_system(S.model)
    model_workspace = get_param(S.model, 'ModelWorkspace');
end
tic

for i=1:2*N
    if i==100 || i==0.2*2*N || i==0.4*2*N || i==0.6*2*N || i==0.8*2*N
        perc = i/S.Nsim;
        [h1,m1,s1] = hms(seconds(toc*(1-perc)/perc));
        [h2,m2,s2] = hms(datetime(datetime('now')) + seconds(toc*(1-perc)/perc));
        disp(['Elapsed percent: ' num2str(perc*100,'%04.1f') '%  Remaining computing time: ' num2str(h1,'%02.0f') 'h:' num2str(m1,'%02.0f')...
            'm:' num2str(s1,'%02.0f') 's' '  Estimated stop time: ' num2str(h2,'%02.0f') ':' num2str(m2,'%02.0f') ':' num2str(s2,'%02.0f')])
    end
    switch S.model_type
        case 'Simulink model'
            assignin(model_workspace,'x', M(i,:));
            simout = sim(S.model,'ReturnWorkspaceOutputs','on');
            y = simout.yout.signals.values';
        case 'M-file model'
            [y,~] = eval([S.model '(M(i,:))']);
    end
    Y(i,:) = y; 
    Ys(i) = eval(S.scalar_characteristic);
end
V = var(Y,1);
Vs = var(Ys, 1);
A = M(1:N,:);
B = M(N+1:2*N,:);
YA = Y(1:N,:);
YB = Y(N+1:2*N,:);
YAs = Ys(1:N,:);
YBs = Ys(N+1:2*N,:);
S1 = zeros(Np,Nt); ST = zeros(Np,Nt); S1s = zeros(Np,1); STs = zeros(Np,Nt);
for i=1:Np
    perc = ((i-1)*N+2*N)/S.Nsim;
    [h1,m1,s1] = hms(seconds(toc*(1-perc)/perc));
    [h2,m2,s2] = hms( datetime(datetime('now')) + seconds(toc*(1-perc)/perc));
    disp(['Elapsed percent: ' num2str(perc*100,'%04.1f') '%  Remaining computing time: ' num2str(h1,'%02.0f') 'h:' num2str(m1,'%02.0f')...
        'm:' num2str(s1,'%02.0f') 's' '  Estimated stop time: ' num2str(h2,'%02.0f') ':' num2str(m2,'%02.0f') ':' num2str(s2,'%02.0f')])
    ABi = A;
    ABi(:,i) = B(:,i);
    YABi = zeros(N,Nt); YABis = zeros(N,1);
    for j = 1:N
        switch S.model_type
            case 'Simulink model'
                assignin(model_workspace,'x', ABi(j,:));
                simout = sim(S.model,'ReturnWorkspaceOutputs','on');
                y = simout.yout.signals.values';
            case 'M-file model'
                y = eval([S.model '(ABi(j,:))']);
        end
        YABi(j,:) = y;
        YABis(j) = eval(S.scalar_characteristic);
    end
    % First-order sensitivity indices
    S1(i,:) = sum(YB.*(YABi - YA),1)./(N*V);
    S1s(i) = sum(YBs.*(YABis - YAs),1)./(N*Vs);
    % Total sensitivity indices
    ST(i,:) = sum((YA - YABi).^2,1)./(2*N*V);
    STs(i) = sum((YAs - YABis).^2,1)./(2*N*Vs);
end
% Normalization of total sensitivity indices
S.ST = ST./sum(ST,1);
S.STs = STs/sum(STs);
S.S1 = S1;
S.S1s = S1s;
S.Y = Y;
S.Ys = Ys;
[th,tm,ts] = hms(seconds(toc)); 
S.tsim = [th tm ts];
disp('Elapsed percent: 100%') 

end