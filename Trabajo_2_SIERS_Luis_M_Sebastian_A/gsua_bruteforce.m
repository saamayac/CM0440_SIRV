function S = gsua_bruteforce(S)
% S = gsua_bruteforce(S)
%
% S                         Structure array with input and output fields
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
ynom = S.ynom; % It is used intrinsically in "eval(S.scalar_characteristic)"
t = S.t;
N = S.N;
Np = S.Np;
S.Nsim = N^2*Np;
Yi_nom = repmat(ynom,N,1);
Nt = size(t,2);
V1 = zeros(Np,Nt);
VT = zeros(Np,Nt);
V1s = zeros(Np,1);
VTs = zeros(Np,1);
if strcmp(S.model_type,'Simulink model') == 1
    load_system(S.model)
    model_workspace = get_param(S.model, 'ModelWorkspace');
end
tic

for i=1:Np
    Yi = zeros(N^2,Nt);
    Yis = zeros(N^2,1);
    Mi = repmat(M,N,1);
    for j=1:N
        Mi((j-1)*N+1:j*N,i) = M(j,i)*ones(N,1);
    end
    for j = 1:N^2
        if j==10 || j==0.3*N^2 || j==0.7*N^2
            perc = ((i-1)*N^2+j)/S.Nsim;
            [h1,m1,s1] = hms(seconds(toc*(1-perc)/perc));
            [h2,m2,s2] = hms(datetime(datetime('now')) + seconds(toc*(1-perc)/perc));
            disp(['Elapsed percent: ' num2str(perc*100,'%04.1f') '%  Remaining computing time: ' num2str(h1,'%02.0f') 'h:' num2str(m1,'%02.0f')...
                'm:' num2str(s1,'%02.0f') 's' '  Estimated stop time: ' num2str(h2,'%02.0f') ':' num2str(m2,'%02.0f') ':' num2str(s2,'%02.0f')])
        end

        switch S.model_type
            case 'Simulink model'
                assignin(model_workspace,'x', Mi(j,:));
                simout = sim(S.model,'ReturnWorkspaceOutputs','on');
                y = simout.yout.signals.values';
            case 'M-file model'
                [y,~] = eval([S.model '(Mi(j,:))']);
        end
        Yi(j,:) = y;
        Yis(j) = eval(S.scalar_characteristic);
    end

    EYij = zeros(N,Nt); % Vectorial case
    EYijs = zeros(N,1); % Scalar case
    for j=1:N
        EYij(j,:) = mean(Yi((j-1)*N+1:j*N,:),1);
        EYijs(j,:) = mean(Yis((j-1)*N+1:j*N,:));
    end
    V1(i,:) = var(EYij,1);
    V1s(i) = var(EYijs);

    YTi = zeros(N^2,Nt);
    YTis = zeros(N^2,1);
    MTi = zeros(N^2,Np);
    for j=1:N
        MTi((j-1)*N+1:j*N,:) = Mi(j:N:N^2,:); % For verification
        YTi((j-1)*N+1:j*N,:) = Yi(j:N:N^2,:);
        YTis((j-1)*N+1:j*N) = sum( (YTi((j-1)*N+1:j*N,:)-Yi_nom).^2,2 );
    end
    EYTij = zeros(N,Nt);
    EYTijs = zeros(N,1);
    for j=1:N
        EYTij(j,:) = var(YTi((j-1)*N+1:j*N,:),1);
        EYTijs(j,:) = var(YTis((j-1)*N+1:j*N,:));
    end
    VT(i,:) = mean(EYTij,1);
    VTs(i) = mean(EYTijs);
end

Y = zeros(N,Nt);
Ys = zeros(N,1);
for j=1:N
    Y(j,:) = Yi((j-1)*N+j,:);
    Ys(j) = Yis((j-1)*N+j);
end
V = var(Y,1);
Vs = var(Ys);
S.S1 = V1./V;
S.ST = (VT./V)./sum((VT./V),1);
S.S1s = V1s/Vs;
S.STs = (VTs/Vs)/sum(VTs/Vs);
[h,m,s] = hms(seconds(toc)); S.tsim = [h m s];
S.Y = Y;
S.Ys = Ys;
disp('Elapsed percent: 100%') 

end