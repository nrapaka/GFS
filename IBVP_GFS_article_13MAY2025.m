% solves an initial boundary value problem
function IBVP()
    clear all; clc;
%     figure(31); clf;
%     figure(32); clf;
%     figure(3); clf;
%     global error;
%     error = zeros(1,0);
    format longE;
    global amp_cutoff mode_cutoff mean_boost alg Tag_merge_modes_eps mode_eps m_monomials; 
    global jump_method r fig_save check_rank xfine ufine u1fine ampS ampC eps_sinc;
    global pseudoInv TOL_PINV TOL_MASK mask_zero_amp x0 hw eps_Jmat regularize_Jmat;
    regularize_Jmat = 1;
    eps_Jmat  = 10^-15;
    eps_sinc  = 10^-15;
    pseudoInv = 0;
    mask_zero_amp = 0;
%     mask_zero_amp = 1; TOL_MASK  = 10^-100;

    fig_save = 0; check_rank = 0;
%     fig_u = '~/Dropbox/RapakaS_ASM/Figures/u_Nk30_n6_N64.jpg'; fig_save = 1;
%     fig_up = '~/Dropbox/RapakaS_ASM/Figures/up_Nk30_n6_N64.jpg'; fig_save = 1;
%     fig_u  = '~/Dropbox/RapakaS_ASM/Figures/u_log_n3_N32.jpg'; fig_save = 1;
%     fig_up = '~/Dropbox/RapakaS_ASM/Figures/up_log_n3_N32.jpg'; fig_save = 1;
%     fig_u  = '~/Dropbox/RapakaS_ASM/Figures/u_eax_n2_N64.jpg'; fig_save = 1;
%     fig_up = '~/Dropbox/RapakaS_ASM/Figures/up_eax_n2_N64.jpg'; fig_save = 1;
%     fig_u  = '~/Dropbox/RapakaS_ASM/Figures/u_x3_n2_N64.jpg'; fig_save = 1;
%     fig_up = '~/Dropbox/RapakaS_ASM/Figures/up_x3_n2_N64.jpg'; fig_save = 1;
%     fig_u  = '~/Dropbox/RapakaS_ASM/Figures/u_x_n1_N64.jpg'; fig_save = 1;
%     fig_up = '~/Dropbox/RapakaS_ASM/Figures/up_x_n1_N64.jpg'; fig_save = 1;
%     fig_u  = '~/Dropbox/RapakaS_ASM/Figures/u_gauss_n3_N64.jpg'; fig_save = 1;
%     fig_up = '~/Dropbox/RapakaS_ASM/Figures/up_guass_n3_N64.jpg'; fig_save = 1;

    %% ============== input parameters ================
    NN = [ 40 80 160 320  ];    % mesh size
    NN = 64;
    Nk = 30; %NN/8; % # of non-periodic modes in the multimode sinusoidal function (IC_profile=4)
    nc=2; ns=nc; % num. of non-periodic modes while decomposing a function into periodic and aperiodic functions
    nmax = max(nc,ns);
    n_apmodes = max(nc,ns);% num. of non-periodic modes while decomposing a function into periodic and aperiodic functions
    n_disc = 1;  % number of discontinuous locations in Eckhoff's method
    order_disc = nmax*4-1;  % order of discontinuities in Eckhoff's & Roache's method
    nmodes_Prony = NN/2;  % number of modes in Prony method
    
    jump_method = 1; % 1-analytical, 2 - FD, 3 - HFD, 4-spline, otherwise - high resolution schemes
    r = 6;      % order of accuracy for jump calculation via FD
    r_fd = 6;   % order of accuracy for derivative computation using FD method
    TOL_DetJ = 10^-14; check_rank = 0; % limits the number of modes based on rank of Jmat
    pseudoInv = 0; TOL_PINV=10^-15;  % use pseudoinverse to compute inverse of J (possibly rank deficit)
    xoff = .5; % test case 10: offset for u = log(x+pi+xoff) profile

    mean_boost = 0.95; alg = 2; 
    Tag_merge_modes_eps = 0; mode_eps = 0.1;
    amp_cutoff = 10^-8; % remove modes with amplitude below 'amp_cutoff'
    L=2*pi;    
    T= .75;   c=1; 
    % 9-eax, 10-log, 4-multimode, 5-monomials, 11-Gaussian
    IC_profile = 4;  m_monomials=3; filter_kN = 0; fig_Ref = 2; fig_Ref2 = 13; 
    time_DM = 1; CFL = 1; dt=10^-2; max_steps=2;
    % 1->FD_Spectral, 2->Linders et al., 3->Classic CD
    stencil = 2;  
    sw=4;   % stencil width (for stencil = 2 or 3)
    % sw=4; p=1; tau=6;     % Linder et al.       
    % sw=4; p=sw; tau=sw;   % classic CD schemes
%     scheme = 'FD';
    scheme = 'FD45';
%     scheme = 'FD113';
%     scheme = 'FFT';
    scheme = 'FFT45';
%     scheme = 'FFT45B';
%     scheme = 'FFTB';
    scheme = 'ASM45';
%     scheme = 'ASM45B';
    scheme = 'ASMB';  % solve Burger's eqn
    % filter_kN = 1;
    RTol = 10^-12; ATol = RTol; 
    Opts=odeset('RelTol',RTol,'AbsTol',ATol);
%     Opts=odeset('RelTol',RTol,'AbsTol',ATol,'InitialStep',dt,'MaxStep',dt);
    % odeset('RelTol',RTol,'AbsTol',ATol,'OutputFcn',@odeplot);
%     leg = [];
    %% ================================================
    fprintf('Numerical scheme : %s\n', scheme);
    fprintf('Temporal scheme  : RK%d (if custom)\n\n', time_DM);
    nMesh = size(NN,2);
    e_N = zeros(nMesh,2);
    e0_u_N = zeros(nMesh,2);
    e0_dudx_N = zeros(nMesh,2);
for nn=1:nMesh   % loop over meshes
    N=NN(nn); 
    mode_cutoff = 10;%floor(N/2)-1;
    dx=L/N; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% mesh generation and initial condition
    xe = linspace(0,L,N+1).';  % xe - column vector
    xe = linspace(-L/2,L/2,N+1).';  % xe - column vector
    xfine = linspace(-L/2,L/2,1000+1).';  % xe - column vector
    ufine = xfine;
    u1fine = xfine;
    
    x  = xe;  % Finite difference formulation
%     x  = xe(1:N);  % Finite difference formulation
    % x  = xe(1:N)+dx/2;  % Finite volume formuation
    xx = x-c*T;    % used for exact solution
    k = ([1.2 3.4 5.6 5.8 5.9 10.6 15.9 ]) * 2*pi/L;   % single mode
    k = ([1.2 2.4 3.3 4.7 5.4 6.5]) * 2*pi/L;   % 6.5 - saddle point
    k = linspace(0.5,Nk-0.5,Nk);
    k = linspace(0,Nk-1,Nk);
    eps = rand(0,Nk-1);
    eps = [
           0.521135830804001
           0.231594386708524
           0.488897743920167
           0.624060088173690
           0.679135540865748
           0.395515215668593
           0.367436648544477
           0.987982003161633
           0.037738866239552
           0.885168008202475
           0.913286827639239
           0.796183873585212
           0.098712278655574
           0.261871183870716
           0.335356839962797
           0.679727951377338
           0.136553137355370
           0.721227498581740
           0.106761861607241
           0.653757348668560
           0.494173936639270
           0.779051723231275
           0.715037078400694
           0.903720560556316
           0.890922504330789
           0.334163052737496
           0.698745832334794
           0.197809826685929
           0.030540946304637
           0.744074260367462
           0.500022435590201
           0.479922141146060    ].';
%     eps = round(eps*100)/100;
%     for ii=1:5
%         if((Nk/10^ii)<1)
%             Nd = 10^ii;
%             break;
%         end
%     end
%     Nd
%     eps = linspace(1,Nk,Nk)/Nd;
    eps = linspace(1,Nk-1,Nk)/(Nk);
    k = k+eps(1:Nk);
    for ii=1:Nk
        k(ii) = ii-1+(ii-1)/Nk*(Nk-2)/(Nk-1)+1/Nk;
    end
    modes = k.'
    ampS = ones(size(k));
    ampC = ones(size(k));


%     ampS = ([0.7 1]);
%     ampC = ([0 0]);
% %     k = ([5 12 ]) * 2*pi/L;  % single mode
%     k = ([5.3 12.4 ]) * 2*pi/L;  % single mode
%     k = ([3 4 6 8 9.3 10.6 12.7 ]) * 2*pi/L;  % single mode
%     k = ([1.67 2.3 3.3 4.2 5.8 6.7 7.3 8.5 ]) * 2*pi/L;           % single mode
%     k = ([  1.3675 ]) * 2*pi/L;   % multi-mode
%     k = ([0.4 1.65 2.3 3.6 4.8 5.5 6.4 7.3 8.7]) * 2*pi/L;   % multi-mode
%     k = ([1.65 2.3 5.6 8.4 9.3 12.5 13.7]) * 2*pi/L;   % multi-mode
%     k = ([1:N/2-1]+0.5) * 2*pi/L;   % multi-mode
    %%% INITIAL CONDITION %%%
    switch(IC_profile)
      case(1) % Gaussian
        x0 = 0;  hw = 1; xn = (x-x0)/hw;
        uM = sqrt(pi)/(2*pi)
        u0 = uM + exp(-((x-x0)/hw).^2);    % initial condition
        u1 = -2*(x-x0)/hw^2 .* exp(-((x-x0)/hw).^2);  % exact first derivative of u0
        u_ex = uM + exp(-((xx-x0)/hw).^2); % EXACT SOLUTION at t=T
        % solve non-linear eqn: u = sin(x-ut)
        for ii=1:1000
            u_ex = uM + exp(-((x-u_ex*T-x0)/hw).^2);   % EXACT SOLUTION at t=T
        end
%         xn = (x-u_ex*T-x0)/hw;
        % u0 = u_ex;
        % u1 = 2*xn.*exp(-(xn).^2)./( 2*T*xn.*exp(-xn.^2)-hw^2 );
        % % u0 = -u0; u1 = -u1; u_ex = -u_ex;
      case(2) % multi-mode sin wave
        u0 = sum(ampS.*sin(k.*x),2);      % initial condition
        u1 = sum(ampS.*k.*cos(k.*x),2);    % exact first derivative of u0
        u_ex = sum(ampS.*sin(k.*xx),2);   % EXACT SOLUTION at t=T
      case(3) % multi-mode cos wave
        u0 = sum(ampC.*cos(k.*x),2);      % initial condition
        u1 = sum(-ampC.*k.*sin(k.*x),2);    % exact first derivative of u0
        u_ex = sum(ampC.*cos(k.*xx),2);   % EXACT SOLUTION at t=T
      case(4) % multi-mode wave
        u0 = sum(ampS.*sin(k.*x)+ampC.*cos(k.*x),2);      % initial condition
        u1 = sum(ampS.*k.*cos(k.*x)-ampC.*k.*sin(k.*x),2);    % exact first derivative of u0
        u_ex = sum(ampS.*sin(k.*xx)+ampC.*cos(k.*xx),2);   % EXACT SOLUTION at t=T
        ufine = sum(ampS.*sin(k.*xfine)+ampC.*cos(k.*xfine),2);      % initial condition
        u1fine = sum(ampS.*k.*cos(k.*xfine)-ampC.*k.*sin(k.*xfine),2);    % exact first derivative of u0
      case(5) % saw-tooth wave / ramp / polynomial (y=x^m)
        m = m_monomials;
        u0 =   x.^(m  );      % initial condition
        u1 = m*x.^(m-1);    % exact first derivative of u0
        u_ex = xx.^(m);   % EXACT SOLUTION at t=T
        k = 2*pi/L * (0:N/4)+.1;
      case(6) % sin wave for Burger's Eqn
        % u0 = 0.25 + 0.5*sin(pi*x);      % initial condition
        pm = -1; % pm = +/- ==> u = pm*sin(x) 
        u0 = pm*sin(x);      % initial condition
        u1 = pm*cos(x);
        xx = x-u0*T;         % characteristic variable
        u_ex = pm*sin(xx);   % guess for solution at t=T
        % solve non-linear eqn: u = sin(x-ut)
        for ii=1:1000
            u_ex = pm*sin(x-u_ex*T);   % EXACT SOLUTION at t=T
        end
        % u0 = u_ex;
        % u1 = pm * cos(x-u_ex*T)./(1 + pm* T*cos(x-u_ex*T));
      case(7) % multi-mode wave
        u0 = exp(sin(65.5*x-27.7*pi)-cos(20.6*x));      % initial condition
        u1 = exp(sin(65.5*x-27.7*pi)-cos(20.6*x))*(65.5*cos(65.5*x-27.7*pi)+20.6*sin(20.6*x));    % exact first derivative of u0
        u_ex = exp(sin(65.5*xx-27.7*pi)-cos(20.6*xx));   % EXACT SOLUTION at t=T
      case(8) % f(x) = e^x-cosh(1)
        u0 = exp(x)-cosh(1);      % initial condition
        u1 = exp(x);    % exact first derivative of u0
        u_ex = exp(xx);   % EXACT SOLUTION at t=T
      case(9) % f(x) = e^(ax)sin(bx), 
        th = pi; a = -1/pi; b= 3/4;% * 2*pi;
        u0 = exp(a*(x+th)).*sin(b*(x+th));      % initial condition
        u1 = a*exp(a*(x+th)).*sin(b*(x+th)) + b*exp(a*(x+th)).*cos(b*(x+th));    % exact first derivative of u0
        u_ex = exp(a*(xx+th)).*sin(b*(xx+th));   % EXACT SOLUTION at t=T
      case(10) % f(x) = log(x+pi+1/2), 
        th = pi+xoff;
        u0 = log(x+th);      % initial condition
        u1 = 1./abs(x+th);    % exact first derivative of u0
        u_ex = log(xx+th);   % EXACT SOLUTION at t=T
      case(11) % f(x) = log(x+pi+1/2), 
        x0 = pi*3/4;
        hw  = 1;
        u0 = exp(-((x-x0)/hw).^2);      % initial condition
        u1 = exp(-((x-x0)/hw).^2).*(-2*(x-x0)/hw^2);    % exact first derivative of u0
        u_ex = exp(-(xx-x0).^2/hw);   % EXACT SOLUTION at t=T
    end
% k.'
    % Remove 2*dx waves from the initial and exact solutions
    if(filter_kN)
        u0 = fft(ifftshift(u0));   % FFT: ifftshift rearranges x=0:L-dx to x=-L/2:L/2-dx
        u0(N/2+1)=0;    % nullify odd-ball wave number    
        u0 = real(fftshift(ifft(u0)));  % iFFT

        u1 = fft(ifftshift(u1));   % FFT: ifftshift rearranges x=0:L-dx to x=-L/2:L/2-dx
        u1(N/2+1)=0;    % nullify odd-ball wave number    
        u1 = real(fftshift(ifft(u1)));  % iFFT

        u_ex = fft(ifftshift(u_ex));   % FFT: ifftshift rearranges x=0:L-dx to x=-L/2:L/2-dx
        u_ex(N/2+1)=0;    % nullify odd-ball wave number    
        u_ex = real(fftshift(ifft(u_ex)));  % iFFT
    end
    
    %% ============== stencil weights ================
    % sw=1; wts = [-1 0 1]/2; % must be a row vector
    % sw=2; wts = [1 -8 0 8 -1]/12; % must be a row vector
    % sw=3; wts = [-1 9 -45 0 45 -9 1]/60; % must be a row vector
    % ==============================
    switch(stencil)
        case(1) % Fast Finite Difference schemes
            [sw,wts] = fd_stencil(N);
            l1=strcat('FDS');
        case(2) % Linders et al., JCP, (2020) schemes            
            % sw: half stencil width, 2p: accuracy, tau: HSW of reference scheme
            p=1; tau=6;       
            wts = optimalStencil( p, sw, tau, u0 ).';  % n=sw 
            wts = [-flip(wts),0,wts];
            l1=strcat('Lo',int2str(2*sw+1),int2str(2*p));
            % l1=strcat('LO',int2str(2*sw+1),int2str(2*p));
        case(3) % classic CD schemes
            % sw: half stencil width, 2p: accuracy, tau: HSW of reference scheme
            p=sw; tau=sw;     % classic CD schemes
            wts = optimalStencil( p, sw, tau, u0 ).';  % n=sw 
            wts = [-flip(wts),0,wts];
            l1=strcat('c',int2str(2*sw+1),int2str(2*p));
            l1=strcat('CD',int2str(2*p));
    end
    % wts.'
    % ==============================
    %% spatial discretization error test
    % FDM
    % ux = dudx(u0,dx,sw,wts);  % periodic domains
    ux = dudx_FD(u0,r_fd,dx); % (r_fd)^th order finite difference scheme with one-sided stencils near boundary
    e_dudxFD = abs(ux-u1);
    error_dudxFD = max(e_dudxFD);

    % FFT
    ux = dudx_FFT(u0,L,1);
    e_dudxFFT = abs(ux-u1);
    error_dudxFFT = max(e_dudxFFT);
    ux_FFT = ux;

    % Roache
    ux_Roache = u1;
    ux_Roache = dudx_Roache(x,u0,order_disc,1,ux_FFT);
    e_dudx_Roache = abs(ux_Roache-u1);
    error_dudx_Roache = max(e_dudx_Roache);

    % Eckhoff
    ux_Eckhoff = u1;
    uj_E = compute_derivative_jumps(u0,order_disc+1);
    [ux_Eckhoff,uxp_E,uxa_E,up_E,ua_E] = ddx_Eckhoff(x,u0,uj_E,order_disc,n_disc,u1);
    e_dudx_Eckhoff = abs(ux_Eckhoff-u1);
    error_dudx_Eckhoff = max(e_dudx_Eckhoff);

    % Prony
    ux_Prony = u1;
    [ux_Prony,u_Prony] = ddx_Prony(x-x(1),u0,nmodes_Prony,u1);
    skip_nodes = 1;%-1;
    ux_Prony(end-skip_nodes:end) = u1(end-skip_nodes:end);
    e_dudx_Prony = abs(ux_Prony-u1);
    error_dudx_Prony = max(e_dudx_Prony);

    ux_compare = cat(2,ux_FFT,ux_Roache,ux_Eckhoff,ux_Prony);
    method = ["FFT";"Roache";"Eckhoff";"Prony"];

%     ux_compare = ux_FFT; method = 'FFT';
%     ux_compare = ux_Roache; method = 'Roache';
%     ux_compare = ux_Eckhoff; method = 'Eckhoff';
%     ux_compare = ux_Prony; method = 'Prony';
    
    % GFS
    [ux,kn] = dudx_GFS(u0,1,ux_compare,method);    % find modes in 'u0' for ASM
    e_GFS = abs(ux-u1);
    error_dudx_GFS = max(e_GFS);
    % kn = find_modes(u0,1);    % find modes in 'u0' for GFS
    % kn = find_modes(u0.^2/2,L,1);    % find modes in 'u0' for GFS
    % ux = dudx_ASM(u0,L,kn,0);
    
    % modified GFS: periodic part is forced to be zero at near-boundary nodes
    % % ux = dudx_modASM_poly(u0,1);    % find modes in 'u0' for GFS
    % % ux = dudx_modASM_randModes(u0,1);    % find modes in 'u0' for GFS
    % ux = dudx_modASM(u0,1);    % find modes in 'u0' for GFS    
    % e_modGFS = abs(ux-u1);
    % error_dudx_modGFS = max(e_modGFS);

    fprintf('||e_GFS||_1   = %.2e\n',norm(e_GFS,1)*dx);
    fprintf('||e_GFS||_2   = %.2e\n',norm(e_GFS,2)*dx^(1/2));
    fprintf('||e_GFS||_Inf = %.2e\n\n',norm(e_GFS,Inf));
    fprintf("N=%d:\t max.|e'| \t ||e'||_2 \n",N);
    fprintf('e_FD  =\t %.2e    &\t %.2e\n',error_dudxFD,norm(e_dudxFD,2)*dx^(1/2));
    fprintf('e_FFT =\t %.2e    &\t %.2e\n',error_dudxFFT,norm(e_dudxFFT,2)*dx^(1/2));
    fprintf('e_GFS =\t %.2e    &\t %.2e\n',error_dudx_GFS,norm(e_GFS,2)*dx^(1/2));
    % fprintf('e_mGFS=\t %.2e    &\t %.2e\n',error_dudx_modGFS,norm(e_modGFS,2)*dx^(1/2));
    fprintf('e_Roa=\t %.2e    &\t %.2e\n',error_dudx_Roache,norm(e_dudx_Roache,2)*dx^(1/2));
    fprintf('e_Eck=\t %.2e    &\t %.2e\n',error_dudx_Eckhoff,norm(e_dudx_Eckhoff,2)*dx^(1/2));
    fprintf('e_Pro=\t %.2e    &\t %.2e\n\n',error_dudx_Prony,norm(e_dudx_Prony,2)*dx^(1/2));
%     fprintf('||e_mGFS||_1   = %.2e\n',norm(e_modGFS,1)*dx);
%     fprintf('||e_mGFS||_2   = %.2e\n',norm(e_modGFS,2)*dx^(1/2));
%     fprintf('||e_mGFS||_Inf = %.2e\n',norm(e_modGFS,Inf));
    e0_dudx_N(nn,1) = error_dudx_GFS;
    e0_dudx_N(nn,2) = error_dudxFFT;  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('%.2e & %.2e & %.2e &  %.2e & %.2e &  %.2e & %.2e & %.2e & %.2e & %.2e & %.2e & %.2e\n\n',error_dudx_GFS,norm(e_GFS,2)*dx^(1/2),error_dudxFD,norm(e_dudxFD,2)*dx^(1/2),error_dudxFFT,norm(e_dudxFFT,2)*dx^(1/2),error_dudx_Roache,norm(e_dudx_Roache,2)*dx^(1/2),error_dudx_Eckhoff,norm(e_dudx_Eckhoff,2)*dx^(1/2),error_dudx_Prony,norm(e_dudx_Prony,2)*dx^(1/2));


  figure(fig_Ref);
  figure(fig_Ref2);
%   figure(44);
%   break;
  return;
    %% time marching
    switch(scheme)
        % custom ODE integration: supports Euler, RK3, RK4
        % use 'ode_solve2()' for full spectral difference in space
        case('FD')  % FD schemes + Euler/RK3/RK4 
        [t,u] = ode_solve(@(u)-c*dudx(u,dx,sw,wts),T,u0,time_DM,dx,CFL,max_steps);
        case('FFT') % FFT + Euler/RK3/RK4 
        [t,u] = ode_solve(@(u)-c*dudx_FFT(u,L,0),T,u0,time_DM,dx,CFL,max_steps);  l1 = 'FFT';
        case('ASM') % FFT + Euler/RK3/RK4
        [t,u] = ode_solve(@(u)-c*dudx_ASM(u,L,kn,0),T,u0,time_DM,dx,CFL,max_steps);  l1 = 'ASM';

        case('FFTB') % FFT + Euler/RK3/RK4 - Burger's eqn
        [t,u] = ode_solve(@(u)-dudx_FFT(u.^2/2,L,0),T,u0,time_DM,dx,CFL,max_steps);  l1 = 'FFT';
        case('ASMB') % FFT + Euler/RK3/RK4 - Burger's eqn
        % [t,u] = ode_solve(@(u)-dudx_ASM(u.^2/2,L,kn,0),T,u0,time_DM,dx,CFL,max_steps);  l1 = 'ASM';        
        [t,u] = ode_solve2(@(uk,kn)uux_ASM(uk,kn),T,u0,time_DM,dx,CFL,max_steps);  l1 = 'ASM';

        % MATLAB in-built ODE integration
        case('FD45')
        [t,ut] = ode45(@(t,u)-c*dudx(u,dx,sw,wts),[0 T],u0,Opts); u = ut(end,:).';
        case('FD113')
        [t,ut] = ode113(@(t,u)-c*dudx(u,dx,sw,wts),[0 T],u0,Opts); u = ut(end,:).';
        case('FFT45')
        [t,ut] = ode45(@(t,u)-c*dudx_FFT(u,L,0),[0 T],u0,Opts); u = ut(end,:).';  l1 = 'FFT';
        case('ASM45')
        [t,ut] = ode45(@(t,u)-c*dudx_ASM(u,L,kn,0),[0 T],u0,Opts); u = ut(end,:).';  l1 = 'ASM';

        case('FFT45B')  % Burgers eqn
        [t,ut] = ode45(@(t,u)-dudx_FFT(u.^2/2,L,0),[0 T],u0,Opts); u = ut(end,:).'; l1 = 'FFT';
        case('ASM45B')  % Burgers eqn
        [t,ut] = ode45(@(t,u)-dudx_ASM(u.^2/2,L,kn,0),[0 T],u0,Opts); u = ut(end,:).'; l1 = 'ASM';
    end
    dtv = [t(1); t(2:end)-t(1:end-1)];

    error = u-u_ex;
    % e_N(nn,1) = max(abs(error(1:N)));
    e_N(nn,1) = max(abs(error));
%     e_N(nn,2) = error_dudxFFT;
    fprintf('nsteps = %d\nt_end  = %.2f\n',numel(t),t(end));
    fprintf('dt_min = %.2e\ndt_max = %.2e\n',min(dtv(2:end)),max(dtv));
    fprintf('-------------------\n');
    fprintf('max|e| = %.2e\n',max(abs(error)));
    fprintf('-------------------\n\n');
    plot_sol(u);
end
    figure(3); %clf; 
    loglog(NN,e_N(:,1),'-o','DisplayName',l1); hold on;
    xlabel('N'); ylabel('max |e|');
    grid on; grid minor;
    axis square;
%     legend('FD','FFT');
    set(gca,'FontSize',20);
%     set(gca,'xscale','log','yscale','log')
%     legend({'FFT','o92','c98','FFT, $\Delta t=2\times10^{-5}$','o72','c76','o52','c54','c32'},'Interpreter','Latex')
%     legend({'FDS','FFT','o92','c98','o72','c76','o52','c54','c32'},'Interpreter','Latex')
%     legend({'FDS, $\Delta t\sim 1.6\times 10^{-4}$','FFT, $\Delta t\sim 1.6\times 10^{-4}$','FFT, $\Delta t\sim 5.8\times 10^{-5}$','FDS, $\Delta t\sim 5.8\times 10^{-5}$','o72','c76','o52','c54','c32'},'Interpreter','Latex')
    lgd=legend('Location','NorthWest');
%     lgd.NumColumns=1;
    lgd.FontSize = 14;
    % lgd.Title.String = 'Scheme';    
    set(findall(gca, 'Type', 'Line'),'LineWidth',1);
%     saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/error_ode45_4mode.jpg');
%     saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/ode45_Gaussian4.jpg');

    figure(31); %clf; 
    loglog(NN,e0_dudx_N(:,1),'-o','DisplayName','ASM'); hold on;
    loglog(NN,e0_dudx_N(:,2),'-^','DisplayName','FFT'); hold on;
    xlabel('N'); ylabel('max |e|');
    grid on; grid minor;    axis square;
    set(gca,'FontSize',20);
    lgd=legend('Location','NorthWest');
    lgd.FontSize = 14;
    lgd.Title.String = 'e in du/dx';
    set(findall(gca, 'Type', 'Line'),'LineWidth',1);
%     saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/error_dudx_4mode_t0.jpg');

    figure(32); %clf; 
    loglog(NN,e0_u_N(:,1),'-o','DisplayName','ASM'); hold on;
    loglog(NN,e0_u_N(:,2),'-^','DisplayName','FFT'); hold on;
    xlabel('N'); ylabel('max |e|');
    grid on; grid minor;    axis square;
    set(gca,'FontSize',20);
    lgd=legend('Location','NorthWest');
    lgd.FontSize = 14;
    lgd.Title.String = 'e in u';
%     lgd.Title.String = 'e in du/dx';
    set(findall(gca, 'Type', 'Line'),'LineWidth',1);
%     saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/error_u_4mode_t0.jpg');

% kn = find_modes(u.^2/2,1);
figure(fig_Ref);  
figure(fig_Ref2);
return;

    %% FUNCTIONS USED IN THE MAIN
    function [uk,kn] = amp_filter(uk,kn,amp_cutoff)
            ind = abs(uk) >= amp_cutoff;
            uk = uk(ind);
            kn = kn(ind);        
    end

    function [uk,kn] = mode_filter(uk,kn,mode_cutoff)
            ind = abs(kn) <= mode_cutoff;
            uk = uk(ind);
            kn = kn(ind);        
    end

    function [uk,kn] = sort_amp(uk,kn)
            [uk,I] = sort(uk,'descend','ComparisonMethod','abs');
            kn = kn(I);
    end

    function [uk,kn] = sort_modes(uk,kn)
            [kn,I] = sort(kn,'ascend','ComparisonMethod','auto');
            uk = uk(I);
    end

    function [uk,kn] = sort_modes_abs(uk,kn)
            [kn,I] = sort(kn,'ascend','ComparisonMethod','abs');
            uk = uk(I);
    end

    function [tn,u] = ode_solve2(func,T,u0,time_DM,dx,CFL,max_steps)
        %CFL = 1; %max_steps = 2000;
        % max_steps = 1;
        % M = 2*floor(N/2);

        isperiodic = false;        
        kn = find_modes(u0,0);    % find modes in 'u0' for ASM
        uk = nmFFT(u0,kn,isperiodic);   % FFT through non-integer basis

%         [uk,I] = sort(uk,'descend','ComparisonMethod','abs');
%         kn = kn(I);
        fprintf('  modes\t amplitudes\n');
        for i=1:numel(kn)
            fprintf('%8.4f  %8.2e\n', kn(i),abs(uk(i)));
        end
        % [uux_k, kn_uu] = uux_ASM(uk,kn);

        dt = abs(CFL*dx/max(abs(u0)));
        t=0; n=0; tn = zeros(max_steps,1);
        while(t<=T && n<max_steps)
            % dt = abs(CFL*dx/max(abs(u)));
            dt = min(dt,10^-4);
            dt = 0.01;
            n=n+1;        
            t=t+dt;
            tn(n) = t;
            fprintf('n = %d\t dt = %.1e  t = %.4f\n', n, dt, t);
            switch(time_DM)
                case(1) % Euler explicit time marching scheme
                    % un = u + dt*func(u);
                    [uux_k,kn_uu] = func(uk,kn);
                    [uk_n,kn_n] = merge_modes([uk;dt*uux_k],[kn;kn_uu]);

                case(3) % RK3 explicit time marching scheme
                    % RK step 1
                    % rh = func(u);
                    % un = u + 8/15*dt*rh;
                    [rh,kn1] = func(uk,kn);
                    [uk_n,kn_n] = merge_modes([uk;8/15*dt*rh],[kn;kn1]);

                    % RK step 2
                    % rh2 = func(un);
                    % un = u + dt*( 1/4*rh + 5/12*rh2 );
                    [rh2,kn2] = func(uk_n,kn_n);
                    [uk_n,kn_n] = merge_modes([uk;1/4*dt*rh;5/12*dt*rh2],[kn;kn1;kn2]);

                    % RK step 3
                    % rh2 = func(un);
                    % un = u + dt*( 1/4*rh + 3/4*rh2 );
                    [rh2,kn2] = func(uk_n,kn_n);
                    [uk_n,kn_n] = merge_modes([uk;1/4*dt*rh;3/4*dt*rh2],[kn;kn1;kn2]);
                    

                case(4) % RK4 explicit time marching scheme
                    % RK step 1
                    % rh = func(u);
                    % un = u + 1/2*dt*rh;
                    % utemp = -u + un;
                    [rh,kn1] = func(uk,kn);
                    [uk_n,kn_n] = merge_modes([uk;1/2*dt*rh],[kn;kn1]);
                    [uk_t,kn_t] = merge_modes([-uk;uk_n],[kn;kn_n]);

                    % RK step 2
                    % rh = func(un);
                    % un = u + 1/2*dt*rh;
                    % utemp = utemp + 2*un;
                    [rh,kn1] = func(uk_n,kn_n);
                    [uk_n,kn_n] = merge_modes([uk;1/2*dt*rh],[kn;kn1]);
                    [uk_t,kn_t] = merge_modes([uk_t;2*uk_n],[kn_t;kn_n]);

                    % RK step 3
                    % rh = func(un);
                    % un = u + dt*rh;
                    % utemp = utemp + un;
                    [rh,kn1] = func(uk_n,kn_n);
                    [uk_n,kn_n] = merge_modes([uk;dt*rh],[kn;kn1]);
                    [uk_t,kn_t] = merge_modes([uk_t;uk_n],[kn_t;kn_n]);

                    % RK step 4
                    % rh = func(un);
                    % un = 1/3*utemp + 1/6*dt*rh;
                    [rh,kn1] = func(uk_n,kn_n);
                    [uk_n,kn_n] = merge_modes([1/3*uk_t;1/6*dt*rh],[kn_t;kn1]);
            end                        
            uk = uk_n;
            kn = kn_n;

            [uk,kn] = amp_filter(uk,kn,amp_cutoff);
%             [uk,kn] = mode_filter(uk,kn,mode_cutoff);
            [uk,kn] = mode_redistribute(uk,kn,x);
            
            % [uk,kn] = sort_amp(uk,kn);
            [uk,kn] = sort_modes(uk,kn);
            if(Tag_merge_modes_eps)
              [uk,kn] = merge_modes_eps(uk,kn,mode_eps);
            end
        end
        fprintf('  modes\t amplitudes\n');
        for i=1:numel(kn)
            fprintf('%8.4f  %8.2e\n', kn(i),abs(uk(i)));
        end
        u = nmiFFT(uk,kn,isperiodic);  % iFFT
        u = real(u);
        tn(n+1:end) = [];
    end

    function pseudo_inverse = my_pinv(A,eps)        
        % Perform Singular Value Decomposition
        [U, S, V] = svd(A);
        
        % Invert the singular values in S
        S_inv = zeros(size(S')); % Transpose of S to match dimensions for pseudo-inverse
        rank_full = min(size(S));
        rank = 0;
        for i = 1:min(size(S))
            if S(i, i) > eps  % Check if singular value is non-zero (avoid division by zero)
                S_inv(i, i) = 1 / S(i, i);
                rank = rank + 1;
            else
                break;
            end
        end
        fprintf("rank : %d / %d\n", rank, rank_full);
        
        % Compute pseudo-inverse
        pseudo_inverse = V * S_inv * U';
%         pseudo_inverse = V * inv(S) * U';
%         pseudo_inverse = inv(A);
    end

    %% non-modal FFT routines
    function uk = nmFFT(u,kn,isperiodic)   % kn - row vec
        if(size(kn,1)>1)    % if kn is a col vector
            kn = kn.';      % make sure kn is a row vector
        end        
%         km = kn.';  % col vector
%         H = (kn-km);    % km-col vec, km-row vec
%         H = arrayfun(@(y)sum(exp(1i*y*x)),H);%,'UniformOutput',false);
%         % E = arrayfun(@(y)exp(-1i*y),km*x.');
%         E = exp(-1i*km*x.'); % km-col vec, x.' - row vec
%         g = E*u;

        if(isperiodic)
            E = exp(1i*kn.*x(1:N)); % periodic case
            g = E'*u(1:N);  % periodic case
        else
            % aperiodic case: i=1:N+1 (includes both end points +/- pi)
            E = exp(1i*kn.*x); % kn->row vec, x->col vec
            % normalize the basis vectors after including 'N+1' node
            % E(1,:) = E(1,:)*sqrt(.5);
            % E(end,:) = E(end,:)*sqrt(.5);
            g = E'*u;
        end

        H = E'*E;        
%         cond_H = cond(H)
%         H/N
%         iH = inv(H)*N
%         max_iH = max(max(iH))
%         S = E*inv(H)*E';
%         for i=1:4
%             norm(S(:,i),2)
%         end
%         S2 = E'*inv(E*E')*E;
        
        uk = H\g;        
    end

    function u = nmiFFT(uk,kn,isperiodic)
        % if(size(kn,1)>1)    % if kn is a col vector
        if(iscolumn(kn))    % if kn is a col vector
            kn = kn.';      % make sure kn is a row vector
        end
        
        if(isperiodic)
            E = exp(1i*kn.*x(1:N)); % periodic case
            u = E*uk;
            u = [u;u(1)];   % periodic case: x(1:N)
        else
            % aperiodic case: i=1:N+1 (includes both end points +/- pi)
            E = exp(1i*kn.*x); % kn->row vec, x->col vec
            % normalize the basis vectors after including 'N+1' node
            % E(1,:) = E(1,:)*sqrt(.5);
            % E(end,:) = E(end,:)*sqrt(.5);
            u = E*uk;
        end
    end

    function e = error_nmFFT(u,kn)
        % if(size(kn,1)>1)    % if kn is a col vector
        if(iscolumn(kn))    % if kn is a col vector
            kn = kn.';      % make sure kn is a row vector
        end

%         fprintf('kn= %.4f\n', kn);
        kn = [-flip(kn),0,kn];
        % aperiodic case: i=1:N+1 (includes both end points +/- pi)
        E = exp(1i*kn.*x); % kn->row vec, x->col vec
        % normalize the basis vectors after including 'N+1' node
        % E(1,:) = E(1,:)*sqrt(.5);
        % E(end,:) = E(end,:)*sqrt(.5);

        H = E'*E;        
%         M = E*inv(H)*E';
% %         eig_H = eig(H)
% %         eig_M = eig(M)
% %         eig_I_M = eig(eye(size(M))-M)
%         e = u'*u - u'*M*u;
        e = u'*u - u'*(E*(inv(H)*(E'*u)));
        e = real(e);
%         error = [error,e];
    end

    function Grad_E2 = grad_e(kl,x,u)
        % both x and u must be col vectors
        Grad_E2 = u'* sin(kl*(x-x.')) *(x.*u);
    end

    %% Adaptive spectral method: find modes in a signal
    % spectral derivative through non-integer modes
    function kn = find_modes(u,plotTag)
        plotTag=1;
        M = floor(N/2)*2;
%         kn = 2*pi/L * (-M/2:M/2-1).'; % modes
        dim  = numel(k); 
%         dim  = M/2-1;
%         M    = 2*dim;
        base = 0.25;
                
%         display('Computing error matrix');
        mode_dim = 1; %numel(k);
        if(numel(k)==2) mode_dim = 2; end
        plot_error = 1;
%         tic;
        if(plot_error)
            if(mode_dim==1)
                dk = 0.01;
                % single-mode evaluation
%                 kvec = [0:dk:min(2*k,N/2)].';
%                 kvec = [k-4:dk:k+4].';
%                 kvec = [dk:dk:N/2-dk].';
                kvec = [dk:dk:N/2-.25].';
                Nk = size(kvec,1);
                error = zeros(Nk,1);
                for i=1:Nk
                    kn = kvec(i);
                    error(i) = (error_nmFFT(u,kn));
                end
                
                e_max = max(error);
                e_min = min(error);
                e_mean = mean(error);
                e_std  = std(error);
                e_skew = skewness(error);
                e_kurt = kurtosis(error);
                delta = e_max-e_min;
                fprintf('\nerror stat:\n----------- \nmax  = %.1f \nmin  = %.1f\n', e_max,e_min);
                fprintf('mean = %5.1f std  = %5.1f \nskew = %5.1f kurt = %5.1f \n',e_mean,e_std,e_skew,e_kurt);
                fprintf('range= %5.1f range/std = %.1f\n',delta,delta/e_std);
%                 e_fluc = error - e_mean;
%                 e_fluc = e_max - error;
%                 max_fluc = max(e_fluc)
%                 e_cutoff = std_cutoff_frac*e_std;
                e_cutoff = e_mean;
                e_cutoff = e_mean + mean_boost*(e_max-e_mean);
            if(alg==1)
                % Alg 1: guess initial values of kn
                % cut-off error above standard deviation
%                 ind = find(abs(e_fluc) > e_cutoff);
                ind = find(error < e_cutoff);
                kl = kvec( ind );
                el = error( ind ).';
                Nind = numel(ind);
                Nmodes = dim; %numel(k);
                dN = ceil(Nind/Nmodes);
                dN2 = ceil(dN/2);
                kn0 = zeros(Nmodes,1);
                for i=1:Nmodes
                    kn0(i) = kl(dN2 + (i-1)*dN);
                end
            else
                % Alg 2: guess initial values of kn
                % find all the local minimum & sort by error magnitude
                mode = 0;
                i = 1;
                    if(error(i)<=error(i+1))
                        if( error(i) < e_cutoff )
                            mode = mode + 1;
                            kn0(mode) = kvec(i);
                            local_min(mode) = error(i);
                        end
                    end
                
                for i=2:Nk-1
                    if(error(i)<error(i-1) && error(i)<=error(i+1))
                        % if( abs(e_fluc(i)) > e_cutoff )
                        if( error(i) < e_cutoff )
                            mode = mode + 1;
                            kn0(mode) = kvec(i);
                            local_min(mode) = error(i);
                        end
                    end
                end
                i=Nk;
                    if(error(i)<error(i-1))
                        if( error(i) < e_cutoff )
                            mode = mode + 1;
                            kn0(mode) = kvec(i);
                            local_min(mode) = error(i);
                        end
                    end
%                 fprintf('No of local min = %d\n',mode);
                % local_Min = local_min.'
                kn0 = kn0(1:end-0);
                kn0 = kn0.';
            end  
            kn0
            elseif(mode_dim==2)
                dk = 0.1;
                % two-mode evaluation
%                 kvec1 = [dk:dk:min(2*k(1),N/2)].';
%                 kvec2 = [dk:dk:min(2*k(2),N/2)].';
                kvec1 = [dk:dk:N/2-dk].';
                kvec2 = [dk:dk:N/2-dk].';
%                 kvec1 = [k(1)-3:dk:k(1)+3].';
%                 kvec2 = [k(2)-3:dk:k(2)+3].';
                Nk1 = size(kvec1,1);
                Nk2 = size(kvec2,1);
                error = zeros(Nk2,Nk1);
                kn = zeros(2,1);
                for i=1:Nk1
                    kn(1) = kvec1(i);
                for j=1:Nk2
                    kn(2) = kvec2(j);
                    % kn = [kvec1(i),kvec2(j)].';
                    if(i==j)
                        error(j,i) = 0;% 100;
                    else
                        error(j,i) = (error_nmFFT(u,kn));
                    end
                end
                end
                e_max = max(max(abs(error)))
                e_min = min(min(abs(error)))
                for i=1:Nk1
                for j=1:Nk2
                    if(i<=j) error(j,i) = e_max;    end
                end
                end
                [M,I] = min(abs(error));
                [M,J] = min(M);
                I = I(J);
                kn(1) = kvec1(I);
                kn(2) = kvec2(J);
                kn0 = kn
            end
        end
%         toc;
%         kn = 2*pi/L * ([1:dim].'+base);
        kn = kn0;
%         kn = 2*pi/L * ([1.57 2.4 3.3 3.6 4.5]).'; % modes
%         kn = 2*pi/L * ([-(M/2+base-1):0, 0, base:M/2]).'; % modes
%         e = error_nmFFT(u,kn);
        TOL = 1e-15;
        MaxIter = 2000;
%         kl = roots(@(kl)grad_e(kl,x,u), [1,10])
%         kl = roots(chebfun(@(kl)grad_e(kl,x,u), [0,19]))

%         options = optimset('TolFun',TOL,'TolX',TOL,'MaxFunEvals',MaxIter,'MaxIter',MaxIter);
% %         options = optimset(options,'Display','iter-detailed');%'notify','iter','final'
%         kn = sort(fminsearch(@(kn)error_nmFFT(u,kn),kn,options));

        options = optimoptions(@fminunc,'Display','none','FunctionTolerance',TOL,'MaxIterations',MaxIter,'SpecifyObjectiveGradient',false);
%         options = optimoptions(@fminunc,'Display','iter-detailed','FunctionTolerance',TOL,'MaxIterations',MaxIter,'SpecifyObjectiveGradient',false);
%         options = optimoptions(options,'Algorithm','trust-region');
        options = optimoptions(options,'Algorithm','quasi-newton');
        options = optimoptions(options,'OptimalityTolerance',TOL,'StepTolerance',TOL);
        options = optimoptions(options,'MaxFunctionEvaluations',10000);
%         options = optimoptions(options,'Display','iter-detailed');
        kn = sort(fminunc(@(kn)error_nmFFT(u,kn),kn,options));
        e_global = error_nmFFT(u,kn);
        kn 
%         step = 1;
%         cnt = 0;
%         tic;
%         % fprintf('kn  = %.4f, e_global = %.4f\n', kn,e_global);
%         while(e_global > 10^-6 && cnt<30)
%           fprintf('kn  = %.4f, e_global = %.4f\n', kn(1),e_global);
%           fprintf('kn  = %.4f, e_global = %.4f\n', kn(2),e_global);
%             cnt = cnt+1;
%             kn1 = min(kn+step,N/2-1);
%             kn1 = sort(fminunc(@(kn1)error_nmFFT(u,kn1),kn1,options));
%             e_new = error_nmFFT(u,kn1);
%             if(e_new < e_global)
%                 kn = kn1;
%                 e_global = e_new;
%             else
%                 step = -step;
%             end
%           fprintf('kn1 = %.4f, e_new    = %.4f\n', kn1(1),e_new);
%           fprintf('kn1 = %.4f, e_new    = %.4f\n', kn1(2),e_new);
%         end
%         toc;

%         error
%         options
% kn = k.';
% e_kn = kn-k.'
        kn = [-flip(kn);0;kn];
%         kn = 2*pi/L * (-M/2:M/2-1).'; % modes
%         kn = 2*pi/L * ([-M/2+0.5:0, 0, 0.5:M/2]).'; % modes

        % use this for even derivatives
        % kn = 2*pi/L * [0:M/2-1,-M/2,-M/2+1:-1].'; 
        % use this for odd derivatives: nullify odd-ball wave number
        % kn = 2*pi/L * [0:M/2-1, 0  ,-M/2+1:-1].'; 

%         mzero = M/2+1;
%         kn(mzero+19) =  2*pi/L * 19.75;
%         kn(mzero-19) = -2*pi/L * 19.75;
% %         kr = round(k).';
%         kr = fix(k).';
%         kn(mzero+kr) = k *2*pi/L;
%         kn(mzero-kr) =-k *2*pi/L;
% modes = kn(mzero:end)
        %mode_range = [kn(1); kn(end)]
%         fprintf('\nmode range = [%.2f : %.2f]\n', kn(1),kn(end));

        isperiodic = true;
        isperiodic = false;
%         fprintf('isperiodic = %d\n', isperiodic);
        uk = nmFFT(u,kn,isperiodic);   % FFT through non-integer basis
        [uk,kn] = amp_filter(uk,kn,amp_cutoff);
        % uk(M/2+1)=0;    % nullify odd-ball wave number    
        ux = nmiFFT(1i*kn.*uk,kn,isperiodic);  % iFFT
        ux = real(ux);
        % copy of the original function via iFFT
        uF = nmiFFT(uk,kn,isperiodic);  % iFFT
        uF = real(uF);
        eF = uF-u;

        if(plotTag)
            if(plot_error)
                if(mode_dim==1)
                    figure(2); clf;
                    semilogy(kvec,error,'r.-'); grid; hold on;
                    semilogy(kvec,e_mean*ones(size(kvec)),'--'); 
%                     semilogy(kvec,abs(e_fluc),':'); 
%                     semilogy(kvec,(e_cutoff)*ones(size(kvec)),':'); 
%                     semilogy(kvec,(e_mean-e_std)*ones(size(kvec)),':'); 
                    set(gca,'FontSize',24);
                    xlabel('k_n');
                    ylabel('e(k_n)');
                    axis square tight;
    %                 ylim([1,100]);
% saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/error_k1_5.3.jpg');
                elseif(mode_dim==2)
                    figure(2); clf;
                    surf(kvec1,kvec2,error); grid;
%                     semilogy(kvec1,error(2,:),'r.-'); hold on;
%                     semilogy(kvec2,error(:,40),'b.-'); grid;
                    set(gca,'FontSize',24);
                    xlabel('k_{1}');
                    ylabel('k_{2}');
                    zlabel('e(k_n)');
                    axis square tight; shading FLAT; %FACETED;
                    box on; grid on;
                    colorbar('LOCATION','NorthOutside');
                    view(3)
% saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/error_k1_5.3_k2_12.7.jpg');
                end
            end
            
            lw=1;
            e0_u_N(nn,1) = max(abs(eF));
            fprintf('max|u-u_ASM|: %.2e\n\n',max(abs(eF)));
%             fprintf("Nyq_comp: \nu : |  uk| = %10.2e \nu': |k*uk| = %10.2e\n\n",abs(uk(1)),abs(kn(1)*uk(1)));
            
            figure(13); clf;    
            subplot(3,1,1);
            plot((kn),real((uk)),'.-','DisplayName','real'); hold on;
            plot((kn),imag((uk)),'-','DisplayName','imag');
%             max_uk_im = max(abs(imag(uk)))
            xlabel('$k_n$','Interpreter','Latex'); ylabel('$\hat{u}_n$','Interpreter','Latex');
            legend('Location','NorthEast'); grid on; grid minor;
            set(gca,'FontSize',16); axis tight;
            xlim([-N/2,N/2]);
            set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
            % set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 1);

            subplot(3,1,2);
            plot(x,u,'ro','DisplayName','exact'); hold on;
            plot(x,real(uF),'b-','DisplayName','ASM');
            plot(x,eF,'k:','DisplayName','error');
            %plot(x,imag(uF),'--','DisplayName','error_im');
            xlabel('x');ylabel('u(x)'); grid on; grid minor;
            set(gca,'FontSize',16);
            set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
            % legend('exact','nm\_spec','nm\_error');  
            axis tight; %legend('Location','Best');
            
            subplot(3,1,3);
            plot(x,u1,'ro','DisplayName','exact'); hold on;
            plot(x,real(ux),'b-','DisplayName','ASM');
            plot(x,real(ux)-u1,'k:','DisplayName','error');
            %plot(x,imag(ux),'--','DisplayName','im(u\_ASM)');
            xlabel('x');ylabel('du/dx'); grid on; grid minor;
            set(gca,'FontSize',16);
            set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
            % legend('exact','nm\_spec','nm\_error'); 
            axis tight; %legend('Location','Best');
%             saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/ASM_mode1.67.jpg');
%             saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/ASM_multimode.jpg');
%             saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/ASM_5modes.jpg');
%             saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/ASM_multimodeFullRank.jpg');
        end
    end

    function uj = compute_derivative_jumps(u,njumps)
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % fdi -> i-th derivative
%         % forward difference stencils with second order accuracy
%         fd1 = [-3/2    2   -1/2 ]/dx;
%         fd2 = [  2    -5     4  -1 ]/dx^2;
%         fd3 = [-5/2    9   -12   7   -3/2 ]/dx^3;
%         fd4 = [  3   -14    26 -24     11  -2 ]/dx^4;
% 
%         % forward difference stencils with fourth order accuracy
%         fd1 = [-25/12   4     -3       4/3    -1/4]/dx;
%         fd2 = [15/4    -77/6   107/6  -13      61/12  -5/6]/dx^2;
%         fd3 = [-49/8    29    -461/8   62     -307/8   13      -15/8]/dx^3;
%         fd4 = [28/3    -111/2  142    -1219/6  176    -185/2    82/3    -7/2]/dx^4;
%         
%         % forward difference stencils with sixth order accuracy
%         fd1 = [-49/20   6     -15/2     20/3    -15/4   6/5 -1/6    ]/dx;
%         fd2 = [469/90  -223/10 879/20  -949/18     41  -201/10  1019/180    -7/10   ]/dx^2;
%         fd3 = [-801/80  349/6   -18353/120  2391/10 -1457/6     4891/30     -561/8  527/30  -469/240 ]/dx^3;
%         % fd4 = []/dx^4;
%         
%         % backward difference stencils
%         bd1 = -flip(fd1);
%         bd2 =  flip(fd2);
%         bd3 = -flip(fd3);
%         % bd4 =  flip(fd4);
% 
%         stw=numel(fd1);
%         uj(1) = u(end)-u(1);
%         uj(2) = bd1*u((end-(stw-1)):end) - fd1*u(1:(stw+0));
%         uj(3) = bd2*u((end-(stw+0)):end) - fd2*u(1:(stw+1));
%         uj(4) = bd3*u((end-(stw+1)):end) - fd3*u(1:(stw+2));
% %         uj(5) = bd4*u((end-(stw+2)):end) - fd4*u(1:(stw+3));
% %         uj(6) = bd5*u((end-(stw+3)):end) - fd5*u(1:(stw+4));
% %         uj(7) = bd6*u((end-(stw+4)):end) - fd6*u(1:(stw+5));
% %         uj(8) = bd7*u((end-(stw+5)):end) - fd7*u(1:(stw+6));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute jumps for the function & derivatives
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        uj = zeros(1,njumps);
        switch(jump_method)
          case(1) % analytical expressions for jumps      
          switch(IC_profile)
              case(5)
                m = m_monomials;
                % exact jump conditions for the monomial u=x^m
%                 uj(1) = ( (pi)^(m  ) - (-pi)^(m  ) );
%                 uj(2) = ( (pi)^(m-1) - (-pi)^(m-1) )*m;
%                 uj(3) = ( (pi)^(m-2) - (-pi)^(m-2) )*m*(m-1);
%                 uj(4) = ( (pi)^(m-3) - (-pi)^(m-3) )*m*(m-1)*(m-2);
                mcoeff = 1;
                for i=1:njumps
                    uj(i) = ( (pi)^(m-i+1) - (-pi)^(m-i+1) )*mcoeff;
                    mcoeff= mcoeff*(m-i+1);
                end

%                 %% alternate approach
%                 uj(1) = pi^m - (-pi)^m;
%                 jumps = polynomial_derivatives(m,njumps-1);
%                 for i = 2:njumps
%                     uj(i) = jumps(i-1);
%                 end
                
              case(8)
                for i=1:njumps
                    uj(i) = exp(pi)-exp(-pi);
                end
              case(9)
%                 th = 0; a = -2; b= 3/4; % * 2*pi;
                pmat = pascal(njumps);
                for i=0:njumps-1
                    tmp = 0;
                    for j=0:i
%                         pmat = pascal(i);
%                         minordiag = sort(find(fliplr(eye(size(pmat)))), 'descend');
%                         iCj = pmat(minordiag);  
                        iCj = pmat(i-j+1,j+1);
                        fb  = a^(i-j)*exp(a*( pi+th));
                        fa  = a^(i-j)*exp(a*(-pi+th));
                        if(mod(j,2)==0) % j even
                            gb  = (-1)^(j/2)*b^j*sin(b*( pi+th));
                            ga  = (-1)^(j/2)*b^j*sin(b*(-pi+th));
                        else            % j odd
                            gb  = (-1)^((j-1)/2)*b^j*cos(b*( pi+th));
                            ga  = (-1)^((j-1)/2)*b^j*cos(b*(-pi+th));
                        end
                        tmp = tmp + iCj*(fb*gb-fa*ga);
                    end
                    uj(i+1) = tmp;
                end
              case(10)
                uj(1) = log(pi+th)-log(-pi+th);
                uj(2) = (1./(pi+th) - 1./(-pi+th));
                fact  = 1;
                for i = 3:njumps
                    fact  = fact * -(i-2);
                    uj(i) = (1./(pi+th).^(i-1) - 1./(-pi+th).^(i-1)) * fact;
                end
  
              case(11)
%                 u0 = exp(-((x-x0)/hw).^2);      % initial condition
%                 u1 = exp(-((x-x0)/hw).^2).*(-2*(x-x0)/hw^2);    % exact first derivative of u0
                uj(1) = exp(-((pi-x0)/hw).^2)-exp(-((-pi-x0)/hw).^2);
                jumps = gaussian_derivatives(x0, hw, njumps-1);
                for i = 2:njumps
                    uj(i) = jumps(i-1);
                end
  
              case(2)
                for i=1:njumps
                    if(mod(i,4)>=2)
                        sn = -1;
                    else
                        sn = 1;
                    end
                    uj(i) = 2*sum( sn* k.^(i-1).*ampS.*sin(k.*pi),2);
                    if(IC_profile==2 && mod(i,2)==0) % pure sin waves
                        uj(i) = 0; % odd jumps are zero
                    end
                    if(IC_profile==3 && mod(i,2)==1) % pure cos waves
                        uj(i) = 0; % even jumps are zero
                    end
                end
                Jm = uj.';
                max_Jm = max(abs(Jm))
                  
              otherwise
%                 uj(1) = 2*sum(       sin(k.*pi),2);
%                 uj(2) = 2*sum(-k   .*sin(k.*pi),2);
%                 uj(3) = 2*sum(-k.^2.*sin(k.*pi),2);
%                 uj(4) = 2*sum( k.^3.*sin(k.*pi),2);
%                 uj(5) = 2*sum( k.^4.*sin(k.*pi),2);
%                 uj(6) = 2*sum(-k.^5.*sin(k.*pi),2);
%                 uj(7) = 2*sum(-k.^6.*sin(k.*pi),2);
%                 uj(8) = 2*sum( k.^7.*sin(k.*pi),2);
                for i=1:njumps
                    if(mod(i,4)>=2)
                        sn = -1;
                    else
                        sn = 1;
                    end
                    if(mod(i,2)==1) 
                        amp = ampS;
                    else
                        amp = ampC;
                    end
                    uj(i) = 2*sum( sn* k.^(i-1).*amp.*sin(k.*pi),2);
                    if(IC_profile==2 && mod(i,2)==0) % pure sin waves
                        uj(i) = 0; % odd jumps are zero
                    end
                    if(IC_profile==3 && mod(i,2)==1) % pure cos waves
                        uj(i) = 0; % even jumps are zero
                    end
                end
                Jm = uj.';
                max_Jm = max(abs(Jm))
          end
                
                
          case(2) % finite difference approximation for jumps
            % r = 10;  % accuracy of each derivative
            stw_min = 1 + r;
            stw_max = njumps-1 + r;
            fprintf('highest derivative is order %d\n',njumps-1);
            fprintf('min stencil width = %d\nmaximum stencil width = %d\n',stw_min,stw_max);
            for i=1:njumps
                q = i-1;    % 'q'th derivative
                stw = r+q;
                % stw = stw_max;
                
                % weights for q_th derivative with O(dx^r) accuracy
                fdq = classic_fd_stencil(0,stw-1,q).'; % 1:stw
                bdq = classic_fd_stencil(stw-1,0,q).'; % end-stw+1:end
%                 uj(i) = (bdq*u((end-(stw-1)):end) - fdq*u(1:stw))/dx^q;
                if(mod(i,2)==0)
                    uj(i) = bdq * ( u((end-(stw-1)):end) + flip(u(1:stw)) )/dx^q;
                else
                    uj(i) = bdq * ( u((end-(stw-1)):end) - flip(u(1:stw)) )/dx^q;
                end
            end
            uj(1) = u(end)-u(1);     
            
          case(3) % finite difference approximation for jumps with fixed stencil width
            % r = 10;  % accuracy of each derivative
            % u = c_0+c_1*x+c_2*x^2+...+c^r*x^r
            % u(0) = c_0 
            % u(1) = c_0+c_1*(  dx)+c_2*(  dx)^2+...+c^r*(  dx)^r
            % u(2) = c_0+c_1*(2*dx)+c_2*(2*dx)^2+...+c^r*(2*dx)^r
            % ...
            % u(w) = c_0+c_1*(w*dx)+c_2*(w*dx)^2+...+c^r*(w*dx)^r
            % w=stw
            stw_min = 1 + r;
            stw_max = njumps-1 + r;
            fprintf('highest derivative is order %d\n',4*n_apmodes-1);
            fprintf('min stencil width = %d\nmaximum stencil width = %d\n',stw_min,stw_max);
            fprintf('maximum stencil width %d is chosen\n',stw_max);
            for i=1:njumps
                q = i-1;    % 'q'th derivative
                stw = r+1; %r+q; %q; %2;
%                 stw = min(stw,stw_max);
%                 stw = max(stw,stw_min);
                stw = stw_max;
%                 stw = stw_min;

                % fprintf('%d th derivative with stencil width %d \n',q,stw)
                % weights for q_th derivative with O(dx^r) accuracy
                fdq = classic_fd_stencil(0,stw-1,q).'; % 1:stw
                bdq = classic_fd_stencil(stw-1,0,q).'; % end-stw+1:end
%                 fwd = fdq.';
%                 bwd = bdq.';
%                 fdw = [fwd,flip(bwd)]
% 
%                 fdq2 = fd_weights(stw, q, 'forward');
%                 bdq2 = fd_weights(stw, q, 'backward');
%                 fdw2 = [fdq2,flip(bdq2)]
                
%                 uj(i) = (bdq*u((end-(stw-1)):end) - fdq*u(1:stw))/dx^q;                
                if(mod(i,2)==0)
                    uj(i) = bdq * ( u((end-(stw-1)):end) + flip(u(1:stw)) )/dx^q;
                else
                    uj(i) = bdq * ( u((end-(stw-1)):end) - flip(u(1:stw)) )/dx^q;
                end
            end
            uj(1) = u(end)-u(1);  
            %uj(10:end) = 0;

          case(4) % spline/polynomial approximation
            stw_min = 1 + r;
            stw_max = njumps-1 + r;
            stw = stw_max;
%             fprintf('highest derivative is order %d\n',njumps-1);
%             fprintf('min stencil width = %d\nmaximum stencil width = %d\n',stw_min,stw_max);
%             fprintf('maximum stencil width %d is chosen\n',stw_max);
            jumps = compute_jumps(x,u,njumps-1,stw);
            uj(1) = u(end)-u(1);
            for i = 2:njumps
                uj(i) = jumps(i-1);
            end

          case(33333)  % comment this out to compute error in the derivatives
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % compute the error in all jumps analytically
            uj_ex = zeros(1,njumps);        
            for i=1:njumps
                if(mod(i,4)>=2)
                    sn = -1;
                else
                    sn = 1;
                end
                uj_ex(i) = 2*sum( sn* k.^(i-1).*sin(k.*pi),2);
                if(IC_profile==2 && mod(i,2)==0) % pure sin waves
                    uj_ex(i) = 0; % odd jumps are zero
                end
                if(IC_profile==3 && mod(i,2)==1) % pure cos waves
                    uj_ex(i) = 0; % even jumps are zero
                end
            end

                % exact jump condition for IC_Profile=10 (log(x))
                uj_ex(1) = log(pi+th)-log(-pi+th);
                uj_ex(2) = (1./(pi+th) - 1./(-pi+th));
                fact  = 1;
                for i = 3:njumps
                    fact  = fact * -(i-2);
                    uj_ex(i) = (1./(pi+th).^(i-1) - 1./(-pi+th).^(i-1)) * fact;
                end
% %             uj = uj_ex;
            Jm = uj_ex.';
            max_Jm = max(abs(Jm))
            (uj - uj_ex).'
            for i = 1:length(uj_ex)
                fprintf('i = %d  uj = %e\tuj_ex = %e\tej = %e\n',i,uj(i),uj_ex(i),uj(i)-uj_ex(i));
            end
%             compare = [uj; uj_ex; uj-uj_ex].'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
          case(4) % Hierarchical finite difference approximation for jumps
            % r = 10;  % accuracy of each derivative
            stw_min = 1 + r;
            stw_max = njumps-1 + r;
            fprintf('highest derivative is order %d\n',4*n_apmodes-1);
            fprintf('min stencil width = %d\nmaximum stencil width = %d\n',stw_min,stw_max);
            for i=1:njumps
                q = i-1;
                stw = r+q; %q; %2;
                stw = stw_max;

                % fprintf('%d th derivative with stencil width %d \n',q,stw)
                % weights for q_th derivative with O(dx^r) accuracy
                fdq = classic_fd_stencil(0,stw-1,q).'; % 1:stw
                bdq = classic_fd_stencil(stw-1,0,q).'; % end-stw+1:end
%                 fwd = fdq.'
%                 bwd = bdq.'
%                 uj(i) = (bdq*u((end-(stw-1)):end) - fdq*u(1:stw))/dx^q;                
                if(mod(i,2)==0)
                    uj(i) = bdq * ( u((end-(stw-1)):end) + flip(u(1:stw)) )/dx^q;
                else
                    uj(i) = bdq * ( u((end-(stw-1)):end) - flip(u(1:stw)) )/dx^q;
                end
            end
            uj(1) = u(end)-u(1);     
            
          otherwise % high resolution FD approximation for jumps
            stw = length(u);
            for i=1:njumps
                q = i-1;
                % weights for q_th derivative with O(dx^r) accuracy
                fdq = high_resolution_stencil(0,stw-1,q,r).'; % 1:stw
                bdq = high_resolution_stencil(stw-1,0,q,r).'; % end-stw+1:end
                uj(i) = (bdq*u((end-(stw-1)):end) - fdq*u(1:stw))/dx^q;
            end
            uj(1) = u(end)-u(1);        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if(regularize_Jmat)
            uj = uj + eps_Jmat;
        end        
    end

    function ux = dudx_Roache(x,u,order_disc,plotTag,ux_FFT)
        njumps = order_disc + 1;
        uj_R = compute_derivative_jumps(u,njumps);
        x_a = x(1);
        x_b = x(end);
        Lx = x(end)-x(1);
        delx = x(2)-x(1);
%         x_scaled = (x-x(1))/Lx; % domain scaled to [0,1]
        %%% Roache's method to find non-periodic part using Jumps %%%
%         % Compute boundary derivatives analytically (up to N)
%         uj = zeros(1,njumps); % uj0 to ujN
%         syms xs
%         f_sym = exp(-2*xs) * sin(2*pi*xs);
%         for n = 0:njumps-1
%             uj(n+1) = double(subs(diff(f_sym, xs, n), 1)) - double(subs(diff(f_sym, xs, n), 0));
%         end
        
        % Compute polynomial coefficients for g(x)
        a_R = zeros(1,njumps);
        a_R(njumps) = 1/(Lx)*uj_R(njumps) / factorial(njumps); % Highest order coefficient
        for n = njumps-1:-1:1
            coeff_sum = 0;
            for j = n+1:njumps%-1
                coeff_sum = coeff_sum + a_R(j) * nchoosek(j,n-1) * (x_b^(j-n+1)-x_a^(j-n+1));
            end
            a_R(n) = (1/Lx) * (uj_R(n)/factorial(n) - coeff_sum/n);
        end
        % a_R(1) = 0; % double(subs(f_sym, 0)); % Constant term
        
        % Evaluate ua(x) and its derivative analytically
        ua = zeros(size(x));
        uax = zeros(size(x));
        for i = 1:njumps
            ua  = ua  + a_R(i)      * x.^(i);
            uax = uax + a_R(i) *(i) * x.^(i-1);
        end
        
        % Residual function (smooth, periodic)
        up = u-ua;
        [upx,upF] = dudx_FFT(up,L,0); % standard MATLAB FFT
        eF = up - upF;
%         % size_up = size(up)
%         % M = floor(size(up)/2)*2;
%         M = floor(N/2)*2;
%         kvec = 2*pi/L * ifftshift(-M/2:M/2-1).';
%         % FFT-based derivative of the periodic part
%         upk = fft(ifftshift(up(1:end-1)));
%         upx = real(fftshift(ifft(1i * kvec .* upk)));
        
        % Total derivative
        ux = upx + uax;

        figure(224); clf;
        lw = 2;
%         plot(xfine,ufine,'k-','DisplayName','$u(x)$, original'); hold on;
        plot(x,u,'k*-','DisplayName','$u(x)$, original'); hold on;
        plot(x,ua,'r^-','DisplayName','$u_a(x)$, aperiodic');
        plot(x,up,'bo-','DisplayName','$u_p(x)$, periodic');
        %plot(x,imag(uF),'--','DisplayName','error_im');
        xlabel('$x$','Interpreter','latex');ylabel('$u(x)$','Interpreter','latex'); grid on; grid minor;
        set(gca,'FontSize',16); legend('Location','Best');
        title('Decomposition with Roache method','FontSize',12);
        set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
        % legend('exact','nm\_spec','nm\_error');  
        axis tight; 
        hl = legend('Location','NorthWest'); set(hl,'Interpreter','latex');
        if(fig_save)
            saveas(gcf,fig_u);
        end
%         saveas(gcf,'~/Dropbox/RapakaS_ASM/Figures/u_eax_2mode.jpg');
%         saveas(gcf,'~/Dropbox/RapakaS_ASM/Figures/u_Nk32_n3.jpg');
                

        fs = 15;
        lw = 2;

        if(plotTag)            
            lw=2;
            e0_u_N(nn,1) = max(abs(eF));
            fprintf('||u-u_ASM||_2  : %.2e\n',norm(abs(eF),2)*dx^(1/2));
            fprintf('||u-u_ASM||_Inf: %.2e\n\n',max(abs(eF)));
%             fprintf("Nyq_comp: \nu : |  uk| = %10.2e \nu': |k*uk| = %10.2e\n\n",abs(uk(1)),abs(kn(1)*uk(1)));
            
            figure(130); clf;    
%             subplot(3,1,1);
%             plot((kn),real((uk)),'.-','DisplayName','real'); hold on;
%             plot((kn),imag((uk)),'-','DisplayName','imag');
% %             max_uk_im = max(abs(imag(uk)))
%             xlabel('$k_n$','Interpreter','Latex'); ylabel('$\hat{u}_{p,n}$','Interpreter','Latex');
%             legend('Location','Best'); grid on; grid minor;
%             set(gca,'FontSize',16); axis tight;
%             xlim([-N/2,N/2]);
%             set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
%             % set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 1);

%             subplot(2,1,1);
%             yyaxis left;
%             plot(x,u,'r-','DisplayName','exact'); hold on;
%             plot(x,real(uF),'bs','DisplayName','GFS');
%             ylabel('$u(x)$','Interpreter','latex');             
%             yyaxis right;
%             set(gca, 'YScale', 'log');
%             semilogy(x,eF,'k:','DisplayName','error');
%             %plot(x,imag(uF),'--','DisplayName','error_im');
%             ylim([10^-16,10^-14]);
%             xlabel('$x$','Interpreter','latex');
%             ylabel('$|error(x)|$','Interpreter','latex'); 
%             grid on; grid minor;
%             set(gca,'FontSize',16);
%             set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
%             % legend('exact','nm\_spec','nm\_error');  
%             axis tight; legend('Location','Best');
%             
%             subplot(2,1,2);
            yyaxis left;
%             plot(xfine,u1fine,'m','DisplayName','analytical'); hold on;
            plot(x,u1,'m','DisplayName','analytical'); hold on;
            plot(x,ux_FFT,'bs','DisplayName','FFT'); hold on;
            plot(x,real(ux),'-o','DisplayName','Roache');
            %plot(x,imag(ux),'--','DisplayName','im(u\_ASM)');
            ylabel("$u^{'}(x)$",'Interpreter','latex'); 
            yyaxis right;
            set(gca, 'YScale', 'log');
            semilogy(x,abs(ux_FFT-u1),':','DisplayName','error (FFT)');
            semilogy(x,abs(real(ux)-u1),'--','DisplayName','error (Roache)');
            ylim([10^-10,10^2]);
            xlabel('$x$','Interpreter','latex');
            ylabel("$|error(x)|$",'Interpreter','latex'); 
            grid on; grid minor;
            set(gca,'FontSize',16);
            set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
            % legend('exact','nm\_spec','nm\_error'); 
            axis tight; legend('Location','Best');
            if(fig_save)
                saveas(gcf,fig_up);
            end
%             saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/ASM_mode1.67.jpg');
%             saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/ASM_multimode.jpg');
%             saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/ASM_5modes.jpg');
%             saveas(gcf,'~/Dropbox/RapakaS_ASM/Figures/up_ramp_4.jpg');
%             saveas(gcf,'~/Dropbox/RapakaS_ASM/Figures/up_eax_2mode.jpg');
%             saveas(gcf,'~/Dropbox/RapakaS_ASM/Figures/up_Nk32_n3.jpg');
        end
    end

    function [ux,kn] = dudx_GFS(u,plotTag,ux_compare,method)
        uj = compute_derivative_jumps(u,4*nmax);
        switch(n_apmodes)
        case(-1)
            fprintf('[u ] = %.2e\n[u1] = %.2e\n[u2] = %.2e\n[u3] = %.2e\n',uj(1),uj(2),uj(3),uj(4));
            fprintf('[u2]/[u ] = %.2e\n[u3]/[u1] = %.2e\n',-uj(3)/uj(1),-uj(4)/uj(2));

            if( abs(uj(1))>=10^-15 )  
                k_i = sqrt(-uj(3)/uj(1)); 
                if( k_i~=0 )
                    u_i =  uj(1)/(2*sin(k_i*pi));
                else
                    u_i = 0;
                end
            else
                k_i = 0; 
                u_i = 0;
            end
            if( abs(uj(2))>=10^-15 ) 
                k_r = sqrt(-uj(4)/uj(2));
                if(k_r~=0)
                    u_r = -uj(2)/(2*k_r*sin(k_r*pi));
                else
                    u_r = 0;
                end
            else
                k_r = 0; 
                u_r = 0;
            end
            fprintf('\nk_r = %.4e\nk_i = %.4e\n',k_r,k_i);
            fprintf('u_r = %.4e\nu_i = %.4e\n\n',u_r,u_i);

            ua   =        u_r*cos(k_r*x) +       u_i*sin(k_i*x);
            uax  = -k_r  *u_r*sin(k_r*x) + k_i  *u_i*cos(k_i*x);
            uaxx = -k_r^2*u_r*cos(k_r*x) - k_i^2*u_i*sin(k_i*x);
            up   = u-ua;
        
%             % verify jumps
%             uj(1) = ua(end)-ua(1);
%             uj(2) = -2*k_r  *u_r*sin(k_r*pi);
%             uj(3) = -2*k_i^2*u_i*sin(k_i*pi);
%             uj(4) =  2*k_r^3*u_r*sin(k_r*pi);
%             fprintf('[u ] = %.2e\n[u1] = %.2e\n[u2] = %.2e\n[u3] = %.2e\n',uj(1),uj(2),uj(3),uj(4));
%             fprintf('[u2]/[u ] = %.2e\n[u3]/[u1] = %.2e\n',-uj(3)/uj(1),-uj(4)/uj(2));
        
        case(-2)
            uj.'
            A = (uj(1)*uj(7)-uj(3)*uj(5))/(uj(3)^2-uj(1)*uj(5));
            B = (uj(5)^2    -uj(3)*uj(7))/(uj(3)^2-uj(1)*uj(5));
            D = sqrt(A^2-4*B);
            ki(1) = sqrt((A+D)/2);
            ki(2) = sqrt((A-D)/2);

            A
            B
            D

            A = (uj(2)*uj(8)-uj(4)*uj(6))/(uj(4)^2-uj(2)*uj(6));
            B = (uj(6)^2    -uj(4)*uj(8))/(uj(4)^2-uj(2)*uj(6));
            D = sqrt(A^2-4*B);
            kr(1) = sqrt((A+D)/2);
            kr(2) = sqrt((A-D)/2);

            A
            B
            D
            
            kr 
            ki 

            ui(1) = -( uj(3) + ki(2)^2*uj(1) )/( 2*(ki(1)^2-ki(2)^2)*sin(ki(1)*pi) );
            ui(2) =  ( uj(3) + ki(1)^2*uj(1) )/( 2*(ki(1)^2-ki(2)^2)*sin(ki(2)*pi) );
            ur(1) =  ( uj(4) + kr(2)^2*uj(2) )/( 2*(kr(1)^2-kr(2)^2)*sin(kr(1)*pi)*kr(1) );
            ur(2) = -( uj(4) + kr(1)^2*uj(2) )/( 2*(kr(1)^2-kr(2)^2)*sin(kr(2)*pi)*kr(2) );

            ur
            ui
            ua   = sum(        ur.*cos(x*kr),2) + sum(        ui.*sin(x*ki),2);
            uax  = sum(-kr   .*ur.*sin(x*kr),2) + sum( ki   .*ui.*cos(x*ki),2);
            uaxx = sum(-kr.^2.*ur.*cos(x*kr),2) + sum(-ki.^2.*ui.*sin(x*ki),2);
            up   = u-ua;

        otherwise
%             kmat = ones(n_apmodes);
%             Jmat = zeros(n_apmodes);
%             Jvec = zeros(n_apmodes,1);
%             Jvec1= zeros(n_apmodes,1);
%             eind = 2*(n_apmodes-1);
            
            %============ sin modes ==================
            cnt = 1;
            while(1)
                kmat = ones(ns);
                Jmat = zeros(ns);
                Jvec = zeros(ns,1);
                Jvec1= zeros(ns,1);
                eind = 2*(ns-1);
                for i=1:ns
                    sind = 2*i-1;
                    Jmat(i,:) = uj(sind+(0:2:eind));
                    Jvec(i)   = uj(sind+eind+2);
                    Jvec1(i)  = uj(sind);
                end
%                 if(regularize_Jmat)
%                     Jmat = Jmat + eps_Jmat*eye(size(Jmat));
%                 end
%                 Jmat_sin = Jmat
                rankJ = rank(Jmat);
                detJ = det(Jmat);
                fprintf('|J|_sin = %e\n',detJ);
                if( (abs(detJ)<TOL_DetJ || rankJ<ns) && check_rank==1 )
                    fprintf('rank(J) = %d<%d, ', rankJ,ns);                    
                    if(cnt>1)
                        ns = ns-1;
                    else
                        ns = rankJ; %ns-1;
                    end
                    fprintf('ns=%d (#sin modes reduced)\n', ns);
                    cnt = cnt + 1;
                else
                    break;
                end
            end
        if(ns==0 || all(Jvec==0))
            ki = 0;
            ui = 0;
        else
            evec = - (Jmat\Jvec);
%             evec = - (inv(Jmat) * Jvec);
            if(pseudoInv == 1)
                disp('evec (sin) with J^-1')
                evec
%                 evec = - (pinv(Jmat,TOL_PINV) * Jvec);
                evec = - (my_pinv(Jmat,TOL_PINV) * Jvec);
                disp('evec (sin) with psuedoinverse of J')
                evec
            end
            evec
            Jmat
            Jvec
            evec = [1; flip(evec)] .* ((-1).^(0:ns)).';
            ki = sqrt(roots(evec)).';   %row vector
            for i=1:ns
                kmat(i,:) = (-1)^(i+1) * ki.^(2*(i-1));
            end
            ui = 0.5* (kmat\Jvec1).' ./ sin(ki*pi);   %row vector 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % regularize '0' wave number with eps_sinc
            num_zeros = sum(ki == 0);
            if(all(ki==0))
                disp('All elements in k_sin are zero.');
                ki = 0;
                ui = 0;
            end
%             if(num_zeros==1)
%                 ki(ki == 0) = eps_sinc; 
%             end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        end
            fprintf('\n');
%             ki
%             ki_simple = sqrt((Jvec./Jvec1).^(1/ns))
%             Jvec1
%             Jvec2 = Jvec./((-ki.^2).^ns).'
            %============ cos modes ==================
            cnt = 1;
            while(1)
                kmat = ones(nc);
                Jmat = zeros(nc);
                Jvec = zeros(nc,1);
                Jvec1= zeros(nc,1);
                eind = 2*(nc-1);
                for i=1:nc
                    sind = 2*i;
                    Jmat(i,:) = uj(sind+(0:2:eind));
                    Jvec(i)   = uj(sind+eind+2);
                    Jvec1(i)  = uj(sind);
                end
%                 if(regularize_Jmat)
%                     Jmat = Jmat + eps_Jmat*eye(size(Jmat));
%                 end
%                 Jmat_cos = Jmat
                rankJ = rank(Jmat);
                detJ = det(Jmat);
                fprintf('|J|_cos = %e\n',detJ);
                if( (abs(detJ)<TOL_DetJ || rankJ<nc) && check_rank==1 )
                    fprintf('rank(J) = %d<%d, ', rankJ,nc);
                    if(cnt>1)
                        nc = nc-1;
                    else
                        nc = rankJ; %nc-1;
                    end
                    fprintf('nc=%d (#cos modes reduced)\n', nc);
                    cnt = cnt + 1;
                else
                    break;
                end
            end
            
            fprintf('Jump vector (J_m*dx^m): \n');
            dxnvec = ones(length(uj.'),1);
            for i=1:length(uj.')
                dxnvec(i) = dx^(i-1);
            end
%             dxnvec
            uj.'
%             uj.' .* dxnvec

        if(nc==0 || all(Jvec==0))
            kr = 0;
            ur = 0;
        else
            evec = - (Jmat\Jvec);
%             evec = - (inv(Jmat) * Jvec);
            if(pseudoInv == 1)
                disp('evec (cos) with J^-1')
                evec
%                 evec = - (pinv(Jmat,TOL_PINV) * Jvec);
                evec = - (my_pinv(Jmat,TOL_PINV) * Jvec);
                disp('evec (cos) with psuedoinverse of J')
                evec
            end
            evec
            Jmat
            Jvec
            evec = [1; flip(evec)] .* ((-1).^(0:nc)).';
            kr = sqrt(roots(evec)).';   %row vector
            for i=1:nc
                kmat(i,:) = (-1)^(i+1) * kr.^(2*(i-1));
            end
            ur = -0.5* (kmat\Jvec1).' ./ (kr.*sin(kr*pi)); %row vector
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % regularize '0' wave number with eps_sinc
            num_zeros = sum(kr == 0);
            if(all(kr==0))
                disp('All elements in k_cos are zero.');
                kr = 0;
                ur = 0;
            end
%             if(num_zeros==1)
%                 kr(kr == 0) = eps_sinc; 
%             end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
            %=========================================

            figure(33)
            scatter(real(kr),imag(kr),'rs');
            xlabel('$Re(k_c)$','Interpreter','latex');
            ylabel('$Im(k_c)$','Interpreter','latex'); grid on; grid minor;
            set(gca,'FontSize',16); %legend;
            figure(34)
            scatter(real(ki),imag(ki),'ro');
            xlabel('$Re(k_s)$','Interpreter','latex');
            ylabel('$Im(k_s)$','Interpreter','latex'); grid on; grid minor;
            set(gca,'FontSize',16); %legend;
            
            if(~isempty(ki))
%                 k_sin = ki.' 
%                 uhat_sin = ui.'
                fprintf('\nk_sin = ');
                fprintf('\n  %4.2e + %4.2ei',[real(ki);imag(ki)]);
                fprintf('\n\nuhat_sin = ');
                fprintf('\n  %4.2e + %4.2ei',[real(ui);imag(ui)]);
                fprintf('\n\n');
            end
            if(~isempty(kr)) 
%                 k_cos = kr.' 
%                 uhat_cos = ur.'
                fprintf('\nk_cos = ');
                fprintf('\n  %4.2e + %4.2ei',[real(kr);imag(kr)]); % fprintf goes down columns for printing next number
                fprintf('\n\nuhat_cos = ');
                fprintf('\n  %4.2e + %4.2ei',[real(ur);imag(ur)]);
                fprintf('\n\n');
            end
%             sxk = sin(x*ki);
%             ind = isinf(imag(sxk));
%             sxk(ind) = real(sxk(ind))+1i*imag(sxk(ind));
%             sxk
%             test   = sum(  real(ui).*real(sin(x*ki))-imag(ui).*imag(sin(x*ki))  ,2)
            if(mask_zero_amp==1)
                mask_cos = abs(ur)<TOL_MASK
                ur(mask_cos) = 0.0;
                kr(mask_cos) = 0.0;
    
                mask_sin = abs(ui)<TOL_MASK 
                ui(mask_sin) = 0.0;
                ki(mask_sin) = 0.0;
            end

            ua   = sum(        ur.*cos(x*kr),2) + sum(        ui.*sin(x*ki),2);
            uax  = sum(-kr   .*ur.*sin(x*kr),2) + sum( ki   .*ui.*cos(x*ki),2);
            uaxx = sum(-kr.^2.*ur.*cos(x*kr),2) + sum(-ki.^2.*ui.*sin(x*ki),2);
            ua   = real(ua);
            uax  = real(uax);
            uaxx = real(uaxx);
            up   = u-ua;

%             ua(1:10)
%             up(1:10)

%             ui
%             ki
% %             sum(        ur.*cos(x*kr),2)
%             sum(        sin(x*ki),2)

%             size(ur.*cos(x*kr))
%             ur.*cos(x*kr)
%             ui.*sin(x*ki)

%             ur
%             kr
%             ui
%             ki
%             ua
%             ureal = real(ua)
%             uimag = imag(ua)
        end

        figure(2); clf;
        lw = 2;
%         plot(xfine,ufine,'k-','DisplayName','$u(x)$, original'); hold on;
        plot(x,u,'k*-','DisplayName','$u(x)$, original'); hold on;
        plot(x,ua,'r^-','DisplayName','$u_a(x)$, aperiodic');
        plot(x,up,'bo-','DisplayName','$u_p(x)$, periodic');
        %plot(x,imag(uF),'--','DisplayName','error_im');
        xlabel('$x$','Interpreter','latex');ylabel('$u(x)$','Interpreter','latex'); grid on; grid minor;
        set(gca,'FontSize',16); legend;
        set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
        % legend('exact','nm\_spec','nm\_error');  
        axis tight; 
        hl = legend('Location','NorthWest'); set(hl,'Interpreter','latex');
        if(fig_save)
            saveas(gcf,fig_u);
        end
%         saveas(gcf,'~/Dropbox/RapakaS_ASM/Figures/u_eax_2mode.jpg');
%         saveas(gcf,'~/Dropbox/RapakaS_ASM/Figures/u_Nk32_n3.jpg');
        
        M = floor(N/2)*2;
        kn = 2*pi/L * (-M/2:M/2-1).'; % modes        
%         [upx,uF] = dudx_FFT(up,L,0);
        
        isperiodic = true;
%         isperiodic = false;
%         fprintf('isperiodic = %d\n', isperiodic);
        uk = nmFFT(up,kn,isperiodic);   % FFT through non-integer basis
        % uk(M/2+1)=0;    % nullify odd-ball wave number    
        upx = nmiFFT(1i*kn.*uk,kn,isperiodic);  % iFFT
        %%%%%%%%%%%%%%%%%%%%%%%%%x
        upx = dudx_FFT(up,L,0); % standard MATLAB FFT
        %%%%%%%%%%%%%%%%%%%%%%%%%x
        upx = real(upx);
        % copy of the original function via iFFT
        uF = nmiFFT(uk,kn,isperiodic);  % iFFT
        uF = real(uF)+ua;
        eF = uF-u;
        
        ux = uax + upx;

        fs = 15;
        lw = 2;
%         % Plot the original signal
%         figure(44);clf;
%         subplot(2, 1, 1);
%         h = plot(x, u,'o-','LineWidth',lw);
% %         title('$f(t)=0.7\sin(5.3 t)+1.0\sin(12.4 t)$','Interpreter','latex');
%         % title('Time-domain signal','Interpreter','latex');
%         xlabel('time, t (s)','Interpreter','latex');
%         ylabel('$f(t)$','Interpreter','latex');
%         set(gca,'FontSize',fs);
%         axis tight;
%         
%         grid on;
%         
%         % Plot the Fourier spectrum
%         subplot(2, 1, 2);
%         h = plot(abs(ki), abs(ui),'o','LineWidth',lw, 'DisplayName', '$\tilde{u}$'); hold on;
%         h = plot(abs(kr), abs(ur),'^-','LineWidth',lw,'DisplayName', '$\hat{u}$');
%         for i = 1:length(ki)
%             line([abs(ki(i)), abs(ki(i))], [0, abs(ui(i))],'LineWidth',lw, 'HandleVisibility', 'off'); % Vertical line
%             line([abs(kr(i)), abs(kr(i))], [0, abs(ur(i))],'LineWidth',lw, 'HandleVisibility', 'off'); % Vertical line
%         end
% %         title('Single-sided amplitude spectrum of f(t)','Interpreter','latex');
%         xlabel('frequency, $\omega (Hz)$','Interpreter','latex');
%         ylabel('$|\hat{f}(\omega)|$','Interpreter','latex');
%         set(gca,'FontSize',fs);
%         legend('Interpreter', 'latex','Location','Best');
%         grid on;
% %         axis tight;
% %         saveas(gcf,'~/Dropbox/RapakaS_ASM/Figures/f_asm_integer.jpg');
        

        if(plotTag)            
            lw=2;
            e0_u_N(nn,1) = max(abs(eF));
            fprintf('||u-u_GFS||_2  : %.2e\n',norm(abs(eF),2)*dx^(1/2));
            fprintf('||u-u_GFS||_Inf: %.2e\n\n',max(abs(eF)));
%             fprintf("Nyq_comp: \nu : |  uk| = %10.2e \nu': |k*uk| = %10.2e\n\n",abs(uk(1)),abs(kn(1)*uk(1)));
            
            figure(13); clf;    
%             subplot(3,1,1);
%             plot((kn),real((uk)),'.-','DisplayName','real'); hold on;
%             plot((kn),imag((uk)),'-','DisplayName','imag');
% %             max_uk_im = max(abs(imag(uk)))
%             xlabel('$k_n$','Interpreter','Latex'); ylabel('$\hat{u}_{p,n}$','Interpreter','Latex');
%             legend('Location','Best'); grid on; grid minor;
%             set(gca,'FontSize',16); axis tight;
%             xlim([-N/2,N/2]);
%             set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
%             % set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 1);

%             subplot(2,1,1);
%             yyaxis left;
%             plot(x,u,'r-','DisplayName','exact'); hold on;
%             plot(x,real(uF),'bs','DisplayName','ASM');
%             ylabel('$u(x)$','Interpreter','latex');             
%             yyaxis right;
%             set(gca, 'YScale', 'log');
%             semilogy(x,eF,'k:','DisplayName','error');
%             %plot(x,imag(uF),'--','DisplayName','error_im');
%             ylim([10^-16,10^-14]);
%             xlabel('$x$','Interpreter','latex');
%             ylabel('$|error(x)|$','Interpreter','latex'); 
%             grid on; grid minor;
%             set(gca,'FontSize',16);
%             set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
%             % legend('exact','nm\_spec','nm\_error');  
%             axis tight; legend('Location','Best');
%             
%             subplot(2,1,2);
            yyaxis left;
%             plot(xfine,u1fine,'m','DisplayName','analytical'); hold on;
            plot(x,u1,'m','DisplayName','analytical'); hold on;
            plot(x,real(ux),'bh','DisplayName','GFS');
            %plot(x,imag(ux),'--','DisplayName','im(u\_ASM)');

            % Style options
            colors = lines(7);  % Pick a colormap with distinct colors
            markers = {'o', 's', '^', 'd', '*', 'x', '+', 'v', '>', '<'};
            lineStyles = {'-', '--', ':', '-.'};
            nitems = size(ux_compare,2);
            for i=1:nitems
                colorIdx = mod(i-1, size(colors, 1)) + 1;
                markerIdx = mod(i-1, length(markers)) + 1;
                styleIdx = mod(i-1, length(lineStyles)) + 1;
                plot(x,ux_compare(:,i),'Color', colors(colorIdx, :), ...
                     'Marker', markers{markerIdx}, ...
                     'LineStyle', lineStyles{styleIdx}, ...
                     'DisplayName',method(i)); hold on;
            end
            ylabel("$u^{'}(x)$",'Interpreter','latex'); 
            yyaxis right;
            set(gca, 'YScale', 'log');
            semilogy(x,abs(real(ux)-u1),'b:','DisplayName','error (GFS)');
            for i=1:nitems
                colorIdx = mod(i-1, size(colors, 1)) + 1;
                markerIdx = mod(i-1, length(markers)) + 1;
                styleIdx = mod(i-1, length(lineStyles)) + 1;
                str1 = sprintf('error (%s)',method(i));
                semilogy(x,abs(ux_compare(:,i)-u1),'Color', colors(colorIdx, :), ...
                      'Marker', 'none', ... %markers{markerIdx}, ...
                     'LineStyle', lineStyles{styleIdx}, ...
                     'DisplayName',str1);
            end
            ylim([10^-10,10^2]);
            xlabel('$x$','Interpreter','latex');
            ylabel("$|error(x)|$",'Interpreter','latex'); 
            grid on; grid minor;
            set(gca,'FontSize',16);
            set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
            % legend('exact','nm\_spec','nm\_error'); 
            axis tight; legend('Location','Best');
            if(fig_save)
                saveas(gcf,fig_up);
            end
%             saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/ASM_mode1.67.jpg');
%             saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/ASM_multimode.jpg');
%             saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/ASM_5modes.jpg');
%             saveas(gcf,'~/Dropbox/RapakaS_ASM/Figures/up_ramp_4.jpg');
%             saveas(gcf,'~/Dropbox/RapakaS_ASM/Figures/up_eax_2mode.jpg');
%             saveas(gcf,'~/Dropbox/RapakaS_ASM/Figures/up_Nk32_n3.jpg');
        end
    end

    function [ux,kn] = dudx_modASM(u,plotTag)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute jumps for the function & derivatives
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        uj = zeros(1,4*n_apmodes);
        switch(jump_method)
          case(1) % analytical expressions for jumps      
          switch(IC_profile)
              case(5)
                m = m_monomials;
                % exact jump conditions for the monomial u=x^m
                mcoeff = 1;
                for i=1:4*n_apmodes
                    uj(i) = ( (pi)^(m-i+1) - (-pi)^(m-i+1) )*mcoeff;
                    mcoeff= mcoeff*(m-i+1);
                end
              case(8)
                for i=1:4*n_apmodes
                    uj(i) = exp(pi)-exp(-pi);
                end
              case(9)
%                 th = 0; a = -2; b= 3/4; % * 2*pi;
                pmat = pascal(4*n_apmodes);
                for i=0:4*n_apmodes-1
                    tmp = 0;
                    for j=0:i
%                         pmat = pascal(i);
%                         minordiag = sort(find(fliplr(eye(size(pmat)))), 'descend');
%                         iCj = pmat(minordiag);  
                        iCj = pmat(i-j+1,j+1);
                        fb  = a^(i-j)*exp(a*( pi+th));
                        fa  = a^(i-j)*exp(a*(-pi+th));
                        if(mod(j,2)==0) % j even
                            gb  = (-1)^(j/2)*b^j*sin(b*( pi+th));
                            ga  = (-1)^(j/2)*b^j*sin(b*(-pi+th));
                        else            % j odd
                            gb  = (-1)^((j-1)/2)*b^j*cos(b*( pi+th));
                            ga  = (-1)^((j-1)/2)*b^j*cos(b*(-pi+th));
                        end
                        tmp = tmp + iCj*(fb*gb-fa*ga);
                    end
                    uj(i+1) = tmp;
                end
              case(10)
                uj(1) = log(pi+th)-log(-pi+th);
                uj(2) = (1./(pi+th) - 1./(-pi+th));
                fact  = 1;
                for i = 3:4*n_apmodes
                    fact  = fact * -(i-2);
                    uj(i) = (1./(pi+th).^(i-1) - 1./(-pi+th).^(i-1)) * fact;
                end
  
              otherwise
                for i=1:4*n_apmodes
                    if(mod(i,4)>=2)
                        sn = -1;
                    else
                        sn = 1;
                    end
                    uj(i) = 2*sum( sn* k.^(i-1).*sin(k.*pi),2);
                    if(IC_profile==2 && mod(i,2)==0) % pure sin waves
                        uj(i) = 0; % odd jumps are zero
                    end
                    if(IC_profile==3 && mod(i,2)==1) % pure cos waves
                        uj(i) = 0; % even jumps are zero
                    end
                end
                Jm = uj.';
                max_Jm = max(abs(Jm))
          end
                
                
          case(2) % finite difference approximation for jumps
            % r = 10;  % accuracy of each derivative
            for i=1:4*n_apmodes
                q = i-1;
                stw = r+q;
                % weights for q_th derivative with O(dx^r) accuracy
                fdq = classic_fd_stencil(0,stw-1,q).'; % 1:stw
                bdq = classic_fd_stencil(stw-1,0,q).'; % end-stw+1:end
                uj(i) = (bdq*u((end-(stw-1)):end) - fdq*u(1:stw))/dx^q;
%                 uj(i) = bdq * ( u((end-(stw-1)):end) + flip(u(1:stw)) )/dx^q;
            end
            uj(1) = u(end)-u(1);     
            
          otherwise % high resolution FD approximation for jumps
            stw = length(u);
            for i=1:4*n_apmodes
                q = i-1;
                % weights for q_th derivative with O(dx^r) accuracy
                fdq = high_resolution_stencil(0,stw-1,q,r).'; % 1:stw
                bdq = high_resolution_stencil(stw-1,0,q,r).'; % end-stw+1:end
                uj(i) = (bdq*u((end-(stw-1)):end) - fdq*u(1:stw))/dx^q;
            end
            uj(1) = u(end)-u(1);        
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        switch(n_apmodes)
        otherwise            
            %============ sin modes ==================
            cnt = 1;
            while(1)
                kmat = ones(ns);
                Jmat = zeros(ns);
                Jvec = zeros(ns,1);
                Jvec1= zeros(ns,1);
                eind = 2*(ns-1);
                for i=1:ns
                    sind = 2*i-1;
                    Jmat(i,:) = uj(sind+(0:2:eind));
                    Jvec(i)   = uj(sind+eind+2);
                    Jvec1(i)  = uj(sind);
                end
%                 Jmat_sin = Jmat
                rankJ = rank(Jmat);
                detJ = det(Jmat);
                fprintf('|J|_sin = %e\n',detJ);
                if( (abs(detJ)<TOL_DetJ || rankJ<ns) && check_rank==1 )
                    fprintf('rank(J) = %d<%d, ', rankJ,ns);                    
                    if(cnt>1)
                        ns = ns-1;
                    else
                        ns = rankJ; %ns-1;
                    end
                    fprintf('ns=%d (#sin modes reduced)\n', ns);
                    cnt = cnt + 1;
                else
                    break;
                end
            end
            evec = - (Jmat\Jvec);
            evec = [1; flip(evec)] .* ((-1).^(0:ns)).';
            ki = sqrt(roots(evec)).';   %row vector
            for i=1:ns
                kmat(i,:) = (-1)^(i+1) * ki.^(2*(i-1));
            end
            ui = 0.5* (kmat\Jvec1).' ./ sin(ki*pi);   %row vector 
            fprintf('\n');
%             ki
%             ki_simple = sqrt((Jvec./Jvec1).^(1/ns))
%             Jvec1
%             Jvec2 = Jvec./((-ki.^2).^ns).'
            %============ cos modes ==================
            cnt = 1;
            while(1)
                kmat = ones(nc);
                Jmat = zeros(nc);
                Jvec = zeros(nc,1);
                Jvec1= zeros(nc,1);
                eind = 2*(nc-1);
                for i=1:nc
                    sind = 2*i;
                    Jmat(i,:) = uj(sind+(0:2:eind));
                    Jvec(i)   = uj(sind+eind+2);
                    Jvec1(i)  = uj(sind);
                end
%                 Jmat_cos = Jmat
                rankJ = rank(Jmat);
                detJ = det(Jmat);
                fprintf('|J|_cos = %e\n',detJ);
                if( (abs(detJ)<TOL_DetJ || rankJ<nc) && check_rank==1 )
                    fprintf('rank(J) = %d<%d, ', rankJ,nc);
                    if(cnt>1)
                        nc = nc-1;
                    else
                        nc = rankJ; %nc-1;
                    end
                    fprintf('nc=%d (#cos modes reduced)\n', nc);
                    cnt = cnt + 1;
                else
                    break;
                end
            end
            fprintf('Jump vector: \n');
            uj.'
            evec = - (Jmat\Jvec);
            evec = [1; flip(evec)] .* ((-1).^(0:nc)).';
            kr = sqrt(roots(evec)).';   %row vector
            for i=1:nc
                kmat(i,:) = (-1)^(i+1) * kr.^(2*(i-1));
            end
            ur = -0.5* (kmat\Jvec1).' ./ (kr.*sin(kr*pi)); %row vector
            %=========================================
            
            ue = 0.5 * ( u + flip(u) );
            uo = 0.5 * ( u - flip(u) );
            ks = ki % row vector
            kc = kr % row vector

            %============ odd part ==================
            Vmat = ks.*x(1:ns);  % ks->row vector, x->column vector
            uvec = uo(1:ns);    % uvec->column vector
            Vmat = sin(Vmat);        
            ao   = Vmat\uvec   % ao->column vector
            cond_Vmat_sin = cond(Vmat)

            uao  = real( sin(ks.*x)*ao );
            uaxo = real( cos(ks.*x)*(ao.*ks.') );
            %============ even part ==================
            Vmat = kc.*x(1:nc);  % ks->row vector, x->column vector
            uvec = ue(1:nc);    % uvec->column vector
            Vmat = cos(Vmat);       
            ae   = Vmat\uvec   % ae->column vector
            cond_Vmat_cos = cond(Vmat)

            uae  = real( cos(kc.*x)*ae );
            uaxe = real(-sin(kc.*x)*(ae.*kc.') );
            %=========================================
            ua   = uae  + uao;
            uax  = uaxe + uaxo;
            up   = u-ua;
            
        end
        figure(3); clf; lw=1;
        plot(x,u,'r--','DisplayName','u'); hold on;
        plot(x,ue,'b.-','DisplayName','u_{even}'); hold on;
        plot(x,uae,'bo','DisplayName','u_{a,even}'); hold on;
        plot(x,uo,'c.-','DisplayName','u_{odd}'); hold on;
        plot(x,uao,'co','DisplayName','u_{a,odd}'); hold on;
        plot(x,ua,'m.-','DisplayName','u_{a}'); hold on;
        plot(x,up,'k.-','DisplayName','u_{p}'); hold on;
        xlabel('x');ylabel('u(x)'); grid on; grid minor;
        set(gca,'FontSize',16);
        set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
        axis tight; legend('Location','Best');

        figure(4); clf;
        lw = 1;
        plot(x,u,'m:','DisplayName','$u(x)$, original'); hold on;
        plot(x,ua,'r.-','DisplayName','$u_a(x)$, aperiodic');
        plot(x,up,'b.-','DisplayName','$u_p(x)$, periodic');
%         plot(x,ua+up,'g.-','DisplayName','$u_a(x)+u_p(x)$');
        %plot(x,imag(uF),'--','DisplayName','error_im');
        xlabel('$x$','Interpreter','latex');ylabel('$u(x)$','Interpreter','latex'); grid on; grid minor;
        set(gca,'FontSize',16); legend;
        set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
        % legend('exact','nm\_spec','nm\_error');  
        axis tight; 
        hl = legend('Location','NorthWest'); set(hl,'Interpreter','latex');
%         if(fig_save)
%             saveas(gcf,fig_u);
%         end
        
        M = floor(N/2)*2;
        kn = 2*pi/L * (-M/2:M/2-1).'; % modes        
%         [upx,uF] = dudx_FFT(up,L,0);
        
        isperiodic = true;
%         isperiodic = false;
%         fprintf('isperiodic = %d\n', isperiodic);
        uk = nmFFT(up,kn,isperiodic);   % FFT through non-integer basis
        % uk(M/2+1)=0;    % nullify odd-ball wave number    
        upx = nmiFFT(1i*kn.*uk,kn,isperiodic);  % iFFT
        upx = real(upx);
        % copy of the original function via iFFT
        uF = nmiFFT(uk,kn,isperiodic);  % iFFT
        uF = real(uF)+ua;
        eF = uF-u;
        
        ux = uax + upx;        

        if(plotTag)            
            lw=1;
            e0_u_N(nn,1) = max(abs(eF));
            fprintf('||u-u_mASM||_2  : %.2e\n',norm(abs(eF),2)*dx^(1/2));
            fprintf('||u-u_mASM||_Inf: %.2e\n\n',max(abs(eF)));
%             fprintf("Nyq_comp: \nu : |  uk| = %10.2e \nu': |k*uk| = %10.2e\n\n",abs(uk(1)),abs(kn(1)*uk(1)));
            
            figure(14); clf;    
%             subplot(3,1,1);
%             plot((kn),real((uk)),'.-','DisplayName','real'); hold on;
%             plot((kn),imag((uk)),'-','DisplayName','imag');
% %             max_uk_im = max(abs(imag(uk)))
%             xlabel('$k_n$','Interpreter','Latex'); ylabel('$\hat{u}_{p,n}$','Interpreter','Latex');
%             legend('Location','Best'); grid on; grid minor;
%             set(gca,'FontSize',16); axis tight;
%             xlim([-N/2,N/2]);
%             set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
%             % set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 1);

            subplot(2,1,1);
            plot(x,u,'ro','DisplayName','exact'); hold on;
            plot(x,real(uF),'b-','DisplayName','mASM');
            plot(x,eF,'k:','DisplayName','error');
            %plot(x,imag(uF),'--','DisplayName','error_im');
            xlabel('$x$','Interpreter','latex');ylabel('$u(x)$','Interpreter','latex'); grid on; grid minor;
            set(gca,'FontSize',16);
            set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
            % legend('exact','nm\_spec','nm\_error');  
            axis tight; legend('Location','Best');
            
            subplot(2,1,2);
            plot(x,u1,'ro','DisplayName','exact'); hold on;
            plot(x,real(ux),'b-','DisplayName','mASM');
            plot(x,real(ux)-u1,'k:','DisplayName','error');
            %plot(x,imag(ux),'--','DisplayName','im(u\_ASM)');
            xlabel('$x$','Interpreter','latex');ylabel("$u^{'}(x)$",'Interpreter','latex'); grid on; grid minor;
            set(gca,'FontSize',16);
            set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
            % legend('exact','nm\_spec','nm\_error'); 
            axis tight; legend('Location','Best');
%             if(fig_save)
%                 saveas(gcf,fig_up);
%             end
        end

    end

    function ux = dudx_modASM_randModes(u,plotTag)
        ue = 0.5 * ( u + flip(u) );
        uo = 0.5 * ( u - flip(u) );
        nc = n_apmodes;
        ns = n_apmodes;
        ks = linspace(1,ns,ns)+0.15; % row vector
        kc = linspace(1,nc,nc)+0.15; % row vector
        ks = rand(1,ns); % row vector
        kc = rand(1,nc); % row vector
        
        %============ odd part ==================
        Vmat = ks.*x(1:ns);  % ks->row vector, x->column vector
        uvec = uo(1:ns);    % uvec->column vector
        Vmat = sin(Vmat);        
        ao   = Vmat\uvec;   % ao->column vector
        cond_Vmat_sin = cond(Vmat)
        
        uao  = real( sin(ks.*x)*ao );
        uaxo = real( cos(ks.*x)*(ao.*ks.') );
        %============ even part ==================
        Vmat = kc.*x(1:nc);  % ks->row vector, x->column vector
        uvec = ue(1:nc);    % uvec->column vector
        Vmat = cos(Vmat);       
        ae   = Vmat\uvec;   % ae->column vector
        cond_Vmat_cos = cond(Vmat)
        ks.'
        kc.'
        ao
        ae

        uae  = real( cos(kc.*x)*ae );
        uaxe = real(-sin(kc.*x)*(ae.*kc.') );
        %=========================================
        ua   = uae  + uao;
        uax  = uaxe + uaxo;
        up   = u-ua;

        figure(3); clf; lw=1;
        plot(x,u,'r--','DisplayName','u'); hold on;
        plot(x,ue,'b.-','DisplayName','u_{even}'); hold on;
        plot(x,uae,'bo','DisplayName','u_{a,even}'); hold on;
        plot(x,uo,'c.-','DisplayName','u_{odd}'); hold on;
        plot(x,uao,'co','DisplayName','u_{a,odd}'); hold on;
        plot(x,ua,'m.-','DisplayName','u_{a}'); hold on;
        plot(x,up,'k.-','DisplayName','u_{p}'); hold on;
        xlabel('x');ylabel('u(x)'); grid on; grid minor;
        set(gca,'FontSize',16);
        set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
        axis tight; legend('Location','Best');

        figure(4); clf;
        lw = 1;
        plot(x,u,'m:','DisplayName','$u(x)$, original'); hold on;
        plot(x,ua,'r.-','DisplayName','$u_a(x)$, aperiodic');
        plot(x,up,'b.-','DisplayName','$u_p(x)$, periodic');
        %plot(x,imag(uF),'--','DisplayName','error_im');
        xlabel('$x$','Interpreter','latex');ylabel('$u(x)$','Interpreter','latex'); grid on; grid minor;
        set(gca,'FontSize',16); legend;
        set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
        % legend('exact','nm\_spec','nm\_error');  
        axis tight; 
        hl = legend('Location','NorthWest'); set(hl,'Interpreter','latex');
%         if(fig_save)
%             saveas(gcf,fig_u);
%         end
        
        M = floor(N/2)*2;
        kn = 2*pi/L * (-M/2:M/2-1).'; % modes        
%         [upx,uF] = dudx_FFT(up,L,0);
        
        isperiodic = true;
%         isperiodic = false;
%         fprintf('isperiodic = %d\n', isperiodic);
        uk = nmFFT(up,kn,isperiodic);   % FFT through non-integer basis
        % uk(M/2+1)=0;    % nullify odd-ball wave number    
        upx = nmiFFT(1i*kn.*uk,kn,isperiodic);  % iFFT
        upx = real(upx);
        % copy of the original function via iFFT
        uF = nmiFFT(uk,kn,isperiodic);  % iFFT
        uF = real(uF)+ua;
        eF = uF-u;
        
        ux = uax + upx;
        

        if(plotTag)            
            lw=1;
            e0_u_N(nn,1) = max(abs(eF));
            fprintf('||u-u_mASM||_2  : %.2e\n',norm(abs(eF),2)*dx^(1/2));
            fprintf('||u-u_mASM||_Inf: %.2e\n\n',max(abs(eF)));
%             fprintf("Nyq_comp: \nu : |  uk| = %10.2e \nu': |k*uk| = %10.2e\n\n",abs(uk(1)),abs(kn(1)*uk(1)));
            
            figure(14); clf;    
%             subplot(3,1,1);
%             plot((kn),real((uk)),'.-','DisplayName','real'); hold on;
%             plot((kn),imag((uk)),'-','DisplayName','imag');
% %             max_uk_im = max(abs(imag(uk)))
%             xlabel('$k_n$','Interpreter','Latex'); ylabel('$\hat{u}_{p,n}$','Interpreter','Latex');
%             legend('Location','Best'); grid on; grid minor;
%             set(gca,'FontSize',16); axis tight;
%             xlim([-N/2,N/2]);
%             set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
%             % set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 1);

            subplot(2,1,1);
            plot(x,u,'ro','DisplayName','exact'); hold on;
            plot(x,real(uF),'b-','DisplayName','mASM');
            plot(x,eF,'k:','DisplayName','error');
            %plot(x,imag(uF),'--','DisplayName','error_im');
            xlabel('$x$','Interpreter','latex');ylabel('$u(x)$','Interpreter','latex'); grid on; grid minor;
            set(gca,'FontSize',16);
            set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
            % legend('exact','nm\_spec','nm\_error');  
            axis tight; legend('Location','Best');
            
            subplot(2,1,2);
            plot(x,u1,'ro','DisplayName','exact'); hold on;
            plot(x,real(ux),'b-','DisplayName','mASM');
            plot(x,real(ux)-u1,'k:','DisplayName','error');
            %plot(x,imag(ux),'--','DisplayName','im(u\_ASM)');
            xlabel('$x$','Interpreter','latex');ylabel("$u^{'}(x)$",'Interpreter','latex'); grid on; grid minor;
            set(gca,'FontSize',16);
            set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
            % legend('exact','nm\_spec','nm\_error'); 
            axis tight; legend('Location','Best');
%             if(fig_save)
%                 saveas(gcf,fig_up);
%             end
        end
    end
    function ux = dudx_modASM_poly(u,plotTag)
        ue = 0.5 * ( u + flip(u) );
        uo = 0.5 * ( u - flip(u) );
        nc = n_apmodes;
        ns = n_apmodes;
        %============ odd part ==================
        Vmat = ones(ns);
        uvec = zeros(ns,1);
        for i=1:ns
            Vmat(i,1) = x(i);
            uvec(i)   = uo(i);
        end
        for i=2:ns
            Vmat(:,i) = Vmat(:,i-1) .* (Vmat(:,1).^2);
        end
        ao = Vmat\uvec;
        uao   = zeros(size(x));
        uaxo  = zeros(size(x));
        for i=1:ns
            uao   = uao   + ao(i) * x.^(2*i-1);
            uaxo  = uaxo  + ao(i) * x.^(2*i-2) * (2*i-1);
        end
        %============ even part ==================
        Vmat = ones(nc);
        uvec = zeros(nc,1);
        for i=1:nc
            Vmat(i,2) = x(i)^2;
            uvec(i)   = ue(i);
        end
        for i=3:nc
            Vmat(:,i) = Vmat(:,i-1) .* Vmat(:,2);
        end
        ae = Vmat\uvec;

        uae   = zeros(size(x));
        uaxe  = zeros(size(x));
        for i=1:nc
            uae   = uae + ae(i) * x.^(2*i-2);
        end
        for i=2:nc
            uaxe  = uaxe + ae(i) * x.^(2*i-3) * (2*i-2);
        end
        %=========================================
        ua   = uae  + uao;
        uax  = uaxe + uaxo;
        up   = u-ua;

        figure(3); clf; lw=1;
        plot(x,u,'r--','DisplayName','u'); hold on;
        plot(x,ue,'b.-','DisplayName','u_{even}'); hold on;
        plot(x,uae,'bo','DisplayName','u_{a,even}'); hold on;
        plot(x,uo,'c.-','DisplayName','u_{odd}'); hold on;
        plot(x,uao,'co','DisplayName','u_{a,odd}'); hold on;
        plot(x,ua,'m.-','DisplayName','u_{a}'); hold on;
        plot(x,up,'k.-','DisplayName','u_{p}'); hold on;
        xlabel('x');ylabel('u(x)'); grid on; grid minor;
        set(gca,'FontSize',16);
        set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
        axis tight; legend('Location','Best');

        figure(4); clf;
        lw = 1;
        plot(x,u,'m:','DisplayName','$u(x)$, original'); hold on;
        plot(x,ua,'r.-','DisplayName','$u_a(x)$, aperiodic');
        plot(x,up,'b.-','DisplayName','$u_p(x)$, periodic');
        %plot(x,imag(uF),'--','DisplayName','error_im');
        xlabel('$x$','Interpreter','latex');ylabel('$u(x)$','Interpreter','latex'); grid on; grid minor;
        set(gca,'FontSize',16); legend;
        set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
        % legend('exact','nm\_spec','nm\_error');  
        axis tight; 
        hl = legend('Location','NorthWest'); set(hl,'Interpreter','latex');
%         if(fig_save)
%             saveas(gcf,fig_u);
%         end
        
        M = floor(N/2)*2;
        kn = 2*pi/L * (-M/2:M/2-1).'; % modes        
%         [upx,uF] = dudx_FFT(up,L,0);
        
        isperiodic = true;
%         isperiodic = false;
%         fprintf('isperiodic = %d\n', isperiodic);
        uk = nmFFT(up,kn,isperiodic);   % FFT through non-integer basis
        % uk(M/2+1)=0;    % nullify odd-ball wave number    
        upx = nmiFFT(1i*kn.*uk,kn,isperiodic);  % iFFT
        upx = real(upx);
        % copy of the original function via iFFT
        uF = nmiFFT(uk,kn,isperiodic);  % iFFT
        uF = real(uF)+ua;
        eF = uF-u;
        
        ux = uax + upx;
        

        if(plotTag)            
            lw=1;
            e0_u_N(nn,1) = max(abs(eF));
            fprintf('||u-u_mASM||_2  : %.2e\n',norm(abs(eF),2)*dx^(1/2));
            fprintf('||u-u_mASM||_Inf: %.2e\n\n',max(abs(eF)));
%             fprintf("Nyq_comp: \nu : |  uk| = %10.2e \nu': |k*uk| = %10.2e\n\n",abs(uk(1)),abs(kn(1)*uk(1)));
            
            figure(14); clf;    
%             subplot(3,1,1);
%             plot((kn),real((uk)),'.-','DisplayName','real'); hold on;
%             plot((kn),imag((uk)),'-','DisplayName','imag');
% %             max_uk_im = max(abs(imag(uk)))
%             xlabel('$k_n$','Interpreter','Latex'); ylabel('$\hat{u}_{p,n}$','Interpreter','Latex');
%             legend('Location','Best'); grid on; grid minor;
%             set(gca,'FontSize',16); axis tight;
%             xlim([-N/2,N/2]);
%             set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
%             % set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 1);

            subplot(2,1,1);
            plot(x,u,'ro','DisplayName','exact'); hold on;
            plot(x,real(uF),'b-','DisplayName','mASM');
            plot(x,eF,'k:','DisplayName','error');
            %plot(x,imag(uF),'--','DisplayName','error_im');
            xlabel('$x$','Interpreter','latex');ylabel('$u(x)$','Interpreter','latex'); grid on; grid minor;
            set(gca,'FontSize',16);
            set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
            % legend('exact','nm\_spec','nm\_error');  
            axis tight; legend('Location','Best');
            
            subplot(2,1,2);
            plot(x,u1,'ro','DisplayName','exact'); hold on;
            plot(x,real(ux),'b-','DisplayName','mASM');
            plot(x,real(ux)-u1,'k:','DisplayName','error');
            %plot(x,imag(ux),'--','DisplayName','im(u\_ASM)');
            xlabel('$x$','Interpreter','latex');ylabel("$u^{'}(x)$",'Interpreter','latex'); grid on; grid minor;
            set(gca,'FontSize',16);
            set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
            % legend('exact','nm\_spec','nm\_error'); 
            axis tight; legend('Location','Best');
%             if(fig_save)
%                 saveas(gcf,fig_up);
%             end
        end
    end


    %% routines to compute derivatives
    function derivatives = dudx_poly(x, y, m, order)
        %{
        Compute derivatives of a function using polynomial approximation.
        Arguments:
            x: Vector of x-values (independent variable).
            y: Vector of y-values (dependent variable).
            m: Degree of the polynomial to fit.
            order: Order of the derivative to compute.
    
        Returns:
            derivatives: Vector of computed derivatives at x-values.
        %}
        
        % Validate input sizes
        if length(x) ~= length(y)
            error('x and y must be vectors of the same length.');
        end
        
        % Fit a polynomial of degree m to the data
        poly_coeff = polyfit(x, y, m);
        
        % Compute the derivative coefficients for the desired order
        for i = 1:order
            poly_coeff = polyder(poly_coeff);
        end
        
        % Evaluate the derivative at the x points
        derivatives = polyval(poly_coeff, x);
    end
    
    % first derivative using finite difference (CD for periodic domain)
    function f = dudx(u,dx,sw,wts)
        % ML=1; wts = 0.5*[-1 0 1]; % must be a row vector
        M = size(u,1);        
        M = N;
        f = zeros(M,1);
        ind = circshift(1:M,sw);
        for w = wts
            f = f + w*u(ind);
            ind = circshift(ind,-1);
        end
        f = 1/dx*f;
        f = [f; f(1)];
    end

    % standard finite difference: qth derivative with O(dx^r) accuracy 
    % central difference in the interior and one sided stencil at the boundaries
    function f = dudx_FD(u,r,dx)
        q = 1;
        stw = r+q;  % stw = M1+M2+1 for the first derivative (q=1)

        M = size(u,1);        
        f = zeros(M,1);
        % weights for q_th derivative with O(dx^r) accuracy
        M1 = floor((stw-1)/2); M2 = M1;                 % central differences
        wcd = classic_fd_stencil(M1,M2,q).';     % -(stw-1)/2:0:(stw-1)/2
        for i=M1+1:M-M2
            f(i) = wcd*u(i-M1:i+M2); 
        end
        
        for i=1:M1
            M1 = i-1;         M2 = stw-i;        % forward differences
            wfd = classic_fd_stencil(M1,M2,q).'; % 1:stw
            f(i) = wfd*u(i-M1:i+M2);
        end
        
        for i=M-M2+1:M
            M2 = M-i;         M1 = stw-1-M2;     % backward differences
            wbd = classic_fd_stencil(M1,M2,q).'; % end-stw+1:end
            f(i) = wbd*u(i-M1:i+M2);
        end
        
        f = f/dx^q;
    end

    % spectral derivative
    function [ux,uF] = dudx_FFT(u,L,plotTag)
        M = floor(N/2)*2;
        kn = 2*pi/L * ifftshift(-M/2:M/2-1).'; % modes
%         kn = 2*pi/L * ifftshift(-M/2:M/2).'; % modes
%         kn(1) = 0;
        % use this for even derivatives
        % kn = 2*pi/L * [0:M/2-1,-M/2,-M/2+1:-1].'; 
        % use this for odd derivatives: nullify odd-ball wave number
        % kn = 2*pi/L * [0:M/2-1, 0  ,-M/2+1:-1].'; 

        uk = fft(ifftshift(u(1:M)));   % FFT: ifftshift rearranges x=-L/2:L/2-dx to x=0:L-dx
        % uk(M/2+1)=0;    % nullify odd-ball wave number
        ux = fftshift(ifft(1i*kn.*uk));  % iFFT
        ux = [ux;ux(1)];
        ux = real(ux);
        uF = fftshift(ifft(uk));  % iFFT
        uF = [uF;uF(1)];
        
        if(plotTag)
            % copy of the original function via iFFT
            uF = fftshift(ifft(uk));  % iFFT
            uF = [uF;uF(1)];
            uF = real(uF);
            eF = uF-u;
%         u
%         uF
%             kn2 = -M/2:M/2-1;
%             mzero = M/2+1;
%             kn2(mzero+1) = 1.35;
%             kn2(mzero-1) = -kn2(mzero+1);
%             %kn2(M) = kn2(M)-0.5;
%             uk2 = nmFFT(u,kn2);
%             uF2 = nmiFFT(uk2,kn2);
%             eF2 = uF2-u;
            e0_u_N(nn,2) = max(abs(eF(1:N)));
            fprintf('max|u-u_FFT|: %.2e\n',max(abs(eF(1:N+1))));
%             fprintf('max|u-u_nmFFT|: %.2e\n\n',max(abs(eF2)));
            figure(12); clf;    
            subplot(3,1,1);
            plot(fftshift(kn),real(fftshift(uk))/M,'.-','DisplayName','real'); hold on;
            plot(fftshift(kn),imag(fftshift(uk))/M,'-','DisplayName','imag');
%             plot((kn2),real((uk2)),'.-','DisplayName','nm_r'); hold on;
%             plot((kn2),imag((uk2)),'-','DisplayName','nm_{im}');
            xlabel('$k_n$','Interpreter','Latex'); ylabel('$\hat{u}_n$','Interpreter','Latex');
%             legend('real','imag','nm_r','nm_{im}'); 
            legend('Location','NorthEast'); grid on; grid minor;
            set(gca,'FontSize',16); axis tight;
            set(findall(gca, 'Type', 'Line'),'LineWidth',1);
            % set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 1);

            subplot(3,1,2);
            plot(x,u,'ro','DisplayName','exact'); hold on;
            plot(x,real(uF),'b-','DisplayName','FFT');
            plot(x,eF,'k:','DisplayName','error');
%             plot(x,real(uF2),'^-','DisplayName','nm\_spec');
%             plot(x,eF2,'.-','DisplayName','nm\_error');
            %plot(x,imag(uF),'--','DisplayName','error_imag');
            xlabel('x');ylabel('u(x)'); grid on; grid minor;
            set(gca,'FontSize',16);
            set(findall(gca, 'Type', 'Line'),'LineWidth',1);
            legend('exact','FFT','error');  
%             legend('exact','FFT','error','nm\_spec','nm\_err');  
            axis tight; %legend('Location','Best');
            
            subplot(3,1,3);
            plot(x,u1,'ro','DisplayName','exact'); hold on;
            plot(x,real(ux),'b-','DisplayName','FFT');
            plot(x,real(ux)-u1,'k:','DisplayName','error');
            %plot(x,imag(ux),'--');
            xlabel('x');ylabel('du/dx'); grid on; grid minor;
            set(gca,'FontSize',16);
            set(findall(gca, 'Type', 'Line'),'LineWidth',1);
            legend('exact','FFT','error'); 
            axis tight; %legend('Location','Best');
%             saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/FFT_mode1.67.jpg');
%             saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/FFT_multimode.jpg');
%             saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/FFT_5modes.jpg');
        end
    end

    % spectral derivative through non-integer modes
    function ux = dudx_ASM(u,L,kn,plotTag)
        % kn = find_modes(u,0);
        % M = floor(N/2)*2;
        % kn = 2*pi/L * (-M/2:M/2-1).'; % modes
                        
        % use this for even derivatives
        % kn = 2*pi/L * [0:M/2-1,-M/2,-M/2+1:-1].'; 
        % use this for odd derivatives: nullify odd-ball wave number
        % kn = 2*pi/L * [0:M/2-1, 0  ,-M/2+1:-1].'; 

        isperiodic = true;
        isperiodic = false;
        % fprintf('isperiodic = %d\n', isperiodic);
        uk = nmFFT(u,kn,isperiodic);   % FFT through non-integer basis
        % uk(M/2+1)=0;    % nullify odd-ball wave number    
        ux = nmiFFT(1i*kn.*uk,kn,isperiodic);  % iFFT
        ux = real(ux);
        % copy of the original function via iFFT
        uF = nmiFFT(uk,kn,isperiodic);  % iFFT
        uF = real(uF);
        eF = uF-u;

        if(plotTag)          
            lw=1;
            e0_u_N(nn,1) = max(abs(eF));
            fprintf('max|u-u_ASM|: %.2e\n\n',max(abs(eF)));
%             fprintf("Nyq_comp: \nu : |  uk| = %10.2e \nu': |k*uk| = %10.2e\n\n",abs(uk(1)),abs(kn(1)*uk(1)));
            
            figure(13); clf;    
            subplot(3,1,1);
            plot((kn),real((uk)),'.-'); hold on;
            plot((kn),imag((uk)),'-');
%             max_uk_im = max(abs(imag(uk)))
            xlabel('$k_n$','Interpreter','Latex'); ylabel('$\hat{u}_n$','Interpreter','Latex');
            legend('real','imag','nm_r','nm_{im}'); grid on; grid minor;
            set(gca,'FontSize',16); axis tight;
            xlim([-N/2,N/2]);
            set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
            % set(findobj(gca, 'Type', 'Line', 'Linestyle', '--'), 'LineWidth', 1);

            subplot(3,1,2);
            plot(x,u,'ro','DisplayName','exact'); hold on;
            plot(x,real(uF),'b-','DisplayName','ASM');
            plot(x,eF,'k:','DisplayName','error');
            %plot(x,imag(uF),'--','DisplayName','error_im');
            xlabel('x');ylabel('u(x)'); grid on; grid minor;
            set(gca,'FontSize',16);
            set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
            % legend('exact','nm\_spec','nm\_error');  
            axis tight; legend('Location','Best');
            
            subplot(3,1,3);
            plot(x,u1,'ro','DisplayName','exact'); hold on;
            plot(x,real(ux),'b-','DisplayName','ASM');
            plot(x,real(ux)-u1,'k:','DisplayName','error');
            %plot(x,imag(ux),'--','DisplayName','im(u\_ASM)');
            xlabel('x');ylabel('du/dx'); grid on; grid minor;
            set(gca,'FontSize',16);
            set(findall(gca, 'Type', 'Line'),'LineWidth',lw);
            % legend('exact','nm\_spec','nm\_error'); 
            axis tight; legend('Location','Best');
%             saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/ASM_mode1.67.jpg');
%             saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/ASM_multimode.jpg');
%             saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/ASM_2modes.jpg');
%             saveas(gcf,'~/Dropbox/RapakaS_FFDSR/Figures/ASM_multimodeFullRank.jpg');
        end
    end

    %% fully spectral derivative of the non-linear term
    function [uux_k, kn_uu] = uux_ASM(uk,kn)
        M = floor(N/2)*2;
        Nk = numel(kn);
                
        Nuuk = min(nchoosek(Nk,1)*2,M);
        uux_k  = zeros(Nuuk,1);
        kn_uu= zeros(size(uux_k));
        % uk(M/2+1)=0;    % nullify odd-ball wave number    
        mode = 0;
        TOL = 10^-15;
        for i=1:Nk         % kn must be col vector
            for j=i:Nk
                ksum = kn(i)+kn(j);
%                 if( abs(ksum) <= M/2 )
                if( -M/2<=ksum && ksum < M/2 )
                    u2xk = -1i*ksum*0.5* uk(i)*uk(j);
                    if( abs(u2xk) >= TOL)
                        mode = mode+1;
                        kn_uu(mode) = ksum;
                        if(i==j)
                            uux_k(mode) = u2xk;
                        else
                            uux_k(mode) = u2xk * 2;
                        end
                    end
                end
            end
        end    
        if(isrow(kn_uu))
            kn_uu = kn_uu.';
        end
        if(isrow(uux_k))
            uux_k = uux_k.';
        end
        
%         % find unique modes (results in sorted order)
%         nmodes = numel(kn_uu);  % num of original modes
%         [kn_uu,IA,IC] = unique(kn_uu);
%         % gather amplitudes of unique modes
%         uu_x_k = zeros(size(kn_uu));
%         for i=1:nmodes
%             uu_x_k(IC(i)) = uu_x_k(IC(i)) + uux_k(i);
%         end
%         uux_k = uu_x_k;
        [uux_k, kn_uu] = merge_modes(uux_k,kn_uu);
        
%         kn_uu
%         uux_k        
    end

%     function [uk,kn] = merge_modes(u1,k1,u2,k2)
%         uk = [u1;u2];   % u1, u2 must be col vec
%         kn = [k1;k2];   % k1, k2 must be col vec        
    function [uk,kn] = merge_modes(uk,kn)
        % % sort the modes
        % [kn,I] = sort(kn);
        % uk = uk(I);

        % find unique modes (results in sorted order)
        nmodes = numel(kn);  % num of original modes
        [kn,IA,IC] = unique(kn);
        % gather amplitudes of unique modes
        uk2 = zeros(size(kn));
        for i=1:nmodes
            uk2(IC(i)) = uk2(IC(i)) + uk(i);
        end
        uk = uk2;
        
%       if(Tag_merge_modes_eps)
%           [uk,kn] = merge_modes_eps(uk,kn,mode_eps);
%       end
    end

    function [uk,kn] = merge_modes_eps(uk,kn,mode_eps)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % merge modes within 'mode_eps' distance
        if(~issorted(kn))
            [uk,kn] = sort_modes(uk,kn);
        end
        nmodes = numel(kn);  % num of original modes
        ind_k1 = ceil(nmodes/2)+1;
        zero_mode = floor(nmodes/2)+1;
        if( kn(zero_mode)==0 )
            uk0 = uk(zero_mode);
        else
            uk0 = inf;
        end
        i = ind_k1;
        while(i < nmodes)
            cnt = 0;
            for j=i+1:nmodes
                if( abs(kn(j)-kn(i)) <= mode_eps && (kn(i)~=0 && kn(j)~=0) )
                    cnt = cnt+1;
                else
                    break;
                end
            end
            if(cnt > 0)
                ksum = mean(kn(i:i+cnt));
                uksum = sum(uk(i:i+cnt));
                kn(i) = ksum;
                uk(i) = uksum;
                kn(i+1:i+cnt) = [];
                uk(i+1:i+cnt) = [];
                nmodes = nmodes-cnt;
            end
            i=i+1;
            % nmodes = numel(kn);
        end
        % fill negative modes: complex conjugate of positive modes
        if( isfinite(uk0) )
            kn = [-flip(kn(ind_k1:end)); 0; kn(ind_k1:end)];
            uk = [ conj(flip(uk(ind_k1:end))); uk0; uk(ind_k1:end)];
        else
            kn = [-flip(kn(ind_k1:end)); kn(ind_k1:end)];
            uk = [ conj(flip(uk(ind_k1:end))); uk(ind_k1:end)];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    function [uk,kn] = mode_redistribute(uk,kn,x)
        Nx = numel(x);
        Nk = numel(kn);
        if(Nk <= Nx)
            return;
        end
        
%         if(~issorted(kn))
%             [uk,kn] = sort_amp(uk,kn);
            [uk,kn] = sort_modes(uk,kn);
%             [uk,kn] = sort_modes_abs(uk,kn);
%         end

kn
        % special mode selection
        ind = -ones(Nx,1);
        indC = -ones(Nk-Nx,1);
        kn1 = zeros(Nx,1);
        kn2 = zeros(Nk-Nx,1);
        uk1 = zeros(Nx,1);
        uk2 = zeros(Nk-Nx,1);
        
        k_Nyq = floor(Nx/2);
        kref = -k_Nyq;
        cnt = 1;
        mode_found = false;
        delta = 1;
        for j=1:Nk
%             if( kref==0 )
                if( kn(j)==0 )
                    ind(cnt) = j;
                    kn1(cnt) = kn(j);
                    uk1(cnt) = uk(j);
                    kref = kref+delta;
                    cnt = cnt+1;
                    mode_found = false;
%                 elseif( kref-0.5<=kn(j) && kn(j)<=kref+0.5)
%                     ind(cnt) = j;
%                     kn1(cnt) = kn(j);
%                     uk1(cnt) = uk(j);
%                     % kref = kref+1;
%                     cnt = cnt+1;
%                     mode_found = false;
%                 end
            elseif( kref<=kn(j) && kn(j)<kref+delta)
                if(mode_found)
                    if( abs(uk(j)) > abs(uk(ind(cnt))) )
                        ind(cnt) = j;
                        kn1(cnt) = kn(j);
                        uk1(cnt) = uk(j);
                    end
                else
                    ind(cnt) = j;
                    kn1(cnt) = kn(j);
                    uk1(cnt) = uk(j);
                    mode_found = true;
                end
            else
                if(mode_found==false)
                    kn1(cnt) = kref;
                    uk1(cnt) = 0;
                end
                if( -1<= kref && kref< 1)
                    delta = 0.5;
                else
                    delta = 1;
                end
                kref = kref+delta;
                cnt = cnt+1;
                mode_found = false;
            end
            if( kref >= k_Nyq )
                kref
                break;
            end
        end
        
        for kref = floor(kn(j)):k_Nyq-1
            kn1(cnt) = kref;
            uk1(cnt) = 0;
            cnt = cnt+1;
        end
%kn(ind(ind~=-1))
kn1

        if( isrow(kn) )
            kn = kn.';
        end
%         if(isperiodic)
%             E = exp(1i* kn.' .*x(1:N)); % periodic case
%         else
            % aperiodic case: i=1:N+1 (includes both end points +/- pi)
            E = exp(1i* kn(1:Nx).' .*x); % kn.'->row vec, x->col vec
            E1 = exp(1i* kn(Nx+1:end).' .*x); % kn.'->row vec, x->col vec
%         end
        uk1 = E \ (E1*uk(Nx+1:end));
        uk = uk(1:Nx)+uk1;
        kn = kn(1:Nx);
        kn
    end

    % plot solution and error
    function plot_sol(u)
        figure(1); clf;
        lw=1;
        subplot(2,1,1); hold on;
        plot(x,u0,'k--','LineWidth',lw);
        plot(x,u_ex,'ro','LineWidth',lw); 
        plot(x,u,'b.-','LineWidth',lw);
        legend('\phi_0','\phi_{ex}','\phi');
        xlabel('$x$','Interpreter','Latex'); ylabel('$\phi(x,t)$','Interpreter','Latex');
        axis tight;
%         xlim([0 L]);
        grid on;    grid minor; 
        set(gca,'fontsize',16);

        subplot(2,1,2);
        plot(x,error,'k.-'); 
        xlabel('$x$','Interpreter','Latex'); ylabel('$e(x,t)$','Interpreter','Latex');
        xlim([-L/2 L/2]);
        grid on;    grid minor; 
        set(gca,'fontsize',16);
    end    
end

%% ODE solvers: time_DM = 1-Euler, 3-RK3, 4-RK4 
function [tn,u] = ode_solve(func,T,u,time_DM,dx,CFL,max_steps)
    %CFL = 1; %max_steps = 2000;
    t=0; n=0; tn = zeros(max_steps,1);
    while(t<=T && n<max_steps)
        dt = abs(CFL*dx/max(abs(u)));
        dt = min(dt,10^-4);
        dt = 0.01;
        n=n+1;        
        t=t+dt;
        tn(n) = t;
        fprintf('n = %d\t dt = %f\t t = %f\n', n, dt, t);
        switch(time_DM)
            case(1) % Euler explicit time marching scheme
                un = u + dt*func(u);

            case(3) % RK3 explicit time marching scheme
                % RK step 1
                rh = func(u);
                un = u + 8/15*dt*rh;

                % RK step 2
                rh2 = func(un);
                un = u + dt*( 1/4*rh + 5/12*rh2 );

                % RK step 3
                rh2 = func(un);
                un = u + dt*( 1/4*rh + 3/4*rh2 );

            case(4) % RK4 explicit time marching scheme
                % RK step 1
                rh = func(u);
                un = u + 1/2*dt*rh;
                utemp = -u + un;

                % RK step 2
                rh = func(un);
                un = u + 1/2*dt*rh;
                utemp = utemp + 2*un;

                % RK step 3
                rh = func(un);
                un = u + dt*rh;
                utemp = utemp + un;

                % RK step 4
                rh = func(un);
                un = 1/3*utemp + 1/6*dt*rh;
        end
        u = un;
    end
    tn(n+1:end) = [];
end

function f = dudx_CD2(u,dx)
    N = size(u,1);        
    % CD2
    f = zeros(size(u));
%     for i=1:N
%         iR = i+1;
%         iL = i-1;
%         if(i==1) 
%             iL=iL+N;
%         elseif(i==N)
%             iR=iR-N;
%         end
%         f(i) = u(iR)-u(iL);
%     end
    f(1) = u(2)-u(N);
    f(2:N-1) = u(3:N)-u(1:N-2);
    f(N) = u(1)-u(N-1);
    f = 1/(2*dx)*f;
end

