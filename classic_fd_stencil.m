function [w] = classic_fd_stencil(M1,M2,q)
format short;
format long;
% k: formal order of accuracy
% M1=40; M2=20; k=5; kc = 2.95; %.75*pi; 
% M1=30; M2=M1-15; k=4; kc = 0.9*pi; %.75*pi; 
% % M1=40; M2=20; k=4; kc = 3.05; %.75*pi; 
% % M1=3; M2=3; k=4; kc = 0.5*pi; %.75*pi;  % Tam and Webb, JCP 1993
% M1=9; M2=8; k=4; kc = 2/3*pi; 
% M1=0; M2=6; k=2; kc = pi; 

    sw   = M1+M2+1;   % total stencil width and number of unknowns
    kmax = M1+M2;     % Max. possible accuracy
%     p    = (M1+M2+1)-(k+1)    % Number of free parameters
    
    % classic FD schemes
    kc = pi;
    p = 0;
    k = M1+M2;
    r = k-q+1; % order of accuracy
%     fprintf('Derivative = %d\nOrder of accuracy = %d\n', q,r);

    C = zeros(sw,sw);
    A1 = zeros(k+1,k+1);
    A2 = zeros(k+1,p);
    dqdn = zeros(sw,p);
    delta_qn = eye(p);

    nn=0;
    for n=0:k
        nn=nn+1;
        j=0; 
        for m=-M1:M2
            j = j+1;
            C(nn,j) = m^n/factorial(n);   % Regular finite difference
            if(n==q)
                b(nn) = 1;
            else
                b(nn) = 0;
            end
        end
    end
    
    A1 = C(1:k+1,1:k+1);
    A2 = -C(1:k+1,k+2:M1+M2+1);

    dqdn(1:k+1,1:p) = A1\(A2*delta_qn);
%     dqdn(1:k+1,1:p) = inv(A1)*A2*delta_qn;
    dqdn(k+2:M1+M2+1,1:p) = delta_qn;
    dqdn = dqdn';

    % spectral optimization
%          nnn = M2;
%    for n=k+1:kmax
    nnn = 0;
   for n=M2-p+1:M2  % derivative wrt free parameter 'a_n', n=M2-p+1:M2
        nn=nn+1; nnn = nnn+1;
        
        j=0; tmp2 = 0;
        for m=-M1:M2
            j = j+1;
            
            % LHS of eqn. for optimum resolution
            l = 0; tmp=0;
            for q=-M1:M2
                l = l+1;
                if(m==q)
                    tmp = tmp + dqdn(nnn,l)*kc;
                else
                    tmp = tmp + dqdn(nnn,l) * sin( kc*(m-q) )/(m-q);
                end
            end
            C(nn,j) = tmp;
            
            % RHS of eqn. for optimum resolution
            if(m~=0)
                tmp2 = tmp2 + dqdn(nnn,j) * (sin(kc*m)-kc*m*cos(kc*m))/(m^2);
            end
        end
        b(nn) = tmp2;
   end

   
%     fprintf('\nWeights: \n');
    w = C\b';
%     w = inv(C)*b';
    return;
    w_abs = (w)
%     fprintf('%6.8e  ', w'); fprintf('\n');
    
    % compute 'keq/k' for DRP test
    im=sqrt(-1);
    kk = 0;
    dk = 0.01;
%     dk = pi/(M1+M2+1);
    for kn=0:dk:pi
        kk = kk+1;  kdx(kk) = kn;
        mm = 0;     ksum = 0;
        for m=-M1:M2
            mm = mm+1;
            ksum = ksum + w(mm)*exp(im*kn*m);
        end       
        keq(kk) = -im*ksum;
        
        % scale 'keq' with 'k'
        if(kn==0)
            keq(kk) = 1;
        else
%             keq(kk) = keq(kk)/kn;
            keq(kk) = real(keq(kk))/kn + im*imag(keq(kk));
        end
    end
        
    k_num = keq.';
%     ind = linspace(-M1,M2,size(w,1));
    ind = -M1:1:M2;
    fs = 16;
    figure(12); clf; hold on;
%     subplot(2,2,1); hold on;
    subplot(3,1,1); hold on;
        plot(kdx, real(keq)); max_kreal = max(real(keq))
    line([0 pi],[1 1],'LineStyle','--');
    xlabel('k \Delta x'); ylabel('[k_{eq}/k]_{real}'); 
    set(gca,'FontSize',fs); box; grid; axis tight; 
%     ylim([0,1.5]);

%     figure(13); clf; hold on;
%     subplot(2,2,2); hold on;
    subplot(3,1,2);hold on;
        plot(kdx, imag(keq)); max_kimag = max(imag(keq))
    line([0 pi],[0 0],'LineStyle','--');
%     xlabel('k \Delta x'); ylabel('[k_{eq}/k]_{imag}');
    xlabel('k \Delta x'); ylabel('[k_{eq}]_{imag}');
    set(gca,'FontSize',fs); box; grid; axis tight; 
%     ylim([-1,.5]);

%     figure(14);clf;
%     subplot(2,2,3);
    subplot(3,1,3);
    plot(ind,abs(w));
    xlabel('m'); ylabel('|a_m|');
    set(gca,'FontSize',fs); box on; grid minor; 
%      ylim([0,1.5]);
   
%     saveas(gcf,'~/Dropbox/Compact_JCP/Figures/compact/stencil_40_20_IVorder.jpg');

end

