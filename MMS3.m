%% method of manufactured solutions :: homogeneous solutions
function [p, s] = MMS3(xx,yy)   % used for homogeneous solution p^h
    test_case = 8;
    rr = sqrt(xx^2+yy^2); 
%     th=atan2(yy,xx)*180/(pi)
    th=atan2(yy,xx);%+pi/32;
    
    switch(test_case)
        case(1)
            % test case 1
            rk(1) = sqrt(xx^2+(yy-0.45)^2);
            rk(2) = sqrt(xx^2+(yy+0.45)^2);
            rk(3) = sqrt((xx-0.45)^2+yy^2);
            rk(4) = sqrt((xx+0.45)^2+yy^2);
            temp  = 0;
            for k=1:4
                temp = temp + 40*(10*rk(k)^2-1)*exp(-10*rk(k)^2);
            end
            s = 12*(10*rr^2-1)*exp(-10*rr^2) + temp;
            temp  = 0;
            for k=1:4
                temp = temp + exp(-10*rk(k)^2);
            end
            p = 0.3*exp(-10*rr^2) + temp;
    
        case(2)
            % test case 2
            s = 7*rr^2*cos(3*th);
            p = rr^4*cos(3*th);

        case(3)
            % test case 3
            k = 2*pi*1;
            s = -2*k*sin(k*xx)*cos(k*yy);
            p =      sin(k*xx)*cos(k*yy);
        case(4)
            % test case 4
            Ri = 0.25; 
            Ro = 0.75;
            q = 1;  %BC: dp/dr=q @ r=Ri and p=0 @ r=Ro
            
            s = 1;
            p = (rr^2-Ro^2)/4 + Ri*(q-Ri/2)*log(rr/Ro);
%             s = 0;
%             p = Ri*q*log(rr/Ro);
        case(5)
            % test case 5
            Ri = 0.25; 
            Ro = 0.75;
            k = 2;
            q = 1;  %BC: dp/dr=q @ r=Ri and p=0 @ r=Ro
            
            s = (4-k^2) * sin(k*th);
            p = ( rr^2 - Ro^2 + Ri*(q-2*Ri)*log(rr/Ro) )*sin(k*th); 
        case(6)
            % test case: Homogeneous solution
            k = 8;            
            s = 0;
            p = ( rr^k )*sin(k*th); 
        case(7)
            % test case: Homogenous solution 
            R = 0.5;
            q = .1;
            s = xx*yy;
            p = xx*yy* ( -2*R^2)/12 + R*q*log(rr); 
        case(8)
            % test case: Homogenous solution 
            R = 0.5;
            q = 1;
            k = 1;
            n = -5;
            Ro = sqrt(2);
            D = Ro^(max(k,n));
%             D = Ro^n;

            s = (k^2-n^2)*rr^(k-2)*sin(n*th)/D;
            p = -(k/n)*(rr/R)^(n-k) * (rr)^k*sin(n*th)/D + R*q*log(rr); 
        case(9)
            % test case: Particular + Homogenous solution             
            s = -200*sin(10*(xx+yy)) + 2;
            p =      sin(10*(xx+yy)) + xx^2 - 3*yy + 8;
        case(10)
            % test case 4
            s = 1;
            p = rr^2/4;
    end
    
    return;
end


