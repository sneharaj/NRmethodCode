% Program loadflow_nr
% THIS IS THE NEWTON-RAPHSON POWER FLOW PROGRAM

clear all

d2r=pi/180;w=100*pi;

% The Ybus matrix is

[yb,ych]=ybus;

g=real(ybus);b=imag(ybus);

% The given parameters and initial conditions are

p=[0;-8;4.4;0;0];
q=[0;-2.8;0;0;0];
mv=[1;1;1.05;1;1];
th=[1e-3;1e-3;1e-3;1e-3;1e-3];
 
del=1;indx=0;

% The Newton-Raphson iterations starts here

%while (del>1e-6)
   for i=1:5
      P(1,i)=0;
      Q(1,i)=0;
      for k=1:5
          if (i==k)
              P(1,i)=P(1,i);
              Q(1,i)=Q(1,i);
          else
              P(1,i)=P(1,i)+abs(mv(i)*mv(k)*yb(i,k))*cos(th(k)-th(i));
              Q(1,i)=Q(1,i)-abs(mv(i)*mv(k)*yb(i,k))*sin(th(k)-th(i));
          end
      end
      P(1,i)=P(1,i)+abs(mv(i))^2*g(i,i);
      Q(1,i)=Q(1,i)-abs(mv(i))^2*b(i,i);
   end

% The mismatches

   Delta_P=p-P';
   Delta_Q=q-Q';

% The Jacobian matrix

%J11
    for i=1:4
       for k=1:4
           if(i==k)
               J_11(i,k)=-1*Q(i)-abs(mv(i)*mv(i))*b(i,i);
           else
               J_11(i,k)=-1*abs(mv(i)*mv(k)*yb(i,k))*sin(th(k)-th(i));
           end
      end
    end
%J21
   for i=1:3
       for k=1:4
           if(i==k)
               J_21(i,k)=P(i)-abs(mv(i)*mv(i))*g(i,i);
           else
               J_21(i,k)=-1*abs(mv(i)*mv(k)*yb(i,k))*cos(th(k)-th(i));
           end
       end
   end

%J12
 for i=1:4
       for k=1:3
           if(i==k)
               J_12(i,k)=2*abs(mv(i)*mv(i))*g(i,i)+J_21(i,i);
           else
               J_12(i,k)=-1*J_21(k,i);
           end
       end
   end   
%J22
 for i=1:3
       for k=1:3
           if(i==k)
               J_22(i,k)=-1*abs(mv(i)*mv(i))*b(i,i)-J_11(i,i);
           else
               J_22(i,k)=J_11(i,k);
           end
       end
   end
   J=[J_11 J_12;J_21 J_22];  

   Delta_PQ=[Delta_P(2:5);Delta_Q(2:4)];  

   Corr=inv(J)*Delta_PQ;  

   th=th+[0;Corr(1:4)];

   mv=mv+[0;mv(2:4).*Corr(5:7);0];  

   del=max(abs(Delta_PQ));

   indx=indx+1;
   index=1;
   for N= 1:0.005:1.08
       mv = [N; 1; 1; 1; 1]
        for i=1:5
      P(1,i)=0;
      Q(1,i)=0;
      for k=1:5
          if (i==k)
              P(1,i)=P(1,i);
              Q(1,i)=Q(1,i);
          else
              P(1,i)=P(1,i)+abs(mv(i)*mv(k)*yb(i,k))*cos(th(k)-th(i));
              Q(1,i)=Q(1,i)-abs(mv(i)*mv(k)*yb(i,k))*sin(th(k)-th(i));
          end
      end
      P(1,i)=P(1,i)+abs(mv(i))^2*g(i,i);
      Q(1,i)=Q(1,i)-abs(mv(i))^2*b(i,i);
        end
       Q_N(index)=Q(1,1)+ Q(1,2)+Q(1,3)+Q(1,4)+Q(1,5);
       mv_N(index)=mv(1);
       index=index+1;
   end
   plot(mv_N,Q_N);
   
       
% preal=(pcal+[0 0 0 0 0.24])*100;
% 
% preac=(qcal+[0 0 0 0 0.11])*100;
% 
% % Power flow calculations
% 
% for i=1:5
%    v(i)=mv(i)*exp(1i*th(i));
% end
% 
% for i=1:4
%    for k=i+1:5
%       if (ybus(i,k)==0)
%          s(i,k)=0;s(k,i)=0;
%          c(i,k)=0;c(k,i)=0;
%          q(i,k)=0;q(k,i)=0;
%          cur(i,k)=0;cur(k,i)=0;
%       else
%          cu=-(v(i)-v(k))*ybus(i,k);
%          s(i,k)=-v(i)*cu'*100;
%          s(k,i)=v(k)*cu'*100;
%          c(i,k)=100*abs(ych(i,k))*abs(v(i))^2;
%          c(k,i)=100*abs(ych(k,i))*abs(v(k))^2;
%          cur(i,k)=cu;cur(k,i)=-cur(i,k);
%       end
%    end
% end
% 
% pwr=real(s);
% qwr=imag(s);
% 
% q=qwr-c;
% 
% % Power loss
% 
% ilin=abs(cur);
% 
% for i=1:4
%    for k=i+1:5
%       if (ybus(i,k)==0)
%          pl(i,k)=0;pl(k,i)=0;
%          ql(i,k)=0;ql(k,i)=0;
%       else
%          z=-1/ybus(i,k);
%          r=real(z);
%          x=imag(z);
%          pl(i,k)=100*r*ilin(i,k)^2;pl(k,i)=pl(i,k);                 ql(i,k)=100*x*ilin(i,k)^2;ql(k,i)=ql(i,k);
%       end
%    end
% end