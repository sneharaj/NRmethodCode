% Function ybus
% THIS IS THE PROGRAM FOR CREATING Ybus MATRIX.
%Just a test comment

function [yb,ych]=ybus

 % The line impedances are

 zz=[0 0 0 0 0.0015+0.02i
    0 0 0 0.009+.1i 0.0045+0.05i
    0 0 0 0.00075+0.01i 0
    0 0.009+.1i 0.00075+0.01i 0 0.002251+0.025i
    0.0015+0.02i 0.0045+0.05i 0 0.002251+0.025i 0];


 % The line chargings are

 ych=[0 0 0 0 0
    0 0 0 0.86i 0.44i
    0 0 0 0 0
    0 0.86i 0 0 0.22i
    0 0.44i 0 0.22i 0];

  % The Ybus matrix is formed here

  for i=1:5
     for j=1:5
        if zz(i,j) == 0
           yb(i,j)=0;
        else
           yb(i,j)=-1/zz(i,j);
        end
     end
  end

  for i=1:5
     ysum=0;
     csum=0;
     for j=1:5
      ysum=ysum+yb(i,j);
      csum=csum+ych(i,j);
   end
   yb(i,i)=csum-ysum;
  end
  
 