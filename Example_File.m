function Example_File
%Eample File
%Run the Sub_Fig_Divider after running this file.
close all;
X1=1:100;Y1=sind(X1);subplot(2,2,1);plot(X1,Y1);legend('L1');title('First');xlabel('x');ylabel('f(x)');
X2=100:200;Y2=cosd(X2);subplot(2,2,2);plot(X2,Y2);legend('L2');title('Second');xlabel('x');ylabel('f(x)');
X3=200:300;Y3=2*cosd(X3);subplot(2,2,3);plot(X3,Y3);legend('L3');title('Third');xlabel('x');ylabel('f(x)');
X4=300:400;Y4=2*sind(X4);subplot(2,2,4);plot(X4,Y4);legend('L4');title('Fourth');xlabel('x');ylabel('f(x)');
end