%�ϑw�p�x�A���a�A�������w�肵�Ēf�ʓ񎟃��[�����g���v�Z
function I=calc_m2(deg,thick,dia)
p=((thick+dia).^4-dia.^4)/32;
%-----Ix
I(:,1)=p.*(deg.*pi./180-(sin(2.*deg.*pi./180)-sin(-2.*deg.*pi./180))/4);
%-----Iz
I(:,2)=p.*(deg.*pi./180+(sin(2.*deg.*pi./180)-sin(-2.*deg.*pi./180))/4);

endfunction
