%積層角度、内径、厚さを指定して断面二次モーメントを計算
function I=calc_m2(deg,thick,dia)
p=((thick+dia).^4-dia.^4)/32;
%-----Ix
I(:,1)=p.*(deg.*pi./180-(sin(2.*deg.*pi./180)-sin(-2.*deg.*pi./180))/4);
%-----Iz
I(:,2)=p.*(deg.*pi./180+(sin(2.*deg.*pi./180)-sin(-2.*deg.*pi./180))/4);

endfunction
