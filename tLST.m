function out=tLST(wing,vel,mangle,nu,aero_data)
%-----状況を入力して揚力係数を返す
%-----wingは[スパン方向位置，翼弦長]の順


%-----各レイノルズ数にて零揚力角と揚力傾斜を求める
%foil=1;
%for i=1:12
	%poly1(i,:)=polyfit(aero_data{foil,i}(1:25,1),aero_data{foil,i}(1:25,2),1);
	%Cl0(i,1)=roots(poly1(i,:));%deg
	%relist(i,1)=50000*i;
%end
	%Clalpha(:,1)=poly1(:,1);%1/deg
	%-----最初の迎え角,レイノルズ数
	aoa(1:rows(wing),1)=mangle;
	re(1:rows(wing),1)=vel.*wing(:,2)./nu;
	
disp("calculating gamma")

for c=1:2
	
	c
	%Cl(:,1)=spline(relist,Clalpha,re).*(aoa-spline(relist,Cl0,re));
	%安心と信頼のRe_spline
	for i=1:rows(wing)
		%printf('%.3d ',i)
		Cl(i,1)=Re_spline(re(i,1),1,1,aoa(i,1),2,aero_data);
	end
		gamma(:,1)=1./2.*vel.*wing(:,2).*Cl;
		
		dgamma_dy(:,1)=gradient(gamma,wing(:,1));
		
		for i=1:rows(wing)
			for n=1:rows(wing)
				if n==i
					dwh(n,1)=dgamma_dy(n,1)/2;
				else
					dwh(n,1)=dgamma_dy(n,1)/abs(wing(n,1)-wing(i,1))/pi/(-4);
				end
				
				%対称翼のための誘導速度見積もり
				%dwh(n+rows(wing),1)=dgamma_dy(n,1)/abs(wing(n,1)+wing(i,1))/pi/(-4);
				
			end
			dw(i,1)=trapz(wing(:,1),dwh);
		end
		
		aoa(:,1)=mangle-atan(dw./vel).*180./pi;
		re(:,1)=sqrt(vel.^2+dw.^2).*wing(:,2)./nu;
				
end		

for i=1:rows(wing)
	%printf('%.3d ',i)
	Cl(i,1)=Re_spline(re(i,1),1,1,aoa(i,1),2,aero_data);
end
			
out=[Cl,re,aoa,dw];
			
endfunction