%任意混合率における任意clを満たすαを探索する関数
function a_def=aoa_def(Re,cl,f_mix,adata)
mix_list={1;0.8;0.6;0.4;0.2;0};
for n=1:6
	n
	a_dex(n,1)=mix_list{n};
	a_dex(n,2)=Re_spline(Re,n,2,cl,1,adata);
	a_dex(n,3)=Re_spline(Re,n,2,cl,2,adata);
	a_dex(n,4)=Re_spline(Re,n,2,cl,3,adata);
	a_dex(n,5)=Re_spline(Re,n,2,cl,5,adata);
	a_dex(n,6)=Re_spline(Re,n,2,cl,10,adata);
end
%clを基準に各空力データを取得する
for n=1:6
	n
	a_def(1,n)=spline(a_dex(:,1),a_dex(:,n),f_mix);
end
%adefはmix,aoa,cl,cd,cm,xcpの順

endfunction
