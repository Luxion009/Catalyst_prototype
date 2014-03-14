%-----—ƒŒ^‚ğ—‘zCL‚É‚ ‚í‚¹‚Ä¬‡‚·‚éŠÖ”
function mix=foil_mix(Re,cl,aoa,adata)
mix_list2={1;0.8;0.6;0.4;0.2;0};
for n=1:6
	mix_cl(n,1)=mix_list2{n};
	mix_cl(n,2)=Re_spline(Re,n,1,aoa,1,adata);
	mix_cl(n,3)=Re_spline(Re,n,1,aoa,2,adata);
	mix_cl(n,4)=Re_spline(Re,n,1,aoa,3,adata);
	mix_cl(n,5)=Re_spline(Re,n,1,aoa,5,adata);
	mix_cl(n,6)=Re_spline(Re,n,1,aoa,10,adata);
end
for n=1:6
	mix(1,n)=spline(mix_cl(:,3),mix_cl(:,n),cl);
end

endfunction
