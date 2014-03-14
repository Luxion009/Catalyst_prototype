%-----Re数と任意一つの空力データ指定、その値のスプラインを出力
%-----Re:欲しいレイノルズ数,foil_num:翼型番号,ad:基準になる空力データ番号,ad_val:基準の値,want:欲しい空力データ番号,adata:adata
function Re_sp=Re_spline(Re,foil_num,ad,ad_val,want,adata)

Relist={50000;100000;150000;200000;250000;300000;350000;400000;450000;500000;550000;600000};
for n=1:12
	for m=1:10
		Re_s_list(n,1)=Relist{n};
		Re_s_list(n,m+1)=spline(adata{foil_num,n}(:,ad),adata{foil_num,n}(:,m),ad_val);
	endfor
endfor
Re_sp=spline(Re_s_list(:,1),Re_s_list(:,want+1),Re);
endfunction