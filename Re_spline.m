%-----Re���ƔC�ӈ�̋�̓f�[�^�w��A���̒l�̃X�v���C�����o��
%-----Re:�~�������C�m���Y��,foil_num:���^�ԍ�,ad:��ɂȂ��̓f�[�^�ԍ�,ad_val:��̒l,want:�~������̓f�[�^�ԍ�,adata:adata
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