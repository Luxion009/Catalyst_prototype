%-----混合率と迎え角から空力データ(揚力係数、抗力係数、モーメント、風圧中心)を引っ張ってくる関数
function lstd=lst_exp(Re,mix,aoa,adata)
mix_list={1;0.8;0.6;0.4;0.2;0};
for n=1:6
	m_lstd(n,1)=mix_list{n};
	m_lstd(n,2)=aoa;
	m_lstd(n,3)=Re_spline(Re,n,1,aoa,2,adata);
	m_lstd(n,4)=Re_spline(Re,n,1,aoa,3,adata);
	m_lstd(n,5)=Re_spline(Re,n,1,aoa,5,adata);
	m_lstd(n,6)=Re_spline(Re,n,1,aoa,10,adata);
endfor
	lstd(1,1)=spline(m_lstd(:,1),m_lstd(:,3),mix);
	lstd(1,2)=spline(m_lstd(:,1),m_lstd(:,4),mix);
	lstd(1,3)=spline(m_lstd(:,1),m_lstd(:,5),mix);
	lstd(1,4)=spline(m_lstd(:,1),m_lstd(:,6),mix);
endfunction