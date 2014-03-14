%重量と座標から慣性テンソルを求める
function ine_res=calc_ine(mass,xd,yd,zd)

	ined_xx=mass*(yd^2+zd^2);
	ined_yy=mass*(xd^2+zd^2);
	ined_zz=mass*(xd^2+yd^2);

	ined_xy=-1*mass*xd*yd;
	ined_xz=-1*mass*xd*zd;
	ined_yz=-1*mass*yd*zd;

	ine_res=[ined_xx ined_xy ined_xz;ined_xy ined_yy ined_yz;ined_xz ined_yz ined_zz];

endfunction