%Ï‘w\¬‚ğ•\¦‚·‚é
function ply_data=ply_disp(ply_data,cw)
	printf("-------------%dwing ply data-------------\n",cw);
	printf("/num\t/dia\t/CFRP\t/ang\t/type\t/sp\t/ep\t/wid\t\n");
	for pn=1:rows(ply_data(:,2))
		printf("/%d\t",pn);
		for dn=1:6
			printf("/%d\t",ply_data(pn,dn));
		endfor
		wid=pi*ply_data(pn,3)*ply_data(pn,1)/180;
		printf("/%.2f\t",wid);
		printf("\n");
	endfor
	
endfunction