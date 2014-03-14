function ply_data=ply_ex(ply_data,ply_num,ply_ang,C_T,ply_end)


ply_len=ply_data(rows(ply_data),6);
ply_data(rows(ply_data),:)=[];
for n=1:abs(ply_num)
	cp=rows(ply_data)+1;
	ply_data(cp,2)=40;
	ply_data(cp,3)=ply_ang;
	ply_data(cp,4)=0;
	ply_data(cp,5)=0;
	ply_data(cp,6)=ply_end;
end
cpe=rows(ply_data)+1;
ply_data(cpe,2)=24;
ply_data(cpe,3)=90;
ply_data(cpe,4)=90;
ply_data(cpe,5)=0;
ply_data(cpe,6)=ply_len;	

for j=2:rows(ply_data)
	ply_data(j,1)=ply_data(j-1,1)+C_T(lookup(C_T(:,1),ply_data(j,2)),2)*2;
end

endfunction