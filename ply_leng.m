%決められた桁の決められた０度層の長さを変更する
%積層の層番号と延長終了位置座標を入力
function ply_data=ply_leng(ply_data,cp,ply_end)

ply_data(cp,6)=ply_end;

endfunction