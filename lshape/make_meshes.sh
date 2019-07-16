cd ../
make -j



echo writing meshes 
for value in {0..8}
do
	./lshape/make_lmesh /home/kraemer/Programmieren/txts/uniform_lshape/mesh/ $value
done


echo All meshes written to file
