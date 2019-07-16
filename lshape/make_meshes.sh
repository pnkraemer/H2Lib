cd ../
make -j



echo writing meshes 
for value in {0..8}
do
	./lshape/make_lmesh $1 $value
done


echo All meshes written to file
