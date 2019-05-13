

#ifndef HELPFCTS_H
#define HELPFCTS_H
void 
lookatcluster(pcluster t)
{
	printf("\n\nSize of the cluster: %d indices", t->size);
	printf("\nIndex set of cluster: ");
	for (uint i = 0; i < t->size; i++){
		printf("%d ", t->idx[i]);
	}
	printf("\nNumber of Sons: %d", t->sons);
	printf("\nDimensionality of bdg box: %d", t->dim);
	printf("\nMinimal coordinate(s) of bdg box:");
	for (uint i = 0; i < t->dim; i++){
		printf("%f ", t->bmin[i]);
	}
	printf("\nMaximal coordinate(s) of bdg box:");
	for (uint i = 0; i < t->dim; i++){
		printf("%f ", t->bmax[i]);
	}
	printf("\nNumber of descendants %d", t->desc);
	printf("\nType of the cluster: %d", t->type);
	printf("\nHow deep is the cluster?: %d level\n\n", getmindepth_cluster(t));
}



pavector clone_avector(pavector vec) {
	uint dim = vec->dim;
	pavector newvec = new_zero_avector(dim);
	copy_avector(vec, newvec);
	return newvec;
}
#endif
