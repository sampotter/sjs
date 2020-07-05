struct mylist {
	struct mylist *next; // pointer to the next node in the list
	struct mylist *previous; // pointer to the previous node in the list
	int ind; // index of node in the list
	int ibucket; // index of bucket where the point is currently located
};

struct mybucket {
	struct mylist *list; // pointer to the list of nodes in this bucket
	double minval; // minimal possible value in the bucket
	int count; // the number of points in the bucket
};

struct bucket_sort_stuff {
	double gap; // the minimal difference between the value at child and the value at parent 
	struct mylist *list; // list is associated with every mesh point
	int Nbuckets; // the number of buckets
	struct mybucket *bucket;
	int bcount; // the number of boundary points
	int *bdry; // indices of boundary points
	double *blist; // list of values of boundary points
	int jbdry; // the first index of  boundary point with no assigned bucket
	double Bmax; // boundary points with values less than Bmax should be assigned to buckets 
	int ibcurrent; // the index of the current bucket
};


void dial_list_init(struct mylist *list,int ind);
void dial_bucket_init(struct mybucket *bucket,int iskip,double gap);
void print_buckets(int Nbuckets,struct mybucket *bucket);
int adjust_bucket(int ind,double newval,double g,int Nbuckets,struct mybucket *bucket,struct mylist *list);
int find_bucket(double utemp,double g);
void myfree(struct bucket_sort_stuff  *BB);
void start_filling_buckets(struct bucket_sort_stuff  *BB,int Nbuckets,struct mybucket *bucket,
		struct mylist *list,double gap,int *bdry,double *blist,int bcount);
int find_number_of_buckets(double gap,double maxgap);		