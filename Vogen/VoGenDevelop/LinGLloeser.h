#define INDEX(i,j) [j+(i)*coln]

#define DOFOR(i,to) for (i=0;i<to;i++)
#define DFOR(i,from,to) for(i=from-1;i<to;i++)
#define DOBY(i,from,to,by) for(i=from-1;i<to;i+=by)
#define DOBYY(i,from,to,by) for(i=from;i<to;i+=by)
#define DOBYYY(i,from,to) for(i=from;i<to;i++)
#define DOV(i,to,by) for(i=0;i<to;i+=by)

#define INDEX1(i,j) [j-1+(i-1)*coln]
#define INDEXC(i,j) [j-1+(i-1)*rown]

#define VECTOR(i) [i-1]

#define minLinGLloeser(a,b) (((a)<(b))?(a):(b))
#define maxLinGLloeser(a,b) (((a)<(b))?(b):(a))
#define abs(x) (((x) > 0.)?(x):-(x))

void LinGlLoesung (float *a, float *b,int n);
