double norm2(const double v[3]);
double dot(const double v1[3],const double v2[3]);
void cross(const double v1[3],const double v2[3],double result[3]);
void vsum(const double v1[3],const double v2[3],double result[3]);
void scalar(const double v1[3],const double v2,double result[3]);
void rotateX(const double v1[3],const double alpha, double result[3]);


/***************************************************************************/
double norm2(const double v[3])
{
return (sqrt(dot(v,v)));
}
/***************************************************************************/
double dot(const double v1[3],const double v2[3])
{
int i;
double sum;

sum = 0.0;
for (i = 0; i < 3; i++)
        sum += v1[i]*v2[i];
return sum;
}
/***************************************************************************/
void cross(const double v1[3],const double v2[3],double result[3])
{
result[0] = v1[1]*v2[2]-v1[2]*v2[1];
result[1] = v1[2]*v2[0]-v1[0]*v2[2];
result[2] = v1[0]*v2[1]-v1[1]*v2[0];
}
/***************************************************************************/

void vsum(const double v1[3],const double v2[3],double result[3])
{
result[0] = v1[0]+v2[0];
result[1] = v1[1]+v2[1];
result[2] = v1[2]+v2[2];
}
/***************************************************************************/
void scalar(const double v1[3],const double v2,double result[3])
{
result[0] = v1[0]*v2;
result[1] = v1[1]*v2;
result[2] = v1[2]*v2;
}
/***************************************************************************/
void rotateX(const double v1[3],const double alpha, double result[3])
{
//transforms vector v1 in new co-system that is rotated by alpha
//around x-axis (RHS!)
result[0] = v1[0];
result[1] = v1[1]*cos(alpha)+v1[2]*sin(alpha);
result[2] = -v1[1]*sin(alpha)+v1[2]*cos(alpha);
}
/***************************************************************************/
