void Glob_Star(const double [3],const double,\
			   const double,const double, double [3]);
void Star_Glob(const double [3],const double,\
			   const double,const double, double [3]);

/***************************************************************************/
void Glob_Star(const double x[3],const double nu,\
			   const double eps,const double psi, double xsi[3])
{
//transformation from x-reference frame to xsi_star-reference frame
//first rotation:  x   -> x', rotation about x-axis by nu
//second rotation: x'  -> xsi, rotation about y'-axis by eps
//third rotation:  xsi -> xsi*, rotation about zeta-axis by psi
//all rotations follow right-hand rule
//
double cnu = cos(nu), snu=sin(nu);
double ceps = cos(eps), seps=sin(eps);
double cpsi = cos(psi), spsi=sin(psi);
double rot[3][3];		//rotation matrix

//the rotation matrix:
	rot[0][0] = cpsi*ceps;
				rot[0][1] = cpsi*seps*snu+spsi*cnu;
							rot[0][2] =-cpsi*seps*cnu+spsi*snu;
	rot[1][0] =-spsi*ceps;
				rot[1][1] =-spsi*seps*snu+cpsi*cnu;
 							rot[1][2] = spsi*seps*cnu+cpsi*snu;
	rot[2][0] = 	 seps;
				rot[2][1] = 	-ceps*snu;
			 					rot[2][2] =		 ceps*cnu;
//printf("\n rot0 %2.2lf %2.3lf %2.3lf ",rot[0][0],rot[0][1],rot[0][2]);
//printf("\n rot1 %2.2lf %2.3lf %2.3lf ",rot[1][0],rot[1][1],rot[1][2]);
//printf("\n rot2 %2.2lf %2.3lf %2.3lf ",rot[2][0],rot[2][1],rot[2][2]);

	//transforming x into xsi*
	xsi[0] = x[0]*rot[0][0] + x[1]*rot[0][1] + x[2]*rot[0][2];
	xsi[1] = x[0]*rot[1][0] + x[1]*rot[1][1] + x[2]*rot[1][2];
	xsi[2] = x[0]*rot[2][0] + x[1]*rot[2][1] + x[2]*rot[2][2];

}
/***************************************************************************/

/***************************************************************************/
void Star_Glob(const double xsi[3],const double nu,\
			   const double eps,const double psi, double x[3])
{
//transformation from xsi*-reference frame to x-reference frame
//third rotation:  xsi* -> xsi, rotation about zeta*-axis by -psi
//second rotation: xsi  -> x', rotation about eta-axis by -eps
//first rotation:  x'   -> x, rotation about x'-axis by -nu
//all rotations follow right-hand rule
//
double cnu = cos(nu), snu=sin(nu);
double ceps = cos(eps), seps=sin(eps);
double cpsi = cos(psi), spsi=sin(psi);
double rot[3][3];		//rotation matrix

//the rotation matrix:
	rot[0][0] = cpsi*ceps;
				rot[0][1] = cpsi*seps*snu+spsi*cnu;
							rot[0][2] =-cpsi*seps*cnu+spsi*snu;
	rot[1][0] =-spsi*ceps;
				rot[1][1] =-spsi*seps*snu+cpsi*cnu;
 							rot[1][2] = spsi*seps*cnu+cpsi*snu;
	rot[2][0] = 	 seps;
				rot[2][1] = 	-ceps*snu;
			 					rot[2][2] =		 ceps*cnu;

	//transforming x into xsi*
	x[0] = xsi[0]*rot[0][0] + xsi[1]*rot[1][0] + xsi[2]*rot[2][0];
	x[1] = xsi[0]*rot[0][1] + xsi[1]*rot[1][1] + xsi[2]*rot[2][1];
	x[2] = xsi[0]*rot[0][2] + xsi[1]*rot[1][2] + xsi[2]*rot[2][2];

}
/***************************************************************************/
