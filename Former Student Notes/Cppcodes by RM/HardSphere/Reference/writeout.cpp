	FILE *writeout;
	char buf[1000];
	sprintf(buf, "%s%s", "initialconfig",name);
	writeout = fopen(buf, "w+");
	fprintf(writeout, "MODEL%*d\n", 9, 0);
	fprintf(writeout, "CRYST1%*.3f%*.3f%*.3f%*.2f%*.2f%*.2f%*s%*d\n", 9, 100.0, 9, 100.0, 9, 100.0, 7, 90.00, 7, 90.00, 7, 90.00, 11, "P", 4, 1);   
	writePeriodicCoordinates(writeout, beadinfo); 
	fprintf(writeout, "ENDMDL\n");
	
	
	/*	writePeriodicCoordinates
	*		This function writes out the periodic coordiantes to the specified file in the format necessary for vmd pdb files
 */
void writePeriodicCoordinates(FILE *PDBFile, double *rx, double *ry, double *rz)
{
	for(int i = 0; i < n; i++)
	{
		double rxi = rx[i];
		double ryi = ry[i];
		double rzi = rz[i];
		double boxlx, boxly, boxlz;
		boxlx = 1.0;
		boxly = 1.0;
		boxlz = 1.0;

		rxi = rxi - (boxlx*(double)((int)(rxi/(double)(boxlx) + (0.5*copysign(1.0,rxi)))));
		ryi = ryi - (boxly*(double)((int)(ryi/(double)(boxly) + (0.5*copysign(1.0,ryi)))));
		rzi = rzi - (boxlz*(double)((int)(rzi/(double)(boxlz) + (0.5*copysign(1.0,rzi))))); 
		fprintf(PDBFile, "ATOM %*d  N%*d", 6, i, 6, beadinfo[i].beadtype);
		fprintf(PDBFile, "%*.3f%*.3f%*.3f\n", 17, 100.0*rxi, 8, 100.0*ryi, 8, 100.0*rzi);
	}
}