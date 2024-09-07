#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

#define X1 0.000
#define Y1 0.000
#define Z1 0.000
#define X2 1.256
#define Y2 -0.715
#define Z2 -0.486
#define DX 2.587
#define DY 2.380
#define DZ 3.500

typedef struct position
{
	int molType;
	int atomType;
	double x, y, z;
} BEAD_POSITIONS;

typedef struct simulationBounds
{
	double xlo, xhi, ylo, yhi, zlo, zhi;
} SIMULATION_BOUNDS;

typedef struct coordinates
{
	double x, y, z;
	int chainID1, chainID2, monomerID, atomType, molType;
	bool endGroup, penultimateGroup;
} COORDINATES;

typedef struct bonds
{
	int sino, bondType, batom1, batom2;
} BONDS;

void countNChainsNAtomsFromGro (const char *filename, int *nAtomsPerChain, int *nChains)
{
	FILE *inputGRO;
	inputGRO = fopen (filename, "r");

	char lineString[5000];

	(*nChains) = 0;

	while (fgets (lineString, 5000, inputGRO) != NULL)
	{
		if (strstr (lineString, "50LPE"))
		{
			(*nChains)++;
		}
	}

	(*nChains) /= 2;

	rewind (inputGRO);

	(*nAtomsPerChain) = 0;
	int atomID_previous = 0, atomID_current = 0, nAtoms_current = 0;

	while (fgets (lineString, 5000, inputGRO) != NULL)
	{
		sscanf (lineString, "%d\n", &atomID_current);
		nAtoms_current++;

		if (atomID_current < atomID_previous)
		{
			if (nAtoms_current > (*nAtomsPerChain))
			{
				(*nAtomsPerChain) = nAtoms_current;
			}

			nAtoms_current = 0;
		}

		atomID_previous = atomID_current;
	}

	fclose (inputGRO);
}

BEAD_POSITIONS rotateInXY (BEAD_POSITIONS nextPosition, double theta_xy)
{
	double tempX, tempY, tempZ;
	tempX = nextPosition.x;
	tempY = nextPosition.y;
	tempZ = nextPosition.z;

	nextPosition.x = tempX * cosf (theta_xy) - tempY * sinf (theta_xy);
	nextPosition.y = tempX * sinf (theta_xy) + tempY * cosf (theta_xy);
	nextPosition.z = tempZ;

	return nextPosition;
}

BEAD_POSITIONS rotateInYZ (BEAD_POSITIONS nextPosition, double theta_yz)
{
	double tempX, tempY, tempZ;
	tempX = nextPosition.x;
	tempY = nextPosition.y;
	tempZ = nextPosition.z;

	nextPosition.x = tempX;
	nextPosition.y = tempY * cosf (theta_yz) - tempZ * sinf (theta_yz);
	nextPosition.z = tempY * sinf (theta_yz) + tempZ * cosf (theta_yz);

	return nextPosition;
}

BEAD_POSITIONS rotateInXZ (BEAD_POSITIONS nextPosition, double theta_xz)
{
	double tempX, tempY, tempZ;
	tempX = nextPosition.x;
	tempY = nextPosition.y;
	tempZ = nextPosition.z;

	nextPosition.x = tempX * cosf (theta_xz) + tempZ * sinf (theta_xz);
	nextPosition.y = tempY;
	nextPosition.z = - tempX * sinf (theta_xz) + tempZ * cosf (theta_xz);

	return nextPosition;
}

BEAD_POSITIONS nextConstrainedWalkPosition (BEAD_POSITIONS previousPosition1, BEAD_POSITIONS previousPosition2, BEAD_POSITIONS nextPosition, double bondLength, double bondAngle, int j)
{
	// srand (time(NULL) + j * time(NULL));

	BEAD_POSITIONS angleOfPreviousBond;

	bondAngle = 180 - bondAngle;
	bondAngle /= 57.2958;

	/*
	beads[currentChain][i].x = beads[currentChain][i - 1].x + cos (randomAngle2) * (double)bondLength * cos (randomAngle1);
	beads[currentChain][i].y = beads[currentChain][i - 1].y + cos (randomAngle2) * (double)bondLength * sin (randomAngle1);
	beads[currentChain][i].z = beads[currentChain][i - 1].z + (double)bondLength * sin (randomAngle2);
	*/

	double randomAngle;
	randomAngle = rand() % 360;

	nextPosition.x = cos (bondAngle) * bondLength;
	nextPosition.y = sin (bondAngle) * bondLength;// * cos (randomAngle);
	// nextPosition.z = sin (randomAngle) * bondLength;

	double tempRandY = nextPosition.y * 2;
	double randomY = ((double)rand () / (double)RAND_MAX);
	randomY *= tempRandY;
	randomY -= nextPosition.y;
	// printf("%f %f\n", nextPosition.y, randomY);
	nextPosition.z = sqrt (nextPosition.y * nextPosition.y - randomY * randomY);
	nextPosition.y = randomY;

	// Find the angle of the projection of the previous bond on XY, YZ, and XZ planes
	double dy, dx, dz;
	// XY plane
	double slope_xy, theta_xy;

	dy = previousPosition2.y - previousPosition1.y;
	dx = previousPosition2.x - previousPosition1.x;

	if (dx != 0)
	{
		slope_xy = dy / dx;
		theta_xy = atanf (slope_xy);
	}
	else
	{
		theta_xy = 90.0 / 57.2958;
	}

	nextPosition = rotateInXY (nextPosition, theta_xy);

	// YZ plane
	double slope_yz, theta_yz;

	dz = previousPosition2.z - previousPosition1.z;
	dy = previousPosition2.y - previousPosition1.y;

	if (dy != 0)
	{
		slope_yz = dz / dy;
		theta_yz = atanf (slope_yz);
	}
	else
	{
		theta_yz = 90.0 / 57.2958;
	}

	nextPosition = rotateInYZ (nextPosition, theta_yz);

	// XZ plane
	double slope_xz, theta_xz;

	dz = previousPosition2.z - previousPosition1.z;
	dx = previousPosition2.x - previousPosition1.x;

	if (dx != 0)
	{
		slope_xz = dz / dx;
		theta_xz = atanf (slope_xz);
	}
	else
	{
		theta_xz = 90.0 / 57.2958;
	}

	nextPosition = rotateInXZ (nextPosition, theta_xz);

	// nextPosition = rotateNewBond (nextPosition, previousPosition1, previousPosition2, bondLength);

	nextPosition.x += previousPosition2.x;
	nextPosition.y += previousPosition2.y;
	nextPosition.z += previousPosition2.z;

	// printf("%lf %lf %lf\n", nextPosition.x, nextPosition.y, nextPosition.z);
	// usleep (100000);

	return nextPosition;
}

COORDINATES **resetAtomTypes (COORDINATES **polymerBeads, int i, int nAtomsPerChain)
{
	for (int j = 0; j < nAtomsPerChain; ++j)
	{
		polymerBeads[i][j].atomType = 0;
	}

	return polymerBeads;
}

COORDINATES **setAtomTypes (COORDINATES **polymerBeads, int nChains, int nAtomsPerChain)
{
	polymerBeads[nChains][0].atomType = 1;
	polymerBeads[nChains][1].atomType = 2;

	for (int j = 2; j < nAtomsPerChain; ++j)
	{
		if (polymerBeads[nChains][j - 1].atomType == 1 && polymerBeads[nChains][j].atomType == 0)
		{
			polymerBeads[nChains][j].atomType = 2;
		}
		if (polymerBeads[nChains][j - 1].atomType == 2 && polymerBeads[nChains][j].atomType == 0)
		{
			polymerBeads[nChains][j].atomType = 1;
		}
		if (polymerBeads[nChains][j - 1].atomType == 3 && polymerBeads[nChains][j].atomType == 0)
		{
			polymerBeads[nChains][j].atomType = 1;
		}
	}

	return polymerBeads;
}

COORDINATES **generateAtomTypes (COORDINATES **polymerBeads, int nChains, int nAtomsPerChain, int nMonomers, float nSideGroupsPerChain)
{
	// nPossibleSideGroupLocations = excludes the first and last monomers of every chain
	// int nSideGroupsLeft = ceil (nSideGroupsPerChain * (float)nChains), nPossibleSideGroupLocations = (nAtomsPerChain * nChains) - (4 * nChains);

	float nSideGroupsLeft = 0;
	int nPossibleSideGroupLocations = ((nMonomers - 2) * 2), sideGroupLocation;

	for (int i = 0; i < nChains; ++i)
	{
		nSideGroupsLeft += nSideGroupsPerChain;

		polymerBeads = resetAtomTypes (polymerBeads, i, nAtomsPerChain);

		if (nSideGroupsLeft > 0)
		{
			nSideGroupsLeft -= 1;

			// generating even random numbers by setting the span as half and multiplying by 2.
			// Multiplication by 2 always gives an even number.
			srand (time(NULL) + i);
			sideGroupLocation = floor (rand() % (int)ceil (((float)nPossibleSideGroupLocations / 2)));
			sideGroupLocation = (sideGroupLocation * 2) + 2;

			polymerBeads[i][sideGroupLocation].atomType = 3;

			polymerBeads = setAtomTypes (polymerBeads, i, nAtomsPerChain);
		}
		else
		{
			polymerBeads = setAtomTypes (polymerBeads, i, nAtomsPerChain - 1);
		}
	}

	return polymerBeads;
}

COORDINATES **generateAtomTypesFromGro (const char *filename, int nSkipLines, COORDINATES **polymerBeads, int nChains, int nAtomsPerChain)
{
	FILE *inputGRO;
	inputGRO = fopen (filename, "r");
	char lineString[6000];

	int atomID_previous = 0, atomID_current = 0, nAtoms_current = 0, currentChainID = 0, currentAtomID = 0, currentAtomType;

	polymerBeads = resetAtomTypes (polymerBeads, currentChainID, nAtomsPerChain);

	for (int i = 0; i < nSkipLines; ++i)
	{
		fgets (lineString, 5000, inputGRO);
	}

	while (fgets (lineString, 5000, inputGRO) != NULL)
	{
		sscanf (lineString, "%d\n", &atomID_current);

		if (atomID_current == 1 && atomID_previous == 50)
		{
			currentChainID++;
			currentAtomID = 0;

			polymerBeads = resetAtomTypes (polymerBeads, currentChainID, nAtomsPerChain);
		}

		currentAtomType = 0;

		if (strstr (lineString, "C1"))
		{
			currentAtomType = 1;
		}
		else if (strstr (lineString, "C2"))
		{
			currentAtomType = 2;
		}
		else if (strstr (lineString, "C3"))
		{
			currentAtomType = 3;
		}

		if (currentAtomType == 0)
		{
			goto breakThisLoop;
		}

		// printf("currentChainID: %d; currentAtomID: %d\n", currentChainID, currentAtomID);
		polymerBeads[currentChainID][currentAtomID].atomType = currentAtomType;

/*		if (currentAtomID > nAtomsPerChain || currentChainID > nChains || 1)
		{
			printf("%d %d [%d] => %s\n", currentChainID, currentAtomID, polymerBeads[currentChainID][currentAtomID].atomType, lineString);
			fflush (stdout);
			usleep (100000);
		}
*/
		currentAtomID++;
		atomID_previous = atomID_current;
	}

	breakThisLoop:;

	return polymerBeads;
}

int main(int argc, char const *argv[])
{
	if (argc != 20)
	{
		printf("INCORRECT ARGUMENTS PASSED:\n~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n{~} argv[0]  = ./program\n{~} argv[1]  = (double) crystal xlo\n{~} argv[2]  = (double) crystal xhi\n{~} argv[3]  = (double) crystal ylo\n{~} argv[4]  = (double) crystal yhi\n{~} argv[5]  = (double) crystal zlo\n{~} argv[6]  = (double) crystal zhi\n{~} argv[7]  = (double) overall xlo\n{~} argv[8]  = (double) overall xhi\n{~} argv[9]  = (double) overall ylo\n{~} argv[10] = (double) overall yhi\n{~} argv[11] = (double) overall zlo\n{~} argv[12] = (double) overall zhi\n{~} argv[13] = (int) Number of tight loops\n{~} argv[14] = (int) Number of loose loops\n{~} argv[15] = (int) Number of dangles\n{~} argv[16] = (double) Lamellae tilt angle\n{~} argv[17] = (double) Bond length\n{~} argv[18] = (int) Number of monomers\n{~} argv[19] = (int) Number of side groups per chain.\n\n");
		exit (1);
	}

	FILE *output, *outputData;
	output = fopen ("chainsAndCrystal.xyz", "w");
	outputData = fopen ("chainsAndCrystal.data", "w");

	SIMULATION_BOUNDS crystalBounds, overallBounds, tempBounds;
	crystalBounds.xlo = atof (argv[1]);
	crystalBounds.xhi = atof (argv[2]);
	crystalBounds.ylo = atof (argv[3]);
	crystalBounds.yhi = atof (argv[4]);
	crystalBounds.zlo = atof (argv[5]);
	crystalBounds.zhi = atof (argv[6]);

	overallBounds.xlo = atof (argv[7]);
	overallBounds.xhi = atof (argv[8]);
	overallBounds.ylo = atof (argv[9]);
	overallBounds.yhi = atof (argv[10]);
	overallBounds.zlo = atof (argv[11]);
	overallBounds.zhi = atof (argv[12]);

	int nTightLoops, nLooseLoops, nDangles;
	nTightLoops = atoi (argv[13]);
	nLooseLoops = atoi (argv[14]);
	nDangles = atoi (argv[15]);

	tempBounds.xlo = overallBounds.xlo * 2.0;
	tempBounds.xhi = overallBounds.xhi * 2.0;
	tempBounds.ylo = overallBounds.ylo * 2.0;
	tempBounds.yhi = overallBounds.yhi * 2.0;
	tempBounds.zlo = overallBounds.zlo * 2.0;
	tempBounds.zhi = overallBounds.zhi * 2.0;

	double tiltAngle = atof (argv[16]) / 57.2958, bondLength = atof (argv[17]);

	printf("\nTilt angle of the crystalline phase: %lf\n", atof (argv[16]));

	int nMonomers = atoi (argv[18]), nSideGroupsPerChain = atoi (argv[19]);

	int nAtomsPerChain = (nMonomers * 2) + ceil (nSideGroupsPerChain);

	printf("Number of monomers: %d\nNumber of atoms per chain: %d\n", nMonomers, nAtomsPerChain);

	// Create a crystal structure larger than the boundary
	//		Calculate the number of repetitions
	//		Replicate the structure in 3D
	double temp_xLength = (tempBounds.xhi - tempBounds.xlo), temp_yLength = (tempBounds.yhi - tempBounds.ylo), temp_zLength = (tempBounds.zhi - tempBounds.zlo);
	int nReps_x = ceil (temp_xLength / (double)DX) + 1, nReps_y = ceil (temp_yLength / (double)DY) + 1, nReps_z = ceil (temp_zLength / (double)DZ) + 1;
	int nAtoms = 2 * nReps_x * nReps_y * nReps_z;

	COORDINATES *crystalAtoms, crystalCOM, *trimmedCrystalAtoms, **addedPolymer;
	crystalAtoms = (COORDINATES *) malloc (nAtoms * sizeof (COORDINATES));

	crystalCOM.x = 0;
	crystalCOM.y = 0;
	crystalCOM.z = 0;

	int currentAtom = 0;

	// the crystal atoms are created in a lattice
	// and they are oriented based on the input tilt angle
	// center of mass is also calculated within these 'for' loops
	for (int k = 1; k < nReps_z; ++k)
	{
		for (int j = 1; j < nReps_y; ++j)
		{
			for (int i = 1; i < nReps_x; ++i)
			{
				printf("Creating crystal atoms: %d/%d                   \r", currentAtom, nAtoms);
				fflush (stdout);

				crystalAtoms[currentAtom].x = (double)X1 + (double)i * (double)DX;
				crystalAtoms[currentAtom].y = (double)Y1 + (double)j * (double)DY;
				crystalAtoms[currentAtom].z = (double)Z1 + (double)k * (double)DZ;

				crystalAtoms[currentAtom].x = crystalAtoms[currentAtom].x * cos (tiltAngle) - crystalAtoms[currentAtom].y * sin (tiltAngle);
				crystalAtoms[currentAtom].y = crystalAtoms[currentAtom].x * sin (tiltAngle) + crystalAtoms[currentAtom].y * cos (tiltAngle);

				crystalCOM.x += crystalAtoms[currentAtom].x;
				crystalCOM.y += crystalAtoms[currentAtom].y;
				crystalCOM.z += crystalAtoms[currentAtom].z;

				crystalAtoms[currentAtom].monomerID = i;
				crystalAtoms[currentAtom].chainID1 = k;
				crystalAtoms[currentAtom].chainID2 = j;

				currentAtom++;

				printf("Creating crystal atoms: %d/%d                   \r", currentAtom, nAtoms);
				fflush (stdout);

				crystalAtoms[currentAtom].x = (double)X2 + (double)i * (double)DX;
				crystalAtoms[currentAtom].y = (double)Y2 + (double)j * (double)DY;
				crystalAtoms[currentAtom].z = (double)Z2 + (double)k * (double)DZ;

				crystalAtoms[currentAtom].x = crystalAtoms[currentAtom].x * cos (tiltAngle) - crystalAtoms[currentAtom].y * sin (tiltAngle);
				crystalAtoms[currentAtom].y = crystalAtoms[currentAtom].x * sin (tiltAngle) + crystalAtoms[currentAtom].y * cos (tiltAngle);

				crystalCOM.x += crystalAtoms[currentAtom].x;
				crystalCOM.y += crystalAtoms[currentAtom].y;
				crystalCOM.z += crystalAtoms[currentAtom].z;

				crystalAtoms[currentAtom].monomerID = i;
				crystalAtoms[currentAtom].chainID1 = k;
				crystalAtoms[currentAtom].chainID2 = j;

				currentAtom++;

				// printf("%lf %lf %lf; ", (double)X1 + (double)i * (double)DX, (double)Y1 + (double)j * (double)DY, (double)Z1 + (double)k * (double)DZ);
				// printf("%lf %lf %lf\n", (double)X2 + (double)i * (double)DX, (double)Y2 + (double)j * (double)DY, (double)Z2 + (double)k * (double)DZ);
				// usleep (100000);
			}
		}
	}

	crystalCOM.x /= currentAtom;
	crystalCOM.y /= currentAtom;
	crystalCOM.z /= currentAtom;

	currentAtom = 0;

	printf("\n");

	// recenter the crystal
	for (int k = 1; k < nReps_z; ++k)
	{
		for (int j = 1; j < nReps_y; ++j)
		{
			for (int i = 1; i < nReps_x; ++i)
			{
				printf("Recentering crystal atoms: %d/%d                   \r", currentAtom, nAtoms);
				fflush (stdout);

				crystalAtoms[currentAtom].x -= crystalCOM.x;
				crystalAtoms[currentAtom].y -= crystalCOM.y;
				crystalAtoms[currentAtom].z -= crystalCOM.z;

				currentAtom++;

				printf("Recentering crystal atoms: %d/%d                   \r", currentAtom, nAtoms);
				fflush (stdout);

				crystalAtoms[currentAtom].x -= crystalCOM.x;
				crystalAtoms[currentAtom].y -= crystalCOM.y;
				crystalAtoms[currentAtom].z -= crystalCOM.z;

				currentAtom++;
			}
		}
	}

	currentAtom = 0;
	int nAtomsWithinBounds = 0;

	// Find the number of crystal atoms within the selected region
	for (int k = 1; k < nReps_z; ++k)
	{
		for (int j = 1; j < nReps_y; ++j)
		{
			for (int i = 1; i < nReps_x; ++i)
			{
				if (crystalAtoms[currentAtom].x > crystalBounds.xlo && 
					crystalAtoms[currentAtom].x < crystalBounds.xhi && 
					crystalAtoms[currentAtom].y > crystalBounds.ylo && 
					crystalAtoms[currentAtom].y < crystalBounds.yhi && 
					crystalAtoms[currentAtom].z > crystalBounds.zlo && 
					crystalAtoms[currentAtom].z < crystalBounds.zhi)
				{
					nAtomsWithinBounds++;
				}

				currentAtom++;

				if (crystalAtoms[currentAtom].x > crystalBounds.xlo && 
					crystalAtoms[currentAtom].x < crystalBounds.xhi && 
					crystalAtoms[currentAtom].y > crystalBounds.ylo && 
					crystalAtoms[currentAtom].y < crystalBounds.yhi && 
					crystalAtoms[currentAtom].z > crystalBounds.zlo && 
					crystalAtoms[currentAtom].z < crystalBounds.zhi)
				{
					nAtomsWithinBounds++;
				}

				currentAtom++;
			}
		}
	}

	int nAtomsWithinBounds_prev = nAtomsWithinBounds;
	nAtomsWithinBounds = 0;
	currentAtom = 0;
	trimmedCrystalAtoms = (COORDINATES *) malloc (nAtomsWithinBounds_prev * sizeof (COORDINATES));

	printf("\n");

	// Store the crystal atoms within the selected region
	for (int k = 1; k < nReps_z; ++k)
	{
		for (int j = 1; j < nReps_y; ++j)
		{
			for (int i = 1; i < nReps_x; ++i)
			{
				if (crystalAtoms[currentAtom].x > crystalBounds.xlo && 
					crystalAtoms[currentAtom].x < crystalBounds.xhi && 
					crystalAtoms[currentAtom].y > crystalBounds.ylo && 
					crystalAtoms[currentAtom].y < crystalBounds.yhi && 
					crystalAtoms[currentAtom].z > crystalBounds.zlo && 
					crystalAtoms[currentAtom].z < crystalBounds.zhi)
				{
					printf("Trimming the crystal: %d/%d (%d/%d)                   \r", currentAtom, nAtoms, nAtomsWithinBounds, nAtomsWithinBounds_prev);
					fflush (stdout);

					trimmedCrystalAtoms[nAtomsWithinBounds].x = crystalAtoms[currentAtom].x;
					trimmedCrystalAtoms[nAtomsWithinBounds].y = crystalAtoms[currentAtom].y;
					trimmedCrystalAtoms[nAtomsWithinBounds].z = crystalAtoms[currentAtom].z;
					trimmedCrystalAtoms[nAtomsWithinBounds].monomerID = crystalAtoms[currentAtom].monomerID;
					trimmedCrystalAtoms[nAtomsWithinBounds].chainID1 = crystalAtoms[currentAtom].chainID1;
					trimmedCrystalAtoms[nAtomsWithinBounds].chainID2 = crystalAtoms[currentAtom].chainID2;
					trimmedCrystalAtoms[nAtomsWithinBounds].endGroup = 0;
					trimmedCrystalAtoms[nAtomsWithinBounds].penultimateGroup = 0;
					nAtomsWithinBounds++;
				}

				currentAtom++;

				if (crystalAtoms[currentAtom].x > crystalBounds.xlo && 
					crystalAtoms[currentAtom].x < crystalBounds.xhi && 
					crystalAtoms[currentAtom].y > crystalBounds.ylo && 
					crystalAtoms[currentAtom].y < crystalBounds.yhi && 
					crystalAtoms[currentAtom].z > crystalBounds.zlo && 
					crystalAtoms[currentAtom].z < crystalBounds.zhi)
				{
					printf("Trimming the crystal: %d/%d (%d/%d)                   \r", currentAtom, nAtoms, nAtomsWithinBounds, nAtomsWithinBounds_prev);
					fflush (stdout);

					trimmedCrystalAtoms[nAtomsWithinBounds].x = crystalAtoms[currentAtom].x;
					trimmedCrystalAtoms[nAtomsWithinBounds].y = crystalAtoms[currentAtom].y;
					trimmedCrystalAtoms[nAtomsWithinBounds].z = crystalAtoms[currentAtom].z;
					trimmedCrystalAtoms[nAtomsWithinBounds].monomerID = crystalAtoms[currentAtom].monomerID;
					trimmedCrystalAtoms[nAtomsWithinBounds].chainID1 = crystalAtoms[currentAtom].chainID1;
					trimmedCrystalAtoms[nAtomsWithinBounds].chainID2 = crystalAtoms[currentAtom].chainID2;
					trimmedCrystalAtoms[nAtomsWithinBounds].endGroup = 0;
					trimmedCrystalAtoms[nAtomsWithinBounds].penultimateGroup = 0;
					nAtomsWithinBounds++;
				}

				currentAtom++;
			}
		}
	}

	printf("\n");

	currentAtom = 0;
	nAtomsWithinBounds = 0;
	bool chainEnd;
	chainEnd = 1;

	int nChainEnds = 0, nPenultimateAtoms = 0;

	// Assigning the end groups
	for (int i = 0; i < nAtomsWithinBounds_prev; ++i)
	{
		if (i == 0)
		{
			trimmedCrystalAtoms[i].endGroup = 1;
			trimmedCrystalAtoms[i + 1].penultimateGroup = 1;

			nChainEnds++;
			nPenultimateAtoms++;
		}
		else
		{
			if (trimmedCrystalAtoms[i].chainID1 != trimmedCrystalAtoms[i - 1].chainID1 ||
				trimmedCrystalAtoms[i].chainID2 != trimmedCrystalAtoms[i - 1].chainID2)
			{
				trimmedCrystalAtoms[i - 1].endGroup = 1;
				trimmedCrystalAtoms[i].endGroup = 1;

				nChainEnds++;
				nChainEnds++;

				if ((i - 1) > 0)
				{
					trimmedCrystalAtoms[i - 2].penultimateGroup = 1;
					nPenultimateAtoms++;
				}

				if ((i + 1) < nAtomsWithinBounds_prev)
				{
					trimmedCrystalAtoms[i + 1].penultimateGroup = 1;
					nPenultimateAtoms++;
				}
			}
		}
	}

	// Checking the trimmed crystal
/*	for (int i = 0; i < nAtomsWithinBounds_prev; ++i)
	{
		printf("%lf %lf %lf %d %d %d %d %d\n", trimmedCrystalAtoms[i].x, trimmedCrystalAtoms[i].y, trimmedCrystalAtoms[i].z, trimmedCrystalAtoms[i].monomerID, trimmedCrystalAtoms[i].chainID1, trimmedCrystalAtoms[i].chainID2, trimmedCrystalAtoms[i].endGroup, trimmedCrystalAtoms[i].penultimateGroup);
		usleep (100000);
	}
*/

	// Count the number of bonds in the trimmed crystal
	int nBondsInTrimmedCrystal;
	nBondsInTrimmedCrystal = nAtomsWithinBounds_prev - ceil((double)nChainEnds / 2);

	BONDS *crystalBonds, **polymerBonds;
	crystalBonds = (BONDS *) malloc (nBondsInTrimmedCrystal * sizeof (BONDS));

	int bondIndex = 0;

	// int sino, bondType, batom1, batom2;
	for (int i = 0; i < (nAtomsWithinBounds_prev - 1); ++i)
	{
		if (trimmedCrystalAtoms[i].chainID1 == trimmedCrystalAtoms[i + 1].chainID1 && trimmedCrystalAtoms[i].chainID2 == trimmedCrystalAtoms[i + 1].chainID2)
		{
			crystalBonds[bondIndex].sino = bondIndex + 1;
			crystalBonds[bondIndex].bondType = 1;
			crystalBonds[bondIndex].batom1 = i + 1;
			crystalBonds[bondIndex].batom2 = i + 2;

			bondIndex++;
		}
	}

	// Checking the added bonds
/*	for (int i = 0; i < nBondsInTrimmedCrystal; ++i)
	{
		printf("%d %d %d %d\n", crystalBonds[i].sino, crystalBonds[i].bondType, crystalBonds[i].batom1, crystalBonds[i].batom2);
		usleep (100000);
	}
*/

	printf("\nNumber of chain ends in the trimmed crystal: %d\nNumber of penultimate groups in the trimmed crystal: %d\n\n", nChainEnds, nPenultimateAtoms);

	BEAD_POSITIONS previousPosition1, previousPosition2, nextPosition;

	int nChains;
	// countNChainsNAtomsFromGro ("../after_npt1.gro", &nAtomsPerChain, &nChains);

	addedPolymer = (COORDINATES **) malloc (nChainEnds * sizeof (COORDINATES *));

	for (int i = 0; i < nChainEnds; ++i)
	{
		addedPolymer[i] = (COORDINATES *) malloc (nAtomsPerChain * sizeof (COORDINATES));
	}

	polymerBonds = (BONDS **) malloc (nChainEnds * sizeof (BONDS *));

	for (int i = 0; i < nChainEnds; ++i)
	{
		polymerBonds[i] = (BONDS *) malloc ((nAtomsPerChain - 1) * sizeof (BONDS));
	}

	int currentChain = -1;
	printf("Allocating memory for %d bonds in crystal and %d bonds in mobile polymer chains...\n", nBondsInTrimmedCrystal, nChainEnds * (nAtomsPerChain - 1));
	printf("\n");

	int nAtomsPerChain_gro, nChains_gro;
	// countNChainsNAtomsFromGro ("../after_npt1.gro", &nAtomsPerChain, &nChains);
	printf("nAtomsPerChain: %d\nnChains: %d\n", nAtomsPerChain, nChains);

	addedPolymer = generateAtomTypes (addedPolymer, nChainEnds, nAtomsPerChain, nMonomers, nSideGroupsPerChain);
	// addedPolymer = generateAtomTypesFromGro ("../after_npt1.gro", 7202, addedPolymer, nChains_gro, nAtomsPerChain_gro);

	// Generating a constrained walk chains from all the end groups
	for (int i = 0; i < nAtomsWithinBounds_prev; ++i)
	{
		if (trimmedCrystalAtoms[i].endGroup == 1)
		{
			previousPosition2.x = trimmedCrystalAtoms[i].x;
			previousPosition2.y = trimmedCrystalAtoms[i].y;
			previousPosition2.z = trimmedCrystalAtoms[i].z;

			currentChain++;

			printf("Adding free chains...%d                                             \r", currentChain);
			fflush (stdout);
		}

		if (trimmedCrystalAtoms[i - 1].penultimateGroup == 1 && 
			trimmedCrystalAtoms[i - 1].chainID1 == trimmedCrystalAtoms[i].chainID1 && 
			trimmedCrystalAtoms[i - 1].chainID2 == trimmedCrystalAtoms[i].chainID2)
		{
			previousPosition1.x = trimmedCrystalAtoms[i - 1].x;
			previousPosition1.y = trimmedCrystalAtoms[i - 1].y;
			previousPosition1.z = trimmedCrystalAtoms[i - 1].z;
		}
		else
		{
			previousPosition1.x = previousPosition2.x + 1.53;
			previousPosition1.y = previousPosition2.y;
			previousPosition1.z = previousPosition2.z;
		}

		if (trimmedCrystalAtoms[i + 1].penultimateGroup == 1 && 
			trimmedCrystalAtoms[i + 1].chainID1 == trimmedCrystalAtoms[i].chainID1 && 
			trimmedCrystalAtoms[i + 1].chainID2 == trimmedCrystalAtoms[i].chainID2)
		{
			previousPosition1.x = trimmedCrystalAtoms[i + 1].x;
			previousPosition1.y = trimmedCrystalAtoms[i + 1].y;
			previousPosition1.z = trimmedCrystalAtoms[i + 1].z;
		}
		else
		{
			previousPosition1.x = previousPosition2.x + 1.53;
			previousPosition1.y = previousPosition2.y;
			previousPosition1.z = previousPosition2.z;
		}

		if (previousPosition1.x != 0 && 
			previousPosition1.y != 0 && 
			previousPosition1.z != 0 && 
			previousPosition2.x != 0 && 
			previousPosition2.y != 0 && 
			previousPosition2.z != 0)
		{
			for (int j = 0; j < nAtomsPerChain; ++j)
			{
				if (addedPolymer[currentChain][j].atomType == 3)
				{
					nextPosition = nextConstrainedWalkPosition (previousPosition1, previousPosition2, nextPosition, bondLength, 90, j + i);
				}
				else
				{
					nextPosition = nextConstrainedWalkPosition (previousPosition1, previousPosition2, nextPosition, bondLength, 120, j + i);

					previousPosition1 = previousPosition2;
					previousPosition2 = nextPosition;
				}

				addedPolymer[currentChain][j].x = nextPosition.x;
				addedPolymer[currentChain][j].y = nextPosition.y;
				addedPolymer[currentChain][j].z = nextPosition.z;
				addedPolymer[currentChain][j].chainID1 = trimmedCrystalAtoms[i].chainID1;
				addedPolymer[currentChain][j].chainID2 = trimmedCrystalAtoms[i].chainID2;

				if (j > 0)
				{
					if (addedPolymer[currentChain][j].atomType == 3)
					{
						addedPolymer[currentChain][j - 1].atomType = 1;
					}
					else
					{
						addedPolymer[currentChain][j - 1].atomType = 2;
					}
				}

				// printf("%lf %lf %lf\n", nextPosition.x, nextPosition.y, nextPosition.z);
				// printf("%d/%d: %lf %lf %lf\n", currentChain, j, addedPolymer[currentChain][j].x, addedPolymer[currentChain][j].y, addedPolymer[currentChain][j].z);
				// usleep (100000);
			}
		}

		previousPosition1.x = 0; previousPosition1.y = 0; previousPosition1.z = 0; previousPosition2.x = 0; previousPosition2.y = 0; previousPosition2.z = 0;
	}

	// Adding bonds for chains
	int nBondsInAmorphousChains = 0;

	for (int i = 0; i < nChainEnds; ++i)
	{
		for (int j = 0; j < (nAtomsPerChain - 1); ++j)
		{
			if (addedPolymer[i][j].x != 0 && addedPolymer[i][j].y != 0 && addedPolymer[i][j].z != 0 && addedPolymer[i][j + 1].x != 0 && addedPolymer[i][j + 1].y != 0 && addedPolymer[i][j + 1].z != 0)
			{
				polymerBonds[i][j].sino = bondIndex + nBondsInAmorphousChains;
				polymerBonds[i][j].bondType = 2;
				polymerBonds[i][j].batom1 = nAtomsWithinBounds_prev + (i * nAtomsPerChain) + j + 1;
				polymerBonds[i][j].batom2 = nAtomsWithinBounds_prev + (i * nAtomsPerChain) + j + 2;

				nBondsInAmorphousChains++;
			}
		}
	}

	int currentAtomNumber = 1;

	// Calculate the boundary
	SIMULATION_BOUNDS finalSimulationBoundary;

	finalSimulationBoundary.xlo = trimmedCrystalAtoms[0].x;
	finalSimulationBoundary.xhi = trimmedCrystalAtoms[0].x;
	finalSimulationBoundary.ylo = trimmedCrystalAtoms[0].y;
	finalSimulationBoundary.yhi = trimmedCrystalAtoms[0].y;
	finalSimulationBoundary.zlo = trimmedCrystalAtoms[0].z;
	finalSimulationBoundary.zhi = trimmedCrystalAtoms[0].z;

	for (int i = 1; i < nAtomsWithinBounds_prev; ++i)
	{
		if (trimmedCrystalAtoms[i].x < finalSimulationBoundary.xlo) {
			finalSimulationBoundary.xlo = trimmedCrystalAtoms[i].x; }
		else if (trimmedCrystalAtoms[i].x > finalSimulationBoundary.xhi) {
			finalSimulationBoundary.xhi = trimmedCrystalAtoms[i].x; }

		if (trimmedCrystalAtoms[i].y < finalSimulationBoundary.ylo) {
			finalSimulationBoundary.ylo = trimmedCrystalAtoms[i].y; }
		else if (trimmedCrystalAtoms[i].y > finalSimulationBoundary.yhi) {
			finalSimulationBoundary.yhi = trimmedCrystalAtoms[i].y; }

		if (trimmedCrystalAtoms[i].z < finalSimulationBoundary.zlo) {
			finalSimulationBoundary.zlo = trimmedCrystalAtoms[i].z; }
		else if (trimmedCrystalAtoms[i].z > finalSimulationBoundary.zhi) {
			finalSimulationBoundary.zhi = trimmedCrystalAtoms[i].z; }
	}

	for (int i = 0; i < nChainEnds; ++i)
	{
		for (int j = 0; j < nAtomsPerChain; ++j)
		{
			if (addedPolymer[i][j].x != 0 && addedPolymer[i][j].y != 0 && addedPolymer[i][j].z != 0)
			{
				if (addedPolymer[i][j].x < finalSimulationBoundary.xlo) {
					finalSimulationBoundary.xlo = addedPolymer[i][j].x; }
				else if (addedPolymer[i][j].x > finalSimulationBoundary.xhi) {
					finalSimulationBoundary.xhi = addedPolymer[i][j].x; }

				if (addedPolymer[i][j].y < finalSimulationBoundary.ylo) {
					finalSimulationBoundary.ylo = addedPolymer[i][j].y; }
				else if (addedPolymer[i][j].y > finalSimulationBoundary.yhi) {
					finalSimulationBoundary.yhi = addedPolymer[i][j].y; }

				if (addedPolymer[i][j].z < finalSimulationBoundary.zlo) {
					finalSimulationBoundary.zlo = addedPolymer[i][j].z; }
				else if (addedPolymer[i][j].z > finalSimulationBoundary.zhi) {
					finalSimulationBoundary.zhi = addedPolymer[i][j].z; }
			}
		}
	}

	// Calculate the number of assigned atoms in amorphous chains
	int nAtomsInAmorphousChains = 0;

	for (int i = 0; i < nChainEnds; ++i)
	{
		for (int j = 0; j < nAtomsPerChain; ++j)
		{
			if (addedPolymer[i][j].x != 0 && addedPolymer[i][j].y != 0 && addedPolymer[i][j].z != 0)
			{
				nAtomsInAmorphousChains++;
			}
		}
	}

	// Printing the header information in data file
	int totalAtomsInDataFile, totalBondsInDataFile;
	totalAtomsInDataFile = nAtomsInAmorphousChains + nAtomsWithinBounds_prev;
	totalBondsInDataFile = nBondsInAmorphousChains + bondIndex - 1;

	// Printing number of atoms, bonds, then types
	fprintf(outputData, "%d atoms\n%d atom types\n%d bonds\n%d bond types\n\n", totalAtomsInDataFile, 5, totalBondsInDataFile, 1);

	// Printing the boundary information
	fprintf(outputData, "%lf %lf xlo xhi\n%lf %lf ylo yhi\n%lf %lf zlo zhi\n", finalSimulationBoundary.xlo, finalSimulationBoundary.xhi, finalSimulationBoundary.ylo, finalSimulationBoundary.yhi, finalSimulationBoundary.zlo, finalSimulationBoundary.zhi);

	fprintf(outputData, "\nMasses\n\n1 14.027\n2 14.027\n3 13.019\n4 14.027\n5 15.035\n");

	fprintf(outputData, "\nAtoms\n\n");

	// Changing the atom types of border crystalline stems
	double finalXLength_crystal = (crystalBounds.xhi - crystalBounds.xlo), finalZLength_crystal = (crystalBounds.zhi - crystalBounds.zlo), finalYLength_crystal = (crystalBounds.yhi - crystalBounds.ylo), yBuffer, xBuffer, zBuffer;
	yBuffer = finalYLength_crystal * 0.1;
	xBuffer = finalXLength_crystal * 0.1;
	zBuffer = finalZLength_crystal * 0.1;

	for (int i = 0; i < nAtomsWithinBounds_prev; ++i) {
		trimmedCrystalAtoms[i].atomType = 1; }

	printf("\nSimulation box length\n~~~~~~~~~~~~~~~~~~~~~\n\nX length: %lf\nY length: %lf\nZ length: %lf\n\n", finalXLength_crystal, finalYLength_crystal, finalZLength_crystal);

	// Check this loop. Atom type 2 is not being assigned properly.
	for (int i = 0; i < nAtomsWithinBounds_prev; ++i)
	{
		//if (finalYLength_crystal > finalXLength_crystal && finalYLength_crystal > finalZLength_crystal)
		//{
			if (trimmedCrystalAtoms[i].y < (crystalBounds.ylo + yBuffer) || trimmedCrystalAtoms[i].y > (crystalBounds.yhi - yBuffer))
			{
				trimmedCrystalAtoms[i].atomType = 2;

				for (int j = 0; j < nAtomsWithinBounds_prev; ++j)
				{
					if (trimmedCrystalAtoms[j].chainID1 == trimmedCrystalAtoms[i].chainID1 && trimmedCrystalAtoms[j].chainID2 == trimmedCrystalAtoms[i].chainID2) {
						trimmedCrystalAtoms[j].atomType = 2; }
				}
			}
		//}
/*		else if (finalXLength_crystal > finalYLength_crystal && finalXLength_crystal > finalZLength_crystal)
		{
			if (trimmedCrystalAtoms[i].x < (crystalBounds.xlo + xBuffer) || trimmedCrystalAtoms[i].x > (crystalBounds.xhi - xBuffer))
			{
				trimmedCrystalAtoms[i].atomType = 2;

				for (int j = 0; j < nAtomsWithinBounds_prev; ++j)
				{
					if (trimmedCrystalAtoms[j].atomType == trimmedCrystalAtoms[i].atomType) {
						trimmedCrystalAtoms[j].atomType = 2; }
				}
			}
		}
		else if (finalZLength_crystal > finalXLength_crystal & finalZLength_crystal > finalYLength_crystal)
		{
			if (trimmedCrystalAtoms[i].z < (crystalBounds.zlo + zBuffer) || trimmedCrystalAtoms[i].z > (crystalBounds.zhi - zBuffer))
			{
				trimmedCrystalAtoms[i].atomType = 2;

				for (int j = 0; j < nAtomsWithinBounds_prev; ++j)
				{
					if (trimmedCrystalAtoms[j].atomType == trimmedCrystalAtoms[i].atomType) {
						trimmedCrystalAtoms[j].atomType = 2; }
				}
			}
		}
*/	}

	// Printing crystal coordinates
	for (int i = 0; i < nAtomsWithinBounds_prev; ++i)
	{
		if (trimmedCrystalAtoms[i].atomType == 1)
		{
			fprintf(output, "C %lf %lf %lf\n", trimmedCrystalAtoms[i].x, trimmedCrystalAtoms[i].y, trimmedCrystalAtoms[i].z);
		}
		else if (trimmedCrystalAtoms[i].atomType == 2)
		{
			fprintf(output, "N %lf %lf %lf\n", trimmedCrystalAtoms[i].x, trimmedCrystalAtoms[i].y, trimmedCrystalAtoms[i].z);
		}

		fprintf(outputData, "%d %d %d %lf %lf %lf\n", currentAtomNumber, 1, trimmedCrystalAtoms[i].atomType, trimmedCrystalAtoms[i].x, trimmedCrystalAtoms[i].y, trimmedCrystalAtoms[i].z);
		currentAtomNumber++;
	}

	// Printing chain coordinates
	for (int i = 0; i < nChainEnds; ++i)
	{
		for (int j = 0; j < nAtomsPerChain; ++j)
		{
			if (addedPolymer[i][j].x != 0 && addedPolymer[i][j].y != 0 && addedPolymer[i][j].z != 0)
			{
				if (addedPolymer[i][j].atomType == 3) {
					fprintf(output, "F %lf %lf %lf\n", addedPolymer[i][j].x, addedPolymer[i][j].y, addedPolymer[i][j].z); }
				else {
					fprintf(output, "O %lf %lf %lf\n", addedPolymer[i][j].x, addedPolymer[i][j].y, addedPolymer[i][j].z); }

				fprintf(outputData, "%d %d %d %lf %lf %lf\n", currentAtomNumber, 2, addedPolymer[i][j].atomType + 2, addedPolymer[i][j].x, addedPolymer[i][j].y, addedPolymer[i][j].z);

				currentAtomNumber++;
			}
		}
	}

	fprintf(outputData, "\nBonds\n\n");

	// Printing bonds in crystalline region
	for (int i = 0; i < bondIndex; ++i)
	{
		fprintf(outputData, "%d %d %d %d\n", crystalBonds[bondIndex].sino, crystalBonds[bondIndex].bondType, crystalBonds[bondIndex].batom1, crystalBonds[bondIndex].batom2);
	}

	// Printing bonds in amorphous region
	for (int i = 0; i < nChainEnds; ++i)
	{
		for (int j = 0; j < (nAtomsPerChain - 1); ++j)
		{
			if (addedPolymer[i][j].x != 0 && addedPolymer[i][j].y != 0 && addedPolymer[i][j].z != 0 && addedPolymer[i][j + 1].x != 0 && addedPolymer[i][j + 1].y != 0 && addedPolymer[i][j + 1].z != 0)
			{
				fprintf(outputData, "%d %d %d %d\n", polymerBonds[i][j].sino, polymerBonds[i][j].bondType, polymerBonds[i][j].batom1, polymerBonds[i][j].batom2);
			}
		}
	}

	// Store the number of tight loops, loose loops, free chains in variables
	//		Generate constrained walk chain from the free points
	//		Generate a single monomer connection for tight loops
	//		Generate a 'n' monomer connection for loose loops

	fclose (output);
	fclose (outputData);
	return 0;
}