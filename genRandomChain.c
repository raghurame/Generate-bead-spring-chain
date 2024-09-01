#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <stdbool.h>

typedef struct position
{
	int molType;
	int atomType;
	double x, y, z;
} BEAD_POSITIONS;

typedef struct bonds
{
	int id, atom1, atom2, bondType;
} BONDS;

double computeDistance (double x1, double y1, double z1, double x2, double y2, double z2)
{
	double distance;

	distance = sqrt (
		(x2 - x1) * (x2 - x1) +
		(y2 - y1) * (y2 - y1) +
		(z2 - z1) * (z2 - z1)
		);

	return distance;
}

BEAD_POSITIONS nextRandomWalkPosition (BEAD_POSITIONS previousPosition, BEAD_POSITIONS nextPosition, float stepLength)
{
	float randomAngle1, randomAngle2;

	srand (time(NULL));

	randomAngle1 = rand() % 360;
	randomAngle2 = rand() % 360;
	randomAngle1 = randomAngle1 / 57.2958;
	randomAngle2 = randomAngle2 / 57.2958;

	nextPosition.x = previousPosition.x + cos (randomAngle2) * (double)stepLength * cos (randomAngle1);
	nextPosition.y = previousPosition.y + cos (randomAngle2) * (double)stepLength * sin (randomAngle1);
	nextPosition.z = previousPosition.z + (double)stepLength * sin (randomAngle2);

	return nextPosition;
}

BEAD_POSITIONS generateRandomInitialPositions (BEAD_POSITIONS initPosition, float xlo, float xhi, float ylo, float yhi, float zlo, float zhi)
{
	float boxLength_x = xhi - xlo, boxLength_y = yhi - ylo, boxLength_z = zhi - zlo;

	sleep (1);
	srand (time(NULL));

	initPosition.x = fmod ((double)rand (), (double)boxLength_x);
	initPosition.y = fmod ((double)rand (), (double)boxLength_y);
	initPosition.z = fmod ((double)rand (), (double)boxLength_z);

	initPosition.x += xlo;
	initPosition.y += ylo;
	initPosition.z += zlo;

	return initPosition;
}

BEAD_POSITIONS **resetAtomTypes (BEAD_POSITIONS **polymerBeads, int i, int nAtomsPerChain)
{
	for (int j = 0; j < nAtomsPerChain; ++j)
	{
		polymerBeads[i][j].atomType = 0;
	}

	return polymerBeads;
}

BEAD_POSITIONS **setAtomTypes (BEAD_POSITIONS **polymerBeads, int nChains, int nAtomsPerChain)
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

BEAD_POSITIONS **generateAtomTypes (BEAD_POSITIONS **polymerBeads, int nChains, int nAtomsPerChain, int nMonomers, float nSideGroupsPerChain)
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

BONDS **generateBonds (BONDS **polymerBonds, BEAD_POSITIONS **polymerBeads, int nChains, int nAtomsPerChain, int nBondsPerChain)
{
	// id, atom1, atom2, bondType;
	int currentBondID = 0, currentAtomID = 0;

	for (int i = 0; i < nChains; ++i)
	{
		currentBondID = 0;

		for (int j = 0; j < nAtomsPerChain; ++j)
		{
			if (polymerBeads[i][j].atomType > 0) {
				currentAtomID++; }

			if (polymerBeads[i][j].atomType == 1 && j != 0)
			{
				polymerBonds[i][currentBondID].atom1 = currentAtomID - 1;
				polymerBonds[i][currentBondID].atom2 = currentAtomID;
				polymerBonds[i][currentBondID].bondType = 1;
				currentBondID++;
			}

			if (polymerBeads[i][j].atomType == 2) {
				polymerBonds[i][currentBondID].atom1 = currentAtomID - 1;
				polymerBonds[i][currentBondID].atom2 = currentAtomID;
				polymerBonds[i][currentBondID].bondType = 1;
				currentBondID++; }

			else if (polymerBeads[i][j].atomType == 3) {
				polymerBonds[i][currentBondID].atom1 = currentAtomID - 2;
				polymerBonds[i][currentBondID].atom2 = currentAtomID;
				polymerBonds[i][currentBondID].bondType = 2;
				currentBondID++; }
		}
	}

	return polymerBonds;
}

BONDS **resetPolymerBonds (BONDS **polymerBonds, int nChains, int nBondsPerChain)
{
	for (int i = 0; i < nChains; ++i)
	{
		for (int j = 0; j < nBondsPerChain; ++j)
		{
			polymerBonds[i][j].atom1 = 0;
			polymerBonds[i][j].atom2 = 0;
			polymerBonds[i][j].bondType = 0;
		}
	}

	return polymerBonds;
}

BEAD_POSITIONS rotateInXY (BEAD_POSITIONS nextPosition, float theta_xy)
{
	float tempX, tempY, tempZ;
	tempX = nextPosition.x;
	tempY = nextPosition.y;
	tempZ = nextPosition.z;

	nextPosition.x = tempX * cosf (theta_xy) - tempY * sinf (theta_xy);
	nextPosition.y = tempX * sinf (theta_xy) + tempY * cosf (theta_xy);
	nextPosition.z = tempZ;

	return nextPosition;
}

BEAD_POSITIONS rotateInYZ (BEAD_POSITIONS nextPosition, float theta_yz)
{
	float tempX, tempY, tempZ;
	tempX = nextPosition.x;
	tempY = nextPosition.y;
	tempZ = nextPosition.z;

	nextPosition.x = tempX;
	nextPosition.y = tempY * cosf (theta_yz) - tempZ * sinf (theta_yz);
	nextPosition.z = tempY * sinf (theta_yz) + tempZ * cosf (theta_yz);

	return nextPosition;
}

BEAD_POSITIONS rotateInXZ (BEAD_POSITIONS nextPosition, float theta_xz)
{
	float tempX, tempY, tempZ;
	tempX = nextPosition.x;
	tempY = nextPosition.y;
	tempZ = nextPosition.z;

	nextPosition.x = tempX * cosf (theta_xz) + tempZ * sinf (theta_xz);
	nextPosition.y = tempY;
	nextPosition.z = - tempX * sinf (theta_xz) + tempZ * cosf (theta_xz);

	return nextPosition;
}

BEAD_POSITIONS rotateNewBond (BEAD_POSITIONS nextPosition, BEAD_POSITIONS previousPosition1, BEAD_POSITIONS previousPosition2, float bondLength)
{
	float unitX, unitY, unitZ;

	unitX = previousPosition2.x - previousPosition1.x;
	unitY = previousPosition2.y - previousPosition1.y;
	unitZ = previousPosition2.z - previousPosition1.z;

	unitX /= bondLength;
	unitY /= bondLength;
	unitZ /= bondLength;

	float norm1X = - unitY * cosf (10 / 57.2958), norm1Y = unitX * cosf (10 / 57.2958), norm1Z = unitZ * sinf (10 / 57.2958);
	float norm2X = unitY * unitZ * sinf (10 / 57.2958) - unitZ * unitX * cosf (10 / 57.2958), norm2Y = -1 * (unitX * unitZ * sinf (10 / 57.2958) + unitZ * unitY * cosf (10 / 57.2958)), norm2Z = unitX * unitX * cosf (10 / 57.2958) + unitY * unitY * cosf (10 / 57.2958);

	float tempX, tempY, tempZ;
	tempX = nextPosition.x;
	tempY = nextPosition.y;
	tempZ = nextPosition.z;

	nextPosition.x = unitX * tempX + norm1X * tempY + norm2X * tempZ;
	nextPosition.y = unitY * tempX + norm1Y * tempY + norm2Y * tempZ;
	nextPosition.z = unitZ * tempX + norm1Z * tempY + norm2Z * tempZ;

	return nextPosition;
}

BEAD_POSITIONS nextConstrainedWalkPosition (BEAD_POSITIONS previousPosition1, BEAD_POSITIONS previousPosition2, BEAD_POSITIONS nextPosition, float bondLength, float bondAngle, int j)
{
	srand (time(NULL) + j);

	BEAD_POSITIONS angleOfPreviousBond;

	bondAngle = 180 - bondAngle;
	bondAngle /= 57.2958;

	/*
	beads[currentChain][i].x = beads[currentChain][i - 1].x + cos (randomAngle2) * (double)bondLength * cos (randomAngle1);
	beads[currentChain][i].y = beads[currentChain][i - 1].y + cos (randomAngle2) * (double)bondLength * sin (randomAngle1);
	beads[currentChain][i].z = beads[currentChain][i - 1].z + (double)bondLength * sin (randomAngle2);
	*/

	float randomAngle;
	randomAngle = rand() % 360;

	nextPosition.x = cos (bondAngle) * bondLength;
	nextPosition.y = sin (bondAngle) * bondLength;// * cos (randomAngle);
	// nextPosition.z = sin (randomAngle) * bondLength;

	float tempRandY = nextPosition.y * 2;
	float randomY = ((float)rand () / (float)RAND_MAX);
	randomY *= tempRandY;
	randomY -= nextPosition.y;
	// printf("%f %f\n", nextPosition.y, randomY);
	nextPosition.z = sqrt (nextPosition.y * nextPosition.y - randomY * randomY);
	nextPosition.y = randomY;

	// Find the angle of the projection of the previous bond on XY, YZ, and XZ planes
	float dy, dx, dz;
	// XY plane
	float slope_xy, theta_xy;

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
	float slope_yz, theta_yz;

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
	float slope_xz, theta_xz;

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

	return nextPosition;
}

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

BEAD_POSITIONS **generateAtomTypesFromGro (const char *filename, int nSkipLines, BEAD_POSITIONS **polymerBeads, int nChains, int nAtomsPerChain)
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
	/*
	to do here:
		364 chains, with 182 methylene C3 groups randomly placed, 0.5 per chain
		box length = x: 5 to 95; y: 5 to 65; z: 37 to 160
	*/
	if (argc != 13)
	{
		printf("REQUIRED ARGUMENTS:\n~~~~~~~~~~~~~~~~~~~\n\n {~} argv[0] = program\n {~} argv[1] = (int)Number of monomers\n {~} argv[2] = (int)Number of chains\n {~} argv[3] = (float)bond length (for all bonds)\n {~} argv[4] = (int)Number of side groups per chain\n {~} argv[5] = (float)xlo\n {~} argv[6] = (float)xhi\n {~} argv[7] = (float)ylo\n {~} argv[8] = (float)yhi\n {~} argv[9] = (float)zlo\n {~} argv[10] = (float)zhi\n {~} argv[11] = (float)average bond angle\n {~} argv[12] = (int) spread in bond angle (this is not stdev or stderr)\n\n");
		exit (1);
	}

	FILE *output, *outputXYZ;
	output = fopen ("chain.data", "w");
	outputXYZ = fopen ("chain.xyz", "w");

	int nMonomers = atoi (argv[1]), nChains = atoi (argv[2]);
	float bondLength = atof (argv[3]), nSideGroupsPerChain = atof (argv[4]), xlo = atof (argv[5]), xhi = atof (argv[6]), ylo = atof (argv[7]), yhi = atof (argv[8]), zlo = atof (argv[9]), zhi = atof (argv[10]), averageAngle = atof (argv[11]);

	int spread = atoi (argv[12]);

	int nAtomsPerChain = (nMonomers * 2) + ceil (nSideGroupsPerChain);

	countNChainsNAtomsFromGro ("../after_npt1.gro", &nAtomsPerChain, &nChains);
	printf("nAtomsPerChain: %d\nnChains: %d\n", nAtomsPerChain, nChains);

	BEAD_POSITIONS **polymerBeads;
	polymerBeads = (BEAD_POSITIONS **) malloc (nChains * sizeof (BEAD_POSITIONS *));

	for (int i = 0; i < nChains; ++i)
	{
		polymerBeads[i] = (BEAD_POSITIONS *) malloc (nAtomsPerChain * sizeof (BEAD_POSITIONS));
	}

	// polymerBeads = generateAtomTypes (polymerBeads, nChains, nAtomsPerChain, nMonomers, nSideGroupsPerChain);
	polymerBeads = generateAtomTypesFromGro ("../after_npt1.gro", 7202, polymerBeads, nChains, nAtomsPerChain);

	int nBondsPerChain = ((nMonomers * 2) - 1) + (nSideGroupsPerChain * nChains);
	BONDS **polymerBonds;

	polymerBonds = (BONDS **) malloc (nChains * sizeof (BONDS *));

	for (int i = 0; i < nChains; ++i)
	{
		polymerBonds[i] = (BONDS *) malloc (nBondsPerChain * sizeof (BONDS));
	}

	polymerBonds = resetPolymerBonds (polymerBonds, nChains, nBondsPerChain);

	BEAD_POSITIONS previousPosition1, previousPosition2, nextPosition;

	for (int i = 0; i < nChains; ++i)
	{
		printf("Generating coordinates for chain %d...                    \r", i + 1);
		fflush (stdout);

		previousPosition1 = generateRandomInitialPositions (previousPosition1, xlo + 10, xhi - 10, ylo + 10, yhi - 10, zlo + 10, zhi - 10);
		previousPosition2 = nextRandomWalkPosition (previousPosition1, previousPosition2, bondLength);
		fprintf(outputXYZ, "C %f %f %f\n", previousPosition2.x, previousPosition2.y, previousPosition2.z);

		polymerBeads[i][0].x = previousPosition1.x;
		polymerBeads[i][0].y = previousPosition1.y;
		polymerBeads[i][0].z = previousPosition1.z;

		polymerBeads[i][1].x = previousPosition2.x;
		polymerBeads[i][1].y = previousPosition2.y;
		polymerBeads[i][1].z = previousPosition2.z;

		for (int j = 1; j < nAtomsPerChain; ++j)
		{
			if (polymerBeads[i][j].atomType > 0)
			{
				if (polymerBeads[i][j].atomType == 3)
				{
					// nextPosition = nextRandomWalkPosition (previousPosition2, nextPosition, bondLength);
					nextPosition = nextConstrainedWalkPosition (previousPosition1, previousPosition2, nextPosition, bondLength, 90, j + i);
					fprintf(outputXYZ, "N %f %f %f\n", j, nextPosition.x, nextPosition.y, nextPosition.z);
					polymerBeads[i][j].x = nextPosition.x;
					polymerBeads[i][j].y = nextPosition.y;
					polymerBeads[i][j].z = nextPosition.z;
				}
				else
				{
					// nextPosition = nextRandomWalkPosition (previousPosition2, nextPosition, bondLength);
					nextPosition = nextConstrainedWalkPosition (previousPosition1, previousPosition2, nextPosition, bondLength, 120, j + i);
					fprintf(outputXYZ, "C %f %f %f\n", j, nextPosition.x, nextPosition.y, nextPosition.z);
					polymerBeads[i][j].x = nextPosition.x;
					polymerBeads[i][j].y = nextPosition.y;
					polymerBeads[i][j].z = nextPosition.z;

					previousPosition1 = previousPosition2;
					previousPosition2 = nextPosition;
				}
			}
		}
	}

/*	polymerBonds = generateBonds (polymerBonds, polymerBeads, nChains, nAtomsPerChain, nBondsPerChain);

	for (int i = 0; i < nChains; ++i)
	{
		for (int j = 0; j < nAtomsPerChain; ++j)
		{
			if (polymerBeads[i][j].atomType > 0)
			{
				printf("%f %f %f\n", polymerBeads[i][j].x, polymerBeads[i][j].y, polymerBeads[i][j].z);
				usleep (10000);
			}
		}
	}
	exit (1); */

	fclose (output);
	fclose (outputXYZ);
	return 0;
}

/*2024-05-14 2024-05-15 2024-05-16 2024-05-17 2024-05-18 2024-05-19 2024-05-20 */

/*

7.8T    coe-rlarson-daily_2024-05-14_01-30
8.4T    coe-rlarson-daily_2024-05-15_01-30
8.4T    coe-rlarson-daily_2024-05-16_01-30
8.4T    coe-rlarson-daily_2024-05-17_01-30
8.9T    coe-rlarson-daily_2024-05-18_01-30
5.5T    coe-rlarson-daily_2024-05-19_01-30
5.5T    coe-rlarson-daily_2024-05-20_01-30

*/

/*mkdir cumv1 cumv2 cumv3 cumv4; 

cat r1/viscosity.output >> cumv1/viscosity.output; cat r2/viscosity.output >> cumv1/viscosity.output; cat r3/viscosity.output >> cumv1/viscosity.output; cat r4/viscosity.output >> cumv1/viscosity.output; cat r5/viscosity.output >> cumv1/viscosity.output; cat r6/viscosity.output >> cumv1/viscosity.output; cat r7/viscosity.output >> cumv1/viscosity.output; cat r8/viscosity.output >> cumv1/viscosity.output; cat r9/viscosity.output >> cumv1/viscosity.output; cat r10/viscosity.output >> cumv1/viscosity.output; 

cat r1/viscosity2.output >> cumv2/viscosity2.output; cat r2/viscosity2.output >> cumv2/viscosity2.output; cat r3/viscosity2.output >> cumv2/viscosity2.output; cat r4/viscosity2.output >> cumv2/viscosity2.output; cat r5/viscosity2.output >> cumv2/viscosity2.output; cat r6/viscosity2.output >> cumv2/viscosity2.output; cat r7/viscosity2.output >> cumv2/viscosity2.output; cat r8/viscosity2.output >> cumv2/viscosity2.output; cat r9/viscosity2.output >> cumv2/viscosity2.output; cat r10/viscosity2.output >> cumv2/viscosity2.output; 

cat r1/viscosity3.output >> cumv3/viscosity3.output; cat r2/viscosity3.output >> cumv3/viscosity3.output; cat r3/viscosity3.output >> cumv3/viscosity3.output; cat r4/viscosity3.output >> cumv3/viscosity3.output; cat r5/viscosity3.output >> cumv3/viscosity3.output; cat r6/viscosity3.output >> cumv3/viscosity3.output; cat r7/viscosity3.output >> cumv3/viscosity3.output; cat r8/viscosity3.output >> cumv3/viscosity3.output; cat r9/viscosity3.output >> cumv3/viscosity3.output; cat r10/viscosity3.output >> cumv3/viscosity3.output; 

cat r1/viscosity4.output >> cumv4/viscosity4.output; cat r2/viscosity4.output >> cumv4/viscosity4.output; cat r3/viscosity4.output >> cumv4/viscosity4.output; cat r4/viscosity4.output >> cumv4/viscosity4.output; cat r5/viscosity4.output >> cumv4/viscosity4.output; cat r6/viscosity4.output >> cumv4/viscosity4.output; cat r7/viscosity4.output >> cumv4/viscosity4.output; cat r8/viscosity4.output >> cumv4/viscosity4.output; cat r9/viscosity4.output >> cumv4/viscosity4.output; cat r10/viscosity4.output >> cumv4/viscosity4.output; 
*/
