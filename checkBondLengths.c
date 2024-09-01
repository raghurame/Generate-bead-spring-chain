#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef struct inputdata
{
	float x, y, z, vx, vy, vz;
	int sino, monomerID;
	char atomType[10], resIndex[10];
	int atomType2;
} GRO_DATA;

float calculateDistance (float x1, float y1, float z1, float x2, float y2, float z2, float distance)
{
	distance = sqrt ((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
	return distance;
}

int main(int argc, char const *argv[])
{
	FILE *input;
	input = fopen ("../after_npt1.gro", "r");

	char lineString[5000];
	int nAtoms;

	for (int i = 0; i < 7202; ++i)
	{
		fgets (lineString, 3000, input);

		if (i == 1) {
			sscanf (lineString, "%d\n", &nAtoms); }
	}

	printf("Number of atoms: %d\n", nAtoms);

	int currentLine = 0;

	float distance;

	GRO_DATA *inputData;
	inputData = (GRO_DATA *) malloc (nAtoms * sizeof (GRO_DATA));

	while (fgets (lineString, 4000, input) != NULL)
	{
		sscanf (lineString, "%d%s %s %d %f %f %f %f %f %f\n", &inputData[currentLine].monomerID, &inputData[currentLine].resIndex, &inputData[currentLine].atomType, &inputData[currentLine].sino, &inputData[currentLine].x, &inputData[currentLine].y, &inputData[currentLine].z, &inputData[currentLine].vx, &inputData[currentLine].vy, &inputData[currentLine].vz);

		if (strstr (inputData[currentLine].atomType, "C1"))
		{
			inputData[currentLine].atomType2 = 1;
		}
		else if (strstr (inputData[currentLine].atomType, "C2"))
		{
			inputData[currentLine].atomType2 = 2;
		}
		else if (strstr (inputData[currentLine].atomType, "C3"))
		{
			inputData[currentLine].atomType2 = 3;
		}

		currentLine++;
	}

	printf("currentLine: %d\n", currentLine);

	for (int i = 1; i < currentLine; ++i)
	{
		if (inputData[i].monomerID == inputData[i - 1].monomerID && inputData[i - 1].atomType2 == 1 && inputData[i].atomType2 == 2)
		{
			// printf("%d %d %f %f %f\n", inputData[i - 1].monomerID, inputData[i - 1].resIndex, inputData[i - 1].atomType2, inputData[i - 1].sino, inputData[i - 1].x, inputData[i - 1].y, inputData[i - 1].z);
			// printf("%d %d %f %f %f\n", inputData[i].monomerID, inputData[i].resIndex, inputData[i].atomType2, inputData[i].sino, inputData[i].x, inputData[i].y, inputData[i].z);

			distance = calculateDistance (inputData[i - 1].x, inputData[i - 1].y, inputData[i - 1].z, inputData[i].x, inputData[i].y, inputData[i].z, distance);
			printf("%f\n", distance);

			// printf("\n\n\n");
			usleep (100000);
		}
	}

	fclose (input);
	return 0;
}