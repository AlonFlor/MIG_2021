// MIGBasicSim.h : Include file for standard system include files,
// or project specific include files.

#pragma once

#include <iostream>



float* diff(float X1[], float X2[])
{
	float ans[2];
	ans[0] = X2[0] - X1[0];
	ans[1] = X2[1] - X1[1];
	return ans;
}

float dist(float diffs[])
{
	if (!((diffs[0] * diffs[0]) + (diffs[1] * diffs[1]) >= 0)) {
		printf("Less than 0\n");
		printf("%f\n", (diffs[0] * diffs[0]) + (diffs[1] * diffs[1]));
		exit(1);
	}
	return sqrt((diffs[0] * diffs[0]) + (diffs[1] * diffs[1]));
}
