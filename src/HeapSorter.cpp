// Copyright 2014, University of Freiburg,
// Algorithmen und Datenstrukturen
// Author: Daniel Lachner <dl115>

#include <vector>
#include <iostream>
#include <assert.h>

#include "./HeapSorter.h"

/**
 * Sort the given array using heapSort.
 */
void Sorter::heapSort(std::vector<double> &array) {
	int n = array.size();
	int i;
	double temp;

	Sorter::heapify(array);
	for (i = n - 1; i > 0; i--) {
		temp = array[0];
		array[0] = array[i];
		array[i] = temp;
		Sorter::repairHeap(array, i);
	}
}

void Sorter::heapify(std::vector<double> &array) {
	int n = array.size();
	int k, child, parent;
	double item;

	for (k = 1; k < n; k++) {
		item = array[k];
		child = k;
		parent = (child - 1) / 2;

		while ((child > 0) && (item > array[parent])) {
			array[child] = array[parent];
			child = parent;
			parent = (child - 1) / 2;
		}
		array[child] = item;
	}
}

void Sorter::repairHeap(std::vector<double> &array, int arrayLen) {
	int i, j;
	double item;

	if (arrayLen > 0) {
		j = 0;
		item = array[j];
		i = 2 * j + 1;

		while (i <= arrayLen - 1) {
			if (i + 1 <= arrayLen - 1)
				if (array[i] < array[i + 1])
					i++;
			if (item < array[i]) {
				array[j] = array[i];
				j = i;
				i = 2 * j + 1;
			}
			else {
				break;
			}
		}
		array[j] = item;
	}
}


// Sort vectorA into ascending order using HeapSort, while making the corresponding rearrangement of vectorB.
void Sorter::heapSort2(std::vector<double> &vecA, std::vector<double> &vecB) {
	assert(vecA.size() == vecB.size());
	int n = vecA.size();
	int i;
	double tempA, tempB;

	Sorter::heapify2(vecA, vecB);
	for (i = n - 1; i > 0; i--) {
		tempA = vecA[0];
		tempB = vecB[0];
		vecA[0] = vecA[i];
		vecB[0] = vecB[i];

		vecA[i] = tempA;
		vecB[i] = tempB;

		Sorter::repairHeap2(vecA, vecB, i);
	}
}

void Sorter::heapify2(std::vector<double> &vecA, std::vector<double> &vecB) {
	int n = vecA.size();
	int k, child, parent;
	double itemA, itemB;

	for (k = 1; k < n; k++) {
		itemA = vecA[k];
		itemB = vecB[k];
		child = k;
		parent = (child - 1) / 2;

		while ((child > 0) && (itemA > vecA[parent])) {
			vecA[child] = vecA[parent];
			vecB[child] = vecB[parent];
			child = parent;
			parent = (child - 1) / 2;
		}
		vecA[child] = itemA;
		vecB[child] = itemB;
	}
}

void Sorter::repairHeap2(std::vector<double> &vecA, std::vector<double> &vecB, int arrayLen) {
	int i, j;
	double itemA, itemB;

	if (arrayLen > 0) {
		j = 0;
		itemA = vecA[j];
		itemB = vecB[j];
		i = 2 * j + 1;

			while (i <= arrayLen - 1) {
				if (i + 1 <= arrayLen - 1)
					if (vecA[i] < vecA[i + 1])
						i++;
				if (itemA < vecA[i]) {
					vecA[j] = vecA[i];
					vecB[j] = vecB[i];
					j = i;
					i = 2 * j + 1;
				}
				else {
					break;
				}
			}
		vecA[j] = itemA;
		vecB[j] = itemB;
	}
}