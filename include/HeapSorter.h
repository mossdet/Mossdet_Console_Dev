// Copyright 2014, University of Freiburg,
// Algorithmen und Datenstrukturen
// Author: Daniel Lachner <dl115>

#ifndef UEBUNGSBLATT_02_HEAPSORTER_H_
#define UEBUNGSBLATT_02_HEAPSORTER_H_

#include <vector>

/**
 * Class for sorting algorithms
 */
class Sorter {
  public:
  // Sort the given array using heapSort.
  static void heapSort(std::vector<double> &array);
  static void heapify(std::vector<double> &array);
  static void repairHeap(std::vector<double> &array, int n);

  static void heapSort2(std::vector<double> &vecA, std::vector<double> &vecB);
  static void heapify2(std::vector<double> &vecA, std::vector<double> &vecB);
  static void repairHeap2(std::vector<double> &vecA, std::vector<double> &vecB, int arrayLen);
};

#endif  // UEBUNGSBLATT_02_HEAPSORTER_H_
