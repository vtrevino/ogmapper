/* 
 * File:   SimpleStatisticsCalculator.cpp
 * Author: victortrevino
 * 
 * Created on July 24, 2022, 3:08 PM
 */

#include <float.h>
#include <math.h>
#include "SimpleStatisticsCalculator.hpp"

SimpleStatisticsCalculator::SimpleStatisticsCalculator() {
    reset();
}

SimpleStatisticsCalculator::~SimpleStatisticsCalculator() {
}

void SimpleStatisticsCalculator::reset() {
    count = 0;
    sum = 0;
    squareSum = 0;
    max = DBL_MIN;
    min = DBL_MAX;
}

void SimpleStatisticsCalculator::enter(double num) {
      // Add the number to the dataset.
   count++;
   sum += num;
   squareSum += num*num;
   if (num > max)
      max = num;
   if (num < min)
      min = num;
}


int SimpleStatisticsCalculator::getCount() {
      // Return number of items that have been entered.
   return count;
}

double SimpleStatisticsCalculator::getSum() {
      // Return the sum of all the items that have been entered.
   return sum;
}

double SimpleStatisticsCalculator::getMean() {
      // Return average of all the items that have been entered.
      // Value is Double.NaN if count == 0.
   return sum / count;
}

double SimpleStatisticsCalculator::getStandardDeviation() {
    if (count < 2) return 0;
     // Return standard deviation of all the items that have been entered.
     // Value will be Double.NaN if count == 0.
   double mean = getMean();
   // return Math.sqrt( squareSum/count - mean*mean );  //// this is biased
   return sqrt( squareSum/(count-1) - count*mean*mean/(count-1) );
}

double SimpleStatisticsCalculator::getVariance() {
    double sd = getStandardDeviation();
    return sd * sd;
}

double SimpleStatisticsCalculator::getMin() {
     // Return the smallest item that has been entered.
     // Value will be infinity if no items have been entered.
   return min;
}

double SimpleStatisticsCalculator::getMax() {
     // Return the largest item that has been entered.
     // Value will be -infinity if no items have been entered.
   return max;
}
