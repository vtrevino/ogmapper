/* 
 * File:   SimpleStatisticsCalculator.hpp
 * Author: victortrevino
 *
 * Created on July 24, 2022, 3:08 PM
 */

#ifndef SIMPLESTATISTICSCALCULATOR_HPP
#define SIMPLESTATISTICSCALCULATOR_HPP

class SimpleStatisticsCalculator {
     int count;   // Number of numbers that have been entered.
     double sum;  // The sum of all the items that have been entered.
     double squareSum;  // The sum of the squares of all the items.
     double max;  // Largest item seen.
     double min;  // Smallest item seen.

public:

    SimpleStatisticsCalculator();
    virtual ~SimpleStatisticsCalculator();
    void        reset();
    void        enter(double num);
    int         getCount();
    double      getSum();
    double      getMean();
    double      getStandardDeviation();
    double      getVariance();
    double      getMin();
    double      getMax();

private:

};


#endif /* SIMPLESTATISTICSCALCULATOR_HPP */
