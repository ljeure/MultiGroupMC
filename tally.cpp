/** 
 @file      tally.cpp
 @brief     contains functions for the Tally class
 @author    Luke Eure
 @date      January 8 2016
*/

#include "tally.h"

/**
 @brief constructor for Tally class
*/
Tally::Tally() {}

/**
 @brief deconstructor
*/
Tally::~Tally() {}

/**
 @brief add an amount to the tally
 @param amt an amount to be added
*/
void Tally::add(double amt) {
    _tally_count += amt;
}

/**
 @brief clears the tally
*/
void Tally::clear() {
    _tally_count = 0;
}
 
/**
  @brief gets the current tally count
  @return a double: the number stored in the tally
*/
double Tally::getCount() {
    return _tally_count;
}
