#include "rosic_PrimeNumbers.h"
#include "rosic_PrimeArray.h"     // this initializes our primeArray-member
using namespace rosic;

PrimeNumbers::PrimeNumbers()
{

}

PrimeNumbers::~PrimeNumbers()
{

}

bool PrimeNumbers::isPrime(int someNumber)
{
 // initializations:
 int lowIndex  = 0;
 int highIndex = primeArrayLength-1;
 int index     = (lowIndex+highIndex)/2;
 int prime     = primeArray[index];

 // binary search:
 while( lowIndex <= highIndex && prime != someNumber )
 {
  index = (lowIndex+highIndex)/2;
  prime = primeArray[index];
  if( prime > someNumber )
   highIndex = index-1; // check only the lower portion in the next iteration
  else if( prime < someNumber )
   lowIndex = index+1;  // check only the upper portion in the next iteration
 }

 if( prime == someNumber )
  return true;
 else
  return false;
}

int PrimeNumbers::findClosestLowerPrimeIndex(int someNumber)
{
 // initializations:
 int lowIndex  = 0;
 int highIndex = primeArrayLength-1;
 int index     = (lowIndex+highIndex)/2;
 int prime     = primeArray[index];
 int result    = 0;

 // binary search:
 while( lowIndex <= highIndex && prime != someNumber )
 {
  index = (lowIndex+highIndex)/2;
  prime = primeArray[index];
  if( prime > someNumber )
   highIndex = index-1; // check only the lower portion in the next iteration
  else if( prime < someNumber )
   lowIndex = index+1;  // check only the upper portion in the next iteration
 }

 if( prime <= someNumber )
  return index;
 else
  return index-1;
}

int PrimeNumbers::findClosestLowerPrime(int someNumber)
{
  return primeArray[findClosestLowerPrimeIndex(someNumber)];
}

int PrimeNumbers::findClosestPrime(int someNumber)
{
  int iLow  = findClosestLowerPrimeIndex(someNumber);
  int pLow  = primeArray[iLow];
  int pHigh = primeArray[iLow+1];
  if( abs(pHigh-someNumber) < abs(someNumber-pLow) )
    return pHigh;
  else
    return pLow;
}