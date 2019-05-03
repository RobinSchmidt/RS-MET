#include "PrimeNumbers.h"
#include "PrimeArray.h"  // this initializes our primeArray-member

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

int PrimeNumbers::findClosestLowerPrime(int someNumber)
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
  return prime;
 else
  return primeArray[index-1];
}