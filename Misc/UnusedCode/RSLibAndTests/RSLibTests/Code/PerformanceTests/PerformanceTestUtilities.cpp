#include "PerformanceTestUtilities.h"

#include <time.h>

#define S(x) std::string(x);

static const int nameColumnWidth = 32;

// maybe move this to RSCore (and use an rsString there and convert here to std::string)
std::string getDateAndTimeString()
{
  time_t t;
  time(&t);
  struct tm *pTime;
  pTime = localtime(&t); // returns pointer to statically allocated memory 
                         // -> we shall not deallocate it via free()
  char tmp[5];

  std::string s = S("Date: ");
  sprintf(tmp, "%d", pTime->tm_year + 1900);
  s += S(tmp);
  s += "/";
  sprintf(tmp, "%02d", pTime->tm_mon + 1);
  s += S(tmp);
  s += "/";
  sprintf(tmp, "%02d", pTime->tm_mday);
  s += S(tmp);

  s += S(", Time: ");
  sprintf(tmp, "%02d", pTime->tm_hour);
  s += S(tmp);
  s += ":";
  sprintf(tmp, "%02d", pTime->tm_min);
  s += S(tmp);
  s += ":";
  sprintf(tmp, "%02d", pTime->tm_sec);
  s += S(tmp);

  return s;
}

std::string createTableHead()
{
  //int nameColumnWidth = 24;

  rsString tmp("Name");
  tmp.padToLength(nameColumnWidth);

  tmp += "Cycles";

  std::string s  = tmp.getAsStdString(); 
  return s;
}

std::string createPerformanceTestHeader()
{
  std::string s = "Performance Test Results\n";

  s += getDateAndTimeString();
  s += "\n";

  // todo:
  // s += getCompilerString();
  // s += getConfigurationString();
  // s += getOperatingSytemString();

  //
  s += "\n";
  s += createTableHead();
  s += "\n";

  return s;
}

void appendResultToReport(std::string &reportString, const std::string &nameOfTest, 
                          const double result)
{
  rsString tmp(nameOfTest);
  tmp.padToLength(nameColumnWidth, '.');

  tmp += rsString(result); // let this constructor have an optional 2nd parameter that determines
                           // the number of digits

  tmp += "\n";
  reportString += tmp.getAsStdString();
}

