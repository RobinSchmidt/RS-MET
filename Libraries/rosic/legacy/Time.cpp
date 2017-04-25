#include "Time.h"

Time::Time()
{
 year        = 0;
 month       = 0;
 day         = 0;
 hour        = 0;
 minute      = 0;
 second      = 0;
 millisecond = 0;
}

Time::~Time()
{

}

int Time::getYear()
{
 return year;
}

void Time::setYear(int newYear)
{
 year = newYear;
}

int Time::getMonth()
{
 return month;
}

void Time::setMonth(int newMonth)
{
 if( newMonth >= 1 && newMonth <= 12 )
  month = newMonth;
 else
 {
  // throw exception ...
 }
}

int Time::getDay()
{
 return day;
}

void Time::setDay(int newDay)
{
 if( newDay >= 1 && newDay <= 7 )
  day = newDay;
 else
 {
  // throw exception ...
 }
}

int Time::getHour()
{
 return hour;
}

void Time::setHour(int newHour)
{
 if( newHour >= 0 && newHour <= 23 )
  hour = newHour;
 else
 {
  // throw exception ...
 }
}

int Time::getMinute()
{
 return minute;
}

void Time::setMinute(int newMinute)
{
 if( newMinute >= 0 && newMinute <= 59 )
  minute = newMinute;
 else
 {
  // throw exception ...
 }
}

int Time::getSecond()
{
 return second;
}

void Time::setSecond(int newSecond)
{
 if( newSecond >= 0 && newSecond <= 59 )
  second = newSecond;
 else
 {
  // throw exception ...
 }
}

int Time::getMillisecond()
{
 return millisecond;
}

void Time::setMillisecond(int newMillisecond)
{
 if( newMillisecond >= 0 && newMillisecond <= 999 )
  millisecond = newMillisecond;
 else
 {
  // throw exception ...
 }
}