#ifndef Time_h
#define Time_h

/**

This is a class for conveniently storing time information.

*/

class Time  
{

public:
 
 //---------------------------------------------------------------------------
 // construction/destruction:

	Time();          ///< Default constructor.

	virtual ~Time(); ///< Destructor.

 //---------------------------------------------------------------------------
 // public access-functions:

 int  getYear();
 void setYear(int newYear);
 int  getMonth();
 void setMonth(int newMonth);
 int  getDay();
 void setDay(int newDay);
 int  getHour();
 void setHour(int newHour);
 int  getMinute();
 void setMinute(int newMinute);
 int  getSecond();
 void setSecond(int newSecond);
 int  getMillisecond();
 void setMillisecond(int newMillisecond);

protected:

 // data members:
 int year;
 int month;
 int day;
 int hour;
 int minute;
 int second;
 int millisecond;

};

#endif // Time_h
