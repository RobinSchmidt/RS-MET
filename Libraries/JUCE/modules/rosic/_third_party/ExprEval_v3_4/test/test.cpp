// File:    test.cpp
// Author:  Brian Allen Vanderburg II
// Purpose: ExprEval 3.x test program
//------------------------------------------------------------------------------

// Includes
#include <iostream>
#include <ctime>

#include "../expreval.h"

using namespace std;
using namespace ExprEval;

// Test data stuff
class MyData : public Data
    {
    public:
        MyData(int x)
            {
            cout << "Made\n";
            m_x = x;
            }
            
        ~MyData()
            {
            cout << "Killed\n";
            }    
            
        void Print()
            {
            cout << "The data is " << m_x << endl;
            }
            
        long GetType() const
            {
            return CreateType('D', 'A', 'T', 'A');
            }
            
    private:
        int m_x;            
    };
    
    
class TestNode : public FunctionNode
    {
    public:
        TestNode(Expression *expr) : FunctionNode(expr)
            {
            SetArgumentCount(0, 0, 0, 0, 1, 1);
            }
            
        double DoEvaluate()
            {
            Data *d = m_data[0]->GetData();
            
            if(d)
                {
                if(d->GetType() != Data::CreateType('D', 'A', 'T', 'A'))
                    throw(InvalidDataException(GetName()));
                else
                    ((MyData*)d)->Print();    
                }
            else
                {
                m_data[0]->SetData(new MyData(5));
                }    
                
            return 0.0;    
            }   
    };
    
class TestFactory : public FunctionFactory
    {
    public:
        string GetName() const
            {
            return "test";
            }
            
        FunctionNode *DoCreate(Expression *expr)
            {
            return new TestNode(expr);
            }    
    };    
    
class TestNode2 : public FunctionNode
    {
    public:
        TestNode2(Expression *expr) : FunctionNode(expr)
            {
            SetArgumentCount(0, 0, 0, 0, 0, 2);
            }
            
        double DoEvaluate()
            {
            return 123.456;    
            }   
    };
    
class TestFactory2 : public FunctionFactory
    {
    public:
        string GetName() const
            {
            return "test2";
            }
            
        FunctionNode *DoCreate(Expression *expr)
            {
            return new TestNode2(expr);
            }    
    };        
    
void HandleException(Exception &e)
    {
    switch(e.GetType())
        {
        case Exception::Type_SyntaxException:
            cout << "Syntax Error\n";
            break;
            
        case Exception::Type_EmptyExpressionException:
            cout << "Empty Expression\n";
            break;
            
        default:
            cout << "Other Error\n";
            break;        
        }
    }

int main()
    {
    string expr;
    long pos, count;
    
    cout << "Enter expression: ";
    cin >> expr;
    
    
    cout << "Enter count: ";
    cin >> count;
    
    try
        {
        ValueList v;
        FunctionList f;
        DataList d;
        Expression e;
        
        v.AddDefaultValues();
        f.AddDefaultFunctions();
        f.Add(new TestFactory());
        f.Add(new TestFactory2());        
        
        e.SetValueList(&v);
        e.SetFunctionList(&f);
        e.SetDataList(&d);
        
        e.Parse(expr);
        
        double result;
        
        // Start a timer
        time_t start = time(NULL);
        
        for(pos = 0; pos < count; pos++)
            {
            result = e.Evaluate();
            }
            
        // Determine total time
        double total = difftime(time(NULL), start);
        
        cout << "Total time: " << total << endl;
        if(total != 0.0)
            {
            cout << "Average/sec: " << (double)count / total << endl;
            }
            
        // Variable dump
        cout << "Value dump" << endl;
        cout << "-------------------------" << endl;
        
        ValueList::size_type i, c;
        c = v.Count();
        for(i = 0; i < c; i++)
            {
            string name;
            double value;
            
            v.Item(i, &name, &value);
            cout << name.c_str() << " = " << value << endl;
            }
            
        cout << endl;
        }
    catch(Exception &e)
        {
        HandleException(e);
        }
    }
