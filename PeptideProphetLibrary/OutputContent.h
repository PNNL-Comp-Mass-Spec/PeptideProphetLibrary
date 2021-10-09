#ifndef OUTPUTCONTENT_H
#define OUTPUTCONTENT_H

class OutputContent
{
public:
    OutputContent()
    {
        HitNum = 0 ;
        fscore = 0.0 ;
        prob = 0.0 ;
        negonly = 0 ;
    } ;

    ~OutputContent() {};

    int HitNum ;
    float fscore ;
    float prob ;
    int negonly ;
} ;

#endif