#ifndef DATASETNUMMAP_H
#define DATASETNUMMAP_H

class DatasetNumMap
{
public:
    DatasetNumMap() 
    {
        dataset_num_start = 0 ;
        dataset_num_other = 0 ;
    } ;

    ~DatasetNumMap() {};

    int dataset_num_start ;
    int dataset_num_other ;
} ;

#endif