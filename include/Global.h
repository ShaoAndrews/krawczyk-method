/*
 * 全局配置文件
*/

#pragma once
extern double PRECISION;
extern double MINISIZE;
extern int newton_kernel;
extern int bezier_kernel;
extern int resultant_kernel;
extern double bezier_time;
extern double resultant_time;
extern double newton_time;
extern double init_time;
extern int unknown_region;
extern double create_tree_time;
extern double connected_time;
#define DEBUG
#ifdef DEBUG
#define PRINT(var) do{\
    std::cout<< "Location: "<<__FILE__<<", Function: "<< __FUNCTION__<<", Line: "<<__LINE__<<std::endl;\
    std::cout <<"Variable name: "<<#var << std::endl;\
    std::cout << var << std::endl;\
    std::cout << ""<<std::endl;\
}while(0)
#else 
#define PRINT(var)
#endif


#ifdef DEBUG
#define PRINTLIST(var) do{\
    std::cout<< "Location: "<<__FILE__<<", Function: "<< __FUNCTION__<<", Line: "<<__LINE__<<std::endl;\
    std::cout <<"Variable name: "<<#var << std::endl;\
    for(int i = 0; i<var.size();++i){\
        std::cout << "id: "<< i<< "    value: "<<var[i]<< std::endl;\
    }\
}while(0)
#else 
#define PRINTLIST(var)
#endif