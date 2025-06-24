/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/file.h to edit this template
 */

/* 
 * File:   ogScheduler.hpp
 * Author: victortrevino
 *
 * Created on April 11, 2023, 7:47 PM
 */

#ifndef OGSCHEDULER_HPP


#include "ogReadKeyMapping.hpp"

typedef struct ogScheduler {
    uint16_t           (*funcSchedule[MAX_SCHEDULE_SIZE])(ogReadKeyMapping *pMapParam); // ogReadsMapper *pMapParam
    uint16_t             nScheduleFunc;
    char                 funcKey[MAX_SCHEDULE_SIZE];
    uint64_t             nCalls[MAX_SCHEDULE_SIZE];
    uint64_t             nMatches[MAX_SCHEDULE_SIZE];
    uint64_t             nFoundFwd[MAX_SCHEDULE_SIZE];
    uint64_t             nFoundRev[MAX_SCHEDULE_SIZE];
    uint64_t             elapsed[MAX_SCHEDULE_SIZE];
    //std::mutex           mtx;
} ogScheduler;

#define OGSCHEDULER_HPP



#endif /* OGSCHEDULER_HPP */

