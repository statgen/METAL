////////////////////////////////////////////////////////////////////// 
// metal/Main.h 
// (c) 2000-2011 Goncalo Abecasis
// 
// This file is distributed as part of the Goncalo source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile Goncalo.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Friday March 25, 2011
// 
 
#ifndef __METAL_H__
#define __METAL_H__

#include "StringArray.h"
#include "StringHash.h"
#include "MathVector.h"
#include "MathStats.h"
#include "InputFile.h"
#include "IntArray.h"
#include <vector>


class FileSummary {
public:
    FileSummary(FileSummary * pointer = NULL) {
        next = pointer;
    }

    virtual ~FileSummary() {
        if (next != NULL) {
            delete next;
        }
    }

    String filename;
    String header;
    String separators;

    double weight;
    double genomicControl;
    bool   useStrand;

    int markerColumn;
    int weightColumn;
    int pvalueColumn;
    int effectColumn;
    int firstColumn;
    int secondColumn;
    int stderrColumn;
    int freqColumn;
    int strandColumn;
    int chromosomeColumn;
    int positionColumn;
    int minColumns;
    int expectedColumns;

    bool strictColumnCounting;
    bool logTransform;

    StringArray filterLabel;
    IntArray    filterColumn;
    IntArray    filterCondition;
    Vector      filterValue;
    StringArray filterAlternate;
    StringHash  filterSets;
    IntArray    filterCounts;

    int processedMarkers;

    FileSummary * next;
};

// Crash and Control-C handlers
void UserBreak(int);
void OutOfMemory(int);
void SetupCrashHandlers();

// Basic setup
void ClearAll();

// Custom input filtering
#define LESS_THAN              1
#define LESS_THAN_OR_EQUAL     2
#define EQUAL_TO               3
#define GREATER_THAN_OR_EQUAL  4
#define GREATER_THAN           5
#define NOT_EQUAL_TO           6
#define STRING_MATCH           (EQUAL_TO | 0x80)
#define STRING_MISMATCH        (NOT_EQUAL_TO | 0x80)
#define STRING_IN_SET          (0x87)

void AddFilter(StringArray & filter);
int  TranslateCondition(String & condition);
void SetupFilters(StringArray & header);
bool ApplyFilter(StringArray & row);
void FilterSummary();
void ClearFilters();

// Processing of Input Files
void NumbersToLetters(String & al);
void FlipAllele(String & al);
bool FlipAlleles(String & al1, String & al2, double & effect, double & freq);
bool GuessSecondAllele(int marker, String & al1, String & al2);

// Workhorses
void Analyze(bool heterogeneity);
void ProcessFile(String & filename, FileSummary * history);
bool ReProcessFile(FileSummary * history);

// Help!
void ShowHelp(bool startup = false);

// General run control
void RunScript(FILE * file);


#endif

 
