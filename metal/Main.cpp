////////////////////////////////////////////////////////////////////// 
// metal/Main.cpp 
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
 
#include "Main.h"

#include <math.h>
#include <ctype.h>
#include <signal.h>
#include <MathGold.h>
#include <map>
#include "../version/VersionInfo.h"

StringIntHash markerLookup;

Vector statistics;
Vector weights;
Vector samples;
Vector frequencies;
Vector frequencies2;
Vector minFrequencies;
Vector maxFrequencies;

Vector *    custom = NULL;
StringArray customVariables;
StringArray customLabels;

Vector    hetStatistic;
IntArray  hetDegreesOfFreedom;

StringArray allele1;
StringArray allele2;
String      original_flipped;

StringArray chromosomes;
StringArray positions;

StringArray filenames;
StringArray directions;
FileSummary * processedFiles = NULL;

double weight = 1.0;
double minweight = 1.0;

String markerLabel  = "MARKER";
String weightLabel  = "N";
String pvalueLabel  = "PVALUE";
String effectLabel  = "EFFECT";
String stderrLabel  = "STDERR";
String strandLabel  = "STRAND";
String frequencyLabel = "FREQ";
String firstAllele  = "ALLELE1";
String secondAllele = "ALLELE2";
String chromosomeLabel = "CHROMOSOME";
String positionLabel = "POSITION";

String separators  = " \t";

bool   useStandardErrors = false;
bool   useStrand = false;
bool   averageFrequencies = false;
bool   minMaxFrequencies = false;
bool   genomicControl = false;
bool   strictColumnCounting = true;
bool   verbose = false;
bool   logPValue = false;
bool   trackPositions = false;
int    effectPrintPrecision = 4;
int    stderrPrintPrecision = 4;

String genomicControlFilter;
double genomicControlLambda = 0.0;
int    genomicControlCount = 0;

int    maxWarnings = 20;

StringArray filterLabel;
IntArray    filterColumn;
IntArray    filterCondition;
Vector      filterValue;
StringArray filterAlternate;
StringHash  filterSets;
IntArray    filterCounts;

bool studyOverlap = false;
double zCutoff = 1.0;
unsigned int studyOverlapStep = 1000;
Vector z1, z2;

void UserBreak(int)
   {
   printf("\n\n## METAL STOPPED BY USER\n\n");

   exit(EXIT_FAILURE);
   }

void OutOfMemory(int)
   {
   printf("\n\nMETAL HAS CRASHED\n\n"
          "The operating system has decided to terminate METAL.\n"
          "This could be due to a memory access problem, or it\n"
          "could be you encountered a bug.\n\n"
          "To help improve this program, please report bugs by\n"
          "e-mailing goncalo@umich.edu. Please include details\n"
          "of how the problem you encountered can be reproduced.\n");

   exit(EXIT_FAILURE);
   }

typedef void (*signal_handler)(int);

void SetupCrashHandlers()
   {
   signal(SIGINT, (signal_handler) UserBreak);
   signal(SIGSEGV, (signal_handler) OutOfMemory);
   signal(SIGABRT, (signal_handler) OutOfMemory);
   }

void PrintablePvalue(String & buffer, double statistic)
   {
   statistic = fabs(statistic);

   if (logPValue)
      {
      double log_pvalue = (logndist(fabs(statistic), true) + log(2.0)) / log(10.0);

      buffer.printf("%.2f", log_pvalue);
      }
   else if (statistic < 20)
      {
      double pvalue = 2.0 * ndist(fabs(statistic), true);

      buffer.printf("%.4g", pvalue);
      }
   else
      {
      double log_pvalue = (logndist(fabs(statistic), true) + log(2.0)) / log(10.0);

      double exponent = floor(log_pvalue);
      double real = pow(10, log_pvalue - exponent);

      if (real >= 9.995)
         {
         real /= 10;
         exponent++;
         }

      buffer.printf("%.2fe%.0f", real, exponent);
      }
   }

int CreateNewMarkerId(const String &markerName) {
    int customColumns = customLabels.Length();
    int marker = 0;

    if (markerLookup.Entries() == 0 && customColumns > 0)
        custom = new Vector[customColumns];

    markerLookup.SetInteger(markerName, marker = markerLookup.Entries());
    statistics.Push(0);
    if (averageFrequencies) frequencies.Push(0);
    if (averageFrequencies) frequencies2.Push(0);
    if (minMaxFrequencies) minFrequencies.Push(_NAN_);
    if (minMaxFrequencies) maxFrequencies.Push(_NAN_);
    weights.Push(0);
    samples.Push(0);
    allele1.Push("");
    allele2.Push("");
    chromosomes.Push("");
    positions.Push("");
    original_flipped += '?';

    for (int i = 0; i < customColumns; i++)
        custom[i].Push(0);

    if (genomicControlFilter.Length())
        genomicControlFilter += '.';

    return marker;
}

int GetMarkerId(const String & markerName)
   {
   int marker = markerLookup.Integer(markerName);

   if (marker >= 0)
      return marker;

   return CreateNewMarkerId(markerName);
   }

void ClearFilters()
   {
   printf("# Clearing user defined filters ...\n");
   filterLabel.Dimension(0);
   filterCondition.Dimension(0);
   filterValue.Dimension(0);
   filterAlternate.Dimension(0);
   filterCounts.Dimension(0);
   filterSets.Clear();
   }

void ClearAll() {
    printf("## Clearing all stored statistics ...\n");
    markerLookup.Clear();
    statistics.Dimension(0);
    weights.Dimension(0);
    samples.Dimension(0);
    allele1.Dimension(0);
    allele2.Dimension(0);
    chromosomes.Dimension(0);
    positions.Dimension(0);
    frequencies.Dimension(0);
    frequencies2.Dimension(0);
    minFrequencies.Dimension(0);
    maxFrequencies.Dimension(0);
    original_flipped.Clear();
    filenames.Clear();
    directions.Clear();
    customVariables.Clear();
    customLabels.Clear();

    if (custom != NULL) {
        delete[] custom;
        custom = NULL;
    }

    ClearFilters();

    if (processedFiles != NULL) {
        delete processedFiles;
        processedFiles = NULL;
    }
}

int TranslateCondition(String &condition) {
    if (condition == ">")
        return GREATER_THAN;
    else if (condition == ">=")
        return GREATER_THAN_OR_EQUAL;
    else if (condition == "<")
        return LESS_THAN;
    else if (condition == "<=")
        return LESS_THAN_OR_EQUAL;
    else if (condition == "==" || condition == "=" || condition == "IS")
        return EQUAL_TO;
    else if (condition == "!=" || condition == "<>" || condition == "ISNOT")
        return NOT_EQUAL_TO;
    else if (condition == "IN")
        return STRING_IN_SET;
    else
        return -1;
}

void AddFilter(StringArray & filter)
   {
   if (filter.Length() != 4)
      {
      printf("# Filter ignored: Command requires exactly 3 arguments separated by whitespace\n"
             "#    The arguments correspond to a column label, condition, and filter value\n");
      return;
      }

   int condition = TranslateCondition(filter[2]);

   if (condition == -1)
      {
      printf("# Filter ignored: the condition '%s' is not recognized\n", (const char *) filter[2]);
      return;
      }

   if (condition == STRING_IN_SET)
      {
      StringArray set;
      String string;

      if (filter[3][0] != '(' || filter[3].Last() != ')')
         {
         printf("# Filter ignored: brackets required around set of valid values, e.g. VALUE IN (ITEM1,ITEM2,...)\n");
         return;
         }

      string = filter[3];
      string.Trim('(');
      string.Trim(')');
      set.AddTokens(string,',');

      string.printf("%d#", filterLabel.Length());
      for (int i = 0; i < set.Length(); i++)
         filterSets.Add(string + set[i]);
      }

   bool isNumber = isdigit(filter[3][0]) || filter[3][0] == '-' || filter[3][0] == '+';

   if (!isNumber && condition < 0x80 /* STRING_ conditions are all >= 0x80 */ )
      {
      printf("# Filter ignored: the condition '%s' requires a numeric right-hand side value", (const char *) filter[2]);
      return;
      }

   if (!isNumber && condition == EQUAL_TO) condition = STRING_MATCH;
   if (!isNumber && condition == NOT_EQUAL_TO) condition = STRING_MISMATCH;

   filterLabel.Push(filter[1]);
   filterCondition.Push(condition);
   filterValue.Push(filter[3]);
   filterAlternate.Push(filter[3]);

   printf("# Added filter based on '%s' column\n", (const char *) filter[1]);
   }

void SetupFilters(StringArray & header)
   {
   filterColumn.Dimension(filterLabel.Length());
   filterCounts.Dimension(filterLabel.Length());
   filterCounts.Zero();

   for (int i = 0; i < filterLabel.Length(); i++)
      {
      filterColumn[i] = header.SlowFind(filterLabel[i]);

      if (filterColumn[i] < 0)
         printf("## WARNING: No '%s' column found -- filter based on this column will be ignored\n", (const char *) filterLabel[i]);
      }
   }

bool ApplyFilter(StringArray & row)
   {
   for (int i = 0; i < filterColumn.Length(); i++)
      if (filterColumn[i] >= 0)
         {
         String & text = row[filterColumn[i]];

         text.Trim();

         bool isValue = isdigit(text[0]) || text[0] == '.' && isdigit(text[1]) ||
            ( text[0] == '-' || text[0] == '+') &&
            ( isdigit(text[1]) || text[1] == '.' && isdigit(text[2]) );

         switch(filterCondition[i])
            {
            case LESS_THAN :
               if (!(row[filterColumn[i]].AsDouble() < filterValue[i]) || !isValue)
                  {
                  filterCounts[i]++;
                  return false;
                  }
               break;
            case LESS_THAN_OR_EQUAL :
               if (!(row[filterColumn[i]].AsDouble() <= filterValue[i]) || !isValue)
                  {
                  filterCounts[i]++;
                  return false;
                  }
               break;
            case GREATER_THAN :
               if (!(row[filterColumn[i]].AsDouble() > filterValue[i]) || !isValue)
                  {
                  filterCounts[i]++;
                  return false;
                  }
               break;
            case GREATER_THAN_OR_EQUAL :
               if (!(row[filterColumn[i]].AsDouble() >= filterValue[i]) || !isValue)
                  {
                  filterCounts[i]++;
                  return false;
                  }
               break;
            case EQUAL_TO :
               if (!(row[filterColumn[i]].AsDouble() == filterValue[i]) || !isValue)
                  {
                  filterCounts[i]++;
                  return false;
                  }
               break;
            case NOT_EQUAL_TO :
               if (!(row[filterColumn[i]].AsDouble() != filterValue[i]) || !isValue)
                  {
                  filterCounts[i]++;
                  return false;
                  }
               break;
            case STRING_MATCH :
               if (!(row[filterColumn[i]] == filterAlternate[i]))
                  {
                  filterCounts[i]++;
                  return false;
                  }
               break;
            case STRING_MISMATCH :
               if (!(row[filterColumn[i]] != filterAlternate[i]))
                  {
                  filterCounts[i]++;
                  return false;
                  }
               break;
            case STRING_IN_SET:
               {
               String key;
               key.printf("%d#%s", i, (const char *) row[filterColumn[i]]);

               if (filterSets.Find(key) < 0)
                  {
                  filterCounts[i]++;
                  return false;
                  }
               }
               break;
            }
         }

   return true;
   }

void FilterSummary()
   {
   if (filterColumn.Length() == 0) return;

   if (filterCounts.Sum() == 0)
      {
      printf("## Although a filter is defined, no lines were excluded because of it\n");
      return;
      }

   const char * conditions[7] = {"<", "<=", "==", ">=", ">", "!=", "IN"};

   if (filterColumn.Length() > 1);
      printf("# %8d lines filtered, for the following reasons:\n", filterCounts.Sum());

   for (int i = 0; i < filterColumn.Length(); i++)
      if (filterCounts[i])
         printf("# %8d lines filtered due to constraint %s %s %s\n",
                filterCounts[i],
                (const char *) filterLabel[i],
                (const char *) conditions[(filterCondition[i] & 0x7F) - 1],
                (const char *) filterAlternate[i]);
   }

void NumbersToLetters(String & al)
   {
   if (al == "1") al = "a"; else
   if (al == "2") al = "c"; else
   if (al == "3") al = "g"; else
   if (al == "4") al = "t";
   }

void FlipAllele(String & al)
   {
   if (al == "A" || al == "a") al = "t"; else
   if (al == "C" || al == "c") al = "g"; else
   if (al == "G" || al == "g") al = "c"; else
   if (al == "T" || al == "t") al = "a";
   }

bool FlipAlleles(String & al1, String & al2, double & effect, double & freq)
   {
   al1.ToLower();
   al2.ToLower();

   if (al1 > al2)
      {
      al1.Swap(al2);
      effect *= -1.0;
      freq = 1.0 - freq;
      }

   if (al1 == "a" || al1 == "c" && al2 == "g")
      return false;

   FlipAllele(al1);
   FlipAllele(al2);

   if (al1 > al2)
      {
      al1.Swap(al2);
      effect *= -1.0;
      freq = 1.0 - freq;
      }

   return true;
   }

bool GuessSecondAllele(int marker, String & al1, String & al2)
   {
   al1.ToLower();

   FlipAllele(al2 = al1);

   if (allele1[marker] == al1 && allele2[marker] == al2 ||
       allele1[marker] == al2 && allele2[marker] == al1)
      return false;

   if (allele1[marker] == al1)
      {
      al2 = allele2[marker];
      return true;
      }
   else if (allele2[marker] == al1)
      {
      al2 = allele1[marker];
      return true;
      }
   else if (allele1[marker] == al2)
      {
      FlipAllele(al2 = allele2[marker]);
      return true;
      }
   else if (allele2[marker] == al2)
      {
      FlipAllele(al2 = allele1[marker]);
      return true;
      }

   return false;
   }

int    outfileround = 1;
String outfile("METAANALYSIS%d.TBL");

void Analyze(bool heterogeneity) {
    if (heterogeneity) {
        printf("###########################################################################\n");
        printf("## Running second pass analysis to evaluate heterogeneity...\n");

        hetDegreesOfFreedom.Dimension(statistics.Length());
        hetDegreesOfFreedom.Zero();

        hetStatistic.Dimension(statistics.Length());
        hetStatistic.Zero();

        FileSummary *pointer = processedFiles;

        while (pointer != NULL) {
            if (!ReProcessFile(pointer)) {
                printf("## WARNING: Heterogeneity analysis failed...\n");
                heterogeneity = false;
                break;
            }

            pointer = pointer->next;
        }
        printf("\n");
    }

    String filename;
    filename.printf(outfile, outfileround++);

    printf("###########################################################################\n");
    printf("## Executing meta-analysis ...\n");
    printf("## Complete results will be stored in file '%s'\n", (const char *) filename);
    printf("## Column descriptions will be stored in file '%s.info'\n", (const char *) filename);

    FILE *f = fopen(filename, "wt");

    if (f == NULL) {
        printf("ERROR: Failed to open output file '%s' ... \n",
               (const char *) filename);
        return;
    }

    fprintf(f, "%sMarkerName\tAllele1\tAllele2\t%s%s%s\t%s%s\t%s\tDirection%s%s",
            trackPositions ? "Chromosome\tPosition\t" : "",
            averageFrequencies ? "Freq1\tFreqSE\t" : "",
            minMaxFrequencies ? "MinFreq\tMaxFreq\t" : "",
            useStandardErrors ? "Effect" : "Weight",
            useStandardErrors ? "StdErr" : "Zscore",
            studyOverlap ? "\tN" : "",
            logPValue ? "log(P)" : "P-value",
            heterogeneity ? "\tHetISq\tHetChiSq\tHetDf\t" : "",
            heterogeneity ? (logPValue ? "logHetP" : "HetPVal") : "");

    for (int i = 0; i < customVariables.Length(); i++)
        fprintf(f, "\t%s", (const char *) customVariables[i]);

    fprintf(f, "\n");

    double largeStatistic = 0.0;
    String smallp("1.0");
    String smallm("nothing tested");

    String al1, al2;

    String direction;

    int count = 0, skipped = 0;
    for (int i = 0; i < markerLookup.Capacity(); i++)
        if (markerLookup.SlotInUse(i)) {
            int marker = markerLookup.Integer(i);

            if (weights[marker] == 0.0) continue;

            if (!useStandardErrors && weights[marker] < minweight) {
                skipped++;
                continue;
            }

            double statistic = statistics[marker] / sqrt(weights[marker]);
            double frequency = averageFrequencies ? frequencies[marker] / weights[marker] : 0.5;
            double frequency2 = averageFrequencies ? frequencies2[marker] / weights[marker] : 0.5;
            double minFrequency = minMaxFrequencies ? minFrequencies[marker] : 0.5;
            double maxFrequency = minMaxFrequencies ? maxFrequencies[marker] : 0.5;


            // double pvalue = 2.0 * ndist(fabs(statistic), true);
            String pvalue;
            PrintablePvalue(pvalue, statistic);

            frequency2 = frequency2 - frequency * frequency;
            frequency2 = frequency2 > 0.0 ? sqrt(frequency2) : 0.0;

            al1 = allele1[marker];
            al2 = allele2[marker];

            if (original_flipped[marker] == 'Y')
                FlipAllele(al1), FlipAllele(al2);

            direction.Clear();
            for (int j = 0; j < filenames.Length(); j++)
                direction += marker < directions[j].Length() ? directions[j][marker] : '?';

            if (trackPositions) {
                fprintf(f, "%s\t%s\t", (const char*) chromosomes[marker], (const char*) positions[marker]);
            }

            fprintf(f, "%s\t%s\t%s\t",
                    (const char *) markerLookup[i],
                    (const char *) al1,
                    (const char *) al2);

            if (averageFrequencies)
                fprintf(f, "%.4f\t%.4f\t", frequency, frequency2);

            if (minMaxFrequencies)
                fprintf(f, "%.4f\t%.4f\t", minFrequency, maxFrequency);

            fprintf(f, "%.*f\t%.*f\t",
                    useStandardErrors ? effectPrintPrecision : 2,
                    useStandardErrors ? statistics[marker] / weights[marker] : weights[marker],
                    useStandardErrors ? stderrPrintPrecision : 3,
                    useStandardErrors ? sqrt(1.0 / weights[marker]) : statistic);

            if (studyOverlap) {
                fprintf(f, "%.2f\t", samples[marker]);
            }

            fprintf(f, "%s\t%s",
                    (const char *) pvalue,
                    (const char *) direction);
            count++;

            if (heterogeneity) {
                double p =
                        (hetStatistic[marker] < 1e-7 || hetDegreesOfFreedom[marker] <= 1) ?
                        1.0 : chidist(hetStatistic[marker], hetDegreesOfFreedom[marker] - 1);
                double I2 =
                        (hetStatistic[marker] <= hetDegreesOfFreedom[marker] - 1) || hetDegreesOfFreedom[marker] <= 1 ?
                        0.0 : (hetStatistic[marker] - hetDegreesOfFreedom[marker] + 1) / hetStatistic[marker] * 100.;

                if (logPValue) p = (p < 1.0) ? log(p) / log(10.0) : 0.0;

                fprintf(f, "\t%.1f\t%.3f\t%d\t%.4g",
                        I2, hetStatistic[marker], hetDegreesOfFreedom[marker] - 1, p);
            }

            for (int j = 0; j < customVariables.Length(); j++)
                fprintf(f, "\t%g", custom[j][marker]);

            fprintf(f, "\n");

            if (fabs(statistic) > largeStatistic) {
                largeStatistic = fabs(statistic);
                smallp = pvalue;
                smallm = markerLookup[i];
            }
        }

    printf("## Completed meta-analysis for %d markers!\n", count);

    if (largeStatistic > 0.0)
        printf("## Smallest p-value is %s at marker '%s'\n", (const char *) smallp, (const char *) smallm);

    if (skipped)
        printf("## A total of %d markers with less than MINWEIGHT = %g were skipped\n", skipped, minweight);

    printf("\n");

    fclose(f);

    f = fopen(filename + ".info", "wt");

    if (f == NULL) {
        printf("ERROR: Failed to open output file '%s.info' ... \n",
               (const char *) filename);
        return;
    }

    fprintf(f, "# This file contains a short description of the columns in the\n"
                    "# meta-analysis summary file, named '%s'\n\n"
                    "%s"
                    "# Marker      - this is the marker name\n"
                    "# Allele1     - the first allele for this marker in the first file where it occurs\n"
                    "# Allele2     - the second allele for this marker in the first file where it occurs\n"
                    "%s"
                    "%s"
                    "%s"
                    "%s"
                    "# %s meta-analysis p-value\n"
                    "# Direction   - summary of effect direction for each study, with one '+' or '-' per study\n"
                    "%s%s",
            (const char *) filename,
            !trackPositions ? "" :
            "# Chromosome  - chromosome name\n"
            "# Position    - chromosomal position\n" ,
            !averageFrequencies ? "" :
            "# Freq1       - weighted average of frequency for allele 1 across all studies\n"
            "# FreqSE      - corresponding standard error for allele frequency estimate\n",
            !minMaxFrequencies ? "" :
            "# MinFreq     - minimum frequency for allele 1 across all studies\n"
            "# MaxFreq     - maximum frequency for allele 1 across all studies\n",
            useStandardErrors ?
            "# Effect      - overall estimated effect size for allele1\n"
            "# StdErr      - overall standard error for effect size estimate\n" :
            "# Weight      - the sum of the individual study weights (typically, N) for this marker\n"
            "# Z-score     - the combined z-statistic for this marker\n",
            !studyOverlap ? "" :
            "# N           - sample size corrected for overlap between studies\n",
            logPValue ? "log(P)      - log of" : "P-value     -",
            !heterogeneity ? "" :
            "# HetISq      - I^2 statistic which measures heterogeneity on scale of 0-100%\n"
            "# HetChiSq    - chi-squared statistic in simple test of heterogeneity\n"
            "# df          - degrees of freedom for heterogeneity statistic\n",
            !heterogeneity ? "" :
            logPValue ?
            "# logHetP     - log of p-value for heterogeneity statistic\n" :
            "# HetPVal     - P-value for heterogeneity statistic\n");

    for (int i = 0; i < customVariables.Length(); i++)
        fprintf(f, "# %-9s - custom variable %d\n", (const char *) customVariables[i], i + 1);

    fprintf(f, "\n# Input for this meta-analysis was stored in the files:\n");
    for (int i = 0; i < filenames.Length(); i++)
        fprintf(f, "# --> Input File %d : %s\n", i + 1, (const char *) filenames[i]);
    fprintf(f, "\n");

    fclose(f);
}

bool flip = false;


double trunc_norm(double rho) {
    double sum_d = 0.0;
    for (int i = 0; i < z1.Length(); ++i) {
        sum_d += log(binormp(z1[i], z2[i], rho));
    }
    double l = -sum_d + z1.Length() * log(binormq(-1.0 * zCutoff, -1.0 * zCutoff, rho) - binormq(-1.0 * zCutoff, zCutoff, rho) - binormq(zCutoff, -1.0 * zCutoff, rho) + binormq(zCutoff, zCutoff, rho));
    return l;
}


void ProcessFile(String & filename, FileSummary * history) {
    IFILE f = ifopen(history->filename = filename, "rb");
    history->processedMarkers = 0;

    if (f == NULL) {
        printf("## Failed to open file '%s'\n", (const char *) filename);
        return;
    }

    printf("###########################################################################\n");
    printf("## Processing file '%s'\n", (const char *) filename);

    String input;
    StringArray tokens;

    tokens.ReplaceTokens(input.ReadLine(f), separators);
    history->header = input;

    if (tokens.Length() == 0) {
        printf("## ERROR: Header line is empty\n\n");
        ifclose(f);
        return;
    }

    int markerColumn = history->markerColumn = tokens.SlowFind(markerLabel);
    if (markerColumn < 0) {
        printf("## ERROR: No '%s' column found\n\n", (const char *) markerLabel);
        ifclose(f);
        return;
    }

    int pvalueColumn = history->pvalueColumn = tokens.SlowFind(pvalueLabel);
    if (pvalueColumn < 0 && !useStandardErrors) {
        printf("## ERROR: No '%s' column found\n\n", (const char *) pvalueLabel);
        ifclose(f);
        return;
    }

    int effectColumn = history->effectColumn = tokens.SlowFind(effectLabel);
    bool logTransform = history->logTransform = false;
    if (effectColumn < 0 && effectLabel.Left(4) == "log(" && effectLabel.Last() == ')') {
        logTransform = history->logTransform = true;
        effectColumn = history->effectColumn = tokens.SlowFind(effectLabel.Mid(4, effectLabel.Length() - 2).Trim());
    }

    if (effectColumn < 0) {
        if (!useStandardErrors) {
            printf("## WARNING: No '%s' column found -- assuming all effects are positive!\n", (const char *) effectLabel);
        } else {
            printf("## ERROR: No '%s' column found\n\n", (const char *) effectLabel);
            ifclose(f);
            return;
        }
    }

    int weightColumn = history->weightColumn = tokens.SlowFind(weightLabel);
    if (weightColumn < 0 && !useStandardErrors) {
        printf("## WARNING: No '%s' column found -- using DEFAULTWEIGHT = %g\n", (const char *) weightLabel, weight);
    }

    int firstColumn = history->firstColumn = tokens.SlowFind(firstAllele);
    int secondColumn = history->secondColumn = tokens.SlowFind(secondAllele);
    bool guessAlleles = false;

    if (firstColumn >= 0 && secondColumn < 0) {
        printf("## WARNING: Column '%s' not found, will try to guess second allele\n", (const char *) secondAllele);
        guessAlleles = true;
    }

    if (firstColumn < 0) {
        printf("## WARNING: Column '%s', '%s' or both not found -- assuming all effects refer to same allele\n", (const char *) firstAllele, (const char *) secondAllele, firstColumn = 0);
    }

    int stderrColumn = history->stderrColumn = tokens.SlowFind(stderrLabel);
    if (stderrColumn < 0 && useStandardErrors) {
        printf("## ERROR: Analysis based on standard errors requested but no '%s' column found\n\n", (const char *) stderrLabel);
        ifclose(f);
        return;
    }

    bool useFrequencies = minMaxFrequencies || averageFrequencies;
    int freqColumn = history->freqColumn = tokens.SlowFind(frequencyLabel);
    if (freqColumn < 0 && useFrequencies) {
        printf("## ERROR: Averaging of allele frequencies requested, but no '%s' column found\n\n", (const char *) frequencyLabel);
        ifclose(f);
        return;
    }

    int strandColumn = history->strandColumn = tokens.SlowFind(strandLabel);
    if (strandColumn < 0 && useStrand) {
        printf("## ERROR: Strand column labeleled '%s' not found\n\n", (const char *) strandLabel);
        ifclose(f);
        return;
    }

    int chromosomeColumn = history->chromosomeColumn = tokens.SlowFind(chromosomeLabel);
    if (chromosomeColumn < 0 && trackPositions) {
        printf("## ERROR: Chromosome column labeled '%s' not found\n\n", (const char *) strandLabel);
        ifclose(f);
        return;
    }

    int positionColumn = history->positionColumn = tokens.SlowFind(positionLabel);
    if (positionColumn < 0 && trackPositions) {
        printf("## ERROR: Position column labeled '%s' not found\n\n", (const char *) positionLabel);
        ifclose(f);
        return;
    }

    SetupFilters(tokens);
    history->filterLabel = filterLabel;
    history->filterColumn = filterColumn;
    history->filterCondition = filterCondition;
    history->filterValue = filterValue;
    history->filterAlternate = filterAlternate;
    history->filterSets = filterSets;
    history->filterCounts = filterCounts;

    int minColumns = markerColumn + 1;
    if (weightColumn >= minColumns && !useStandardErrors) minColumns = weightColumn + 1;
    if (effectColumn >= minColumns) minColumns = effectColumn + 1;
    if (pvalueColumn >= minColumns && !useStandardErrors) minColumns = pvalueColumn + 1;
    if (firstColumn >= minColumns) minColumns = firstColumn + 1;
    if (secondColumn >= minColumns) minColumns = secondColumn + 1;
    if (stderrColumn >= minColumns && useStandardErrors) minColumns = stderrColumn + 1;
    if (strandColumn >= minColumns && useStrand) minColumns = strandColumn + 1;
    if (freqColumn >= minColumns && useFrequencies) minColumns = freqColumn + 1;
    if (chromosomeColumn >= minColumns && trackPositions) minColumns = chromosomeColumn + 1;
    if (positionColumn >= minColumns && trackPositions) minColumns = positionColumn + 1;
    if (filterColumn.Max() >= minColumns) minColumns = filterColumn.Max() + 1;

    IntArray customColumns(customLabels.Length());
    for (int i = 0; i < customLabels.Length(); i++) {
        customColumns[i] = tokens.SlowFind(customLabels[i]);
        if (customColumns[i] < 0) {
            printf("## ERROR: Required column '%s' not found\n\n", (const char *) customLabels[i]);
            ifclose(f);
            return;
        }
    }
    if (customColumns.Max() > minColumns)
        minColumns = customColumns.Max();

    history->minColumns = minColumns;

    int expectedColumns = tokens.Length();
    history->expectedColumns = expectedColumns;
    history->strictColumnCounting = strictColumnCounting;

    Vector backupStatistics, backupWeights, backupSamples, backupFrequencies, backupFrequencies2;
    Vector chiSquareds;

    if (genomicControl) {
        backupStatistics = statistics;
        statistics.Zero();
    }

    if (genomicControl && useStandardErrors) {
        backupFrequencies2 = frequencies2;
        backupFrequencies = frequencies;
        backupWeights = weights;

        frequencies2.Zero();
        frequencies.Zero();
        weights.Zero();
    }

    // DT
    if (studyOverlap) {
        if (!genomicControl) {
            backupStatistics = statistics;
            statistics.Zero();
        }
        if (!useStandardErrors) { // currently always true for this method
            backupWeights = weights;
            backupSamples = samples;
            weights.Zero();
            samples.Zero();
        }
    }

    String direction;
    direction.Fill('?', allele1.Length());

    int invalid = 0;
    int invalidEffect = 0;
    int badAlleles = 0;
    int badGuesses = 0;
    int duplicates = 0;
    int blindGuesses = 0;
    int wrongColumnCount = 0;
    int shortColumnCount = 0;
    int badChromosome = 0;
    int badPosition = 0;

    history->weight = weight;

    filenames.Push(filename);

    if (verbose) {
        printf("# MARKER\tAL1\tAL2\t");

        if (!useStandardErrors)
            printf("N\tZ\tPVAL\t");
        else
            printf("EFF\tSTDERR\tPVAL\t");

        if (minMaxFrequencies || averageFrequencies)
            printf("FREQ\t");

        printf("SOURCE\n");
    }

    history->separators = separators;
    history->useStrand = useStrand;
    while (!ifeof(f)) {
        if (separators.Length() != 1)
            tokens.ReplaceTokens(input.ReadLine(f).Trim(), separators);
        else
            tokens.ReplaceColumns(input.ReadLine(f), separators[0]);

        if (input[0] == '#') continue;
        if (tokens.Length() == 0) continue;

        if (tokens.Length() != expectedColumns) {
            wrongColumnCount++;
            if (strictColumnCounting) continue;
        }

        if (tokens.Length() < minColumns) {
            shortColumnCount++;
            continue;
        }

        if (!ApplyFilter(tokens)) continue;

        int marker = markerLookup.Integer(tokens[markerColumn]);

        if (marker < 0 && guessAlleles) {
            if (++blindGuesses > maxWarnings)
                continue;
            printf("## WARNING: Can't guess alleles for '%s', marker not seen in previous files\n", (const char *) tokens[markerColumn]);
            continue;
        }

        if (marker < 0) {
            marker = CreateNewMarkerId(tokens[markerColumn]);
            direction += '?';
        }

        double w, z;

        if (!useStandardErrors) {
            long double p = tokens[pvalueColumn].AsLongDouble();

            if (p <= 0 || p > 1.0) {
                if (++invalid > maxWarnings) continue;
                printf("## WARNING: Invalid pvalue for marker %s, ignored\n", (const char *) tokens[markerColumn]);
                continue;
            }

            z = -ninv(p * 0.5);
            w = weightColumn >= 0 ? tokens[weightColumn].AsDouble() : weight;

            if (!logTransform && effectColumn >= 0 && tokens[effectColumn][0] == '-')
                z *= -1;

            if (logTransform) {
                double eff = tokens[effectColumn].AsDouble();

                if (eff <= 0.0) {
                    if (++invalidEffect > maxWarnings) continue;
                    printf("## WARNING: Invalid effect %s for marker %s, ignored\n",
                           (const char *) effectLabel, (const char *) tokens[markerColumn]);
                    continue;
                }

                if (eff < 1.0)
                    z *= -1;
            }
        } else {
            double eff = tokens[effectColumn].AsDouble();
            double sd = tokens[stderrColumn].AsDouble();

            if (logTransform) {
#if  defined(BSD_SOURCE) || defined(_SVID_SOURCE) || defined(_ISOC99_SOURCE) || XOPEN_SOURCE > 600
                if (eff <= 0.0 || isnan(eff) || isinf(eff))
#else
                if (eff <= 0.0)
#endif
                {
                    if (++invalidEffect > maxWarnings) continue;
                    printf("## WARNING: Invalid log(effect) for marker %s, ignored\n",
                           (const char *) tokens[markerColumn]);
                    continue;
                }

                eff = log(eff);
            }

            if (sd <= 0) {
                if (++invalid > maxWarnings) continue;
                printf("## WARNING: Invalid standard error for marker %s, ignored\n",
                       (const char *) tokens[markerColumn]);
                continue;
            }

            z = eff;
            w = 1.0 / (sd * sd);
        }

        double freq = 0.0;
        if (useFrequencies)
            freq = tokens[freqColumn].AsDouble();

        if (flip)
            z *= -1;

        if (firstColumn >= 0 && secondColumn >= 0) {
            NumbersToLetters(tokens[firstColumn]);
            NumbersToLetters(tokens[secondColumn]);

            if (useStrand && tokens[strandColumn] == "-") {
                FlipAllele(tokens[firstColumn]);
                FlipAllele(tokens[secondColumn]);
            }

            bool flip = FlipAlleles(tokens[firstColumn], tokens[secondColumn], z, freq);

            if (allele1[marker] == "")
                allele1[marker] = tokens[firstColumn],
                allele2[marker] = tokens[secondColumn],
                original_flipped[marker] = flip ? 'Y' : 'N';
            else if (tokens[firstColumn] != allele1[marker] ||
                     tokens[secondColumn] != allele2[marker]) {
                if (++badAlleles > maxWarnings) continue;
                printf("## WARNING: Bad alleles for marker '%s', expecting '%s/%s' found '%s/%s'\n",
                       (const char *) tokens[markerColumn],
                       (const char *) allele1[marker],
                       (const char *) allele2[marker],
                       (const char *) tokens[firstColumn],
                       (const char *) tokens[secondColumn]);
                continue;
            }
        } else if (guessAlleles) {
            tokens.Push("");
            NumbersToLetters(tokens[firstColumn]);
            if (!GuessSecondAllele(marker, tokens[firstColumn], tokens.Last())) {
                if (++badGuesses <= maxWarnings)
                    printf("## WARNING: Bad allele or ambiguous strand for marker '%s', expecting '%s/%s' found '%s'\n",
                           (const char *) tokens[markerColumn],
                           (const char *) allele1[marker],
                           (const char *) allele2[marker],
                           (const char *) tokens[firstColumn]);
                continue;
            }
            FlipAlleles(tokens[firstColumn], tokens.Last(), z, freq);
        }

        if (direction[marker] != '?') {
            if (++duplicates > maxWarnings) continue;
            printf("## WARNING: Marker '%s' duplicated in input, first instance used, others skipped\n", (const char *) tokens[markerColumn]);
            continue;
        }

        if (trackPositions) {
            if (chromosomes[marker] == "") {
                chromosomes[marker] = tokens[chromosomeColumn];
            } else if (tokens[chromosomeColumn] != chromosomes[marker]) {
                if (++badChromosome > maxWarnings) {
                    continue;
                }
                printf("## WARNING: Bad chromosome for marker '%s', expecting '%s' found '%s'\n", (const char *) tokens[markerColumn], (const char *) chromosomes[marker], (const char *) tokens[chromosomeColumn]);
                continue;
            }
            if (positions[marker] == "") {
                positions[marker] = tokens[positionColumn];
            } else if (tokens[positionColumn] != positions[marker]) {
                if (++badPosition > maxWarnings) {
                    continue;
                }
                printf("## WARNING: Bad position for marker '%s', expecting '%s' found '%s'\n", (const char *) tokens[markerColumn], (const char *) positions[marker], (const char *) tokens[positionColumn]);
                continue;
            }
        }

        direction[marker] = z == 0.0 ? '0' : (z > 0.0 ? '+' : '-');

        if (verbose) {
            String al1 = allele1[marker];
            String al2 = allele2[marker];

            if (original_flipped[marker] == 'Y')
                FlipAllele(al1), FlipAllele(al2);

            printf("# %s\t%s\t%s\t", (const char *) tokens[markerColumn],
                   (const char *) al1,
                   (const char *) al2);

            if (!useStandardErrors)
                printf("%g\t%+.22f\t%.4g\t", w, z, ndist(fabs(z), true) * 2.0);
            else
                printf("%.3f\t%.22f\t%.4g\t", z, w > 0.0 ? sqrt(1.0 / w) : 0.0, chidist(z * z * w, 1));

            if (minMaxFrequencies || averageFrequencies)
                printf("%.3f\t", freq);

            printf("%s\n", (const char *) filenames.Last());
        }

        if (studyOverlap) {
            samples[marker] += w;
        }

        if (!useStandardErrors) {
            statistics[marker] += sqrt(w) * z;
            weights[marker] += w;
            if (genomicControl)
                if (genomicControlFilter.Length() == 0 || genomicControlFilter[marker] == 'Y')
                    chiSquareds.Push(z * z);
        } else {
            statistics[marker] += w * z;
            weights[marker] += w;
            if (genomicControl)
                if (genomicControlFilter.Length() == 0 || genomicControlFilter[marker] == 'Y')
                    chiSquareds.Push(z * z * w);
        }

        if (averageFrequencies) {
            frequencies[marker] += w * freq;
            frequencies2[marker] += w * freq * freq;
        }

        if (minMaxFrequencies) {
            if (minFrequencies[marker] > freq || minFrequencies[marker] == _NAN_)
                minFrequencies[marker] = freq;
            if (maxFrequencies[marker] < freq || maxFrequencies[marker] == _NAN_)
                maxFrequencies[marker] = freq;
        }

        for (int i = 0; i < customColumns.Length(); i++)
            custom[i][marker] += tokens[customColumns[i]].AsDouble();

        history->processedMarkers++;
    }

    if (invalid > maxWarnings)
        printf("## WARNING: Invalid %s for %d other markers also ignored\n",
               useStandardErrors ? "standard errors" : "p-values", invalid - maxWarnings);

    if (invalidEffect > maxWarnings)
        printf("## WARNING: Invalid log(effect) for %d other markers also ignored\n", invalidEffect - maxWarnings);

    if (badGuesses > maxWarnings)
        printf("## WARNING: Failed to guess second allele for %d other markers\n", badGuesses - maxWarnings);

    if (duplicates > maxWarnings)
        printf("## WARNING: An additional %d rows with duplicate marker names were ignored\n", duplicates - maxWarnings);

    if (badAlleles > maxWarnings)
        printf("## WARNING: Allele names don't match previous occurences at %d additional markers\n", badAlleles - maxWarnings);

    if (blindGuesses > maxWarnings)
        printf("## WARNING: An additional %d markers whose alleles couldn't be guessed were skipped\n", blindGuesses - maxWarnings);

    if (wrongColumnCount) {
        if (strictColumnCounting)
            printf("## WARNING: %d input lines were skipped because they did not include exactly %d entries\n"
                           "##          For lenient handling of trailing tokens, use the COLUMNCOUNTING LENIENT command\n",
                   wrongColumnCount, expectedColumns);
        else {
            printf("## WARNING: %d input lines did not include exactly %d columns as in the header line\n", wrongColumnCount, expectedColumns);
            printf("## WARNING: %d input lines included less than the %d required columns and were skipped\n", shortColumnCount, minColumns);
        }
    }

    if (trackPositions) {
        if (badChromosome > maxWarnings) {
            printf("## WARNING: An additional %d markers whose chromosomes didn't match\n", badChromosome - maxWarnings);
        }
        if (badPosition > maxWarnings) {
            printf("## WARNING: An additional %d markers whose positions didn't match\n", badPosition - maxWarnings);
        }
    }

    FilterSummary();

    directions.Push(direction);

    history->genomicControl = 1.0;
    if (genomicControl) {
        if (chiSquareds.Length() == 0 && genomicControlLambda == 0.0 && genomicControlFilter.Length() == 0) {
            printf("## WARNING: Genomic control parameter cannot be calculated, no valid input\n");

            statistics.Swap(backupStatistics);

            if (useStandardErrors) {
                backupFrequencies2.Swap(frequencies2);
                backupFrequencies.Swap(frequencies);
                backupWeights.Swap(weights);
            } else {
                if (studyOverlap) {
                    backupSamples.Swap(samples);
                    backupWeights.Swap(weights);
                }
            }
        } else {
            bool warning = false;
            double gc = genomicControlLambda;

            if (gc == 0.0)
                if (chiSquareds.Length()) {
                    chiSquareds.Sort();
                    gc = chiSquareds[0.5] / 0.4549364;
                } else {
                    warning = true;
                    printf("## WARNING: Genomic control parameter cannot be calculated, no null markers\n");
                    gc = 1.0;
                }

            if (gc <= 1.0) {
                if (!warning)
                    printf("## Genomic control parameter is %.3f, no adjustment made\n", gc);
                gc = 1.0;
            } else
                printf("## Genomic control parameter is %.3f, adjusting test statistics\n", gc);

            history->genomicControl = gc;

            if (!useStandardErrors) {
                statistics *= 1.0 / sqrt(gc);
            } else {
                statistics *= 1.0 / gc;
                weights *= 1.0 / gc;
                if (averageFrequencies) frequencies *= 1.0 / gc;
                if (averageFrequencies) frequencies2 *= 1.0 / gc;
            }
        }
    }

    std::map<const int, double> rho;
    std::map<const int, double>::iterator rho_it;

    if (studyOverlap & (backupStatistics.Length() > 0)) {
        double z1_value = 0.0, z2_value = 0.0;

        for (int i = 0; i < backupStatistics.Length(); i++) {
            if (direction[i] == '?') {
                continue;
            }
            rho.insert(std::pair<const int, double>((((int) (backupSamples[i] + samples[i])) / studyOverlapStep), 0.0));
        }

        for (rho_it = rho.begin(); rho_it != rho.end(); ++rho_it) {
            z1.Clear();
            z2.Clear();
            for (int i = 0; i < backupStatistics.Length(); ++i) {
                if (direction[i] == '?') {
                    continue;
                }
                if (((int) (backupSamples[i] + samples[i])) / studyOverlapStep != rho_it->first) {
                    continue;
                }
                z1_value = backupStatistics[i] / sqrt(backupWeights[i]);
                z2_value = statistics[i] / sqrt(weights[i]);
                if ((fabs(z1_value) < zCutoff) && (fabs(z2_value) < zCutoff)) {
                    z1.Push(z1_value);
                    z2.Push(z2_value);
                }
            }

            ScalarMinimizer s;
            s.func = trunc_norm;
            s.a = -1.0;
            s.b = 1.0;
            s.c = s.b - (s.b - s.a) / (GOLD + 1.0);
            s.fa = s.f(s.a);
            s.fb = s.f(s.b);
            s.fc = s.f(s.c);
            s.Brent();
            if (s.min < 0.0) {
                rho_it->second = 0.0;
            } else {
                rho_it->second = s.min;
            }
            printf("## Samples overlap for %d marker(s) with %d<=WEIGHT<%d: %.10f\n", z1.Length(), rho_it->first * studyOverlapStep,
                   rho_it->first * studyOverlapStep + studyOverlapStep, rho_it->second);
        }
    }

    if (studyOverlap) {
        double n12 = 0.0;
        for (int i = 0; i < backupStatistics.Length(); ++i) {
            n12 = rho[((int) (samples[i] + backupSamples[i])) / studyOverlapStep] * sqrt(weights[i] * backupWeights[i]);
            samples[i] += backupSamples[i] - n12;
            weights[i] += backupWeights[i] + 2 * n12;
            statistics[i] += backupStatistics[i];
        }
    } else if (genomicControl) {
        for (int i = 0; i < backupStatistics.Length(); i++)
            statistics[i] += backupStatistics[i];

        if (useStandardErrors) {
            for (int i = 0; i < backupWeights.Length(); i++)
                weights[i] += backupWeights[i];

            for (int i = 0; i < backupFrequencies.Length(); i++) {
                frequencies2[i] += backupFrequencies2[i];
                frequencies[i] += backupFrequencies[i];
            }
        }
    }

//    history->genomicControl = 1.0;
//    if (genomicControl) {
//        if (chiSquareds.Length() == 0 && genomicControlLambda == 0.0 && genomicControlFilter.Length() == 0) {
//            printf("## WARNING: Genomic control parameter cannot be calculated, no valid input\n");
//
//            statistics.Swap(backupStatistics);
//
//            if (useStandardErrors) {
//                backupFrequencies2.Swap(frequencies2);
//                backupFrequencies.Swap(frequencies);
//                backupWeights.Swap(weights);
//            }
//        } else {
//            bool warning = false;
//            double gc = genomicControlLambda;
//
//            if (gc == 0.0)
//                if (chiSquareds.Length()) {
//                    chiSquareds.Sort();
//                    gc = chiSquareds[0.5] / 0.4549364;
//                } else {
//                    warning = true;
//                    printf("## WARNING: Genomic control parameter cannot be calculated, no null markers\n");
//                    gc = 1.0;
//                }
//
//            if (gc <= 1.0) {
//                if (!warning)
//                    printf("## Genomic control parameter is %.3f, no adjustment made\n", gc);
//                gc = 1.0;
//            } else
//                printf("## Genomic control parameter is %.3f, adjusting test statistics\n", gc);
//
//            history->genomicControl = gc;
//
//            if (!useStandardErrors)
//                statistics *= 1.0 / sqrt(gc);
//            else
//                statistics *= 1.0 / gc;
//
//            for (int i = 0; i < backupStatistics.Length(); i++)
//                statistics[i] += backupStatistics[i];
//
//            if (useStandardErrors) {
//                weights *= 1.0 / gc;
//
//                if (averageFrequencies) frequencies *= 1.0 / gc;
//                if (averageFrequencies) frequencies2 *= 1.0 / gc;
//
//                for (int i = 0; i < backupWeights.Length(); i++)
//                    weights[i] += backupWeights[i];
//
//                for (int i = 0; i < backupFrequencies.Length(); i++) {
//                    frequencies2[i] += backupFrequencies2[i];
//                    frequencies[i] += backupFrequencies[i];
//                }
//            }
//        }
//    }

    printf("## Processed %d markers ...\n\n", history->processedMarkers);

    ifclose(f);
}





bool ReProcessFile(FileSummary * history) {
    if (history->processedMarkers == 0)
        return true;

    IFILE f = ifopen(history->filename, "rb");

    if (f == NULL) {
        printf("## Failed to open file '%s'\n", (const char *) history->filename);
        return false;
    }

    printf("## Processing file '%s'\n", (const char *) history->filename);

    String input;
    StringArray tokens;

    input.ReadLine(f);

    if (input != history->header) {
        printf("## ERROR: Input file has changed since analysis started\n\n");
        ifclose(f);
        return false;
    }

    int markerColumn = history->markerColumn;
    int pvalueColumn = history->pvalueColumn;
    int effectColumn = history->effectColumn;
    int weightColumn = history->weightColumn;
    int firstColumn = history->firstColumn;
    int secondColumn = history->secondColumn;
    int stderrColumn = history->stderrColumn;
    int freqColumn = history->freqColumn;
    int strandColumn = history->strandColumn;
    int expectedColumns = history->expectedColumns;

    bool strictColumnCounting = history->strictColumnCounting;
    bool useStrand = history->useStrand;

    if (firstColumn >= 0 && secondColumn < 0) {
        printf("## ERROR: Heterogeneity analysis requires both allele labels\n\n");
        ifclose(f);
        return false;
    }

    bool useFrequencies = minMaxFrequencies || averageFrequencies;

    history->filterLabel.Swap(filterLabel);
    history->filterColumn.Swap(filterColumn);
    history->filterCondition.Swap(filterCondition);
    history->filterValue.Swap(filterValue);
    history->filterAlternate.Swap(filterAlternate);
    history->filterSets.Swap(filterSets);
    history->filterCounts.Swap(filterCounts);

    int minColumns = history->minColumns;
    int processedMarkers = 0;

    String direction;
    direction.Fill('?', allele1.Length());

    while (!ifeof(f)) {
        if (history->separators.Length() != 1)
            tokens.ReplaceTokens(input.ReadLine(f).Trim(), history->separators);
        else
            tokens.ReplaceColumns(input.ReadLine(f), history->separators[0]);

        if (input[0] == '#') continue;

        if (tokens.Length() != expectedColumns)
            if (strictColumnCounting)
                continue;

        if (tokens.Length() < minColumns)
            continue;

        if (!ApplyFilter(tokens)) continue;

        int marker = markerLookup.Integer(tokens[markerColumn]);

        if (marker < 0)
            break;

        double w, z;

        if (!useStandardErrors) {
            long double p = tokens[pvalueColumn].AsLongDouble();

            if (p <= 0 || p > 1.0)
                continue;

            z = -ninv(p * 0.5);
            w = weightColumn >= 0 ? tokens[weightColumn].AsDouble() : history->weight;

            if (!history->logTransform && effectColumn >= 0 && tokens[effectColumn][0] == '-')
                z *= -1;

            if (history->logTransform) {
                double eff = tokens[effectColumn].AsDouble();

                if (eff <= 0.0)
                    continue;

                if (eff < 1.0)
                    z *= -1;
            }
        } else {
            double eff = tokens[effectColumn].AsDouble();
            double sd = tokens[stderrColumn].AsDouble();

            if (history->logTransform) {
                if (eff <= 0.0)
                    continue;

                eff = log(eff);
            }

            if (sd <= 0)
                continue;

            z = eff;
            w = 1.0 / (sd * sd);
        }

        double freq = 0.0;
        if (useFrequencies)
            freq = tokens[freqColumn].AsDouble();

        if (flip)
            z *= -1;

        if (firstColumn >= 0 && secondColumn >= 0) {
            NumbersToLetters(tokens[firstColumn]);
            NumbersToLetters(tokens[secondColumn]);

            if (useStrand && tokens[strandColumn] == "-") {
                FlipAllele(tokens[firstColumn]);
                FlipAllele(tokens[secondColumn]);
            }

            FlipAlleles(tokens[firstColumn], tokens[secondColumn], z, freq);

            if (allele1[marker] == "")
                break;
            else if (tokens[firstColumn] != allele1[marker] ||
                     tokens[secondColumn] != allele2[marker])
                continue;
        }

        if (direction[marker] != '?')
            continue;

        direction[marker] = z == 0.0 ? '0' : (z > 0.0 ? '+' : '-');

        if (!useStandardErrors) {
            if (weights[marker] == 0.0) continue;

            double ez = sqrt(w) * statistics[marker] / weights[marker];

            z /= sqrt(history->genomicControl);

            hetStatistic[marker] += (z - ez) * (z - ez);
            hetDegreesOfFreedom[marker]++;
        } else {
            double e = statistics[marker] / weights[marker];

            hetStatistic[marker] += (z - e) * (z - e) * w / history->genomicControl;
            hetDegreesOfFreedom[marker]++;
        }

        processedMarkers++;
    }

    history->filterLabel.Swap(filterLabel);
    history->filterColumn.Swap(filterColumn);
    history->filterCondition.Swap(filterCondition);
    history->filterValue.Swap(filterValue);
    history->filterAlternate.Swap(filterAlternate);
    history->filterSets.Swap(filterSets);
    history->filterCounts.Swap(filterCounts);

    ifclose(f);

    if (processedMarkers != history->processedMarkers) {
        printf("## ERROR: Input file has changed since analysis started\n");
        return false;
    }

    return true;
}


void ShowHelp(bool startup)
{
    const char * setting = startup ? "default" : "setting";

    String gc;

    if (genomicControl && genomicControlLambda != 0.0)
        gc.printf("%.3f", genomicControlLambda);
    else if (genomicControl && genomicControlFilter.Length())
        gc.printf("LIST WITH %d MARKERS", genomicControlCount);
    else
        gc = genomicControl ? "ON" : "OFF";

    printf("# This program faciliates meta-analysis of genome-wide association studies.\n"
                   "# Commonly used commands are listed below:\n"
                   "#\n"
                   "# Options for describing input files ...\n"
                   "#   SEPARATOR        [WHITESPACE|COMMA|BOTH|TAB] (default = WHITESPACE)\n"
                   "#   COLUMNCOUNTING   [STRICT|LENIENT]            (%s = '%s')\n"
                   "#   MARKERLABEL      [LABEL]                     (%s = '%s')\n"
                   "#   ALLELELABELS     [LABEL1 LABEL2]             (%s = '%s','%s')\n"
                   "#   EFFECTLABEL      [LABEL|log(LABEL)]          (%s = '%s')\n"
                   "#   FLIP\n"
                   "#\n"
                   "# Options for filtering input files ...\n"
                   "#   ADDFILTER        [LABEL CONDITION VALUE]     (example = ADDFILTER N > 10)\n"
                   "#                    (available conditions are <, >, <=, >=, =, !=, IN)\n"
                   "#   REMOVEFILTERS\n"
                   "#\n"
                   "# Options for sample size weighted meta-analysis ...\n"
                   "#   WEIGHTLABEL      [LABEL]                     (%s = '%s')\n"
                   "#   PVALUELABEL      [LABEL]                     (%s = '%s')\n"
                   "#   DEFAULTWEIGHT    [NUMBER]                    (%s = %.1f)\n"
                   "#   MINWEIGHT        [NUMBER]                    (%s = %.1f)\n"
                   "#\n"
                   "# Options for inverse variance weighted meta-analysis ...\n"
                   "#   STDERRLABEL      [LABEL]                     (%s = '%s')\n"
                   "#   SCHEME           [SAMPLESIZE|STDERR]         (%s = %s)\n"
                   "#\n"
                   "# Options to enable tracking of allele frequencies ...\n"
                   "#   AVERAGEFREQ      [ON|OFF]                    (%s = %s)\n"
                   "#   MINMAXFREQ       [ON|OFF]                    (%s = %s)\n"
                   "#   FREQLABEL        [LABEL]                     (%s = '%s')\n"
                   "#\n"
                   "# Options to enable tracking of user defined variables ...\n"
                   "#   CUSTOMVARIABLE   [VARNAME]\n"
                   "#   LABEL            [VARNAME] AS [HEADER]\n"
                   "#\n"
                   "# Options to enable tracking of chromosomes and positions ...\n"
                   "#   TRACKPOSITIONS   [ON|OFF]                    (%s = %s\n"
                   "#   CHROMOSOMELABEL  [LABEL]                     (%s = '%s')\n"
                   "#   POSITIONLABEL    [LABEL]                     (%s = '%s')\n"
                   "#\n"
                   "# Options to enable explicit strand information ...\n"
                   "#   USESTRAND        [ON|OFF]                    (%s = %s)\n"
                   "#   STRANDLABEL      [LABEL]                     (%s = '%s')\n"
                   "#\n"
                   "# Automatic genomic control correction of input statistics ...\n"
                   "#   GENOMICCONTROL   [ON|OFF|VALUE|LIST snps.txt](%s = %s)\n"
                   "#\n"
                   "# Options to account for samples overlap ...\n"
                   "#   OVERLAP          [ON|OFF]                    (%s = %s)\n"
                   "#   ZCUTOFF          [NUMBER]                    (%s = %.1f)\n"
                   "#\n"
                   "# Options for general analysis control ...\n"
                   "#   PROCESSFILE            [FILENAME]\n"
                   "#   OUTFILE                [PREFIX SUFFIX]       (default = 'METAANALYSIS','.TBL')\n"
                   "#   MAXWARNINGS            [NUMBER]              (%s = %d)\n"
                   "#   VERBOSE                [ON|OFF]              (%s = '%s')\n"
                   "#   LOGPVALUE              [ON|OFF]              (%s = '%s')\n"
                   "#   EFFECT_PRINT_PRECISION [NUMBER]              (%s = '%d')\n"
                   "#   STDERR_PRINT_PRECISION [NUMBER]              (%s = '%d')\n"
                   "#   ANALYZE                [HETEROGENEITY]\n"
                   "#   CLEAR\n\n"
                   "# Options for general run control ...\n"
                   "#   SOURCE           [SCRIPTFILE]\n"
                   "#   RETURN\n"
                   "#   QUIT\n\n",
           setting, (const char *) strictColumnCounting ? "STRICT" : "LENIENT",
           setting, (const char *) markerLabel,
           setting, (const char *) firstAllele, (const char *) secondAllele,
           setting, (const char *) effectLabel,
           setting, (const char *) weightLabel,
           setting, (const char *) pvalueLabel,
           setting, weight,
           setting, minweight,
           setting, (const char *) stderrLabel,
           setting, useStandardErrors ? "STDERR" : "SAMPLESIZE",
           setting, averageFrequencies ? "ON" : "OFF",
           setting, minMaxFrequencies ? "ON" : "OFF",
           setting, (const char *) frequencyLabel,
           setting, trackPositions ? "ON" : "OFF",
           setting, (const char *) chromosomeLabel,
           setting, (const char *) positionLabel,
           setting, useStrand ? "ON" : "OFF",
           setting, (const char *) strandLabel,
           setting, (const char *) gc,
           setting, "OFF",
           setting, zCutoff,
           setting, maxWarnings,
           setting, verbose ? "ON" : "OFF",
           setting, logPValue ? "ON" : "OFF",
           setting, effectPrintPrecision,
           setting, stderrPrintPrecision);
}

void RunScript(FILE * file)
{
    String input;
    StringArray tokens;

    while (!feof(file)) {
        if (file == stdin)
            input.ReadLine().Trim();
        else
            input.ReadLine(file).Trim();

        tokens.ReplaceTokens(input);

        if (tokens.Length() == 0)
            continue;

        if (input[0] == '#') continue;

        if ((tokens[0].MatchesBeginningOf("ANALYZE") == 0 || tokens[0].MatchesBeginningOf("ANALYSE") == 0) && tokens[0].Length() > 1) {
            Analyze(tokens.Length() > 1 && tokens[1].MatchesBeginningOf("HETEROGENEITY") == 0);
            continue;
        }

        if (tokens[0].MatchesBeginningOf("CLEAR") == 0 && tokens[0].Length() > 1) {
            ClearAll();
            continue;
        }

        if (tokens[0].MatchesBeginningOf("REMOVEFILTERS") == 0 && tokens[0].Length() > 2) {
            ClearFilters();
            continue;
        }

        if (tokens[0].MatchesBeginningOf("FLIP") == 0 && tokens[0].Length() > 1) {
            flip = !flip;
            printf("## All effects will %s flipped\n", flip ? "be" : "not be");
            continue;
        }

        if (tokens[0].MatchesBeginningOf("HELP") == 0) {
            ShowHelp();
            continue;
        }

        if (tokens[0].MatchesBeginningOf("QUIT") == 0 ||
            tokens[0].MatchesBeginningOf("EXIT") == 0 && tokens[0].Length() > 1) {
            ClearAll();
            exit(0);
        }

        if (tokens[0].MatchesBeginningOf("RETURN") == 0 && tokens[0].Length() > 2)
            break;

        if (tokens.Length() == 1) {
            printf("##  ERROR: The command you issued could not be processed ...\n");
            continue;
        }

        if (tokens[0].MatchesBeginningOf("COLUMNCOUNTING") == 0 && tokens[0].Length() > 1) {
            if (tokens[1].MatchesBeginningOf("STRICT") == 0) {
                strictColumnCounting = true;
                printf("## STRICT MODE: Every input line must have exactly the same number of columns\n");
                continue;
            }
            if (tokens[1].MatchesBeginningOf("LENIENT") == 0) {
                strictColumnCounting = false;
                printf("## LENIENT MODE: Input lines can include extra columns, which will be ignored\n");
                continue;
            }
        }

        if (tokens[0].MatchesBeginningOf("CUSTOMVARIABLE") == 0 && tokens[0].Length() > 1) {
            if (markerLookup.Entries()) {
                printf("## ERROR: Meta-analysis in progress - before creating custom variables, use CLEAR command\n");
                continue;
            }

            if (customVariables.SlowFind(tokens[1]) >= 0)
                printf("## Variable '%s' already defined\n", (const char *) tokens[1]);
            else {
                printf("## Created custom variable '%s'\n", (const char *) tokens[1]);
                customVariables.Push(tokens[1]);
                customLabels.Push(tokens[1]);
            }
            continue;
        }

        if (tokens[0].MatchesBeginningOf("DEFAULTWEIGHT") == 0) {
            weight = tokens[1].AsDouble();
            printf("## Set default weight to %.2f ...\n", weight);
            continue;
        }

        if (tokens[0].MatchesBeginningOf("WEIGHTLABEL") == 0) {
            weightLabel = tokens[1];
            printf("## Set weight header to %s ...\n", (const char *) weightLabel);
            continue;
        }

        if (tokens[0].MatchesBeginningOf("GENOMICCONTROL") == 0) {
            if (tokens[1] == "ON") {
                genomicControl = true;
                genomicControlFilter.Clear();
                genomicControlLambda = 0.0;
                printf("## Genomic control correction of input statistics enabled\n");
                continue;
            } else if (tokens[1].MatchesBeginningOf("OFF") == 0 && tokens[1].Length() > 1) {
                genomicControl = false;
                genomicControlFilter.Clear();
                printf("## Genomic control correction of input statistics disabled\n");
                continue;
            } else if (tokens[1].MatchesBeginningOf("LIST") == 0 && tokens.Length() == 3) {
                genomicControl = true;
                genomicControlCount = 0;

                String label;
                StringArray markerList;
                markerList.Read(tokens[2]);

                for (int i = 0; i < markerList.Length(); i++) {
                    label = markerList[i].Trim();

                    if (label.Length() == 0) continue;

                    GetMarkerId(label);
                    genomicControlCount++;
                }

                if (genomicControlCount == 0) {
                    genomicControl = false;
                    genomicControlFilter.Clear();

                    printf("## Genomic control disabled, no markers listed in [%s]\n", (const char *) tokens[2]);
                    continue;
                }

                genomicControlFilter.Fill('.', markerLookup.Entries());

                for (int i = 0; i < markerList.Length(); i++) {
                    label = markerList[i].Trim();

                    if (label.Length() == 0) continue;

                    genomicControlFilter[GetMarkerId(label)] = 'Y';
                }

                printf("## Genomic control correction will be based on %d markers listed in [%s]\n",
                       genomicControlCount, (const char *) tokens[2]);
                continue;
            } else if (tokens[1].AsDouble() > 0.0) {
                genomicControl = true;
                genomicControlLambda = tokens[1].AsDouble();
                genomicControlFilter.Clear();
                printf("## Genomic control parameter set to %.3f\n", genomicControlLambda);
                continue;
            }
        }

        if (tokens[0].MatchesBeginningOf("LOGPVALUE") == 0) {
            if (tokens[1] == "ON") {
                logPValue = true;
                printf("## Log(p-value) will be output\n");
                continue;
            } else if (tokens[1].MatchesBeginningOf("OFF") == 0 && tokens[1].Length() > 1) {
                logPValue = false;
                printf("## Untransformed p-value will be output\n");
                continue;
            }
        }

        if (tokens[0].MatchesBeginningOf("USESTRAND") == 0 && tokens[0].Length() > 3) {
            if (tokens[1] == "ON") {
                useStrand = true;
                printf("## Strand flipping according to strand column enabled\n");
                continue;
            } else if (tokens[1].MatchesBeginningOf("OFF") == 0 && tokens[1].Length() > 1) {
                useStrand = false;
                printf("## Strand flipping according to strand column disabled\n");
                continue;
            }
        }

        if (tokens[0].MatchesBeginningOf("VERBOSE") == 0) {
            if (tokens[1] == "ON") {
                verbose = true;
                printf("## Verbose output during meta-analysis enabled\n");
                continue;
            } else if (tokens[1].MatchesBeginningOf("OFF") == 0 && tokens[1].Length() > 1) {
                verbose = false;
                printf("## Verbose output during meta-analysis disabled\n");
                continue;
            }
        }

        if (tokens[0].MatchesBeginningOf("AVERAGEFREQUENCY") == 0 && tokens[0].Length() > 1) {
            if (markerLookup.Entries()) {
                printf("## ERROR: Meta-analysis in progress - before issuing this command, use CLEAR command\n");
                continue;
            }

            if (tokens[1] == "ON") {
                averageFrequencies = true;
                printf("## Averaging of allele frequencies enabled\n");
                continue;
            } else if (tokens[1].MatchesBeginningOf("OFF") == 0 && tokens[1].Length() > 1) {
                averageFrequencies = false;
                printf("## Averaging of allele frequencies disabled\n");
                continue;
            }
        }

        if (tokens[0].MatchesBeginningOf("MINMAXFREQUENCY") == 0 && tokens[0].Length() > 3) {
            if (markerLookup.Entries()) {
                printf("## ERROR: Meta-analysis in progress - before issuing this command, use CLEAR command\n");
                continue;
            }

            if (tokens[1] == "ON") {
                minMaxFrequencies = true;
                printf("## Tracking of extreme allele frequencies enabled\n");
                continue;
            } else if (tokens[1].MatchesBeginningOf("OFF") == 0 && tokens[1].Length() > 1) {
                minMaxFrequencies = false;
                printf("## Tracking of extreme allele frequencies disabled\n");
                continue;
            }
        }

        if (tokens[0].MatchesBeginningOf("MARKERLABEL") == 0 &&
            tokens[0].Length() > 2) {
            markerLabel = tokens[1];
            printf("## Set marker header to %s ...\n", (const char *) markerLabel);
            continue;
        }

        if (tokens[0].MatchesBeginningOf("MAXWARNINGS") == 0 &&
            tokens[0].Length() > 2) {
            maxWarnings = tokens[1].AsInteger();
            printf("## Set the maximum number of warnings per file to %d ...\n", maxWarnings);
            continue;
        }

        if (tokens[0].MatchesBeginningOf("PVALUELABEL") == 0 && tokens[0].Length() > 1) {
            pvalueLabel = tokens[1];
            printf("## Set p-value header to %s ...\n", (const char *) pvalueLabel);
            continue;
        }

        if (tokens[0].MatchesBeginningOf("EFFECTLABEL") == 0 && tokens[0].Length() > 1) {
            effectLabel = tokens[1];
            printf("## Set effect header to %s ...\n", (const char *) effectLabel);
            continue;
        }

        if (tokens[0].MatchesBeginningOf("FREQLABEL") == 0 && tokens[0].Length() > 1) {
            frequencyLabel = tokens[1];
            printf("## Set frequency header to %s ...\n", (const char *) frequencyLabel);
            if (!averageFrequencies)
                printf("## If you want frequencies to be averaged, issue the 'AVERAGEFREQ ON' command\n");
            continue;
        }

        if (tokens[0].MatchesBeginningOf("SOURCE") == 0 && tokens[0].Length() > 1) {
            FILE *ifile = fopen(tokens[1], "rt");

            if (ifile != NULL) {
                printf("# Processing commands in %s ...\n", (const char *) tokens[1]);

                RunScript(ifile);
                fclose(ifile);
                continue;
            } else {
                printf("# ERROR: Failed to open file %s ...\n", (const char *) tokens[1]);
                continue;
            }
        }

        if (tokens[0].MatchesBeginningOf("STDERRLABEL") == 0 && tokens[0].Length() > 2) {
            stderrLabel = tokens[1];
            printf("## Set standard error header to %s ...\n", (const char *) stderrLabel);
            continue;
        }

        if (tokens[0].MatchesBeginningOf("STRANDLABEL") == 0 && tokens[0].Length() > 2) {
            strandLabel = tokens[1];
            printf("## Set strand header to %s ...\n", (const char *) strandLabel);
            if (!useStrand) printf("## If you want strand information to be used, issue the 'USESTRAND ON' command\n");
            continue;
        }

        if (tokens[0].MatchesBeginningOf("CHROMOSOMELABEL") == 0 && tokens[0].Length() > 2) {
            chromosomeLabel = tokens[1];
            printf("## Set chromosome header to %s ...\n", (const char*) chromosomeLabel);
            continue;
        }

        if (tokens[0].MatchesBeginningOf("POSITIONLABEL") == 0 && tokens[0].Length() > 2) {
            positionLabel = tokens[1];
            printf("## Set position header to %s ...\n", (const char*) positionLabel);
            continue;
        }

        if (tokens[0].MatchesBeginningOf("SCHEME") == 0 &&
            tokens[0].Length() > 1)
            if (markerLookup.Entries()) {
                printf("## ERROR: Meta-analysis in progress - before changing analysis scheme, use CLEAR command\n");
                continue;
            } else if (tokens[1].MatchesBeginningOf("SAMPLESIZE") == 0 &&
                       tokens[1].Length() > 1) {
                useStandardErrors = false;
                printf("## Meta-analysis will be based on sample sizes, p-values and direction of effect ...\n");
                continue;
            } else if (tokens[1].MatchesBeginningOf("STDERR") == 0 &&
                       tokens[1].Length() > 1) {
                if (studyOverlap) {
                    printf("## ERROR: STDERR can't be used together with enabled OVERLAP\n");
                    continue;
                }
                useStandardErrors = true;
                printf("## Meta-analysis will be based on effect sizes and their standard errors ...\n");
                continue;
            }

        if (tokens[0].MatchesBeginningOf("SEPARATOR") == 0 &&
            tokens[0].Length() > 1)
            if (tokens[1].MatchesBeginningOf("WHITESPACE") == 0) {
                separators = " \t";
                printf("## Set column separator to WHITESPACE ...\n");
                continue;
            } else if (tokens[1].MatchesBeginningOf("COMMAS") == 0) {
                separators = ",";
                printf("## Set column separator to COMMAS ...\n");
                continue;
            } else if (tokens[1].MatchesBeginningOf("TABS") == 0) {
                separators = "\t";
                printf("## Set column separator to TAB ...\n");
                continue;
            } else if (tokens[1].MatchesBeginningOf("BOTH") == 0) {
                separators = " \t,";
                printf("## Set column separator to BOTH ...\n");
                continue;
            }

        if (tokens[0].MatchesBeginningOf("PROCESSFILE") == 0 && tokens[0].Length() > 1) {
            ProcessFile(tokens[1], processedFiles = new FileSummary(processedFiles));
            continue;
        }

        if (tokens[0].MatchesBeginningOf("MINWEIGHT") == 0 &&
            tokens[0].Length() > 3) {
            minweight = tokens[1].AsDouble();
            printf("## Set minimum weight for meta-analysis to %.2f ...\n", minweight);
            continue;
        }

        if (tokens[0].MatchesBeginningOf("OVERLAP") == 0) {
            if (tokens[1] == "ON") {
                if (useStandardErrors) {
                    printf("## ERROR: OVERLAP can't be used if meta-analysis is based on effect sizes and their standard errors\n");
                    continue;
                }
                studyOverlap = true;
                printf("## Correction for samples overlap enabled\n");
                continue;
            } else if (tokens[1].MatchesBeginningOf("OFF") == 0 && tokens[1].Length() > 1) {
                studyOverlap = false;
                printf("## Correction for samples overlap disabled\n");
                continue;
            }
        }

        if (tokens[0].MatchesBeginningOf("TRACKPOSITIONS") == 0) {
            if (markerLookup.Entries()) {
                printf("## ERROR: Meta-analysis in progress - before turning on/off tracking of chromosomes and positions, use CLEAR command\n");
                continue;
            } else if (tokens[1] == "ON"){
                trackPositions = true;
                printf("## Tracking of chromosomes and positions is enabled\n");
                continue;
            } else if (tokens[1].MatchesBeginningOf("OFF") == 0 && tokens[1].Length() > 1) {
                trackPositions = false;
                printf("## Tracking of chromosomes and positions is disabled\n");
                continue;
            }
        }

        if ((tokens[0].MatchesBeginningOf("EFFECT_PRINT_PRECISION") == 0) && (tokens.Length() > 1)) {
            effectPrintPrecision = tokens[1].AsInteger();
            printf("## Set print pecision for Effect to %d ...\n", effectPrintPrecision);
            continue;
        }

        if ((tokens[0].MatchesBeginningOf("STDERR_PRINT_PRECISION") == 0) && (tokens.Length() > 1)) {
            stderrPrintPrecision = tokens[1].AsInteger();
            printf("## Set print pecision for StdErr to %d ...\n", stderrPrintPrecision);
            continue;
        }

        if (tokens[0].MatchesBeginningOf("ZCUTOFF") == 0) {
            zCutoff = tokens[1].AsDouble();
            printf("## Set Z cutoff to %.2f ...\n", zCutoff);
            continue;
        }

        if (tokens.Length() == 2) {
            printf("## ERROR: The command you issued could not be processed ...\n");
            continue;
        }

        if (tokens[0].MatchesBeginningOf("ALLELELABELS") == 0 &&
            tokens[0].Length() > 1) {
            firstAllele = tokens[1];
            secondAllele = tokens[2];
            printf("## Set allele headers to %s and %s ...\n",
                   (const char *) firstAllele, (const char *) secondAllele);
            continue;
        }

        if (tokens[0].MatchesBeginningOf("OUTFILE") == 0) {
            outfileround = 1;
            outfile = tokens[1] + "%d" + tokens[2];
            printf("## Set output file prefix and suffix to %s and %s ...\n",
                   (const char *) tokens[1], (const char *) tokens[2]);
            continue;
        }

        if (tokens[0].MatchesBeginningOf("ADDFILTER") == 0 &&
            tokens[0].Length() > 2) {
            AddFilter(tokens);
            continue;
        }

        if (tokens[0].MatchesBeginningOf("LABEL") == 0 &&
            tokens[2].MatchesBeginningOf("AS") == 0 &&
            tokens.Length() > 3) {
            int customId = customVariables.SlowFind(tokens[1]);

            if (customId < 0) {
                printf("## ERROR: Custom variable '%s' is undefined\n", (const char *) tokens[0]);
                continue;
            }

            customLabels[customId] = tokens[3];
            printf("## Set header for '%s' to '%s'\n", (const char *) tokens[1], (const char *) tokens[3]);
            continue;
        }

        printf("## Command not recognized - type HELP to list available commands\n");
    }
}

int main(int argc, char ** argv)
{
   // suggested by Albert Vernon Smith to facilitate monitoring of background jobs
   setvbuf(stdout, NULL, _IONBF, 0);
   setvbuf(stderr, NULL, _IONBF, 0);

   SetupCrashHandlers();

   String::caseSensitive = false;

   printf("MetaAnalysis Helper - (c) 2007 - 2009 Goncalo Abecasis\n");

#ifdef VERSION
   printf("This version released on %s\n", VERSION);
#endif

   printf("\n");

   ShowHelp(true);

   if (argc <= 1)
      RunScript(stdin);
   else
      for (int i = 1; i < argc; i++)
      {
         FILE * ifile = fopen(argv[i], "rt");

         if (ifile != NULL)
         {
            printf("# Processing commands in %s ...\n", (const char *) argv[i]);

            RunScript(ifile);
            fclose(ifile);
            continue;
         }
         else
         {
            printf("# Error opening file %s ...\n", (const char *) argv[i]);
            continue;
         }
      }

   ClearAll();
}
 
