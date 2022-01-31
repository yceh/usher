/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   Interface.cpp
 * CREATED ON:  04 August 2011, 11:33
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 *
 * DESCRIPTION: See "Interface.h"
 *
 * HISTORY:     Version     Date            Description
 *              1.0         2011-08-04      Created
 *
 * VERSION:     1.0
 * LAST EDIT:   04 August 2011
 */

#include "Interface.h"
#include "util/Util.h"
#include "stat/PTableFile.h"


////////////////////////////////////////////////////////////////////////////////
// Static variables
////////////////////////////////////////////////////////////////////////////////

#ifdef __CYGWIN__
const string Interface::SEPARATOR
        = "--------------------------------------------------------------------------------";
const string Interface::BOLD_SEPARATOR
        = "--------------------------------------------------------------------------------";
#else
const string Interface::SEPARATOR
        = "────────────────────────────────────────────────────────────────────────────────";

const string Interface::BOLD_SEPARATOR
        = "════════════════════════════════════════════════════════════════════════════════";
#endif /*__CYGWIN__*/

const string Interface::DEFAULT_INDENT = "  "; // 2 spaces

const string Interface::YES_STR = " - YES";
const string Interface::NO_STR = " - NO";
const string Interface::WARNING_STR = "WARNING: ";
const string Interface::ERROR_STR = "ERROR: ";


/* Colours      black   red     green   brown   blue    purple  cyan    gray
 * Foreground   30      31      32      33      34      35      36      37
 * Background   40      41      42      43      44      45      46      47
 *
 * Style    plain   bold    underline   strike-through
 *          0       1       4           9
 *
 * Reversed 7
 */
const string Interface::ERROR_FORMAT = "31";
const string Interface::WARNING_FORMAT = "33";
const string Interface::OUTPUT_FORMAT = "1;32";


////////////////////////////////////////////////////////////////////////////////

string Interface::getProgramDir(void) const {
    string result = AppInfo::instance()->getExecutablePath();
    
    if (result.length() > 0) {
        unsigned long index = result.length();
        while (result[index - 1] != '/' && index > 0) index--;
        
        if (index < result.length()) {
            result = result.substr(0, index);
        }
    }
    
    if (result.length() <= 0) result = "./";
    
    return result;
}

void Interface::initLogFile(string logFileName) {
    TextFile * newLogFile = new TextFile(logFileName);
    
    if (newLogFile->exists()) {
        logFile = nullptr; // temporarily suppress writing the log file
        
        clear();
        (*this) << "The file \"" << newLogFile->getPath() << "\" already exists.\n"
                << "Do you want to overwrite it?";
        if (!showWarning(false)) {
            (*this) << "No log file will be created.";
            logFile = nullptr;
            
        } else {
            logFile = newLogFile;
            logFile->removeFile();
        }
        
        clear();
        
    } else {
        logFile = newLogFile;
    }
    
    if (logFile != nullptr) {
        logFile->openToWrite();
        logFile->writeLine(BOLD_SEPARATOR);
        logFile->writeLine(runDescription);
        logFile->writeLine(BOLD_SEPARATOR + "\n");
        
        if (runCommand.length() > 0) {
            logFile->writeLine("COMMAND:");
            logFile->writeLine(DEFAULT_INDENT + runCommand + "\n");
        }
    }
}

const string Interface::align(string firstIndent) {
    string alignedStr = "";
    string nextIndent = "";
    
    for (unsigned long i = 0; i < firstIndent.length(); i++) {
        if (firstIndent[i] == '\n' || firstIndent[i] == '\r') {
            nextIndent = "";
        } else if (firstIndent[i] == '\t') {
            nextIndent += "\t";
        } else {
            nextIndent += " "; // add a white space
        }
    }
    
    stringstream tmpStream;
    tmpStream.str(this->str());
    string newLine = "";
    int lineCount = 0;
    
    do {
        std::getline(tmpStream, newLine);
        
        if (lineCount > 0) {
            alignedStr += "\n";
        }
        
        if (newLine.length() > 0) {
            if (lineCount <= 0) {
                alignedStr += firstIndent + newLine;
            } else {
                alignedStr += nextIndent + newLine;
            }
        }
        
        lineCount++;
    } while (!tmpStream.eof() && !tmpStream.fail());
    
    return alignedStr;
}

void Interface::startProgram(string startMessage) {
    clear();
    
    runDescription = AppInfo::instance()->getApplicationName();
    if (startMessage.length() > 0) {
        runDescription += " -- " + startMessage;
    }
    
    switch (mode) {
        case CONSOLE:
            (*this) << endl
                    << BOLD_SEPARATOR << endl
                    << runDescription << endl
                    << BOLD_SEPARATOR << endl;
            cout << this->align(DEFAULT_INDENT) << endl;
            break;
        
        
        case GUI:
            //TODO: Show output on GUI.
            assert(false);
            break;
        
        
        case SILENT:
            // Show nothing
            break;
    }
    
    clear();
    cout.flush();
}

void Interface::showPTableLink() {
    string lastUsedPTableFile = "";
    
    try {
        lastUsedPTableFile = Util::getConfig(PTableFile::getConfigKey());
        
        // Quickly validate the file
        fstream pTableFile(lastUsedPTableFile.c_str(), ios::in | ios::binary);
        if (!pTableFile.is_open() || pTableFile.fail() || pTableFile.eof()) {
            // File not valid
            lastUsedPTableFile = "";
        }
        pTableFile.close();
        
    } catch (...) {
        lastUsedPTableFile = "";
    }
    
    clear();
    
    switch (mode) {
        case CONSOLE:
            if (lastUsedPTableFile.length() > 0) {
                (*this) << "The current p-value table is located at: \""
                        << lastUsedPTableFile << "\""
                        << endl;
            } else {
                (*this) << "Currently, there is no p-value table associated with "
                        << AppInfo::instance()->getApplicationName() << "."
                        << endl;
            }
            
            cout << this->align(DEFAULT_INDENT) << endl;
            break;
        
        
        case GUI:
            //TODO: Show output on GUI.
            assert(false);
            break;
        
        
        case SILENT:
            // Show nothing
            break;
    }
    
    clear();
    cout.flush();
}

void Interface::throwExitSignal(int status) {
    try {
        clear();
        
        switch (mode) {
            case CONSOLE:
                (*this) << BOLD_SEPARATOR << endl;
                cout << this->align(DEFAULT_INDENT) << endl << endl;
                break;
            
            
            case GUI:
                //TODO: Show output on GUI.
                assert(false);
                break;
            
            
            case SILENT:
                // Show nothing
                break;
        }
        
        if (logFile != nullptr) {
            time_t currentTime;
            struct tm * timeInfo;
            string timeStr;
            time(&currentTime);
            timeInfo = localtime(&currentTime);
            
            timeStr = String::formatInt(timeInfo->tm_year + 1900, 4) + "/"
                      + String::formatInt(timeInfo->tm_mon + 1, 2) + "/"
                      + String::formatInt(timeInfo->tm_mday, 2)
                      + " "
                      + String::formatInt(timeInfo->tm_hour, 2) + ":"
                      + String::formatInt(timeInfo->tm_min, 2) + ":"
                      + String::formatInt(timeInfo->tm_sec, 2);
            
            logFile->writeLine(BOLD_SEPARATOR + "\n\n");
            logFile->writeLine(
                    "                                            ╔══════════════════════════════════╗");
            logFile->writeLine(
                    "                                            ║ Log saved at " + timeStr + " ║");
            logFile->writeLine(
                    "                                            ╚══════════════════════════════════╝");
            logFile->writeLine();
            logFile->close();
        }
        
        clear();
        cout.flush();
        
        
    } catch (...) {
        /* Do nothing.
         * The reason for using try-catch is about to make sure that the
         * exitSignal will be throw.
         */
    }
    
    
    if (status == 0) {
        throw SUCCESS;
    } else {
        throw FAIL;
    }
}

void Interface::showOutput(bool autoEndl) {
    assert(!isCounting);
    
    switch (mode) {
        case CONSOLE: // fall through
        case SILENT:
            setOutputStreamFormat(OUTPUT_FORMAT);
            cout << this->align(DEFAULT_INDENT);
            if (autoEndl) {
                cout << endl;
            }
            cleanOutputStreamFormat();
            break;
        
        
        case GUI:
            //TODO: Show output on GUI.
            assert(false);
            break;
    }
    
    if (logFile != nullptr) {
        if (autoEndl)
            logFile->writeLine(this->str());
        else
            logFile->write(this->str());
    }
    
    clear();
    cout.flush();
}

void Interface::showLog(bool autoEndl) {
    assert(!isCounting);
    
    if (isLogStreamEnabled) {
        switch (mode) {
            case CONSOLE:
                cout << this->align(DEFAULT_INDENT);
                if (autoEndl) {
                    cout << endl;
                }
                break;
            
            
            case GUI:
                //TODO: Show log message on GUI.
                assert(false);
                break;
            
            
            case SILENT:
                // no log
                break;
        }
    }
    
    if (logFile != nullptr && isLogFileEnabled) {
        if (autoEndl)
            logFile->writeLine(this->str());
        else
            logFile->write(this->str());
    }
    
    clear();
    cout.flush();
}

const bool Interface::showWarning(bool ignoreAndCont /* = true */) {
    assert(!isCounting);
    
    bool option = true;
    string answerStr = "";
    
    switch (mode) {
        case CONSOLE:
            setOutputStreamFormat(WARNING_FORMAT);
            if (ignoreAndCont) {
                cout << this->align(DEFAULT_INDENT + WARNING_STR) << endl << endl;
            } else {
                (*this) << " (Y/N)";
                cout << this->align(DEFAULT_INDENT + WARNING_STR);
                
                int answerKey = 0;
                do {
                    answerKey = Util::getKey();
                    if (answerKey > 0) {
                        answerKey = toupper(answerKey);
                    }
                } while (answerKey != 'Y' && answerKey != 'N');
                
                if (answerKey == 'Y') {
                    answerStr = YES_STR;
                    option = true;
                } else {
                    answerStr = NO_STR;
                    option = false;
                }
                
                cout << answerStr << endl << endl;
            }
            
            cleanOutputStreamFormat();
            break;
        
        
        case GUI:
            //TODO: Show warning message on GUI.
            break;
        
        
        case SILENT:
            // no warning
            break;
    }
    
    if (logFile != nullptr) {
        logFile->writeLine(this->align(WARNING_STR) + answerStr + "\n");
    }
    
    clear();
    cout.flush();
    
    return option;
}

void Interface::showError(bool autoEndl, bool exitAfterShowing) {
    if (isCounting) {
        count(counter.currentValue, true); // Break counting when error occurs
    }
    
    switch (mode) {
        case CONSOLE: // fall through
        case SILENT:
            setErrorStreamFormat(ERROR_FORMAT);
            cerr << this->align(DEFAULT_INDENT + ERROR_STR);
            if (autoEndl) {
                cerr << endl;
            }
            cleanErrorStreamFormat();
            break;
        
        
        case GUI:
            //TODO: Show error on GUI.
            assert(false);
            break;
    }
    
    if (logFile != nullptr) {
        if (autoEndl)
            logFile->writeLine(this->align(ERROR_STR));
        else
            logFile->write(this->align(ERROR_STR));
    }
    
    clear();
    cerr.flush();
    
    if (exitAfterShowing) {
        throwExitSignal(1); // exit on error
    }
}

bool Interface::yesNoQuestion(ThreeStateOpt defaultAnswer) {
    assert(!isCounting);
    bool answer = false;
    string answerStr = "";
    
    
    if (mode == CONSOLE ||
        (mode == SILENT && defaultAnswer == UNKNOWN)) {
        (*this) << " (Y/N)";
        cout << this->align(DEFAULT_INDENT);
        
        int answerKey = 0;
        do {
            answerKey = Util::getKey();
            if (answerKey > 0) {
                answerKey = toupper(answerKey);
            }
        } while (answerKey != 'Y' && answerKey != 'N');
        
        if (answerKey == 'Y') {
            answerStr = YES_STR;
            answer = true;
        } else {
            answerStr = NO_STR;
            answer = false;
        }
        cout << answerStr << endl << endl;
        
        
    } else if (mode == SILENT) {
        answer = (defaultAnswer == YES);
        
    } else if (mode == GUI) {
        //TODO: Show question message on GUI.
        assert(false);
    }
    
    if (logFile != nullptr) {
        logFile->writeLine(this->str() + answerStr + "\n");
    }
    
    clear();
    cout.flush();
    
    return answer;
}

void Interface::count(double currentValue, bool isLastValue /*= false*/) {
    assert(isCounting);
    
    counter.currentValue = currentValue;
    string message = counter.showValue();
    
    if (counter.processingText.length() > 0) {
        message = counter.processingText + " :  " + message;
    }
    if (counter.isPercentCounter) {
        message += "%";
    }
    if (this->str().length() > 0) {
        message += " " + this->str();
    }
    if (isLastValue) {
        message += "\n";
    } else {
        message += "\r";
    }
    this->str(message);
    
    
    switch (mode) {
        case CONSOLE:
            cout << this->align(DEFAULT_INDENT);
            
            //if (isLastValue) {
            //    cout << endl;
            //}
            break;
        
        
        case GUI:
            //TODO: Show process counter on GUI.
            assert(false);
            break;
        
        
        case SILENT:
            // No counting message.
            break;
    }
    
    if (isLastValue) {
        isCounting = false;
        if (logFile != nullptr) {
            logFile->write(this->str());
        }
    }
    
    clear();
    cout.flush();
}

void Interface::initCounter(string processingText, double minValue,
                            double maxValue, bool percentCounter /*= true*/) {
    
    counter.initialize(processingText, minValue, maxValue, percentCounter);
    isCounting = true;
}

void Interface::finishCounting() {
    count(counter.maxValue, true);
}

void Interface::catchUsageModes(void) {
    std::string exeName = AppInfo::instance()->getExecutablePath();
    
    switch (mode) {
        case CONSOLE:
            (*this) << "USAGE:\n"
                    << "    " << exeName << " -<run_mode>  <run_options>\n"
                    << endl
                    << "There are currently 10 run modes:\n"
                    << "    " << exeName << " -help\n"
                    << endl
                    << "    " << exeName << " -version\n"
                    << endl
                    << "    " << exeName << " -gen-p    <ptable_file>  <ptable_size>\n"
                    << endl
                    << "    " << exeName << " -check-p  [ptable_file]\n"
                    << endl
                    << "    " << exeName << " -write    <seq_file>                [options]\n"
                    << "    " << exeName << " -write    <seq_file> <output_file>  [options]\n"
                    << endl
                    << "    " << exeName << " -info     <seq_file>                [options]\n"
                    << endl
                    << "    " << exeName
                    << " -match    <seq_file> <seq_name> <minDistance> <maxDistance> [options]\n"
                    << endl
                    << "    " << exeName
                    << " -triplet  <seq_file>  <P_name>    <Q_name>      <C_name>    [options]\n"
                    << endl
                    << "    " << exeName << " -single   <seq_file>                 [options]\n"
                    << "    " << exeName << " -single   <parent_file> <child_file> [options]\n"
                    << endl
                    << "    " << exeName
                    << " -full     <seq_file>                 [-ptable ptable_file]  [options]\n"
                    << "    " << exeName
                    << " -full     <parent_file> <child_file> [-ptable ptable_file]  [options]\n";
            break;
        
        
        case GUI:
            //TODO: catch usage modes to show on GUI.
            assert(false);
            break;
        
        
        case SILENT:
            // No message.
            break;
    }
}

void Interface::catchVersion(void) {
    switch (mode) {
        case CONSOLE:
            (*this) << "Version " << AppInfo::instance()->getVersionFullString()
                    << endl << endl
                    //<< AppInfo::instance()->getCopyright() << endl
                    << AppInfo::instance()->getCitation();
            break;
        
        
        case GUI:
            //TODO: catch version info to shown on GUI.
            assert(false);
            break;
        
        
        case SILENT:
            // No message.
            break;
    }
}

void Interface::catchOptions(void) {
    std::string exeName = AppInfo::instance()->getExecutablePath();
    
    switch (mode) {
        case CONSOLE:
            (*this) << "OPTIONS:\n"
                    << endl
                    << "    -a       use all sites rather than just polymorphic sites; default is off.\n"
                    << endl
                    << "    -b -e    beginning sequence / end sequence; e.g. -b13 -e17 tests sequences\n"
                    << "             #13 through #17 inclusive (each as the child sequence) to test if\n"
                    << "             they are a recombinant of the remaining sequences in the data set.\n"
                    << endl
                    << "             This option can be used to parallelize the algorithm (e.g. -b0 -e9\n"
                    << "             on processor #1; -b10 -e19 on processor #2).\n"
                    << endl
                    << "             To test a single query sequence for recombination, simply use -b19\n"
                    << "             -e19.\n"
                    << endl
                    << "             Default is to test all sequences for recombination.\n"
                    << endl
                    << "    -d       distinct sequences only; removes sequences that are identical to\n"
                    << "             other sequences; default is off.\n"
                    << endl
                    << "    -f -l    first and last nucleotide to be analyzed, e.g. -f100 -l200\n"
                    << "             considers nucleotide positions 100-200 inclusive (first position is\n"
                    << "             1, not 0); default is all positions.\n"
                    << endl
                    << "    -12po    first and second positions only; strips sequences down to positions\n"
                    << "             1, 2, 4, 5, 7, 8 etc. Note that if your start codon is at positions\n"
                    << "             20-22, then you must use the option -f20 (and perhaps -l as well)\n"
                    << "             to ensure that you are in the right reading frame and that you are\n"
                    << "             looking at a coding region.\n"
                    << endl
                    << "    -3po     third positions only; strips sequences down to positions 3, 6, 9,\n"
                    << "             etc. Note that if your start codon is at positions 20-22, then you\n"
                    << "             must use the option -f20 (and perhaps -l as well) to ensure that\n"
                    << "             you are in the right reading frame and that you are looking at a\n"
                    << "             coding region.\n"
                    << endl
                    << "    -L       set the minimum length to count a segment as recombinant;\n"
                    << "             e.g. -L150 sets minimum length to 150 nucleotides; default is 100.\n"
                    << endl
                    << "    -r       record skipped computation to file \""
                    << config::DEFAULT_SKIPPED_TRIPLETS_FILE_NAME
                    << "\"; default is off.\n"
                    << endl
                    << "    -t       rejection threshold, e.g. -t0.01, -t1e-6; default is -t0.05.\n"
                    << endl
                    << "    -x       cut; this option can be used with the -f and -l switches to cut out\n"
                    << "             a portion of the sequence and use the remainder of the sequence for\n"
                    << "             analysis; example:\n"
                    << endl
                    << "                  " << exeName << " -full seq_file -f100 -l200 -x\n"
                    << endl
                    << "    -y       indicate YesNo-mode; algorithm stops once a significant triple has\n"
                    << "             been found; this is off by default.\n"
                    << endl
                    << "    -#       indicate that for multiple comparisons corrections, the number of\n"
                    << "             comparisons to be used is the actual number performed in the run;\n"
                    << "             when this option is off, the number of comparisons used in the\n"
                    << "             correction is N*(N-1)*(N-2) where N is the total number of\n"
                    << "             sequences; this options is off by default.\n"
                    << endl
                    << "    -nohs    suppress using Hogan-Siegmund approximations (Adv. Appl. Math.,\n"
                    << "             1986), when the exact p-values cannot be computed; default is off.\n"
                    << endl
                    << "    -fasta   format output as FASTA; default is PHYLIP format; used in write\n"
                    << "             mode only.\n"
                    << endl
                    << "    -nexus   format output as NEXUS; default is PHYLIP format; used in write\n"
                    << "             mode only.\n"
                    << endl
                    << "    -nr      suppress writing to \""
                    << config::DEFAULT_RECOMBINANTS_FILE_NAME << "\" file.\n"
                    << endl
                    << "    -id      give unique identifier for this run; e.g. when running the\n"
                    << "             algorithm in parallel, you can specify\n"
                    << endl
                    << "                 " << exeName << " -full seq_file -b0  -e9  -id run01\n"
                    << "                 " << exeName << " -full seq_file -b10 -e19 -id run02\n"
                    << endl
                    << "             and this will create output files called run01."
                    << config::DEFAULT_RECOMBINANTS_FILE_NAME
                    << ",\n"
                    << "             run02." << config::DEFAULT_RECOMBINANTS_FILE_NAME << ", run01."
                    << config::DEFAULT_PVAL_HISTOGRAM_FILE_NAME << ", run02."
                    << config::DEFAULT_PVAL_HISTOGRAM_FILE_NAME << ", etc.\n"
                    << endl
                    << "    -id-auto          automatically give a unique identifier for this run based\n"
                    << "                      on the run-time.\n"
                    << endl
                    << "    -bp-all           calculate all breakpoint positions for all candidate\n"
                    << "                      recombinants; this will significantly slow down the\n"
                    << "                      algorithm when the sequence data are highly recombinant;\n"
                    << "                      default is to calculate only the best breakpoint positions\n"
                    << "                      per candidate recombinant.\n"
                    << endl
                    << "    -bp-none          do not calculate any breakpoint positions; this mode will\n"
                    << "                      generate a file \""
                    << config::DEFAULT_RECOMBINANTS_FILE_NAME
                    << "\" with all candidate recombinant\n"
                    << "                      triplets but with no breakpoint information; default is to\n"
                    << "                      calculate the best breakpoint positions per candidate\n"
                    << "                      recombinant and include only these triplets in \""
                    << config::DEFAULT_RECOMBINANTS_FILE_NAME << "\".\n"
                    << endl
                    << "    -subset           designate subset of sequences to be analyzed; the subset\n"
                    << "                      file must contain sequence accessions, each separated by\n"
                    << "                      an end-line character; command-line usage is\n"
                    << endl
                    << "                          " << exeName
                    << " -full seq_file -subset subset_file\n"
                    << endl
                    << "                      where the file \"myfile.subset\" is text file with accession\n"
                    << "                      numbers separated by whitespace.\n"
                    << endl
                    << "    -subset-remove   as above, except that the sequences in the subset file\n"
                    << "                     \"subset_file\" are removed from the dataset.\n"
                    << endl
                    << "    -rand            specify the number of random child sequences that will be\n"
                    << "                     sub-sampled from the original dataset and used in a recombi\n"
                    << "                     -nation sensitivity analysis.  If the parent and the\n"
                    << "                     child datasets are taken from different files, this option\n"
                    << "                     will not affect the parent dataset.  In case the parent and\n"
                    << "                     child sequences are from the same file, only the sub-\n"
                    << "                     sampled sequences will be used as both the parent and the\n"
                    << "                     child datasets for the sensitivity analysis.\n"
                    << endl
                    << "    -n               specify the number of iterations of the recombination\n"
                    << "                     analysis.  For example, the following command:\n"
                    << endl
                    << "                          " << exeName << "-full seq_file -rand 20 -n 100\n"
                    << endl
                    << "                     will test the recombination sensitivity of the given\n"
                    << "                     dataset 100 times, each time considers 20 randomly sub-\n"
                    << "                     sampled sequences.\n";
            break;
        
        
        case GUI:
            //TODO: catch all command line options to shown on GUI.
            assert(false);
            break;
        
        
        case SILENT:
            // No message.
            break;
    }
}

void Interface::enableLog(const bool &isLogStreamEnabled, const bool &isLogFileEnabled) {
    this->isLogStreamEnabled = isLogStreamEnabled;
    this->isLogFileEnabled = isLogFileEnabled;
}

void Interface::enableWarning(const bool &isWarningEnabled) {
    this->isWarningEnabled = isWarningEnabled;
}




////////////////////////////////////////////////////////////////////////////////
//  COUNTER
////////////////////////////////////////////////////////////////////////////////

const int Interface::Counter::DECIMAL_PLACES = 3;

void Interface::Counter::initialize(string newProcessingText, double newMinValue,
                                    double newMaxValue, bool percentCounter) {
    
    processingText = newProcessingText;
    minValue = newMinValue;
    maxValue = newMaxValue;
    isPercentCounter = percentCounter;
    
    currentValue = minValue;
    startTime = time(nullptr);
}

const string Interface::Counter::showValue() const {
    stringstream tmpStream;
    
    if (isPercentCounter) {
        double percentDone = (currentValue - minValue) / (maxValue - minValue);
        percentDone *= 100.0;
        
        tmpStream.setf(ios::fixed);
        tmpStream.setf(ios::showpoint);
        tmpStream.precision(DECIMAL_PLACES);
        tmpStream.width(4 + DECIMAL_PLACES);
        tmpStream << percentDone;
        
    } else {
        tmpStream << currentValue;
    }
    
    return tmpStream.str();
}

const long Interface::Counter::getElapsedSeconds(void) const {
    time_t currentTime = time(nullptr);
    return currentTime - startTime;
}