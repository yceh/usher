/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   FullRun.cpp
 * CREATED ON:  17 August 2011
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 *
 * DESCRIPTION: See "FullRun.h"
 *
 * HISTORY:     Version     Date        Description
 *              1.0         2011-08-17  created
 *
 * VERSION:     1.0
 * LAST EDIT:   17 August 2011
 */

#include "FullRun.h"
#include "../data/bio/Triplet.h"

FullRun::FullRun(int argc, char ** argv) : Run(argc, argv) {
    Interface::instance().startProgram("Full Run");
    
    fileSkippedTriplets = nullptr;
    fileRecombinants = nullptr;
    fileLongRecombinants = nullptr;
}

FullRun::~FullRun() {
    if (fileSkippedTriplets != nullptr) delete fileSkippedTriplets;
    if (fileRecombinants != nullptr) delete fileRecombinants;
    if (fileLongRecombinants != nullptr) delete fileLongRecombinants;
}

void FullRun::parseCmdLine() {
    int argNum = getRunArgsNum();
    
    if (argNum < 1) {
        Interface::instance() << "Not enough parameter for full run.\n";
        Interface::instance().showError(true, true);
    }
    
    string filePath = argVector[2];
    SequenceFile parentFile(filePath, SequenceFile::UNKNOWN);
    if (!parentFile.exists()) {
        Interface::instance() << "File \"" << filePath << "\" not found.\n";
        Interface::instance().showError(true, true);
    }
    parentFile.detectFileType();
    parentDataset = new GenomeSet();
    parentFile.read(parentDataset, "");
    
    if (argNum >= 2) {
        argNum = 2;
        filePath = argVector[3];
        SequenceFile childFile(filePath, SequenceFile::UNKNOWN);
        if (!childFile.exists()) {
            Interface::instance() << "File \"" << filePath << "\" not found.\n";
            Interface::instance().showError(true, true);
        }
        childFile.detectFileType();
        childDataset = new GenomeSet();
        childFile.read(childDataset, "");
    }
    
    isRandomMode = false;
    randomSeqNum = 0;
    randomLoopNum = 1;
    pTableFile = nullptr;
    
    for (int i = 2 + argNum; i < argCount - 1; i++) {
        string arg = argVector[i];
        if (arg == "-ptable" || arg == "-p") {
            string pTableFilePath = argVector[i + 1];
            pTableFile = new PTableFile(pTableFilePath);
            argVector[i] = "";
            argVector[i + 1] = "";
            
        } else if (arg == "-rand") {
            errno = 0;
            randomSeqNum = strtoul(argVector[i + 1].c_str(), nullptr, 0);
            if (errno || randomSeqNum <= 0) {
                Interface::instance() << "Cannot convert \"" << argVector[i + 1]
                                      << "\" to a positive number\n";
                Interface::instance().showError(true, true);
            }
            argVector[i] = "";
            argVector[i + 1] = "";
            isRandomMode = true;
            
        } else if (arg == "-n") {
            errno = 0;
            randomLoopNum = strtoul(argVector[i + 1].c_str(), nullptr, 0);
            if (errno || randomLoopNum <= 0) {
                Interface::instance() << "Cannot convert \"" << argVector[i + 1]
                                      << "\" to a positive number\n";
                Interface::instance().showError(true, true);
            }
            argVector[i] = "";
            argVector[i + 1] = "";
            isRandomMode = true;
        }
    }
    
    if (pTableFile != nullptr && !pTableFile->exists()) {
        Interface::instance() << "The file \"" << pTableFile->getFilePath()
                              << "\" does not exist.\n";
        Interface::instance().showError(true, true);
    }
    
    deleteCmdLineArgs(argNum + 1); // --fullrun inputFile ...
    Run::parseCmdLine(); // Get common run options
}

void FullRun::perform() {
    Run::perform(); // process common data
    
    verifyData();
    setup();
    loadPTable(pTableFile);
    
    
    Interface::instance() << Interface::SEPARATOR << endl;
    Interface::instance().showLog(true);
    
    progress();
    
    Interface::instance() << Interface::SEPARATOR << endl;
    Interface::instance().showLog(true);
    
    if (!cloAllBpCalculated && !cloNoBpCalculated && !cloNo3sRecFile) {
        /* We calculate breakpoints for the best triplets */
        for (unsigned long i = 0; i < childDataset->getSize(); i++) {
            Genome * child = childDataset->getGenome(i);
            Triplet * bestRecTriplet = child->getBestRecombinantTriplet();
            if (bestRecTriplet != nullptr) {
                recordRecombinantTriplet(bestRecTriplet);
            }
        }
    }
    
    
    if (isRandomMode && randomLoopNum > 1) {
        Interface::instance() << "Recombination sensitivity: "
                              << recombinationSensitivity << "/"
                              << randomLoopNum << "\n";
        Interface::instance().showOutput(true);
        
    } else {
        displayResult();
        Interface::instance() << Interface::SEPARATOR << endl;
        Interface::instance().showLog(true);
        savePValHistogram();
    }
    
    
    /* Close all files */
    if (fileSkippedTriplets)
        fileSkippedTriplets->close();
    if (fileRecombinants)
        fileRecombinants->close();
    if (fileLongRecombinants != nullptr)
        fileLongRecombinants->close();
}

void FullRun::verifyData(void) {
    if (isRandomMode) {
        cloStopAtFirstRecombinant = true;
        
        if (cloBeginSequence > 0 || cloEndSequence < parentDataset->getSize()) {
            Interface::instance()
                    << "-b and -e options cannot be used with random sub-sampling.\n";
            Interface::instance().showError(true, true);
        }
        
        if (isChildParentInDiffFiles) {
            Interface::instance()
                    << "Parent and child sequences are from different datasets.\n"
                    << "Only the child dataset will be sub-sampled.\n";
            Interface::instance().showWarning(true);
        }
        
        if (cloIsAllTriplesUsed) {
            Interface::instance() << ""
                    "Auto-enable -# option for random sub-sampling. Actual number of\n"
                    "comparisons will be used in corrections.";
            Interface::instance().showWarning(true);
            cloIsAllTriplesUsed = false;
        }
        
        childDataset->randomlyActivate(randomSeqNum);
        calculateTripletNum();
    }
    
    /* Calculate the number of actual comparisons */
    double totalTripletNum = getNumTotalTriplets();
    
    if (totalTripletNum <= 0) {
        Interface::instance()
                << "The number of sequences is not enough (at least 3 sequences needed).\n";
        Interface::instance().showError(true, true);
    } else if (totalTripletNum > MAX_UL) {
        Interface::instance()
                << "The number of comparisons you are requesting is larger than LONG_MAX for an\n"
                << "\"unsigned long\" type on your system. The P-values below should be correct,\n"
                << "but the numbers of comparisons that were approximated/computed/skipped might\n"
                << "not be the exact values.";
        Interface::instance().showWarning();
    }
    
    /* Show some logs */
    Interface::instance()
            << "Using " << parentDataset->getNumActiveSeq() << " sequences as parents.\n"
            << "Using " << childDataset->getNumActiveSeq() << " sequences as children.\n";
    Interface::instance().showLog(true);
    
    Interface::instance() << "Need a p-value of "
                          << cloRejectThreshold / getNumTripletsForCorrection()
                          << " to survive multiple comparisons correction.\n";
    Interface::instance().showLog(true);
}

void FullRun::setup(void) {
    if (!isRandomMode) {
        randomLoopNum = 1;
    } else {
        cloNoBpCalculated = true;
        cloNo3sRecFile = true;
    }
    
    /* Prepare file names */
    string recombinantsFileName = config::DEFAULT_RECOMBINANTS_FILE_NAME;
    string longRecombinantsFileName = config::DEFAULT_LONG_RECOMBINANTS_FILE_NAME;
    string skippedTripletFileName = config::DEFAULT_SKIPPED_TRIPLETS_FILE_NAME;
    if (id.length() > 0) {
        const string period = ".";
        recombinantsFileName = id + period + recombinantsFileName;
        longRecombinantsFileName = id + period + longRecombinantsFileName;
        skippedTripletFileName = id + period + skippedTripletFileName;
    }
    
    /* 3s.skipped */
    fileSkippedTriplets = nullptr;
    if (cloIsSkippedTriplesRecorded) {
        fileSkippedTriplets = new TextFile(skippedTripletFileName);
        
        Interface::instance()
                << "Skipped triplets will be recorded to the file \""
                << fileSkippedTriplets->getPath() << "\".\n";
        
        if (fileSkippedTriplets->exists()) {
            Interface::instance()
                    << "This file already exists. Do you want to overwrite it?";
            if (!Interface::instance().showWarning(false)) {
                /* User says "No" */
                delete fileSkippedTriplets;
                fileSkippedTriplets = nullptr;
                Interface::instance()
                        << "Skipped triplets will not be saved.\n";
                Interface::instance().showLog(true);
            } else {
                fileSkippedTriplets->removeFile();
            }
            
        } else {
            Interface::instance().showLog(true);
        }
    }
    
    /* 3s.rec */
    fileRecombinants = nullptr;
    if (!cloNo3sRecFile) {
        fileRecombinants = new TextFile(recombinantsFileName);
        
        Interface::instance()
                << "All recombinant triplets will be recorded to the file \""
                << fileRecombinants->getPath() << "\".\n";
        
        if (fileRecombinants->exists()) {
            Interface::instance()
                    << "This file already exists. Do you want to overwrite it?";
            if (!Interface::instance().showWarning(false)) {
                /* User says "No" */
                delete fileRecombinants;
                fileRecombinants = nullptr;
                Interface::instance()
                        << "Recombinant triplets will not be saved.\n";
                Interface::instance().showLog(true);
            } else {
                fileRecombinants->removeFile();
            }
            
        } else {
            Interface::instance().showLog(true);
        }
        
        if (fileRecombinants != nullptr) {
            /** Write the header */
            fileRecombinants->writeLine(
                    "P_ACCNUM\tQ_ACCNUM\tC_ACCNUM\tm\tn\tk\t"
                            "p\tHS?\tlog(p)\tDS(p)\tDS(p)\tmin_rec_length\tbreakpoints");
        }
    }
    
    /* 3s.longRec */
    fileLongRecombinants = new TextFile(longRecombinantsFileName);
    Interface::instance()
            << "Long recombinants will be recorded to the file \""
            << fileLongRecombinants->getPath() << "\".\n";
    
    if (fileLongRecombinants->exists()) {
        Interface::instance()
                << "This file already exists. Do you want to overwrite it?";
        if (!Interface::instance().showWarning(false)) {
            /* User says "No" */
            delete fileLongRecombinants;
            fileLongRecombinants = nullptr;
            Interface::instance() << "Long Recombinants will not be saved.\n";
            Interface::instance().showLog(true);
        } else {
            fileLongRecombinants->removeFile();
        }
        
    } else {
        Interface::instance().showLog(true);
    }
}

void FullRun::showProgress(double currentLoop, bool isFinish) const {
    char strBuf[200];
    
    if (cloAllBpCalculated) {
        sprintf(strBuf,
                "     %s       %e       %10.0lf            %lu",
                Interface::instance().getElapsedTime().c_str(),
                minPVal, numRecombinantTriplets,
                longestRecombinantSegment);
    } else if (cloStopAtFirstRecombinant) {
        std::string recFound = "NO";
        if (numRecombinantTriplets > 0) {
            recFound = "YES";
        }
        sprintf(strBuf,
                "     %s       %e       %6s",
                Interface::instance().getElapsedTime().c_str(),
                minPVal, recFound.c_str());
    } else {
        sprintf(strBuf,
                "     %s       %e       %10.0lf",
                Interface::instance().getElapsedTime().c_str(),
                minPVal, numRecombinantTriplets);
    }
    Interface::instance() << strBuf;
    
    if (isFinish) {
        Interface::instance().finishCounting();
    } else {
        Interface::instance().count(currentLoop);
    }
}

void FullRun::progress(void) {
    if (cloAllBpCalculated) {
        Interface::instance()
                << "                                                Recombinant      Longest Recombinant\n"
                << "Progress    Time Elapsed    Minimum P-Value    Triplets Found          Segment\n";
    } else {
        Interface::instance()
                << "                                                Recombinant\n"
                << "Progress    Time Elapsed    Minimum P-Value    Triplets Found\n";
    }
    Interface::instance().showLog(false);
    
    recombinationSensitivity = 0;
    
    for (unsigned long randLoopCounter = 0;
         randLoopCounter < randomLoopNum; randLoopCounter++) {
        
        if (isRandomMode) {
            childDataset->randomlyActivate(randomSeqNum);
            calculateTripletNum();
        }
        
        numRecombinantTriplets = 0.0;
        numComputedExactly = 0.0;
        numApproximated = 0.0;
        numSkipped = getNumTotalTriplets();
        
        minPVal = 1.0;
        longestRecombinantSegment = 0;
        
        /* Initialise progress counter */
        unsigned long activeChildNum = childDataset->getNumActiveSeq();
        unsigned long activeParentNum = parentDataset->getNumActiveSeq();
        
        double totalLoop = static_cast<double> (activeChildNum)
                           * static_cast<double> (activeParentNum);
        Interface::instance().initCounter("", 0, totalLoop);
        
        
        /* This flag indicates if the progressing should be stopped */
        bool isStoped = false;
        
        /* Progressing begins */
        unsigned long activeChildCounter = 0;
        unsigned long activeDadCounter = 0;
        time_t lastTime = time(nullptr) - config::PROGRESS_MONITOR_UPDATE_RATE - 1;
        
        Genome * child = nullptr;
        Genome * dad = nullptr;
        Genome * mum = nullptr;
        Triplet * triplet = nullptr;
        
        for (unsigned long childIndex = 0;
             !isStoped && childIndex < childDataset->getSize(); childIndex++) {
            child = childDataset->getGenome(childIndex);
            if (child->isActive()) {
                activeChildCounter++;
            } else {
                continue;
            }
            
            activeDadCounter = 0;
            for (unsigned long dadIndex = 0;
                 !isStoped && dadIndex < parentDataset->getSize(); dadIndex++) {
                
                dad = parentDataset->getGenome(dadIndex);
                if (dad->isActive()) {
                    activeDadCounter++;
                }
                if (!dad->isActive() || dad == child) {
                    continue;
                }
                
                /* Show progress counter */
                time_t currentTime = time(nullptr);
                if (currentTime - lastTime >=
                    config::PROGRESS_MONITOR_UPDATE_RATE) {
                    double currentLoop =
                            static_cast<double> (activeChildCounter - 1)
                            * static_cast<double> (activeParentNum)
                            + static_cast<double> (activeDadCounter - 1);
                    showProgress(currentLoop, false);
                    lastTime = currentTime;
                }
                
                for (unsigned long mumIndex = 0;
                     !isStoped && mumIndex < parentDataset->getSize();
                     mumIndex++) {
                    
                    mum = parentDataset->getGenome(mumIndex);
                    if (!mum->isActive() || mum == child || mum == dad)
                        continue;
                    
                    triplet = new Triplet(dad, mum, child);
                    
                    /* Let's see if we can calculate the exact P-value */
                    bool hasPValue = triplet->calculatePVal(
                            !cloNoSiegmundApprox);
                    if (!hasPValue) {
                        if (fileSkippedTriplets) {
                            /* No need to open file, writeLine() will do that automatically */
                            fileSkippedTriplets->writeLine(triplet->toString());
                        }
                        delete triplet;
                        continue;
                    }
                    numSkipped -= 1.0;
                    
                    if (triplet->isPValApproximated()) {
                        numApproximated += 1.0;
                    } else {
                        numComputedExactly += 1.0;
                    }
                    
                    double pValue = triplet->getPVal();
                    addPValIntoHistogram(pValue);
                    
                    if (pValue < minPVal)
                        minPVal = pValue;
                    
                    /* Let's see if the P-value is significant */
                    if (stat::correction::dunnSidak(pValue, getNumTripletsForCorrection())
                        < cloRejectThreshold) {
                        
                        numRecombinantTriplets += 1.0;
                        if (child->getRecombinantType() == Genome::Na_REC) {
                            child->setRecombinantType(Genome::SHORT_REC);
                        }
                        
                        /* Go into this loop either if we're doing ALL breakpoint or
                         * NO breakpoint if we are doing no breakpoint, then the 
                         * breakpoint calculations is skipped. */
                        if (cloAllBpCalculated || cloNoBpCalculated) {
                            recordRecombinantTriplet(triplet);
                        }
                        
                        /* Now check if this is the *best* recombinant (meaning lowest p-value)
                         * and if it is, record it */
                        if (child->getBestRecombinantTriplet() == nullptr
                            || pValue <
                               child->getBestRecombinantTriplet()->getPVal()) {
                            child->setBestRecombinantTriplet(triplet);
                        } else {
                            delete triplet;
                        }
                        
                        if (cloStopAtFirstRecombinant) {
                            /* Finish progressing early */
                            isStoped = true;
                        }
                        
                    } else {
                        delete triplet;
                    }
                }
            }
        }
        
        /* Progressing finished */
        showProgress(0.0, true);
        
        if (numRecombinantTriplets > 0) {
            recombinationSensitivity++;
        }
    }
    
    Interface::instance().showLog(true);
}

void FullRun::recordRecombinantTriplet(Triplet * triplet) {
    if (!cloNoBpCalculated) {
        if (triplet->getMinRecLength() >= cloMinLongRecombinantLength)
            triplet->getChild()->setRecombinantType(Genome::LONG_REC);
        
        if (triplet->getMinRecLength() > longestRecombinantSegment)
            longestRecombinantSegment = triplet->getMinRecLength();
    }
    
    /* Begin writing 3s.rec file */
    if (fileRecombinants == nullptr) return;
    
    double pValue = triplet->getPVal();
    long double dunnSidak = stat::correction::dunnSidak(pValue,
                                                        getNumTripletsForCorrection());
    
    char buff[1000];
    sprintf(buff,
            "%s\t%s\t%s\t%lu\t%lu\t%lu\t%1.12f\t%d\t%2.4f\t%1.5Lf\t%Le",
            triplet->getDad()->getAccession().c_str(),
            triplet->getMum()->getAccession().c_str(),
            triplet->getChild()->getAccession().c_str(),
            triplet->getUpStep(), triplet->getDownStep(),
            triplet->getMaxDecent(),
            pValue, triplet->isPValApproximated(), log10(pValue),
            dunnSidak, dunnSidak);
    fileRecombinants->write(buff);
    
    if (!cloNoBpCalculated) {
        sprintf(buff, "\t%lu", triplet->getMinRecLength());
        fileRecombinants->write(buff);
        
        vector<Triplet::BreakPointPair> breakPointPairs = triplet->getBreakPointPairs();
        
        for (unsigned long i = 0; i < breakPointPairs.size(); i++) {
            Triplet::BreakPointPair bpPair = breakPointPairs[i];
            fileRecombinants->write("\t " + bpPair.toString());
        }
    }
    
    fileRecombinants->writeLine();
    /* Finish writing 3s.rec file */
}

void FullRun::displayResult(void) {
    unsigned long numRecombinant = 0;
    unsigned long numLongRec = 0;
    
    
    for (unsigned long i = 0; i < childDataset->getSize(); i++) {
        Genome * child = childDataset->getGenome(i);
        
        if (child->getRecombinantType() != Genome::Na_REC) {
            numRecombinant++;
        }
        
        if (child->getRecombinantType() == Genome::LONG_REC) {
            numLongRec++;
            if (fileLongRecombinants != nullptr) {
                fileLongRecombinants->writeLine(child->getAccession());
            }
        }
    }
    
    
    Interface::instance()
            << "Number of triples tested :              "
            << getNumTotalTriplets() << "\n"
            << "Number of p-values computed exactly :   " << numComputedExactly
            << "\n"
            << "Number of p-values approximated (HS) :  " << numApproximated
            << "\n"
            << "Number of p-values not computed :       " << numSkipped << "\n"
            << endl
            << "Number of recombinant triplets :                               \t"
            << numRecombinantTriplets << "\n"
            << "Number of distinct recombinant sequences :                     \t"
            << numRecombinant << "\n";
    
    if (cloStopAtFirstRecombinant) {
        Interface::instance() << "[[ search stopped when first one found ]]\n"
                              << endl;
    }
    
    if (!cloNoBpCalculated) {
        Interface::instance()
                << "Number of distinct recombinant sequences at least "
                << cloMinLongRecombinantLength << "nt long : \t" << numLongRec
                << "\n"
                << "Longest of short recombinant segments :                        \t"
                << longestRecombinantSegment << "nt\n";
    }
    Interface::instance().showOutput(true);
    
    Interface::instance() << Interface::SEPARATOR << endl;
    Interface::instance().showLog(true);
    
    char formatedPVal[20];
    sprintf(formatedPVal,
            "%1.3e",
            static_cast<double> (stat::correction::dunnSidak(minPVal,
                                                             getNumTripletsForCorrection()))
    );
    Interface::instance()
            << "Rejection of the null hypothesis of clonal evolution at p = "
            << stat::correction::dunnSidak(minPVal,
                                           getNumTripletsForCorrection())
            << "\n"
            << "                                                        p = "
            << formatedPVal << "\n"
            << "                                            Uncorrected p = "
            << minPVal << "\n"
            << "                                            Bonferroni  p = "
            << stat::correction::bonferroni(minPVal,
                                            getNumTripletsForCorrection())
            << "\n";
    Interface::instance().showOutput(true);
}


//void FullRun::outputForPython(void) const {
//    printf("%1.4f\t%d\t%d\t%d\t%d\t%d",
//            correction::dunnSidak(minPVal, numForCorrection),
//            numComputedExactly,
//            numApproximated,
//            numSkipped,
//            parentDataset->getSize(),
//            numRecombinantTriplets
//            );
//}
