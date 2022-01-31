/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 *
 * FILE NAME:   config.h
 * CREATED ON:  15 June 2011
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 *
 * DESCRIPTION: This file contains the declaration of all default values of
 *              variables in the 3SEQ program.
 *
 * HISTORY:     Version    Date          Description
 *              1.0        2011-06-15    Created.
 *
 * VERSION:     1.0
 * LAST EDIT:   15 June 2011
 */

#ifndef CONFIG_H
#define CONFIG_H

#include <string>

namespace config {
    
    /* Default file names */
    const std::string DEFAULT_LOG_FILE_NAME = "3s.log";
    const std::string DEFAULT_RECOMBINANTS_FILE_NAME = "3s.rec";
    const std::string DEFAULT_LONG_RECOMBINANTS_FILE_NAME = "3s.longRec";
    const std::string DEFAULT_PVAL_HISTOGRAM_FILE_NAME = "3s.pvalHist";
    const std::string DEFAULT_SKIPPED_TRIPLETS_FILE_NAME = "3s.skipped";
    
    const std::string DEFAULT_MATCH_FILE_NAME = "3s.match"; // used in match-run
    const std::string DEFAULT_TRIPLET_PS_FILE_EXT = "triplet.eps"; // used in triplet-run
    
    
    /* Default analysis settings */
    const unsigned long DEFAULT_MIN_RECOMBINANT_LENGTH = 100;
    const double DEFAULT_REJECT_THRESHOLD = 0.05;
    
    
    /** The size of the histogram of P-values */
    const int PVAL_HISTOGRAM_SIZE = 41; // 0 to 40
    
    
    /** The directory where all the program data are stored */
#ifdef _WIN32
    const std::string DEFAULT_CONF_DIR = ".3seq";   // not tested yet
#elif __APPLE__
    const std::string DEFAULT_CONF_DIR = "Library/3seq";
#elif __linux__
    /* Linux */
    const std::string DEFAULT_CONF_DIR = ".config/3seq";
#else
    /* Unknown OS */
    const std::string DEFAULT_CONF_DIR = ".3seq";
#endif
    
    /**
     * The name of the file where all the program configurations are stored
     */
    const std::string DEFAULT_CONF_FILE = "3seq.conf";
    
    /**
     * The program will only ask for permission when it needs to allocate more than this amount
     * (in MB) of RAM.
     */
    const int MAX_AUTO_MEM_ALLOT = 100;
    
    /**
     * Whenever the program needs to allocate more than this amount of RAM (in MB), it will show a
     * warning telling the user to wait for memory allocation.
     */
    const int MEM_ALLOT_WARNING_SIZE = 100;
    
    /**
     * The rate (in seconds) at which the progress counter is updated.
     */
    const long PROGRESS_MONITOR_UPDATE_RATE = 2;
}

#endif	/* CONFIG_H */
