/**
 * Created by Ha Minh Lam on 2017-01-19.
 */


#include "AppInfo.h"
#include <sstream>

// /////////////////////////////////////////////////////////////////////////////
// Set static constants
// /////////////////////////////////////////////////////////////////////////////

const std::string AppInfo::APP_NAME = "3SEQ";
const std::string AppInfo::DESCRIPTION =
        "Software For Identifying Recombination In Sequence Data";

const std::string AppInfo::AUTHOR_NAME = "Maciej F. Boni, Ha Minh Lam";
const std::string AppInfo::AUTHOR_CONTACT =
        "maciej.boni@ndm.ox.ac.uk, mboni@oucru.org";

const int AppInfo::VERSION_MAJOR = 1;
const int AppInfo::VERSION_MINOR = 7;

/**
 * This date can be set to a fixed value of the following format:
 *     <i>"Jan 19 2017"</i>.    <br/>
 * If this string is empty, the build date will be read from the __DATE__ macro.
 */
const std::string AppInfo::BUILD_DATE = "Jun 12 2017";
//const std::string AppInfo::BUILD_DATE = "";

const int AppInfo::VERSION_DETAIL_LEVEL = 3;
const std::string AppInfo::DEVELOPMENT_STATUS = "beta 3";


const std::string AppInfo::COPYRIGHT = ""
        "Copyright © 2006-10 Maciej F. Boni. All Rights Reserved.\n"
        "Licensed for non-commercial use only.\n";

const std::string AppInfo::CITATION = ""
        "When using this software software, please cite\n"
        "\n"
        "    Lam HM, Ratmann O, Boni MF.  Improved space complexity for 3SEQ\n"
        "    recombination detection algorithm.  In preparation, 2017.\n"
        "\n"
        "When referring to a core part of the statistics used,you can cite\n"
        "\n"
        "    Boni MF, Posada D, Feldman MW. An exact nonparametric method for inferring\n"
        "    mosaic structure in sequence triplets.  Genetics, 176:1035-1047, 2007.\n"
        "\n"
        "If you make extensive use of the Hogan-Siegmund approximations in the results,\n"
        "please also cite\n"
        "\n"
        "    Hogan ML, Siegmund D.  Large deviations for the maxima of some random\n"
        "    fields.  Advances in Applied Mathematics, 7:2-22, 1986.\n";


// /////////////////////////////////////////////////////////////////////////////
// Singleton implementation
// /////////////////////////////////////////////////////////////////////////////

AppInfo * AppInfo::instancePtr = nullptr;

AppInfo * AppInfo::instance() {
    if (AppInfo::instancePtr == nullptr) {
        AppInfo::instancePtr = new AppInfo();
    }
    return AppInfo::instancePtr;
}



// /////////////////////////////////////////////////////////////////////////////
// Definition of private methods
// /////////////////////////////////////////////////////////////////////////////

AppInfo::AppInfo(void) {
    appName = AppInfo::APP_NAME;
    description = AppInfo::DESCRIPTION;
    authorName = AppInfo::AUTHOR_NAME;
    authorContact = AppInfo::AUTHOR_CONTACT;
    
    copyright = AppInfo::COPYRIGHT;
    citation = AppInfo::CITATION;
    
    /* Example of __DATE__ string: "Jan 19 2017" */
    buildDate = AppInfo::BUILD_DATE;
    if (AppInfo::BUILD_DATE.length() <= 0) {
        buildDate = __DATE__;
    }
    
    /* Example of __TIME__ string: "21:06:19" */
    buildTime = __TIME__;
    
    /* Presume that the executable path is the same as the application name.
     * This should be changed after the program starts. */
    executablePath = appName;
    
    /* Initialize version strings */
    std::stringstream shortVersionStream;
    shortVersionStream << AppInfo::VERSION_MAJOR << "."
                       << AppInfo::VERSION_MINOR;
    if (AppInfo::DEVELOPMENT_STATUS.length() > 0) {
        shortVersionStream << " " << AppInfo::DEVELOPMENT_STATUS;
    }
    shortVersionString = shortVersionStream.str();
    
    std::stringstream fullVersionStream;
    fullVersionStream << shortVersionString;
    if (AppInfo::VERSION_DETAIL_LEVEL > 0) {
        fullVersionStream << " " << getBuildDetails();
    }
    fullVersionString = fullVersionStream.str();
}


const std::string AppInfo::getBuildDetails() const {
    std::string buildDetails = "";
    
    if (AppInfo::VERSION_DETAIL_LEVEL > 0) {
        buildDetails += "build " + getBuildYear();
    }
    if (AppInfo::VERSION_DETAIL_LEVEL > 1) {
        buildDetails += getBuildMonth();
    }
    if (AppInfo::VERSION_DETAIL_LEVEL > 2) {
        buildDetails += getBuildDay();
    }
    if (AppInfo::VERSION_DETAIL_LEVEL > 3) {
        buildDetails += "." + getBuildHour();
    }
    if (AppInfo::VERSION_DETAIL_LEVEL > 4) {
        buildDetails += getBuildMinute();
    }
    if (AppInfo::VERSION_DETAIL_LEVEL > 5) {
        buildDetails += getBuildSecond();
    }
    
    return buildDetails;
}


// /////////////////////////////////////////////////////////////////////////////
// Definition of public methods
// /////////////////////////////////////////////////////////////////////////////

const std::string &AppInfo::getApplicationName(void) const {
    return appName;
}

const std::string &AppInfo::getDescription(void) const {
    return description;
}

const std::string &AppInfo::getAuthorName(void) const {
    return authorName;
}

const std::string &AppInfo::getAuthorContact(void) const {
    return authorContact;
}

const std::string &AppInfo::getCopyright(void) const {
    return copyright;
}

const std::string &AppInfo::getCitation(void) const {
    return citation;
}

const std::string &AppInfo::getVersionShortString(void) const {
    return shortVersionString;
}

const std::string &AppInfo::getVersionFullString(void) const {
    return fullVersionString;
}

const std::string &AppInfo::getExecutablePath() const {
    return executablePath;
}

const std::string &AppInfo::getBuildDate() const {
    return buildDate;
}

const std::string &AppInfo::getBuildTime() const {
    return buildTime;
}

std::string AppInfo::getBuildYear() const {
    std::string buildYear = getBuildDate().substr(9, 2);
    return buildYear;
}

std::string AppInfo::getBuildMonth() const {
    std::string buildMonth = getBuildDate().substr(0, 3);
    buildMonth =
            buildMonth == "Jan" ? "01" :
            buildMonth == "Feb" ? "02" :
            buildMonth == "Mar" ? "03" :
            buildMonth == "Apr" ? "04" :
            buildMonth == "May" ? "05" :
            buildMonth == "Jun" ? "06" :
            buildMonth == "Jul" ? "07" :
            buildMonth == "Aug" ? "08" :
            buildMonth == "Sep" ? "09" :
            buildMonth == "Oct" ? "10" :
            buildMonth == "Nov" ? "11" :
            buildMonth == "Dec" ? "12" : "??";
    return buildMonth;
}

std::string AppInfo::getBuildDay() const {
    std::string buildDay = getBuildDate().substr(4, 2);
    if (buildDay[0] == ' ') {
        buildDay[0] = '0';
    }
    return buildDay;
}

std::string AppInfo::getBuildHour() const {
    std::string buildHour = getBuildTime().substr(0, 2);
    return buildHour;
}

std::string AppInfo::getBuildMinute() const {
    std::string buildMinute = getBuildTime().substr(3, 2);
    return buildMinute;
}

std::string AppInfo::getBuildSecond() const {
    std::string buildSecond = getBuildTime().substr(6, 2);
    return buildSecond;
}

void AppInfo::setExecutablePath(const std::string &executablePath) {
    if (executablePath.length() > 0 && executablePath[0] != '\0') {
        this->executablePath = executablePath;
    }
}

void AppInfo::setExecutablePath(const char * executablePath) {
    if (executablePath != nullptr && executablePath[0] != '\0') {
        setExecutablePath(std::string(executablePath));
    }
}
