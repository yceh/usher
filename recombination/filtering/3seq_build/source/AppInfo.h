/**
 * Created by Ha Minh Lam on 2017-01-19.
 */

#ifndef INC_APPINFO_H
#define INC_APPINFO_H


/* Include string library for the implementation of the AppInfo class */
#include <string>


/**
 * This is a singleton to manage all the basic information of the program.
 */
class AppInfo {
public:
    static AppInfo * instance();
    
    /**
     * Retrieve the name of the application (this may not be the name of the
     * executable file).
     * @return The application name.
     */
    const std::string &getApplicationName(void) const;
    
    /**
     * Retrieve the description of the application.
     * @return The application description.
     */
    const std::string &getDescription(void) const;
    
    /**
     * Retrieve the name(s) of the author(s).
     * @return Author(s)' name(s).
     */
    const std::string &getAuthorName(void) const;
    
    /**
     * Retrieve the contact(s) of the author(s).
     * @return The author(s)' contact(s).
     */
    const std::string &getAuthorContact(void) const;
    
    /**
     * Retrieve the copyright of the application.
     * @return The copyright.
     */
    const std::string &getCopyright(void) const;
    
    /**
     * Retrieve the citation text for the application.
     * @return The citation text.
     */
    const std::string &getCitation(void) const;
    
    /**
     * Retrieve the version of the application in short format.
     * @return The version of the application in short format.
     */
    const std::string &getVersionShortString(void) const;
    
    /**
     * Retrieve the version of the application with all the build details.
     * @return The version of the application in full format.
     */
    const std::string &getVersionFullString(void) const;
    
    /**
     * Retrieve the path (including the file name) to the executable file.
     * @return The path to the executable file.
     */
    const std::string &getExecutablePath() const;
    
    /**
     * Set the path (including the file name) to the executable file.
     * @param executablePath
     */
    void setExecutablePath(const std::string &executablePath);
    
    void setExecutablePath(const char * executablePath);
    
    const std::string & getBuildDate() const;
    
    const std::string &getBuildTime() const;
    
    std::string getBuildYear() const;
    
    std::string getBuildMonth() const;
    
    std::string getBuildDay() const;
    
    std::string getBuildHour() const;
    
    std::string getBuildMinute() const;
    
    std::string getBuildSecond() const;

private:
    AppInfo(void);
    
    const std::string getBuildDetails() const;
    
    static AppInfo * instancePtr;
    
    std::string appName;
    std::string executablePath;
    std::string description;
    std::string authorName;
    std::string authorContact;
    
    std::string copyright;
    std::string citation;

    std::string shortVersionString;
    std::string fullVersionString;
    
    std::string buildDate;
    std::string buildTime;
    
    
    static const std::string APP_NAME;
    static const std::string DESCRIPTION;
    static const std::string AUTHOR_NAME;
    static const std::string AUTHOR_CONTACT;
    
    static const std::string COPYRIGHT;
    static const std::string CITATION;
    
    
    static const int VERSION_MAJOR;
    static const int VERSION_MINOR;
    static const std::string BUILD_DATE;
    
    /**
     * The status of the current version (e.g. beta, alpha). If the development is completed,
     * this string should be empty ("").
     */
    static const std::string DEVELOPMENT_STATUS;
    
    /**
     * Specify whether the full version string includes the build time.
     * The value means:
     *      0 : Build time is not included at all;
     *      1 : Include build year ;
     *      2 : Include build year and month;
     *      3 : Include build year, month, and day;
     *      4 : Include build year, month, day, and hour;
     *      5 : Include build year, month, day, hour, and minute;
     *      over 5 : Include build year, month, day, hour, minute, and second;
     */
    static const int VERSION_DETAIL_LEVEL;
    
    
};

#endif //INC_APPINFO_H
