/**
 * FILE NAME:   util.cpp
 * CREATED ON:  26 September 2011, 15:28
 * AUTHOR:      lamhm
 *
 * DESCRIPTION: See "util.h"
 * 
 * HISTORY:     Version     Date            Description
 *              1.0         26 September 2011      Created
 *
 * VERSION:     1.0
 * LAST EDIT:   26 September 2011
 */

#include "Util.h"

#include <termios.h>
#include <unistd.h>
#include <sys/stat.h>
//#include <cstdlib>

#include <vector>
#include "../data/file/TextFile.h"

namespace Util {

    /// DECLARATION OF PRIVATE METHODS
    namespace {
        const std::string getConfigDir(void);
    }

    int getKey(void) {
        static const long MAX_CODE_LENGTH = 3;
        static const int MAX_CONTROL_CODE = 31;

        struct termios oldSettings;
        struct termios newSettings;
        char keycodes[MAX_CODE_LENGTH];

        int codeLength;
        int keyCode;

        fflush(stdout);
        fflush(stderr);

        tcgetattr(STDIN_FILENO, &oldSettings);
        newSettings = oldSettings;

        newSettings.c_cc[VTIME] = 1;
        newSettings.c_cc[VMIN] = MAX_CODE_LENGTH;
        newSettings.c_iflag &= ~(IXOFF);
        newSettings.c_lflag &= ~(ECHO | ICANON);
        tcsetattr(STDIN_FILENO, TCSANOW, &newSettings);

        codeLength = read(fileno(stdin), (void *) keycodes, MAX_CODE_LENGTH);
        keyCode = static_cast<int> (keycodes[0]);
        if (keyCode <= MAX_CONTROL_CODE && codeLength > 1) {
            keyCode = -(static_cast<int> (keycodes[codeLength - 1]));
        }

        tcsetattr(STDIN_FILENO, TCSANOW, &oldSettings);

        return keyCode;
    }

    const std::string getCurrentDir(void) {
        char stackBuffer[10240];
        if (getcwd(stackBuffer, sizeof(stackBuffer)) == nullptr) {
            throw "Cannot get the current directory.";
        }
        std::string currentDir = stackBuffer;

        return currentDir;
    }

    bool saveConfig(std::string key, std::string value) {
        try {
            unsigned long keyLen = key.length();
            string confDir = getConfigDir();
            std::vector<std::string> oldConfigList;

            TextFile configFile = TextFile(confDir + "/" + config::DEFAULT_CONF_FILE);

            if (configFile.exists()) {
                configFile.openToRead();
                oldConfigList = configFile.readAllLines();
                configFile.close();
            }

            configFile.openToWrite();

            for (unsigned long i = 0; i < oldConfigList.size(); i++) {
                std::string oldConf = oldConfigList[i];

                if (oldConf.length() > keyLen && oldConf.substr(0, keyLen) == key && oldConf[keyLen] == '=') {
                    continue;
                } else {
                    configFile.writeLine(oldConf);
                }
            }

            configFile.writeLine(key + "=" + value); // write new config
            configFile.close();

            return true;

        } catch (...) {
            return false;
        }
    }

    std::string getConfig(std::string key) {
        unsigned long keyLen = key.length();
        string confDir = getConfigDir();

        TextFile configFile = TextFile(confDir + "/" + config::DEFAULT_CONF_FILE);

        if (!configFile.exists()) {
            /* No configuration found */
            // throw "No config file.";
            return std::string("\0");
        }

        configFile.openToRead();
        std::vector<std::string> configList = configFile.readAllLines();
        configFile.close();

        for (unsigned long i = 0; i < configList.size(); i++) {
            std::string conf = configList[i];

            if (conf.length() > keyLen && conf.substr(0, keyLen) == key && conf[keyLen] == '=') {
                return conf.substr(keyLen + 1);
            }
        }

        /* No configuration found */
        return std::string("\0");
    }

}   /* Util */


/// DEFINITION OF PRIVATE METHODS
namespace Util {
    namespace {
        const std::string getConfigDir(void) {
            /* Get the user's home directory */
            char *homeDir = getenv("HOME");
            if (homeDir == nullptr) {
                throw "Cannot detect user's home directory.";
            }

            std::string confDir = homeDir;
            confDir += "/" + config::DEFAULT_CONF_DIR;

            /* Create a folder (if it does not exist) to store all the configurations of 3SEQ */
            struct stat st;
            if (stat(confDir.c_str(), &st) != 0) {
                if (mkdir(confDir.c_str(), S_IRWXU | S_IRWXG) != 0) {
                    throw "Fail to create config directory.";
                }
            }

            return confDir;
        }
    }
}   /* Util */
