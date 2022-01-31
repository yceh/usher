/**
 * Copyright (c) ___. All Rights Reserved.
 *
 * FILE NAME:   util.h
 * CREATED ON:  26 September 2011, 15:28
 * AUTHOR:      lamhm
 *
 * DESCRIPTION: 
 *
 * HISTORY:     Version     Date            Description
 *              1.0         26 September 2011      Created
 *
 * VERSION:     1.0
 * LAST EDIT:   26 September 2011
 */

#ifndef UTIL_H
#define    UTIL_H

#include <string>
#include "../config.h"

namespace Util {


    /**
     * Get a key from standard input without showing anything on the output stream.
     * @return  The code of the key.
     */
    int getKey(void);

    const std::string getCurrentDir(void);

    bool saveConfig(std::string key, std::string value);

    /**
     * Get the configuration value of the given key
     * @param key The configuration key
     * @return If a value is found, return it. If not, return "\0".
     * @notes If the given key has more than one value, only the first value (found from the
     *        configuration file) is returned.
     */
    std::string getConfig(std::string key);

};

#endif	/* UTIL_H */

