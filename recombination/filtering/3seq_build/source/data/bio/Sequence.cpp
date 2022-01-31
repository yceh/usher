/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 * 
 * FILE NAME:   Sequence.cpp
 * CREATED ON:  27 July 2011, 16:46
 * AUTHOR:      Maciej F. Boni, Ha Minh Lam
 * 
 * DESCRIPTION: 
 * 
 * HISTORY:     Version     Date        Description
 *              1.0         2011-07-27  Created
 * 
 * VERSION:     1.0
 * LAST EDIT:   27 July 2011
 */

#include "Sequence.h"


//Sequence::Sequence(const string& dnaStr, const string& newName)
//: name(newName) {
//
//    try {
//        nucleotides.reserve(dnaStr.length());
//
//        for (unsigned long i = 0; i < dnaStr.length(); i++) {
//            Nucleotide * nucleotidePtr = new Nucleotide(dnaStr[i], i, this);
//            nucleotides.push_back(nucleotidePtr);
//        }
//
//    } catch (bad_alloc) {
//        Interface::instancePtr() << "Not enough memory to create new nucleotide object.\n";
//        Interface::instancePtr().showError(true, true);
//    } catch (length_error) {
//        Interface::instancePtr() << "Not enough memory to store the sequence " << name << "\n";
//        Interface::instancePtr().showError(true, true);
//    }
//}

//Sequence::Sequence(const string& newName)
//: name(newName) {
//}

//Sequence::~Sequence() {
//    for (unsigned long i = 0; i < nucleotides.size(); i++) {
//        if (nucleotides[i]) {
//            delete nucleotides[i];
//        }
//    }
//}

const string Sequence::toString(void) const {
    string sequenceStr = "";
    sequenceStr.reserve(nucleotides.size());
    
    for (auto nu = nucleotides.begin(), nuEnd = nucleotides.end(); nu != nuEnd; nu++) {
        sequenceStr += (*nu)->getCode();
    }
    
    return sequenceStr;
}

const unsigned long Sequence::ntDistanceTo(Sequence * const &another, bool const &ignoreGap) const {
    assert(another && this->getLength() == another->getLength());
    
    unsigned long distance = 0;
    
    Nucleotide * thisNu;
    Nucleotide * anotherNu;
    unsigned long seqLen = this->nucleotides.size();
    for (unsigned long i = 0; i < seqLen; i++) {
        thisNu = this->nucleotides[i];
        anotherNu = another->getNuAt(i);
        if (!ignoreGap || (thisNu->isNotGap() && anotherNu->isNotGap())) {
            if (thisNu->getCode() != anotherNu->getCode()) {
                distance++;
            }
        }
    }
    
    return distance;
}

bool Sequence::isSubsetOf(Sequence * const &another) const {
    assert(another && this->getLength() == another->getLength());
    
    Nucleotide * thisNu;
    Nucleotide * anotherNu;
    unsigned long seqLen = this->nucleotides.size();
    
    for (unsigned long i = 0; i < seqLen; i++) {
        thisNu = this->nucleotides[i];
        anotherNu = another->getNuAt(i);
        if (thisNu->getCode() != anotherNu->getCode() && thisNu->isNotGap()) {
            return false;
        }
    }
    
    return true;
}
