// ***************************************************************************
// SamHeader.h (c) 2010 Derek Barnett
// Marth Lab, Department of Biology, Boston College
// All rights reserved.
// ---------------------------------------------------------------------------
// Last modified: 18 April 2011 (DB)
// ---------------------------------------------------------------------------
// Provides direct read/write access to the SAM header data fields.
// ***************************************************************************

#ifndef SAM_HEADER_H
#define SAM_HEADER_H

#include "api_global.h"
#include "SamProgramChain.h"
#include "SamReadGroupDictionary.h"
#include "SamSequenceDictionary.h"
#include <string>
#include <vector>

namespace BamTools {

struct API_EXPORT SamHeader {

    // ctor & dtor
    SamHeader(const std::string& headerText = "");
    SamHeader(const SamHeader& other);
    ~SamHeader(void);

    // query/modify entire SamHeader
    void Clear(void);                                   // clears all header contents
    bool IsValid(bool verbose = false) const;           // returns true if SAM header is well-formed
    void SetHeaderText(const std::string& headerText);  // replaces data fields with contents of SAM-formatted text
    std::string ToString(void) const;                   // returns the printable, SAM-formatted header text

    // convenience query methods
    bool HasVersion(void) const;            // returns true if header contains format version entry
    bool HasSortOrder(void) const;          // returns true if header contains sort order entry
    bool HasGroupOrder(void) const;         // returns true if header contains group order entry
    bool HasSequences(void) const;          // returns true if header contains any sequence entries
    bool HasReadGroups(void) const;         // returns true if header contains any read group entries
    bool HasPrograms(void) const;           // returns true if header contains any program record entries
    bool HasComments(void) const;           // returns true if header contains comments

    // --------------
    // data members
    // --------------

    // header metadata (@HD line)
    std::string Version;                    // VN:<Version>  *Required for valid SAM header, if @HD record is present*
    std::string SortOrder;                  // SO:<SortOrder>
    std::string GroupOrder;                 // GO:<GroupOrder>

    // header sequences (@SQ entries)
    SamSequenceDictionary Sequences;

    // header read groups (@RG entries)
    SamReadGroupDictionary ReadGroups;

    // header program data (@PG entries)
    SamProgramChain Programs;

    // header comments (@CO entries)
    std::vector<std::string> Comments;
};

} // namespace BamTools

#endif // SAM_HEADER_H
