/**
 * UGENE - Integrated Bioinformatics Tools.
 * Copyright (C) 2008-2026 UniPro <ugene@unipro.ru>
 * http://ugene.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */

#include "HostMethylation.h"

#include <cstring>

#include <QHash>

namespace U2 {

namespace {

struct BlockRule {
    const char* enzyme;
    const char* pattern;
    int siteOffset;
    int siteLength;
    int flags;
};

/**
 * Patterns are the recognition site plus the flanking bases NEB requires for an overlapping
 * Dam (GmATC) or Dcm (CmCWGG) site. siteOffset is the recognition site inside the pattern.
 * A hit on the opposite strand is tested with the reverse complement of the pattern.
 * Palindromic enzymes are blocked from either end unless the pattern already includes both
 * flanks (Acc65I, BsaHI, BslI), which is the NEB "multiple overlaps required" case.
 *
 * Containing GATC or CCWGG is not enough. BamHI and Sau3AI cut Dam-methylated DNA.
 * DpnI requires Dam methylation. BstNI cuts Dcm-methylated DNA. Those enzymes are absent here.
 * Impaired-only enzymes (BcgI, and SfiI's internal overlap) are absent too.
 *
 * Sources: NEB "Dam and Dcm Methylases of E. coli" (table archived 2021-01-21) and the REBASE
 * isoschizomer lists for enzymes blocked by a complete GATC or CCWGG site (EcoRII, NdeII, and
 * the other GATC cutters that share MboI's block).
 */
const BlockRule BLOCK_RULES[] = {
    {"AlwI", "GGATC", 0, 5, HostMethylationDam},
    {"BclI", "TGATCA", 0, 6, HostMethylationDam},
    {"BsaBI", "GATCNNNATC", 0, 10, HostMethylationDam},
    {"BscFI", "GATC", 0, 4, HostMethylationDam},
    {"Bsp105I", "GATC", 0, 4, HostMethylationDam},
    {"BspDI", "ATCGATC", 0, 6, HostMethylationDam},
    {"BspEI", "TCCGGATC", 0, 6, HostMethylationDam},
    {"BspHI", "TCATGATC", 0, 6, HostMethylationDam},
    {"BspKT6I", "GATC", 0, 4, HostMethylationDam},
    {"BstEIII", "GATC", 0, 4, HostMethylationDam},
    {"BstKTI", "GATC", 0, 4, HostMethylationDam},
    {"CacI", "GATC", 0, 4, HostMethylationDam},
    {"ClaI", "ATCGATC", 0, 6, HostMethylationDam},
    {"CpfAI", "GATC", 0, 4, HostMethylationDam},
    {"CtyI", "GATC", 0, 4, HostMethylationDam},
    {"CviAI", "GATC", 0, 4, HostMethylationDam},
    {"DpnII", "GATC", 0, 4, HostMethylationDam},
    {"Gst1588II", "GATC", 0, 4, HostMethylationDam},
    {"HphI", "GGTGATC", 0, 5, HostMethylationDam},
    {"Hpy188I", "TCNGATC", 0, 5, HostMethylationDam},
    {"Hpy188III", "TCNNGATC", 0, 6, HostMethylationDam},
    {"HpyAIII", "GATC", 0, 4, HostMethylationDam},
    {"LlaAI", "GATC", 0, 4, HostMethylationDam},
    {"MboI", "GATC", 0, 4, HostMethylationDam},
    {"MboII", "GAAGATC", 0, 5, HostMethylationDam},
    {"Mel2TI", "GATC", 0, 4, HostMethylationDam},
    {"Mel3JI", "GATC", 0, 4, HostMethylationDam},
    {"Mel4OI", "GATC", 0, 4, HostMethylationDam},
    {"Mel5JI", "GATC", 0, 4, HostMethylationDam},
    {"Mel5OI", "GATC", 0, 4, HostMethylationDam},
    {"Mel5TI", "GATC", 0, 4, HostMethylationDam},
    {"Mel7JI", "GATC", 0, 4, HostMethylationDam},
    {"MmeII", "GATC", 0, 4, HostMethylationDam},
    {"NciAI", "GATC", 0, 4, HostMethylationDam},
    {"NdeII", "GATC", 0, 4, HostMethylationDam},
    {"NruI", "TCGCGATC", 0, 6, HostMethylationDam},
    {"RalF40I", "GATC", 0, 4, HostMethylationDam},
    {"TaqI", "TCGATC", 0, 4, HostMethylationDam},
    {"TrsKTI", "GATC", 0, 4, HostMethylationDam},
    {"TrsSI", "GATC", 0, 4, HostMethylationDam},
    {"TrsTI", "GATC", 0, 4, HostMethylationDam},
    {"XbaI", "TCTAGATC", 0, 6, HostMethylationDam},

    {"Acc65I", "CCWGGTACCWGG", 3, 6, HostMethylationDcm},
    {"AlwNI", "CAGNNCCTGG", 0, 9, HostMethylationDcm},
    {"ApaI", "GGGCCCWGG", 0, 6, HostMethylationDcm},
    {"AvaII", "GGWCCWGG", 0, 5, HostMethylationDcm},
    {"BanI", "GGYRCCWGG", 0, 6, HostMethylationDcm},
    {"BsaHI", "CCWGGCGCCWGG", 3, 6, HostMethylationDcm},
    {"BsaI", "GAGACCWGG", 0, 6, HostMethylationDcm},
    {"BslI", "CCWGGNCCWGG", 0, 11, HostMethylationDcm},
    {"BsmFI", "GGGACT", 0, 5, HostMethylationDcm},
    {"BssKI", "CCWGG", 0, 5, HostMethylationDcm},
    {"BstXI", "CCAGGNNNNTGG", 0, 12, HostMethylationDcm},
    {"CbrI", "CCWGG", 0, 5, HostMethylationDcm},
    {"EaeI", "YGGCCAGG", 0, 6, HostMethylationDcm},
    {"EcoO109I", "RGGNCCTGG", 0, 7, HostMethylationDcm},
    {"EcoRII", "CCWGG", 0, 5, HostMethylationDcm},
    {"FseI", "GGCCGGCCWGG", 0, 8, HostMethylationDcm},
    {"MscI", "TGGCCAGG", 0, 6, HostMethylationDcm},
    {"NlaIV", "GGNNCCWGG", 0, 6, HostMethylationDcm},
    {"PflMI", "CCAGGNNNTGG", 0, 11, HostMethylationDcm},
    {"PpuMI", "RGGWCCTGG", 0, 7, HostMethylationDcm},
    {"PspGI", "CCWGG", 0, 5, HostMethylationDcm},
    {"PspOMI", "GGGCCCWGG", 0, 6, HostMethylationDcm},
    {"Sau96I", "GGNCCWGG", 0, 5, HostMethylationDcm},
    {"ScrFI", "CCWGG", 0, 5, HostMethylationDcm},
    {"SexAI", "ACCWGGT", 0, 7, HostMethylationDcm},
    {"SfiI", "GGCCNNNNNGGCCWGG", 0, 13, HostMethylationDcm},
    {"SfoI", "GGCGCCWGG", 0, 6, HostMethylationDcm},
    {"StuI", "AGGCCTGG", 0, 6, HostMethylationDcm},
};

const QHash<QString, QVector<BlockRule>>& rulesByEnzyme() {
    static const QHash<QString, QVector<BlockRule>> rules = []() {
        QHash<QString, QVector<BlockRule>> map;
        for (const BlockRule& rule : BLOCK_RULES) {
            map[QString::fromLatin1(rule.enzyme)].append(rule);
        }
        return map;
    }();
    return rules;
}

char normalizedBase(char base) {
    const char upper = (base >= 'a' && base <= 'z') ? char(base - 'a' + 'A') : base;
    return upper == 'U' ? 'T' : upper;
}

char complementBase(char base) {
    switch (normalizedBase(base)) {
        case 'A':
            return 'T';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'T':
            return 'A';
        case 'R':
            return 'Y';
        case 'Y':
            return 'R';
        case 'M':
            return 'K';
        case 'K':
            return 'M';
        case 'B':
            return 'V';
        case 'D':
            return 'H';
        case 'H':
            return 'D';
        case 'V':
            return 'B';
        case 'W':
            return 'W';
        case 'S':
            return 'S';
        case 'N':
            return 'N';
        default:
            return 'N';
    }
}

QByteArray reverseComplementPattern(const char* pattern) {
    const int length = int(std::strlen(pattern));
    QByteArray reversed(length, '\0');
    for (int i = 0; i < length; ++i) {
        reversed[length - 1 - i] = complementBase(pattern[i]);
    }
    return reversed;
}

bool iupacMatch(char sequenceBase, char patternBase) {
    const char base = normalizedBase(sequenceBase);
    const char pattern = normalizedBase(patternBase);
    switch (pattern) {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
            return base == pattern;
        case 'R':
            return base == 'A' || base == 'G';
        case 'Y':
            return base == 'C' || base == 'T';
        case 'W':
            return base == 'A' || base == 'T';
        case 'S':
            return base == 'C' || base == 'G';
        case 'M':
            return base == 'A' || base == 'C';
        case 'K':
            return base == 'G' || base == 'T';
        case 'B':
            return base == 'C' || base == 'G' || base == 'T';
        case 'D':
            return base == 'A' || base == 'G' || base == 'T';
        case 'H':
            return base == 'A' || base == 'C' || base == 'T';
        case 'V':
            return base == 'A' || base == 'C' || base == 'G';
        case 'N':
            return base == 'A' || base == 'C' || base == 'G' || base == 'T';
        default:
            return false;
    }
}

bool baseAtSiteOffset(const QByteArray& sequence, bool circular, const QVector<U2Region>& regions, qint64 siteLength, qint64 offset, char* base) {
    if (regions.isEmpty() || siteLength <= 0) {
        return false;
    }
    qint64 pos = 0;
    if (offset < 0) {
        pos = regions.first().startPos + offset;
    } else if (offset >= siteLength) {
        pos = regions.last().endPos() + (offset - siteLength);
    } else {
        qint64 remaining = offset;
        pos = regions.last().endPos();
        for (const U2Region& region : regions) {
            if (remaining < region.length) {
                pos = region.startPos + remaining;
                break;
            }
            remaining -= region.length;
        }
    }

    const qint64 length = sequence.size();
    if (circular) {
        if (length <= 0) {
            return false;
        }
        pos %= length;
        if (pos < 0) {
            pos += length;
        }
    } else if (pos < 0 || pos >= length) {
        return false;
    }
    *base = sequence.at(int(pos));
    return true;
}

bool patternMatches(const QByteArray& sequence, bool circular, const QVector<U2Region>& regions, qint64 siteLength, const QByteArray& pattern, int siteOffset) {
    if (siteOffset < 0 || siteOffset + siteLength > pattern.size()) {
        return false;
    }
    for (int i = 0; i < pattern.size(); ++i) {
        char base = 0;
        if (!baseAtSiteOffset(sequence, circular, regions, siteLength, qint64(i - siteOffset), &base)) {
            return false;
        }
        if (!iupacMatch(base, pattern.at(i))) {
            return false;
        }
    }
    return true;
}

}  // namespace

bool enzymeMayBeBlockedByHostMethylation(const QString& enzymeName) {
    return rulesByEnzyme().contains(enzymeName);
}

int hostMethylationBlock(const QString& enzymeName, const QByteArray& sequence, bool circular, const QVector<U2Region>& regions) {
    if (regions.isEmpty() || sequence.isEmpty()) {
        return HostMethylationNone;
    }
    qint64 siteLength = 0;
    for (const U2Region& region : regions) {
        siteLength += region.length;
    }
    if (siteLength <= 0) {
        return HostMethylationNone;
    }

    const QVector<BlockRule> rules = rulesByEnzyme().value(enzymeName);
    int flags = HostMethylationNone;
    for (const BlockRule& rule : rules) {
        if (rule.siteLength != siteLength) {
            continue;
        }
        const QByteArray pattern = QByteArray(rule.pattern);
        const bool direct = patternMatches(sequence, circular, regions, siteLength, pattern, rule.siteOffset);
        const QByteArray reversed = reverseComplementPattern(rule.pattern);
        const int reversedOffset = pattern.size() - rule.siteOffset - rule.siteLength;
        const bool complement = patternMatches(sequence, circular, regions, siteLength, reversed, reversedOffset);
        if (direct || complement) {
            flags |= rule.flags;
        }
    }
    return flags;
}

}  // namespace U2
